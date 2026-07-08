// Exact GF(p) kernel for the observability analytic path of symmetryDetection().
// Builds the observability-identifiability matrix by a Taylor-mode construction
// over a finite field (Lie derivatives as truncated power series in time,
// parameter gradients as forward-mode dual numbers), reduces it, and supports
// the rational reconstruction used by the R orchestration. Primes are < 2^31 so
// residue products stay below 2^62 in uint64; the CRT product of four primes
// stays below 2^124 in unsigned __int128.

#include <Rcpp.h>

#include <vector>
#include <string>
#include <cstdint>
#include <algorithm>
#include <map>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

typedef uint64_t u64;
// __int128 is a GCC/Clang extension (not ISO C++), used here for a 128-bit
// accumulator in modular multiply / CRT. __extension__ marks the use as a
// deliberate extension so -Wpedantic stays quiet without blanket-suppressing it.
__extension__ typedef __int128 i128;
__extension__ typedef unsigned __int128 u128;

namespace {

const int OP_CONST = 0, OP_ADD = 1, OP_MUL = 2, OP_INV = 3;

inline u64 addmod(u64 a, u64 b, u64 p) { u64 s = a + b; return s >= p ? s - p : s; }
inline u64 submod(u64 a, u64 b, u64 p) { return a >= b ? a - b : a + p - b; }
inline u64 mulmod(u64 a, u64 b, u64 p) { return (u64)((u128)a * b % p); }

u64 powmod(u64 a, u64 e, u64 p) {
  u64 r = 1;
  a %= p;
  while (e) {
    if (e & 1) r = mulmod(r, a, p);
    a = mulmod(a, a, p);
    e >>= 1;
  }
  return r;
}

inline u64 invmod(u64 a, u64 p) { return powmod(a % p, p - 2, p); }

// Reduce a (possibly negative) R integer to its residue in [0, p).
inline u64 red(long long x, u64 p) {
  i128 v = (i128)x % (i128)p; if (v < 0) v += p; return (u64)v;
}

inline std::vector<int> to_ivec(const IntegerVector& v) {
  return std::vector<int>(v.begin(), v.end());
}

inline std::vector<std::string> to_svec(const CharacterVector& v) {
  std::vector<std::string> o(v.size());
  for (int i = 0; i < v.size(); ++i) o[i] = as<std::string>(v[i]);
  return o;
}

// Parse a signed decimal string to its residue modulo p (handles arbitrary
// length without overflow).
u64 parse_mod(const std::string& s, u64 p) {
  size_t i = 0;
  bool neg = false;
  if (i < s.size() && (s[i] == '-' || s[i] == '+')) { neg = s[i] == '-'; ++i; }
  u64 r = 0;
  for (; i < s.size(); ++i) {
    if (s[i] < '0' || s[i] > '9') continue;
    r = (r * 10 + (u64)(s[i] - '0')) % p;
  }
  return neg ? (p - r) % p : r;
}

// Residue of the rational constant num/den (decimal strings) modulo p.
inline u64 reduce_rational(const std::string& num, const std::string& den, u64 p) {
  return mulmod(parse_mod(num, p), invmod(parse_mod(den, p), p), p);
}

// Accumulate the dual product of duals X and Y into o (width w: value at 0,
// partial derivatives at 1..w-1).
inline void dual_mul_acc(u64* o, const u64* X, const u64* Y, int w, u64 p) {
  u64 xv = X[0], yv = Y[0];
  o[0] = addmod(o[0], mulmod(xv, yv, p), p);
  for (int c = 1; c < w; ++c) {
    u64 t = addmod(mulmod(xv, Y[c], p), mulmod(yv, X[c], p), p);
    o[c] = addmod(o[c], t, p);
  }
}

// ---- Multivariate truncated power series over K formal gap lengths ----
// A "polyseries" is a flat array of nMono monomial coefficients, each a width-w
// dual number [value, d/dz_1, ...]. Monomials are the multi-indices mu in N^K with
// total degree |mu| <= Mtot, ordered by increasing total degree. It carries the
// exact dependence of a propagated state on the (generic) inter-event gap lengths
// Delta t_i; the rank over GF(p) later treats each monomial as a separate row (a
// polynomial in the gaps vanishes for generic gaps iff every coefficient does).
struct PolyBasis {
  int K, Mtot, nMono, w;
  std::vector<std::vector<int> > idx;          // nMono multi-indices, degree-ordered
  std::vector<std::vector<int> > sumSlot;      // sumSlot[a][b] = slot(idx[a]+idx[b]) or -1
  std::vector<std::vector<std::pair<int,int> > > decomp;  // decomp[s] = {(a,b): a+b == s}
  std::map<std::vector<int>, int> slotOf;

  // slot of idx[a] with `axis` incremented by k (promoting a time-order into a
  // gap axis), or -1 if it overflows the total-degree budget
  int shiftSlot(int a, int axis, int k) const {
    std::vector<int> mu = idx[a];
    mu[axis] += k;
    int td = 0; for (int i = 0; i < K; ++i) td += mu[i];
    if (td > Mtot) return -1;
    std::map<std::vector<int>, int>::const_iterator it = slotOf.find(mu);
    return it == slotOf.end() ? -1 : it->second;
  }

  PolyBasis(int K_, int Mtot_, int w_) : K(K_), Mtot(Mtot_), w(w_) {
    std::vector<int> mu(K, 0);
    for (int d = 0; d <= Mtot; ++d) enumerate(d, 0, d, mu);
    nMono = (int)idx.size();
    for (int s = 0; s < nMono; ++s) slotOf[idx[s]] = s;
    sumSlot.assign(nMono, std::vector<int>(nMono, -1));
    decomp.assign(nMono, std::vector<std::pair<int,int> >());
    std::vector<int> sum(K);
    for (int a = 0; a < nMono; ++a)
      for (int b = 0; b < nMono; ++b) {
        int td = 0;
        for (int i = 0; i < K; ++i) { sum[i] = idx[a][i] + idx[b][i]; td += sum[i]; }
        if (td > Mtot) continue;
        int s = slotOf[sum];
        sumSlot[a][b] = s;
        decomp[s].push_back(std::make_pair(a, b));
      }
  }
  // recursively emit all K-compositions of `remaining` (total degree d) at axis pos
  void enumerate(int remaining, int pos, int d, std::vector<int>& mu) {
    if (pos == K - 1) { mu[pos] = remaining; idx.push_back(mu); return; }
    for (int v = 0; v <= remaining; ++v) {
      mu[pos] = v;
      enumerate(remaining - v, pos + 1, d, mu);
    }
  }
};

// out (nMono x w) += A * B, convolution over monomials with the dual product rule
inline void poly_mul_acc(u64* out, const u64* A, const u64* B,
                         const PolyBasis& pb, u64 p) {
  int w = pb.w;
  for (int a = 0; a < pb.nMono; ++a) {
    const u64* Aa = A + (size_t)a * w;
    bool nz = false;
    for (int c = 0; c < w; ++c) if (Aa[c]) { nz = true; break; }
    if (!nz) continue;
    const std::vector<int>& srow = pb.sumSlot[a];
    for (int b = 0; b < pb.nMono; ++b) {
      int s = srow[b];
      if (s < 0) continue;
      dual_mul_acc(out + (size_t)s * w, Aa, B + (size_t)b * w, w, p);
    }
  }
}

// B = 1 / A as a polyseries of duals (A's constant monomial must have nonzero
// value). Computed monomial by monomial in increasing total degree.
inline bool poly_inv(u64* B, const u64* A, const PolyBasis& pb, u64 p) {
  int w = pb.w;
  std::fill(B, B + (size_t)pb.nMono * w, (u64)0);
  u64 a0 = A[0] % p;                       // slot 0 is the degree-0 monomial
  if (a0 == 0) return false;
  u64 vi = invmod(a0, p), vi2 = mulmod(vi, vi, p);
  std::vector<u64> invA0(w, 0);
  invA0[0] = vi;
  for (int c = 1; c < w; ++c) invA0[c] = submod(0, mulmod(A[c], vi2, p), p);
  for (int c = 0; c < w; ++c) B[c] = invA0[c];
  std::vector<u64> s(w, 0);
  for (int m = 1; m < pb.nMono; ++m) {     // monomials are degree-ordered
    std::fill(s.begin(), s.end(), (u64)0);
    const std::vector<std::pair<int,int> >& dc = pb.decomp[m];
    for (size_t t = 0; t < dc.size(); ++t) {
      int a = dc[t].first, b = dc[t].second;
      if (a == 0) continue;                // exclude the A[0]*B[m] term
      dual_mul_acc(s.data(), A + (size_t)a * w, B + (size_t)b * w, w, p);
    }
    for (int c = 0; c < w; ++c) s[c] = submod(0, s[c], p);   // negate
    dual_mul_acc(B + (size_t)m * w, invA0.data(), s.data(), w, p);
  }
  return true;
}

// In-place Gauss-Jordan over GF(p); returns the pivot columns.
std::vector<int> rref_mod(std::vector<std::vector<u64> >& A, u64 p) {
  std::vector<int> pivots;
  if (A.empty()) return pivots;
  int nrows = (int)A.size(), ncols = (int)A[0].size(), r = 0;
  for (int c = 0; c < ncols && r < nrows; ++c) {
    int piv = -1;
    for (int i = r; i < nrows; ++i)
      if (A[i][c] % p != 0) { piv = i; break; }
    if (piv < 0) continue;
    std::swap(A[r], A[piv]);
    u64 inv = invmod(A[r][c], p);
    for (int j = 0; j < ncols; ++j) A[r][j] = mulmod(A[r][j], inv, p);
    for (int i = 0; i < nrows; ++i) {
      if (i == r) continue;
      u64 f = A[i][c] % p;
      if (f == 0) continue;
      for (int j = 0; j < ncols; ++j)
        A[i][j] = submod(A[i][j], mulmod(f, A[r][j], p), p);
    }
    pivots.push_back(c);
    ++r;
  }
  return pivots;
}

u128 u128_isqrt(u128 n) {
  if (n == 0) return 0;
  int bits = 0;
  for (u128 t = n; t; t >>= 1) ++bits;
  u128 x = (u128)1 << ((bits + 1) / 2);
  for (int it = 0; it < 300; ++it) {
    u128 y = (x + n / x) / 2;
    if (y >= x) break;
    x = y;
  }
  while (x > 0 && x * x > n) --x;
  while ((x + 1) * (x + 1) <= n) ++x;
  return x;
}

std::string i128_to_string(i128 v) {
  if (v == 0) return "0";
  bool neg = v < 0;
  u128 u = neg ? (u128)(-(v + 1)) + 1 : (u128)v;
  std::string s;
  while (u > 0) { s += (char)('0' + (int)(u % 10)); u /= 10; }
  if (neg) s += '-';
  std::reverse(s.begin(), s.end());
  return s;
}

// Recover n/d with x = n * d^{-1} (mod M), |n|, |d| <= sqrt(M/2); returns false
// if no such bounded rational exists.
bool rational_reconstruct(u128 x, u128 M, i128& num, i128& den) {
  x %= M;
  if (x == 0) { num = 0; den = 1; return true; }
  u128 bound = u128_isqrt(M / 2);
  i128 r0 = (i128)M, r1 = (i128)x, s0 = 0, s1 = 1;
  while (r1 > (i128)bound) {
    i128 q = r0 / r1;
    i128 r2 = r0 - q * r1; r0 = r1; r1 = r2;
    i128 s2 = s0 - q * s1; s0 = s1; s1 = s2;
  }
  if (s1 == 0) return false;
  num = r1; den = s1;
  if (den < 0) { num = -num; den = -den; }
  u128 an = (u128)(num < 0 ? -num : num);
  if (an > bound || (u128)den > bound) return false;
  return true;
}

// Evaluate the initial-condition tape at time order 0 only. icVal has one width-w
// dual vector per slot; leaves [0, nLeaves) are pre-filled with their value and
// dual, instruction i writes slot nLeaves + i. Returns false on a zero reciprocal.
bool eval_ic_order0(std::vector<std::vector<u64> >& icVal,
                    const std::vector<int>& icOp, const std::vector<int>& icA,
                    const std::vector<int>& icB, const std::vector<u64>& icCval,
                    int nLeaves, int w, u64 p) {
  int n = (int)icOp.size();
  for (int i = 0; i < n; ++i) {
    u64* o = icVal[nLeaves + i].data();
    switch (icOp[i]) {
      case OP_CONST:
        o[0] = icCval[i];
        break;
      case OP_ADD: {
        const u64* pa = icVal[icA[i]].data();
        const u64* pb = icVal[icB[i]].data();
        for (int c = 0; c < w; ++c) o[c] = addmod(pa[c], pb[c], p);
        break;
      }
      case OP_MUL: {
        const u64* pa = icVal[icA[i]].data();
        const u64* pb = icVal[icB[i]].data();
        for (int c = 0; c < w; ++c) o[c] = 0;
        dual_mul_acc(o, pa, pb, w, p);
        break;
      }
      case OP_INV: {
        const u64* pa = icVal[icA[i]].data();
        u64 a0 = pa[0] % p;
        if (a0 == 0) return false;
        u64 vi = invmod(a0, p), vi2 = mulmod(vi, vi, p);
        o[0] = vi;
        for (int c = 1; c < w; ++c) o[c] = submod(0, mulmod(pa[c], vi2, p), p);
        break;
      }
    }
  }
  return true;
}

bool build_obs_rows(std::vector<std::vector<u64> >& val,
                    const std::vector<int>& op, const std::vector<int>& a,
                    const std::vector<int>& b, const std::vector<u64>& cval,
                    int instrBase, const std::vector<int>& stateSlots,
                    const std::vector<int>& fOut, const std::vector<int>& gOut,
                    int nz, int w, int Nt, u64 p,
                    std::vector<std::vector<u64> >& rows) {
  int nInstr = (int)op.size();
  int m = (int)stateSlots.size();
  std::vector<std::vector<u64> > invA0(nInstr);
  bool fail = false;
  for (int k = 0; k <= Nt && !fail; ++k) {
    for (int i = 0; i < nInstr; ++i) {
      u64* o = val[instrBase + i].data() + (size_t)k * w;
      switch (op[i]) {
        case OP_CONST:
          if (k == 0) o[0] = cval[i];
          break;
        case OP_ADD: {
          const u64* pa = val[a[i]].data() + (size_t)k * w;
          const u64* pb = val[b[i]].data() + (size_t)k * w;
          for (int c = 0; c < w; ++c) o[c] = addmod(pa[c], pb[c], p);
          break;
        }
        case OP_MUL: {
          const u64* pa = val[a[i]].data();
          const u64* pb = val[b[i]].data();
          for (int ii = 0; ii <= k; ++ii)
            dual_mul_acc(o, pa + (size_t)ii * w, pb + (size_t)(k - ii) * w, w, p);
          break;
        }
        case OP_INV: {
          const u64* pa = val[a[i]].data();
          if (k == 0) {
            u64 a0 = pa[0] % p;
            if (a0 == 0) { fail = true; break; }
            u64 vi = invmod(a0, p), vi2 = mulmod(vi, vi, p);
            invA0[i].assign(w, 0);
            invA0[i][0] = vi;
            for (int c = 1; c < w; ++c) invA0[i][c] = submod(0, mulmod(pa[c], vi2, p), p);
            for (int c = 0; c < w; ++c) o[c] = invA0[i][c];
          } else {
            std::vector<u64> s(w, 0);
            const u64* self = val[instrBase + i].data();
            for (int j = 1; j <= k; ++j)
              dual_mul_acc(s.data(), pa + (size_t)j * w,
                           self + (size_t)(k - j) * w, w, p);
            for (int c = 0; c < w; ++c) s[c] = submod(0, s[c], p);
            dual_mul_acc(o, invA0[i].data(), s.data(), w, p);
          }
          break;
        }
      }
      if (fail) break;
    }
    if (fail || k == Nt) break;
    u64 invk = invmod((u64)(k + 1), p);
    for (int i = 0; i < m; ++i) {
      const u64* fk = val[fOut[i]].data() + (size_t)k * w;
      u64* st = val[stateSlots[i]].data() + (size_t)(k + 1) * w;
      for (int c = 0; c < w; ++c) st[c] = mulmod(fk[c], invk, p);
    }
  }
  if (fail) return false;
  for (int gi = 0; gi < (int)gOut.size(); ++gi) {
    const u64* g = val[gOut[gi]].data();
    for (int k = 0; k <= Nt; ++k) {
      std::vector<u64> row(nz);
      for (int c = 0; c < nz; ++c) row[c] = g[(size_t)k * w + 1 + c];
      rows.push_back(row);
    }
  }
  return true;
}

// Polyseries generalisation of build_obs_rows: every value is a width-(nMono*w)
// block (a truncated power series in the gap lengths, each monomial a dual). The
// time-Taylor recurrence mirrors build_obs_rows in structure; the scalar dual
// operations become their polyseries versions. Computes the state and output coefficients
// up to time order Nrun, emitting z-gradient rows for time orders 0..Ntemit and all
// gap monomials. `val` slots must be pre-sized to (Nrun+1)*blk and the state slots
// pre-seeded at time order 0.
bool build_obs_rows_poly(std::vector<std::vector<u64> >& val,
                         const std::vector<int>& op, const std::vector<int>& a,
                         const std::vector<int>& b, const std::vector<u64>& cval,
                         int instrBase, const std::vector<int>& stateSlots,
                         const std::vector<int>& fOut, const std::vector<int>& gOut,
                         const PolyBasis& pb, int Nrun, int Ntemit, u64 p,
                         std::vector<std::vector<u64> >& rows) {
  int nInstr = (int)op.size();
  int m = (int)stateSlots.size();
  int w = pb.w, nM = pb.nMono, blk = nM * w, nz = w - 1;
  std::vector<std::vector<u64> > invA0(nInstr);
  bool fail = false;
  for (int k = 0; k <= Nrun && !fail; ++k) {
    for (int i = 0; i < nInstr; ++i) {
      u64* o = val[instrBase + i].data() + (size_t)k * blk;
      switch (op[i]) {
        case OP_CONST:
          if (k == 0) o[0] = cval[i];       // value at monomial 0, dual 0
          break;
        case OP_ADD: {
          const u64* pa = val[a[i]].data() + (size_t)k * blk;
          const u64* pbb = val[b[i]].data() + (size_t)k * blk;
          for (int c = 0; c < blk; ++c) o[c] = addmod(pa[c], pbb[c], p);
          break;
        }
        case OP_MUL: {
          const u64* pa = val[a[i]].data();
          const u64* pbb = val[b[i]].data();
          for (int ii = 0; ii <= k; ++ii)
            poly_mul_acc(o, pa + (size_t)ii * blk, pbb + (size_t)(k - ii) * blk, pb, p);
          break;
        }
        case OP_INV: {
          const u64* pa = val[a[i]].data();
          if (k == 0) {
            invA0[i].assign(blk, 0);
            if (!poly_inv(invA0[i].data(), pa, pb, p)) { fail = true; break; }
            for (int c = 0; c < blk; ++c) o[c] = invA0[i][c];
          } else {
            std::vector<u64> s(blk, 0);
            const u64* self = val[instrBase + i].data();
            for (int j = 1; j <= k; ++j)
              poly_mul_acc(s.data(), pa + (size_t)j * blk,
                           self + (size_t)(k - j) * blk, pb, p);
            for (int c = 0; c < blk; ++c) s[c] = submod(0, s[c], p);
            poly_mul_acc(o, invA0[i].data(), s.data(), pb, p);
          }
          break;
        }
      }
      if (fail) break;
    }
    if (fail || k == Nrun) break;
    u64 invk = invmod((u64)(k + 1), p);
    for (int i = 0; i < m; ++i) {
      const u64* fk = val[fOut[i]].data() + (size_t)k * blk;
      u64* st = val[stateSlots[i]].data() + (size_t)(k + 1) * blk;
      for (int c = 0; c < blk; ++c) st[c] = mulmod(fk[c], invk, p);
    }
  }
  if (fail) return false;
  for (int gi = 0; gi < (int)gOut.size(); ++gi) {
    const u64* g = val[gOut[gi]].data();
    for (int k = 0; k <= Ntemit; ++k)
      for (int mo = 0; mo < nM; ++mo) {
        std::vector<u64> row(nz);
        const u64* d = g + (size_t)k * blk + (size_t)mo * w;
        for (int c = 0; c < nz; ++c) row[c] = d[1 + c];
        rows.push_back(row);
      }
  }
  return true;
}

// One experimental condition's tape in thread-safe plain C++: the instruction
// stream and slot wiring as ints, the rational constants of the tape and the
// initial-condition map kept as decimal strings so they can be reduced against
// any prime, and an optional steady-state seed icSeed (constraint mode). Every
// numeric constant is held raw (strings, or unreduced ints for icSeed) and
// reduced per prime in build_one_condition, so a single extraction serves every
// prime -- which the batch path relies on.
struct CondRaw {
  std::vector<int> op, a, b, stateSlots, fOut, gOut, icLeaf, icOp, icA, icB, icOut;
  std::vector<std::string> cnum, cden, icCnum, icCden, icNum, icDen;
  bool hasIcTape, hasIcSeed;
  std::vector<std::vector<int> > icSeedRaw;
};

CondRaw extract_cond_raw(List tp, int nStates, int w) {
  CondRaw cd;
  cd.op = to_ivec(tp["op"]); cd.a = to_ivec(tp["a"]); cd.b = to_ivec(tp["b"]);
  cd.stateSlots = to_ivec(tp["stateSlots"]);
  cd.fOut = to_ivec(tp["fOut"]); cd.gOut = to_ivec(tp["gOut"]);
  cd.icLeaf = to_ivec(tp["icLeaf"]);
  cd.cnum = to_svec(tp["cnum"]); cd.cden = to_svec(tp["cden"]);
  cd.hasIcTape = tp.containsElementNamed("icOp");
  cd.hasIcSeed = tp.containsElementNamed("icSeed");
  if (cd.hasIcSeed) {
    IntegerMatrix icSeed = as<IntegerMatrix>(tp["icSeed"]);
    cd.icSeedRaw.assign(nStates, std::vector<int>(w, 0));
    for (int i = 0; i < nStates; ++i)
      for (int c = 0; c < w; ++c) cd.icSeedRaw[i][c] = icSeed(i, c);
  } else if (cd.hasIcTape) {
    cd.icOp = to_ivec(tp["icOp"]); cd.icA = to_ivec(tp["icA"]);
    cd.icB = to_ivec(tp["icB"]); cd.icOut = to_ivec(tp["icOut"]);
    cd.icCnum = to_svec(tp["icCnum"]); cd.icCden = to_svec(tp["icCden"]);
  } else {
    cd.icNum = to_svec(tp["icNum"]); cd.icDen = to_svec(tp["icDen"]);
  }
  return cd;
}

// Build one condition's observability rows at a single point/prime. leafVal holds
// each leaf's value and dual already reduced mod p; the tape and IC constants are
// reduced from their strings here, so the same CondRaw serves any prime. Returns
// false on a vanishing denominator (tape reciprocal or reciprocal IC).
bool build_one_condition(const CondRaw& cd,
                         const std::vector<std::vector<u64> >& leafVal,
                         int nLeaves, int nStates, int nz, int w, int Nt, u64 p,
                         std::vector<std::vector<u64> >& outRows) {
  int nInstr = (int)cd.op.size();
  int instrBase = nLeaves + nStates;
  int S = instrBase + nInstr;

  std::vector<u64> cval(nInstr, 0);
  for (int i = 0; i < nInstr; ++i)
    if (cd.op[i] == OP_CONST)
      cval[i] = reduce_rational(cd.cnum[i], cd.cden[i], p);

  std::vector<std::vector<u64> > icVal;
  if (cd.hasIcTape && !cd.hasIcSeed) {
    int nIc = (int)cd.icOp.size();
    std::vector<u64> icCval(nIc, 0);
    for (int i = 0; i < nIc; ++i)
      if (cd.icOp[i] == OP_CONST)
        icCval[i] = reduce_rational(cd.icCnum[i], cd.icCden[i], p);
    icVal.assign(nLeaves + nIc, std::vector<u64>(w, 0));
    for (int L = 0; L < nLeaves; ++L)
      for (int c = 0; c < w; ++c) icVal[L][c] = leafVal[L][c];
    if (!eval_ic_order0(icVal, cd.icOp, cd.icA, cd.icB, icCval, nLeaves, w, p))
      return false;
  }

  std::vector<u64> icConst;
  if (!cd.hasIcTape && !cd.hasIcSeed) {
    icConst.assign(nStates, 0);
    for (int i = 0; i < nStates; ++i)
      icConst[i] = reduce_rational(cd.icNum[i], cd.icDen[i], p);
  }

  std::vector<std::vector<u64> > val(S, std::vector<u64>((size_t)(Nt + 1) * w, 0));
  for (int L = 0; L < nLeaves; ++L)
    for (int c = 0; c < w; ++c) val[L][c] = leafVal[L][c];
  for (int i = 0; i < (int)cd.stateSlots.size(); ++i) {
    int slot = cd.stateSlots[i];
    if (cd.hasIcSeed) {
      for (int c = 0; c < w; ++c) val[slot][c] = red(cd.icSeedRaw[i][c], p);
      continue;
    }
    int src = cd.hasIcTape ? cd.icOut[i] : (cd.icLeaf[i] >= 0 ? cd.icLeaf[i] : -1);
    if (src >= 0) {
      const std::vector<u64>& s = src < nLeaves ? leafVal[src] : icVal[src];
      for (int c = 0; c < w; ++c) val[slot][c] = s[c];
    } else {
      val[slot][0] = icConst[i];
    }
  }

  return build_obs_rows(val, cd.op, cd.a, cd.b, cval, instrBase, cd.stateSlots,
                        cd.fOut, cd.gOut, nz, w, Nt, p, outRows);
}

// One chain segment in thread-safe plain C++: the regime-substituted tape, its
// reduced constants, and the optional first-segment steady-state seed, IC tape
// and state-dose event map. Constants are reduced against the single call prime.
struct SegRaw {
  std::vector<int> op, a, b, stateSlots, fOut, gOut;
  std::vector<u64> cval;
  bool hasIcSeed, hasIcTape, hasEv;
  std::vector<std::vector<u64> > icSeed;
  std::vector<int> icOp, icA, icB, icOut;
  std::vector<u64> icCval;
  std::vector<int> evVarIdx, evMethod, evOp, evA, evB, evOut;
  std::vector<u64> evCval;
};

SegRaw extract_seg_raw(List seg, int nStates, int w, u64 p) {
  SegRaw s;
  s.op = to_ivec(seg["op"]); s.a = to_ivec(seg["a"]); s.b = to_ivec(seg["b"]);
  s.stateSlots = to_ivec(seg["stateSlots"]);
  s.fOut = to_ivec(seg["fOut"]); s.gOut = to_ivec(seg["gOut"]);
  std::vector<std::string> cnum = to_svec(seg["cnum"]), cden = to_svec(seg["cden"]);
  int nInstr = (int)s.op.size();
  s.cval.assign(nInstr, 0);
  for (int i = 0; i < nInstr; ++i)
    if (s.op[i] == OP_CONST)
      s.cval[i] = reduce_rational(cnum[i], cden[i], p);
  s.hasIcSeed = seg.containsElementNamed("icSeed");
  if (s.hasIcSeed) {
    IntegerMatrix icSeed = as<IntegerMatrix>(seg["icSeed"]);
    s.icSeed.assign(nStates, std::vector<u64>(w, 0));
    for (int i = 0; i < nStates; ++i)
      for (int c = 0; c < w; ++c) s.icSeed[i][c] = red(icSeed(i, c), p);
  }
  s.hasIcTape = seg.containsElementNamed("icOp");
  if (s.hasIcTape) {
    s.icOp = to_ivec(seg["icOp"]); s.icA = to_ivec(seg["icA"]);
    s.icB = to_ivec(seg["icB"]); s.icOut = to_ivec(seg["icOut"]);
    std::vector<std::string> icCnum = to_svec(seg["icCnum"]), icCden = to_svec(seg["icCden"]);
    int nIc = (int)s.icOp.size();
    s.icCval.assign(nIc, 0);
    for (int i = 0; i < nIc; ++i)
      if (s.icOp[i] == OP_CONST)
        s.icCval[i] = reduce_rational(icCnum[i], icCden[i], p);
  }
  s.hasEv = seg.containsElementNamed("evVarIdx");
  if (s.hasEv) {
    s.evVarIdx = to_ivec(seg["evVarIdx"]); s.evMethod = to_ivec(seg["evMethod"]);
    s.evOp = to_ivec(seg["evOp"]); s.evA = to_ivec(seg["evA"]); s.evB = to_ivec(seg["evB"]);
    s.evOut = to_ivec(seg["evOut"]);
    std::vector<std::string> ecn = to_svec(seg["evCnum"]), ecd = to_svec(seg["evCden"]);
    int nEv = (int)s.evOp.size();
    s.evCval.assign(nEv, 0);
    for (int i = 0; i < nEv; ++i)
      if (s.evOp[i] == OP_CONST)
        s.evCval[i] = reduce_rational(ecn[i], ecd[i], p);
  }
  return s;
}

// Build one condition's chain (all segments) at a single point/prime, propagating
// the state across each gap as a formal power series and applying state-dose
// events on the carry. leafPt holds the leaf residues; the z-leaves seed a dual
// unit. Returns false on a vanishing denominator. Touches no R object.
bool build_one_chain(const std::vector<SegRaw>& segs, int nLeaves, int nStates,
                     int nz, int w, const std::vector<int>& dualCol,
                     const std::vector<int>& leafPt, int Nt, int Mtot, u64 p,
                     std::vector<std::vector<u64> >& rows) {
  int instrBase = nLeaves + nStates;
  int nSeg = (int)segs.size();
  int K = nSeg - 1; if (K < 1) K = 1;
  PolyBasis pb(K, Mtot, w);
  int nM = pb.nMono, blk = nM * w;
  int Nrun = Nt > Mtot ? Nt : Mtot;

  std::vector<std::vector<u64> > leafVal(nLeaves, std::vector<u64>(blk, 0));
  for (int L = 0; L < nLeaves; ++L) {
    i128 pv = (i128)leafPt[L] % (i128)p; if (pv < 0) pv += p;
    leafVal[L][0] = (u64)pv;
    if (dualCol[L] >= 0) leafVal[L][1 + dualCol[L]] = 1;
  }

  std::vector<std::vector<u64> > carry;
  for (int sj = 0; sj < nSeg; ++sj) {
    const SegRaw& seg = segs[sj];
    int nInstr = (int)seg.op.size();
    int Sn = instrBase + nInstr;
    std::vector<std::vector<u64> > val(Sn, std::vector<u64>((size_t)(Nrun + 1) * blk, 0));
    for (int L = 0; L < nLeaves; ++L)
      for (int c = 0; c < blk; ++c) val[L][c] = leafVal[L][c];

    std::vector<std::vector<u64> > icVal;
    bool useIcTape = (sj == 0) && !seg.hasIcSeed && seg.hasIcTape;
    if (useIcTape) {
      int nIc = (int)seg.icOp.size();
      icVal.assign(nLeaves + nIc, std::vector<u64>(w, 0));
      for (int L = 0; L < nLeaves; ++L)
        for (int c = 0; c < w; ++c) icVal[L][c] = leafVal[L][c];
      if (!eval_ic_order0(icVal, seg.icOp, seg.icA, seg.icB, seg.icCval, nLeaves, w, p))
        return false;
    }
    for (int i = 0; i < nStates; ++i) {
      u64* st = val[nLeaves + i].data();
      if (sj == 0 && seg.hasIcSeed) {
        for (int c = 0; c < w; ++c) st[c] = seg.icSeed[i][c];
      } else if (useIcTape) {
        const std::vector<u64>& s = seg.icOut[i] < nLeaves ? leafVal[seg.icOut[i]]
                                                           : icVal[seg.icOut[i]];
        for (int c = 0; c < w; ++c) st[c] = s[c];
      } else {
        for (int c = 0; c < blk; ++c) st[c] = carry[i][c];
      }
    }
    if (!build_obs_rows_poly(val, seg.op, seg.a, seg.b, seg.cval, instrBase,
                             seg.stateSlots, seg.fOut, seg.gOut, pb, Nrun, Nt, p, rows))
      return false;
    if (sj < nSeg - 1) {
      int axis = sj;
      carry.assign(nStates, std::vector<u64>(blk, 0));
      for (int i = 0; i < nStates; ++i) {
        const u64* xv = val[seg.stateSlots[i]].data();
        for (int k = 0; k <= Mtot; ++k)
          for (int aMono = 0; aMono < nM; ++aMono) {
            int dst = pb.shiftSlot(aMono, axis, k);
            if (dst < 0) continue;
            const u64* src = xv + (size_t)k * blk + (size_t)aMono * w;
            u64* d = carry[i].data() + (size_t)dst * w;
            for (int c = 0; c < w; ++c) d[c] = addmod(d[c], src[c], p);
          }
      }
      const SegRaw& nxt = segs[sj + 1];
      if (nxt.hasEv) {
        int nEv = (int)nxt.evOp.size();
        std::vector<std::vector<u64> > evVal(nLeaves + nEv, std::vector<u64>(w, 0));
        for (int L = 0; L < nLeaves; ++L)
          for (int c = 0; c < w; ++c) evVal[L][c] = leafVal[L][c];
        if (!eval_ic_order0(evVal, nxt.evOp, nxt.evA, nxt.evB, nxt.evCval, nLeaves, w, p))
          return false;
        for (int e = 0; e < (int)nxt.evVarIdx.size(); ++e) {
          int s = nxt.evVarIdx[e], meth = nxt.evMethod[e];
          const std::vector<u64>& vv = nxt.evOut[e] < nLeaves ? leafVal[nxt.evOut[e]]
                                                              : evVal[nxt.evOut[e]];
          u64* cv = carry[s].data();
          if (meth == 0) {
            std::fill(carry[s].begin(), carry[s].end(), (u64)0);
            for (int c = 0; c < w; ++c) cv[c] = vv[c] % p;
          } else if (meth == 1) {
            for (int c = 0; c < w; ++c) cv[c] = addmod(cv[c], vv[c], p);
          } else {
            std::vector<u64> scaled((size_t)nM * w, 0);
            for (int mo = 0; mo < nM; ++mo)
              dual_mul_acc(scaled.data() + (size_t)mo * w,
                           cv + (size_t)mo * w, vv.data(), w, p);
            carry[s] = scaled;
          }
        }
      }
    }
  }
  return true;
}

// Value of every output slot of a straight-line value tape (opcodes
// CONST/ADD/MUL/INV) at the integer leaves mod p. Returns false on a zero
// reciprocal.
[[maybe_unused]] bool seed_tape_value(const std::vector<int>& op, const std::vector<int>& a,
                     const std::vector<int>& b, const std::vector<u64>& cval,
                     const std::vector<u64>& leaf, int nLeaves, u64 p,
                     std::vector<u64>& val) {
  int n = (int)op.size();
  val.assign(nLeaves + n, 0);
  for (int L = 0; L < nLeaves; ++L) val[L] = leaf[L] % p;
  for (int i = 0; i < n; ++i) {
    int s = nLeaves + i;
    if (op[i] == OP_CONST) val[s] = cval[i] % p;
    else if (op[i] == OP_ADD) val[s] = addmod(val[a[i]], val[b[i]], p);
    else if (op[i] == OP_MUL) val[s] = mulmod(val[a[i]], val[b[i]], p);
    else { u64 av = val[a[i]] % p; if (av == 0) return false; val[s] = invmod(av, p); }
  }
  return true;
}

// Value and forward-mode duals (dual column j is d/d(leaf j), w = nLeaves + 1) of
// every output slot of a straight-line tape. Returns false on a zero reciprocal.
[[maybe_unused]] bool seed_tape_dual(const std::vector<int>& op, const std::vector<int>& a,
                    const std::vector<int>& b, const std::vector<u64>& cval,
                    const std::vector<u64>& leaf, int nLeaves, u64 p,
                    std::vector<std::vector<u64> >& val) {
  int n = (int)op.size(), w = nLeaves + 1;
  val.assign(nLeaves + n, std::vector<u64>(w, 0));
  for (int L = 0; L < nLeaves; ++L) { val[L][0] = leaf[L] % p; val[L][1 + L] = 1; }
  for (int i = 0; i < n; ++i) {
    int s = nLeaves + i;
    if (op[i] == OP_CONST) {
      val[s][0] = cval[i] % p;
    } else if (op[i] == OP_ADD) {
      const std::vector<u64>& A = val[a[i]]; const std::vector<u64>& B = val[b[i]];
      for (int c = 0; c < w; ++c) val[s][c] = addmod(A[c], B[c], p);
    } else if (op[i] == OP_MUL) {
      const std::vector<u64>& A = val[a[i]]; const std::vector<u64>& B = val[b[i]];
      val[s][0] = mulmod(A[0], B[0], p);
      for (int c = 1; c < w; ++c)
        val[s][c] = addmod(mulmod(A[0], B[c], p), mulmod(A[c], B[0], p), p);
    } else {
      const std::vector<u64>& A = val[a[i]];
      u64 a0 = A[0] % p; if (a0 == 0) return false;
      u64 vi = invmod(a0, p), vi2 = mulmod(vi, vi, p);
      val[s][0] = vi;
      for (int c = 1; c < w; ++c) val[s][c] = submod(0, mulmod(A[c], vi2, p), p);
    }
  }
  return true;
}

// Solve A X = B over GF(p), A n x n, B n x k. Returns false when A is singular.
[[maybe_unused]] bool seed_solve_mod(const std::vector<std::vector<u64> >& A,
                    const std::vector<std::vector<u64> >& B, u64 p,
                    std::vector<std::vector<u64> >& X) {
  int n = (int)A.size();
  int k = (n > 0 && !B.empty()) ? (int)B[0].size() : 0;
  std::vector<std::vector<u64> > M(n, std::vector<u64>(n + k, 0));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) M[i][j] = A[i][j] % p;
    for (int j = 0; j < k; ++j) M[i][n + j] = B[i][j] % p;
  }
  for (int col = 0; col < n; ++col) {
    int piv = -1;
    for (int r = col; r < n; ++r) if (M[r][col] % p) { piv = r; break; }
    if (piv < 0) return false;
    std::swap(M[col], M[piv]);
    u64 inv = invmod(M[col][col] % p, p);
    for (int j = 0; j < n + k; ++j) M[col][j] = mulmod(M[col][j], inv, p);
    for (int r = 0; r < n; ++r) if (r != col && M[r][col] % p) {
      u64 f = M[r][col] % p;
      for (int j = 0; j < n + k; ++j)
        M[r][j] = submod(M[r][j], mulmod(f, M[col][j], p), p);
    }
  }
  X.assign(n, std::vector<u64>(k, 0));
  for (int i = 0; i < n; ++i) for (int j = 0; j < k; ++j) X[i][j] = M[i][n + j];
  return true;
}

}  // namespace

// Multi-condition observability over a shared coordinate space. Each element of
// `tapes` is one experimental condition compiled against the same slot layout:
// leaves [0, nLeaves), states [nLeaves, nLeaves + nStates), then that
// condition's own instructions. A condition supplies op/a/b/cnum/cden, its
// integrated stateSlots (with aligned fOut), its gOut, and per-state initial
// conditions as either a leaf slot (icLeaf >= 0, carrying that leaf's value and
// dual) or a rational constant (icLeaf < 0, value icNum/icDen). The leaf point
// and the dual seeding are shared across conditions, so a parameter is one
// unknown everywhere; the observability rows of all conditions are stacked and
// reduced once, giving the intersection nullspace: a direction is
// non-identifiable only when every condition leaves it unobservable. ok is
// FALSE when any condition hits a vanishing denominator at the point. The
// per-condition Taylor build is independent across conditions, so it is run in
// parallel over `cores` OpenMP threads; the rows are merged and reduced once.
// [[Rcpp::export]]
List symObsNullMulti(List tapes, int nLeaves, int nStates,
                     IntegerVector zSlots, IntegerVector point,
                     double pIn, int Nt, int cores = 1) {
  u64 p = (u64)pIn;
  int nz = zSlots.size(), w = nz + 1;
  int T = tapes.size();

  std::vector<int> dualCol(nLeaves, -1);
  for (int c = 0; c < nz; ++c) dualCol[zSlots[c]] = c;

  std::vector<std::vector<u64> > leafVal(nLeaves, std::vector<u64>(w, 0));
  for (int L = 0; L < nLeaves; ++L) {
    i128 pv = (i128)point[L] % (i128)p;
    if (pv < 0) pv += p;
    leafVal[L][0] = (u64)pv;
    if (dualCol[L] >= 0) leafVal[L][1 + dualCol[L]] = 1;
  }

  // serial pre-pass: extract every condition's tape into thread-safe plain C++,
  // so the parallel region below touches no R object (Rcpp/R is not thread-safe).
  // A condition seeded directly from a numeric interior steady state carries its
  // modular seed and IFT parameter-duals in icSeed (constraint mode).
  std::vector<CondRaw> td(T);
  for (int t = 0; t < T; ++t) td[t] = extract_cond_raw(tapes[t], nStates, w);

  // parallel per-condition build: each thread owns its row block; a vanishing
  // denominator in any condition marks failure (no early return inside OpenMP).
  std::vector<std::vector<std::vector<u64> > > perRows(T);
  std::vector<char> okFlag(T, 1);
  #pragma omp parallel for num_threads(cores > 0 ? cores : 1) schedule(dynamic) \
          if (cores > 1 && T > 1)
  for (int t = 0; t < T; ++t)
    if (!build_one_condition(td[t], leafVal, nLeaves, nStates, nz, w, Nt, p, perRows[t]))
      okFlag[t] = 0;

  for (int t = 0; t < T; ++t)
    if (!okFlag[t]) return List::create(_["ok"] = false);

  std::vector<std::vector<u64> > rows;
  for (int t = 0; t < T; ++t)
    for (size_t i = 0; i < perRows[t].size(); ++i)
      rows.push_back(std::move(perRows[t][i]));

  std::vector<int> pivots = rref_mod(rows, p);
  int rank = (int)pivots.size();
  IntegerMatrix R(rank, nz);
  for (int i = 0; i < rank; ++i)
    for (int c = 0; c < nz; ++c) R(i, c) = (int)rows[i][c];
  IntegerVector piv(rank);
  for (int i = 0; i < rank; ++i) piv[i] = pivots[i];

  return List::create(_["ok"] = true, _["R"] = R, _["pivots"] = piv,
                      _["rank"] = rank, _["dim"] = nz);
}

// Batched single-segment observability: evaluate the shared condition tapes at
// many (point, prime) pairs in one call, one OpenMP task per pair, so the rational
// reconstruction's sample bank is built with full core occupancy. `points` is
// nB x nLeaves and `primes` has length nB, each row carrying its own prime;
// conditions within a pair are built serially and parallelism is over the batch.
// The tapes carry an IC tape or constant ICs, seeded per prime from the leaves.
// Returns a list of nB results, each shaped like symObsNullMulti (ok, R, pivots,
// rank, dim); a pair with a vanishing denominator yields ok=FALSE.
// [[Rcpp::export]]
List symObsNullBatch(List tapes, int nLeaves, int nStates, IntegerVector zSlots,
                     IntegerMatrix points, NumericVector primes, int Nt,
                     int cores = 1) {
  int nz = zSlots.size(), w = nz + 1;
  int T = tapes.size();
  int nB = points.nrow();

  std::vector<int> dualCol(nLeaves, -1);
  for (int c = 0; c < nz; ++c) dualCol[zSlots[c]] = c;

  // copy points and primes out of R into plain C++ so the parallel region below
  // reads no R object; every constant (tape rationals and icSeed) is reduced per
  // prime inside the build, so a single tape extraction serves every pair.
  std::vector<int> pts((size_t)nB * nLeaves);
  for (int bi = 0; bi < nB; ++bi)
    for (int L = 0; L < nLeaves; ++L) pts[(size_t)bi * nLeaves + L] = points(bi, L);
  std::vector<u64> pr(nB);
  for (int bi = 0; bi < nB; ++bi) pr[bi] = (u64)primes[bi];

  std::vector<CondRaw> td(T);
  for (int t = 0; t < T; ++t) td[t] = extract_cond_raw(tapes[t], nStates, w);

  std::vector<char> okFlag(nB, 1);
  std::vector<std::vector<std::vector<u64> > > redRows(nB);
  std::vector<std::vector<int> > redPiv(nB);
  #pragma omp parallel for num_threads(cores > 0 ? cores : 1) schedule(dynamic) \
          if (cores > 1 && nB > 1)
  for (int bi = 0; bi < nB; ++bi) {
    u64 p = pr[bi];
    std::vector<std::vector<u64> > leafVal(nLeaves, std::vector<u64>(w, 0));
    for (int L = 0; L < nLeaves; ++L) {
      i128 pv = (i128)pts[(size_t)bi * nLeaves + L] % (i128)p; if (pv < 0) pv += p;
      leafVal[L][0] = (u64)pv;
      if (dualCol[L] >= 0) leafVal[L][1 + dualCol[L]] = 1;
    }
    std::vector<std::vector<u64> > rows;
    bool ok = true;
    for (int t = 0; t < T && ok; ++t)
      if (!build_one_condition(td[t], leafVal, nLeaves, nStates, nz, w, Nt, p, rows))
        ok = false;
    if (!ok) { okFlag[bi] = 0; continue; }
    std::vector<int> pivots = rref_mod(rows, p);
    int rank = (int)pivots.size();
    redRows[bi].assign(rank, std::vector<u64>(nz));
    for (int i = 0; i < rank; ++i)
      for (int c = 0; c < nz; ++c) redRows[bi][i][c] = rows[i][c];
    redPiv[bi] = pivots;
  }

  List out(nB);
  for (int bi = 0; bi < nB; ++bi) {
    if (!okFlag[bi]) { out[bi] = List::create(_["ok"] = false); continue; }
    int rank = (int)redPiv[bi].size();
    IntegerMatrix R(rank, nz);
    for (int i = 0; i < rank; ++i)
      for (int c = 0; c < nz; ++c) R(i, c) = (int)redRows[bi][i][c];
    IntegerVector piv(rank);
    for (int i = 0; i < rank; ++i) piv[i] = redPiv[bi][i];
    out[bi] = List::create(_["ok"] = true, _["R"] = R, _["pivots"] = piv,
                           _["rank"] = rank, _["dim"] = nz);
  }
  return out;
}

// Multi-segment observability with exact generic-timing propagation across events.
// `chains` is one entry per condition: an ordered list of segments. Each segment
// supplies its regime-substituted tape (op/a/b/cnum/cden, stateSlots, fOut, gOut);
// the first additionally supplies an icSeed (nStates x w) anchoring the state at
// the first event time. A condition with nSeg segments has K = nSeg-1 inter-event
// gaps, each a formal length Delta t_i; the state is propagated across a gap by
// promoting the segment's local-time Taylor coefficients into that gap's axis (an
// identity event map: the regime change leaves the state continuous). Every
// segment's output jet, expanded to total gap-degree Mtot, contributes its
// z-gradient rows for every gap monomial; a direction is non-identifiable iff it is
// annihilated at every monomial (a polynomial vanishing for generic gaps), so the
// stacked rows are reduced once over GF(p) exactly as in the single-segment path.
// [[Rcpp::export]]
List symObsNullChain(List chains, int nLeaves, int nStates, IntegerVector zSlots,
                     IntegerVector point, double pIn, int Nt, int Mtot, int cores = 1) {
  u64 p = (u64)pIn;
  int nz = zSlots.size(), w = nz + 1;
  std::vector<int> dualCol(nLeaves, -1);
  for (int c = 0; c < nz; ++c) dualCol[zSlots[c]] = c;
  int T = chains.size();

  // serial pre-pass: copy the leaf point and extract every condition's segments
  // into thread-safe plain C++, so the parallel build below touches no R object.
  std::vector<int> leafPt(nLeaves);
  for (int L = 0; L < nLeaves; ++L) leafPt[L] = point[L];
  std::vector<std::vector<SegRaw> > chainsRaw(T);
  for (int t = 0; t < T; ++t) {
    List segs = chains[t];
    int nSeg = segs.size();
    chainsRaw[t].reserve(nSeg);
    for (int sj = 0; sj < nSeg; ++sj)
      chainsRaw[t].push_back(extract_seg_raw(segs[sj], nStates, w, p));
  }

  // parallel per-condition build: each thread owns its row block; a vanishing
  // denominator in any condition marks failure (no early return inside OpenMP).
  std::vector<std::vector<std::vector<u64> > > perRows(T);
  std::vector<char> okFlag(T, 1);
  #pragma omp parallel for num_threads(cores > 0 ? cores : 1) schedule(dynamic) \
          if (cores > 1 && T > 1)
  for (int t = 0; t < T; ++t)
    if (!build_one_chain(chainsRaw[t], nLeaves, nStates, nz, w, dualCol, leafPt,
                         Nt, Mtot, p, perRows[t]))
      okFlag[t] = 0;

  for (int t = 0; t < T; ++t)
    if (!okFlag[t]) return List::create(_["ok"] = false);

  std::vector<std::vector<u64> > rows;
  for (int t = 0; t < T; ++t)
    for (size_t i = 0; i < perRows[t].size(); ++i)
      rows.push_back(std::move(perRows[t][i]));

  std::vector<int> pivots = rref_mod(rows, p);
  int rank = (int)pivots.size();
  IntegerMatrix R(rank, nz);
  for (int i = 0; i < rank; ++i)
    for (int c = 0; c < nz; ++c) R(i, c) = (int)rows[i][c];
  IntegerVector piv(rank);
  for (int i = 0; i < rank; ++i) piv[i] = pivots[i];
  return List::create(_["ok"] = true, _["R"] = R, _["pivots"] = piv,
                      _["rank"] = rank, _["dim"] = nz);
}

// Batched seeded single-chain observability: evaluate, for many (chain, seed, prime)
// triples, one condition's multi-segment observability rows, one OpenMP task per
// triple. This is the batched twin of the joint/equilibrate path's per-condition,
// per-point symObsNullChain calls (which run serially, cores=1, one condition each):
// the reconstruction sample bank supplies every (sample point, condition) pair with
// its pre-solved steady-state seed, and they are evaluated in parallel here. `chains`
// is a list of conditions' segment lists; evalChain[e] selects the chain for eval e;
// `seeds` is nB x nLeaves (the seed leaf residues, states already written in); each
// eval carries its own prime. Segments are extracted once per (chain, distinct prime)
// in a serial pre-pass (extract_seg_raw reduces constants per prime and touches R).
// Returns a list of nB results, each shaped like symObsNullChain (ok, R = rank x nz
// reduced rows, pivots, rank, dim); a vanishing denominator yields ok = FALSE.
// [[Rcpp::export]]
List symObsNullChainSeedBatch(List chains, IntegerVector evalChain,
                              IntegerMatrix seeds, NumericVector primes,
                              int nLeaves, int nStates, IntegerVector zSlots,
                              int Nt, int Mtot, int cores = 1) {
  int nz = zSlots.size(), w = nz + 1;
  int T = chains.size();
  int nB = evalChain.size();
  std::vector<int> dualCol(nLeaves, -1);
  for (int c = 0; c < nz; ++c) dualCol[zSlots[c]] = c;

  // primes and their distinct set (extraction is per distinct prime)
  std::vector<u64> pr(nB);
  for (int e = 0; e < nB; ++e) pr[e] = (u64)primes[e];
  std::vector<u64> distinct;
  std::vector<int> primeIdx(nB);
  for (int e = 0; e < nB; ++e) {
    int idx = -1;
    for (size_t k = 0; k < distinct.size(); ++k) if (distinct[k] == pr[e]) { idx = (int)k; break; }
    if (idx < 0) { idx = (int)distinct.size(); distinct.push_back(pr[e]); }
    primeIdx[e] = idx;
  }
  int nPr = (int)distinct.size();

  // serial pre-pass: extract every chain's segments per distinct prime into
  // thread-safe plain C++, so the parallel build below touches no R object.
  std::vector<std::vector<std::vector<SegRaw> > >
      segsRaw(nPr, std::vector<std::vector<SegRaw> >(T));
  for (int t = 0; t < T; ++t) {
    List segs = chains[t];
    int nSeg = segs.size();
    for (int pi = 0; pi < nPr; ++pi) {
      segsRaw[pi][t].reserve(nSeg);
      for (int sj = 0; sj < nSeg; ++sj)
        segsRaw[pi][t].push_back(extract_seg_raw(segs[sj], nStates, w, distinct[pi]));
    }
  }

  std::vector<int> ec(nB);
  for (int e = 0; e < nB; ++e) ec[e] = evalChain[e];
  std::vector<int> sd((size_t)nB * nLeaves);
  for (int e = 0; e < nB; ++e)
    for (int L = 0; L < nLeaves; ++L) sd[(size_t)e * nLeaves + L] = seeds(e, L);

  std::vector<char> okFlag(nB, 1);
  std::vector<std::vector<std::vector<u64> > > redRows(nB);
  std::vector<std::vector<int> > redPiv(nB);
  #pragma omp parallel for num_threads(cores > 0 ? cores : 1) schedule(dynamic) \
          if (cores > 1 && nB > 1)
  for (int e = 0; e < nB; ++e) {
    u64 p = pr[e];
    std::vector<int> leafPt(nLeaves);
    for (int L = 0; L < nLeaves; ++L) leafPt[L] = sd[(size_t)e * nLeaves + L];
    std::vector<std::vector<u64> > rows;
    if (!build_one_chain(segsRaw[primeIdx[e]][ec[e]], nLeaves, nStates, nz, w,
                         dualCol, leafPt, Nt, Mtot, p, rows)) { okFlag[e] = 0; continue; }
    std::vector<int> pivots = rref_mod(rows, p);
    int rank = (int)pivots.size();
    redRows[e].assign(rank, std::vector<u64>(nz));
    for (int i = 0; i < rank; ++i)
      for (int c = 0; c < nz; ++c) redRows[e][i][c] = rows[i][c];
    redPiv[e] = pivots;
  }

  List out(nB);
  for (int e = 0; e < nB; ++e) {
    if (!okFlag[e]) { out[e] = List::create(_["ok"] = false); continue; }
    int rank = (int)redPiv[e].size();
    IntegerMatrix R(rank, nz);
    for (int i = 0; i < rank; ++i)
      for (int c = 0; c < nz; ++c) R(i, c) = (int)redRows[e][i][c];
    IntegerVector piv(rank);
    for (int i = 0; i < rank; ++i) piv[i] = redPiv[e][i];
    out[e] = List::create(_["ok"] = true, _["R"] = R, _["pivots"] = piv,
                          _["rank"] = rank, _["dim"] = nz);
  }
  return out;
}

// Solve A x = b over GF(p); returns the solution (free variables zero) or
// R_NilValue when the system is inconsistent.
// [[Rcpp::export]]
SEXP symSolveMod(IntegerMatrix A, IntegerVector b, double pIn) {
  u64 p = (u64)pIn;
  int nr = A.nrow(), nc = A.ncol();
  std::vector<std::vector<u64> > aug(nr, std::vector<u64>(nc + 1));
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) aug[i][j] = red(A(i, j), p);
    aug[i][nc] = red(b[i], p);
  }
  std::vector<int> pivots = rref_mod(aug, p);
  for (size_t i = 0; i < pivots.size(); ++i)
    if (pivots[i] == nc) return R_NilValue;
  IntegerVector x(nc);
  for (size_t i = 0; i < pivots.size(); ++i)
    if (pivots[i] < nc) x[pivots[i]] = (int)aug[i][nc];
  return x;
}

// Reduced row echelon form of a stacked GF(p) matrix. Entries arrive as doubles
// already in [0, p) (they are log-normalised residues assembled in R); returns the
// pivot (rank) rows in pivot-column order, the 0-based pivot columns, and the rank.
// This is the compiled twin of R's .sym_rref_modp for the joint/equilibrate path,
// where the final reduction runs per accepted sample point.
// [[Rcpp::export]]
List symRrefMod(NumericMatrix M, double pIn) {
  u64 p = (u64)pIn;
  int nr = M.nrow(), nc = M.ncol();
  std::vector<std::vector<u64> > A(nr, std::vector<u64>(nc));
  for (int i = 0; i < nr; ++i)
    for (int j = 0; j < nc; ++j) {
      i128 v = (i128)(long long)M(i, j) % (i128)p;
      if (v < 0) v += p;
      A[i][j] = (u64)v;
    }
  std::vector<int> pivots = rref_mod(A, p);
  int rank = (int)pivots.size();
  NumericMatrix R(rank, nc);
  for (int i = 0; i < rank; ++i)
    for (int j = 0; j < nc; ++j) R(i, j) = (double)A[i][j];
  IntegerVector piv(rank);
  for (int i = 0; i < rank; ++i) piv[i] = pivots[i];   // 0-based, ascending
  return List::create(_["R"] = R, _["piv"] = piv, _["rank"] = rank);
}

// Fit a rational function num/den to samples of one nullspace entry over GF(p).
// sampleU holds the values of the relevant variables at each sample; mons holds
// the monomial exponents (row 0 is the constant monomial); rvals holds the entry
// value at each sample. Monomials are formed modulo p here to avoid overflow.
// The relation num - r*den = 0 is homogeneous in the 2*nMon coefficients, so its
// solution is the kernel of the sample matrix; num and den share one scale that
// cancels (this also fits denominators with no constant term). status is
// "inconsistent" when only the trivial solution exists (raise the degree),
// "ambiguous" when the kernel is more than one dimensional (degree too high),
// and "ok" with the canonical kernel vector (free coefficient set to one) and
// the free column index, which must agree across primes.
// [[Rcpp::export]]
List symFitRational(IntegerMatrix sampleU, IntegerMatrix mons,
                    IntegerVector rvals, double pIn) {
  u64 p = (u64)pIn;
  int nS = sampleU.nrow(), nrel = sampleU.ncol(), nMon = mons.nrow();
  int ncols = 2 * nMon;
  std::vector<std::vector<u64> > A(nS, std::vector<u64>(ncols, 0));
  for (int s = 0; s < nS; ++s) {
    i128 rv = (i128)rvals[s] % (i128)p; if (rv < 0) rv += p;
    u64 r = (u64)rv;
    std::vector<u64> uu(nrel);
    for (int k = 0; k < nrel; ++k) {
      i128 uv = (i128)sampleU(s, k) % (i128)p; if (uv < 0) uv += p;
      uu[k] = (u64)uv;
    }
    for (int j = 0; j < nMon; ++j) {
      u64 mon = 1;
      for (int k = 0; k < nrel; ++k) mon = mulmod(mon, powmod(uu[k], (u64)mons(j, k), p), p);
      A[s][j] = mon;
      A[s][nMon + j] = submod(0, mulmod(r, mon, p), p);
    }
  }
  std::vector<int> pivots = rref_mod(A, p);
  int rank = (int)pivots.size();
  int nullity = ncols - rank;
  if (nullity == 0) return List::create(_["status"] = "inconsistent");
  if (nullity > 1) return List::create(_["status"] = "ambiguous");
  std::vector<char> isPiv(ncols, 0);
  for (int c : pivots) isPiv[c] = 1;
  int freeCol = 0;
  for (int c = 0; c < ncols; ++c) if (!isPiv[c]) { freeCol = c; break; }
  IntegerVector coeffs(ncols);
  coeffs[freeCol] = 1;
  for (int ri = 0; ri < rank; ++ri)
    coeffs[pivots[ri]] = (int)submod(0, A[ri][freeCol], p);
  return List::create(_["status"] = "ok", _["coeffs"] = coeffs,
                      _["freeCol"] = freeCol);
}

// Per-row CRT over the given primes followed by rational reconstruction.
// Returns num/den as decimal strings; den is "0" when reconstruction fails.
// [[Rcpp::export]]
List symRatRecon(IntegerMatrix residues, IntegerVector primes) {
  int k = residues.nrow(), nprime = residues.ncol();
  // the CRT accumulator holds the running prime product in u128; r + M*t stays
  // below 2*product, so the product must fit under 2^127 to avoid silent
  // overflow. Primes near 2^31 allow up to 4; guard exactly by summed bit length.
  int prodBits = 0;
  for (int j = 0; j < nprime; ++j) {
    u64 pj = (u64)primes[j];
    prodBits += pj ? (64 - __builtin_clzll(pj)) : 0;
  }
  if (prodBits >= 127)
    stop("symRatRecon: prime product exceeds the u128 CRT capacity (2^127); "
         "reduce the number or size of primes in .symPrimes");
  CharacterVector num(k), den(k);
  for (int row = 0; row < k; ++row) {
    u128 r = 0, M = 1;
    bool first = true;
    for (int j = 0; j < nprime; ++j) {
      u64 pj = (u64)primes[j];
      i128 rv = (i128)residues(row, j) % (i128)pj; if (rv < 0) rv += pj;
      u64 res = (u64)rv;
      if (first) { r = res; M = pj; first = false; continue; }
      u64 Mmod = (u64)(M % pj);
      u64 diff = submod(res, (u64)(r % pj), pj);
      u64 t = mulmod(diff, invmod(Mmod, pj), pj);
      r = r + M * (u128)t;
      M = M * (u128)pj;
    }
    i128 n, d;
    if (rational_reconstruct(r, M, n, d)) {
      num[row] = i128_to_string(n);
      den[row] = i128_to_string(d);
    } else {
      num[row] = "0";
      den[row] = "0";
    }
  }
  return List::create(_["num"] = num, _["den"] = den);
}


// Berlekamp-Massey: minimal connection polynomial C (C[0]=1) of a GF(p) sequence.
// The sequence satisfies s[n] = -sum_{i>=1} C[i] s[n-i]; deg C is the LFSR length.
static std::vector<u64> berlekamp_massey(const std::vector<u64>& s, u64 p) {
  std::vector<u64> C(1, 1), B(1, 1);
  int L = 0, m = 1;
  u64 b = 1;
  for (int n = 0; n < (int)s.size(); ++n) {
    u64 d = s[n] % p;
    for (int i = 1; i <= L; ++i) d = addmod(d, mulmod(C[i], s[n - i], p), p);
    if (d == 0) { ++m; continue; }
    u64 coef = mulmod(d, invmod(b, p), p);
    if ((int)C.size() < (int)B.size() + m) C.resize(B.size() + m, 0);
    if (2 * L <= n) {
      std::vector<u64> T = C;
      for (int i = 0; i < (int)B.size(); ++i)
        C[i + m] = submod(C[i + m], mulmod(coef, B[i], p), p);
      L = n + 1 - L; B = T; b = d; m = 1;
    } else {
      for (int i = 0; i < (int)B.size(); ++i)
        C[i + m] = submod(C[i + m], mulmod(coef, B[i], p), p);
      ++m;
    }
  }
  C.resize(L + 1, 0);
  return C;
}

// Sparse multivariate polynomial interpolation (Ben-Or-Tiwari) from a geometric
// evaluation sequence seq[k] = f(b_1^k, ..., b_n^k) mod p. monoTab lists candidate
// exponent vectors (rows) and monoRes their monomial values prod b_j^{e_j} mod p;
// the recovered term monomials are the candidates whose value is a root of the
// Berlekamp-Massey locator, and their coefficients come from a transposed
// Vandermonde solve. status is "ok", "needmore" (sequence too short for the
// recovered order) or "noroots" (locator roots are not legal monomial values).
// [[Rcpp::export]]
List symSparsePoly(IntegerVector seq, IntegerMatrix monoTab, IntegerVector monoRes,
                   double pIn) {
  u64 p = (u64) pIn;
  int len = seq.size();
  int nvar = monoTab.ncol();
  std::vector<u64> s(len);
  for (int i = 0; i < len; ++i) s[i] = red(seq[i], p);

  std::vector<u64> C = berlekamp_massey(s, p);
  int L = (int)C.size() - 1;
  if (L == 0)
    return List::create(_["status"] = "ok", _["nterms"] = 0,
                        _["exps"] = IntegerMatrix(0, nvar),
                        _["coeffs"] = IntegerVector(0));
  if (2 * L > len) return List::create(_["status"] = "needmore");

  // roots of the characteristic polynomial chi(x) = sum_i C[i] x^{L-i}
  int ncand = monoRes.size();
  std::vector<int> rootIdx;
  for (int c = 0; c < ncand; ++c) {
    u64 m = red(monoRes[c], p);
    u64 r = 0;
    for (int i = 0; i <= L; ++i) r = addmod(mulmod(r, m, p), C[i] % p, p);
    if (r == 0) rootIdx.push_back(c);
  }
  if ((int)rootIdx.size() != L) return List::create(_["status"] = "noroots");

  // transposed Vandermonde solve: sum_j coeff_j * nodes_j^k = s[k], k = 0..L-1
  int t = L;
  std::vector<u64> nodes(t);
  for (int j = 0; j < t; ++j)
    nodes[j] = red(monoRes[rootIdx[j]], p);
  std::vector<std::vector<u64> > A(t, std::vector<u64>(t + 1, 0));
  for (int k = 0; k < t; ++k) {
    for (int j = 0; j < t; ++j) A[k][j] = powmod(nodes[j], (u64)k, p);
    A[k][t] = s[k];
  }
  std::vector<int> piv = rref_mod(A, p);
  if ((int)piv.size() != t) return List::create(_["status"] = "noroots");

  IntegerMatrix exps(t, nvar);
  IntegerVector coeffs(t);
  for (int j = 0; j < t; ++j) {
    for (int v = 0; v < nvar; ++v) exps(j, v) = monoTab(rootIdx[j], v);
    coeffs[j] = (int)A[j][t];
  }
  return List::create(_["status"] = "ok", _["nterms"] = t,
                      _["exps"] = exps, _["coeffs"] = coeffs);
}


// Modular values of the monomials prod_j bases_j^{expts(i,j)} mod p, one per row
// of `expts`. Negative exponents use the modular inverse, so Laurent monomials are
// supported. Builds the candidate-monomial residues for sparse interpolation.
// [[Rcpp::export]]
IntegerVector symMonoResidues(IntegerMatrix expts, IntegerVector bases, double pIn) {
  u64 p = (u64) pIn;
  int nr = expts.nrow(), nc = expts.ncol();
  IntegerVector out(nr);
  for (int i = 0; i < nr; ++i) {
    u64 v = 1;
    for (int j = 0; j < nc; ++j) {
      int e = expts(i, j);
      u64 b = red(bases[j], p);
      if (e >= 0) v = mulmod(v, powmod(b, (u64)e, p), p);
      else v = mulmod(v, invmod(powmod(b, (u64)(-e), p), p), p);
    }
    out[i] = (int)v;
  }
  return out;
}


// Berlekamp-Massey order (LFSR length) of a GF(p) sequence: the number of terms
// the sequence requires, used to grow the sparse sampling until it stabilises.
// [[Rcpp::export]]
int symBMorder(IntegerVector seq, double pIn) {
  u64 p = (u64) pIn;
  int len = seq.size();
  std::vector<u64> s(len);
  for (int i = 0; i < len; ++i) s[i] = red(seq[i], p);
  return (int)berlekamp_massey(s, p).size() - 1;
}


// Univariate Cauchy (rational) interpolation A(t)/B(t) = r(t) over GF(p), with
// deg A = dN, deg B = dD, from samples (tnodes, rvals); normalised so B(0) = 1.
// Returns the values A(1) and B(1) (= N(point)/D(s) and D(point)/D(s) when the
// caller samples r along the ray s + t*(point - s)). status is "ok",
// "ambiguous" (the (dN, dD) guess does not give a one-dimensional fit) or
// "badshift" (B(0) = 0 at this shift). Used by the general sparse-rational path.
// [[Rcpp::export]]
List symCauchyEval(IntegerVector tnodes, IntegerVector rvals, int dN, int dD,
                   double pIn) {
  u64 p = (u64) pIn;
  int ns = tnodes.size(), nc = dN + dD + 2;
  std::vector<std::vector<u64> > M(ns, std::vector<u64>(nc, 0));
  for (int s = 0; s < ns; ++s) {
    u64 t = red(tnodes[s], p);
    u64 r = red(rvals[s], p);
    u64 tp = 1;
    for (int d = 0; d <= dN; ++d) { M[s][d] = tp; tp = mulmod(tp, t, p); }
    tp = 1;
    for (int e = 0; e <= dD; ++e) {
      M[s][dN + 1 + e] = (p - mulmod(r, tp, p)) % p;
      tp = mulmod(tp, t, p);
    }
  }
  std::vector<int> piv = rref_mod(M, p);
  std::vector<bool> isPiv(nc, false);
  for (size_t i = 0; i < piv.size(); ++i) isPiv[piv[i]] = true;
  int fcol = -1;
  for (int c = 0; c < nc; ++c) if (!isPiv[c]) { if (fcol >= 0) { fcol = -2; break; } fcol = c; }
  if (fcol < 0) return List::create(_["status"] = "ambiguous");
  std::vector<u64> v(nc, 0);
  v[fcol] = 1;
  for (size_t ri = 0; ri < piv.size(); ++ri) v[piv[ri]] = (p - M[ri][fcol]) % p;
  u64 b0 = v[dN + 1] % p;
  if (b0 == 0) return List::create(_["status"] = "badshift");
  u64 ib0 = invmod(b0, p), A1 = 0, B1 = 0;
  for (int d = 0; d <= dN; ++d) A1 = addmod(A1, v[d], p);
  for (int e = 0; e <= dD; ++e) B1 = addmod(B1, v[dN + 1 + e], p);
  return List::create(_["status"] = "ok",
                      _["N"] = (double) mulmod(A1, ib0, p),
                      _["D"] = (double) mulmod(B1, ib0, p));
}
