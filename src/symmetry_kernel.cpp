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

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

typedef uint64_t u64;
typedef __int128 i128;
typedef unsigned __int128 u128;

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

inline std::vector<int> to_ivec(const IntegerVector& v) {
  return std::vector<int>(v.begin(), v.end());
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

// Run the Taylor recurrence on a value table whose leaves and state order-0
// coefficients are already seeded, then collect the observability rows (the
// dual parts of every time coefficient of every observable). Instruction i
// writes slot instrBase + i. Returns false when a reciprocal hits a zero
// leading coefficient at the evaluation point.
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
  for (int gi = 0; gi < gOut.size(); ++gi) {
    const u64* g = val[gOut[gi]].data();
    for (int k = 0; k <= Nt; ++k) {
      std::vector<u64> row(nz);
      for (int c = 0; c < nz; ++c) row[c] = g[(size_t)k * w + 1 + c];
      rows.push_back(row);
    }
  }
  return true;
}

}  // namespace

// Build the observability matrix from a tape at one integer point modulo p and
// return its reduced row echelon form (pivot rows only) together with the pivot
// columns. ok is FALSE when a denominator vanishes at the point.
// [[Rcpp::export]]
List symObsNull(IntegerVector op, IntegerVector a, IntegerVector b,
                CharacterVector cnum, CharacterVector cden,
                int nLeaves, IntegerVector stateSlots, IntegerVector fOut,
                IntegerVector zSlots, IntegerVector gOut,
                IntegerVector point, double pIn, int Nt) {
  u64 p = (u64)pIn;
  // stateSlots: all integrated states (aligned with fOut). zSlots: the
  // dual-carrying coordinates z (free states then free parameters); a fixed
  // state is integrated but absent from zSlots, so it carries no dual.
  int nz = zSlots.size(), w = nz + 1;
  int nInstr = op.size();
  int S = nLeaves + nInstr;

  std::vector<int> dualCol(nLeaves, -1);
  for (int c = 0; c < nz; ++c) dualCol[zSlots[c]] = c;

  std::vector<u64> cval(nInstr, 0);
  for (int i = 0; i < nInstr; ++i)
    if (op[i] == OP_CONST)
      cval[i] = mulmod(parse_mod(as<std::string>(cnum[i]), p),
                       invmod(parse_mod(as<std::string>(cden[i]), p), p), p);

  std::vector<std::vector<u64> > val(S, std::vector<u64>((size_t)(Nt + 1) * w, 0));
  for (int L = 0; L < nLeaves; ++L) {
    i128 pv = (i128)point[L] % (i128)p;
    if (pv < 0) pv += p;
    val[L][0] = (u64)pv;
    if (dualCol[L] >= 0) val[L][1 + dualCol[L]] = 1;
  }

  std::vector<std::vector<u64> > rows;
  if (!build_obs_rows(val, to_ivec(op), to_ivec(a), to_ivec(b), cval, nLeaves,
                      to_ivec(stateSlots), to_ivec(fOut), to_ivec(gOut),
                      nz, w, Nt, p, rows))
    return List::create(_["ok"] = false);

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
  int instrBase = nLeaves + nStates;
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

  // serial pre-pass: extract every condition's tape into plain-C++ data and
  // precompute the modular constants, so the parallel region below touches no R
  // objects (Rcpp/R is not thread-safe).
  struct CondTape {
    std::vector<int> op, a, b, stateSlots, fOut, gOut, icLeaf, icOp, icA, icB, icOut;
    std::vector<u64> cval, icCval, icConst;
    bool hasIcTape, hasIcSeed;
    std::vector<std::vector<u64> > icSeed;
  };
  std::vector<CondTape> td(T);
  for (int t = 0; t < T; ++t) {
    List tp = tapes[t];
    CondTape& cd = td[t];
    cd.op = to_ivec(tp["op"]); cd.a = to_ivec(tp["a"]); cd.b = to_ivec(tp["b"]);
    cd.stateSlots = to_ivec(tp["stateSlots"]);
    cd.fOut = to_ivec(tp["fOut"]); cd.gOut = to_ivec(tp["gOut"]);
    cd.icLeaf = to_ivec(tp["icLeaf"]);
    cd.hasIcTape = tp.containsElementNamed("icOp");
    // constraint mode: states seeded directly from a numeric interior
    // steady-state point and its IFT sensitivities (icSeed: nStates x w),
    // recomputed by R per evaluation point and prime; overrides the IC tape.
    cd.hasIcSeed = tp.containsElementNamed("icSeed");
    CharacterVector cnum = tp["cnum"], cden = tp["cden"];
    int nInstr = cd.op.size();
    cd.cval.assign(nInstr, 0);
    for (int i = 0; i < nInstr; ++i)
      if (cd.op[i] == OP_CONST)
        cd.cval[i] = mulmod(parse_mod(as<std::string>(cnum[i]), p),
                            invmod(parse_mod(as<std::string>(cden[i]), p), p), p);
    if (cd.hasIcSeed) {
      IntegerMatrix icSeed = as<IntegerMatrix>(tp["icSeed"]);
      cd.icSeed.assign(nStates, std::vector<u64>(w, 0));
      for (int i = 0; i < nStates; ++i)
        for (int c = 0; c < w; ++c) {
          i128 v = (i128)icSeed(i, c) % (i128)p; if (v < 0) v += p;
          cd.icSeed[i][c] = (u64)v;
        }
    } else if (cd.hasIcTape) {
      cd.icOp = to_ivec(tp["icOp"]); cd.icA = to_ivec(tp["icA"]);
      cd.icB = to_ivec(tp["icB"]); cd.icOut = to_ivec(tp["icOut"]);
      CharacterVector icCnum = tp["icCnum"], icCden = tp["icCden"];
      int nIc = cd.icOp.size();
      cd.icCval.assign(nIc, 0);
      for (int i = 0; i < nIc; ++i)
        if (cd.icOp[i] == OP_CONST)
          cd.icCval[i] = mulmod(parse_mod(as<std::string>(icCnum[i]), p),
                                invmod(parse_mod(as<std::string>(icCden[i]), p), p), p);
    } else {
      CharacterVector icNum = tp["icNum"], icDen = tp["icDen"];
      cd.icConst.assign(nStates, 0);
      for (int i = 0; i < nStates; ++i)
        cd.icConst[i] = mulmod(parse_mod(as<std::string>(icNum[i]), p),
                               invmod(parse_mod(as<std::string>(icDen[i]), p), p), p);
    }
  }

  // parallel per-condition build: each thread owns its row block; a vanishing
  // denominator in any condition marks failure (no early return inside OpenMP).
  std::vector<std::vector<std::vector<u64> > > perRows(T);
  std::vector<char> okFlag(T, 1);
  #pragma omp parallel for num_threads(cores > 0 ? cores : 1) schedule(dynamic) \
          if (cores > 1 && T > 1)
  for (int t = 0; t < T; ++t) {
    const CondTape& cd = td[t];
    int nInstr = cd.op.size();
    int S = instrBase + nInstr;

    std::vector<std::vector<u64> > icVal;
    if (cd.hasIcTape && !cd.hasIcSeed) {
      int nIc = cd.icOp.size();
      icVal.assign(nLeaves + nIc, std::vector<u64>(w, 0));
      for (int L = 0; L < nLeaves; ++L)
        for (int c = 0; c < w; ++c) icVal[L][c] = leafVal[L][c];
      if (!eval_ic_order0(icVal, cd.icOp, cd.icA, cd.icB, cd.icCval, nLeaves, w, p)) {
        okFlag[t] = 0; continue;
      }
    }

    std::vector<std::vector<u64> > val(S, std::vector<u64>((size_t)(Nt + 1) * w, 0));
    for (int L = 0; L < nLeaves; ++L)
      for (int c = 0; c < w; ++c) val[L][c] = leafVal[L][c];
    for (int i = 0; i < (int)cd.stateSlots.size(); ++i) {
      int slot = cd.stateSlots[i];
      if (cd.hasIcSeed) {
        for (int c = 0; c < w; ++c) val[slot][c] = cd.icSeed[i][c];
        continue;
      }
      int src = cd.hasIcTape ? cd.icOut[i] : (cd.icLeaf[i] >= 0 ? cd.icLeaf[i] : -1);
      if (src >= 0) {
        const std::vector<u64>& s = src < nLeaves ? leafVal[src] : icVal[src];
        for (int c = 0; c < w; ++c) val[slot][c] = s[c];
      } else {
        val[slot][0] = cd.icConst[i];
      }
    }

    if (!build_obs_rows(val, cd.op, cd.a, cd.b, cd.cval, instrBase, cd.stateSlots,
                        cd.fOut, cd.gOut, nz, w, Nt, p, perRows[t]))
      okFlag[t] = 0;
  }

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

// Solve A x = b over GF(p); returns the solution (free variables zero) or
// R_NilValue when the system is inconsistent.
// [[Rcpp::export]]
SEXP symSolveMod(IntegerMatrix A, IntegerVector b, double pIn) {
  u64 p = (u64)pIn;
  int nr = A.nrow(), nc = A.ncol();
  std::vector<std::vector<u64> > aug(nr, std::vector<u64>(nc + 1));
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {
      i128 v = (i128)A(i, j) % (i128)p; if (v < 0) v += p;
      aug[i][j] = (u64)v;
    }
    i128 v = (i128)b[i] % (i128)p; if (v < 0) v += p;
    aug[i][nc] = (u64)v;
  }
  std::vector<int> pivots = rref_mod(aug, p);
  for (size_t i = 0; i < pivots.size(); ++i)
    if (pivots[i] == nc) return R_NilValue;
  IntegerVector x(nc);
  for (size_t i = 0; i < pivots.size(); ++i)
    if (pivots[i] < nc) x[pivots[i]] = (int)aug[i][nc];
  return x;
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
  for (int i = 0; i < len; ++i) s[i] = (u64)(((i128)seq[i] % (i128)p + p) % p);

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
    u64 m = (u64)(((i128)monoRes[c] % (i128)p + p) % p);
    u64 r = 0;
    for (int i = 0; i <= L; ++i) r = addmod(mulmod(r, m, p), C[i] % p, p);
    if (r == 0) rootIdx.push_back(c);
  }
  if ((int)rootIdx.size() != L) return List::create(_["status"] = "noroots");

  // transposed Vandermonde solve: sum_j coeff_j * nodes_j^k = s[k], k = 0..L-1
  int t = L;
  std::vector<u64> nodes(t);
  for (int j = 0; j < t; ++j)
    nodes[j] = (u64)(((i128)monoRes[rootIdx[j]] % (i128)p + p) % p);
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
      u64 b = (u64)(((i128)bases[j] % (i128)p + p) % p);
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
  for (int i = 0; i < len; ++i) s[i] = (u64)(((i128)seq[i] % (i128)p + p) % p);
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
    u64 t = (u64)(((i128)tnodes[s] % (i128)p + p) % p);
    u64 r = (u64)(((i128)rvals[s] % (i128)p + p) % p);
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
