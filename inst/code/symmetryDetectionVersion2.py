# Structural identifiability of ODE models. Self-contained module driven from
# R/symmetryDetection.R via reticulate. Three exact engines, selected by the
# R-side `method` argument:
#
#   "observability": local identifiability from the rank of the
#     observability-identifiability matrix O = d/dz [g, L_f g, ...]. For a
#     rational model O is built by a Taylor-mode construction over GF(p) at a
#     rational point (no symbolic O), the rank certified across several primes;
#     this scales to large, deep systems. Non-rational observables route to a
#     symbolic build whose rank is certified at several generic points, with an
#     optional closed-form nullspace reconstructed only in the relevant
#     variables.
#   "liesym": the polynomial Lie-symmetry ansatz of Merkt et al. 2015
#     (PRE 92, 012920); the determining system X(g)=0, [f,X]=0 (and the Lie
#     chain X(L_f^k g)=0 and the steady-state field [f_ss,X]=0) is solved in
#     exact rational arithmetic (modular GF(p) + CRT + rational reconstruction,
#     validated, sympy fallback; exact=False selects a floating-point path) and
#     every generator is verified symbolically.
#   "scaling": scaling (toric) symmetries from the integer kernel of the
#     monomial-exponent conditions; pure integer linear algebra.
#
# Observation functions may be non-rational with rational gradient (log10);
# gradients are taken before any polynomial is formed. A clean parse_expr
# symbol table keeps state and parameter names unmangled, and an optional
# symengine backend accelerates the symbolic build.

import io
import sys
import math
import tokenize

import numpy as np
import sympy as spy
from sympy.parsing.sympy_parser import parse_expr

try:
    import symengine as _seng
    _HAVE_SYMENGINE = True
except Exception:
    _HAVE_SYMENGINE = False

# module state set per run
_EXACT = True
_BACKEND = "sympy"
_warned_symengine = False

# primes < 2**31 so int64 products stay < 2**62 during modular elimination
_PRIMES = [2147483647, 2147483629, 2147483587, 2147483579,
           2147483563, 2147483549, 2147483543, 2147483497]


def _select_backend(name):
    global _BACKEND, _warned_symengine
    if name == "symengine" and not _HAVE_SYMENGINE:
        if not _warned_symengine:
            sys.stdout.write("symengine not available; using sympy.\n")
            sys.stdout.flush()
            _warned_symengine = True
        _BACKEND = "sympy"
    else:
        _BACKEND = name


def _diff(expr, var):
    """Differentiate with the selected backend, return a sympy expression."""
    if _BACKEND == "symengine" and _HAVE_SYMENGINE:
        try:
            d = _seng.diff(_seng.sympify(expr), _seng.sympify(var))
            return d._sympy_()
        except Exception:
            return spy.diff(expr, var)
    return spy.diff(expr, var)


###########################################################################
#####################     input parsing (clean)     ######################
###########################################################################

def _as_list(x):
    if x is None:
        return []
    if isinstance(x, str):
        return [x]
    return list(x)


def _clean(line):
    return line.replace('"', '').replace(',', '').replace('^', '**').strip()


def _build_symbol_table(all_lines):
    """Map every non-function identifier to a sympy Symbol, so names like
    E, I, S, O, Q, N are taken as model variables, not sympy constants."""
    syms, funcs = set(), set()
    for raw in all_lines:
        line = _clean(raw)
        if not line:
            continue
        try:
            toks = list(tokenize.generate_tokens(io.StringIO(line).readline))
        except (tokenize.TokenError, IndentationError):
            toks = []
        for i, t in enumerate(toks):
            if t.type == tokenize.NAME:
                nxt = toks[i + 1] if i + 1 < len(toks) else None
                if nxt is not None and nxt.type == tokenize.OP and nxt.string == '(':
                    funcs.add(t.string)
                else:
                    syms.add(t.string)
    syms -= funcs
    return {s: spy.Symbol(s) for s in syms}


def _function_aliases():
    """Functions used in dMod that sympy does not provide by name; mapped to
    differentiable sympy expressions so their gradients are rational."""
    x = spy.Symbol('_x_')
    return {
        'log10': spy.Lambda(x, spy.log(x) / spy.log(10)),
        'log2': spy.Lambda(x, spy.log(x) / spy.log(2)),
    }


def _make_parse(local_dict):
    def parse(rhs):
        return parse_expr(_clean(rhs), local_dict=local_dict, evaluate=True)
    return parse


def _make_local_parse(all_lines):
    """Symbol table (model names as symbols, not sympy constants) with the dMod
    function aliases, plus a parser bound to it."""
    local = _build_symbol_table(all_lines)
    local.update(_function_aliases())
    return local, _make_parse(local)


def _read_equations(lines, parse):
    variables, functions = [], []
    for raw in lines:
        line = _clean(raw)
        if '=' not in line:
            continue
        lhs, rhs = line.split('=', 1)
        variables.append(spy.Symbol(lhs.strip()))
        functions.append(parse(rhs))
    allsyms = set()
    for f in functions:
        allsyms |= set(spy.sympify(f).free_symbols)
    parameters = sorted(allsyms - set(variables), key=spy.default_sort_key)
    return variables, functions, parameters


###########################################################################
#####################     exact polynomial class     #####################
###########################################################################

class Apoly:
    """Sparse multivariate polynomial. Coefficients are kept exact (sympy
    Rational) when _EXACT, else float64 (the floating-point path)."""

    def __init__(self, expr, variables, rs):
        self.vars = variables
        self.rs = rs
        if expr is None:
            self.coefs = []
            self.exps = []
            return
        poly = spy.Poly(spy.sympify(expr), variables).as_dict()
        if rs is None:
            self.coefs = list(poly.values())
        else:
            coefsTmp = list(poly.values())
            self.coefs = [None] * len(coefsTmp)
            for i in range(len(coefsTmp)):
                if _EXACT:
                    row = np.empty(len(rs), dtype=object)
                    for j, r in enumerate(rs):
                        row[j] = spy.diff(coefsTmp[i], r) if coefsTmp[i].has(r) \
                            else spy.Integer(0)
                else:
                    row = np.zeros(len(rs))
                    for j, r in enumerate(rs):
                        if coefsTmp[i].has(r):
                            row[j] = float(spy.diff(coefsTmp[i], r))
                self.coefs[i] = row
        self.exps = [np.array(e) for e in poly.keys()]

    def __repr__(self):
        return str(self.coefs) + '\n' + str(self.exps)

    def getCopy(self):
        newPoly = Apoly(None, self.vars, self.rs)
        newPoly.coefs = [c.copy() if hasattr(c, "copy") else c for c in self.coefs]
        newPoly.exps = [e.copy() for e in self.exps]
        return newPoly

    def add(self, otherPoly):
        for i in range(len(otherPoly.exps)):
            for j in range(len(self.exps)):
                if np.array_equal(otherPoly.exps[i], self.exps[j]):
                    self.coefs[j] = self.coefs[j] + otherPoly.coefs[i]
                    if not np.any(self.coefs[j]):
                        self.coefs.pop(j)
                        self.exps.pop(j)
                    break
            else:
                self.coefs.append(otherPoly.coefs[i])
                self.exps.append(otherPoly.exps[i])

    def sub(self, otherPoly):
        for i in range(len(otherPoly.exps)):
            for j in range(len(self.exps)):
                if np.array_equal(otherPoly.exps[i], self.exps[j]):
                    self.coefs[j] = self.coefs[j] - otherPoly.coefs[i]
                    if not np.any(self.coefs[j]):
                        self.coefs.pop(j)
                        self.exps.pop(j)
                    break
            else:
                self.coefs.append(-1 * otherPoly.coefs[i])
                self.exps.append(otherPoly.exps[i])

    def mul(self, otherPoly):
        newPoly = Apoly(None, self.vars, self.rs)
        n = len(self.coefs) * len(otherPoly.coefs)
        newPoly.coefs = [0] * n
        newPoly.exps = [0] * n
        k = 0
        for i in range(len(otherPoly.exps)):
            for j in range(len(self.exps)):
                # works only because at most one factor carries rs
                newPoly.coefs[k] = otherPoly.coefs[i] * self.coefs[j]
                newPoly.exps[k] = otherPoly.exps[i] + self.exps[j]
                k += 1
        i = 0
        while i < len(newPoly.coefs):
            j = i + 1
            while j < len(newPoly.coefs):
                if np.array_equal(newPoly.exps[i], newPoly.exps[j]):
                    newPoly.exps.pop(j)
                    newPoly.coefs[i] = newPoly.coefs[i] + newPoly.coefs.pop(j)
                else:
                    j += 1
            i += 1
        return newPoly

    def diff(self, j):
        newPoly = self.getCopy()
        i = 0
        while i < len(newPoly.exps):
            if newPoly.exps[i][j] != 0:
                newPoly.coefs[i] = newPoly.coefs[i] * int(newPoly.exps[i][j])
                newPoly.exps[i][j] -= 1
                i += 1
            else:
                newPoly.coefs.pop(i)
                newPoly.exps.pop(i)
        return newPoly

    def as_expr(self):
        expr = 0
        for i in range(len(self.coefs)):
            fact = 1
            for j in range(len(self.vars)):
                fact = fact * self.vars[j] ** int(self.exps[i][j])
            if self.rs is None:
                expr += self.coefs[i] * fact
            else:
                coef = 0
                for j in range(len(self.rs)):
                    coef += self.rs[j] * self.coefs[i][j]
                expr += coef * fact
        return spy.nsimplify(expr)


###########################################################################
#####################     infinitesimal ansatz     #######################
###########################################################################

def _sym(name):
    return spy.Symbol(name)


def giveDegree(vars, i, p, summand, poly, num, k, rs):
    if i == len(vars) - 1:
        rs.append(_sym('r_' + str(vars[k]) + '_' + str(num)))
        poly += rs[-1] * summand * vars[i] ** p
        return poly, num + 1
    else:
        for j in range(p + 1):
            poly, num = giveDegree(vars, i + 1, p - j, summand * vars[i] ** j,
                                   poly, num, k, rs)
    return poly, num


def makeAnsatz(ansatz, allVariables, m, q, pMax, fixed):
    n = len(allVariables)
    rs = []
    infis = []

    if ansatz == 'uni':
        for k in range(n):
            infis.append(spy.sympify(0))
            if allVariables[k] in fixed:
                continue
            for p in range(pMax + 1):
                rs.append(_sym('r_' + str(allVariables[k]) + '_' + str(p)))
                infis[-1] += rs[-1] * allVariables[k] ** p
        diffInfis = [[0] * n]
        for i in range(n):
            diffInfis[0][i] = spy.diff(infis[i], allVariables[i])

    elif ansatz == 'par':
        for k in range(n):
            infis.append(spy.sympify(0))
            if allVariables[k] in fixed:
                continue
            num = 0
            for p in range(pMax + 1):
                vari = list(allVariables[m + q:])
                if k < (m + q):
                    vari.append(allVariables[k])
                    kp = len(vari) - 1
                else:
                    kp = k - (m + q)
                degree, num = giveDegree(vari, 0, p, 1, 0, num, kp, rs)
                infis[-1] += degree
        diffInfis = [[0] * n]
        for i in range(n):
            diffInfis[0][i] = spy.diff(infis[i], allVariables[i])

    elif ansatz == 'multi':
        for k in range(n):
            infis.append(spy.sympify(0))
            if allVariables[k] in fixed:
                continue
            num = 0
            for p in range(pMax + 1):
                if k < m:
                    vari = list(allVariables[:m]) + list(allVariables[m + q:])
                    kp = k
                elif k < m + q:
                    vari = list(allVariables[:])
                    kp = k
                else:
                    vari = list(allVariables[m + q:])
                    kp = k - (m + q)
                degree, num = giveDegree(vari, 0, p, 1, 0, num, kp, rs)
                infis[-1] += degree
        diffInfis = [0] * n
        for i in range(n):
            diffInfis[i] = [0] * n
        for i in range(n):
            for j in range(n):
                diffInfis[i][j] = spy.diff(infis[i], allVariables[j])
    else:
        raise UserWarning("ansatz must be one of 'uni', 'par', 'multi'")

    return infis, diffInfis, rs


def transformInfisToPoly(infis, diffInfis, allVariables, rs):
    n = len(allVariables)
    k = len(diffInfis)
    infisPoly = [Apoly(infis[i], allVariables, rs) for i in range(n)]
    diffInfisPoly = [[0] * n for _ in range(k)]
    for a in range(k):
        for i in range(n):
            diffInfisPoly[a][i] = Apoly(diffInfis[a][i], allVariables, rs)
    return infisPoly, diffInfisPoly


###########################################################################
#####################     determining equations     ######################
###########################################################################

def _quotient_derivatives(numerators, denominators, allVariables):
    m = len(numerators)
    n = len(allVariables)
    derivativesNum = [[None] * n for _ in range(m)]
    for k in range(m):
        for l in range(n):
            d = Apoly(None, allVariables, None)
            d.add(numerators[k].diff(l).mul(denominators[k]))
            d.sub(numerators[k].mul(denominators[k].diff(l)))
            derivativesNum[k][l] = d
    return derivativesNum


def doEquation(k, numerators, denominators, derivativesNum, infis, diffInfis,
               allVariables, rs, ansatz):
    n = len(allVariables)
    m = len(numerators)
    polynomial = Apoly(None, allVariables, rs)
    if ansatz in ('uni', 'par'):
        polynomial.add(diffInfis[0][k].mul(denominators[k]).mul(numerators[k]))
        for i in range(n):
            polynomial.sub(infis[i].mul(derivativesNum[k][i]))
    else:  # multi
        for j in range(m):
            summand = diffInfis[k][j].mul(denominators[k]).mul(numerators[j])
            for l in range(m):
                if l != j:
                    summand = summand.mul(denominators[l])
            polynomial.add(summand)
        for i in range(n):
            summand = infis[i].mul(derivativesNum[k][i])
            for l in range(m):
                if l != k:
                    summand = summand.mul(denominators[l])
            polynomial.sub(summand)
    return list(polynomial.coefs)


def _obs_rows(obsExpr, infisSym, allVariables, rs):
    # X(g) = sum_l xi_l dg/dvar_l = 0. Only the gradient of g enters, so g may
    # be non-rational (e.g. log10) as long as every dg/dvar is rational; the
    # numerator of together(X(g)) is then polynomial in the variables.
    n = len(allVariables)
    Xg = spy.Integer(0)
    for l in range(n):
        d = _diff(obsExpr, allVariables[l])
        if d != 0:
            Xg += infisSym[l] * d
    num, _ = spy.fraction(spy.together(Xg))
    num = spy.expand(num)
    try:
        poly = Apoly(num, allVariables, rs)
    except spy.PolynomialError:
        raise UserWarning(
            "Observable gradient is not rational (e.g. exp/sqrt mixed with "
            "other terms). The liesym engine needs rational gradients; strip "
            "an invertible outer function g = phi(h) to its argument h, or use "
            "method='observability'.")
    return list(poly.coefs)


###########################################################################
#####################     exact / float solvers     ######################
###########################################################################

def _rref_mod_p(A, p):
    """Vectorized int64 Gauss-Jordan over GF(p). Returns (R, pivots)."""
    A = (np.asarray(A, dtype=np.int64) % p)
    nrows, ncols = A.shape
    pivots = []
    r = 0
    for c in range(ncols):
        piv = -1
        for i in range(r, nrows):
            if A[i, c] % p != 0:
                piv = i
                break
        if piv == -1:
            continue
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        for i in range(nrows):
            if i != r and A[i, c] % p != 0:
                A[i] = (A[i] - A[i, c] * A[r]) % p
        pivots.append(c)
        r += 1
        if r == nrows:
            break
    return A, pivots


def _rational_reconstruct(a, m):
    """Recover n/d with a ≡ n*d^{-1} (mod m), |n|,|d| bounded by sqrt(m/2)."""
    a %= m
    if a == 0:
        return spy.Integer(0)
    bound = math.isqrt(m // 2)
    r0, r1 = m, a
    s0, s1 = 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return spy.Rational(r1, s1)


def symRatReconBig(residues, primes):
    """Arbitrary-precision per-row CRT + rational reconstruction, without the u128
    cap of the C++ symRatRecon (which limits the product of ~2**31 primes to 4). Each
    row of `residues` is one coefficient's residues across `primes`; returns
    {'num': [...], 'den': [...]} as decimal strings, with den '0' when a coefficient
    does not rationally reconstruct at the given prime product. Used by the
    gauge-robust reconstruction path, whose log-carrying coefficients can have a
    height beyond the 4-prime (~4.6e18) bound."""
    from sympy.ntheory.modular import crt
    mods = [int(p) for p in primes]
    M = 1
    for m in mods:
        M *= m
    rows = [[int(x) for x in row] for row in np.asarray(residues, dtype=object)]
    num, den = [], []
    for row in rows:
        res = [r % m for r, m in zip(row, mods)]
        x, _ = crt(mods, res)
        val = _rational_reconstruct(int(x), int(M))
        if val is None:
            num.append("0"); den.append("0")
        else:
            val = spy.Rational(val)
            num.append(str(int(val.p))); den.append(str(int(val.q)))
    return {'num': num, 'den': den}


def _modular_nullspace(M):
    """Exact nullspace of sympy Matrix M via multi-prime GF(p) + CRT +
    rational reconstruction. Returns list of sympy column vectors, or None
    if reconstruction/validation fails (caller falls back to sympy)."""
    nrows, ncols = M.shape
    if nrows == 0:
        return [spy.Matrix([1 if j == c else 0 for j in range(ncols)])
                for c in range(ncols)]
    Aint = []
    for i in range(nrows):
        row = [spy.Rational(M[i, j]) for j in range(ncols)]
        L = 1
        for r in row:
            L = spy.ilcm(L, r.q)
        Aint.append([int(r * L) for r in row])

    ref_pivots = None
    free = None
    residues = {}
    mods = []
    for p in _PRIMES:
        Ap = [[x % p for x in row] for row in Aint]
        R, pivots = _rref_mod_p(Ap, p)
        if ref_pivots is None:
            ref_pivots = pivots
            free = [c for c in range(ncols) if c not in pivots]
        elif pivots != ref_pivots:
            continue
        mods.append(p)
        for ki in range(len(ref_pivots)):
            for f in free:
                residues.setdefault((ki, f), []).append(int(R[ki, f]) % p)
        if len(mods) >= 4:
            break

    if not mods:
        return None

    from sympy.ntheory.modular import crt
    exact = {}
    for key, res in residues.items():
        x, Mmod = crt(mods, res)
        val = _rational_reconstruct(int(x), int(Mmod))
        if val is None:
            return None
        exact[key] = val

    basis = []
    for f in free:
        v = [spy.Integer(0)] * ncols
        v[f] = spy.Integer(1)
        for ki, c in enumerate(ref_pivots):
            v[c] = -exact[(ki, f)]
        basis.append(spy.Matrix(v))

    # validate exactly
    for v in basis:
        if not (M * v).is_zero_matrix:
            return None
    return basis


def exactNullspace(rows, ncols):
    """Return basis vectors (list of sympy column vectors)."""
    if not rows:
        return [spy.Matrix([1 if j == c else 0 for j in range(ncols)])
                for c in range(ncols)]
    M = spy.Matrix([[spy.sympify(x) for x in row] for row in rows])
    # the modular solver needs a rational matrix; transcendental constants
    # (e.g. log(10) from a log10 observable) route to the exact sympy nullspace
    if all(bool(e.is_rational) for e in M):
        basis = _modular_nullspace(M)
        if basis is not None:
            return basis
    return M.nullspace()


# ---- floating-point path (exact=False) ----

def legacyNullspace(rows, ncols):
    # Float SVD-based nullspace, valid for rectangular and underdetermined systems.
    from scipy.linalg import null_space
    if not rows:
        base = np.eye(ncols)
    else:
        A = np.array([[float(x) for x in row] for row in rows])
        base = null_space(A)
    return [spy.Matrix([spy.nsimplify(base[i, l], rational=True)
                        for i in range(ncols)])
            for l in range(base.shape[1])]


###########################################################################
#####################     post-processing / report     ###################
###########################################################################

def checkForCommonFactor(infisTmp, allVariables, m):
    factors = []
    for i in range(len(allVariables)):
        if infisTmp[i] != 0:
            fac = spy.factor(infisTmp[i])
            if fac.is_Add:
                factors = [infisTmp[i]]
            elif fac.is_Mul:
                factors = list(fac.args)
            else:
                factors = [fac]
            break
    i = 0
    while i < len(factors):
        if factors[i].is_number:
            factors.pop(i)
        elif factors[i] in allVariables[:m]:
            factors.pop(i)
        elif factors[i].is_Add:
            factors.pop(i)
        elif factors[i].is_Pow:
            if not factors[i].args[0].is_Add:
                factors[i] = factors[i].args[0]
                i += 1
            else:
                factors.pop(i)
        else:
            i += 1
    for i in range(1, len(infisTmp)):
        if infisTmp[i] == 0:
            continue
        fac = spy.factor(infisTmp[i])
        if fac.is_Mul:
            factorsTmp = list(fac.args)
        else:
            factorsTmp = [fac]
        j = 0
        while j < len(factors):
            k = 0
            while k < len(factorsTmp):
                if factorsTmp[k].is_number:
                    factorsTmp.pop(k)
                elif factorsTmp[k] in allVariables[:m]:
                    factorsTmp.pop(k)
                elif factorsTmp[k].is_Add:
                    factorsTmp.pop(k)
                elif factorsTmp[k].is_Pow:
                    if not factorsTmp[k].args[0].is_Add:
                        factorsTmp[k] = factorsTmp[k].args[0]
                        k += 1
                    else:
                        factorsTmp.pop(k)
                else:
                    k += 1
            if factors[j] in factorsTmp:
                j += 1
            else:
                factors.pop(j)
        if len(factors) != 0:
            continue
        else:
            break
    return len(factors) != 0


def buildTransformation(infis, allVariables):
    n = len(allVariables)
    epsilon = spy.Symbol('epsilon')
    transformations = [0] * n
    tType = [False] * 6
    for i in range(n):
        if infis[i] == 0:
            transformations[i] = allVariables[i]
        else:
            poly = spy.Poly(infis[i], allVariables).as_dict()
            monomials = list(poly.keys())
            coefs = list(poly.values())
            if len(monomials) == 1:
                p = None
                broke = False
                for j in range(n):
                    if monomials[0][j] != 0:
                        if j == i and p is None:
                            p = monomials[0][i]
                        elif p is None and monomials[0][j] == 1:
                            p = -1 - j
                        else:
                            transformations[i] = '-?-'
                            tType[0] = True
                            broke = True
                            break
                if not broke:
                    if p is None:
                        transformations[i] = allVariables[i] + epsilon * coefs[0]
                        tType[2] = True
                    elif p <= 0:
                        transformations[i] = allVariables[i] + epsilon * coefs[0] * allVariables[-p - 1]
                        tType[5] = True
                    elif p == 1:
                        transformations[i] = spy.exp(epsilon * coefs[0]) * allVariables[i]
                        tType[1] = True
                    else:
                        transformations[i] = spy.simplify(
                            allVariables[i] / (1 - (p - 1) * epsilon * allVariables[i] ** (p - 1)) ** (spy.sympify(1) / (p - 1)))
                        if p == 2:
                            tType[3] = True
                        else:
                            tType[4] = True
            else:
                transformations[i] = '-?-'
                tType[0] = True

    labels = ['unknown', 'scaling', 'translation', 'MM-like', 'p>2', 'gen. translation']
    parts = [labels[i] for i in range(6) if tType[i]]
    return transformations, 'Type: ' + ', '.join(parts)


def printTransformations(infisAll, allVariables, verified):
    n = len(allVariables)
    print('\n\n' + str(len(infisAll)) + ' transformation(s) found:')
    for l in range(len(infisAll)):
        for i in range(n):
            infisAll[l][i] = spy.nsimplify(infisAll[l][i])
        trans, ttype = buildTransformation(infisAll[l], allVariables)
        flag = ''
        if verified is not None:
            flag = '  [verified]' if verified[l] else '  [UNVERIFIED]'
        print('-' * 60)
        print('#' + str(l + 1) + ': ' + ttype + flag)
        for i in range(n):
            if infisAll[l][i] != 0:
                print('  {0:>14s} : {1:s}  ->  {2:s}'.format(
                    str(allVariables[i]), str(infisAll[l][i]), str(trans[i])))


###########################################################################
#####################     verification     ###############################
###########################################################################

def _is_zero(expr):
    if expr == 0:
        return True
    e = spy.together(spy.sympify(expr))
    num, _ = spy.fraction(e)
    return spy.simplify(spy.expand(num)) == 0


def _verify_generator(infis, allVariables, diffEquations, obsExprs):
    n = len(allVariables)
    # observation / Lie-chain invariance: X(o) = 0
    for o in obsExprs:
        s = 0
        for l in range(n):
            if infis[l] != 0:
                s += infis[l] * spy.diff(o, allVariables[l])
        if not _is_zero(s):
            return False
    # flow conditions [f, X] = 0 for each dynamic field. The field has a
    # component only for the dynamic states (index < len(fields)); parameter
    # and input directions evolve trivially (component 0).
    fields = diffEquations
    mf = len(fields)
    for k in range(mf):
        fk = fields[k]
        bracket = 0
        for l in range(n):
            fl = fields[l] if l < mf else 0
            if fl != 0 and infis[k] != 0:
                bracket += fl * spy.diff(infis[k], allVariables[l])
            if infis[l] != 0 and fk != 0:
                bracket -= infis[l] * spy.diff(fk, allVariables[l])
        if not _is_zero(bracket):
            return False
    return True


###########################################################################
#####################     core driver     ################################
###########################################################################

def _rationalize(functions, allVariables):
    """Split each expression into numerator/denominator Apoly pairs."""
    nums, dens = [], []
    for fexpr in functions:
        rational = spy.together(spy.sympify(fexpr))
        nums.append(Apoly(spy.numer(rational), allVariables, None))
        dens.append(Apoly(spy.denom(rational), allVariables, None))
    return nums, dens


def symmetryDetection(allVariables, diffEquations, obsFunctions,
                      ansatz='uni', pMax=2, inputs=(), fixed=(), lieOrder=0,
                      allTrafos=False, verify=True):
    n = len(allVariables)
    m = len(diffEquations)
    q = len(inputs)

    sys.stdout.write('Preparing equations...')
    sys.stdout.flush()

    infisSym, diffInfis, rs = makeAnsatz(ansatz, allVariables, m, q, pMax, list(fixed))
    infis, diffInfis = transformInfisToPoly(infisSym, diffInfis, allVariables, rs)

    numerators, denominators = _rationalize(diffEquations, allVariables)
    derivativesNum = _quotient_derivatives(numerators, denominators, allVariables)

    # observation invariance plus the Lie-derivative chain X(L_f^k g) = 0
    fieldExprs = [spy.together(spy.sympify(f)) for f in diffEquations]
    obsExprs = []
    for g in obsFunctions:
        expr = spy.sympify(g)
        chain = [expr]
        for _ in range(int(lieOrder)):
            cur = chain[-1]
            lie = 0
            for i in range(m):
                d = _diff(cur, allVariables[i])
                if d != 0:
                    lie += fieldExprs[i] * d
            chain.append(spy.together(lie))
        obsExprs.extend(chain)
    h = len(obsExprs)

    sys.stdout.write('done\nBuilding system...')
    sys.stdout.flush()

    rows = []
    for k in range(m):
        rows.extend(doEquation(k, numerators, denominators, derivativesNum,
                               infis, diffInfis, allVariables, rs, ansatz))
    for k in range(h):
        rows.extend(_obs_rows(obsExprs[k], infisSym, allVariables, rs))

    ncols = len(rs)
    sys.stdout.write('done\nSolving system of size %dx%d (%s)...'
                     % (len(rows), ncols, 'exact' if _EXACT else 'float'))
    sys.stdout.flush()

    if _EXACT:
        basis = exactNullspace(rows, ncols)
    else:
        basis = legacyNullspace(rows, ncols)

    sys.stdout.write('done\nProcessing results...\n')
    sys.stdout.flush()

    infisAll = []
    for v in basis:
        infisTmp = [0] * n
        for i in range(n):
            poly = infis[i].getCopy()
            poly.rs = [v[j] for j in range(ncols)]
            infisTmp[i] = poly.as_expr()
        if allTrafos or not checkForCommonFactor(infisTmp, allVariables, m):
            infisAll.append(infisTmp)

    verified = None
    if verify:
        verified = [_verify_generator(infisTmp, allVariables, fieldExprs, obsExprs)
                    for infisTmp in infisAll]

    printTransformations(infisAll, allVariables, verified)

    # structured return
    result = []
    for l, infisTmp in enumerate(infisAll):
        trans, ttype = buildTransformation(infisTmp, allVariables)
        nz = [i for i in range(n) if infisTmp[i] != 0]
        result.append({
            'infinitesimals': {str(allVariables[i]): str(infisTmp[i]) for i in nz},
            'transformation': {str(allVariables[i]): str(trans[i]) for i in nz},
            'type': ttype,
            'verified': (None if verified is None else bool(verified[l])),
        })
    return result


# raised when an expression is not a rational function of the coordinates
class _NotRational(Exception):
    pass

###########################################################################
#####################   observability tape compiler   ####################
###########################################################################

# Opcodes for the flat tape consumed by the C++ kernel (src/symmetry_kernel.cpp).
# Integer powers are expanded into multiplications, so the kernel needs only
# these four operations.
_OP_CONST, _OP_ADD, _OP_MUL, _OP_INV = 0, 1, 2, 3


def _is_rational_expr(e):
    """True if e is a rational function of its symbols (no transcendental
    function, no non-integer power)."""
    e = spy.sympify(e)
    for a in spy.preorder_traversal(e):
        if a.is_Function:
            return False
        if a.is_Pow and not a.exp.is_Integer:
            return False
    return True


def _emit_tape_shared(fexpr, gexpr, slotOf, base):
    """Emit a straight-line tape over an explicit slot map covering both
    states and leaves; instruction i writes slot base + i. Used by the
    multi-condition compiler, where states occupy their own slots (seeded from
    an initial condition) rather than being leaves."""
    op, a, b, cnum, cden = [], [], [], [], []
    memo = {}

    def emit(opcode, aa, bb, num=0, den=1):
        op.append(int(opcode))
        a.append(int(aa))
        b.append(int(bb))
        cnum.append(str(num))
        cden.append(str(den))
        return base + len(op) - 1

    def build(e):
        s = memo.get(e)
        if s is not None:
            return s
        if e.is_Symbol:
            nm = str(e)
            if nm not in slotOf:
                raise _NotRational()
            memo[e] = slotOf[nm]
            return slotOf[nm]
        if e.is_Integer:
            s = emit(_OP_CONST, 0, 0, int(e), 1)
        elif e.is_Rational:
            s = emit(_OP_CONST, 0, 0, int(e.p), int(e.q))
        elif e.is_Float:
            r = spy.nsimplify(e, rational=True)
            if not r.is_Rational:
                r = spy.Rational(e)
            s = emit(_OP_CONST, 0, 0, int(r.p), int(r.q))
        elif e.is_Add:
            args = list(e.args)
            s = build(args[0])
            for t in args[1:]:
                s = emit(_OP_ADD, s, build(t))
        elif e.is_Mul:
            args = list(e.args)
            s = build(args[0])
            for t in args[1:]:
                s = emit(_OP_MUL, s, build(t))
        elif e.is_Pow and e.exp.is_Integer:
            n = int(e.exp)
            baseS = build(e.base)
            if n == 0:
                s = emit(_OP_CONST, 0, 0, 1, 1)
            else:
                s = baseS
                for _ in range(abs(n) - 1):
                    s = emit(_OP_MUL, s, baseS)
                if n < 0:
                    s = emit(_OP_INV, s, 0)
        else:
            raise _NotRational()
        memo[e] = s
        return s

    repl, red = spy.cse(list(fexpr) + list(gexpr))
    for sym, sub in repl:
        memo[sym] = build(sub)
    outslots = [build(e) for e in red]
    return op, a, b, cnum, cden, outslots


_ssModularCache = {}
# coupled-subsystem solutions keyed by the canonical system mod p, so a call whose
# coupled part is unchanged (e.g. a probe perturbing an unrelated parameter) reuses
# the solution instead of recomputing a Groebner basis
_coupledCache = {}
# compiled t0 event value/derivative term lists keyed by the event tuple and the
# parameter order
_eventCompileCache = {}

# force the symbolic steady-state solve; a test seam to cross-check it against the
# compiled linear plan
_SS_FORCE_SYMPY = False

# force the numeric per-point steady-state seed; a test seam to cross-check it
# against the compiled steady-state IC tape
_FORCE_CONSTRAINT_SEED = False


def setSteadyStateForceSympy(on=True):
    """Force the symbolic steady-state solve and return the previous setting; a test
    seam to cross-check it against the compiled linear plan."""
    global _SS_FORCE_SYMPY
    prev = _SS_FORCE_SYMPY
    _SS_FORCE_SYMPY = bool(on)
    return prev


def setForceConstraintSeed(on=True):
    """Force the numeric per-point steady-state seed and return the previous setting;
    a test seam to cross-check it against the compiled steady-state IC tape."""
    global _FORCE_CONSTRAINT_SEED
    prev = _FORCE_CONSTRAINT_SEED
    _FORCE_CONSTRAINT_SEED = bool(on)
    return prev


def _modp_rational(val, p):
    """Map a sympy Rational (or integer-valued expression) to its residue mod p."""
    r = spy.Rational(val)
    return (int(r.p) % p) * pow(int(r.q) % p, p - 2, p) % p


def _poly_terms(expr, gens):
    """Compile a rational expression in `gens` to (numerator, denominator) term
    lists [(exponent tuple, (p, q))] with rational coefficients as integer pairs.
    Prime-independent, so the compile is cached and only the per-prime reduction
    (in _eval_terms) repeats."""
    num, den = spy.fraction(spy.together(expr))
    def terms(e):
        P = spy.Poly(e, *gens)
        out = []
        for monom, coef in P.terms():
            r = spy.Rational(coef)
            out.append((tuple(int(m) for m in monom), (int(r.p), int(r.q))))
        return out
    return terms(num), terms(den)


def _bipoly(expr, stateGens, paramGens):
    """Compile a polynomial in (states, params) to a list of (state-monomial
    exponent tuple, param-coefficient term list), the latter from _poly_terms, so
    the per-point state-polynomial is rebuilt by evaluating the parameter
    coefficients mod p."""
    P = spy.Poly(expr, *stateGens)
    out = []
    for sm, coef in P.terms():
        ct = _poly_terms(coef, paramGens)
        out.append((tuple(int(m) for m in sm), ct))
    return out


def _eval_bipoly(bip, stateGens, paramvals, p):
    """Rebuild the state-polynomial at numeric parameter values mod p: each
    state-monomial scaled by its parameter coefficient, evaluated via _eval_terms
    at `paramvals` (the coefficients are functions of the parameters only)."""
    terms = []
    for sm, ct in bip:
        c = _eval_terms(ct, paramvals, p)
        if c == 0:
            continue
        mono = spy.Integer(c)
        for k, e in enumerate(sm):
            if e:
                mono = mono * stateGens[k] ** int(e)
        terms.append(mono)
    return spy.Add(*terms) if terms else spy.Integer(0)


def _solve_mod(A, B, p):
    """Solve A X = B over GF(p) by Gauss-Jordan (A: n x n, B: n x k, integer lists
    already reduced mod p). Returns X as a list of rows, or None if A is singular."""
    n = len(A)
    if n == 0:
        return [[] for _ in range(0)]
    k = len(B[0]) if B and B[0] else 0
    M = [list(A[i]) + list(B[i]) for i in range(n)]
    for col in range(n):
        piv = next((r for r in range(col, n) if M[r][col] % p), None)
        if piv is None:
            return None
        M[col], M[piv] = M[piv], M[col]
        inv = pow(M[col][col] % p, p - 2, p)
        M[col] = [(x * inv) % p for x in M[col]]
        for r in range(n):
            if r != col and M[r][col] % p:
                f = M[r][col] % p
                M[r] = [(M[r][c] - f * M[col][c]) % p for c in range(n + k)]
    return [M[i][n:] for i in range(n)]


def _eval_poly_terms(terms, ptvals, p):
    """Value of one term list (from _poly_terms) at the integer point `ptvals`
    (gens order, already mod p) as a residue mod p."""
    acc = 0
    for monom, (cn, cd) in terms:
        c = (cn % p) * pow(cd % p, p - 2, p) % p
        for k, e in enumerate(monom):
            if e:
                c = c * pow(ptvals[k], e, p) % p
        acc = (acc + c) % p
    return acc


def _eval_terms(numden, ptvals, p):
    """Residue num/den mod p of a (num, den) term list at the point `ptvals`."""
    nv = _eval_poly_terms(numden[0], ptvals, p)
    dv = _eval_poly_terms(numden[1], ptvals, p)
    return nv * pow(dv % p, p - 2, p) % p


def _eval_terms_guarded(numden, ptvals, p):
    """Residue num/den mod p, or None when the denominator vanishes mod p. The
    explicit denominator check is required because a Fermat inverse of 0 returns
    0, which would silently corrupt the value."""
    dv = _eval_poly_terms(numden[1], ptvals, p)
    if dv % p == 0:
        return None
    nv = _eval_poly_terms(numden[0], ptvals, p)
    return nv * pow(dv, p - 2, p) % p


def _reduce_modp(e, p):
    """Reduce the integer coefficients of a polynomial expression modulo p."""
    e = spy.sympify(e)
    syms = sorted(e.free_symbols, key=str)
    if not syms:
        return spy.Integer(int(e) % p)
    return spy.Poly(e, *syms, modulus=p).as_expr()


def evalRationalMod(expr, names, vals, q):
    """Value of the rational expression `expr` at the integer point names -> vals,
    reduced modulo the prime q. Returns None if the denominator vanishes mod q.
    Used to certify a reconstructed identifiability direction at a fresh prime."""
    q = int(q)
    if not isinstance(names, (list, tuple)): names = [names]
    if not isinstance(vals, (list, tuple)): vals = [vals]
    local, parse = _make_local_parse([str(expr)] + [str(n) for n in names])
    e = spy.together(spy.sympify(parse(str(expr))).subs(
        {spy.Symbol(str(n)): spy.Integer(int(v))
         for n, v in zip(names, vals)}))
    num, den = spy.fraction(e)
    num, den = int(num), int(den) % q
    if den == 0:
        return None
    return int((num % q) * pow(den, q - 2, q) % q)


_LINSOLVE_OPS_CAP = 2000  # bail to the numeric solver past this elimination size


def _linear_solution(polys, solveStates):
    """Generic linear-elimination solution of polys = 0 for `solveStates`, each a
    rational function of the parameters. Returns {state: expr} or None when the
    system is not generically linear (a coupled residual remains, or an
    intermediate expression exceeds `_LINSOLVE_OPS_CAP` and is left to the numeric
    solver). The order is fixed generically; a point where a pivot denominator
    vanishes mod p is caught at evaluation and routed to the symbolic solve."""
    remSet = set(solveStates)
    remP = [spy.sympify(pl) for pl in polys]
    elim = []
    progress = True
    while progress and remSet:
        progress = False
        for idx in range(len(remP)):
            pl = remP[idx]
            if pl is None or not (pl.free_symbols & remSet):
                continue
            if pl.count_ops() > _LINSOLVE_OPS_CAP:
                return None
            picked = None
            for v in pl.free_symbols & remSet:
                pv = spy.Poly(pl, v)
                if pv.degree() != 1:
                    continue
                a = pv.nth(1)
                if (a.free_symbols & remSet) or a.is_zero:
                    continue
                picked = (v, spy.cancel(-pv.nth(0) / a))
                break
            if picked is None:
                continue
            v, expr = picked
            elim.append((v, expr))
            remP[idx] = None
            for j in range(len(remP)):
                if remP[j] is not None and v in remP[j].free_symbols:
                    sub = remP[j].subs(v, expr)
                    if sub.count_ops() > _LINSOLVE_OPS_CAP:
                        return None
                    remP[j] = spy.fraction(spy.together(sub))[0]
            remSet.discard(v)
            progress = True
    if remSet:
        return None
    sol = {}
    for v, expr in reversed(elim):
        sol[v] = spy.cancel(expr.subs(sol))
    return sol


def _compile_linear_plan(polys, solveStates, paramSyms):
    """Compile the linear-elimination solution of polys = 0 to per-state
    prime-independent (num, den) term lists in `solveStates` order. Returns
    (True, perStateTerms) when every state eliminates linearly, else (False, None)."""
    sol = _linear_solution(polys, solveStates)
    if sol is None:
        return False, None
    return True, [(str(s), _poly_terms(sol[s], paramSyms)) for s in solveStates]


def _coupled_groebner(residPolys, remaining, p):
    """Interior point of the coupled residual system (sympy polynomials in
    `remaining` over GF(p)) by a saturated grevlex Groebner basis, FGLM to a
    triangular lex basis (direct lex as fallback) and interior-root
    back-substitution. Returns (sol-dict, None) or (None, failure-dict)."""
    remaining = list(remaining)
    try:
        ckey = (p, tuple(str(v) for v in remaining), tuple(sorted(
            tuple(sorted((tuple(int(e) for e in m), int(c) % p)
                         for m, c in spy.Poly(pl, *remaining).terms()))
            for pl in residPolys)))
    except Exception:
        ckey = None
    cached = _coupledCache.get(ckey) if ckey is not None else None
    if cached is not None:
        if cached == 'NONE':
            return None, {'ok': False, 'why': 'no consistent interior point'}
        return {spy.Symbol(k): spy.Integer(val) for k, val in cached}, None
    if len(remaining) == 1:
        # one coupled variable: its interior value is a nonzero root of any residual
        # that contains it at which all residuals vanish
        v = remaining[0]
        cand = next((pl for pl in residPolys
                     if v in spy.sympify(pl).free_symbols), None)
        if cand is not None:
            try:
                roots = spy.Poly(cand, v, modulus=p).ground_roots()
            except Exception:
                roots = {}
            residE1 = [spy.sympify(pl) for pl in residPolys]
            for r in roots:
                rr = int(r) % p
                if rr == 0:
                    continue
                if all(int(e.subs(v, rr)) % p == 0 for e in residE1):
                    if ckey is not None:
                        _coupledCache[ckey] = ((str(v), rr),)
                    return {v: spy.Integer(rr)}, None
    w = spy.Symbol('_w_sat_')
    prod = spy.Integer(1)
    for s in remaining:
        prod = prod * s
    sat = residPolys + [w * prod - 1]
    try:
        G = spy.groebner(sat, w, *remaining, order='grevlex', modulus=p)
        basis = [g for g in G.fglm('lex')
                 if w not in spy.sympify(g).free_symbols]
    except Exception:
        try:
            G = spy.groebner(sat, w, *remaining, order='lex', modulus=p)
        except Exception as e:
            return None, {'ok': False, 'why': 'groebner: %s' % e}
        basis = [g for g in G.exprs if w not in g.free_symbols]
    rvars = list(reversed(remaining))
    rset = set(remaining)
    residE = [spy.sympify(pl) for pl in residPolys]

    def backtrack(i, partial):
        if i == len(rvars):
            for e in residE:
                if int(e.subs(partial)) % p != 0:
                    return None
            return partial
        v = rvars[i]
        cand = None
        for g in basis:
            gg = spy.sympify(g).subs(partial)
            if gg.free_symbols & rset == {v}:
                cand = gg
                break
        if cand is None:
            return None
        try:
            roots = spy.Poly(cand, v, modulus=p).ground_roots()
        except Exception:
            return None
        for r in roots:
            rr = int(r) % p
            if rr == 0:
                continue
            nxt = dict(partial)
            nxt[v] = spy.Integer(rr)
            got = backtrack(i + 1, nxt)
            if got is not None:
                return got
        return None

    rsol = backtrack(0, {})
    if rsol is None:
        if ckey is not None:
            _coupledCache[ckey] = 'NONE'
        return None, {'ok': False, 'why': 'no consistent interior point'}
    if ckey is not None:
        _coupledCache[ckey] = tuple((str(k), int(v) % p) for k, v in rsol.items())
    return rsol, None


def _solve_states_modular(polysN, solveStates, p):
    """Interior point of f = 0 over GF(p) from the per-point specialized state
    polynomials `polysN`. Eliminates linearly where possible, then solves the
    coupled residual by a saturated Groebner basis with interior-root
    back-substitution. Returns (sol, None) on success or (None, failure-dict)."""
    sol = {}
    elim = []
    remaining = list(solveStates)
    remSet = set(remaining)
    remPolys = list(polysN)
    progress = True
    while progress and remaining:
        progress = False
        for idx, pl in enumerate(remPolys):
            if pl is None:
                continue
            present = pl.free_symbols & remSet
            if not present:
                continue
            picked = None
            for v in present:
                pv = spy.Poly(pl, v)
                if pv.degree() != 1:
                    continue
                a = pv.nth(1)
                if a.free_symbols & remSet or int(a) % p == 0:
                    continue
                b = pv.nth(0)
                expr = _reduce_modp(-b * pow(int(a) % p, p - 2, p), p)
                picked = (v, expr)
                break
            if picked is None:
                continue
            v, expr = picked
            elim.append((v, expr))
            remPolys[idx] = None
            for j in range(len(remPolys)):
                if remPolys[j] is not None and v in remPolys[j].free_symbols:
                    remPolys[j] = _reduce_modp(remPolys[j].subs(v, expr), p)
            remaining.remove(v)
            remSet.discard(v)
            progress = True

    # coupled residual solved by the saturated Groebner basis
    if remaining:
        residPolys = [pl for pl in remPolys if pl is not None]
        sol, fail = _coupled_groebner(residPolys, remaining, p)
        if fail is not None:
            return None, fail

    # resolve the eliminated states to numbers, latest first
    for v, expr in reversed(elim):
        sol[v] = spy.Integer(_modp_rational(expr.subs(sol), p))
    return sol, None


def _compile_t0events(events, paramNames):
    """Compile each t0 event's value and its parameter derivatives to
    prime-independent term lists over the parameters, so a point evaluates by
    integer arithmetic. Keyed by the event tuple and parameter order."""
    key = (tuple((str(e['var']), str(e['value']), str(e['method'])) for e in events),
           tuple(paramNames))
    compiled = _eventCompileCache.get(key)
    if compiled is not None:
        return compiled
    paramSyms = [spy.Symbol(nm) for nm in paramNames]
    eloc = _build_symbol_table([str(e['value']) for e in events] + list(paramNames))
    eloc.update(_function_aliases())
    eparse = _make_parse(eloc)
    compiled = []
    for evt in events:
        val = spy.sympify(eparse(str(evt['value'])))
        duals = {}
        for nm in paramNames:
            d = spy.diff(val, spy.Symbol(nm))
            duals[nm] = None if d == 0 else _poly_terms(d, paramSyms)
        compiled.append({'var': str(evt['var']), 'method': str(evt['method']),
                         'valT': _poly_terms(val, paramSyms), 'duals': duals})
    _eventCompileCache[key] = compiled
    return compiled


def _ss_compile(model, stateNames, paramNames, forcings):
    """Prime-independent compile of f = 0 for the modular steady-state solve,
    memoised in _ssModularCache. Inputs are the already-normalised model/state/param
    lists and the forcings set. Returns the cached tuple (paramSyms, solveStates,
    polys, Jx, Jt, gens, JxTerms, JtTerms, polyBi, genericLinear, linTerms)."""
    key = (tuple(model), tuple(stateNames), tuple(paramNames), tuple(sorted(forcings)))
    cached = _ssModularCache.get(key)
    if cached is not None:
        return cached
    local, parse = _make_local_parse(model + list(paramNames))
    rhsByName = {}
    for raw in model:
        line = _clean(raw)
        if '=' not in line:
            continue
        lhs, rhs = line.split('=', 1)
        rhsByName[lhs.strip()] = parse(rhs)
    stateSyms = [spy.Symbol(nm) for nm in stateNames]
    paramSyms = [spy.Symbol(nm) for nm in paramNames]
    solveStates = [s for s in stateSyms if str(s) not in forcings]
    forcingSubs = {spy.Symbol(nm): spy.Integer(0) for nm in forcings}
    fSym = [spy.sympify(rhsByName[str(s)]).subs(forcingSubs) for s in solveStates]
    polys = [spy.expand(spy.fraction(spy.together(fi))[0]) for fi in fSym]
    Jx = spy.Matrix([[spy.diff(fi, sj) for sj in solveStates] for fi in fSym])
    Jt = {str(th): spy.Matrix([[spy.diff(fi, th)] for fi in fSym])
          for th in paramSyms}
    gens = list(paramSyms) + list(solveStates)
    nSc = len(solveStates)
    JxTerms = [[_poly_terms(Jx[i, j], gens) for j in range(nSc)]
               for i in range(nSc)]
    JtTerms = {str(th): [_poly_terms(Jt[str(th)][i], gens) for i in range(nSc)]
               for th in paramSyms}
    polyBi = [_bipoly(pl, list(solveStates), list(paramSyms)) for pl in polys]
    genericLinear, linTerms = _compile_linear_plan(polys, solveStates, paramSyms)
    cached = (paramSyms, solveStates, polys, Jx, Jt, gens, JxTerms, JtTerms,
              polyBi, genericLinear, linTerms)
    _ssModularCache[key] = cached
    return cached


def solveSteadyStateModular(model, stateNames, paramNames, paramVals, prime,
                            forcings=None, backend='sympy', t0events=None,
                            recast=None, lVals=None, jointMode=False):
    """Numeric point on the interior component of f = 0 over GF(prime) with its
    implicit-function-theorem parameter sensitivities.

    `model` is a list of "X = rhs" lines; `stateNames`/`paramNames` the state and
    parameter order; `paramVals` a {name: int} map of residues mod `prime`;
    `forcings` are held at 0; `t0events` are dose-style events composed onto the
    point after the solve. Returns {'ok': True, 'xstar': [...], 'dx': {param:
    [...]}} (state-ordered) or {'ok': False, 'why': ...}."""
    _select_backend(backend)
    model = _as_list(model)
    forcings = set(_as_list(forcings))
    stateNames = _as_list(stateNames)
    paramNames = _as_list(paramNames)
    p = int(prime)

    cached = _ss_compile(model, stateNames, paramNames, forcings)
    (paramSyms, solveStates, polys, Jx, Jt, gens, JxTerms, JtTerms, polyBi,
     genericLinear, linTerms) = cached

    paramvals = [int(paramVals.get(str(th), 0)) % p for th in paramSyms]

    # fast path: each state is a compiled rational function of the parameters,
    # evaluated in integer arithmetic; a vanishing pivot denominator mod p (None)
    # routes the point to the symbolic solve
    sol = None
    if genericLinear and not _SS_FORCE_SYMPY:
        sol = {}
        for nm, terms in linTerms:
            val = _eval_terms_guarded(terms, paramvals, p)
            if val is None:
                sol = None
                break
            sol[spy.Symbol(nm)] = spy.Integer(val)
    if sol is None:
        polysN = [_eval_bipoly(bip, list(solveStates), paramvals, p) for bip in polyBi]
        sol, fail = _solve_states_modular(polysN, solveStates, p)
        if fail is not None:
            return fail

    nS = len(solveStates)
    valBy = {str(s): int(sol[s]) % p for s in solveStates}
    sIdx = {str(s): i for i, s in enumerate(solveStates)}
    # point in gens order (params then solve-states), already reduced mod p
    ptvals = [int(paramVals.get(str(th), 0)) % p for th in paramSyms] + \
             [valBy[str(s)] for s in solveStates]
    JtBy = {str(th): [_eval_terms(JtTerms[str(th)][i], ptvals, p) for i in range(nS)]
            for th in paramSyms}
    Jxeff = [[_eval_terms(JxTerms[i][j], ptvals, p) for j in range(nS)]
             for i in range(nS)]

    # raw resting Jacobian rows for the implicit/joint determining system: the
    # constraint df_rest . xi = 0 (tangency to the resting manifold) in the
    # UNELIMINATED (x, theta) coordinates. Snapshot BEFORE the recast folding below
    # mutates Jxeff/JtBy. Row i is d f_rest[i] / d(solveStates, params). State
    # columns are Jxeff, parameter columns are JtBy (which already carry the recast
    # gen coordinate E as a param). Returned so R can stack these rows and run the
    # scaling peel over the enlarged coordinate set (dfStateCols / dfParamCols name
    # the columns).
    dfJx = [list(row) for row in Jxeff]
    dfJt = {nm: list(col) for nm, col in JtBy.items()}
    dfStateCols = [str(sst) for sst in solveStates]
    dfParamCols = [str(th) for th in paramSyms]

    if jointMode:
        # joint/implicit mode uses only the pre-event resting value (valBy) and the
        # raw constraint Jacobian df_rest; the IFT parameter-duals, the recast dual
        # chain and the t0-event composition below are all for the eliminated icSeed
        # and are unused here. Returning now also skips the (occasionally singular)
        # dual solve, so a point that seeds fine is not rejected for a dual failure.
        return {'ok': True, 'stateNames': list(stateNames), 'valBy': dict(valBy),
                'dfJx': dfJx, 'dfJt': dfJt, 'dfStateCols': dfStateCols,
                'dfParamCols': dfParamCols}

    # power/Hill recast. A normal entry holds E = base^exp generic and folds its
    # chain rule into the base and exponent columns, so the IFT duals of the solved
    # base pick up the exponent. An inverted entry has E solved (the balance is
    # linear in E) and holds base and L = log(base) generic; no folding is done and
    # their duals come from dE by the inverse chain rule. base0, E0, L0 are the
    # resting values; the generic partner of each entry is an independent residue.
    recastOut = []
    if recast:
        lVals = lVals or {}
        for rc in recast:
            E, L, base, exp = rc['E'], rc['L'], rc['base'], rc['exp']
            inverted = bool(rc.get('inverted'))
            E0 = valBy[E] if E in valBy else int(paramVals.get(E, 0)) % p
            L0 = int(lVals.get(L, 0)) % p
            base0 = valBy[base] if base in valBy else int(paramVals.get(base, 0)) % p
            expv = int(paramVals.get(exp, 0)) % p
            if not inverted:
                binv = pow(base0, p - 2, p) if base0 % p else 0
                dEbase = expv * E0 % p * binv % p
                dEexp = E0 * L0 % p
                jE = JtBy.get(E, [0] * nS)
                if base in sIdx:
                    bc = sIdx[base]
                    for i in range(nS):
                        Jxeff[i][bc] = (Jxeff[i][bc] + jE[i] * dEbase) % p
                elif base in JtBy:
                    JtBy[base] = [(JtBy[base][i] + jE[i] * dEbase) % p for i in range(nS)]
                if exp in JtBy:
                    JtBy[exp] = [(JtBy[exp][i] + jE[i] * dEexp) % p for i in range(nS)]
            recastOut.append({'E': E, 'L': L, 'base': base, 'exp': exp,
                              'inverted': inverted, 'E0': E0, 'L0': L0,
                              'base0': base0, 'expv': expv})

    # solve Jxeff * dx = -Jt for every parameter column at once over GF(p)
    params_l = [str(th) for th in paramSyms]
    B = [[(-JtBy[params_l[j]][i]) % p for j in range(len(params_l))]
         for i in range(nS)]
    X = _solve_mod(Jxeff, B, p)
    if X is None:
        return {'ok': False, 'why': 'singular jacobian mod p'}
    dxBy = {params_l[j]: [X[i][j] % p for i in range(nS)]
            for j in range(len(params_l))}

    # assemble in the full state order; held states stay at 0 with no duals
    xstar = [valBy.get(nm, 0) for nm in stateNames]
    dx = {}
    for nm in paramNames:
        col = dxBy.get(nm, [0] * nS)
        dx[nm] = [col[sIdx[s]] if s in sIdx else 0 for s in stateNames]

    # parameter-duals for the recast icSeed rows. Each entry seeds its generic
    # partner (gen) and L: a normal entry derives them from dx[base], an inverted
    # one from dx[E]. gen is E for a normal entry and base for an inverted one.
    for rc in recastOut:
        E, L, base, exp = rc['E'], rc['L'], rc['base'], rc['exp']
        E0, L0, base0, expv = rc['E0'], rc['L0'], rc['base0'], rc['expv']
        gD = {}
        lD = {}
        if not rc['inverted']:
            bi = sIdx.get(base)
            binv = pow(base0, p - 2, p) if base0 % p else 0
            dEbase = expv * E0 % p * binv % p
            dEexp = E0 * L0 % p
            for nm in paramNames:
                dxb = dxBy[nm][bi] if (bi is not None and nm in dxBy) else 0
                gD[nm] = (dEbase * dxb + (dEexp if nm == exp else 0)) % p
                lD[nm] = binv * dxb % p
            rc['gen'], rc['gen0'] = E, E0
        else:
            ei = sIdx.get(E)
            xe = pow(expv, p - 2, p) if expv % p else 0
            einv = pow(E0, p - 2, p) if E0 % p else 0
            for nm in paramNames:
                dE = dxBy[nm][ei] if (ei is not None and nm in dxBy) else 0
                gD[nm] = (base0 * xe % p * einv % p * dE
                          - (base0 * L0 % p * xe if nm == exp else 0)) % p
                lD[nm] = (xe * einv % p * dE - (L0 * xe if nm == exp else 0)) % p
            rc['gen'], rc['gen0'] = base, base0
        rc['genDual'] = gD
        rc['lDual'] = lD
        for k in ('E', 'base', 'exp', 'E0', 'base0', 'expv'):
            rc.pop(k, None)

    # compose the t0 events onto the seed (value mod p and its parameter-duals)
    events = list(t0events) if t0events else []
    if events:
        idxOfState = {nm: i for i, nm in enumerate(stateNames)}
        for evt in _compile_t0events(events, list(paramNames)):
            X = evt['var']
            if X not in idxOfState:
                continue
            i = idxOfState[X]
            v0 = _eval_terms(evt['valT'], paramvals, p)
            vdu = {nm: (_eval_terms(evt['duals'][nm], paramvals, p)
                        if evt['duals'][nm] is not None else 0)
                   for nm in paramNames}
            meth = evt['method']
            if meth == 'replace':
                xstar[i] = v0
                for nm in paramNames:
                    dx[nm][i] = vdu[nm]
            elif meth == 'add':
                xstar[i] = (xstar[i] + v0) % p
                for nm in paramNames:
                    dx[nm][i] = (dx[nm][i] + vdu[nm]) % p
            elif meth == 'multiply':
                old = xstar[i]
                oldd = {nm: dx[nm][i] for nm in paramNames}
                xstar[i] = old * v0 % p
                for nm in paramNames:
                    dx[nm][i] = (old * vdu[nm] + v0 * oldd[nm]) % p
    return {'ok': True, 'xstar': xstar, 'dx': dx, 'stateNames': list(stateNames),
            'recast': recastOut,
            'dfJx': dfJx, 'dfJt': dfJt, 'dfStateCols': dfStateCols,
            'dfParamCols': dfParamCols, 'valBy': dict(valBy)}


# ---- forward steady-state solve (choose the resting states, solve for rates) ----------
# The backward solve above finds x*(theta); the FORWARD solve inverts it: CHOOSE the resting
# state values and the free parameters, and solve f = 0 for a turnover subset of the rate
# constants. Because every mass-action rate enters f linearly (and no two rates multiply),
# f = 0 is a LINEAR system in the chosen rates -- so it is solvable at EVERY prime (no
# per-prime Groebner degeneracy), which lets the gauge-robust reconstruction sample a
# direction on a SHARED slice across primes (the states become independent coordinates, so a
# direction whose entries depend on x* -- e.g. via a Hill term C3^n -- reconstructs as a
# rational in (theta, states) instead of an inconsistent per-prime constant).

def _forward_rate_pick(rhsByName, solveStates, paramNames, forcings, keepFree=None):
    """Match each state to one rate constant to solve for: a parameter that enters that state's
    balance LINEARLY (f = 0 stays linear in it), preferring a turnover term (the rate multiplies
    a monomial containing the state, with a negative sign) that is DEDICATED (appears in few
    balances). Parameters in `keepFree` are never solved for -- pass the direction's support so
    the reconstruction can vary them. Assigns the most-constrained states first. Returns the
    rate-name list (one per solve state) or None if no complete matching exists."""
    keepFree = set(keepFree or [])
    # never solve for a recast coordinate (E = base^exp, L = log base): they are independent
    # generic coordinates the reconstruction varies, and choosing one makes f bilinear in the
    # rates (E multiplies a real turnover rate in the recast balance).
    paramset = {spy.Symbol(pn) for pn in paramNames
                if pn not in keepFree and not pn.startswith('_E_') and not pn.startswith('_L_')}
    forc = {spy.Symbol(nm) for nm in forcings}
    cand, appear = {}, {}
    for s in solveStates:
        f = spy.sympify(rhsByName[s]); sSym = spy.Symbol(s); lst = []
        for r in sorted(f.free_symbols & paramset, key=str):
            if r in forc:
                continue
            try:
                pv = spy.Poly(f, r)
            except spy.PolynomialError:
                continue
            if pv.degree() != 1 or r in pv.nth(1).free_symbols:
                continue
            coef = pv.nth(1)
            turnover = sSym in coef.free_symbols
            neg = coef.could_extract_minus_sign()
            lst.append((0 if (turnover and neg) else 1 if turnover else 2, str(r)))
            appear[str(r)] = appear.get(str(r), 0) + 1
        cand[s] = lst
    order = sorted(solveStates, key=lambda s: len(cand[s]))
    used, assign = set(), {}
    for s in order:
        opts = sorted((sc, appear.get(r, 99), r) for sc, r in cand[s] if r not in used)
        if not opts:
            return None
        assign[s] = opts[0][2]; used.add(opts[0][2])
    return [assign[s] for s in solveStates]


_forwardCache = {}


def _forward_compile(model, stateNames, paramNames, forcings, solveRates):
    """Prime-independent compile of the forward solve, cached: the coefficient of each solve
    rate in each state balance and the rate-free constant, as (num, den) term lists over the
    free coords rgens = (non-solve params) + solve states, PLUS the reused steady-state Jacobian
    term lists (JxTerms/JtTerms). So a forward solve is per-point integer arithmetic, not
    symbolic. Returns {'bad': reason} if the rate set is not linear (two rates multiply, a
    concentration chosen, etc.)."""
    key = (tuple(model), tuple(stateNames), tuple(paramNames), tuple(sorted(forcings)),
           tuple(solveRates))
    c = _forwardCache.get(key)
    if c is not None:
        return c
    (paramSyms, solveStates, polys, Jx, Jt, gens0, JxTerms, JtTerms,
     polyBi, genLin, linTerms) = _ss_compile(model, stateNames, paramNames, forcings)
    solveSet = set(solveRates); rSyms = [spy.Symbol(r) for r in solveRates]
    rgens = [th for th in paramSyms if str(th) not in solveSet] + list(solveStates)
    rgenset = set(rgens)
    coefT, constT = [], []
    for i in range(len(solveStates)):
        fexp = spy.sympify(polys[i]); row = []
        for r in rSyms:
            dexp = spy.diff(fexp, r)
            if dexp.free_symbols - rgenset:
                return {'bad': 'not linear in the chosen rates at %s' % str(solveStates[i])}
            row.append(_poly_terms(dexp, rgens))
        cexp = fexp.subs({r: spy.Integer(0) for r in rSyms})
        if cexp.free_symbols - rgenset:
            return {'bad': 'unsubstituted symbol in balance %s' % str(solveStates[i])}
        coefT.append(row); constT.append(_poly_terms(cexp, rgens))
    c = {'paramNames': [str(s) for s in paramSyms], 'solveStates': [str(s) for s in solveStates],
         'rgens': [str(g) for g in rgens], 'coefT': coefT, 'constT': constT,
         'JxTerms': JxTerms, 'JtTerms': JtTerms, 'gens': [str(g) for g in gens0]}
    _forwardCache[key] = c
    return c


def solveForwardModular(model, stateNames, paramNames, stateVals, paramVals, prime,
                        forcings=None, solveRates=None, keepFree=None, backend='sympy'):
    """Solve f = 0 over GF(prime) for a turnover subset of rate constants, given CHOSEN resting
    `stateVals` and `paramVals` (residues mod prime; forcings held at 0). `solveRates` names the
    rates to solve (one per non-forcing state) -- auto-picked (avoiding `keepFree`) if None.
    Uses the cached `_forward_compile`, so it is per-point integer arithmetic. Returns the SAME
    joint-mode payload as the backward solve -- {'valBy','dfJx','dfJt','dfStateCols',
    'dfParamCols'} at the forward point -- plus {'rates','solveRates'}, or {'ok': False,'why'}.
    Linear in the rates, so it solves at every prime unless a pivot vanishes there."""
    _select_backend(backend)
    p = int(prime)
    model = _as_list(model); forcings = set(_as_list(forcings))
    stateNames = _as_list(stateNames); paramNames = _as_list(paramNames)
    svd = {str(k): int(v) % p for k, v in dict(stateVals).items()}
    pvd = {str(k): int(v) % p for k, v in dict(paramVals).items()}
    if solveRates is None:
        (_ps, _ss, _polys, *_r) = _ss_compile(model, stateNames, paramNames, forcings)
        rhsByName = {str(_ss[i]): _polys[i] for i in range(len(_ss))}
        solveRates = _forward_rate_pick(rhsByName, [str(s) for s in _ss], paramNames, forcings, keepFree=keepFree)
        if solveRates is None:
            return {'ok': False, 'why': 'no complete forward rate matching'}
    solveRates = _as_list(solveRates)
    c = _forward_compile(model, stateNames, paramNames, forcings, solveRates)
    if 'bad' in c:
        return {'ok': False, 'why': c['bad']}
    solveStates = c['solveStates']; nS = len(solveStates)
    if len(solveRates) != nS:
        return {'ok': False, 'why': 'rate/state count mismatch (%d rates, %d states)'
                % (len(solveRates), nS)}
    def cvfree(nm):
        return svd[nm] if nm in svd else (pvd[nm] if nm in pvd else 0)
    rgv = [cvfree(nm) for nm in c['rgens']]
    A = [[_eval_terms(c['coefT'][i][j], rgv, p) for j in range(len(solveRates))] for i in range(nS)]
    b = [[(-_eval_terms(c['constT'][i], rgv, p)) % p] for i in range(nS)]
    X = _solve_mod(A, b, p)
    if X is None:
        return {'ok': False, 'why': 'singular forward system mod p'}
    ratesDict = {solveRates[j]: int(X[j][0]) % p for j in range(len(solveRates))}
    def cvfull(nm):
        if nm in ratesDict: return ratesDict[nm]
        return svd[nm] if nm in svd else (pvd[nm] if nm in pvd else 0)
    ptvals = [cvfull(nm) for nm in c['gens']]
    dfJx = [[_eval_terms(c['JxTerms'][i][j], ptvals, p) for j in range(nS)] for i in range(nS)]
    dfJt = {th: [_eval_terms(c['JtTerms'][th][i], ptvals, p) for i in range(nS)]
            for th in c['paramNames']}
    return {'ok': True, 'rates': ratesDict, 'solveRates': list(solveRates),
            'solveStates': list(solveStates), 'valBy': {s: cvfull(s) for s in solveStates},
            'dfJx': dfJx, 'dfJt': dfJt, 'dfStateCols': list(solveStates),
            'dfParamCols': list(c['paramNames'])}


def _detect_power_atoms(perCond):
    """Find base^exp terms with a non-numeric exponent. base must be a single
    symbol and exp = c*n (c rational, n a symbol); returns the unique (base, n)
    pairs, or None if a power is outside this form (then it stays non-rational)."""
    pairs = set()
    for (f_c, g_c, ic_c, f_ss) in perCond:
        for e in list(f_c) + list(g_c) + list(ic_c.values()) + list(f_ss):
            for pw in spy.sympify(e).atoms(spy.Pow):
                if pw.exp.is_number:
                    continue
                if not pw.base.is_Symbol:
                    return None
                c, rest = pw.exp.as_coeff_Mul()
                if not (rest.is_Symbol and c.is_rational):
                    return None
                pairs.add((pw.base, rest))
    return sorted(pairs, key=lambda t: (str(t[0]), str(t[1])))


def _apply_power_recast(pairs, S, perCond):
    """Recast each base^(c*n) as E^c with E = base^n a new state, and add a
    companion L = log(base) per base. E' = n*E*base'/base, L' = base'/base when
    base is a state (0 when it is a parameter). E, L are appended to the state
    list; f_ss keeps only the real states (E is held generic in the solve)."""
    Sset = set(S)
    recast = {}
    Lof = {}
    for (base, n) in pairs:
        recast[(base, n)] = spy.Symbol('_E_%s_%s' % (base, n))
        if base not in Lof:
            Lof[base] = spy.Symbol('_L_%s' % base)

    def sub(expr):
        def repl(pw):
            if pw.exp.is_number or not pw.base.is_Symbol:
                return pw
            c, rest = pw.exp.as_coeff_Mul()
            E = recast.get((pw.base, rest))
            return E ** c if E is not None else pw
        return spy.sympify(expr).replace(lambda x: x.is_Pow, repl)

    extra = [recast[(b, n)] for (b, n) in pairs] + \
            [Lof[b] for b in sorted(Lof, key=str)]
    newPerCond = []
    for (f_c, g_c, ic_c, f_ss) in perCond:
        f_r = [sub(e) for e in f_c]
        g_r = [sub(e) for e in g_c]
        ic_r = {k: sub(v) for k, v in ic_c.items()}
        ss_r = [sub(e) for e in f_ss]
        rhsOf = {str(X): f_r[i] for i, X in enumerate(S)}
        for (b, n) in pairs:
            E = recast[(b, n)]
            brhs = rhsOf.get(str(b), spy.Integer(0))
            ic_r[str(E)] = E
            f_r.append(n * E * brhs / b if brhs != 0 else spy.Integer(0))
        for b in sorted(Lof, key=str):
            L = Lof[b]
            brhs = rhsOf.get(str(b), spy.Integer(0))
            ic_r[str(L)] = L
            f_r.append(brhs / b if brhs != 0 else spy.Integer(0))
        newPerCond.append((f_r, g_r, ic_r, ss_r))
    meta = [{'E': str(recast[(b, n)]), 'L': str(Lof[b]),
             'base': str(b), 'exp': str(n)} for (b, n) in pairs]
    return list(S) + extra, newPerCond, meta


def recastBacksub(expr, eNames, lNames, bases, exps):
    """Substitute the recast coordinates of a reconstructed direction back to
    their meaning, E -> base**exp and L -> log(base), and cancel. `expr` is a
    rational expression string; the name vectors are aligned per recast atom.
    Every name goes through the shared symbol table so model names like E, I, N
    are taken as symbols, not sympy constants."""
    asList = lambda v: list(v) if isinstance(v, (list, tuple)) else [v]
    eNames, lNames, bases, exps = (asList(eNames), asList(lNames),
                                   asList(bases), asList(exps))
    local, parse = _make_local_parse(
        [str(expr)] + [str(x) for x in eNames + lNames + bases + exps])
    e = parse(str(expr))
    subs = {}
    for E, L, base, exp in zip(eNames, lNames, bases, exps):
        b = parse(str(base))
        subs[parse(str(E))] = b ** parse(str(exp))
        subs[parse(str(L))] = spy.log(b)
    return str(spy.cancel(e.subs(subs)))


def compileObservabilityTapeMulti(model, observation, conditionSubs, conditionIC0,
                                  fixed=None, parameters=None, backend='sympy',
                                  equilibrate=False, forcings=None,
                                  segEquilibrate=None, conditionEvents=None,
                                  conditionT0Events=None, jointSteadyState=False,
                                  jointFixedStates=None):
    """Compile one observability tape per experimental condition over a shared
    coordinate space, for the multi-condition observability path.

    The base `model` and `observation` are symbolic. `conditionSubs` is a list
    (one entry per condition) of {symbol: replacement} maps: a numeric
    replacement bakes that symbol to a constant in the condition, a symbol
    replacement renames it (e.g. a knockdown-specific rate), and pre-equilibration
    switches enter here too. `conditionIC0` is a list of {state: expression} maps
    giving the start-point (t0+) initial condition of each state in each
    condition, already composed in R from the steady state / `initial` and the
    t0 events; an entry may be a symbol, a number, or an arbitrary rational
    expression in the parameters. A state with no entry starts free (its own
    unknown initial value).

    Returns the per-condition tapes (each carrying a small IC tape that seeds the
    state initial values, with their parameter-duals, at order 0) plus the shared
    leaf/state layout, the dual-carrying coordinates z = free state initial values
    + free parameters (excluding `fixed`), and their names. A non-rational
    right-hand side, observable or initial condition returns
    {'ok': False, 'nonrational': ...}."""
    _select_backend(backend)
    model = _as_list(model)
    observation = _as_list(observation)
    conditionSubs = conditionSubs or []
    conditionIC0 = conditionIC0 or []
    forcings = set(_as_list(forcings))
    K = len(conditionSubs)
    # one flag per (condition, segment): the first segment of an equilibrated
    # condition is seeded from the resting steady state (icSeed); every later
    # segment, and every segment without equilibrate, is seeded from an IC tape
    # whose initial-value expressions may reference fresh carry coordinates.
    if segEquilibrate is None:
        segEq = [bool(equilibrate)] * K
    else:
        segEq = [bool(x) for x in list(segEquilibrate)]
        segEq += [bool(equilibrate)] * (K - len(segEq))

    extra = []
    for d in conditionSubs:
        for k, v in dict(d).items():
            extra += [str(k), str(v)]
    for d in conditionIC0:
        for k, v in dict(d).items():
            extra += [str(k), str(v)]
    all_lines = model + observation + _as_list(parameters) + extra
    local, parse = _make_local_parse(all_lines)

    variables, diffEquations, _ = _read_equations(model, parse)
    obsVars, obsFunctions, _ = _read_equations(observation, parse)
    fixedNames = set(str(s) for s in
                     [local.get(nm, spy.Symbol(nm)) for nm in _as_list(fixed)])

    isConst = [spy.sympify(e) == 0 for e in diffEquations]
    S = [variables[i] for i in range(len(variables)) if not isConst[i]]
    Srhs = [spy.sympify(diffEquations[i]) for i in range(len(variables))
            if not isConst[i]]
    nS = len(S)

    def pval(s):
        return spy.sympify(parse(str(s)))

    perCond = []
    for c in range(K):
        subsMap = {}
        for k, v in dict(conditionSubs[c]).items():
            subsMap[local.get(str(k), spy.Symbol(str(k)))] = pval(v)
        ic0 = dict(conditionIC0[c]) if c < len(conditionIC0) else {}
        f_c = [e.subs(subsMap) for e in Srhs]
        g_c = [spy.sympify(e).subs(subsMap) for e in obsFunctions]
        # resting-state model for the equilibrate solve: forcings at 0, no events,
        # non-forcing per-condition substitutions baked in
        forcZero = {local.get(nm, spy.Symbol(nm)): spy.Integer(0) for nm in forcings}
        subsMapNF = {k: v for k, v in subsMap.items() if str(k) not in forcings}
        f_ss = [e.subs(subsMapNF).subs(forcZero) for e in Srhs]
        ic_c = {}
        for X in S:
            e = pval(ic0[str(X)]) if str(X) in ic0 else X
            ic_c[str(X)] = spy.sympify(e).subs(subsMap)
        perCond.append((f_c, g_c, ic_c, f_ss))

    # power/Hill recast (equilibrate only): replace base^exp (exp a parameter) by
    # a state E with E' = exp*E*base'/base and a companion L = log(base). Both are
    # held generic in the f = 0 solve and excluded from z, so f stays rational and
    # the log enters only through the generic L coordinate (sound by the algebraic
    # independence of base, log(base) and base^exp).
    nReal = nS
    powerRecast = []
    invSolveName = {}
    if equilibrate:
        pairs = _detect_power_atoms(perCond)
        if pairs is None:
            return {'ok': False, 'nonrational':
                    ['unsupported power form: base must be a symbol and exponent c*param']}
        if pairs:
            S, perCond, powerRecast = _apply_power_recast(pairs, S, perCond)
            nS = len(S)
            # a base with a linear turnover term keeps its bare symbol in the
            # steady-state balance and is solved directly; a base without one has a
            # balance linear in E and is "inverted": E is solved, base and L stay
            # generic with duals from the inverse chain rule.
            realSet = {str(X) for X in S[:nReal]}
            bareSyms = set()
            for (_f, _g, _ic, ss_r) in perCond:
                for e in ss_r[:nReal]:
                    bareSyms |= spy.fraction(spy.together(e))[0].free_symbols
            for rc in powerRecast:
                rc['inverted'] = (rc['base'] in realSet
                                  and spy.Symbol(rc['base']) not in bareSyms)
            byBase = {}
            for rc in powerRecast:
                byBase.setdefault(rc['base'], []).append(rc)
            for base, group in byBase.items():
                if group[0]['inverted'] and len(group) > 1:
                    return {'ok': False, 'nonrational':
                            ['free exponent base %s carries multiple independent '
                             'powers with no linear turnover term' % base]}
            invSolveName = {rc['base']: rc['E']
                            for rc in powerRecast if rc['inverted']}

    nonrational = []
    for (f_c, g_c, ic_c, f_ss) in perCond:
        for X, e in zip(S, f_c):
            if not _is_rational_expr(e):
                nonrational.append('d%s/dt = %s' % (X, e))
        for y, e in zip(obsVars, g_c):
            if not _is_rational_expr(e):
                nonrational.append('%s = %s' % (y, e))
        for X in S:
            if not _is_rational_expr(ic_c[str(X)]):
                nonrational.append('%s(0) = %s' % (X, ic_c[str(X)]))
    if nonrational:
        return {'ok': False, 'nonrational': nonrational}

    # a state carries a free initial-value leaf when its initial condition still
    # depends on the state symbol itself (a bare free value, or a free value with
    # an additive/multiplicative dose composed on top). In constraint mode no
    # state is a free coordinate: every state is seeded numerically from the
    # interior steady-state point per evaluation point and prime (icSeed).
    freeState = {str(X): False for X in S}
    for c, (f_c, g_c, ic_c, f_ss) in enumerate(perCond):
        if segEq[c]:
            # joint/implicit mode: keep every equilibrate state as a free coordinate
            # (its value is seeded on-manifold to x* by R, the steady-state constraint
            # enters as stacked df rows) so the direction is low-degree in (x, theta).
            if jointSteadyState:
                jfs = set(jointFixedStates or [])
                for X in S:
                    if str(X) not in forcings and str(X) not in jfs:
                        freeState[str(X)] = True
            continue  # equilibrate-seeded: no state is a free coordinate (unless joint)
        for X in S:
            if X in spy.sympify(ic_c[str(X)]).free_symbols:
                freeState[str(X)] = True

    # boundary state-dose values (e.g. a second-dose parameter) may appear only in
    # an event, never in f/g/ic; collect their symbols so they become coordinates
    conditionEvents = conditionEvents or []
    eventVals = []
    for c in range(K):
        evs = list(conditionEvents[c]) if c < len(conditionEvents) else []
        eventVals.append([spy.sympify(pval(e['value'])) for e in evs])

    # f_ss is included so parameters the perturbed dynamics drop still count
    paramset = set()
    for (f_c, g_c, ic_c, f_ss) in perCond:
        for e in list(f_c) + list(g_c) + list(ic_c.values()) + list(f_ss):
            paramset |= set(spy.sympify(e).free_symbols)
    for vs in eventVals:
        for e in vs:
            paramset |= set(e.free_symbols)
    paramset -= set(S)
    params = sorted(paramset, key=spy.default_sort_key)
    # A free-exponent base that appears ONLY under the exponent (a Michaelis constant
    # Km in C^n/(Km^n+C^n)) was replaced by its recast atom E=Km^n and so dropped from
    # the coordinates. Re-add it as a parameter coordinate so the recast relation
    # E=base^exp ties it and the pool/Michaelis co-scaling is reported and closes as an
    # exact scaling (rather than a doomed finite-field fit that under-reports Km). A state
    # base (e.g. C3) is already a coordinate; forcings never scale.
    if powerRecast:
        known = set(str(s) for s in S) | set(str(s) for s in params) | set(forcings)
        for rc in powerRecast:
            if rc['base'] not in known:
                params.append(spy.Symbol(rc['base'])); known.add(rc['base'])
        params = sorted(params, key=spy.default_sort_key)

    freeStates = [X for X in S if freeState[str(X)]]
    leafNames = [str(X) for X in freeStates] + [str(s) for s in params]
    leafSlot = {nm: i for i, nm in enumerate(leafNames)}
    nLeaves = len(leafNames)
    slotOf = dict(leafSlot)
    for i in range(nS):
        slotOf[str(S[i])] = nLeaves + i
    base = nLeaves + nS

    tapes = []
    for c, (f_c, g_c, ic_c, f_ss) in enumerate(perCond):
        try:
            op, a, b, cnum, cden, outslots = _emit_tape_shared(f_c, g_c, slotOf, base)
        except _NotRational:
            return {'ok': False}
        fOut = outslots[:nS]
        tape = {
            'op': op, 'a': a, 'b': b, 'cnum': cnum, 'cden': cden,
            'stateSlots': [nLeaves + i for i in range(nS)],
            'fOut': fOut,
            'gOut': outslots[nS:nS + len(g_c)],
            'icLeaf': [-1] * nS, 'icNum': ['0'] * nS, 'icDen': ['1'] * nS,
            # substituted dynamics and observation of this segment, for the
            # multi-condition scaling peel
            'modelLines': ['%s = %s' % (str(S[i]), spy.sympify(f_c[i]))
                           for i in range(nS)],
            'obsLines': ['%s = %s' % (str(obsVars[j]), spy.sympify(g_c[j]))
                         for j in range(len(g_c))],
        }
        if segEq[c] and jointSteadyState:
            # joint/implicit mode: every non-forcing state is a free leaf whose
            # order-0 value R seeds to the on-manifold resting value x* per sample.
            # The IC tape aliases each state slot to its own leaf (identity), so the
            # kernel gives each state an independent dual column. The resting model
            # is kept so R can solve x* and read the df constraint rows [Jx|Jt].
            jfs = set(jointFixedStates or [])
            icMap = {X: (X if (str(X) not in forcings and str(X) not in jfs)
                         else spy.Integer(0)) for X in S}
            # steady state BEFORE the events: the t0 events (a dose at t0) are applied
            # to the state COORDINATE x_ss here, so the observability jet starts at
            # x0 = E(x_ss) while the df constraint is on the pre-event x_ss. R seeds
            # the state leaf with the pre-event value (valBy), and the IC tape carries
            # E with its chain-rule duals w.r.t. the leaves.
            t0evs = (conditionT0Events[c]
                     if conditionT0Events and c < len(conditionT0Events) else [])
            for e in t0evs:
                Xv = spy.Symbol(str(e['var']))
                if Xv not in icMap:
                    continue
                val = spy.sympify(pval(e['value'])).subs(subsMap)
                meth = str(e['method'])
                if meth == 'replace':
                    icMap[Xv] = val
                elif meth == 'add':
                    icMap[Xv] = icMap[Xv] + val
                elif meth == 'multiply':
                    icMap[Xv] = icMap[Xv] * val
            try:
                icOp, icA, icB, icCnum, icCden, icOut = _emit_tape_shared(
                    [icMap[X] for X in S], [], leafSlot, nLeaves)
            except _NotRational:
                return {'ok': False}
            tape.update({'icOp': icOp, 'icA': icA, 'icB': icB, 'icCnum': icCnum,
                         'icCden': icCden, 'icOut': icOut})
            tape['constraintModel'] = [
                '%s = %s' % (invSolveName.get(str(S[i]), str(S[i])),
                             spy.sympify(f_ss[i])) for i in range(nReal)]
        elif segEq[c]:
            # equilibrate-seeded first segment. A generically-linear resting state
            # (no recast) is solved symbolically and emitted as an IC tape, so the
            # kernel seeds each state with its steady-state value and computes the
            # parameter sensitivities by forward-mode duals; dead states fall out as
            # 0 and forcings stay 0. Otherwise the state is seeded numerically per
            # point from the interior steady-state point and its IFT duals (icSeed).
            ssIC = None
            if not powerRecast and not _FORCE_CONSTRAINT_SEED:
                solveS = [S[i] for i in range(nReal) if str(S[i]) not in forcings]
                fssPolys = [spy.expand(spy.fraction(spy.together(spy.sympify(f_ss[i])))[0])
                            for i in range(nReal) if str(S[i]) not in forcings]
                ssSol = _linear_solution(fssPolys, solveS)
                if ssSol is not None:
                    # solved real states take their steady-state value; recast
                    # coordinates E = base^exp and L = log base (states beyond nReal)
                    # stay generic; forcings and dead states stay 0
                    icMap = {X: (ssSol[X] if X in ssSol
                                 else X if i >= nReal else spy.Integer(0))
                             for i, X in enumerate(S)}
                    # t0 events (a dose at t0) compose onto the resting state; their
                    # parameter sensitivities then come from the same forward-mode duals
                    t0evs = (conditionT0Events[c]
                             if conditionT0Events and c < len(conditionT0Events) else [])
                    for e in t0evs:
                        Xv = spy.Symbol(str(e['var']))
                        if Xv not in icMap:
                            continue
                        val = spy.sympify(pval(e['value'])).subs(subsMap)
                        meth = str(e['method'])
                        if meth == 'replace':
                            icMap[Xv] = val
                        elif meth == 'add':
                            icMap[Xv] = icMap[Xv] + val
                        elif meth == 'multiply':
                            icMap[Xv] = icMap[Xv] * val
                    try:
                        icOp, icA, icB, icCnum, icCden, icOut = _emit_tape_shared(
                            [icMap[X] for X in S], [], leafSlot, nLeaves)
                        ssIC = {'icOp': icOp, 'icA': icA, 'icB': icB,
                                'icCnum': icCnum, 'icCden': icCden, 'icOut': icOut}
                    except _NotRational:
                        ssIC = None
            if ssIC is not None:
                tape.update(ssIC)
            else:
                tape['constraintModel'] = [
                    '%s = %s' % (invSolveName.get(str(S[i]), str(S[i])),
                                 spy.sympify(f_ss[i])) for i in range(nReal)]
        else:
            # IC tape: seed each state from its initial-value expression in the
            # leaves (free initial values, carry coordinates, doses, resets).
            try:
                icOp, icA, icB, icCnum, icCden, icOut = _emit_tape_shared(
                    [ic_c[str(X)] for X in S], [], leafSlot, nLeaves)
            except _NotRational:
                return {'ok': False}
            tape.update({'icOp': icOp, 'icA': icA, 'icB': icB,
                         'icCnum': icCnum, 'icCden': icCden, 'icOut': icOut})
        # state-dose event map at this segment's left boundary, applied by the
        # kernel to the propagated state (replace/add/multiply by a parametric
        # value). The values are emitted as a small order-0 tape over the leaves.
        evs = list(conditionEvents[c]) if c < len(conditionEvents) else []
        if evs:
            idxOfState = {str(X): i for i, X in enumerate(S)}
            methodCode = {'replace': 0, 'add': 1, 'multiply': 2}
            keep = [e for e in evs if str(e['var']) in idxOfState]
            if keep:
                vals = [spy.sympify(pval(e['value'])).subs(subsMap) for e in keep]
                try:
                    evOp, evA, evB, evCnum, evCden, evO = _emit_tape_shared(
                        [], vals, leafSlot, nLeaves)
                except _NotRational:
                    return {'ok': False}
                tape.update({
                    'evOp': evOp, 'evA': evA, 'evB': evB, 'evCnum': evCnum,
                    'evCden': evCden, 'evOut': evO,
                    'evVarIdx': [idxOfState[str(e['var'])] for e in keep],
                    'evMethod': [methodCode[str(e['method'])] for e in keep]})
        tapes.append(tape)

    zStateNames = [str(X) for X in freeStates if str(X) not in fixedNames]
    zParamNames = [str(s) for s in params if str(s) not in fixedNames]
    znames = zStateNames + zParamNames
    zSlots = [leafSlot[nm] for nm in znames]

    out = {
        'ok': True,
        'tapes': tapes,
        'nLeaves': nLeaves,
        'nStates': nS,
        'zSlots': zSlots,
        'znames': znames,
        'leafNames': leafNames,
    }
    if equilibrate:
        out['equilibrate'] = True
        out['stateNames'] = [str(X) for X in S]
        out['paramNames'] = [str(s) for s in params]
        out['forcings'] = sorted(forcings)
        out['realStateNames'] = [invSolveName.get(str(X), str(X))
                                 for X in S[:nReal]]
        out['powerRecast'] = powerRecast
        if jointSteadyState:
            out['jointSteadyState'] = True
            # state z-columns (the free states that entered znames), for the R joint
            # branch to seed on-manifold and stack the df constraint rows
            out['zStateNames'] = zStateNames
    return out


###########################################################################
#####################     scaling symmetries     ########################
###########################################################################

def _poly_monomials(expr, zvars):
    """Numerator and denominator monomial-exponent lists of expr over zvars
    (other symbols are treated as weight-zero coefficients). Raises if expr is
    not rational-polynomial in zvars."""
    e = spy.together(spy.sympify(expr))
    p, q = spy.fraction(e)
    pe = spy.Poly(spy.expand(p), *zvars).as_dict()
    qe = spy.Poly(spy.expand(q), *zvars).as_dict()
    return list(pe.keys()), list(qe.keys())


def _scaling_rows(diffEquations, obsFunctions, m, zvars, interOffset):
    """Sparse monomial-exponent rows of one (f, g) system for the scaling kernel:
    weight columns 0..nz-1 are shared over zvars, intermediate columns run from
    interOffset. Returns (rows, ninter, skipped); each row is a {col: coeff} map."""
    nz = len(zvars)
    exprs = []          # (numer monomials, denom monomials, target weight vector)
    skipped = 0
    for g in obsFunctions:
        try:
            pmon, qmon = _poly_monomials(g, zvars)
        except Exception:
            skipped += 1
            continue
        exprs.append((pmon, qmon, [0] * nz))
    for i in range(m):
        try:
            pmon, qmon = _poly_monomials(diffEquations[i], zvars)
        except Exception:
            skipped += 1
            continue
        tvec = [0] * nz
        tvec[i] = 1     # state i is column i of zvars; x_i has weight c_i
        exprs.append((pmon, qmon, tvec))

    rows = []
    for k, (pmon, qmon, tvec) in enumerate(exprs):
        wp, wq = interOffset + 2 * k, interOffset + 2 * k + 1
        for a in pmon:                       # a . c - wp = 0
            row = {j: int(a[j]) for j in range(nz) if a[j]}
            row[wp] = -1
            rows.append(row)
        for b in qmon:                       # b . c - wq = 0
            row = {j: int(b[j]) for j in range(nz) if b[j]}
            row[wq] = -1
            rows.append(row)
        row = {j: -int(tvec[j]) for j in range(nz) if tvec[j]}  # wp - wq - t.c = 0
        row[wp], row[wq] = 1, -1
        rows.append(row)
    return rows, 2 * len(exprs), skipped


def _materialize_rows(rows, ncols):
    """Dense matrix (list of lists) from sparse {col: coeff} rows."""
    out = []
    for r in rows:
        dense = [0] * ncols
        for c, v in r.items():
            dense[c] = v
        out.append(dense)
    return out


def _scaling_gens(rows, ncols, nz):
    """Reduced integer generators (primitive, over the first nz weight columns)
    of the scaling kernel given the determining rows."""
    basis = exactNullspace(rows, ncols) if rows else []
    cvecs = [[v[j] for j in range(nz)] for v in basis
             if any(v[j] != 0 for j in range(nz))]
    gens = []
    if cvecs:
        red, _ = spy.Matrix(cvecs).rref()
        for r in range(red.rows):
            row = [red[r, j] for j in range(nz)]
            if all(x == 0 for x in row):
                continue
            L = 1
            for x in row:
                L = spy.ilcm(L, spy.Rational(x).q)
            ints = [int(spy.Rational(x) * L) for x in row]
            d = 0
            for x in ints:
                d = spy.igcd(d, x)
            if d:
                ints = [x // d for x in ints]
            gens.append(ints)
    return gens


def _scaling_nonid(gens, znames, nz):
    """Scaling entries (support, integer vector, type) from generators."""
    nonId = []
    for c in gens:
        comp = {znames[j]: c[j] for j in range(nz) if c[j] != 0}
        nonId.append({'support': sorted(comp.keys()),
                      'vector': {k: str(v) for k, v in comp.items()},
                      'type': 'scaling'})
    return nonId


def _scaling_gens_recast(gens, znames, nz, recast):
    """Impose the recast relations c_E = exp * c_base on the integer scaling lattice
    span(gens) and return the PHYSICAL scaling generators, whose weights may involve
    the exponent parameter. A free Hill exponent makes c_E = nhill * c_base, which no
    integer weight can satisfy, so the plain integer kernel treats E as free and
    misses (or over-reports) the Hill scaling. Solving span(gens) intersect {relations}
    over Q(exp) recovers it exactly, e.g. xi_kinh = -nhill * kinh, xi_FB = FB -- no
    finite-field sampling, so it is instant at any model size. E = base^exp holds in
    both recast branches, so the same relation applies whether or not `inverted`.
    Returns lists of sympy expressions (possibly symbolic in the exponents)."""
    if not gens:
        return []
    G = spy.Matrix(gens)                      # rows = basis vectors, cols = coords
    idx = {nm: j for j, nm in enumerate(znames)}
    relRows = []
    for rc in recast:
        E, base, exp = str(rc['E']), str(rc['base']), str(rc['exp'])
        if E not in idx or base not in idx:
            continue
        expSym = spy.sympify(exp)
        jE, jB = idx[E], idx[base]
        # (G[:,E] - exp*G[:,base]) . alpha = 0 over the basis coefficients alpha
        relRows.append([G[i, jE] - expSym * G[i, jB] for i in range(G.rows)])
    if not relRows:
        return [[spy.Integer(x) for x in g] for g in gens]
    alphas = spy.Matrix(relRows).nullspace()   # over Q(exp)
    out = []
    for a in alphas:
        v = [spy.together(sum(a[i] * G[i, j] for i in range(G.rows)))
             for j in range(nz)]
        # L = log(base) transforms ADDITIVELY under the scaling: it shifts by c_base
        # (log(lam^{c_base} base) = c_base*log(lam) + L), so its generator weight is
        # c_base. The multiplicative kernel above leaves L at 0; set it here so the
        # generator matches the joint nullspace (which carries the L shift).
        for rc in recast:
            base, L = str(rc['base']), str(rc['L'])
            if base in idx and L in idx:
                v[idx[L]] = v[idx[base]]
        # clear denominators so the weights are polynomial in the exponents
        dens = [spy.denom(x) for x in v if x != 0]
        Lden = spy.Integer(1)
        for d in dens:
            Lden = spy.lcm(Lden, d)
        v = [spy.expand(x * Lden) for x in v]
        if any(x != 0 for x in v):
            out.append(v)
    return out


def scalingSymmetries(allVariables, diffEquations, obsFunctions, m, params,
                      fixed=(), verbose=True):
    """Exact scaling (toric) symmetries z_i -> lam^{c_i} z_i via the integer
    kernel of the monomial-exponent conditions: every observable is invariant
    and every f_i scales with the weight of x_i. Pure linear algebra over the
    integers, so it is exact and scales to large models with no expression
    swell. Returns only the scaling part of the symmetry algebra. Symbols in
    `fixed` are known and may not scale, so their weight is forced to zero."""
    zvars = list(allVariables[:m]) + list(params)
    nz = len(zvars)
    znames = [str(s) for s in zvars]
    fixedset = set(str(s) for s in fixed)

    sparse, ninter, skipped = _scaling_rows(diffEquations, obsFunctions, m, zvars, nz)
    ncols = nz + ninter
    rows = _materialize_rows(sparse, ncols)
    for j in range(nz):                      # a fixed coordinate does not scale
        if znames[j] in fixedset:
            row = [0] * ncols
            row[j] = 1
            rows.append(row)

    nonId = _scaling_nonid(_scaling_gens(rows, ncols, nz), znames, nz)

    if verbose:
        print('-' * 60)
        print('%d scaling symmetry/ies (exact integer kernel)%s:'
              % (len(nonId), '' if not skipped else
                 ', %d non-polynomial term(s) skipped' % skipped))
        for d in nonId:
            print('  ' + ', '.join('%s^(%s)' % (k, d['vector'][k])
                                   for k in d['support']))

    return {
        'method': 'scaling',
        'count': len(nonId),
        'nonIdentifiable': nonId,
    }


def scalingSymmetriesMulti(perCondModel, perCondObs, inputs=None, fixed=None,
                           recast=None):
    """Scaling symmetries common to every condition: the integer kernel of the
    monomial-exponent conditions of all conditions, stacked over a shared weight
    space (each condition contributes its own intermediate columns). The kernel
    of the stack is the intersection of the per-condition scaling lattices, so a
    returned scaling is a symmetry of every condition. Coordinates are the shared
    dynamic states plus the parameters; `inputs` and `fixed` symbols do not scale.
    Returns scalings keyed by coordinate name (states and parameters)."""
    perCondModel = [_as_list(m) for m in perCondModel]
    perCondObs = [_as_list(o) for o in perCondObs]
    inputset = set(_as_list(inputs))
    fixedset = set(_as_list(fixed))
    K = len(perCondModel)
    if K == 0:
        return {'method': 'scaling', 'count': 0, 'nonIdentifiable': []}

    all_lines = [l for lines in perCondModel + perCondObs for l in lines]
    local, parse = _make_local_parse(all_lines)

    stateNames = [_clean(l).split('=', 1)[0].strip() for l in perCondModel[0]]
    stateSyms = [local.get(nm, spy.Symbol(nm)) for nm in stateNames]

    perF, perG, paramset = [], [], set()
    for c in range(K):
        rhs = {}
        for l in perCondModel[c]:
            lhs, expr = _clean(l).split('=', 1)
            rhs[lhs.strip()] = spy.sympify(parse(expr))
        f = [rhs[nm] for nm in stateNames]
        g = [spy.sympify(parse(_clean(l).split('=', 1)[1])) for l in perCondObs[c]]
        perF.append(f)
        perG.append(g)
        for e in f + g:
            paramset |= set(spy.sympify(e).free_symbols)

    paramset -= set(stateSyms)
    paramset -= {local.get(nm, spy.Symbol(nm)) for nm in inputset}
    params = sorted(paramset, key=spy.default_sort_key)
    # re-add eliminated free-exponent bases (see compileObservabilityTapeMulti): a
    # Michaelis constant Km in C^n/(Km^n+C^n) is dropped from the model by the recast
    # but must stay a coordinate so _scaling_gens_recast can impose c_E = exp*c_base and
    # recover the pool/Michaelis co-scaling exactly over Q(exp).
    for rc in _as_list(recast):
        b = local.get(str(rc['base']), spy.Symbol(str(rc['base'])))
        if (b not in set(stateSyms) and b not in params
                and str(b) not in inputset and str(b) not in fixedset):
            params.append(b)
    params = sorted(params, key=spy.default_sort_key)
    zvars = stateSyms + params
    nz = len(zvars)
    znames = [str(s) for s in zvars]
    m = len(stateSyms)

    rows, interOffset, skipped = [], nz, 0
    for c in range(K):
        sparse, ninter, sk = _scaling_rows(perF[c], perG[c], m, zvars, interOffset)
        rows.extend(sparse)
        interOffset += ninter
        skipped += sk
    if skipped:
        print('scalingSymmetriesMulti: %d non-polynomial term(s) skipped '
              '(a scaling they would forbid may be over-reported)' % skipped)
    ncols = interOffset
    dense = _materialize_rows(rows, ncols)
    for j in range(nz):
        if znames[j] in fixedset:
            row = [0] * ncols
            row[j] = 1
            dense.append(row)

    gens = _scaling_gens(dense, ncols, nz)
    # with a power/Hill recast, the integer kernel treats each E = base^exp as a free
    # coordinate; impose c_E = exp*c_base to recover the parameter-weighted (Hill)
    # scalings exactly over Q(exp), instead of leaving them to the rational fit.
    recast = _as_list(recast)
    if recast:
        physical = _scaling_gens_recast(gens, znames, nz, recast)
        nonId = _scaling_nonid(physical, znames, nz)
    else:
        nonId = _scaling_nonid(gens, znames, nz)
    return {'method': 'scaling', 'count': len(nonId), 'nonIdentifiable': nonId}


###########################################################################
#####################  pure-symbolic observability  ######################
###########################################################################

def _lie_deriv(h, states, rhs):
    """Lie derivative L_f h = sum_i (dh/dx_i) f_i along the flow (parameters are
    constants of motion, so they do not contribute)."""
    return sum(spy.diff(h, states[i]) * rhs[i] for i in range(len(states)))


def observabilitySympy(model, observation, fixed=None, parameters=None,
                       inputs=None, equilibrate=False, backend='sympy'):
    """Pure-symbolic observability-identifiability, an INDEPENDENT cross-check of the
    prime-kernel engine for SMALL models. The observability-identifiability matrix
    O = d/dz [g, L_f g, L_f^2 g, ...] over z = (states, unknown parameters) is built
    and rank/nullspace-reduced with sympy over the exact rational-function field --
    NO finite fields and NO power/Hill recast (base^exp stays a symbolic atom, e.g.
    C**nhill is differentiated as nhill*C**(nhill-1)). Because it is exact and
    symbolic the Lie order is carried to the guaranteed saturation bound (# states),
    so there is no premature-truncation risk and no separate saturation guard is
    needed. With equilibrate=True the resting-manifold tangency rows Jf = df/dz are
    appended (the implicit determining system: states stay coordinates, f = 0 enters
    as df.xi = 0). Single condition; events/gaps are not handled here. Returns rank,
    dim, the Lie order used and the non-identifiable directions (a nullspace basis)
    with their support and exact symbolic entries."""
    _select_backend(backend)
    asL = lambda v: list(v) if isinstance(v, (list, tuple)) else ([] if v is None else [v])
    model = [str(l) for l in asL(model)]           # reticulate passes a length-1 vector
    observation = [str(l) for l in asL(observation)]  # as a scalar string; keep it a list
    lines = list(model) + list(observation)
    local, parse = _make_local_parse(lines)
    stateNames = [_clean(l).split('=', 1)[0].strip() for l in model]
    states = [local.get(nm, spy.Symbol(nm)) for nm in stateNames]
    rhs = [spy.sympify(parse(_clean(l).split('=', 1)[1])) for l in model]
    obs = [spy.sympify(parse(_clean(l).split('=', 1)[1])) for l in observation]
    known = {local.get(nm, spy.Symbol(nm)) for nm in (asL(inputs) + asL(fixed))}
    pset = set()
    for e in rhs + obs:
        pset |= spy.sympify(e).free_symbols
    pset -= set(states)
    pset -= known
    params = sorted(pset, key=spy.default_sort_key)
    order = 0
    if equilibrate:
        # resting steady state: solve f = 0 for the states (needs a rational solution;
        # reduceCQ upstream removes conserved moieties so the system is determined) and
        # substitute x*(theta) into the observables. At rest the Lie jet vanishes, so
        # identifiability is over the PARAMETERS through the resting observation g(x*).
        # A non-rational steady state (e.g. a free Hill exponent) is not solvable this
        # way -- that case needs the modular kernel's implicit (tangency) construction.
        try:
            sols = spy.solve([spy.together(fi) for fi in rhs], states, dict=True)
        except Exception:
            sols = []
        if not sols:
            return {'ok': False, 'rank': 0, 'dim': 0, 'lieOrder': 0,
                    'nonIdentifiable': [], 'identifiable': False,
                    'why': 'equilibrate: the steady state is not symbolically solvable '
                           '(e.g. a free Hill exponent); use symEngine="modular"'}
        Gobs = [spy.together(h.subs(sols[0])) for h in obs]
        z = list(params)
        nz = len(z)
        znames = [str(s) for s in z]
        if nz == 0:
            return {'ok': True, 'rank': 0, 'dim': 0, 'lieOrder': 0,
                    'nonIdentifiable': [], 'identifiable': True}
        rows = [[spy.diff(h, zj) for zj in z] for h in Gobs]
    else:
        z = states + params
        nz = len(z)
        znames = [str(s) for s in z]
        if nz == 0:
            return {'ok': True, 'rank': 0, 'dim': 0, 'lieOrder': 0,
                    'nonIdentifiable': [], 'identifiable': True}
        rows = []
        jet = list(obs)
        prev = -1
        flat = 0
        # the extended (states + constant parameters) system has dimension nz, so its
        # observability codistribution is spanned by the Lie derivatives up to order
        # nz - 1; carry the jet to that guaranteed bound (early exit once the rank holds
        # for two consecutive orders) -- exact, so no premature-truncation risk.
        while True:
            for h in jet:
                rows.append([spy.diff(h, zj) for zj in z])
            rank = spy.Matrix(rows).rank()
            flat = flat + 1 if rank == prev else 0
            if rank >= nz or (flat >= 2 and order >= 1) or order >= nz:
                break
            prev = rank
            order += 1
            jet = [spy.expand(_lie_deriv(h, states, rhs)) for h in jet]
    M = spy.Matrix(rows)
    rank = M.rank()
    nonId = []
    if rank < nz:
        for vec in M.nullspace():
            v = [spy.cancel(vec[i]) for i in range(nz)]
            nz_i = [i for i in range(nz) if v[i] != 0]
            nonId.append({'support': sorted(znames[i] for i in nz_i),
                          'vector': {znames[i]: str(v[i]) for i in nz_i},
                          'type': 'general', 'closedForm': True})
    return {'ok': True, 'rank': int(rank), 'dim': int(nz), 'lieOrder': int(order),
            'nonIdentifiable': nonId, 'identifiable': bool(rank >= nz)}


###########################################################################
#####################     R entry point     ##############################
###########################################################################

def symmetryDetectiondMod(model, observation,
                          ansatz='uni', pMax=2, inputs=None, fixed=None,
                          allTrafos=False, lieOrder=0, exact=True,
                          verify=True, backend='sympy', parameters=None,
                          method='liesym'):
    global _EXACT
    _EXACT = bool(exact)
    _select_backend(backend)

    model = _as_list(model)
    observation = _as_list(observation)
    inputNames = _as_list(inputs)
    fixedNames = _as_list(fixed)
    extraParams = _as_list(parameters)

    all_lines = model + observation + inputNames + fixedNames + extraParams
    local, parse = _make_local_parse(all_lines)

    sys.stdout.write('\nReading input...')
    sys.stdout.flush()

    variables, diffEquations, params = _read_equations(model, parse)

    _, obsFunctions, params = _read_observation(
        observation, variables, params, parse)

    # explicit parameters
    for pn in extraParams:
        s = local.get(pn, spy.Symbol(pn))
        if s not in params and s not in variables:
            params.append(s)

    inputSyms = [local.get(nm, spy.Symbol(nm)) for nm in inputNames]
    fixedSyms = [local.get(nm, spy.Symbol(nm)) for nm in fixedNames]
    for s in inputSyms:
        if s in params:
            params.remove(s)

    # a declared input that is also a state (a known switch) is a constant, not
    # a dynamic variable: drop its equation so it appears once, as an input
    keep = [i for i, v in enumerate(variables) if v not in inputSyms]
    variables = [variables[i] for i in keep]
    diffEquations = [diffEquations[i] for i in keep]

    allVariables = list(variables) + list(inputSyms) + list(params)

    sys.stdout.write('done\n')
    sys.stdout.flush()

    if str(method) == 'scaling':
        return scalingSymmetries(
            allVariables, diffEquations, obsFunctions, len(variables), params,
            fixed=fixedSyms)

    return symmetryDetection(
        allVariables, diffEquations, obsFunctions, ansatz=ansatz, pMax=int(pMax),
        inputs=inputSyms, fixed=fixedSyms, lieOrder=int(lieOrder),
        allTrafos=bool(allTrafos), verify=bool(verify))


def _read_observation(observation, variables, parameters, parse):
    obsVars, obsFunctions, obsParameters = _read_equations(observation, parse)
    for var in variables:
        if var in obsParameters:
            obsParameters.remove(var)
    for par in parameters:
        if par in obsParameters:
            obsParameters.remove(par)
    return obsVars, obsFunctions, parameters + obsParameters
