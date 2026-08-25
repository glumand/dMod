## ---------------------------------------------------------------------------
## FIM null-space group regularization -- reference implementation.
##
## Implements the penalty derived in notes/group_regularization_null_space.tex
##
## The idea (see the note for the full derivation): near a reference fit p*,
## the Gauss-Newton Hessian H = S'WS splits parameter space into a "stiff"
## subspace the data pins down and a "sloppy"/null subspace it says nothing
## about (H's near-zero eigenvalues). Instead of an ordinary ridge penalty
## R = lambda*||p-p*||^2, which also pulls on the identified (stiff)
## direction, project the displacement onto the null subspace FIRST and only
## penalize that:
##
##   P      = V_null %*% t(V_null)                    (n x n projector, symmetric, idempotent)
##   R(p)   = lambda * (p-p*)' P (p-p*) = lambda*||V_null' (p-p*)||^2
##   grad R = 2*lambda * P (p-p*)
##   hess R = 2*lambda * P                             (constant, exact -- R is exactly quadratic)
##
## R is IDENTICALLY zero (not just small) for any displacement confined to
## the stiff subspace (note, Section 4.4) -- unlike a diagonal
## constraintL2() prior, which cannot null out a direction that is a mix of
## several raw parameters.
##
#---------------------------------------------------------------------------

## nullSpaceBasis(): split a symmetric Hessian's eigenbasis into "stiff"
## (identified) and "null" (structurally non-identifiable) subspaces by
## thresholding eigenvalues relative to the largest one:
##   is_zero_i = 1[ |s_i| < tol * max_j|s_j| ]
##
## H should be a named (square, symmetric) parameter-by-parameter Hessian,
## e.g. obj(pars)$hessian at a reference fit -- typically restricted to the
## outer/fitted parameters first, as in `H <- o_best$hessian[nm, nm]`.
##
## Returns a list:
##   values  -- eigenvalues, decreasing (as returned by eigen())
##   vectors -- full eigenvector matrix V (n x n), columns = eigenvectors
##   is_zero -- logical, TRUE for columns classified as null/sloppy
##   V_null  -- eigenvector columns with is_zero == TRUE  (n x m, m = #null directions)
##   V_stiff -- eigenvector columns with is_zero == FALSE (n x (n-m))
##

nullSpaceBasis <- function(H, tol = 1e-6) {
  stopifnot(is.matrix(H), nrow(H) == ncol(H))
  nm <- rownames(H)

  eig <- eigen(H, symmetric = TRUE)
  s <- eig$values
  is_zero <- abs(s) < tol * max(abs(s)) # scale tolerance on maximum eigenvalue

  V_null  <- eig$vectors[, is_zero,  drop = FALSE]  # drop = false keeps dimension of matrix
  V_stiff <- eig$vectors[, !is_zero, drop = FALSE]
  if (!is.null(nm)) {
    rownames(V_null)  <- nm
    rownames(V_stiff) <- nm
  }

  list(values = s, vectors = eig$vectors, is_zero = is_zero,
       V_null = V_null, V_stiff = V_stiff)
}

## nullSpaceProjector(): the n x n projector P = V_null %*% t(V_null) onto
## the null subspace found by nullSpaceBasis(). P is
## symmetric and idempotent (P %*% P == P) by construction -- it depends
## only on the null SUBSPACE, never on eigen()'s arbitrary choice of basis
## within it when eigenvalues repeat 
##
## If H carries no dimnames, P comes back without dimnames too -- groupL2()
## below still works as long as p0 is named and matches H's parameter order.
nullSpaceProjector <- function(H, tol = 1e-6) {
  basis <- nullSpaceBasis(H, tol = tol)
  P <- basis$V_null %*% t(basis$V_null)
  if (!is.null(rownames(H))) dimnames(P) <- list(rownames(H), rownames(H))
  P
}

## groupL2(): the regularization itself. Builds a dMod objfn for
##   R(p) = lambda * (p-p0)' P (p-p0)
## centered on a reference point p0 (typically the fit H was computed at --
## note, Section 2, "Why expand around a point at all?"), penalizing only
## the component of p-p0 that falls inside P's span. Composes with any other
## dMod objective via "+" ("+.fn"/"+.objlist" sum value/gradient/hessian by
## parameter NAME, so obj_group <- obj + groupL2(P, best) is all that's
## needed, same pattern as quadraticPenalty()/constraintL2()).
##
## Follows quadraticPenalty()'s calling convention (established in
## scripts/example_toy_mapk_fim_fit.R and vetted there against normL2's
## actual objfn call signature) rather than the note's own simplified
## `list(...)[[1]]` sketch -- pulls the parameter vector out via
## match.fnargs(arglist, "pars") so this also works when callers pass
## `pars = ` by name rather than positionally.
##
## Value/gradient/Hessian are exact and closed-form (R is exactly quadratic
## in p, note Section 4.1) -- no deriv/deriv2 machinery needed, unlike the
## data objective it gets added to.
##
## Arguments:
##   P         n x n projector (from nullSpaceProjector(), or any symmetric
##             idempotent matrix -- e.g. hand-built as in the note's
##             three-parameter example). If P carries dimnames, they are
##             matched against names(p0) BY NAME and P is reordered to
##             p0's order if needed -- deliberately not positional (see the
##             trust()/trust_kernel.cpp align_grad/align_hess fix: silently
##             trusting matrix/vector position over names is exactly the
##             class of bug that fix exists to avoid).
##   p0        Named numeric vector, the reference point p* (e.g. `best`
##             from a trust()/mstrust() fit), at the SAME (outer/trafo)
##             parameter level P/H were computed at.
##   lambda    Overall penalty strength (scalar, default 0.1 -- same default
##             as the note's code listing).
##   condition dMod condition this penalty applies to; NULL (default) makes
##             it condition-unspecific, so it composes with an objective
##             defined over any/all conditions (same convention as
##             quadraticPenalty()).
groupL2 <- function(P, p0, lambda = 0.1, condition = NULL) {
  nm <- names(p0)
  if (is.null(nm)) stop("groupL2: p0 must be a named parameter vector.")

  if (!is.null(dimnames(P))) {
    stopifnot(setequal(rownames(P), nm), setequal(colnames(P), nm))
    P <- P[nm, nm, drop = FALSE]   # align to p0's name order, not position
  } else {
    stopifnot(nrow(P) == length(nm), ncol(P) == length(nm))
    dimnames(P) <- list(nm, nm)
  }
  force(P); force(p0); force(lambda)  # evaluate arguments immediately, to avoid overwriting issues

  myfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                    conditions = condition, env = NULL) {
    arglist <- list(...)
    pars <- arglist[[match.fnargs(arglist, "pars")]]
    p <- pars[nm]
    delta <- p - p0   # deviation from the p0 we regularize to
    Pd <- as.numeric(P %*% delta)   # apply projection on deviation to only regularize uninformed directions

    value    <- lambda * sum(delta * Pd)          # R(p)      = lambda * delta' P delta
    gradient <- setNames(2 * lambda * Pd, nm)      # grad R    = 2*lambda * P delta
    hessian  <- 2 * lambda * P                     # hessian R = 2*lambda * P  (constant)
    dimnames(hessian) <- list(nm, nm)

    out <- objlist(value = value, gradient = gradient, hessian = hessian)
    attr(out, "prior") <- value
    attr(out, "env") <- env
    out
  }
  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- nm
  myfn
}
