# Behavioral tests for the explicit parameter transformation Pexpl() / P().
#
# Verifies:
#   * identity trafo round-trips (value + Jacobian)
#   * log trafo gives Jacobian = diag(exp(theta)) = diag(p)
#   * a mixed nonlinear trafo's Jacobian matches the algebraic derivative
#   * derivMode "symbolic" and "dual" agree on value and Jacobian
#   * getParameters() consistency through composition (Y * Xs * P)
#
# Second-order chain rule is covered by test-deriv2-Pexpl.R.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
}


## ---- Identity transformation -------------------------------------------

test_that("Pexpl identity trafo round-trips and has identity Jacobian", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  outer <- c(A = 1.5, k = 0.3)
  inner <- bench$pfn_id(outer, deriv = TRUE)
  # parlist[[C1]] is a parvec carrying value + "deriv" attr.
  pv <- inner$C1
  expect_equal(as.numeric(pv), as.numeric(outer))
  J  <- attr(pv, "deriv")
  expect_equal(unname(J), diag(2))
  expect_setequal(rownames(J), c("A", "k"))
  expect_setequal(colnames(J), c("A", "k"))
})


## ---- Log transformation -----------------------------------------------

test_that("Pexpl log trafo maps theta -> exp(theta) with Jacobian diag(exp(theta))", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  outer <- c(A_log = log(1.5), k_log = log(0.3))
  inner <- bench$pfn_log(outer, deriv = TRUE)
  pv <- inner$C1
  expect_equal(as.numeric(pv), c(1.5, 0.3))

  # J[i, j] = d (inner_i) / d (outer_j). Diagonal entries are exp(theta) = p.
  J <- attr(pv, "deriv")
  expect_equal(unname(diag(J)), c(1.5, 0.3))
  expect_equal(J[upper.tri(J)], rep(0, sum(upper.tri(J))))
  expect_equal(J[lower.tri(J)], rep(0, sum(lower.tri(J))))
})


## ---- Mixed nonlinear trafo: analytical Jacobian -----------------------

test_that("Pexpl Jacobian on a mixed nonlinear trafo equals the algebraic derivative", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # Mixed trafo: A = a^2, k = a * b
  # Analytical Jacobian: J[1,] = (2a, 0), J[2,] = (b, a).
  pfn <- P(eqnvec(A = "a^2", k = "a * b"),
           condition = "C1",
           modelname = "test_P_mix", compile = TRUE)

  outer <- c(a = 1.3, b = 0.7)
  inner <- pfn(outer, deriv = TRUE)$C1
  J <- attr(inner, "deriv")

  J_ref <- rbind(
    A = c(a = 2 * outer[["a"]], b = 0),
    k = c(a = outer[["b"]],     b = outer[["a"]]))
  expect_equal(unname(J), unname(J_ref), tolerance = 1e-8)
})


## ---- derivMode parity --------------------------------------------------

test_that("Pexpl derivMode 'symbolic' and 'dual' agree on value and Jacobian", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  pfn_sym <- P(eqnvec(A = "exp(a)", k = "exp(b)"), condition = "C1",
               method = "explicit", derivMode = "symbolic",
               modelname = "test_P_dm_sym", compile = TRUE)
  pfn_dual <- P(eqnvec(A = "exp(a)", k = "exp(b)"), condition = "C1",
                method = "explicit", derivMode = "dual",
                modelname = "test_P_dm_dual", compile = TRUE)

  outer <- c(a = 0.3, b = -0.5)
  i_sym  <- pfn_sym (outer, deriv = TRUE)$C1
  i_dual <- pfn_dual(outer, deriv = TRUE)$C1
  expect_equal(as.numeric(i_sym), as.numeric(i_dual), tolerance = 1e-10)
  expect_equal(attr(i_sym, "deriv"), attr(i_dual, "deriv"), tolerance = 1e-10)
})


## ---- Parameter set propagation through composition --------------------

test_that("getParameters(Y * Xs * P) equals getParameters(P) (outer-pars view)", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  expect_setequal(getParameters(bench$prd_id),  getParameters(bench$pfn_id))
  expect_setequal(getParameters(bench$prd_log), getParameters(bench$pfn_log))
})
