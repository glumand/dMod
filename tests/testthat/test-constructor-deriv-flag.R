# Constructor-level `deriv = TRUE/FALSE` gating for P, Pexpl, Pimpl,
# Pequil, and Y. Symmetric with the existing `deriv2` flag: the
# constructor decides whether the artifact carries first-order
# sensitivities, and the runtime call errors out if it asks for
# something the construction didn't produce.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
}


test_that("Pexpl(deriv = FALSE) yields a parvec without deriv attribute", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  trafo <- c(A = "a * x", B = "b + y")
  pf <- Pexpl(trafo, deriv = FALSE,
              modelname = paste0("test_pexpl_nod1_", as.integer(Sys.time())),
              compile = TRUE, derivMode = "symbolic", verbose = FALSE)
  out <- pf(c(a = 2, b = 3, x = 4, y = 5))
  expect_null(attr(out[[1]], "deriv"))
  # Default runtime deriv = TRUE is silently capped by the constructor:
  # no error, just no deriv attribute on the result.
  out2 <- pf(c(a = 2, b = 3, x = 4, y = 5), deriv = TRUE)
  expect_null(attr(out2[[1]], "deriv"))
})


test_that("Pequil(deriv = FALSE) skips the sensitivity model", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  pf <- Pequil(c(A = "k_in - k_out * A"),
               parameters = c("k_in", "k_out"),
               deriv = FALSE,
               modelname = paste0("test_pequil_nod1_", as.integer(Sys.time())),
               compile = TRUE, verbose = FALSE, attach.input = FALSE)
  out <- pf(c(k_in = 1, k_out = 0.5, A = 0.1))
  expect_null(attr(out[[1]], "deriv"))
})


test_that("Pimpl(deriv = FALSE) drops the IFT chain rule from output", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  pf <- Pimpl(c(x = "x - a"), parameters = "a", deriv = FALSE,
              modelname = paste0("test_pimpl_nod1_", as.integer(Sys.time())),
              compile = TRUE, verbose = FALSE)
  out <- pf(c(a = 1.5, x = 0.5))
  expect_null(attr(out[[1]], "deriv"))
})


test_that("Y(deriv = FALSE) produces output without deriv attribute", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  g <- c(obs = "k * A")
  gfn <- Y(g, states = c("A", "time"), parameters = "k", deriv = FALSE,
           modelname = paste0("test_y_nod1_", as.integer(Sys.time())),
           compile = TRUE, derivMode = "symbolic", verbose = FALSE,
           attach.input = FALSE)

  prd <- structure(
    cbind(time = c(0, 1), A = c(1, 2)),
    parameters = structure(c(k = 1.5), fixed = NULL),
    class = c("prdframe", "matrix", "array"))

  res <- gfn(out = prd, pars = c(k = 1.5))
  expect_null(attr(res[[1]], "deriv"))
})


test_that("Constructors reject deriv = FALSE combined with deriv2 = TRUE", {
  expect_error(Pexpl(c(A = "x"), deriv = FALSE, deriv2 = TRUE, compile = FALSE),
               "requires deriv = TRUE")
  expect_error(Pimpl(c(x = "x - a"), parameters = "a", deriv = FALSE, deriv2 = TRUE),
               "requires deriv = TRUE")
  expect_error(Pequil(c(A = "k - A"), parameters = "k", deriv = FALSE, deriv2 = TRUE),
               "requires deriv = TRUE")
  expect_error(Y(c(obs = "A"), states = c("A", "time"),
                 deriv = FALSE, deriv2 = TRUE),
               "requires deriv = TRUE")
  expect_error(P(c(A = "k"), method = "explicit",
                 deriv = FALSE, deriv2 = TRUE),
               "requires deriv = TRUE")
})


test_that("P() dispatcher forwards deriv to each method", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  pf <- P(c(A = "a * x"), method = "explicit", deriv = FALSE,
          modelname = paste0("test_P_nod1_", as.integer(Sys.time())),
          compile = TRUE, verbose = FALSE)
  out <- pf(c(a = 2, x = 3))
  expect_null(attr(out[[1]], "deriv"))
})
