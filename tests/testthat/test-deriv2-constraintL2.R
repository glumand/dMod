test_that("constraintL2 deriv2 adds gi . dP2 chain term after Pexpl", {
  prev <- getwd(); on.exit(setwd(prev), add = TRUE)
  setwd(tempdir())

  # Pexpl: a = exp(la). Then constraintL2(mu = mu_a, sigma = s) on `a`,
  # composed via attr(p, "deriv") and attr(p, "deriv2") from Pexpl.
  pfn <- Pexpl(c(a = "exp(la)"), parameters = NULL,
               modelname = paste0("c2_pexpl_", as.integer(Sys.time())),
               compile = TRUE, deriv2 = TRUE, derivMode = "dual")

  mu <- c(a = 1.0); sg <- 0.5
  cfn <- constraintL2(mu = mu, sigma = sg)

  pars <- c(la = 0.3)
  pinner <- pfn(pars, deriv = TRUE, deriv2 = TRUE)[[1]]

  # The constraint expects `a` in its parameter set. Provide the inner value.
  res_gn <- cfn(pinner, deriv = TRUE, deriv2 = FALSE)
  res_ex <- cfn(pinner, deriv = TRUE, deriv2 = TRUE)

  # Reference: L = ((exp(la) - mu) / sg)^2.
  # dL/dla = 2 * (exp(la) - mu) * exp(la) / sg^2
  # d^2L/dla^2 = 2 * exp(la)^2 / sg^2 + 2 * (exp(la) - mu) * exp(la) / sg^2
  la <- pars["la"]
  ref_grad <- unname(2 * (exp(la) - mu) * exp(la) / sg^2)
  ref_hess_GN <- unname(2 * exp(la)^2 / sg^2)
  ref_hess_exact <- ref_hess_GN + unname(2 * (exp(la) - mu) * exp(la) / sg^2)

  expect_equal(unname(res_gn$value), unname(((exp(la) - mu) / sg)^2),
               tolerance = 1e-10)
  expect_equal(unname(res_gn$gradient["la"]), ref_grad, tolerance = 1e-10)
  expect_equal(unname(res_gn$hessian["la", "la"]), ref_hess_GN,
               tolerance = 1e-10)
  expect_equal(unname(res_ex$hessian["la", "la"]), ref_hess_exact,
               tolerance = 1e-10)
})
