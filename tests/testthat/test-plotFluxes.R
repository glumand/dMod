context("plotFluxes")
test_that("plotFluxes runs end-to-end on a small reaction network", {

  withr::local_dir(tempdir())

  f <- NULL |>
    addReaction("A", "B", "k1*A", "production") |>
    addReaction("B",  "", "k2*B", "degradation")

  model <- odemodel(f, modelname = "pf_model")
  x <- Xs(model)

  pars <- c(A = 1, B = 0, k1 = 0.5, k2 = 0.3)
  times <- seq(0, 10, length.out = 21)

  fluxEquations <- getFluxes(f)$B

  P <- plotFluxes(pars, x, times, fluxEquations, nameFlux = "B fluxes")

  # plotFluxes returns a ggplot
  expect_s3_class(P, "ggplot")

  # The long-format flux table is attached as attribute "out" and should
  # have a row per (time, condition, flux name).
  out <- attr(P, "out")
  expect_true(is.data.frame(out))
  expect_true(all(c("time", "condition", "name", "value") %in% colnames(out)))
  expect_equal(nrow(out), length(times) * length(fluxEquations))

  # Flux values should be finite (the underlying CppODE::funCpp evaluation
  # is what we're really exercising here — the previous cOde::funC0 backend
  # is gone, so this guards against regressions in the new dispatch).
  expect_true(all(is.finite(out$value)))
})
