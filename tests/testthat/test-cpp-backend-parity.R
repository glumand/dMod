# Cross-backend parity smoke for the C++ objective kernels.
#
# Behavioral correctness of both backends is verified in the test-normL2-*
# and test-constraintL2-* files via for_each_backend(). This file adds a
# minimal cross-backend parity check on a fixed input so a future kernel
# refactor that silently breaks the wiring is caught even if the
# behavioural tests pass on one backend only.
#
# One block per objective kernel pair (normL2, constraintL2). The
# residual kernel (src/residual_kernel.{h,cpp}) retains its own dedicated
# behavioural-plus-parity file (test-residual-kernel.R) because that file
# also covers BLOQ-deriv2-exact via finite-difference Hessian validation
# which has no equivalent at the R / dMod-API layer.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("CppODE")
  testthat::skip_on_cran()
}


test_that("normL2 cpp kernel agrees with R reference on the linear-decay fixture", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.05)
  pars  <- bench$outerpars_id

  with_cpp_backend(FALSE, {
    o_R <- normL2(data, bench$prd_id)(pars)
  })
  with_cpp_backend(TRUE, {
    o_C <- normL2(data, bench$prd_id)(pars)
  })
  expect_equal(o_C$value,    o_R$value,    tolerance = 1e-9)
  expect_equal(o_C$gradient, o_R$gradient, tolerance = 1e-8)
  expect_equal(o_C$hessian,  o_R$hessian,  tolerance = 1e-8)
})


test_that("constraintL2 cpp kernel agrees with R reference on a small diagonal case", {
  obj <- constraintL2(mu = c(a = 0.0, b = 1.0, c = -0.5), sigma = c(0.5, 1, 2))
  pars <- c(a = 0.3, b = 1.4, c = -0.6)
  with_cpp_backend(FALSE, { o_R <- obj(pars) })
  with_cpp_backend(TRUE,  { o_C <- obj(pars) })
  expect_equal(o_C$value,    o_R$value,    tolerance = 1e-12)
  expect_equal(o_C$gradient, o_R$gradient, tolerance = 1e-12)
  expect_equal(o_C$hessian,  o_R$hessian,  tolerance = 1e-12)
})
