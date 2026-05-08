context("omega() helper for FOCEI")


test_that("diagonal omega returns K Cholesky parameters", {
  om <- omega(eta = c("eta_Cl", "eta_V"))
  expect_equal(om$K, 2L)
  expect_equal(om$structure, "diag")
  expect_setequal(om$chol_pars, c("omega_Cl_Cl", "omega_V_V"))
  expect_true(all(om$is_diag))
  expect_null(om$subject_etas)
})


test_that("full omega returns K(K+1)/2 Cholesky parameters", {
  om <- omega(eta = c("eta_Cl", "eta_V", "eta_Ka"), structure = "full")
  expect_equal(om$K, 3L)
  expect_equal(length(om$chol_pars), 6L)
  expect_setequal(om$chol_pars,
                  c("omega_Cl_Cl", "omega_V_V", "omega_Ka_Ka",
                    "omega_V_Cl", "omega_Ka_Cl", "omega_Ka_V"))
  expect_equal(sum(om$is_diag), 3L)
})


test_that("selective correlation adds exactly one off-diagonal", {
  om <- omega(eta = c("eta_Cl", "eta_V", "eta_Ka"),
              correlate = list(c("eta_Cl", "eta_V")))
  expect_equal(length(om$chol_pars), 4L)
  expect_true("omega_V_Cl" %in% om$chol_pars)
  expect_false("omega_Ka_Cl" %in% om$chol_pars)
})


test_that("subject expansion produces N x K name matrix", {
  subjects <- paste0("subj", 1:4)
  om <- omega(eta = c("eta_Cl", "eta_V"), subjects = subjects)
  expect_equal(dim(om$subject_etas), c(4L, 2L))
  expect_equal(om$subject_etas[1, "eta_Cl"], "eta_Cl_subj1")
  expect_equal(om$subject_etas[3, "eta_V"],  "eta_V_subj3")
})


test_that("build_L returns lower-triangular matrix with exp(diag)", {
  om <- omega(eta = c("eta_Cl", "eta_V"), structure = "full")
  chol_vec <- c(omega_Cl_Cl = log(0.3), omega_V_V = log(0.2),
                omega_V_Cl  = 0.05)
  L <- om$build_L(chol_vec)
  expect_equal(dim(L), c(2L, 2L))
  expect_equal(L[1, 1], 0.3)
  expect_equal(L[2, 2], 0.2)
  expect_equal(L[2, 1], 0.05)
  expect_equal(L[1, 2], 0)
})


test_that("Omega = L L^T is symmetric positive definite", {
  set.seed(1)
  om <- omega(eta = c("eta_a", "eta_b", "eta_c"), structure = "full")
  chol_vec <- c(omega_a_a = log(0.3), omega_b_b = log(0.4), omega_c_c = log(0.2),
                omega_b_a = 0.1, omega_c_a = -0.05, omega_c_b = 0.07)
  L <- om$build_L(chol_vec)
  Omega_mat <- L %*% t(L)
  expect_equal(Omega_mat, t(Omega_mat))
  ev <- eigen(Omega_mat, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > 0))
})


test_that("missing Cholesky parameter raises an error", {
  om <- omega(eta = c("eta_Cl", "eta_V"))
  expect_error(om$build_L(c(omega_Cl_Cl = 0)), "missing Cholesky")
})


test_that("invalid eta arguments are rejected", {
  expect_error(omega(eta = character(0)), "non-empty")
  expect_error(omega(eta = c("eta_a", "eta_a")), "unique")
})


test_that("correlate with unknown eta name fails clearly", {
  expect_error(omega(eta = c("eta_Cl"), correlate = list(c("eta_Cl", "eta_X"))),
               "unknown eta names")
})


test_that("parnames.omegaSpec returns correct subsets", {
  om <- omega(eta = c("eta_Cl", "eta_V"), subjects = c("s1", "s2"))
  all_names  <- parnames.omegaSpec(om, "all")
  eta_names  <- parnames.omegaSpec(om, "eta")
  chol_names <- parnames.omegaSpec(om, "chol")
  expect_equal(length(eta_names),  4L)
  expect_equal(length(chol_names), 2L)
  expect_equal(length(all_names),  6L)
  expect_setequal(all_names, c(eta_names, chol_names))
})
