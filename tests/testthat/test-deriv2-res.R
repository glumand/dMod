test_that("res() reduces 4D deriv2 attribute to 3D [n_residuals, p, p]", {
  # Build a synthetic prediction matrix and 4D deriv2 attribute, run through
  # res(), and verify the Hessian-of-residuals indexing is correct.

  times <- c(0.0, 0.5, 1.0)
  obs   <- c("y1", "y2")
  pars  <- c("p1", "p2")

  # value matrix [time, name]
  out <- cbind(time = times,
               y1 = c(1.0, 2.0, 3.0),
               y2 = c(0.1, 0.2, 0.3))
  d2 <- array(NA_real_, c(length(times), length(obs), length(pars), length(pars)),
              dimnames = list(NULL, obs, pars, pars))
  # Fill with a unique encoding so we can verify indexing.
  for (i in seq_along(times))
    for (j in seq_along(obs))
      for (k1 in seq_along(pars))
        for (k2 in seq_along(pars))
          d2[i, j, k1, k2] <- 1000 * i + 100 * j + 10 * k1 + k2
  attr(out, "deriv2") <- d2
  attr(out, "deriv") <- array(0, c(length(times), length(obs), length(pars)),
                              dimnames = list(NULL, obs, pars))

  data_df <- data.frame(name = c("y1", "y2", "y1"),
                        time = c(0.0, 0.5, 1.0),
                        value = c(1.05, 0.18, 2.95),
                        sigma = c(0.1, 0.1, 0.1),
                        lloq = -Inf)

  ro <- res(data_df, out)
  rd2 <- attr(ro, "deriv2")
  expect_equal(dim(rd2), c(3, 2, 2))
  # Row 1: y1 at t=0.0 -> i=1, j=1 -> entries 1100 + 10*k1 + k2.
  expect_equal(rd2[1, , ], matrix(c(1111, 1121, 1112, 1122), 2, 2,
                                  dimnames = list(pars, pars)))
  # Row 2: y2 at t=0.5 -> i=2, j=2 -> 2200 + 10*k1 + k2.
  expect_equal(rd2[2, , ], matrix(c(2211, 2221, 2212, 2222), 2, 2,
                                  dimnames = list(pars, pars)))
  # Row 3: y1 at t=1.0 -> i=3, j=1 -> 3100 + 10*k1 + k2.
  expect_equal(rd2[3, , ], matrix(c(3111, 3121, 3112, 3122), 2, 2,
                                  dimnames = list(pars, pars)))
})
