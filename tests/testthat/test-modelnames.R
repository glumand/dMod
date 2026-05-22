context("Conflicting modelnames")
test_that("modelnames behave as expected", {
  
  # What needs to be checked
  # 1. That modelname is what goes in (merged-SO: P() with multi-condition trafo
  #    links every per-condition .cpp into one shared library named after the
  #    `modelname` arg, and modelname() reports that library name).
  # 2. That recompiling the same modelname is rebuild-safe: previously-built
  #    objects keep working because their symbols re-resolve against the
  #    (re-)loaded SO.
  # 3. That compiling the same structural model into a different modelname lets both functions intact
  
  #-!Start example code
  #-! library(conveniencefunctions)
  #-! library(dMod)
  library(dplyr)
  # Run codegen in a temp dir so the .c/.so artefacts don't leak into
  # tests/testthat/. withr::local_dir auto-restores at test_that exit so
  # later tests don't load stale shared libs from this test's tempdir.
  withr::local_dir(tempdir())

  ## Model definition (text-based, scripting part)
  f <- NULL %>%
    addReaction("A", "B", "k1*A", "translation") %>%
    addReaction("B",  "", "k2*B", "degradation") %>%
    as.eqnvec()
  events <- eventlist(var = "A", time = 5, value = "A_add", method = "add")

  # Use a name distinct from the default `odemodel` used in test-odemodel.R —
  # both tests would otherwise share a DLL symbol space across the testthat
  # session and corrupt each other's parms-length contract.
  x1 <- odemodel(f, events = events, modelname = "mn_odemodel", compile = FALSE) %>% Xs
  compile(x1)
  g1 <- Y(c(Bobs = "s1*B"), x1, compile = T, modelname = "obsfn")
  
  conditions <- c("a", "b")
  trafo <-
    getParameters(g1,x1) %>%
    setNames(.,.) %>%
    branch(conditions = conditions) %>%
    insert("x~x_cond", x = "s1", cond = condition) %>%
    insert("x~exp(x)", x = getSymbols(mytrafo[[i]])) %>%
    {.}
  
  p1 <- P(trafo, modelname = "p", compile = T)
  
  parameters <- getParameters(p1)
  pars <- structure(rnorm(length(parameters)), names = parameters)
  (g1*x1*p1)(0:10, pars)
  
  #-!End example code
  # 2. Rerunning the same parts breaks existing objects
  g2 <- Y(c(Bobs = "s1*B"), x1, compile = T, modelname = "obsfn")
  p2 <- P(trafo, modelname = "p", compile = T)
  
  # 3. Compiling the same structural model into a different modelname lets both functions intact
  g3 <- Y(c(Bobs = "s1*B"), x1, compile = T, modelname = "obsfn3")
  x3 <- odemodel(f, events = events, modelname = "mn_odemodel3", compile = FALSE) %>% Xs
  compile(x3)

  # Define your expectations here
  # 1. Modelname is what goes in. P() with compile=TRUE links the per-condition
  #    sources into a single SO named `modelname`, so modelname(p1) is "p"
  #    rather than the per-condition codegen names ("p_a", "p_b").
  expect_equal(modelname(x1), "mn_odemodel")
  expect_equal(modelname(g1), "obsfn")
  expect_equal(modelname(p1), "p")
  # 2. Recompiling the same modelname is rebuild-safe: every combination still
  #    evaluates because compile() re-loads the merged SO and old objects
  #    re-resolve their symbols against it.
  expect_true(!inherits(try((g2*x1*p2)(0:10,pars)), "try-error"))
  expect_true(!inherits(try((g1*x1*p2)(0:10,pars)), "try-error"))
  expect_true(!inherits(try((g2*x1*p1)(0:10,pars)), "try-error"))
  # 3. Compiling the same structural model into a different modelname lets both functions intact
  expect_true(!inherits(try((g2*x3*p2)(0:10,pars)), "try-error"))
  expect_true(!inherits(try((g3*x3*p2)(0:10,pars)), "try-error"))
  expect_true(!inherits(try((g3*x1*p2)(0:10,pars)), "try-error"))
})
