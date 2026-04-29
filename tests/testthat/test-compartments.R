context("compartments")

test_that("back-compat flat volumes auto-translate to compartment IDs", {
  # Flat volumes c(A=V1, B=V1, C=V2) should produce two compartments c1, c2.
  f <- eqnlist(
    smatrix = matrix(c(-1, 1, 0, 0, -1, 1), nrow = 2, byrow = TRUE,
                     dimnames = list(NULL, c("A", "B", "C"))),
    states = c("A", "B", "C"),
    rates = c("k1*A", "k2*B"),
    description = c("r1", "r2"),
    volumes = c(A = "V1", B = "V1", C = "V2")
  )
  expect_true(is.eqnlist(f))
  expect_equal(sort(names(f$compartments)), c("c1", "c2"))
  # States with same volume expression share a compartment.
  expect_equal(unname(f$compartmentOf[["A"]]), unname(f$compartmentOf[["B"]]))
  expect_false(unname(f$compartmentOf[["A"]]) == unname(f$compartmentOf[["C"]]))
  # Derived $volumes view matches.
  expect_equal(unname(f$volumes[["A"]]), "V1")
  expect_equal(unname(f$volumes[["C"]]), "V2")
})

test_that("explicit compartments/compartmentOf are respected", {
  f <- eqnlist(
    smatrix = matrix(c(-1, 1), nrow = 1, dimnames = list(NULL, c("A", "B"))),
    states = c("A", "B"),
    rates = "k*A",
    description = "transport",
    compartments = list(cyt = "V_cyt", nuc = "V_nuc"),
    compartmentOf = c(A = "cyt", B = "nuc")
  )
  expect_true(is.eqnlist(f))
  expect_equal(unname(f$compartmentOf[["A"]]), "cyt")
  expect_equal(unname(f$compartmentOf[["B"]]), "nuc")
  expect_equal(f$compartments$cyt$volume, "V_cyt")
})

test_that("unassigned states fall into an implicit default compartment", {
  f <- NULL
  f <- addReaction(f, "A", "B", "k*A")
  expect_equal(names(f$compartments), "default")
  expect_equal(f$compartments$default$volume, "1")
  expect_true(all(f$compartmentOf %in% "default"))
})

test_that("is.eqnlist rejects broken compartmentOf references", {
  f <- eqnlist(
    smatrix = matrix(c(-1, 1), nrow = 1, dimnames = list(NULL, c("A", "B"))),
    states = c("A", "B"), rates = "k*A", description = "r",
    compartments = list(cyt = "V_cyt"),
    compartmentOf = c(A = "cyt", B = "cyt")
  )
  # Corrupt it: point B at a nonexistent compartment
  f$compartmentOf[["B"]] <- "nonexistent"
  expect_false(is.eqnlist(f))
})

test_that("addReaction respects the compartment argument for new states", {
  f <- addReaction(NULL, "A", "B", "k*A", compartment = "cyt")
  expect_equal(unname(f$compartmentOf[["A"]]), "cyt")
  expect_equal(unname(f$compartmentOf[["B"]]), "cyt")
  f <- addReaction(f, "B", "C", "k2*B", compartment = "nuc")
  expect_equal(unname(f$compartmentOf[["B"]]), "cyt")   # existing state keeps its compartment
  expect_equal(unname(f$compartmentOf[["C"]]), "nuc")   # new state picks up the arg
})

test_that("subset.eqnlist drops unreferenced compartments", {
  f <- NULL
  f <- addReaction(f, "A", "B", "k1*A", "r1")
  f <- addReaction(f, "B", "C", "k2*B", "r2")
  f$compartments <- list(c1 = list(volume = "V1", rule = NULL),
                         c2 = list(volume = "V2", rule = NULL))
  f$compartmentOf <- c(A = "c1", B = "c1", C = "c2")
  f$volumes <- c(A = "V1", B = "V1", C = "V2")

  f_sub <- subset(f, Description == "r1")
  expect_equal(names(f_sub$compartments), "c1")
  expect_equal(sort(names(f_sub$compartmentOf)), c("A", "B"))
})

test_that("c.eqnlist merges on equal-volume and errors on conflict", {
  g1 <- eqnlist(smatrix = matrix(-1, nrow = 1, dimnames = list(NULL, "X")),
                states = "X", rates = "k*X", description = "deg_X",
                compartments = list(cyt = "V_cyt"), compartmentOf = c(X = "cyt"))
  g2 <- eqnlist(smatrix = matrix(-1, nrow = 1, dimnames = list(NULL, "Y")),
                states = "Y", rates = "k*Y", description = "deg_Y",
                compartments = list(nuc = "V_nuc"), compartmentOf = c(Y = "nuc"))
  g <- c(g1, g2)
  expect_equal(sort(names(g$compartments)), c("cyt", "nuc"))

  # Same ID, different volume -> conflict error
  h1 <- eqnlist(smatrix = matrix(-1, nrow = 1, dimnames = list(NULL, "X")),
                states = "X", rates = "k*X", description = "r",
                compartments = list(cyt = "V1"), compartmentOf = c(X = "cyt"))
  h2 <- eqnlist(smatrix = matrix(-1, nrow = 1, dimnames = list(NULL, "Y")),
                states = "Y", rates = "k*Y", description = "r",
                compartments = list(cyt = "V2"), compartmentOf = c(Y = "cyt"))
  expect_error(c(h1, h2), "conflict")
})

test_that("getFluxes produces identical output for legacy single-compartment models", {
  # The README enzyme-kinetics pipeline: sanity check that nothing in the
  # default-compartment code path changes the flux expressions.
  f <- NULL
  f <- addReaction(f, "Enz + Sub", "Compl", "k1*Enz*Sub", "production")
  f <- addReaction(f, "Compl", "Enz + Sub", "k2*Compl", "decay")
  f <- addReaction(f, "Compl", "Enz + Prod", "k3*Compl", "prod_product")
  f <- addReaction(f, "Enz", "", "k4*Enz", "deg_enzyme")
  fl <- getFluxes(f)
  # No volume-ratio factors should appear — every reaction lives in the same
  # (default) compartment.
  flat <- unlist(fl)
  expect_false(any(grepl("\\(1/1\\)", flat)))
  # Structure matches expected states
  expect_equal(sort(names(fl)), sort(c("Enz", "Sub", "Compl", "Prod")))
})

test_that("getFluxes emits cross-compartment volume ratio", {
  f <- eqnlist(
    smatrix = matrix(c(-1, 1), nrow = 1, dimnames = list(NULL, c("A", "B"))),
    states = c("A", "B"),
    rates = "k*A",
    description = "transport",
    compartments = list(cyt = "V_cyt", nuc = "V_nuc"),
    compartmentOf = c(A = "cyt", B = "nuc")
  )
  fl <- getFluxes(f)
  # A (cyt, origin) should not carry a ratio; B (nuc, destin) should carry V_cyt/V_nuc.
  expect_false(grepl("V_cyt/V_nuc", fl$A))
  expect_true(grepl("V_cyt/V_nuc", fl$B))
})

test_that("getFluxes emits a dilution term when a compartment has a rule", {
  f <- eqnlist(
    smatrix = matrix(c(-1, 1), nrow = 1, dimnames = list(NULL, c("A", "B"))),
    states = c("A", "B"),
    rates = "k*A",
    description = "r",
    compartments = list(cyt = list(volume = "V_cyt", rule = "mu*V_cyt")),
    compartmentOf = c(A = "cyt", B = "cyt")
  )
  fl <- getFluxes(f)
  # Dilution term: -(A)*(mu*V_cyt)/(V_cyt) and -(B)*(mu*V_cyt)/(V_cyt)
  expect_true(any(grepl("dilution", names(fl$A))))
  expect_true(any(grepl("dilution", names(fl$B))))
  expect_true(any(grepl("mu\\*V_cyt", fl$A)))
})

test_that("print.eqnlist does not surface the default compartment", {
  f <- NULL
  f <- addReaction(f, "A", "B", "k*A")
  captured <- paste(capture.output(print(f)), collapse = "\n")
  expect_false(grepl("Compartments:", captured))
})

test_that("print.eqnlist shows compartments when meaningful", {
  f <- eqnlist(
    smatrix = matrix(c(-1, 1), nrow = 1, dimnames = list(NULL, c("A", "B"))),
    states = c("A", "B"), rates = "k*A", description = "r",
    compartments = list(cyt = "V_cyt", nuc = "V_nuc"),
    compartmentOf = c(A = "cyt", B = "nuc")
  )
  captured <- paste(capture.output(print(f)), collapse = "\n")
  expect_true(grepl("Compartments:", captured))
  expect_true(grepl("cyt", captured))
  expect_true(grepl("nuc", captured))
})

test_that("conservedQuantities accepts a weight argument", {
  S <- matrix(c(-1, 1, -1, 1), nrow = 2, byrow = TRUE,
              dimnames = list(NULL, c("A", "B")))

  unweighted <- conservedQuantities(S)                  # default "none"
  weighted   <- conservedQuantities(S, weight = "volume",
                                    volumes = c(A = 2, B = 3))
  expect_true(!is.null(unweighted))
  # For A <-> B with volumes (2, 3), mass conservation weights: 2*A + 3*B conserved
  expect_true(!is.null(weighted))

  expect_error(conservedQuantities(S, weight = "volume"), "volumes")
  expect_error(conservedQuantities(S, weight = "volume",
                                    volumes = c(A = "V1", B = 3)),
                "numeric")
})

test_that("rateCompartment enables reactions with educts across compartments", {
  # Membrane binding: L in extraCellular, R in cytosol. The rate expression
  # k*L*R is a concentration-rate in V_ext.
  f <- NULL
  f <- addReaction(f, "",          "L", "k_prod",       compartment = "extraCellular")
  f <- addReaction(f, "",          "R", "k_Rprod",      compartment = "cytosol")
  f <- addReaction(f, "L + R",     "Compl", "k_on*L*R",
                   compartment = "cytosol",
                   rateCompartment = "extraCellular")
  f$compartments$extraCellular$volume <- "V_ext"
  f$compartments$cytosol$volume       <- "V_cyt"

  # reactionCompartment vector stored
  expect_equal(f$reactionCompartment[3], "extraCellular")
  expect_true(is.na(f$reactionCompartment[1]))

  # uniform flux formula: V_ref = V_ext for the binding reaction,
  # so R (cytosol) is scaled by V_ext/V_cyt, and L (ext) has ratio 1.
  fl <- getFluxes(f)
  # L_ext educt: no ratio (ratio is V_ext/V_ext == "")
  expect_true(any(grepl("-1\\*\\(k_on\\*L\\*R\\)$", fl$L)))
  # R cytosol educt: scaled by V_ext/V_cyt
  expect_true(any(grepl("k_on\\*L\\*R\\)\\*\\(V_ext/V_cyt\\)", fl$R)))
  # Compl cytosol product: same scaling
  expect_true(any(grepl("k_on\\*L\\*R\\)\\*\\(V_ext/V_cyt\\)", fl$Compl)))
})

test_that("getFluxes errors with a helpful message when educts span compartments without annotation", {
  f <- eqnlist(
    smatrix = matrix(c(-1, -1, 1), nrow = 1, dimnames = list(NULL, c("L", "R", "Compl"))),
    states = c("L", "R", "Compl"), rates = "k*L*R", description = "bind",
    compartments = list(ext = "V_ext", cyt = "V_cyt"),
    compartmentOf = c(L = "ext", R = "cyt", Compl = "cyt")
  )
  expect_error(getFluxes(f), "reactionCompartment")
})

test_that("subset.eqnlist and c.eqnlist preserve reactionCompartment", {
  f <- NULL
  f <- addReaction(f, "",       "L", "k1",     "produce_L",
                   compartment = "ext")
  f <- addReaction(f, "L + R",  "C", "k2*L*R", "bind",
                   compartment = "cyt", rateCompartment = "ext")
  f$compartments$ext$volume <- "V_ext"
  f$compartments$cyt$volume <- "V_cyt"

  expect_equal(f$reactionCompartment, c(NA, "ext"))

  # subset: keep only the binding reaction
  f_sub <- subset(f, Description == "bind")
  expect_equal(f_sub$reactionCompartment, "ext")

  # c.eqnlist: annotation is preserved across concatenation
  f_first  <- subset(f, Description == "produce_L")
  f_second <- subset(f, Description == "bind")
  f_comb   <- c(f_first, f_second)
  expect_equal(unname(f_comb$reactionCompartment), c(NA, "ext"))
})

test_that("getParameters.eqnlist picks up symbolic compartment volumes", {
  f <- eqnlist(
    smatrix = matrix(c(-1, 1), nrow = 1, dimnames = list(NULL, c("A", "B"))),
    states = c("A", "B"), rates = "k*A", description = "r",
    compartments = list(cyt = "V_cyt", nuc = "V_nuc"),
    compartmentOf = c(A = "cyt", B = "nuc")
  )
  params <- getParameters(f)
  expect_true("V_cyt" %in% params)
  expect_true("V_nuc" %in% params)
  expect_true("k" %in% params)
})
