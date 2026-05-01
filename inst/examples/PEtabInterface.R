\dontrun{
  ## Import a PEtab problem from disk. importPEtab dispatches on the YAML's
  ## `format_version` key, so v1 and v2 manifests are both supported. The
  ## TSVs and the SBML model live next to the YAML; native artefacts
  ## (`.c`, `.cpp`, `.so`) land in the current working directory —
  ## `setwd(tempdir())` keeps them out of the project tree.
  setwd(tempdir())

  ## v2 example (the bundled v2 test suite):
  yamlPath <- system.file("PEtabTests/v2/0001/_0001.yaml", package = "dMod")
  petab <- importPEtab(yamlPath, solver = "deSolve")

  ## v1 still works — same call, different YAML schema:
  ## yamlPath <- system.file("PEtabTests/0001/_0001.yaml", package = "dMod")

  print(petab)

  ## petab is a plain list — every slot is a regular dMod object.
  ## The objective has the PEtab `fixed` parameters baked in, so calling
  ## obj(bestfit) evaluates the likelihood at the nominal estimate:
  petab$obj(petab$bestfit)$value

  ## Predict and plot:
  times <- seq(0, 10, length.out = 51)
  prediction <- petab$prd(times, c(petab$bestfit,
                                   attr(petab, "petab_meta")$fixed))
  plot(prediction, petab$dataList)

  ## Fit:
  fit <- trust(petab$obj, petab$bestfit, rinit = 1, rmax = 10)

  ## Round-trip back to PEtab on disk. exportPEtabObject takes the full
  ## petabProblem list; for a dMod-native problem (no PEtab origin), use
  ## exportPEtab(data, reactions, observables, p, pouter, ...).
  ## formatVersion defaults to "2.0.0"; pass "1" to write the legacy
  ## schema instead.
  outDir <- file.path(tempdir(), "petab_export")
  yamlOut <- exportPEtabObject(petab, outDir, formatVersion = "2.0.0",
                               overwrite = TRUE)
  petab2 <- importPEtab(yamlOut, solver = "deSolve")
}
