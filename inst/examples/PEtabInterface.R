\dontrun{
  ## Import a PEtab v1 problem from disk. The YAML manifest lives next to
  ## the four TSVs and the SBML model file. Native artefacts (`.c`, `.cpp`,
  ## `.so`) land in the current working directory — set tempdir() first
  ## to avoid polluting the project tree.
  setwd(tempdir())

  yaml_path <- system.file("PEtabTests/0001/_0001.yaml", package = "dMod")
  petab <- importPEtab(yaml_path, solver = "deSolve")

  print(petab)

  ## petab is a plain list — every slot is a regular dMod object.
  ## Evaluate the objective at the nominal pouter values:
  petab$obj(petab$pouter, fixed = petab$fixed)$value

  ## Predict and plot:
  times <- seq(0, 10, length.out = 51)
  prediction <- petab$prd(times, c(petab$pouter, petab$fixed))
  plot(prediction, petab$data)

  ## Fit:
  fit <- trust(petab$obj, petab$pouter, fixed = petab$fixed,
               rinit = 1, rmax = 10)

  ## Round-trip back to PEtab on disk. exportPEtabObject takes the full
  ## petabProblem list; for a dMod-native problem (no PEtab origin), use
  ## exportPEtab(data, obj, model, g, x, p, pouter, lower, upper, ...).
  out_dir <- file.path(tempdir(), "petab_export")
  yaml_out <- exportPEtabObject(petab, out_dir, overwrite = TRUE)
  petab2 <- importPEtab(yaml_out, solver = "deSolve")
}
