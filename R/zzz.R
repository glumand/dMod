.onAttach <- function(libname, pkgname) {
  if (is.loaded("_dMod_has_batch_strided") && has_batch_strided()) {
    packageStartupMessage("\nIntel\u00AE oneAPI Math Kernel Library (oneMKL) available")
    Sys.setenv(MKL_NUM_THREADS = "1", MKL_THREADING_LAYER = "SEQUENTIAL")
  }
}