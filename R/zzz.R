.onAttach <- function(libname, pkgname) {
  if (is.loaded("_dMod_has_batch_gemm") && has_batch_gemm()) {
    packageStartupMessage("\nIntel\u00AE oneAPI Math Kernel Library (oneMKL) available")
    Sys.setenv(MKL_NUM_THREADS = "1", MKL_THREADING_LAYER = "SEQUENTIAL")
  } else {
    packageStartupMessage("\nIntel\u00AE oneAPI Math Kernel Library (oneMKL) not available")
    packageStartupMessage("For instructions on enabling the MKL backend for faster batched matrix multiplication, run mklHelp().")
  }
}

#' @export
mklHelp <- function() {
  os <- Sys.info()[["sysname"]]
  if (os == "Windows") {
    cat(
      "  ┌────────────────────────────────────────────────────────────────────────────────────────────┐\n",
      "  │  For maximum performance of batched matrix multiplication:                                 │\n",
      "  │                                                                                            │\n",
      "  │  1. Install Intel oneAPI Math Kernel Library if not already installed:                     │\n",
      "  │       https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html  │\n",
      "  │                                                                                            │\n",
      "  │  2. Check if MKLROOT is set:                                                               │\n",
      "  │       In R, run:  Sys.getenv('MKLROOT')                                                    │\n",
      "  │       - If it prints a valid path, MKLROOT is set.                                         │\n",
      "  │       - If it is empty or invalid, set it:                                                 │\n",
      "  │             Sys.setenv(MKLROOT='C:/Program Files (x86)/Intel/oneAPI/mkl/latest')           │\n",
      "  │                                                                                            │\n",
      "  │  3. Ensure mkl_rt.dll is in your PATH:                                                     │\n",
      "  │       In R, run:                                                                           │\n",
      "  │             Sys.setenv(PATH=paste(Sys.getenv('PATH'),                                      │\n",
      "  │                                  paste0(Sys.getenv('MKLROOT'), '/redist/intel64'),         │\n",
      "  │                                  sep=.Platform$path.sep))                                  │\n",
      "  │                                                                                            │\n",
      "  │  4. Reinstall dMod to enable MKL backend.                                                  │\n",
      "  └────────────────────────────────────────────────────────────────────────────────────────────┘\n"
      , sep = "")
  } else {
    cat(
      "  ┌────────────────────────────────────────────────────────────────────────────────────────────┐\n",
      "  │  For maximum performance of batched matrix multiplication:                                 │\n",
      "  │                                                                                            │\n",
      "  │  1. Install Intel oneAPI Math Kernel Library if not already installed:                     │\n",
      "  │                                                                                            │\n",
      "  │       Ubuntu/Debian:  sudo apt install intel-mkl                                           │\n",
      "  │       Fedora/RHEL:    sudo dnf install intel-oneapi-mkl-devel                              │\n",
      "  │       macOS/conda:    conda install mkl mkl-devel                                          │\n",
      "  │                                                                                            │\n",
      "  │  2. Make sure, that MKL is set as R's BLAS backend so that R links against it:             │\n",
      "  │                                                                                            │\n",
      "  │       Ubuntu/Debian:  sudo update-alternatives --config libblas.so.3-x86_64-linux-gnu      │\n",
      "  │                                                                                            │\n",
      "  │  3. Reinstall dMod to pick up the optimized backend.                                       │\n",
      "  └────────────────────────────────────────────────────────────────────────────────────────────┘\n"
      , sep = "")
  }
}