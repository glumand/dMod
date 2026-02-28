.onAttach <- function(libname, pkgname) {
  blas_info <- extSoftVersion()["BLAS"]
  
  if (is.loaded("_dMod_has_batch_gemm") && has_batch_gemm()) {
    packageStartupMessage("\nBLAS: ", blas_info)
    # Set single-threaded BLAS to avoid nested parallelism
    Sys.setenv(
      OMP_NUM_THREADS       = "1",
      MKL_NUM_THREADS       = "1",
      MKL_THREADING_LAYER   = "SEQUENTIAL",
      OPENBLAS_NUM_THREADS   = "1",
      BLIS_NUM_THREADS       = "1"
    )
  } else {
    packageStartupMessage("\nBLAS: ", blas_info)
    packageStartupMessage("cblas_dgemm_batch not available (using fallback implementation)")
    packageStartupMessage("For instructions on enabling optimized batched GEMM, run blasHelp().")
  }
}


#' Display instructions for enabling optimized batched matrix multiplication
#' Shows platform-specific instructions for installing and configuring a BLAS
#' library that provides \code{cblas_dgemm_batch}.
#'
#' @export
blasHelp <- function() {
  os <- Sys.info()[["sysname"]]
  if (os == "Windows") {
    cat(
      "  ┌────────────────────────────────────────────────────────────────────────────────────────────┐\n",
      "  │  For optimized batched matrix multiplication, a BLAS with cblas_dgemm_batch is needed.     │\n",
      "  │                                                                                            │\n",
      "  │  Recommended: Intel oneAPI Math Kernel Library (MKL)                                       │\n",
      "  │    All MKL versions include cblas_dgemm_batch.                                             │\n",
      "  │                                                                                            │\n",
      "  │  1. Install Intel oneAPI Math Kernel Library:                                              │\n",
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
      "  │  4. Reinstall dMod to enable the optimized backend.                                        │\n",
      "  └────────────────────────────────────────────────────────────────────────────────────────────┘\n"
      , sep = "")
  } else if (os == "Darwin") {
    cat(
      "  ┌────────────────────────────────────────────────────────────────────────────────────────────┐\n",
      "  │  For optimized batched matrix multiplication, a BLAS with cblas_dgemm_batch is needed.     │\n",
      "  │                                                                                            │\n",
      "  │  Note: Apple's Accelerate framework does NOT provide cblas_dgemm_batch.                    │\n",
      "  │                                                                                            │\n",
      "  │  Recommended: Intel oneAPI Math Kernel Library (MKL)                                       │\n",
      "  │    All MKL versions include cblas_dgemm_batch.                                             │\n",
      "  │                                                                                            │\n",
      "  │       conda install mkl mkl-devel  (in conda env)                                         │\n",
      "  │                                                                                            │\n",
      "  │  Alternative: OpenBLAS (>= 0.3.19)                                                         │\n",
      "  │    NOTE: Homebrew's OpenBLAS may not include cblas_dgemm_batch.                            │\n",
      "  │    If not, use MKL or build OpenBLAS from source:                                          │\n",
      "  │    https://github.com/OpenMathLib/OpenBLAS                                                 │\n",
      "  │                                                                                            │\n",
      "  │       brew install openblas                                                                │\n",
      "  │                                                                                            │\n",
      "  │  R links against the BLAS found at install time. Either reinstall R or set                 │\n",
      "  │  DYLD_LIBRARY_PATH to point to your preferred BLAS.                                       │\n",
      "  │                                                                                            │\n",
      "  │  Reinstall dMod to pick up the optimized backend.                                          │\n",
      "  └────────────────────────────────────────────────────────────────────────────────────────────┘\n"
      , sep = "")
  } else {
    cat(
      "  ┌────────────────────────────────────────────────────────────────────────────────────────────┐\n",
      "  │  For optimized batched matrix multiplication, a BLAS with cblas_dgemm_batch is needed.     │\n",
      "  │                                                                                            │\n",
      "  │  Recommended: Intel oneAPI Math Kernel Library (MKL)                                       │\n",
      "  │    All MKL versions include cblas_dgemm_batch.                                             │\n",
      "  │                                                                                            │\n",
      "  │       Ubuntu/Debian:  sudo apt install intel-mkl                                           │\n",
      "  │       Fedora/RHEL:    sudo dnf install intel-oneapi-mkl-devel                              │\n",
      "  │                                                                                            │\n",
      "  │  Alternative: OpenBLAS (>= 0.3.19)                                                         │\n",
      "  │    NOTE: Some Linux distributions (e.g. Ubuntu) ship OpenBLAS without                      │\n",
      "  │    cblas_dgemm_batch even in recent versions. In that case, use MKL or                     │\n",
      "  │    build OpenBLAS from source: https://github.com/OpenMathLib/OpenBLAS                     │\n",
      "  │                                                                                            │\n",
      "  │       Ubuntu/Debian:  sudo apt install libopenblas-dev                                     │\n",
      "  │       Fedora/RHEL:    sudo dnf install openblas-devel                                      │\n",
      "  │                                                                                            │\n",
      "  │  IMPORTANT: On Linux, ALL FOUR alternatives must be switched consistently:                 │\n",
      "  │                                                                                            │\n",
      "  │       sudo update-alternatives --config libblas.so.3-x86_64-linux-gnu                      │\n",
      "  │       sudo update-alternatives --config libblas.so-x86_64-linux-gnu                        │\n",
      "  │       sudo update-alternatives --config liblapack.so.3-x86_64-linux-gnu                    │\n",
      "  │       sudo update-alternatives --config liblapack.so-x86_64-linux-gnu                      │\n",
      "  │                                                                                            │\n",
      "  │  On HPC systems, load the appropriate module (e.g. module load numlib/mkl/2022.2)          │\n",
      "  │  before running R CMD INSTALL.                                                             │\n",
      "  │                                                                                            │\n",
      "  │  Reinstall dMod to pick up the optimized backend.                                          │\n",
      "  └────────────────────────────────────────────────────────────────────────────────────────────┘\n"
      , sep = "")
  }
}