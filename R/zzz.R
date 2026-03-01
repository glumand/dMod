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
#' library that provides \code{cblas_dgemm_batch} and
#' \code{cblas_dgemm_batch_strided}.
#'
#' @export
blasHelp <- function() {
  os <- Sys.info()[["sysname"]]
  if (os == "Windows") {
    cat(
      "  \u250c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2510\n",
      "  \u2502  For optimized batched matrix multiplication, a BLAS with cblas_dgemm_batch and            \u2502\n",
      "  \u2502  cblas_dgemm_batch_strided is needed.                                                      \u2502\n",
      "  \u2502                                                                                            \u2502\n",
      "  \u2502  Recommended: Intel oneAPI Math Kernel Library (MKL)                                       \u2502\n",
      "  \u2502    All MKL versions include both batched BLAS functions.                                   \u2502\n",
      "  \u2502                                                                                            \u2502\n",
      "  \u2502  1. Install Intel oneAPI Math Kernel Library:                                              \u2502\n",
      "  \u2502       https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html  \u2502\n",
      "  \u2502                                                                                            \u2502\n",
      "  \u2502  2. Check if MKLROOT is set:                                                               \u2502\n",
      "  \u2502       In R, run:  Sys.getenv('MKLROOT')                                                    \u2502\n",
      "  \u2502       - If it prints a valid path, MKLROOT is set.                                         \u2502\n",
      "  \u2502       - If it is empty or invalid, set it (adapt if necessary):                            \u2502\n",
      "  \u2502             Sys.setenv(MKLROOT='C:/Program Files (x86)/Intel/oneAPI/mkl/latest')           \u2502\n",
      "  \u2502                                                                                            \u2502\n",
      "  \u2502  3. Ensure mkl_rt.dll is in your PATH:                                                     \u2502\n",
      "  \u2502       In R, run:                                                                           \u2502\n",
      "  \u2502             Sys.setenv(PATH=paste(Sys.getenv('PATH'),                                      \u2502\n",
      "  \u2502                                  paste0(Sys.getenv('MKLROOT'), '/redist/intel64'),         \u2502\n",
      "  \u2502                                  sep=.Platform$path.sep))                                  \u2502\n",
      "  \u2502                                                                                            \u2502\n",
      "  \u2502  4. Reinstall dMod to enable the optimized backend.                                        \u2502\n",
      "  \u2514\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2518\n"
      , sep = "")
  } else if (os == "Darwin") {
    cat(
      "  \u250c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2510\n",
      "  \u2502  For optimized batched matrix multiplication, a BLAS with cblas_dgemm_batch and            \u2502\n",
      "  \u2502  cblas_dgemm_batch_strided is needed.                                                      \u2502\n",
      "  \u2502                                                                                            \u2502\n",
      "  \u2502  Note: Apple's Accelerate framework does NOT provide these functions.                      \u2502\n",
      "  \u2502                                                                                            \u2502\n",
      "  \u2502  Recommended: Intel oneAPI Math Kernel Library (MKL)                                       \u2502\n",
      "  \u2502    All MKL versions include both batched BLAS functions.                                   \u2502\n",
      "  \u2502                                                                                            \u2502\n",
      "  \u2502       conda install mkl mkl-devel  (in conda env)                                          \u2502\n",
      "  \u2502                                                                                            \u2502\n",
      "  \u2502  Alternative: OpenBLAS (>= 0.3.28)                                                         \u2502\n",
      "  \u2502    NOTE: Batched BLAS was added in OpenBLAS 0.3.28 and may be unstable in early            \u2502\n",
      "  \u2502    versions. Homebrew's OpenBLAS may not include cblas_dgemm_batch.                        \u2502\n",
      "  \u2502    If not, use MKL or build OpenBLAS from source:                                          \u2502\n",
      "  \u2502    https://github.com/OpenMathLib/OpenBLAS                                                 \u2502\n",
      "  \u2502                                                                                            \u2502\n",
      "  \u2502       brew install openblas                                                                \u2502\n",
      "  \u2502                                                                                            \u2502\n",
      "  \u2502  R links against the BLAS found at install time. Either reinstall R or set                 \u2502\n",
      "  \u2502  DYLD_LIBRARY_PATH to point to your preferred BLAS.                                        \u2502\n",
      "  \u2502                                                                                            \u2502\n",
      "  \u2502  Reinstall dMod to enable  up the optimized backend.                                       \u2502\n",
      "  \u2514\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2518\n"
      , sep = "")
  } else {
    cat(
      "  \u250c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2510\n",
      "  \u2502  For optimized batched matrix multiplication, a BLAS with cblas_dgemm_batch and            \u2502\n",
      "  \u2502  cblas_dgemm_batch_strided is needed.                                                      \u2502\n",
      "  \u2502                                                                                            \u2502\n",
      "  \u2502  Recommended: Intel oneAPI Math Kernel Library (MKL)                                       \u2502\n",
      "  \u2502                                                                                            \u2502\n",
      "  \u2502       Ubuntu/Debian:  sudo apt install intel-mkl                                           \u2502\n",
      "  \u2502       Fedora/RHEL:    sudo dnf install intel-oneapi-mkl-devel                              \u2502\n",
      "  \u2502                                                                                            \u2502\n",
      "  \u2502  Alternative: OpenBLAS (>= 0.3.28)                                                         \u2502\n",
      "  \u2502    NOTE: Batched BLAS was added in OpenBLAS 0.3.28 and may be unstable in early            \u2502\n",
      "  \u2502    versions. Some Linux distributions ship OpenBLAS without it or with a                   \u2502\n",
      "  \u2502    non-functional implementation.                                                          \u2502\n",
      "  \u2502                                                                                            \u2502\n",
      "  \u2502       Ubuntu/Debian:  sudo apt install libopenblas-dev                                     \u2502\n",
      "  \u2502       Fedora/RHEL:    sudo dnf install openblas-devel                                      \u2502\n",
      "  \u2502                                                                                            \u2502\n",
      "  \u2502  IMPORTANT: On Linux, ALL FOUR alternatives must be switched consistently:                 \u2502\n",
      "  \u2502                                                                                            \u2502\n",
      "  \u2502       sudo update-alternatives --config libblas.so.3-x86_64-linux-gnu                      \u2502\n",
      "  \u2502       sudo update-alternatives --config libblas.so-x86_64-linux-gnu                        \u2502\n",
      "  \u2502       sudo update-alternatives --config liblapack.so.3-x86_64-linux-gnu                    \u2502\n",
      "  \u2502       sudo update-alternatives --config liblapack.so-x86_64-linux-gnu                      \u2502\n",
      "  \u2502                                                                                            \u2502\n",
      "  \u2502  Reinstall dMod to pick up the optimized backend.                                          \u2502\n",
      "  \u2514\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2518\n"
      , sep = "")
  }
}