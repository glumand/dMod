#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>
using namespace Rcpp;

#ifndef FCONE
#define FCONE
#endif

/*
 * ============================================================================
 * Batched Matrix Multiplication (BLAS-optimized, batch-last convention)
 * ============================================================================
 *
 * All arrays use [M,K,B] convention (batch index LAST).
 * Each slice [,,b] is standard column-major -> no copying needed.
 *
 * OPERATIONS:
 *
 *   bmm_lb:  [M,K,B] x [K,N]   -> [M,N,B]   (B dgemm / 1 batched call)
 *   bmm_rb:  [M,K]   x [K,N,B] -> [M,N,B]   (1 dgemm call always!)
 *   bmm_bb:  [M,K,B] x [K,N,B] -> [M,N,B]   (B dgemm / 1 batched call)
 *
 * If R is linked against OpenBLAS >= 0.3.19 or Intel MKL, configure
 * detects cblas_dgemm_batch_strided and defines HAVE_BATCH_STRIDED.
 * In that case bmm_lb and bmm_bb use a single batched BLAS call.
 *
 * ============================================================================
 */


/* --------------------------------------------------------------------------
 * Batched BLAS (optional, detected at package build)
 * -------------------------------------------------------------------------- */

#ifdef HAVE_BATCH_STRIDED

// Try to include proper CBLAS header
#if defined(__has_include)
#if __has_include(<mkl_cblas.h>)
#include <mkl_cblas.h>
#define USING_MKL_HEADER 1
#elif __has_include(<cblas.h>)
#include <cblas.h>
#define USING_CBLAS_HEADER 1
#else
// No header found, use manual declarations
#define NEED_CBLAS_DECL 1
#endif
#else
// Compiler doesn't support __has_include, try cblas.h
#if __cplusplus >= 201703L
#include <cblas.h>
#define USING_CBLAS_HEADER 1
#else
#define NEED_CBLAS_DECL 1
#endif
#endif

// Manual declarations if no header available
#ifdef NEED_CBLAS_DECL
enum CBLAS_LAYOUT { CblasRowMajor=101, CblasColMajor=102 };
enum CBLAS_TRANSPOSE { CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113 };

extern "C" void cblas_dgemm_batch_strided(
    const CBLAS_LAYOUT Layout,
    const CBLAS_TRANSPOSE TransA,
    const CBLAS_TRANSPOSE TransB,
    const int M, const int N, const int K,
    const double alpha,
    const double *A, const int lda, const int stride_a,
    const double *B, const int ldb, const int stride_b,
    const double beta,
    double *C, const int ldc, const int stride_c,
    const int batch_size);
#endif

#endif // HAVE_BATCH_STRIDED


/* --------------------------------------------------------------------------
 * Standard BLAS wrapper: C = A * B
 * -------------------------------------------------------------------------- */

inline void dgemm_nn(int M, int N, int K,
                     const double* A, int lda,
                     const double* B, int ldb,
                     double* C, int ldc) {
  const double one = 1.0, zero = 0.0;
  F77_CALL(dgemm)("N", "N", &M, &N, &K, &one,
           A, &lda, B, &ldb, &zero, C, &ldc FCONE FCONE);
}


/* --------------------------------------------------------------------------
 * bmm_lb: [M,K,B] x [K,N] -> [M,N,B]
 * --------------------------------------------------------------------------
 *
 * Left argument batched: same B multiplied into each batch of A.
 *
 * With HAVE_BATCH_STRIDED: 1 call (stride_B = 0 reuses the same matrix).
 * Fallback: B dgemm calls on contiguous slices.
 */
// [[Rcpp::export]]
NumericVector bmm_lb(NumericVector A, NumericVector B,
                     int Bn, int M, int K, int N) {
  
  // SAFETY CHECKS
  if (A.size() != M * K * Bn) {
    stop("bmm_lb: A size mismatch: got %d, expected %d (M=%d, K=%d, Bn=%d)", 
         A.size(), M*K*Bn, M, K, Bn);
  }
  if (B.size() != K * N) {
    stop("bmm_lb: B size mismatch: got %d, expected %d (K=%d, N=%d)", 
         B.size(), K*N, K, N);
  }
  
  NumericVector C(M * N * Bn);
  
#ifdef HAVE_BATCH_STRIDED
  // DEBUG OUTPUT
  Rprintf("\n=== bmm_lb DEBUG ===\n");
  Rprintf("Using batched BLAS (cblas_dgemm_batch_strided)\n");
#ifdef USING_MKL_HEADER
  Rprintf("Header: mkl_cblas.h\n");
#elif defined(USING_CBLAS_HEADER)
  Rprintf("Header: cblas.h\n");
#else
  Rprintf("Header: manual declaration\n");
#endif
  Rprintf("Parameters:\n");
  Rprintf("  M=%d, N=%d, K=%d, Bn=%d\n", M, N, K, Bn);
  Rprintf("  A: size=%d, lda=%d, stride=%d\n", A.size(), M, M*K);
  Rprintf("  B: size=%d, ldb=%d, stride=%d\n", B.size(), K, 0);
  Rprintf("  C: size=%d, ldc=%d, stride=%d\n", C.size(), M, M*N);
  Rprintf("Calling cblas_dgemm_batch_strided...\n");
  
  cblas_dgemm_batch_strided(
    CblasColMajor, CblasNoTrans, CblasNoTrans,
    M, N, K,
    1.0,
    &A[0], M, M * K,   // A strides between batches
    &B[0], K, 0,        // stride 0: reuse same B
    0.0,
    &C[0], M, M * N,
    Bn);
  
  Rprintf("SUCCESS! Call returned without crash.\n");
  Rprintf("===================\n\n");
#else
  const int MK = M * K;
  const int MN = M * N;
  for (int b = 0; b < Bn; ++b) {
    dgemm_nn(M, N, K,
             &A[b * MK], M,
             &B[0], K,
             &C[b * MN], M);
  }
#endif
  
  C.attr("dim") = IntegerVector::create(M, N, Bn);
  return C;
}


/* --------------------------------------------------------------------------
 * bmm_rb: [M,K] x [K,N,B] -> [M,N,B]
 * --------------------------------------------------------------------------
 *
 * Right argument batched. ONE dgemm call always — no batched BLAS needed!
 *
 * Key insight: [K,N,B] is contiguous as [K, N*B] in memory.
 * So: [M,K] x [K, N*B] = [M, N*B] which IS [M,N,B].
 */
// [[Rcpp::export]]
NumericVector bmm_rb(NumericVector A, NumericVector B,
                     int Bn, int M, int K, int N) {
  
  // SAFETY CHECKS
  if (A.size() != M * K) {
    stop("bmm_rb: A size mismatch: got %d, expected %d (M=%d, K=%d)", 
         A.size(), M*K, M, K);
  }
  if (B.size() != K * N * Bn) {
    stop("bmm_rb: B size mismatch: got %d, expected %d (K=%d, N=%d, Bn=%d)", 
         B.size(), K*N*Bn, K, N, Bn);
  }
  
  NumericVector C(M * N * Bn);
  
  Rprintf("\n=== bmm_rb (no batched BLAS) ===\n");
  Rprintf("M=%d, N=%d, K=%d, Bn=%d\n", M, N, K, Bn);
  Rprintf("Calling single dgemm: M=%d, N*Bn=%d, K=%d\n", M, N*Bn, K);
  
  dgemm_nn(M, N * Bn, K,
           &A[0], M,
           &B[0], K,
           &C[0], M);
           
           Rprintf("SUCCESS!\n===================\n\n");
           
           C.attr("dim") = IntegerVector::create(M, N, Bn);
           return C;
}


/* --------------------------------------------------------------------------
 * bmm_bb: [M,K,B] x [K,N,B] -> [M,N,B]
 * --------------------------------------------------------------------------
 *
 * Both arguments batched. Each slice [,,b] is contiguous column-major,
 * so no copying needed (unlike batch-first convention).
 *
 * With HAVE_BATCH_STRIDED: 1 call.
 * Fallback: B dgemm calls on contiguous slices.
 */
// [[Rcpp::export]]
NumericVector bmm_bb(NumericVector A, NumericVector B,
                     int Bn, int M, int K, int N) {
  
  // SAFETY CHECKS
  if (A.size() != M * K * Bn) {
    stop("bmm_bb: A size mismatch: got %d, expected %d (M=%d, K=%d, Bn=%d)", 
         A.size(), M*K*Bn, M, K, Bn);
  }
  if (B.size() != K * N * Bn) {
    stop("bmm_bb: B size mismatch: got %d, expected %d (K=%d, N=%d, Bn=%d)", 
         B.size(), K*N*Bn, K, N, Bn);
  }
  
  NumericVector C(M * N * Bn);
  
#ifdef HAVE_BATCH_STRIDED
  // DEBUG OUTPUT
  Rprintf("\n=== bmm_bb DEBUG ===\n");
  Rprintf("Using batched BLAS (cblas_dgemm_batch_strided)\n");
#ifdef USING_MKL_HEADER
  Rprintf("Header: mkl_cblas.h\n");
#elif defined(USING_CBLAS_HEADER)
  Rprintf("Header: cblas.h\n");
#else
  Rprintf("Header: manual declaration\n");
#endif
  Rprintf("Parameters:\n");
  Rprintf("  M=%d, N=%d, K=%d, Bn=%d\n", M, N, K, Bn);
  Rprintf("  A: size=%d, lda=%d, stride=%d\n", A.size(), M, M*K);
  Rprintf("  B: size=%d, ldb=%d, stride=%d\n", B.size(), K, K*N);
  Rprintf("  C: size=%d, ldc=%d, stride=%d\n", C.size(), M, M*N);
  Rprintf("Calling cblas_dgemm_batch_strided...\n");
  
  cblas_dgemm_batch_strided(
    CblasColMajor, CblasNoTrans, CblasNoTrans,
    M, N, K,
    1.0,
    &A[0], M, M * K,
    &B[0], K, K * N,
    0.0,
    &C[0], M, M * N,
    Bn);
  
  Rprintf("SUCCESS! Call returned without crash.\n");
  Rprintf("===================\n\n");
#else
  const int MK = M * K;
  const int KN = K * N;
  const int MN = M * N;
  for (int b = 0; b < Bn; ++b) {
    dgemm_nn(M, N, K,
             &A[b * MK], M,
             &B[b * KN], K,
             &C[b * MN], M);
  }
#endif
  
  C.attr("dim") = IntegerVector::create(M, N, Bn);
  return C;
}