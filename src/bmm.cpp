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
 * Batched BLAS (optional, detected at package built)
 * -------------------------------------------------------------------------- */

#ifdef HAVE_BATCH_STRIDED
extern "C" void cblas_dgemm_batch_strided(
    int layout, int transa, int transb,
    int M, int N, int K,
    double alpha,
    const double *A, int lda, int stride_a,
    const double *B, int ldb, int stride_b,
    double beta,
    double *C, int ldc, int stride_c,
    int batch_size);

static const int CBLAS_COL_MAJOR = 102;
static const int CBLAS_NO_TRANS  = 111;
#endif


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
  
  NumericVector C(M * N * Bn);
  
#ifdef HAVE_BATCH_STRIDED
  cblas_dgemm_batch_strided(
    CBLAS_COL_MAJOR, CBLAS_NO_TRANS, CBLAS_NO_TRANS,
    M, N, K,
    1.0,
    &A[0], M, M * K,   // A strides between batches
    &B[0], K, 0,        // stride 0: reuse same B
    0.0,
    &C[0], M, M * N,
    Bn);
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
  
  NumericVector C(M * N * Bn);
  
  dgemm_nn(M, N * Bn, K,
           &A[0], M,
           &B[0], K,
           &C[0], M);
           
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
  
  NumericVector C(M * N * Bn);
  
#ifdef HAVE_BATCH_STRIDED
  cblas_dgemm_batch_strided(
    CBLAS_COL_MAJOR, CBLAS_NO_TRANS, CBLAS_NO_TRANS,
    M, N, K,
    1.0,
    &A[0], M, M * K,
    &B[0], K, K * N,
    0.0,
    &C[0], M, M * N,
    Bn);
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