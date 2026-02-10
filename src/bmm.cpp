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
 *   bmm_lb:  [M,K,B] x [K,N]   -> [M,N,B]   (B dgemm / 1 batched dgemm call)
 *   bmm_rb:  [M,K]   x [K,N,B] -> [M,N,B]   (1 dgemm call)
 *   bmm_bb:  [M,K,B] x [K,N,B] -> [M,N,B]   (B dgemm / 1 batched dgemm call)
 *
 * If R is linked against Intel MKL, configure detects
 * cblas_dgemm_batch_strided and defines HAVE_BATCH_STRIDED.
 * In that case bmm_lb and bmm_bb use a single batched BLAS call.
 *
 * ============================================================================
 */

#include <Rcpp.h>

#ifdef HAVE_BATCH_STRIDED
#include <mkl.h>
#else
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>
#endif

using namespace Rcpp;

#ifndef FCONE
#define FCONE
#endif


/* --------------------------------------------------------------------------
 * Query whether batched BLAS is available
 * -------------------------------------------------------------------------- */

// [[Rcpp::export]]
bool has_batch_strided() {
#ifdef HAVE_BATCH_STRIDED
  return true;
#else
  return false;
#endif
}


#ifndef HAVE_BATCH_STRIDED
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
#endif


/* --------------------------------------------------------------------------
 * bmm_lb: [M,K,B] x [K,N] -> [M,N,B]
 * --------------------------------------------------------------------------
 *
 * Left argument batched: same B multiplied into each batch of A.
 *
 * With HAVE_BATCH_STRIDED: 1 stided dgemm call (stride_B = 0 reuses the same matrix).
 * Fallback: B dgemm calls on contiguous slices.
 */
// [[Rcpp::export]]
NumericVector bmm_lb(NumericVector A, NumericVector B,
                     int Bn, int M, int K, int N) {
  
  NumericVector C((size_t) M * N * Bn);
  
#ifdef HAVE_BATCH_STRIDED
  MKL_INT m = M, n = N, k = K, batch = Bn;
  MKL_INT lda = M, ldb = K, ldc = M;
  MKL_INT strideA = (MKL_INT) M * K;
  MKL_INT strideB = 0;
  MKL_INT strideC = (MKL_INT) M * N;
  
  cblas_dgemm_batch_strided(
    CblasColMajor, CblasNoTrans, CblasNoTrans,
    m, n, k,
    1.0,
    A.begin(), lda, strideA,
    B.begin(), ldb, strideB,
    0.0,
    C.begin(), ldc, strideC,
    batch
  );
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
 * Right argument batched.
 *
 * [K,N,B] is contiguous as [K, N*B] in memory.
 * So: [M,K] x [K, N*B] = [M, N*B] which IS [M,N,B].
 * ONE dgemm call always: no batched BLAS needed!
 * 
 */
// [[Rcpp::export]]
NumericVector bmm_rb(NumericVector A, NumericVector B,
                     int Bn, int M, int K, int N) {
  
  NumericVector C(M * N * Bn);
  
#ifdef HAVE_BATCH_STRIDED
  // MKL: use standard cblas_dgemm
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
              M, N*Bn, K,
              1.0, &A[0], M, &B[0], K, 0.0, &C[0], M);
#else
  // fallback: use our wrapper
  dgemm_nn(M, N * Bn, K,
           &A[0], M,
           &B[0], K,
           &C[0], M);
#endif
           
           C.attr("dim") = IntegerVector::create(M, N, Bn);
           return C;
}


/* --------------------------------------------------------------------------
 * bmm_bb: [M,K,B] x [K,N,B] -> [M,N,B]
 * --------------------------------------------------------------------------
 *
 * Both arguments batched. Each slice [,,b] is contiguous column-major.
 *
 * With HAVE_BATCH_STRIDED: 1 stided dgemm call.
 * Fallback: B dgemm calls on contiguous slices.
 */
// [[Rcpp::export]]
NumericVector bmm_bb(NumericVector A, NumericVector B,
                     int Bn, int M, int K, int N) {
  
  NumericVector C((size_t) M * N * Bn);
  
#ifdef HAVE_BATCH_STRIDED
  MKL_INT m = M, n = N, k = K, batch = Bn;
  MKL_INT lda = M, ldb = K, ldc = M;
  MKL_INT strideA = (MKL_INT) M * K;
  MKL_INT strideB = (MKL_INT) K * N;
  MKL_INT strideC = (MKL_INT) M * N;
  
  cblas_dgemm_batch_strided(
    CblasColMajor, CblasNoTrans, CblasNoTrans,
    m, n, k,
    1.0,
    A.begin(), lda, strideA,
    B.begin(), ldb, strideB,
    0.0,
    C.begin(), ldc, strideC,
    batch
  );
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
