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
 *   bmm_lb:  [M,K,B] x [K,N]   -> [M,N,B]   (1 cblas_dgemm_batch call)
 *   bmm_rb:  [M,K]   x [K,N,B] -> [M,N,B]   (1 dgemm call - reshape trick)
 *   bmm_bb:  [M,K,B] x [K,N,B] -> [M,N,B]   (1 cblas_dgemm_batch_strided call)
 *
 * If configure detects batched BLAS (MKL, FlexiBLAS, etc.), it defines
 * HAVE_BATCH_GEMM. Then:
 *   - bmm_bb uses cblas_dgemm_batch_strided (all strides > 0)
 *   - bmm_lb uses cblas_dgemm_batch (pointer arrays, B pointers all
 *     point to the same matrix for broadcast - avoids stride=0 issue)
 *
 * ============================================================================
 */

#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>
#include <vector>

using namespace Rcpp;

#ifndef FCONE
#define FCONE
#endif

/* CBLAS constants (match MKL LP64 ABI - plain integers) */
#define CBLAS_COL_MAJOR  102
#define CBLAS_NO_TRANS   111


#ifdef HAVE_BATCH_GEMM
/* --------------------------------------------------------------------------
 * CBLAS declarations - no MKL header needed
 * --------------------------------------------------------------------------
 * All parameters are plain int to match MKL LP64 ABI exactly.
 * Using C++ enums can cause ABI mismatches (different underlying type).
 */
extern "C" {
  
  void cblas_dgemm(const int Layout,
                   const int TransA, const int TransB,
                   const int M, const int N, const int K,
                   const double alpha,
                   const double* A, const int lda,
                   const double* B, const int ldb,
                   const double beta,
                   double* C, const int ldc);
  
  void cblas_dgemm_batch(const int Layout,
                         const int* transa_array,
                         const int* transb_array,
                         const int* m_array,
                         const int* n_array,
                         const int* k_array,
                         const double* alpha_array,
                         const double** a_array,
                         const int* lda_array,
                         const double** b_array,
                         const int* ldb_array,
                         const double* beta_array,
                         double** c_array,
                         const int* ldc_array,
                         const int group_count,
                         const int* group_size);
  
  void cblas_dgemm_batch_strided(const int Layout,
                                 const int TransA, const int TransB,
                                 const int M, const int N, const int K,
                                 const double alpha,
                                 const double* A, const int lda,
                                 const int stridea,
                                 const double* B, const int ldb,
                                 const int strideb,
                                 const double beta,
                                 double* C, const int ldc,
                                 const int stridec,
                                 const int batch_size);
  
}  // extern "C"
#endif


/* --------------------------------------------------------------------------
 * Query whether batched BLAS is available
 * -------------------------------------------------------------------------- */

// [[Rcpp::export]]
bool has_batch_gemm() {
#ifdef HAVE_BATCH_GEMM
  return true;
#else
  return false;
#endif
}


/* --------------------------------------------------------------------------
 * Standard BLAS wrapper: C = alpha * A * B + beta * C
 * -------------------------------------------------------------------------- */
inline void dgemm_nn(int M, int N, int K,
                     double alpha,
                     const double* A, int lda,
                     const double* B, int ldb,
                     double beta,
                     double* C, int ldc) {
#ifdef HAVE_BATCH_GEMM
  cblas_dgemm(CBLAS_COL_MAJOR, CBLAS_NO_TRANS, CBLAS_NO_TRANS,
              M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
#else
  F77_CALL(dgemm)("N", "N", &M, &N, &K, &alpha,
           A, &lda, B, &ldb, &beta, C, &ldc FCONE FCONE);
#endif
}


/* --------------------------------------------------------------------------
 * bmm_lb: [M,K,B] x [K,N] -> [M,N,B]
 * --------------------------------------------------------------------------
 *
 * Left argument batched: same B multiplied into each batch of A.
 *
 * With HAVE_BATCH_GEMM: 1 cblas_dgemm_batch call with pointer arrays.
 *   All b_ptrs point to the same matrix (broadcast without stride=0).
 * Fallback: Bn dgemm calls on contiguous slices.
 */
// [[Rcpp::export]]
NumericVector bmm_lb(NumericVector A, NumericVector B,
                     int Bn, int M, int K, int N) {
  
  NumericVector C((size_t) M * N * Bn);
  
#ifdef HAVE_BATCH_GEMM
  const double** a_ptrs = (const double**) R_alloc(Bn, sizeof(const double*));
  const double** b_ptrs = (const double**) R_alloc(Bn, sizeof(const double*));
  double**       c_ptrs = (double**)       R_alloc(Bn, sizeof(double*));
  
  const int MK = M * K;
  const int MN = M * N;
  const double* A_data = A.begin();
  const double* B_data = B.begin();
  double* C_data = C.begin();
  
  for (int b = 0; b < Bn; ++b) {
    a_ptrs[b] = A_data + b * MK;
    b_ptrs[b] = B_data;
    c_ptrs[b] = C_data + b * MN;
  }
  
  int transa = CBLAS_NO_TRANS;
  int transb = CBLAS_NO_TRANS;
  double alpha = 1.0;
  double beta = 0.0;
  int lda = M;
  int ldb = K;
  int ldc = M;
  int group_count = 1;
  int group_size = Bn;
  
  cblas_dgemm_batch(CBLAS_COL_MAJOR,
                    &transa, &transb,
                    &M, &N, &K,
                    &alpha,
                    a_ptrs, &lda,
                    b_ptrs, &ldb,
                    &beta,
                    c_ptrs, &ldc,
                    group_count, &group_size);
#else
  const int MK = M * K;
  const int MN = M * N;
  const double* Bptr = B.begin();
  
  for (int b = 0; b < Bn; ++b) {
    dgemm_nn(M, N, K,
             1.0,
             &A[b * MK], M,
             Bptr, K,
             0.0,
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
 */
// [[Rcpp::export]]
NumericVector bmm_rb(NumericVector A, NumericVector B,
                     int Bn, int M, int K, int N) {
  
  NumericVector C(M * N * Bn);
  
  dgemm_nn(M, N * Bn, K,
           1.0,
           A.begin(), M,
           B.begin(), K,
           0.0,
           C.begin(), M);
  
  C.attr("dim") = IntegerVector::create(M, N, Bn);
  return C;
}


/* --------------------------------------------------------------------------
 * bmm_bb: [M,K,B] x [K,N,B] -> [M,N,B]
 * --------------------------------------------------------------------------
 *
 * Both arguments batched. Each slice [,,b] is contiguous column-major.
 *
 * With HAVE_BATCH_GEMM: 1 cblas_dgemm_batch_strided call.
 * Fallback: Bn dgemm calls on contiguous slices.
 */
// [[Rcpp::export]]
NumericVector bmm_bb(NumericVector A, NumericVector B,
                     int Bn, int M, int K, int N) {
  
  NumericVector C((size_t) M * N * Bn);
  
#ifdef HAVE_BATCH_GEMM
  cblas_dgemm_batch_strided(CBLAS_COL_MAJOR, CBLAS_NO_TRANS, CBLAS_NO_TRANS,
                            M, N, K,
                            1.0,
                            A.begin(), M, M * K,
                            B.begin(), K, K * N,
                            0.0,
                            C.begin(), M, M * N,
                            Bn);
#else
  const int MK = M * K;
  const int KN = K * N;
  const int MN = M * N;
  
  for (int b = 0; b < Bn; ++b) {
    dgemm_nn(M, N, K,
             1.0,
             &A[b * MK], M,
             &B[b * KN], K,
             0.0,
             &C[b * MN], M);
  }
#endif
  
  C.attr("dim") = IntegerVector::create(M, N, Bn);
  return C;
}