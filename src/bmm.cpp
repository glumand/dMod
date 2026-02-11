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
 *   bmm_lb:  [M,K,B] x [K,N]   -> [M,N,B]   (1 batched dgemm call if available)
 *   bmm_rb:  [M,K]   x [K,N,B] -> [M,N,B]   (1 dgemm call - reshape trick)
 *   bmm_bb:  [M,K,B] x [K,N,B] -> [M,N,B]   (1 batched dgemm call if available)
 *
 * If configure detects cblas_dgemm_batch (MKL, FlexiBLAS, etc.),
 * it defines HAVE_BATCH_GEMM. Then bmm_lb and bmm_bb use a single
 * batched BLAS call.
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


#ifdef HAVE_BATCH_GEMM
/* --------------------------------------------------------------------------
 * CBLAS declarations - no MKL header needed
 * --------------------------------------------------------------------------
 * We declare these ourselves to guarantee LP64 ABI (32-bit integers).
 */
extern "C" {
  
  typedef enum { CblasRowMajor = 101, CblasColMajor = 102 } CBLAS_LAYOUT;
  typedef enum { CblasNoTrans = 111, CblasTrans = 112, CblasConjTrans = 113 } CBLAS_TRANSPOSE;
  
  void cblas_dgemm(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB,
                   const int M, const int N, const int K,
                   const double alpha, const double* A, const int lda,
                   const double* B, const int ldb,
                   const double beta, double* C, const int ldc);
  
  /*
   * cblas_dgemm_batch: Group API for batched GEMM
   * 
   * Performs multiple groups of GEMM operations. Each group shares the same
   * parameters (transa, transb, m, n, k, alpha, beta, lda, ldb, ldc) but
   * operates on different matrices specified by pointer arrays.
   */
  void cblas_dgemm_batch(const CBLAS_LAYOUT Layout,
                         const CBLAS_TRANSPOSE* transa_array,
                         const CBLAS_TRANSPOSE* transb_array,
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
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
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
 * With HAVE_BATCH_GEMM: 1 batched dgemm call with pointer arrays.
 *   - All B pointers point to the same matrix (broadcast).
 * Fallback: Bn dgemm calls on contiguous slices.
 */
// [[Rcpp::export]]
NumericVector bmm_lb(NumericVector A, NumericVector B,
                     int Bn, int M, int K, int N) {
  
  NumericVector C((size_t) M * N * Bn);
  
#ifdef HAVE_BATCH_GEMM
  // Build pointer arrays
  std::vector<const double*> a_ptrs(Bn);
  std::vector<const double*> b_ptrs(Bn);
  std::vector<double*> c_ptrs(Bn);
  
  const int MK = M * K;
  const int MN = M * N;
  const double* A_data = A.begin();
  const double* B_data = B.begin();
  double* C_data = C.begin();
  
  for (int b = 0; b < Bn; ++b) {
    a_ptrs[b] = A_data + b * MK;
    b_ptrs[b] = B_data;  // Same B for all batches (broadcast)
    c_ptrs[b] = C_data + b * MN;
  }
  
  // Single group with Bn operations
  CBLAS_TRANSPOSE transa = CblasNoTrans;
  CBLAS_TRANSPOSE transb = CblasNoTrans;
  double alpha = 1.0;
  double beta = 0.0;
  int lda = M;
  int ldb = K;
  int ldc = M;
  int group_count = 1;
  int group_size = Bn;
  
  cblas_dgemm_batch(CblasColMajor,
                    &transa, &transb,
                    &M, &N, &K,
                    &alpha,
                    a_ptrs.data(), &lda,
                    b_ptrs.data(), &ldb,
                    &beta,
                    c_ptrs.data(), &ldc,
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
 * With HAVE_BATCH_GEMM: 1 batched dgemm call with pointer arrays.
 * Fallback: Bn dgemm calls on contiguous slices.
 */
// [[Rcpp::export]]
NumericVector bmm_bb(NumericVector A, NumericVector B,
                     int Bn, int M, int K, int N) {
  
  NumericVector C((size_t) M * N * Bn);
  
#ifdef HAVE_BATCH_GEMM
  // Build pointer arrays
  std::vector<const double*> a_ptrs(Bn);
  std::vector<const double*> b_ptrs(Bn);
  std::vector<double*> c_ptrs(Bn);
  
  const int MK = M * K;
  const int KN = K * N;
  const int MN = M * N;
  const double* A_data = A.begin();
  const double* B_data = B.begin();
  double* C_data = C.begin();
  
  for (int b = 0; b < Bn; ++b) {
    a_ptrs[b] = A_data + b * MK;
    b_ptrs[b] = B_data + b * KN;
    c_ptrs[b] = C_data + b * MN;
  }
  
  // Single group with Bn operations
  CBLAS_TRANSPOSE transa = CblasNoTrans;
  CBLAS_TRANSPOSE transb = CblasNoTrans;
  double alpha = 1.0;
  double beta = 0.0;
  int lda = M;
  int ldb = K;
  int ldc = M;
  int group_count = 1;
  int group_size = Bn;
  
  cblas_dgemm_batch(CblasColMajor,
                    &transa, &transb,
                    &M, &N, &K,
                    &alpha,
                    a_ptrs.data(), &lda,
                    b_ptrs.data(), &ldb,
                    &beta,
                    c_ptrs.data(), &ldc,
                    group_count, &group_size);
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