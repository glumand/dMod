#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>
using namespace Rcpp;

/*
 * ============================================================================
 * Batched Matrix Multiplication (BLAS-optimized)
 * ============================================================================
 *
 * All functions use R's [B,M,K] convention (batch index FIRST) and exploit
 * memory layout for maximum BLAS efficiency.
 *
 * MEMORY LAYOUT INSIGHT:
 *
 *   R's [B,M,K] array stores element [b,m,k] at position: b + B*m + B*M*k
 *
 *   Key observations:
 *   1. [B,M,K] IS contiguous as [B*M, K] matrix  -> bmm_lb uses ONE dgemm!
 *   2. [B,K,N] slice [,,n] is contiguous [B,K]   -> bmm_rb uses N dgemms
 *   3. [B,M,K] slice [b,,] is NOT standard column-major -> bmm_bb copies
 *
 * OPERATIONS:
 *
 *   bmm_lb:  [B,M,K] x [K,N]   -> [B,M,N]   (1 dgemm call!)
 *   bmm_rb:  [M,K]   x [B,K,N] -> [B,M,N]   (N dgemm calls)
 *   bmm_bb:  [B,M,K] x [B,K,N] -> [B,M,N]   (B dgemm calls with copies)
 *
 * ============================================================================
 */


/* --------------------------------------------------------------------------
 * BLAS wrappers
 * -------------------------------------------------------------------------- */

// C = A * B
inline void dgemm_nn(int M, int N, int K,
                     const double* A, int lda,
                     const double* B, int ldb,
                     double* C, int ldc) {
  const double one = 1.0, zero = 0.0;
  F77_CALL(dgemm)("N", "N", &M, &N, &K, &one,
           A, &lda, B, &ldb, &zero, C, &ldc FCONE FCONE);
}

// C = A * B^T
inline void dgemm_nt(int M, int N, int K,
                     const double* A, int lda,
                     const double* B, int ldb,
                     double* C, int ldc) {
  const double one = 1.0, zero = 0.0;
  F77_CALL(dgemm)("N", "T", &M, &N, &K, &one,
           A, &lda, B, &ldb, &zero, C, &ldc FCONE FCONE);
}


/* --------------------------------------------------------------------------
 * bmm_lb: [B,M,K] x [K,N] -> [B,M,N]
 * --------------------------------------------------------------------------
 *
 * ONE dgemm call!
 *
 * Exploits: [B,M,K] is contiguous as [B*M, K] in memory.
 * Compute:  [B*M, K] x [K, N] = [B*M, N] which IS [B,M,N].
 */
// [[Rcpp::export]]
NumericVector bmm_lb(NumericVector A, NumericVector B,
                     int Bn, int M, int K, int N) {
  
  NumericVector C(Bn * M * N);
  const int BM = Bn * M;
  
  dgemm_nn(BM, N, K,
           &A[0], BM,    // [B,M,K] as [B*M, K]
           &B[0], K,
           &C[0], BM);   // [B*M, N] = [B,M,N]
           
           C.attr("dim") = IntegerVector::create(Bn, M, N);
           return C;
}


/* --------------------------------------------------------------------------
 * bmm_rb: [M,K] x [B,K,N] -> [B,M,N]
 * --------------------------------------------------------------------------
 *
 * N dgemm calls, each on contiguous data.
 *
 * For [B,K,N]: slice B[,,n] at offset B*K*n is contiguous [B,K].
 *
 * Want: C[b,m,n] = sum_k A[m,k] * B[b,k,n]
 *     = sum_k B[b,k,n] * A^T[k,m]
 *
 * So: C[,,n] = B[,,n] * A^T  with shapes [B,K] x [K,M] = [B,M]
 */
// [[Rcpp::export]]
NumericVector bmm_rb(NumericVector A, NumericVector B,
                     int Bn, int M, int K, int N) {
  
  NumericVector C(Bn * M * N);
  const int BK = Bn * K;
  const int BM = Bn * M;
  
  for (int n = 0; n < N; ++n) {
    dgemm_nt(Bn, M, K,
             &B[n * BK], Bn,   // B[,,n] as [B,K]
             &A[0], M,         // A^T via 'T' flag
             &C[n * BM], Bn);  // C[,,n] as [B,M]
  }
  
  C.attr("dim") = IntegerVector::create(Bn, M, N);
  return C;
}


/* --------------------------------------------------------------------------
 * bmm_bb: [B,M,K] x [B,K,N] -> [B,M,N]
 * --------------------------------------------------------------------------
 *
 * R's [B,M,K] layout is incompatible with strided BLAS for [b,,] slices:
 * - Row stride = B, Column stride = B*M
 * - But dgemm requires col_stride = lda * nrows (standard column-major)
 * 
 * Solution: Copy slices to contiguous buffers per batch.
 * For typical sizes (M,K,N < 100), copy overhead is negligible.
 */
// [[Rcpp::export]]
NumericVector bmm_bb(NumericVector A, NumericVector B,
                     int Bn, int M, int K, int N) {
  
  NumericVector C(Bn * M * N);
  std::vector<double> Aslice(M * K);
  std::vector<double> Bslice(K * N);
  std::vector<double> Cslice(M * N);
  
  for (int b = 0; b < Bn; ++b) {
    
    // Copy A[b,,] to contiguous Aslice[M,K] (column-major)
    for (int k = 0; k < K; ++k)
      for (int m = 0; m < M; ++m)
        Aslice[m + M*k] = A[b + Bn*m + Bn*M*k];
    
    // Copy B[b,,] to contiguous Bslice[K,N] (column-major)
    for (int n = 0; n < N; ++n)
      for (int k = 0; k < K; ++k)
        Bslice[k + K*n] = B[b + Bn*k + Bn*K*n];
    
    // Cslice[M,N] = Aslice[M,K] * Bslice[K,N]
    dgemm_nn(M, N, K,
             &Aslice[0], M,
             &Bslice[0], K,
             &Cslice[0], M);
             
             // Copy Cslice back to C[b,,]
             for (int n = 0; n < N; ++n)
               for (int m = 0; m < M; ++m)
                 C[b + Bn*m + Bn*M*n] = Cslice[m + M*n];
  }
  
  C.attr("dim") = IntegerVector::create(Bn, M, N);
  return C;
}