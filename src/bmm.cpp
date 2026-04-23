/*
 * ============================================================================
 * Batched Matrix Multiplication (batch-first convention)
 * ============================================================================
 *
 * All 3D arrays use [B, M, N] convention (batch index FIRST) matching
 * CppODE's sensitivity output. In column-major R storage the batch axis
 * varies fastest, so per-batch slices are NOT contiguous. Only bmm_lb
 * gets a single-dgemm fast path (via the [B,M,K] == [B*M, K] memory
 * identity); bmm_rb and bmm_bb loop per batch with scatter-gather copies
 * into contiguous scratch buffers. The inner dgemm goes through R's
 * linked BLAS.
 *
 * OPERATIONS:
 *
 *   bmm_lb:  [B,M,K] x [K,N]   -> [B,M,N]   (1 dgemm: [B*M,K] x [K,N])
 *   bmm_rb:  [M,K]   x [B,K,N] -> [B,M,N]   (Bn dgemm calls on scratch)
 *   bmm_bb:  [B,M,K] x [B,K,N] -> [B,M,N]   (Bn dgemm calls on scratch)
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


/* --------------------------------------------------------------------------
 * Standard BLAS wrapper: C = alpha * A * B + beta * C (column-major)
 * -------------------------------------------------------------------------- */
inline void dgemm_nn(int M, int N, int K,
                     double alpha,
                     const double* A, int lda,
                     const double* B, int ldb,
                     double beta,
                     double* C, int ldc) {
  F77_CALL(dgemm)("N", "N", &M, &N, &K, &alpha,
           A, &lda, B, &ldb, &beta, C, &ldc FCONE FCONE);
}


/* --------------------------------------------------------------------------
 * bmm_lb: [B,M,K] x [K,N] -> [B,M,N]
 * --------------------------------------------------------------------------
 *
 * In column-major storage [B,M,K] and [B*M, K] share memory. One dgemm on
 * the reshaped left operand computes all batches at once with no copies.
 */
// [[Rcpp::export]]
NumericVector bmm_lb(NumericVector A, NumericVector B,
                     int Bn, int M, int K, int N) {

  NumericVector C((size_t) Bn * M * N);

  dgemm_nn(Bn * M, N, K,
           1.0,
           A.begin(), Bn * M,
           B.begin(), K,
           0.0,
           C.begin(), Bn * M);

  C.attr("dim") = IntegerVector::create(Bn, M, N);
  return C;
}


/* --------------------------------------------------------------------------
 * bmm_rb: [M,K] x [B,K,N] -> [B,M,N]
 * --------------------------------------------------------------------------
 *
 * Right argument batched under batch-first: slice B[b,,] has row stride Bn
 * (not 1), so dgemm cannot read it directly. Gather each slice into a
 * contiguous [K,N] scratch, dgemm against the shared A, scatter the [M,N]
 * result into C[b,,].
 */
// [[Rcpp::export]]
NumericVector bmm_rb(NumericVector A, NumericVector B,
                     int Bn, int M, int K, int N) {

  NumericVector C((size_t) Bn * M * N);

  std::vector<double> Bbuf((size_t) K * N);
  std::vector<double> Cbuf((size_t) M * N);

  const double* Bdata = B.begin();
  double* Cdata = C.begin();

  for (int b = 0; b < Bn; ++b) {
    // Gather B[b,,]: element (k,n) at Bdata[b + k*Bn + n*Bn*K]
    for (int n = 0; n < N; ++n) {
      for (int k = 0; k < K; ++k) {
        Bbuf[k + n * K] = Bdata[b + k * Bn + n * (size_t) Bn * K];
      }
    }

    dgemm_nn(M, N, K,
             1.0,
             A.begin(), M,
             Bbuf.data(), K,
             0.0,
             Cbuf.data(), M);

    // Scatter into C[b,,]: element (m,n) at Cdata[b + m*Bn + n*Bn*M]
    for (int n = 0; n < N; ++n) {
      for (int m = 0; m < M; ++m) {
        Cdata[b + m * Bn + n * (size_t) Bn * M] = Cbuf[m + n * M];
      }
    }
  }

  C.attr("dim") = IntegerVector::create(Bn, M, N);
  return C;
}


/* --------------------------------------------------------------------------
 * bmm_bb: [B,M,K] x [B,K,N] -> [B,M,N]
 * --------------------------------------------------------------------------
 *
 * Both operands batched. Per slice: gather A[b,,] and B[b,,] into
 * contiguous scratch, dgemm, scatter into C[b,,].
 */
// [[Rcpp::export]]
NumericVector bmm_bb(NumericVector A, NumericVector B,
                     int Bn, int M, int K, int N) {

  NumericVector C((size_t) Bn * M * N);

  std::vector<double> Abuf((size_t) M * K);
  std::vector<double> Bbuf((size_t) K * N);
  std::vector<double> Cbuf((size_t) M * N);

  const double* Adata = A.begin();
  const double* Bdata = B.begin();
  double* Cdata = C.begin();

  for (int b = 0; b < Bn; ++b) {
    // Gather A[b,,]
    for (int k = 0; k < K; ++k) {
      for (int m = 0; m < M; ++m) {
        Abuf[m + k * M] = Adata[b + m * Bn + k * (size_t) Bn * M];
      }
    }
    // Gather B[b,,]
    for (int n = 0; n < N; ++n) {
      for (int k = 0; k < K; ++k) {
        Bbuf[k + n * K] = Bdata[b + k * Bn + n * (size_t) Bn * K];
      }
    }

    dgemm_nn(M, N, K,
             1.0,
             Abuf.data(), M,
             Bbuf.data(), K,
             0.0,
             Cbuf.data(), M);

    // Scatter into C[b,,]
    for (int n = 0; n < N; ++n) {
      for (int m = 0; m < M; ++m) {
        Cdata[b + m * Bn + n * (size_t) Bn * M] = Cbuf[m + n * M];
      }
    }
  }

  C.attr("dim") = IntegerVector::create(Bn, M, N);
  return C;
}
