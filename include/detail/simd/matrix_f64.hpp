// MIT License
// Copyright 2026 CleverElsie

// AVX2 / AVX512 を用いた double 行列の高速 FMA 実装
// C += A * B をブロック分割 + SIMD で計算する。
//   ブロック: 64x64x64（キャッシュ効率を狙った素直な分割）
//   AVX512: 8 要素 (512bit) を一括で __m512d で処理
//   AVX2  : 4 要素 (256bit) を __m256d で処理
// いずれも数式は同じ: c[j..] = a[i][k] * b[k][j..] + c[j..] （FMA）
#pragma once

#include <immintrin.h>
#include <algorithm>
#include "../_image.hpp"

namespace elsie {

using dmat_t = elsie::matrix<double>;

#ifdef __AVX512F__
// -----------------------------------------------------------------------------
// AVX512 実装（8 double 同時）
// -----------------------------------------------------------------------------
template<>
template<>
__attribute__((target("avx512f")))
__attribute__((optimize("O3")))
inline dmat_t& dmat_t::fma_impl_blocked(const dmat_t& a, const dmat_t& b) {
  [[assume(dim().row == a.dim().row)]];
  [[assume(dim().col == b.dim().col)]];
  [[assume(a.dim().col == b.dim().row)]];

  constexpr size_t BLOCK_SIZE = 64;
  const auto [row, col] = dim();
  const size_t midCR = a.dim().col;

  size_t I, J, K, i, j, k, it, jt, kt;
  constexpr size_t simd = 8; // 8 doubles

  for (I = 0, it = std::min(BLOCK_SIZE, row); I < row; I += BLOCK_SIZE, it = std::min(it + BLOCK_SIZE, row)) {
    for (J = 0, jt = std::min(BLOCK_SIZE, col); J < col; J += BLOCK_SIZE, jt = std::min(jt + BLOCK_SIZE, col)) {
      for (K = 0, kt = std::min(BLOCK_SIZE, midCR); K < midCR; K += BLOCK_SIZE, kt = std::min(kt + BLOCK_SIZE, midCR)) {
        for (i = I; i < it; ++i) {
          const auto as = a[i];
          auto cs = (*this)[i];
          for (k = K; k < kt; ++k) {
            const double ak = as[k];
            const auto b_row = b[k];

            __m512d a8 = _mm512_set1_pd(ak);
            for (j = J; j + simd - 1 < jt; j += simd) {
              __m512d b8 = _mm512_loadu_pd(&b_row[j]);
              __m512d c8 = _mm512_loadu_pd(&cs[j]);
              c8 = _mm512_fmadd_pd(a8, b8, c8); // c += a * b
              _mm512_storeu_pd(&cs[j], c8);
            }
            for (; j < jt; ++j) cs[j] += ak * b_row[j];
          }
        }
      }
    }
  }
  return *this;
}
#elif defined(__AVX2__)
// -----------------------------------------------------------------------------
// AVX2 実装（4 double 同時）
// -----------------------------------------------------------------------------
template<>
template<>
__attribute__((target("avx2")))
__attribute__((optimize("O3")))
inline dmat_t& dmat_t::fma(const dmat_t& a, const dmat_t& b) {
  [[assume(dim().row == a.dim().row)]];
  [[assume(dim().col == b.dim().col)]];
  [[assume(a.dim().col == b.dim().row)]];

  constexpr size_t BLOCK_SIZE = 64;
  const auto [row, col] = dim();
  const size_t midCR = a.dim().col;

  size_t I, J, K, i, j, k, it, jt, kt;
  constexpr size_t simd = 4; // 4 doubles

  for (I = 0, it = std::min(BLOCK_SIZE, row); I < row; I += BLOCK_SIZE, it = std::min(it + BLOCK_SIZE, row)) {
    for (J = 0, jt = std::min(BLOCK_SIZE, col); J < col; J += BLOCK_SIZE, jt = std::min(jt + BLOCK_SIZE, col)) {
      for (K = 0, kt = std::min(BLOCK_SIZE, midCR); K < midCR; K += BLOCK_SIZE, kt = std::min(kt + BLOCK_SIZE, midCR)) {
        for (i = I; i < it; ++i) {
          const auto as = a[i];
          auto cs = (*this)[i];
          for (k = K; k < kt; ++k) {
            const double ak = as[k];
            const auto b_row = b[k];

            __m256d a4 = _mm256_set1_pd(ak);
            for (j = J; j + simd - 1 < jt; j += simd) {
              __m256d b4 = _mm256_loadu_pd(&b_row[j]);
              __m256d c4 = _mm256_loadu_pd(&cs[j]);
              c4 = _mm256_fmadd_pd(a4, b4, c4);
              _mm256_storeu_pd(&cs[j], c4);
            }
            for (; j < jt; ++j) cs[j] += ak * b_row[j];
          }
        }
      }
    }
  }
  return *this;
}
#endif

} // namespace elsie
