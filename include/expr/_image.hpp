#ifndef ELSIE_MATRIX_EXPR_IMAGE_HPP
#define ELSIE_MATRIX_EXPR_IMAGE_HPP
#include <vector>
#include <limits>

#include "../detail/_image.hpp"

namespace elsie{
namespace matrix_mul{
  
  template<class T>
  class half_matrix{
    private:
    size_t n;
    std::vector<T> data;
    static constexpr size_t row_index(size_t i){
      return i*n - (i?(i*(i-1)>>1):0);
    }
    public:
    half_matrix(size_t n): n(n), data(n*(n+1)>>1) {}
    T& operator[](size_t i, size_t j){
      return data[row_index(i)+j];
    }
  };

  template<class T>
  matrix<T> chain(const std::vector<matrix<T>>&matrices){
    constexpr size_t inf=std::numeric_limits<size_t>::max()/3;
    const size_t n=matrices.size();
    half_matrix<size_t> dp(n),s(n);
    for(size_t i=0;i<n;++i){
      // dp[i][i]=0; vectorの初期化で0になるので不要
    }
    for(size_t l=1;l<n;++l){
      for(size_t i=0;i<=n-l;++i){
        const size_t j=i+l-1;
        dp[i,j]=inf;
        for(size_t k=i;k<j;++k){
          size_t cost = dp[i,k] + dp[k+1,j]
            + matrices[i].dim().row * matrices[k].dim().col * matrices[j].dim().col;
          if(cost<dp[i,j]){
            dp[i,j]=cost;
            s[i,j]=k;
          }
        }
      }
    }
  }

}// namespace matrix_mul
}// namespace elsie
#endif