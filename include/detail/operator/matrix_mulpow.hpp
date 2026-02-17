// MIT License
// Copyright 2026 CleverElsie

#ifndef ELSIE_MATRIX_OPERATOR_MATRIX_MULPOW_HPP
#define ELSIE_MATRIX_OPERATOR_MATRIX_MULPOW_HPP

#include "../_image.hpp"
#include <cassert>

namespace elsie{

template<class T>
template<class U, class V>
matrix<T>& matrix<T>::fma_impl_naive(const matrix<U>&a,const matrix<V>&b){
  [[assume(dim().row==a.dim().row)]];
  [[assume(dim().col==b.dim().col)]];
  [[assume(a.dim().col==b.dim().row)]];
  const auto[row,col]=dim();
  const size_t midCR=a.dim().col;
  for(size_t i=0;i<row;++i){
    const auto as=a[i];
    auto cs = (*this)[i];
    for(size_t k=0;k<midCR;++k){
      const auto&ak=as[k];
      const auto bs=b[k];
      for(size_t j=0;j<col;++j)
        cs[j]+=ak*bs[j];
    }
  }
  return *this;
}

template<class T>
template<class U, class V>
matrix<T>& matrix<T>::fma_impl_blocked(const matrix<U>&a,const matrix<V>&b){
  [[assume(dim().row==a.dim().row)]];
  [[assume(dim().col==b.dim().col)]];
  [[assume(a.dim().col==b.dim().row)]];
  constexpr const size_t bs=64; // block size
  const auto[row,col]=dim();
  const size_t midCR=a.dim().col;
  size_t I,J,K,i,j,k,it,jt,kt;
  for(I=0,it=std::min(bs,row);I<row;I+=bs,it=std::min(it+bs,row))
  for(J=0,jt=std::min(bs,col);J<col;J+=bs,jt=std::min(jt+bs,col))
  for(K=0,kt=std::min(bs,midCR);K<midCR;K+=bs,kt=std::min(kt+bs,midCR))
    for(i=I;i<it;++i){
      const auto as=a[i];
      auto cs=(*this)[i];
      for(k=K;k<kt;++k){
        const auto&ak=as[k];
        const auto bs=b[k];
        for(j=J;j<jt;++j)
          cs[j]+=ak*bs[j];
      }
    }
  return *this;
}

template<class T>
template<class U, class V>
matrix<T>& matrix<T>::fma_impl_strassen(const matrix<U>&a,const matrix<V>&b){
  [[assume(dim().row==a.dim().row)]];
  [[assume(dim().col==b.dim().col)]];
  [[assume(a.dim().col==b.dim().row)]];
  const auto[row,col]=dim();
  const size_t midCR=a.dim().col;
  std::array<matrix<U>,4>as=a.split();
  std::array<matrix<V>,4>bs=b.split();
  matrix<U> s2=as[0]+as[1], s9=as[0]-as[2], s5=as[0]+as[3], s3=as[2]+as[3], s7=as[1]-as[3];
  matrix<V> s1=bs[1]-bs[3], s8=bs[2]+bs[3], s6=bs[0]+bs[3], s4=bs[2]-bs[0], sA=bs[0]+bs[1];
  matrix<T> p1=as[0]*s1;
  matrix<T> p2=s2*bs[3];
  matrix<T> p3=s3*bs[0];
  matrix<T> p4=as[3]*s4;
  matrix<T> p5=s5*s6;
  matrix<T> p6=s7*s8;
  matrix<T> p7=s9*sA;

  matrix<T>::merge(*this, 0, 0,
    ((p5+p4)-=p2)+=p6,p1+p2,
    p3+p4, ((p5+p1)-=p3)-=p7
  );
  return *this;
}

template<class T>
template<class U,class V>
matrix<T>& matrix<T>::fma(const matrix<U>&a,const matrix<V>&b){
  //assert((dim().row==a.dim().row) && (dim().col==b.dim().col) && (a.dim().col==b.dim().row));
  const auto[row,col]=dim();
  const size_t midCR=a.dim().col;
  constexpr const size_t block_size=64;
  constexpr const size_t strassen_threshold=256;
  const size_t min=std::min<size_t>({row,col,midCR});
  if(min<block_size) return fma_impl_naive(a,b);
  else return fma_impl_blocked(a,b);
  /*
  else if(min<strassen_threshold) return fma_impl_blocked(a,b);
  else return fma_impl_strassen(a,b);
  */
}

template<class T>
template<class U>
matrix<T>& matrix<T>::operator*=(const matrix<U>&rhs){
  [[assume(dim().col==rhs.dim().row)]];
  return *this=*this*rhs;
}

template<class U,class V,class W=std::common_type_t<U,V>>
matrix<W> operator*(const matrix<U>&lhs,const matrix<V>&rhs){
  [[assume(lhs.dim().col==rhs.dim().row)]];
  matrix<W> ret(lhs.dim().row,rhs.dim().col,W());
  ret.fma(lhs,rhs);
  return ret;
}

template<class T>
matrix<T> matrix<T>::pow_impl(matrix<T>&&k, uint64_t b){
  [[assume(dim().row==dim().col)]];
  matrix<T> ret(k.dim(),T());
  for(size_t i=0;i<dim().row;++i)
    ret[i,i]=T(1);
  for(;b;b>>=1){
    if(b&1) ret=ret*k;
    k=k*k;
  }
}

template<class T>
void matrix<T>::pow(uint64_t b){
  *this=pow_impl(std::move(*this),b);
}

template<class T>
matrix<T> matrix<T>::pow(uint64_t b)const{
  return pow_impl(matrix<T>(this->copy()),b);
}

} // namespace elsie
#endif