// MIT License
// Copyright 2026 CleverElsie

#ifndef ELSIE_MATRIX_OPERATOR_SCALAR_HPP
#define ELSIE_MATRIX_OPERATOR_SCALAR_HPP

#include "../_image.hpp"
#include <type_traits>

namespace elsie{

template<class T>
template<class U>
matrix<T>& matrix<T>::operator*=(const U&rhs){
  for(size_t i=0;i<dim().row;++i)
    for(size_t j=0;j<dim().col;++j)
      (*this)[i,j]*=rhs;
  return*this;
}

template<class T>
template<class U>
matrix<T>& matrix<T>::operator/=(const U&rhs){
  for(size_t i=0;i<dim().row;++i)
    for(size_t j=0;j<dim().col;++j)
      (*this)[i,j]/=rhs;
  return*this;
}

template<class U,class V,class W=std::common_type_t<U,V>>
matrix<W> operator*(const U&lhs,const matrix<V>&rhs){
  matrix<W> ret(rhs.dim());
  for(size_t i=0;i<rhs.dim().row;++i)
    for(size_t j=0;j<rhs.dim().col;++j)
      ret[i,j]=lhs*rhs[i,j];
  return ret;
}

template<class U,class V,class W=std::common_type_t<U,V>>
matrix<W> operator*(const matrix<U>&lhs,const V&rhs){
  matrix<W> ret(lhs.dim());
  for(size_t i=0;i<lhs.dim().row;++i)
    for(size_t j=0;j<lhs.dim().col;++j)
      ret[i,j]=lhs[i,j]*rhs;
  return ret;
}

template<class U,class V,class W=std::common_type_t<U,V>>
matrix<W> operator/(const matrix<U>&lhs,const V&rhs){
  matrix<W> ret(lhs.dim());
  for(size_t i=0;i<lhs.dim().row;++i)
    for(size_t j=0;j<lhs.dim().col;++j)
      ret[i,j]=lhs[i,j]/rhs;
  return ret;
}
} // namespace elsie
#endif