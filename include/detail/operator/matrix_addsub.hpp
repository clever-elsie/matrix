// MIT License
// Copyright 2026 CleverElsie

#ifndef ELSIE_MATRIX_OPERATOR_MATRIX_ADDSUB_HPP
#define ELSIE_MATRIX_OPERATOR_MATRIX_ADDSUB_HPP

#include "../_image.hpp"

namespace elsie{

template<class T>
template<class U>
matrix<T>& matrix<T>::operator+=(const matrix<U>&rhs){
  [[assume(dim().row==rhs.dim().row)]];
  [[assume(dim().col==rhs.dim().col)]];
  for(size_t i=0;i<dim().row;++i)
    for(size_t j=0;j<dim().col;++j)
      (*this)[i,j]+=rhs[i,j];
  return*this;
}

template<class T>
template<class U>
matrix<T>& matrix<T>::operator-=(const matrix<U>&rhs){
  [[assume(dim().row==rhs.dim().row)]];
  [[assume(dim().col==rhs.dim().col)]];
  for(size_t i=0;i<dim().row;++i)
    for(size_t j=0;j<dim().col;++j)
      (*this)[i,j]-=rhs[i,j];
  return*this;
}

template<class U,class V,class W=std::common_type_t<U,V>>
matrix<W> operator+(const matrix<U>&lhs,const matrix<V>&rhs){
  [[assume(lhs.dim().row==rhs.dim().row)]];
  [[assume(lhs.dim().col==rhs.dim().col)]];
  matrix<W> ret(lhs.dim());
  for(size_t i=0;i<lhs.dim().row;++i)
    for(size_t j=0;j<lhs.dim().col;++j)
      ret[i,j]=lhs[i,j]+rhs[i,j];
  return ret;
}

template<class U,class V,class W=std::common_type_t<U,V>>
matrix<W> operator-(const matrix<U>&lhs,const matrix<V>&rhs){
  [[assume(lhs.dim().row==rhs.dim().row)]];
  [[assume(lhs.dim().col==rhs.dim().col)]];
  matrix<W> ret(lhs.dim());
  for(size_t i=0;i<lhs.dim().row;++i)
    for(size_t j=0;j<lhs.dim().col;++j)
      ret[i,j]=lhs[i,j]-rhs[i,j];
  return ret;
}

} // namespace elsie
#endif