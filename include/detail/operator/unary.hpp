// MIT License
// Copyright 2026 CleverElsie

#ifndef ELSIE_MATRIX_OPERATOR_UNARY_HPP
#define ELSIE_MATRIX_OPERATOR_UNARY_HPP

#include "../_image.hpp"

namespace elsie{

template<class T>
matrix<T>::operator std::string()const{
  std::string ret;
  ret.reserve(get_capacity()*8);
  for(size_t i=0;i<dim().row;++i){
    const std::span<const T> row=(*this)[i];
    for(const auto&x:row)
      ret+=static_cast<std::string>(x)+" ";
    ret+="\n";
  }
  return ret;
}

template<class T>
std::string matrix<T>::string()const{
  return static_cast<std::string>(*this);
}

template<class T>
void matrix<T>::transpose(){
  constexpr size_t block_size=64;
  for(size_t i=0;i<dim().row;i+=block_size){
    const size_t Erow=std::min(i+block_size,dim().row);
    for(size_t j=0;j<dim().col;j+=block_size){
      const size_t Ecol=std::min(j+block_size,dim().col);
      for(size_t I=i;I<Erow;++I)
        for(size_t J=j;J<Ecol;++J)
          (*this)[J,I]=(*this)[I,J];
    }
  }
}

template<class T>
matrix<T> matrix<T>::transpose()const{
  matrix<T> ret(dim().col,dim().row);
  constexpr size_t block_size=64;
  for(size_t i=0;i<dim().row;i+=block_size){
    const size_t Erow=std::min(i+block_size,dim().row);
    for(size_t j=0;j<dim().col;j+=block_size){
      const size_t Ecol=std::min(j+block_size,dim().col);
      for(size_t I=i;I<Erow;++I)
        for(size_t J=j;J<Ecol;++J)
          ret[J,I]=(*this)[I,J];
    }
  }
  return ret;
}

template<class T>
void matrix<T>::negate(){
  for(size_t i=0;i<dim_.row;++i)
    for(size_t j=0;j<dim_.col;++j)
      (*this)[i,j]=-(*this)[i,j];
}

template<class T>
matrix<T> matrix<T>::operator-()const{
  matrix<T> ret(dim());
  for(size_t i=0;i<dim_.row;++i)
    for(size_t j=0;j<dim_.col;++j)
      ret[i,j]=-(*this)[i,j];
  return ret;
}

} // namespace elsie
#endif