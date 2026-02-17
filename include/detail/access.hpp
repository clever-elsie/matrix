// MIT License
// Copyright 2026 CleverElsie

#ifndef ELSIE_MATRIX_ACCESS_HPP
#define ELSIE_MATRIX_ACCESS_HPP

#include "./_image.hpp"

namespace elsie{

template<class T>
matrix<T>::dimension matrix<T>::dim()const{
  return dim_;
}

template<class T>
matrix<T>::dimension matrix<T>::leading_dim()const{
  return leading_;
}

template<class T>
T& matrix<T>::operator[](size_t i, size_t j){
  return data_[i*leading_.col+j];
}

template<class T>
const T& matrix<T>::operator[](size_t i, size_t j)const{
  return data_[i*leading_.col+j];
}

template<class T>
T matrix<T>::val(size_t i, size_t j)const{
  return data_[i*leading_.col+j];
}

template<class T>
std::span<T> matrix<T>::operator[](size_t i){
  return {data_+i*leading_.col,data_+i*leading_.col+leading_.col};
}

template<class T>
std::span<const T> matrix<T>::operator[](size_t i)const{
  return {data_+i*leading_.col,data_+i*leading_.col+leading_.col};
}

} // namespace elsie
#endif