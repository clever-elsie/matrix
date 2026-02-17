// MIT License
// Copyright 2026 CleverElsie

#ifndef ELSIE_MATRIX_IO_HPP
#define ELSIE_MATRIX_IO_HPP

#include "./_image.hpp"
#include <charconv>

namespace elsie{

template<class Char, class Traits, class T>
std::basic_ostream<Char, Traits>& operator<<(std::basic_ostream<Char, Traits>& os, const matrix<T>& mat) {
  char element_buf[128];
  typename std::basic_ostream<Char, Traits>::sentry sent(os);
  if(!sent)[[unlikely]] return os;
  const auto [row_sz,col_sz]=mat.dim();
  for(size_t i=0;i<row_sz;++i){
    const auto& row = mat[i];
    for(size_t j=0;j<col_sz;++j){
      const auto& element = row[j];
      os<<element;
      if(j+1<row.size())[[likely]]
        os.rdbuf()->sputc(' ');
    }
    if(i+1<mat.dim().row)[[likely]]
      os.rdbuf()->sputc('\n');
  }
  return os;
}

template<class T>
template<class Char, class Traits>
void matrix<T>::read(std::basic_istream<Char, Traits>& is){
  typename std::basic_istream<Char, Traits>::sentry sent(is);
  if(!sent)[[unlikely]] return;
  const auto [row_sz,col_sz]=dim();
  for(size_t i=0;i<row_sz;++i){
    auto row = (*this)[i];
    for(size_t j=0;j<col_sz;++j)
      is>>row[j];
  }
}

} // namespace elsie

#endif