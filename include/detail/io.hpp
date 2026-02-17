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
      if constexpr(matrix_io::Writable<T>){
        auto [ptr, len] = element.write(element_buf);
        os.rdbuf()->sputn(ptr, len);
      }else if constexpr(matrix_io::ostreamable<Char, Traits, T>)
        os<<element;
      else if constexpr(std::convertible_to<T,std::string>){
        auto element_str = static_cast<std::string>(element);
        os.rdbuf()->sputn(element_str.data(), element_str.size());
      }else static_assert(false, "T is not able to be outputted to ostream.");
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
void matrix<T>::read(std::basic_istream<Char, Traits>& is) requires (!matrix_io::modint_v<T>){
  typename std::basic_istream<Char, Traits>::sentry sent(is);
  if(!sent)[[unlikely]] return;
  const auto [row_sz,col_sz]=dim();
  for(size_t i=0;i<row_sz;++i){
    auto row = (*this)[i];
    for(size_t j=0;j<col_sz;++j)
      if constexpr(!matrix_io::istreamable<Char, Traits, T>)
        static_assert(false, "Tis not able to be inputted from istream.");
      else is>>row[j];
  }
}

template<class T>
template<bool validated,class Char, class Traits>
void matrix<T>::read(std::basic_istream<Char, Traits>& is) requires (matrix_io::modint_v<T>){
  typename std::basic_istream<Char, Traits>::sentry sent(is);
  if(!sent)[[unlikely]] return;
  const auto [row_sz,col_sz]=dim();
  for(size_t i=0;i<row_sz;++i){
    auto row = (*this)[i];
    for(size_t j=0;j<col_sz;++j){
      typename T::value_type value;
      if constexpr(std::same_as<typename T::value_type,uint32_t> || std::same_as<typename T::value_type,uint64_t>){
        static_assert(std::same_as<Char, char>, "fast modint input uses std::from_chars, which requires Char=char.");
        auto* sb = is.rdbuf();
        char c;
        // whitespace skip (locale 非依存)
        do{
          const auto ic = sb->sbumpc();
          if(Traits::eq_int_type(ic, Traits::eof())) return;
          c = Traits::to_char_type(ic);
        }while(c==' ' || c=='\n' || c=='\t' || c=='\r' || c=='\f' || c=='\v');

        char buf[32];
        size_t len = 0;
        for(;;){
          buf[len++] = c;
          const auto ic = sb->sbumpc();
          if(Traits::eq_int_type(ic, Traits::eof())) break;
          c = Traits::to_char_type(ic);
          if(c<'0' || c>'9') break;
        }
        (void)std::from_chars(buf, buf + len, value, 10);
      }else{
        is>>value;
      }
      if constexpr(validated){
        if constexpr(matrix_io::raw_constructible<T>){
          row[j]=T::raw(value);
        }else row[j]=T(value);
      }else row[j]=T(value);
    }
  }
}

} // namespace elsie

#endif