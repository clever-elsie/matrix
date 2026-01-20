#ifndef ELSIE_MATRIX_IO_HPP
#define ELSIE_MATRIX_IO_HPP

#include "./_image.hpp"

namespace elsie{

namespace matrix_io{
  template<class Char, class Traits, class UInt>
  inline bool fast_read_uint(std::basic_istream<Char, Traits>& is, UInt& out)
    requires (std::same_as<UInt, uint32_t> || std::same_as<UInt, uint64_t>)
  {
    using stream_type = std::basic_istream<Char, Traits>;
    using int_type = typename Traits::int_type;
    auto* sb = is.rdbuf();
    if(!sb){
      is.setstate(stream_type::badbit);
      return false;
    }

    int_type ic = sb->sgetc();
    if(Traits::eq_int_type(ic, Traits::eof())){
      is.setstate(stream_type::eofbit | stream_type::failbit);
      return false;
    }

    // whitespace skip (locale 非依存)
    auto is_space = [](Char ch) constexpr {
      return ch==' ' || ch=='\n' || ch=='\t' || ch=='\r' || ch=='\f' || ch=='\v';
    };
    Char c = Traits::to_char_type(ic);
    while(is_space(c)){
      ic = sb->snextc();
      if(Traits::eq_int_type(ic, Traits::eof())){
        is.setstate(stream_type::eofbit | stream_type::failbit);
        return false;
      }
      c = Traits::to_char_type(ic);
    }

    if(c < Char('0') || c > Char('9')){
      is.setstate(stream_type::failbit);
      return false;
    }

    UInt value = 0;
    for(;;){
      value = value * 10 + static_cast<UInt>(c - Char('0'));
      ic = sb->snextc();
      if(Traits::eq_int_type(ic, Traits::eof())){
        is.setstate(stream_type::eofbit);
        break;
      }
      c = Traits::to_char_type(ic);
      if(c < Char('0') || c > Char('9'))
        break;
    }

    out = value;
    return true;
  }
} // namespace matrix_io

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
        if(!matrix_io::fast_read_uint(is, value))[[unlikely]]
          return;
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