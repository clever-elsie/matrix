// MIT License
// Copyright 2026 CleverElsie

#ifndef ELSIE_MATRIX_IMAGE_HPP
#define ELSIE_MATRIX_IMAGE_HPP
#include <cstddef>
#include <cstdint>
#include <span>
#include <array>
#include <memory>
#include <string>
#include <new>
#include <concepts>
#include <type_traits>
#include <iostream>

namespace elsie{

namespace matrix_io{
  template<class T>
  concept Writable=requires(const T&x, char (&buf)[128]){
    {x.write(buf)}->std::convertible_to<std::pair<char*,size_t>>;
  };
  template<class T>
  concept modint_v=requires{
    requires std::integral<typename T::value_type>;
    { T::is_modint }->std::convertible_to<bool>;
  } && T::is_modint;
  template<class T>
  concept raw_constructible=requires{
    { T::raw(0) }->std::convertible_to<T>;
    { T::raw(0LL) }->std::convertible_to<T>;
  };
  template<class Char, class Traits, class T>
  concept istreamable=requires(std::basic_istream<Char, Traits>& is, T&x){
    {is>>x}->std::convertible_to<std::basic_istream<Char, Traits>&>;
  };
  template<class Char, class Traits, class T>
  concept ostreamable=requires(std::basic_ostream<Char, Traits>& os, const T&x){
    {os<<x}->std::convertible_to<std::basic_ostream<Char, Traits>&>;
  };
};

// コピーコンストラクタ，コピー代入はviewの作成
// moveは所有権を持っていればmove，viewの場合は単にコピーしてコピー元のviewを無効化
// deep copyは.copy()を使う
// 所有権の有無はプログラマに委ねる．assert用にowned()を用意する
template<class T>
class matrix{
  public:
  struct dimension{
    size_t row,col;
    dimension():row(0),col(0){}
    dimension(size_t row, size_t col):row(row),col(col){}
    dimension(const dimension& other):row(other.row),col(other.col){}
    dimension& operator=(const dimension& other){row=other.row; col=other.col; return *this;}
    size_t capacity()const{return row*col;}
    bool operator==(const dimension& other)const{return row==other.row&&col==other.col;}
    bool operator!=(const dimension& other)const{return !(*this==other);}
  };

  private:
  dimension dim_,leading_;
  // capacityは本当のオーバーフローをしないために，確保している配列の末尾までの距離を持つ
  size_t capacity_data; // { capacity:63, is_view:1 }
  T* data_;

  private:
  static size_t calc_capacity_data(size_t capacity, bool is_view=false){
    return (capacity<<1) | (is_view?1:0);
  }
  size_t get_capacity()const{ return capacity_data>>1; }
  bool is_view()const{ return capacity_data&1; }

  // make_view用のコンストラクタ
  matrix(const matrix<T>& other, size_t row, size_t col, size_t row_offset, size_t col_offset);
  public:
  // constructor.destructor.operator=.hpp
  matrix();
  matrix(const dimension&dim);
  matrix(const dimension&dim, const T& init);
  matrix(const size_t row, const size_t col);
  matrix(const size_t row, const size_t col, const T& init);
  matrix(const matrix& other);
  matrix(matrix&& other);
  matrix& operator=(const matrix& other);
  matrix& operator=(matrix&& other);
  ~matrix();
  matrix<T> make_view(size_t row, size_t col, size_t row_offset, size_t col_offset)const;
  matrix<T> copy()const;
  matrix<T> copy(size_t row, size_t col, size_t row_offset, size_t col_offset)const;
  private:
  static T* allocate(size_t n){
    if constexpr (alignof(T)<=__STDCPP_DEFAULT_NEW_ALIGNMENT__)
      return static_cast<T*>(::operator new(n*sizeof(T)));
    else return static_cast<T*>(::operator new(n*sizeof(T), std::align_val_t(alignof(T))));
  }
  static void default_construct(T*p,size_t n,T init=T()){
    if constexpr(false==std::is_trivially_default_constructible_v<T>)
      for(size_t i=0;i<n;++i)
        ::new (static_cast<void*>(p+i)) T(init);
  }
  static void free(T*p){
    if constexpr (alignof(T)<=__STDCPP_DEFAULT_NEW_ALIGNMENT__)
      ::operator delete(p);
    else ::operator delete(p, std::align_val_t(alignof(T)));
  }
  public:

  // access.hpp
  dimension dim()const;
  dimension leading_dim()const;
  T& operator[](size_t i, size_t j);
  const T& operator[](size_t i, size_t j)const;
  T val(size_t i, size_t j)const;
  std::span<T> operator[](size_t i);
  std::span<const T> operator[](size_t i)const;
  bool owned()const{ return !is_view(); }
  private:
  T* data(size_t i=0, size_t j=0)const{
    return data_+i*leading_.col+j;
  }
  public:
  
  // split.merge.hpp
  private:
  std::array<matrix,4> split()const;
  static void merge(matrix& dst,
    const size_t row_offset, const size_t col_offset,
    const matrix&m00,const matrix&m01,
    const matrix&m10,const matrix&m11);
  public:

  // operator/unary.hpp
  operator std::string()const;
  std::string string()const;
  void transpose();
  matrix<T> transpose()const;
  void negate();
  matrix<T> operator-()const;

  // operator/scalar.hpp
  template<class U> matrix<T>& operator*=(const U&);
  template<class U> matrix<T>& operator/=(const U&);
  template<class U,class V,class W>
  friend matrix<W> operator*(const U&,const matrix<V>&);
  template<class U,class V,class W>
  friend matrix<W> operator*(const matrix<U>&,const V&);
  template<class U,class V,class W>
  friend matrix<W> operator/(const matrix<U>&,const V&);

  // operator/matrix_addsub.hpp
  template<class U> matrix<T>& operator+=(const matrix<U>&);
  template<class U> matrix<T>& operator-=(const matrix<U>&);
  template<class U,class V,class W>
  friend matrix<W> operator+(const matrix<U>&,const matrix<V>&);
  template<class U,class V,class W>
  friend matrix<W> operator-(const matrix<U>&,const matrix<V>&);

  // operator/matrix_mulpow.hpp
  template<class U,class V> // *this+=A*B
  matrix<T>& fma(const matrix<U>&,const matrix<V>&);
  template<class U> matrix<T>& operator*=(const matrix<U>&);
  template<class U,class V,class W>
  friend matrix<W> operator*(const matrix<U>&,const matrix<V>&);
  void pow(uint64_t);
  matrix<T> pow(uint64_t)const;
  private:
  template<class U, class V>
  matrix<T>&fma_impl_naive(const matrix<U>&,const matrix<V>&);
  template<class U, class V>
  matrix<T>&fma_impl_blocked(const matrix<U>&,const matrix<V>&);
  template<class U, class V>
  matrix<T>&fma_impl_strassen(const matrix<U>&,const matrix<V>&);
  matrix<T> pow_impl(matrix<T>&&k,uint64_t b);
  public:

  // operation/reduction.hpp
  // io.hpp
  template<class Char, class Traits, class U>
  friend std::basic_ostream<Char, Traits>& operator<<(std::basic_ostream<Char, Traits>& os, const matrix<U>& mat);
  template<class Char, class Traits>
  void read(std::basic_istream<Char, Traits>& is) requires (!matrix_io::modint_v<T>);
  template<bool validated=true,class Char, class Traits>
  void read(std::basic_istream<Char, Traits>& is) requires (matrix_io::modint_v<T>);
};

} // namespace elsie
#endif