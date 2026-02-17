#ifndef ELSIE_MATRIX_CONSTRUCTOR_DESTRUCTOR_OPERATOR_EQ_HPP
#define ELSIE_MATRIX_CONSTRUCTOR_DESTRUCTOR_OPERATOR_EQ_HPP

#include "./_image.hpp"
#include <algorithm>

namespace elsie{

template<class T>
matrix<T>::matrix():dim_(),leading_(),capacity_data(0),data_(nullptr){}

template<class T>
matrix<T>::matrix(const dimension&dim)
  :dim_(dim),leading_(ceil_to_alignment(dim.row),ceil_to_alignment(dim.col))
  ,capacity_data(calc_capacity_data(leading_.capacity()))
  ,data_(allocate(get_capacity())){
  default_construct(data_,get_capacity());
}

template<class T>
matrix<T>::matrix(const dimension&dim, const T& init)
  :dim_(dim),leading_(ceil_to_alignment(dim.row),ceil_to_alignment(dim.col))
  ,capacity_data(calc_capacity_data(leading_.capacity()))
  ,data_(allocate(get_capacity())
){
  default_construct(data_,get_capacity(),init);
}

template<class T>
matrix<T>::matrix(const size_t row, const size_t col)
  :matrix(dimension(row,col)){}

template<class T>
matrix<T>::matrix(const size_t row, const size_t col, const T& init)
  :matrix(dimension(row,col),init){}

template<class T>
matrix<T>::matrix(const matrix<T>& other)
  :dim_(other.dim_),leading_(other.leading_)
  ,capacity_data(other.capacity_data|1)
  ,data_(other.data_){}

template<class T>
matrix<T>::matrix(matrix<T>&& other)
  :dim_(other.dim_),leading_(other.leading_)
  ,capacity_data(other.capacity_data)
  ,data_(other.data_){
  other.data_=nullptr;
}

template<class T>
matrix<T>& matrix<T>::operator=(const matrix<T>& other){
  if(owned()) free(data_);
  dim_=other.dim_;
  leading_=other.leading_;
  capacity_data=other.capacity_data|1;
  data_=other.data_;
  return *this;
}

template<class T>
matrix<T>& matrix<T>::operator=(matrix<T>&& other){
  if(owned()) free(data_);
  dim_=other.dim_;
  leading_=other.leading_;
  capacity_data=other.capacity_data;
  data_=other.data_;
  return *this;
}

template<class T>
matrix<T>::~matrix(){
  if(owned()) free(data_);
  dim_=leading_={0,0};
  capacity_data=1; // 破棄済みならviewとして次の代入時の~matrix()を呼ばせない
  data_=nullptr;
}


// private constructor
template<class T>
matrix<T>::matrix(const matrix<T>& other, size_t row, size_t col, size_t row_offset, size_t col_offset)
  :dim_(row,col),leading_(other.leading_)
  ,capacity_data(calc_capacity_data(other.get_capacity()-row*col,true))
  ,data_(other.data(row_offset,col_offset)){}

template<class T>
matrix<T> matrix<T>::make_view(size_t row, size_t col, size_t row_offset, size_t col_offset)const{
  return matrix<T>(*this, row, col, row_offset, col_offset);
}

template<class T>
matrix<T> matrix<T>::copy()const{
  matrix<T> ret(dim());
  for(size_t i=0;i<dim_.row;++i)
    std::copy(data(i,0), data(i+1,0), ret.data(i,0));
  return ret;
}

template<class T>
matrix<T> matrix<T>::copy(size_t row, size_t col, size_t row_offset, size_t col_offset)const{
  return this->make_view(row, col, row_offset, col_offset).copy();
}
} // namespace elsie
#endif