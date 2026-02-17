// MIT License
// Copyright 2026 CleverElsie

#ifndef ELSIE_MATRIX_SPLIT_MERGE_HPP
#define ELSIE_MATRIX_SPLIT_MERGE_HPP

#include "./_image.hpp"
#include <algorithm>

namespace elsie{

template<class T>
std::array<matrix<T>,4>
matrix<T>::split()const{
  // 4分割する
  // ただし，長い方の半分未満の短い方である場合，そちらは行か列が0次元になる
  // 行列の長い方をWとし，左上成分がW/2 x W/2 の行列になるように分割する．
  // ただし，W/2は切り上げ
  const auto[row,col]=dim();
  const ssize_t W=(std::max(row,col)+1)>>1;
  // W x W,   W x C-W
  // R-W x W, R-W x C-W
  const ssize_t RW=std::max<ssize_t>(ssize_t(0),(ssize_t)row-W);
  const ssize_t CW=std::max<ssize_t>(ssize_t(0),(ssize_t)col-W);
  return std::array<matrix<T>,4>{
    make_view(W,W,0,0),  make_view(W,CW,0,W),
    make_view(RW,W,W,0), make_view(RW,CW,W,W)
  };
}

template<class T>
void matrix<T>::merge(
  matrix& dst,
  const size_t row_offset, const size_t col_offset,
  const matrix&m00, const matrix&m01,
  const matrix&m10, const matrix&m11
){
  [[assume(m00.dim().row==m01.dim().row)]];
  [[assume(m00.dim().col==m10.dim().col)]];
  [[assume(m10.dim().row==m11.dim().row)]];
  [[assume(m01.dim().col==m11.dim().col)]];
  [[assume(row_offset+m00.dim().row+m10.dim().row<=dst.dim().row)]];
  [[assume(col_offset+m00.dim().col+m01.dim().col<=dst.dim().col)]];

  const auto[Wr,Wc]=m00.dim();
  const auto[Dr,Dc]=m10.dim();
  for(size_t i=0;i<Wr;++i){
    std::copy(m00.data(i,0), m00.data(i+1,0), dst.data(row_offset+i, col_offset));
    std::copy(m01.data(i,0), m01.data(i+1,0), dst.data(row_offset+i, col_offset+Wc));
  }
  for(size_t i=0;i<Dr;++i){
    std::copy(m10.data(i,0), m10.data(i+1,0), dst.data(row_offset+Wr+i, col_offset));
    std::copy(m11.data(i,0), m11.data(i+1,0), dst.data(row_offset+Wr+i, col_offset+Dc));
  }
}

} // namespace elsie
#endif