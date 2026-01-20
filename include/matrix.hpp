#ifndef ELSIE_MATRIX
#define ELSIE_MATRIX
#include <cassert>
#include <cstddef>
#include <cstdint>

#include <iostream>

#include <bit>
#include <new>
#include <memory>

#include <span>
#include <array>
#include <vector>
#include <utility>

#include <algorithm>
#include <iterator>
#include <concepts>
#include <type_traits>
#include <format>

#if 1
#include "./detail/_image.hpp"
#include "./detail/constructor.destructor.operator=.hpp"
#include "./detail/access.hpp"
#include "./detail/split.merge.hpp"
#include "./detail/operator/unary.hpp"
#include "./detail/operator/scalar.hpp"
#include "./detail/operator/matrix_addsub.hpp"
#include "./detail/operator/matrix_mulpow.hpp"
#include "./detail/operator/reduction.hpp"
#include "./detail/io.hpp"
#else
namespace elsie{

template<class T> class matrix{
  using sz_t = int32_t;
  size_t row,col,capacity;
  std::vector<T> data;
  public:
  matrix()=default;
	~matrix()=default;
  matrix(matrix&&)=default;
  matrix(const matrix&)=default;
	matrix&operator=(matrix&&)=default;
	matrix&operator=(const matrix&)=default;

  matrix(size_t row_,size_t col_)
  :row(row_),col(col_),capacity(row*col),data(capacity){}

  matrix(size_t row_,size_t col_,const T&init)
  :row(row_),col(col_),capacity(row*col),data(capacity,init){}

  template<class S>
  requires std::same_as<std::decay_t<S>,std::vector<T>>
  && (!std::same_as<std::decay_t<S>,T&>)
  matrix(size_t row_,size_t col_,S&&data_)
  :row(row_),col(col_),capacity(row*col),data(std::forward<S>(data_)){ assert(row&&col); }

  template<size_t row_,size_t col_>
  matrix(const std::array<std::array<T,col_>,row_>&data_)
  :row(row_),col(col_),capacity(row*col),data(capacity){
    static_assert(row_&&col_);
    auto dtr=data.begin();
    for(const auto&y:data_)
    for(const auto&x:y) *(dtr++)=x;
  }

  matrix(const std::vector<std::vector<T>>&data_)
  :row(data_.size()),col(data_[0].size()),capacity(row*col),data(capacity){
    assert(row&&col);
    auto dtr=data.begin();
    for(const auto&y:data_){
      assert(y.size()==col);
      for(const auto&x:y) *(dtr++)=x;
    }
  }

  inline std::pair<size_t,size_t>size()const{return {row,col};}
  inline T& operator[](size_t i,size_t j){ return data[i*col+j]; }
  inline const T& operator[](size_t i,size_t j)const{ return data[i*col+j]; }
  inline std::span<T> operator[](size_t i){
    const auto ptr=data.begin()+i*col;
    return {ptr,ptr+col};
  }
  inline std::span<const T> operator[](size_t i)const{
    const auto ptr=data.begin()+i*col;
    return {ptr,ptr+col};
  }

  inline T val(size_t i,size_t j)const{
    if(j>=col||i>=row)return T();
    return data[i*col+j];
  }
  private:
  public:
  template<class Char, class Traits>
  inline void print(std::basic_ostream<Char,Traits>& os=std::cout)const{
    const size_t sz_=row*col;
    for(size_t i=0;i<sz_;i+=col)
      for(size_t j=0;j<col;++j)
        os<<data[i+j].val()<<" \n"[j==col-1];
  }

  void transpose();
  matrix<T> pow(uint64_t b)const;
  private:
  enum class ReduceClass{
    rank,solve_LE,solve_LE_with_basis
  };
  template<ReduceClass reduce_class>
  std::tuple<sz_t,std::vector<T>,matrix> reduce_impl();
  public:
  sz_t rank()const;
  sz_t reduce();
  std::pair<sz_t,std::vector<T>> solve_linear_equation();
  std::tuple<sz_t,std::vector<T>,matrix> solve_linear_equation_with_basis();
  std::pair<sz_t,matrix> inverse()const;

  using iterator=std::vector<T>::iterator;
  using const_iterator=std::vector<T>::const_iterator;
  using reverse_iterator=std::vector<T>::reverse_iterator;
  using const_reverse_iterator=std::vector<T>::const_reverse_iterator;
  inline iterator begin(){return data.begin();}  
  inline iterator end(){return data.end();}
  inline const_iterator begin()const{return data.cbegin();}
  inline const_iterator end()const{return data.cend();}
  inline const_iterator cbegin()const{return data.cbegin();}
  inline const_iterator cend()const{return data.cend();}
  inline reverse_iterator rbegin(){return data.rbegin();}
  inline reverse_iterator rend(){return data.rend();}
  inline const_reverse_iterator rbegin()const{return data.crbegin();}
  inline const_reverse_iterator rend()const{return data.crend();}
  inline const_reverse_iterator crbegin()const{return data.crbegin();}
  inline const_reverse_iterator crend()const{return data.crend();}


  // スカラ倍
  template<class U> matrix& operator*=(const U&);
  template<class U> matrix& operator/=(const U&);
	template<class U,class V,class W> friend matrix<W>& operator*(const U&,const matrix<V>&);
	template<class U,class V,class W> friend matrix<W>& operator*(const matrix<U>&,const V&);
	template<class U,class V,class W> friend matrix<W>& operator/(const matrix<U>&,const V&);

  // 行列演算
  template<class U> matrix& operator+=(const matrix<U>&);
  template<class U> matrix& operator-=(const matrix<U>&);
  template<class U> matrix& operator*=(const matrix<U>&);
  template<class U,class V,class W> friend matrix<W> operator+(const matrix<U>&,const matrix<V>&);
  template<class U,class V,class W> friend matrix<W> operator-(const matrix<U>&,const matrix<V>&);
  template<class U,class V,class W> friend matrix<W> operator*(const matrix<U>&,const matrix<V>&);

	protected:
  /**
   * @brief a:NxN matrix. N is 2^i and N>1.
   * @return {A_00, A_01, A_10, A_11}
   */
  template<bool by2>
  std::array<matrix,4>split(size_t n)const;
  static matrix merge(const matrix&,const matrix&,const matrix&,const matrix&);

  template<size_t naive, size_t base,
    bool q2,class U,class V,class W>
  friend matrix<W> matrix_mul(const matrix<U>&,const matrix<V>&);
  
  template<class U,class V,class W>
  friend matrix<W> matrix_mul_naive_impl(const matrix<U>&,const matrix<V>&,size_t,size_t,size_t);

  template<size_t naive,class U,class V,class W>
  friend matrix<W> matrix_mul_blocked_impl(const matrix<U>&,const matrix<V>&,size_t,size_t,size_t);

  template<size_t naive,size_t base,bool q2,class U,class V,class W>
  friend matrix<W> matrix_mul_strassen_impl(const matrix<U>&,const matrix<V>&);
}; // class matrix
 
///////////////////
///////////////////
///////////////////
/// matrix impl ///
///////////////////
///////////////////
///////////////////

template<class T>
template<bool by2>
std::array<matrix<T>,4> matrix<T>::split(const size_t n)const{
  assert(std::popcount(n)==1);
  const size_t n2=n>>1;
  std::array<matrix<T>,4>r{matrix<T>(n2,n2),matrix<T>(n2,n2),matrix<T>(n2,n2),matrix<T>(n2,n2)};
  if constexpr(by2){
    auto db=data.cbegin();
    auto db0=r[0].begin(),db1=r[1].begin(),db2=r[2].begin(),db3=r[3].begin();
    auto end1=r[0].end();
    for(;db0<end1;db+=n,db0+=n2,db1+=n2){
      std::copy(db,db+n2,db0);
      std::copy(db+n2,db+n,db1);
    }
    auto end2=r[2].end();
    for(;db2<end2;db+=n,db2+=n2,db3+=n2){
      std::copy(db,db+n2,db2);
      std::copy(db+n2,db+n,db3);
    }
  }else{
    if constexpr(true){
      auto d=data.cbegin();
      auto d0=r[0].data.begin();
      auto d1=r[1].data.begin();
      auto d2=r[2].data.begin();
      auto d3=r[3].data.begin();
      if(row<n2){
        if(col<n2){
          //10
          //00
          auto e=data.cend();
          for(;d<e;d+=col,d0+=n2){
            std::copy(d,d+col,d0);
            std::fill(d0+col,d0+n2,T());
          }
        }else{
          //11
          //00
          auto e=data.cend();
          if(col!=n)
          for(;d<e;d+=col,d0+=n2,d1+=n2){
            std::copy(d,d+n2,d0);
            std::copy(d+n2,d+col,d1);
            std::fill(d1+col-n2,d1+n2,T());
          }else for(;d<e;d+=n,d0+=n2,d1+=n2){
            std::copy(d,d+n2,d0);
            std::copy(d+n2,d+n,d1);
          }
          std::fill(d0,r[0].end(),T());
        }
        std::fill(d0,r[0].data.end(),T());
        std::fill(d1,r[1].data.end(),T());
        std::fill(d2,r[2].data.end(),T());
        std::fill(d3,r[3].data.end(),T());
      }else{
        if(col<n2){
          //10
          //10
          auto e0=r[0].data.end(),e=data.end();
          for(;d0<e0;d+=col,d0+=n2){
            std::copy(d,d+col,d0);
            std::fill(d0+col,d0+n2,T());
          }
          for(;d<e;d+=col,d2+=n2){
            std::copy(d,d+col,d2);
            std::fill(d2+col,d2+n2,T());
          }
          std::fill(d2,r[2].data.end(),T());
          std::fill(d1,r[1].data.end(),T());
          std::fill(d3,r[3].data.end(),T());
        }else{
          //11
          //11
          auto e0=r[0].data.end(),e=data.end();
          if(col!=n){
            for(;d0<e0;d+=col,d0+=n2,d1+=n2){
              std::copy(d,d+n2,d0);
              std::copy(d+n2,d+col,d1);
              std::fill(d1+col-n2,d1+n2,T());
            }
            for(;d<e;d+=col,d2+=n2,d3+=n2){
              std::copy(d,d+n2,d2);
              std::copy(d+n2,d+col,d3);
              std::fill(d3+col-n2,d3+n2,T());
            }
          }else{
            for(;d0<e0;d+=n,d0+=n2,d1+=n2){
              std::copy(d,d+n2,d0);
              std::copy(d+n2,d+n,d1);
            }
            for(;d<e;d+=col,d2+=n2,d3+=n2){
              std::copy(d,d+n2,d2);
              std::copy(d+n2,d+n,d3);
            }
          }
        }
      }
    }else{
      auto&ref0=r[0];
      auto&ref1=r[1];
      auto&ref2=r[2];
      auto&ref3=r[3];
      size_t i,j;
      for(i=0;i<n2;++i){
        for(j=0;j<n2;++j) ref0[i,j]=val(i,j);
        for(;j<n;++j) ref1[i,j-n2]=val(i,j);
      }
      for(;i<n;++i){
        for(j=0;j<n2;++j) ref2[i-n2,j]=val(i,j);
        for(;j<n;++j) ref3[i-n2,j-n2]=val(i,j);
      }
    }
  }
  return r;
}

template<class T>
matrix<T> matrix<T>::merge(const matrix<T>&c11,const matrix<T>&c12,const matrix<T>&c21,const matrix<T>&c22){
  // n satisfy n=2^c, n>=256
  // 11 12
  // 21 22
  const size_t n=c11.row+c21.row, m=c11.col+c12.col;
  [[assume(c11.row==c12.row)]];
  [[assume(c21.row==c22.row)]];
  [[assume(c11.col==c21.col)]];
  [[assume(c12.col==c22.col)]];
  matrix<T> c(n,m);
  const size_t upper_capacity=c11.capacity+c12.capacity;
  for(size_t i=0,i11=0,i12=0;i<upper_capacity;i+=m,i11+=c11.col,i12+=c12.col){
    std::copy(c11.data.begin()+i11,c11.data.begin()+i11+c11.col,c.data.begin()+i);
    std::copy(c12.data.begin()+i12,c12.data.begin()+i12+c12.col,c.data.begin()+i+c11.col);
  }
  for(size_t i=upper_capacity,i21=0,i22=0;i<c.capacity;i+=m,i21+=c21.col,i22+=c22.col){
    std::copy(c21.data.begin()+i21,c21.data.begin()+i21+c21.col,c.data.begin()+i);
    std::copy(c22.data.begin()+i22,c22.data.begin()+i22+c22.col,c.data.begin()+i+c21.col);
  }
  return c;
}

template<class U,class V,class W=std::common_type_t<U,V>>
inline __attribute__((always_inline)) __attribute__((optimize("O3"))) matrix<W> matrix_mul_naive_impl(const matrix<U>&a,const matrix<V>&b,size_t N,size_t M,size_t L){
  matrix<W>c(N,L,W());
  uint32_t i,j,k;
  uint32_t x,y,z;
  uint32_t ia,ic,kb;
  constexpr const uint16_t unrool=8-1;
  for(i=0,ic=0,ia=0;i<N;++i,ic+=c.col,ia+=a.col)
  for(k=0,kb=0;k<M;++k,kb+=b.col){
    const auto&ak=a.data[ia+k];
    for(j=0;j+unrool<L;j+=(unrool+1)){
      c.data[ic+j  ]+=ak*b.data[kb+j  ];
      c.data[ic+j+1]+=ak*b.data[kb+j+1];
      c.data[ic+j+2]+=ak*b.data[kb+j+2];
      c.data[ic+j+3]+=ak*b.data[kb+j+3];
      c.data[ic+j+4]+=ak*b.data[kb+j+4];
      c.data[ic+j+5]+=ak*b.data[kb+j+5];
      c.data[ic+j+6]+=ak*b.data[kb+j+6];
      c.data[ic+j+7]+=ak*b.data[kb+j+7];
    }
    for(j=L&(~unrool);j<L;++j)
      c.data[ic+j]+=ak*b.data[kb+j];
  }
  return c;
}

template<size_t naive,class U,class V,class W=std::common_type_t<U,V>>
inline __attribute__((always_inline)) __attribute__((optimize("O3"))) matrix<W> matrix_mul_blocked_impl(const matrix<U>&a,const matrix<V>&b,size_t N,size_t M,size_t L){
  matrix<W>c(N,L,W());
  constexpr static uint16_t bs=naive; // block size
  constexpr static uint16_t unrool=8-1;
  uint16_t I,J,K,i,j,k,it,jt,kt;
  for(I=0,it=std::min<uint16_t>(bs,N);I<N;I+=bs,it=std::min<uint16_t>(it+bs,N))
  for(J=0,jt=std::min<uint16_t>(bs,L);J<L;J+=bs,jt=std::min<uint16_t>(jt+bs,L))
  for(K=0,kt=std::min<uint16_t>(bs,M);K<M;K+=bs,kt=std::min<uint16_t>(kt+bs,M))
    for(i=I;i<it;++i){
      uint16_t ci=i*c.col,ai=i*a.col+K;
      for(k=K;k<kt;++k,++ai){
        uint16_t bk=k*b.col;
        const auto&ak=a.data[ai];
        for(j=J;j+unrool<jt;j+=(unrool+1)){
          c.data[ci+j  ]+=ak*b.data[bk+j  ];
          c.data[ci+j+1]+=ak*b.data[bk+j+1];
          c.data[ci+j+2]+=ak*b.data[bk+j+2];
          c.data[ci+j+3]+=ak*b.data[bk+j+3];
          c.data[ci+j+4]+=ak*b.data[bk+j+4];
          c.data[ci+j+5]+=ak*b.data[bk+j+5];
          c.data[ci+j+6]+=ak*b.data[bk+j+6];
          c.data[ci+j+7]+=ak*b.data[bk+j+7];
        }
        for(j=jt&(~unrool);j<jt;++j)
          c.data[ci+j]+=a.data[ai]*b.data[bk+j];
      }
    }
  return c;
}

template<size_t naive,size_t base,bool q2,class U,class V,class W>
matrix<W> matrix_mul_strassen_impl(const matrix<U>&a,const matrix<V>&b){
  const size_t N=a.row,M=b.row,L=b.col;
  const auto[Q,P]=std::minmax({N,M,L});
  size_t n;
  if constexpr(q2) n=std::bit_ceil(P);
  else n=L;
  auto as=a.template split<!q2>(n);
  auto bs=b.template split<!q2>(n);
  matrix<U> s2=as[0]+as[1], s9=as[0]-as[2], s5=as[0]+as[3], s3=as[2]+as[3], s7=as[1]-as[3];
  matrix<V> s1=bs[1]-bs[3], s8=bs[2]+bs[3], s6=bs[0]+bs[3], s4=bs[2]-bs[0], sA=bs[0]+bs[1];

  matrix<W> p1=matrix_mul<naive,base,false>(as[0],s1);
  matrix<W> p2=matrix_mul<naive,base,false>(s2,bs[3]);
  matrix<W> p3=matrix_mul<naive,base,false>(s3,bs[0]);
  matrix<W> p4=matrix_mul<naive,base,false>(as[3],s4);
  matrix<W> p5=matrix_mul<naive,base,false>(s5,s6);
  matrix<W> p6=matrix_mul<naive,base,false>(s7,s8);
  matrix<W> p7=matrix_mul<naive,base,false>(s9,sA);

  auto r=matrix<W>::merge(
    ((p5+p4)-=p2)+=p6,p1+p2,
    p3+p4, ((p5+p1)-=p3)-=p7
  );
  if constexpr(q2)
    if(N!=r.row||L!=r.col){
      matrix<W> ret(N,L);
      for(size_t i=0;i<N;++i){
        auto reti=ret[i];
        const auto rri=r[i];
        for(size_t j=0;j<L;++j)
          reti[j]=rri[j];
      }
      r=std::move(ret);
    }
  return r;
}

template<size_t naive=16,size_t base=64,bool q2=true,class U,class V,class W=std::common_type_t<U,V>>
matrix<W> __attribute__((optimize("O3"))) matrix_mul(const matrix<U>&a,const matrix<V>&b){
  if(a.col!=b.row) [[unlikely]]
    throw std::invalid_argument("matrix_mul: incompatible dimensions");
  const size_t N=a.row,M=b.row,L=b.col;
  const auto[Q,P]=std::minmax({N,M,L});
  if(P<base){
    if(Q<naive) return matrix_mul_naive_impl(a,b,N,M,L);
    else return matrix_mul_blocked_impl<naive>(a,b,N,M,L);
  }else return matrix_mul_strassen_impl<naive,base,q2,U,V,W>(a,b);
}

template<class T>
template<class U>
matrix<T>& matrix<T>::operator+=(const matrix<U>&rhs){
  assert(row==rhs.row&&col==rhs.col);
  constexpr size_t unrool=8-1;
  for(size_t i=0;i+unrool<data.size();i+=(unrool+1)){
    data[i   ]+=rhs.data[i   ];
    data[i+ 1]+=rhs.data[i+ 1];
    data[i+ 2]+=rhs.data[i+ 2];
    data[i+ 3]+=rhs.data[i+ 3];
    data[i+ 4]+=rhs.data[i+ 4];
    data[i+ 5]+=rhs.data[i+ 5];
    data[i+ 6]+=rhs.data[i+ 6];
    data[i+ 7]+=rhs.data[i+ 7];
  }
  for(size_t i=data.size()&(~unrool);i<data.size();++i)
    data[i]+=rhs.data[i];
  return*this;
}

template<class T>
template<class U>
matrix<T>& matrix<T>::operator-=(const matrix<U>&rhs){
  assert(row==rhs.row&&col==rhs.col);
  constexpr size_t unrool=8-1;
  for(size_t i=0;i+unrool<data.size();i+=(unrool+1)){
    data[i   ]-=rhs.data[i   ];
    data[i+ 1]-=rhs.data[i+ 1];
    data[i+ 2]-=rhs.data[i+ 2];
    data[i+ 3]-=rhs.data[i+ 3];
    data[i+ 4]-=rhs.data[i+ 4];
    data[i+ 5]-=rhs.data[i+ 5];
    data[i+ 6]-=rhs.data[i+ 6];
    data[i+ 7]-=rhs.data[i+ 7];
  }
  for(size_t i=data.size()&(~unrool);i<data.size();++i)
    data[i]-=rhs.data[i];
  return*this;
}

template<class T>
template<class U>
matrix<T>& matrix<T>::operator*=(const U&rhs){
  for(auto&x:data)x*=rhs;
  return*this;
}

template<class U,class V,class W=std::common_type_t<U,V>>
matrix<W>& operator*(const U&lhs,const matrix<V>&rhs){
  auto [row,col]=rhs.size();
  size_t sz=rhs.data.size();
  std::vector<W> data;
  data.reserve(sz);
  for(const auto&x:rhs)
    data.emplace_back(lhs*x);
  return matrix<W>(row,col,std::move(data));
}

template<class U,class V,class W=std::common_type_t<U,V>>
matrix<W>& operator*(const matrix<U>&lhs,const V&rhs){
  auto [row,col]=lhs.size();
  size_t sz=lhs.data.size();
  std::vector<W> data;
  data.reserve(sz);
  for(const auto&x:lhs)
    data.emplace_back(x*rhs);
  return matrix<W>(row,col,std::move(data));
}

template<class T>
template<class U>
matrix<T>& matrix<T>::operator/=(const U&rhs){
  for(auto&x:data)x/=rhs;
  return*this;
}

template<class T>
template<class U>
matrix<T>& matrix<T>::operator*=(const matrix<U>&rhs){
  assert(col==rhs.col&&col==rhs.row);
  *this=matrix_mul(*this,rhs);
  return*this;
}

template<class U,class V,class W=std::common_type_t<U,V>>
matrix<W>& operator/(const matrix<U>&lhs,const V&rhs){
  auto [row,col]=lhs.size();
  size_t sz=lhs.data.size();
  std::vector<W> data;
  data.reserve(sz);
  for(const auto&x:lhs)
    data.emplace_back(x/rhs);
  return matrix<W>(row,col,std::move(data));
}

template<class U,class V,class W=std::common_type_t<U,V>>
matrix<W> operator+(const matrix<U>&lhs,const matrix<V>&rhs){
  assert(lhs.row==rhs.row&&lhs.col==rhs.col);
  matrix<W> res(lhs.row,lhs.col);
  constexpr size_t unrool=8-1;
  for(size_t i=0;i+unrool<res.data.size();i+=(unrool+1)){
    res.data[i   ]=lhs.data[i   ]+rhs.data[i   ];
    res.data[i+ 1]=lhs.data[i+ 1]+rhs.data[i+ 1];
    res.data[i+ 2]=lhs.data[i+ 2]+rhs.data[i+ 2];
    res.data[i+ 3]=lhs.data[i+ 3]+rhs.data[i+ 3];
    res.data[i+ 4]=lhs.data[i+ 4]+rhs.data[i+ 4];
    res.data[i+ 5]=lhs.data[i+ 5]+rhs.data[i+ 5];
    res.data[i+ 6]=lhs.data[i+ 6]+rhs.data[i+ 6];
    res.data[i+ 7]=lhs.data[i+ 7]+rhs.data[i+ 7];
  }
  for(size_t i=res.data.size()&(~unrool);i<res.data.size();++i)
    res.data[i]=lhs.data[i]+rhs.data[i];
  return res;
}

template<class U,class V,class W=std::common_type_t<U,V>>
matrix<W> operator-(const matrix<U>&lhs,const matrix<V>&rhs){
  assert(lhs.row==rhs.row&&lhs.col==rhs.col);
  matrix<W> res(lhs.row,lhs.col);
  constexpr size_t unrool=8-1;
  for(size_t i=0;i+unrool<res.data.size();i+=(unrool+1)){
    res.data[i   ]=lhs.data[i   ]-rhs.data[i   ];
    res.data[i+ 1]=lhs.data[i+ 1]-rhs.data[i+ 1];
    res.data[i+ 2]=lhs.data[i+ 2]-rhs.data[i+ 2];
    res.data[i+ 3]=lhs.data[i+ 3]-rhs.data[i+ 3];
    res.data[i+ 4]=lhs.data[i+ 4]-rhs.data[i+ 4];
    res.data[i+ 5]=lhs.data[i+ 5]-rhs.data[i+ 5];
    res.data[i+ 6]=lhs.data[i+ 6]-rhs.data[i+ 6];
    res.data[i+ 7]=lhs.data[i+ 7]-rhs.data[i+ 7];
  }
  for(size_t i=res.data.size()&(~unrool);i<res.data.size();++i)
    res.data[i]=lhs.data[i]-rhs.data[i];
  return res;
}

template<class U,class V,class W=std::common_type_t<U,V>>
matrix<W> operator*(const matrix<U>&lhs,const matrix<V>&rhs){
  return matrix_mul(lhs,rhs);
}

template<class T>
void matrix<T>::transpose(){
  std::vector<T> res(col*row);
  constexpr const size_t block_size=64;
  for(size_t i=0;i<row;i+=block_size){
    const size_t Erow=std::min(i+block_size,row);
    for(size_t j=0;j<col;j+=block_size){
      const size_t Ecol=std::min(j+block_size,col);
      for(size_t I=i;I<Erow;++I)
        for(size_t J=j;J<Ecol;++J)
          res[J*row+I]=data[I*col+J];
    }
  }
  std::swap(row,col);
  data=std::move(res);
}

template<class T>
matrix<T> matrix<T>::pow(uint64_t b)const{
  assert(b<=1||row==col);
  matrix<T> r(row,row,T());
  for(size_t i=0;i<row;++i) r[i,i]=1;
  matrix<T> k(*this);
  for(;b;b>>=1){
    if(b&1) r=r*k;
    k=k*k;
  }
  return r;
}

template<class T>
template<typename matrix<T>::ReduceClass reduce_class>
auto matrix<T>::reduce_impl() -> std::tuple<matrix<T>::sz_t,std::vector<T>,matrix<T>> {
  if(row==0||col==0) return {0,std::vector<T>(),matrix<T>()};
  const size_t col=this->col;
  std::vector<T> sol(col-1);
  std::vector<sz_t> p(row),q(col-1,-1);
  std::vector<bool> is_main(col-1,0);
  // O(n^3) なので32bitよりも大きいrow,colは考慮しない
  sz_t start=0;
  for(sz_t i=0;i<row;){
    const sz_t r=i*col;
    if(data[r+start]==0){
      // 全0行削除，最後だけ非ゼロで制約矛盾
      // 他行が非ゼロで入れ替え，それがないなら主成分を右へずらす
      sz_t nonzero=-1;
      for(sz_t j=start+1;j<col;++j)
        if(data[r+j]!=0){
          nonzero=j;
          break;
        }
      if(nonzero==-1){
        if(i+1!=row){
          const sz_t tmp_r=(row-1)*col;
          for(sz_t j=0;j<col;++j) data[r+j]=data[tmp_r+j];
        }
        --row;
        data.resize(data.size()-col);
        continue;
      }else if constexpr(reduce_class!=ReduceClass::rank)
        if(nonzero==col-1)
          return {-1,std::vector<T>(),matrix<T>()};
      bool empty=true;
      for(sz_t ii=i+1,rr=r+col;ii<row;++ii,rr+=col){
        if(data[rr+start]==0) continue;
        empty=false;
        auto zitr=data.begin()+r;
        auto itr=data.begin()+rr;
        auto end=itr+col;
        for(;itr<end;++itr,++zitr)
          swap(*zitr,*itr);
        break;
      }
      if(empty){
        ++start;
        continue;
      }
    }
    if constexpr(reduce_class!=ReduceClass::rank)
    if(start+1>=col) return {-1,std::vector<T>(),matrix<T>()};
    if constexpr(reduce_class==ReduceClass::solve_LE_with_basis){
      is_main[start]=true;
      p[i]=start;
    }
    const T divk=1/data[r+start];
    for(sz_t j=start;j<col;++j)
      data[r+j]*=divk;
    for(sz_t ii=i+1,f=r+col;ii<row;++ii,f+=col){
      if(data[f+start]==0) continue;
      const T mulk=data[f+start];
      for(sz_t j=start;j<col;++j)
        data[f+j]-=data[r+j]*mulk;
    }
    ++start,++i;
  }
  for(sz_t i=row-1,r=i*col;i<row;--i,r-=col){
    for(sz_t j=0;j<col;++j)if(data[r+j]!=0){
      start=j;
      break;
    }
    if constexpr(reduce_class != ReduceClass::rank)
      sol[start]=data[r+col-1];
    for(sz_t ii=0,rr=0;ii<i;++ii,rr+=col){
      const T subk=data[rr+start];
      for(sz_t j=start;j<col;++j)
        data[rr+j]-=subk*data[r+j];
    }
  }
  if constexpr(reduce_class==ReduceClass::rank)
    return {row,std::vector<T>(),matrix<T>()};
  else if constexpr(reduce_class==ReduceClass::solve_LE)
    return {row,std::move(sol),matrix<T>()};
  else if constexpr(reduce_class==ReduceClass::solve_LE_with_basis){
    start=0;
    const sz_t dim=col-1-row;
    matrix<T> b(dim,col-1,T(0));
    for(sz_t j=0,cur=-1;j<col-1;++j)
      if(!is_main[j]) b[q[j]=++cur,j]=1;
    for(sz_t i=0,r=0;i<row;++i,r+=col)
      for(sz_t j=0;j<col-1;++j)if(q[j]>=0)
        b[q[j],p[i]]=-data[r+j];
    return {row,std::move(sol),std::move(b)};
  }
}


template<class T>
matrix<T>::sz_t matrix<T>::rank()const{
  matrix<T> tmp(*this);
  return tmp.reduce();
}

template<class T>
matrix<T>::sz_t matrix<T>::reduce(){
  return get<0>(reduce_impl<ReduceClass::rank>());
}

template<class T>
auto matrix<T>::solve_linear_equation() -> std::pair<matrix<T>::sz_t,std::vector<T>> {
  auto [r,sol,b]=reduce_impl<ReduceClass::solve_LE>();
  return {r,std::move(b)};
}

template<class T>
auto matrix<T>::solve_linear_equation_with_basis()-> std::tuple<matrix<T>::sz_t,std::vector<T>,matrix<T>> {
  return reduce_impl<ReduceClass::solve_LE_with_basis>();
}

template<class T>
auto matrix<T>::inverse() const ->std::pair<sz_t,matrix<T>>{
  assert(row==col);
  const sz_t n=row;
  matrix<T> a(n,2*n,0);
  for(sz_t i=0;i<n;++i){
    for(sz_t j=0;j<n;++j)
      a[i,j]=(*this)[i,j];
    a[i,n+i]=1;
  }
  a.reduce();
  for(sz_t i=0;i<n;++i)
    if(a[i,i]!=1) return {-1,matrix<T>()};
  matrix<T> b(n,n);
  for(sz_t i=0;i<n;++i)
    for(sz_t j=0;j<n;++j)
      b[i,j]=a[i,j+n];
  return {n,b};
}

}// namespace elsie
#endif
#endif