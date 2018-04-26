#ifndef MATRIX_H
#define MATRIX_H
/*!
 * @file Matrix.h


 */


#include <limits> // for std::numeric_limits
#include <ostream>
#include <istream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <cassert>
#include <random>

#include <iterator>     // std::iterator, std::input_iterator_tag

#include <iostream>
#include <algorithm>
#include <memory>

#include "mySerializer.h"
#include "mySerializer.h"


inline double logit(double x){return std::log(x/(1.0-x));}



inline std::pair<double,double> logit(double x,double sd){
  return {std::log(x/(1.0-x)),sd/(x*(1.0-x))};}


inline double logistic(double x){return 1.0/(1.0+std::exp(-x));}


template<typename T1, typename T2>
std::pair<T1,T2>& operator+=(std::pair<T1,T2>& x, const std::pair<T1,T2>& other)
{
  x.first+=other.first;
  x.second+=other.second;
  return x;
}

inline double average(double x, double y){return 0.5*(x+y);}

inline double sqr(double x){return x*x;}


namespace Vector_Unary_Index_Function
{

  template<typename T>
  std::size_t i_max(const std::vector<T>& x)
  {
    std::size_t j=0;
    for (std::size_t i=1; i<x.size(); ++i)
      if (x[i]>x[j]) j=i;
    return j;
  }
}

namespace Container_Unary_Functions
{

  template <class V>
  double mean(const V& x)
  {
    double sum=0;
    for (std::size_t i=0; i<x.size(); ++i)
      sum+=x[i];
    return sum/x.size();
  }

  template <class V>
  double stddev(const V& x)
  {
    double sum=0;
    double sum2=0;
    for (std::size_t i=0; i<x.size(); ++i)
      {
        sum+=x[i];sum2+=sqr(x[i]);
      }
    return std::sqrt(sum2/x.size()-sqr(sum/x.size()));
  }
}


namespace Container_Binary_Transformations
{

  namespace additive
  {

    template<template<typename...>class M,typename... Ts>
    M<Ts...>& operator_additive_assigment(M<Ts...>& itself, const M<Ts...>&  x)
    {
      for (size_t i=0; i<itself.size(); i++)
        itself[i]+=x[i];
      return itself;
    }





    /** @name Aritmetic Assigment Operations between a Matrix and a scalar
    (single element)
      */
    //@{
    /**
     Scalar Adition assignment.
     @returns a reference to itself
     @post all the values of the matrix are summed up by the value x
     */

    template<template<typename...>class M,typename T, typename... Ts>
    M<Ts...>& operator_additive_assigment(M<Ts...>& itself, T x)
    {
      for (size_t i=0; i<itself.size(); i++)
        itself[i]+=x;
      return itself;
    }
  }
  template<typename T>
  std::vector<T>& operator+=(std::vector<T>& itself, const std::vector<T>& other)
  {
    return additive::operator_additive_assigment(itself,other);
  }

  template<typename T>
  std::vector<T>& operator+=(std::vector<T>& itself, const T& other)
  {
    return additive::operator_additive_assigment(itself,other);
  }

}






namespace lapack
{
  extern "C" void dgetrf_(int *M,
                          int* N,
                          double *A,
                          int* LDA,
                          int* IPIV,
                          int * INFO );

  extern "C" void dsytrf_( char *UPLO,
                           int *N,
                           double* A,
                           int *LDA,
                           int* IPIV,
                           double *WORK,
                           int *LWORK,
                           int *INFO );

  extern "C" void dgetri_(int* n,
                          double *B,
                          int* dla,
                          int* ipiv,
                          double* work1,
                          int* lwork,
                          int* info);


  extern "C" void dsytri_(char* UPLO,
                          int * N,
                          double */*dimension( lda, * )*/A,
                          int *	LDA,
                          int */* dimension( * ) */ 	IPIV,
                          double * /*dimension( * )*/  	WORK,
                          int *	INFO );

  extern "C" void dgemm_(char * 	TRANSA,
                         char * 	TRANSB,
                         int * 	M,
                         int * 	N,
                         int * 	K,
                         double * ALPHA,
                         double * A,
                         int * 	LDA,
                         double * B,
                         int * 	LDB,
                         double * BETA,
                         double * C,
                         int * 	LDC
                         );

  extern "C" void dsymm_ 	(char*  	SIDE,
                                 char*  	UPLO,
                                 int*  	M,
                                 int*  	N,
                                 double*  	ALPHA,
                                 double * /*, dimension(lda,*)*/  	A,
                                 int*  	LDA,
                                 double */*, dimension(ldb,*)*/  	B,
                                 int*  	LDB,
                                 double *  	BETA,
                                 double * /*, dimension(ldc,*) */ 	C,
                                 int*  	LDC
                                 );
}



/**

   SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
  *
  *  -- LAPACK routine (version 3.3.1) --
  *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  *  -- April 2011                                                      --
  *
  *     .. Scalar Arguments ..
        CHARACTER          UPLO
        INTEGER            INFO, LDA, N
  *     ..
  *     .. Array Arguments ..
        DOUBLE PRECISION   A( LDA, * )
  *     ..
  *
  *  Purpose
  *  =======
  *
  *  DPOTRF computes the Cholesky factorization of a real symmetric
  *  positive definite matrix A.
  *
  *  The factorization has the form
  *     A = U**T * U,  if UPLO = 'U', or
  *     A = L  * L**T,  if UPLO = 'L',
  *  where U is an upper triangular matrix and L is lower triangular.
  *
  *  This is the block version of the algorithm, calling Level 3 BLAS.
  *
  *  Arguments
  *  =========
  *
  *  UPLO    (input) CHARACTER*1
  *          = 'U':  Upper triangle of A is stored;
  *          = 'L':  Lower triangle of A is stored.
  *
  *  N       (input) INTEGER
  *          The order of the matrix A.  N >= 0.
  *
  *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  *          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
  *          N-by-N upper triangular part of A contains the upper
  *          triangular part of the matrix A, and the strictly lower
  *          triangular part of A is not referenced.  If UPLO = 'L', the
  *          leading N-by-N lower triangular part of A contains the lower
  *          triangular part of the matrix A, and the strictly upper
  *          triangular part of A is not referenced.
  *
  *          On exit, if INFO = 0, the factor U or L from the Cholesky
  *          factorization A = U**T*U or A = L*L**T.
  *
  *  LDA     (input) INTEGER
  *          The leading dimension of the array A.  LDA >= max(1,N).
  *
  *  INFO    (output) INTEGER
  *          = 0:  successful exit
  *          < 0:  if INFO = -i, the i-th argument had an illegal value
  *          > 0:  if INFO = i, the leading minor of order i is not
  *                positive definite, and the factorization could not be
  *                completed.
  *
  *  =====================================================================


    */

template<typename T>
class M_Matrix
{
public:

  class MyConstIterator;
  enum ITER_TYPE{FULL_IT,LT_IT,UT_IT,DIAG_IT,SCALAR_IT, ZERO_IT};

  enum TYPE{ZERO,FULL,SYMMETRIC,DIAGONAL,SCALAR_FULL,SCALAR_DIAGONAL};
  static
  ITER_TYPE
  iter_type(TYPE type)
  {
    switch(type)
      {
      case ZERO:
        return ZERO_IT;
      case FULL:
        return FULL_IT;
      case SYMMETRIC:
        return LT_IT;
      case DIAGONAL:
        return DIAG_IT;
      case SCALAR_DIAGONAL:
      case SCALAR_FULL:
        return SCALAR_IT;
      default: return {};
      }
  }
  static
  ITER_TYPE const_iter_type(TYPE type)
  {
    switch(type)
      {
      case FULL:
      case SYMMETRIC:
      case SCALAR_FULL:
        return FULL_IT;
      case DIAGONAL:
      case SCALAR_DIAGONAL:
        return DIAG_IT;
      case ZERO:
      default:
        return ZERO_IT;
      }
  }

  class MyIterator : public std::iterator<std::input_iterator_tag, T>
  {
    M_Matrix<T>& m;
    std::size_t i_;
    std::size_t j_;
    ITER_TYPE type_;
  public:

    friend class MyConstIterator;
    ITER_TYPE type()const {return type_;}
    std::size_t iRow()const{return i_;}
    std::size_t jCol()const{return j_;}

    MyIterator(M_Matrix<T>& x) :m(x),i_(0),j_(0),type_(iter_type(x.type())) {}

    MyIterator(M_Matrix<T>& x, ITER_TYPE t) :m(x),i_(0),j_(0),type_(t) {}

    MyIterator(M_Matrix<T>& x, std::size_t i, std::size_t j)
      :m(x),i_(i),j_(j),type_(iter_type(x.type()))  {}

    MyIterator(const MyIterator& mit)
      :m(mit.m), i_(mit.i_),j_(mit.j_),type_(mit.type_){}
    MyIterator& operator++()
    {
      ++j_;
      switch(type())
        {
        case FULL_IT:
          if (j_>=m.ncols())
            {
              j_=0;
              ++i_;
            }
          break;
        case LT_IT:
          if (j_>i_)
            {
              j_=0;
              ++i_;
            }
          break;
        case UT_IT:
          if (j_>=m.ncols())
            {
              ++i_;
              if (i_>=m.nrows())
                j_=0;
              else
                j_=i_;
            }
          break;


        case DIAG_IT:
          if (j_>=m.ncols())
            {
              j_=0;
              i_=m.nrows();
            }
          else
            {
              i_=j_;
            }
          break;
        case SCALAR_IT:
        case ZERO_IT:
          if (j_>0)
            {
              j_=0;
              i_=m.nrows();
            }
          break;


        }
      return *this;
    }
    MyIterator operator++(int)
    {MyIterator tmp(*this); operator++(); return tmp;}
    bool operator==(const MyIterator& rhs)
    {
      if (i_!=rhs.i_)
        return false;
      else if (j_!=rhs.j_)
        return false;
      else return true;
    }
    bool operator!=(const MyIterator& rhs) {return ! (*this==rhs);}
    T& operator*() {return m(i_,j_);}
    const T& operator*() const {return m(i_,j_);}
  };


  class MyConstIterator : public std::iterator<std::input_iterator_tag, T>
  {
    const M_Matrix<T>& m;
    std::size_t i_;
    std::size_t j_;
    ITER_TYPE type_;
  public:
    std::size_t iRow()const{return i_;}
    std::size_t jCol()const{return j_;}
    ITER_TYPE type()const {return type_;}

    MyConstIterator(const M_Matrix<T>& x)
      :m(x),i_(0),j_(0),type_(const_iter_type(x.type())) {}

    MyConstIterator(const M_Matrix<T>& x, ITER_TYPE t)
      :m(x),i_(0),j_(0),type_(t) {}

    MyConstIterator(const M_Matrix<T>& x, std::size_t i, std::size_t j)
      :m(x),i_(i),j_(j), type_(const_iter_type(x.type())) {}

    MyConstIterator(const MyConstIterator& mit)
      : m(mit.m),i_(mit.i_),j_(mit.j_),
        type_(mit.type_){}

    MyConstIterator(const MyIterator& mit)
      : m(mit.m),i_(mit.i_),j_(mit.j_),type_(mit.type_) {}

    MyConstIterator& operator++()
    {
      ++j_;
      switch(type())
        {
        case FULL_IT:
          if (j_>=m.ncols())
            {
              j_=0;
              ++i_;
            }
          break;
        case LT_IT:
          if (j_>i_)
            {
              j_=0;
              ++i_;
            }
          break;
        case UT_IT:
          if (j_>=m.ncols())
            {
              ++i_;
              if (i_>=m.nrows())
                j_=0;
              else
                j_=i_;
            }
          break;


        case DIAG_IT:
          if (j_>=m.ncols())
            {
              j_=0;
              i_=m.nrows();
            }
          else
            {
              i_=j_;
            }
          break;
        case SCALAR_IT:
        case ZERO_IT:
          if (j_>0)
            {
              j_=0;
              i_=m.nrows();
            }
          break;


        }
      return *this;
    }
    MyConstIterator operator++(int)
    {MyConstIterator tmp(*this); operator++(); return tmp;}
    bool operator==(const MyConstIterator& rhs)
    {
      if (i_!=rhs.i_)
        return false;
      else if (j_!=rhs.j_)
        return false;
      else return true;
    }
    bool operator!=(const MyConstIterator& rhs) {return ! (*this==rhs);}
    const T& operator*() const {return m(i_,j_);}
  };



  typedef  MyIterator iterator;
  typedef  MyConstIterator const_iterator;


  iterator begin()
  {
    if (type()!=ZERO)
      {
        MyIterator out(*this);
        return out;
      }
    else
      return end();
  }


  iterator begin(ITER_TYPE t)
  {
    if ((type()!=ZERO)&&(t!=ZERO_IT))
      {
        MyIterator out(*this,t);
        return out;
      }
    else
      return end();
  }

  iterator end()
  {
    MyIterator out(*this,nrows(),0);
    return out;
  }

  const_iterator begin()const
  {
    if (type()!=ZERO)
      {
        MyConstIterator out(*this);
        return out;
      }
    else
      return end();
  }

  const_iterator begin(ITER_TYPE t)const
  {
    if ((type()!=ZERO)&&(t!=ZERO_IT))
      {
        const_iterator out(*this,t);
        return out;
      }
    else
      return end();
  }
  const_iterator end() const
  {
    MyConstIterator out(*this,nrows(),0);
    return out;
  }


  M_Matrix()=default;

  M_Matrix (const M_Matrix<T> & sample)=default;
  M_Matrix (M_Matrix<T> && sample)=default;


  M_Matrix (std::size_t nrows,std::size_t ncols)
    : type_(FULL),
      _nrows(nrows),
      _ncols(ncols),
      _data(nrows*ncols)
  {
  }


  template<typename S>
  M_Matrix (const M_Matrix<S> & sample)
    : type_(TYPE(sample.type())),
      _nrows(sample.nrows()),
      _ncols(sample.ncols()),
      _data(sample.size())

  {
    for (auto it=begin(); it!=end(); ++it)
      *it = sample(it);

  }
  M_Matrix (const std::vector<std::vector<T>> & sample)
    : type_(FULL),
      _nrows(sample.size()),
      _ncols(sample[0].size()),
      _data(sample.size()*sample[0].size())
  {
    for (std::size_t i = 0; i < sample.size(); ++i)
      for (std::size_t j=0; j<sample[0].size(); ++j)
        (*this)(i,j) = sample[i][j];

  }

  M_Matrix (TYPE t,const std::vector<std::vector<T>> & sample)
    : type_(t),
      _nrows(sample.size()),
      _ncols(sample[0].size()),
      _data(getSize(t,sample.size(),sample[0].size()))
  {
    for (auto it=begin(); it!=end(); ++it)
      (*it)=sample[it.iRow()][it.jCol()];

  }

  M_Matrix (TYPE t,const M_Matrix<T> & sample)
    : type_(t),
      _nrows(sample.nrows()),
      _ncols(sample.ncols()),
      _data(getSize(t,sample.nrows(),sample.ncols()))
  {
    for (auto it=begin(); it!=end(); ++it)
      (*it)=sample(it);

  }

  template<typename S>
  M_Matrix (const std::vector<std::vector<S>> & sample)
    :type_(FULL),
      _nrows(sample.size()),
      _ncols(sample[0].size()),
      _data(sample.size()*sample[0].size())
  {
    for (std::size_t i = 0; i < sample.size(); ++i)
      (*this)[i] = sample[i];

  }

  M_Matrix(std::size_t nrows_,std::size_t ncols_,TYPE t,T data):
    type_(t),
    _nrows(nrows_),
    _ncols(ncols_),
    _data(getSize(t,nrows_,ncols_),data)
  {
  }

  M_Matrix(std::size_t nrows_,std::size_t ncols_,TYPE t):
    type_(t),
    _nrows(nrows_),
    _ncols(ncols_),
    _data(getSize(t,nrows_,ncols_))
  {
  }




  M_Matrix(std::size_t nrows_,std::size_t ncols_, std::vector<T> data):
    type_(FULL),
    _nrows(nrows_),
    _ncols(ncols_),
    _data(data)
  {

  }

  M_Matrix(std::size_t nrows_,std::size_t ncols_, const M_Matrix<T>& data):
    type_(FULL),
    _nrows(nrows_),
    _ncols(ncols_),
    _data(nrows_*ncols_)
  {
    assert(data.nrows()==nrows());
    assert(data.ncols()==ncols());
    for (size_t i=0; i<nrows(); ++i)
      for (std::size_t j=0; j<ncols(); ++j)
        (*this)(i,j)=data(i,j);
  }

  template <typename S>
  M_Matrix(std::size_t nrows_,std::size_t ncols_, const M_Matrix<S>& data):
    type_(FULL),
    _nrows(nrows_),
    _ncols(ncols_),
    _data(nrows_*ncols_)
  //   _data(new T[_ncells])
  {

    assert(data.nrows()==nrows());
    assert(data.ncols()==ncols());
    for (size_t i=0; i<nrows(); ++i)
      for (std::size_t j=0; j<ncols(); ++j)
        (*this)(i,j)=data(i,j);
  }


  template <typename S>
  M_Matrix(std::size_t nrows_,std::size_t ncols_, TYPE type,const M_Matrix<S>& data):
    type_(type),
    _nrows(nrows_),
    _ncols(ncols_),
    _data(nrows_*ncols_)
  //   _data(new T[_ncells])
  {

    assert(data.nrows()==nrows());
    assert(data.ncols()==ncols());
    for (auto it=begin(); it!=end(); ++it)
      (*it)=data(it);
  }


  M_Matrix(std::size_t nrows_,std::size_t ncols_,T data):
    type_(FULL),
    _nrows(nrows_),
    _ncols(ncols_),
    //   _data(new T[_ncells])
    _data(nrows_*ncols_,data)
  {

  }


  M_Matrix& operator= (const M_Matrix<T> & sample)=default;
  M_Matrix& operator= (M_Matrix<T> && sample)=default;


  ~M_Matrix()
  {
    //  if (_data>0)
    //	delete [] _data;
  }



  template<class F>
  M_Matrix<T>
  apply(const F& f)const
  {
    M_Matrix<T> out(nrows(),ncols(),type());
    for (std::size_t i=0; i<size(); ++i)
      out[i]=f((*this)[i]);
    return out;
  }




  M_Matrix<T>& operator=(T X)
  {
    for (std::size_t i=0; i<size(); ++i) _data[i]=X;
    return *this;
  }



  size_t size()const
  {
    return _data.size();
  }

  size_t nrows()const
  {
    return _nrows;
  }

  size_t ncols()const
  {
    return _ncols;
  }



  T& operator[](std::size_t n)
  {
    return _data[n];
  }
  bool empty()const
  {
    return _data.empty();
  }

  const T& operator[](std::size_t n) const
  {
    return _data[n];
  }


  T&  operator() (std::size_t i,std::size_t j)
  {
    assert(i<nrows());
    assert(j<ncols());
    switch(type_){
      case FULL:
        return (*this)[i*_ncols+j];
      case SYMMETRIC:
        if(i>=j)
          return (*this)[(i*(i+1))/2+j]; //stores at lower triangular portion
        else
          return (*this)[(j*(j+1))/2+i];
      case DIAGONAL:
        assert(i==j);
        return (*this)[i];
      case SCALAR_DIAGONAL:
        assert(i==j);
        return (*this)[0];
      case SCALAR_FULL:
        return (*this)[0];
      case ZERO:
        assert(false);
        return (*this)[0];
      default:
        assert(false);
        return (*this)[0];
      }
  }

  T const&  operator() (std::size_t i,std::size_t j) const
  {
    if(i>=nrows())
      assert(true);
    assert(j<ncols());
    switch(type_){
      case FULL:
        return (*this)[i*_ncols+j];
      case SYMMETRIC:
        if(i>=j)
          return (*this)[i*(i+1)/2+j]; //stores at lower triangular portion
        else
          return (*this)[j*(j+1)/2+i];
      case DIAGONAL:
        if(i==j)
          return (*this)[i];
        else
          return zero_;
      case SCALAR_DIAGONAL:
        if(i==j)
          return (*this)[0];
        else
          return zero_;
      case SCALAR_FULL:
        return (*this)[0];
      case ZERO:
        return zero_;
      default:
        return zero_;
      }
  }




  /** @name  Accesing all the values of a Row or Column at once
     */
  //@{

  /**
    Replacement of the ith Row.
    @param iRow the extracted row (0 is the first as always in C)
    @param newValues contains the numbers used to relace iRow
    @param dummy an ignored variable to differentiate from column
    extraction
    @pre iRow is smaller than nrows(*this)
    @pre size(newValues)==ncols(*this)
    @returns *this
    @post if the preconditions are met, the ith row is replaced by
    newValues
    @post newValues is treated as a vector, although is a Matrix. Its
    internal structure (i.e., ncols and nrows) is ignored.
    @post assert the precoditions
    */

  M_Matrix<T>&  operator() (std::size_t iRow,
                            std::string /*dummy*/,
                            const M_Matrix<T>& newValues)
  {
    assert(type_==FULL);
    assert(iRow<nrows());//number of rows
    assert(newValues.size()==ncols()); //number of columns
    for (std::size_t j=0; j<std::min(ncols(),size()); j++)
      this->operator()(iRow,j)=newValues[j];
    return *this;
  }





  /**
    Replacement of the jth Column.
    @param newValues contains the numbers used to relace jth Column
    @pre newValues is treated as a vector, although is a Matrix. Its
    internal structure (i.e., ncols and nrows) is ignored.
    @param jColumn the replaced column (0 is the first as always in C)
    @param dummy an ignored variable to differentiate from column
    extraction
    @pre jColumn is smaller than ncols(*this)
    @pre size(newValues)==nrows(*this)
    \returns *this
    @post if the preconditions are met, the jth Column is replaced by
    newValues
    @post assert the precoditions
    */

  template <typename V>
  M_Matrix<T>&  operator() (const std::string /*dummy*/,
                            std::size_t jColumn,
                            const V& newValues)
  {
    assert(type_==FULL);

    for (std::size_t i=0; i<std::min(nrows(),newValues.size()); i++)
      this->operator()(i,jColumn)=newValues[i];
    //  assert(ndim>1);
    //  assert(i<n[0]);//number of rows
    //  assert(j<n[1]); //number of columns
    return *this;
  }





  /**
    Copy of the ith Row
    @pre iRow is smaller than nrows(*this)
    @param iRow the extracted row (0 is the first as always in C)
    @param dummy an ignored variable to differentiate from column
    extraction
    \returns a 1-row ncols Matrix with the values of the ith row
    */

  M_Matrix<T>  operator() (std::size_t iRow,
                           const std::string /*dummy*/
                           ) const
  {
    assert(type_==FULL);
    M_Matrix<T> out(1,ncols());
    for (std::size_t j=0; j<ncols(); j++)
      out[j]=this->operator()(iRow,j);
    return out;
  }



  template <class Iter>
  T& operator() (const Iter &it)
  {
    return (*this)(it.iRow(), it.jCol());
  }
  template <class Iter>
  const T& operator()(const Iter& it) const
  {
    return (*this)(it.iRow(), it.jCol());
  }


  /**
    Copy of the jth Column
    @pre jColumn is smaller than ncols(*this)
    @param jColumn the extracted column (0 is the first as always in C)
    @param dummy is an ignored const string (like "") to differentiate
    from row extraction
    \returns a nrows 1-column Matrix with the values of the jth column
    */

  M_Matrix<T>  operator() (std::string /*dummy*/,
                           std::size_t jColumn
                           ) const
  {
    assert(type_==FULL);
    M_Matrix<T> out(nrows(),1);
    for (std::size_t i=0; i<nrows(); i++)
      out[i]=(*this)(i,jColumn);
    return out;
  }


  void clear()
  {
    type_=ZERO;
    _nrows=0;
    _ncols=0;
    _data.clear();
  }


  std::vector<T> toVector()const
  {
    return _data;
  }

  M_Matrix<T> toVector_of_Rows()const
  {
    return M_Matrix<T>(size(),1,_data);
  }
  M_Matrix<T> toVector_of_Cols()const
  {
    return M_Matrix<T>(1,size(),_data);
  }


  std::vector<std::vector<T>> toMatrix()const
  {
    std::vector<std::vector<T>> out(nrows(),std::vector<T>(ncols()));
    for (std::size_t i=0;i<nrows();++i)
      for (std::size_t j=0;j<ncols();++j)
        out[i][j]=(*this)(i,j);
    return out;
  }

  /**
     Returns a custom sized Matrix filled with ones
    @post (ones(n,m))(i,j)==T(1)
     */



  static
  M_Matrix
  unpackForLapack(const M_Matrix<double>& packedMatrix, char UPLO='L')
  {
    switch (packedMatrix.type())
      {
      case SYMMETRIC:
        {
          M_Matrix<double> out(packedMatrix.nrows(),packedMatrix.ncols());
          if (UPLO=='L')
            {
              for (std::size_t i=0; i<out.nrows(); ++i)
                for (std::size_t j=0; j<=i; ++j)
                  out(i,j)=packedMatrix(i,j);
              return out;
            }
          else
            {
              for (std::size_t i=0; i<out.nrows(); ++i)
                for (std::size_t j=i; j<=out.ncols(); ++j)
                  out(i,j)=packedMatrix(i,j);
              return out;
            }

        }
      default:
        return M_Matrix<double>(packedMatrix);
      }
  }


  M_Matrix
  full()const
  {
    M_Matrix out(nrows(),ncols());
    for (std::size_t i=0; i<nrows(); i++)
      for (std::size_t j=0; j<ncols(); ++j)
        out(i,j)=(*this)(i,j);
    return out;
  }

  static
  std::size_t getSize(TYPE t,std::size_t nrows,std::size_t ncols)
  {
    switch (t)
      {
      case ZERO: return 0;
      case FULL: return nrows*ncols;
      case SYMMETRIC:
        {
          assert(nrows==ncols);
          return (nrows*(ncols+1))/2;
        }
      case DIAGONAL:
        {
          return std::min(nrows,ncols);
          break;
        }
      case SCALAR_FULL:
      case SCALAR_DIAGONAL:
      default:
        return 1;
      }
  }
  std::size_t getSize(TYPE t,std::size_t nrows)
  {
    switch (t)
      {
      case ZERO: return 0;
      case FULL: return nrows*nrows;
      case SYMMETRIC:
        {
          return (nrows*(nrows+1))/2;
        }
      case DIAGONAL:
        {
          return nrows;
          break;
        }
      case SCALAR_FULL:
      case SCALAR_DIAGONAL:
      default:
        return 1;
      }
  }


  TYPE type()const {return type_;}


  bool isSymmetric()const {
    switch (type_)
      {
      case ZERO:
      case SYMMETRIC:
        return true;
      case SCALAR_DIAGONAL:
      case DIAGONAL:
      case SCALAR_FULL:
        return ncols()==nrows();
      case FULL:
      default:
        return false;
      }
  }

  bool isDiagonal()const {
    switch(type())
      {
      case FULL:
      case SYMMETRIC:
      case SCALAR_FULL:
        return false;
      case DIAGONAL:
      case SCALAR_DIAGONAL:
      case ZERO:
      default:
        return true;
      }
  }

  /**
    Matrix of zeros with the shape of the provided Matrix
     @post nrows(ones(x))==x.nrows()   ncols(ones(x))=x.ncols()

    */

  template<typename S>
  friend
  auto
  TranspMult(const M_Matrix<T>& x,const M_Matrix<S>& y)->
  M_Matrix<decltype(operator*(std::declval<T>(),std::declval<S>()))>;
  M_Matrix(TYPE t, std::size_t nrows, std::size_t ncols, std::vector<T> data):
    type_(t),_nrows(nrows),_ncols(ncols),_data(data){}

  template<typename...Ts , template<typename...>class V>
  M_Matrix(std::size_t nrows, std::size_t ncols,TYPE t,  V<Ts...> data):
    type_(t),_nrows(nrows),_ncols(ncols),_data(data.size()){
    for (std::size_t i=0; i<data.size(); ++i)
      _data[i]=data[i];
  }

  /*!
       Matrix Addition assignment.
       @returns a reference to itself
       @pre  same number of rows and columns
       @post all the values of the matrix are summed up by the corresponing
             values of the other matrix\n\n
             assert(nrows(*this)==nrows(other))&& assert(ncols(*this)==ncols(other))

       */



private:
  TYPE type_;
  std::size_t          _nrows;    /**< number of rows */
  std::size_t          _ncols;    /**< number of columns */

  /** internal data
          @remarks can be a pointer or a vector
          @remarks now is a vector for debugging purposes
          */
  //	T                      *_data; /**< pointer to the data  */
  std::vector<T>        _data;
  T zero_=T{};
};

template<typename T>
M_Matrix<T> operator-(const M_Matrix<T>& x)
{
  M_Matrix<T> out(x.nrows(), x.ncols(), x.type());
  for (auto it=out.begin(); it!=out.end(); ++it)
    *it= -x(it);
  return out;
}

template<typename T>
M_Matrix<T> operator-(M_Matrix<T>&& x)
{
  for (auto it=x.begin(); it!=x.end(); ++it)
    *it=-*it;
  return x;
}

namespace Vector_Binary_Transformations
{

  template <typename T>
  std::vector<T> operator -(const std::vector<T>& x,const std::vector<T>& y)
  {
    std::vector<T> out(x.size());
    for (std::size_t i=0; i<x.size(); ++i)
      out[i]=x[i]-y[i];
    return out;
  }

  template <typename T>
  std::vector<T> operator +(const std::vector<T>& x,const std::vector<T>& y)
  {
    std::vector<T> out(x.size());
    for (std::size_t i=0; i<x.size(); ++i)
      out[i]=x[i]+y[i];
    return out;
  }


  template <typename T>
  std::vector<T> elemMult(const std::vector<T>& x,const std::vector<T>& y)
  {
    std::vector<T> out(x.size());
    for (std::size_t i=0; i<x.size(); ++i)
      out[i]=x[i]*y[i];
    return out;
  }

  template <typename T>
  std::vector<T> elemDiv(const std::vector<T>& x,const std::vector<T>& y)
  {
    std::vector<T> out(x.size());
    for (std::size_t i=0; i<x.size(); ++i)
      out[i]=x[i]/y[i];
    return out;
  }

}
namespace Vector_Binary_Transformations
{

  inline
  std::vector<double> operator*(const std::vector<double>& x,
                                const M_Matrix<double>& y)
  {
    std::vector<double> out(y.ncols(),0);
    for (std::size_t i=0; i<x.size(); ++i)
      for (std::size_t j=0; j<y.ncols(); ++j)
        out[j]+=x[i]*y(i,j);
    return out;
  }



  inline
  double xTSigmaX(const M_Matrix<double> &vector, const M_Matrix<double> &matrix)
  {
    double sum=0;
    if (matrix.isDiagonal())
      {
        for (std::size_t i=0; i<matrix.nrows(); ++i)
          sum+=vector[i]*matrix(i,i)*vector[i];

      }
    else
      {
        for (std::size_t i=0; i<matrix.nrows(); ++i)
          {
            sum+=vector[i]*matrix(i,i)*vector[i];
            for (std::size_t j=i+1; j<matrix.ncols();++j)
              sum+=2*vector[i]*matrix(i,j)*vector[j];
          }
      }
    return sum;
  }

  template<class T>
  double xTSigmaX(const M_Matrix<T> &vector, const M_Matrix<T> &matrix)
  {
    double sum=0;
    if (matrix.isDiagonal())
      {
        for (std::size_t i=0; i<matrix.nrows(); ++i)
          sum+=xTSigmaX(vector[i],matrix(i,i));
      }
    else
      {
        for (std::size_t i=0; i<matrix.nrows(); ++i)
          {
            sum+=xTSigmaX(vector[i],matrix(i,i));
            for (std::size_t j=i+1; j<matrix.ncols();++j)
              sum+=2*xTSigmaX(vector[i],matrix(i,j));
          }
      }
    return sum;
  }



  inline
  double xTSigmaX(const std::vector<double> &v, const M_Matrix<double> &matrix)
  {
    double sum=0;
    for (std::size_t i=0; i<matrix.nrows(); ++i)
      {
        sum+=v[i]*matrix(i,i)*v[i];
        for (std::size_t j=i+1; j<matrix.ncols();++j)
          sum+=2*v[i]*matrix(i,j)*v[j];
      }
    return sum;
  }


}

namespace lapack
{

  extern "C" void dpotrf_(char * 	UPLO,
                          int * N,
                          double * A,
                          int * LDA,
                          int * INFO);


  inline M_Matrix<double> UT(const M_Matrix<double>& x)
  {
    M_Matrix<double> y(x.nrows(),x.ncols());
    for (std::size_t i=0;i<x.nrows();i++)
      {
        for (std::size_t j=0; j<i; j++)
          y(i,j)=0;
        for (std::size_t j=i;j<x.ncols(); j++)
          y(i,j)=x(i,j);

      }
    return y;
  }

  inline M_Matrix<double> LT(const M_Matrix<double>& x)
  {
    M_Matrix<double> y(x.nrows(),x.ncols());
    for (std::size_t i=0;i<x.nrows();i++)
      {
        for (std::size_t j=0;j<i+1; j++)
          y(i,j)=x(i,j);
        for (std::size_t j=i+1;j<x.ncols(); j++)
          y(i,j)=0;
      }
    return y;
  }



}


namespace Matrix_Binary_Transformations
{
  template<typename T>
  M_Matrix<T> operator+(const M_Matrix<T>& x,const M_Matrix<T>& y);
  template<typename T>
  M_Matrix<T> operator-(const M_Matrix<T>& x,const M_Matrix<T>& y);
  template<typename T>
  M_Matrix<T> operator-(M_Matrix<T>&& x,M_Matrix<T>&& y);

  template<typename T, typename S>
  M_Matrix<T>& operator+=(M_Matrix<T>& x,const M_Matrix<S>& y);

  template<typename T, typename S>
  M_Matrix<T>& operator-=(M_Matrix<T>& x,const M_Matrix<S>& y);

  template<typename T, typename S>
  auto
  operator *
  (const M_Matrix<T>& one,
   const M_Matrix<S>& other)
  ->M_Matrix<decltype(std::declval<T>()*std::declval<S>())>;

  template<typename T, typename S>
  auto
  quadraticForm_BT_A_B (const M_Matrix<T>& A, const M_Matrix<S>& B)
  ->M_Matrix<decltype(std::declval<T>()*std::declval<S>())>;

  template<typename T, typename S>
  auto
  quadraticForm_B_A_BT (const M_Matrix<T>& A, const M_Matrix<S>& B)
  ->M_Matrix<decltype(std::declval<T>()*std::declval<S>())>;


  inline
  M_Matrix<double>
  Forward_Sustitution_Ly_b(const M_Matrix<double>& L,
                           const M_Matrix<double>& b);

  template<typename T>
  M_Matrix<T>
  Forward_Sustitution_Ly_b(const M_Matrix<T>& L, const M_Matrix<T>& b);


}

using namespace Matrix_Binary_Transformations;

namespace Matrix_Generators
{
  template<typename T>
  M_Matrix<T>  ones(size_t nrows, size_t ncols)
  {
    return M_Matrix<T>(nrows,ncols,M_Matrix<T>::SCALAR_FULL,T(1));
  }

  template<typename T>
  M_Matrix<T>  zeros(size_t nrows, size_t ncols)
  {
    return M_Matrix<T>(nrows,ncols,M_Matrix<T>::ZERO);
  }

  /**
  Identity Matrix of the specified size
 */
  template<typename T>
  M_Matrix<T> eye(std::size_t n)
  {
    return M_Matrix<T>(n,n,M_Matrix<T>::SCALAR_DIAGONAL,T(1));
  }

  template<typename T>
  M_Matrix<T>  Rand(const M_Matrix<T>& x)
  {
    std::normal_distribution<> normal;
    std::random_device rd;
    std::mt19937_64 sto(rd());
    auto out=M_Matrix<T>(x);
    for (std::size_t i=0; i<out.size(); ++i)
      out[i]=normal(sto);
    return out;
  }

  template<class D>
  M_Matrix<double>  Rand(const M_Matrix<double>& x,D d, std::mt19937_64& sto)
  {
    auto out=M_Matrix<double>(x.nrows(),x.ncols(),x.type());
    for (auto it=out.begin(); it!=out.end(); ++it)
      *it=d(sto);
    return out;
  }


  template<typename T, class D>
  M_Matrix<T>  Rand(const M_Matrix<T>& x,D d, std::mt19937_64& sto)
  {
    auto out=M_Matrix<T>(x.nrows(),x.ncols(),x.type());;
    for (auto it=out.begin(); it!=out.end(); ++it)
      *it=Rand(*it,d,sto);
    return out;
  }
}

namespace Matrix_Stream_Operators
{

  template<typename T>
  std::ostream& operator<<(std::ostream& os,const M_Matrix<T>& x)
  {
    os<<"[";
    for (std::size_t i=0; i<x.nrows(); ++i)
      {
        for (std::size_t j=0; j<x.ncols(); ++j)
          os<<x(i,j)<<" ";
        os<<";";
      }
    os<<"]";
    return os;
  }


  template<typename T>
  std::istream& operator>>(std::istream& is,M_Matrix<T>& x)
  {
    std::vector<T> o;
    std::size_t nrows=0;
    char ch;
    while ((is>>ch)&&(ch!='[')){}
    if(ch!='[')
      return is;
    else
      while (ch!=']')
        {
          std::string s;
          while ((is.get(ch))&&((ch!=']')&&ch!=';'))
            {
              s.push_back(ch);
            }
          std::stringstream ss(s);
          T e;
          std::size_t i=o.size();
          while (ss>>e) o.push_back(e);
          if (o.size()>i) ++nrows;
        }
    std::size_t ncols=o.size()/nrows;
    x=M_Matrix<T>(nrows,ncols,o);
    return is;

  }
}


namespace Matrix_Unary_Predicates
{
  template<typename T, class Predicate>
  bool all(const M_Matrix<T>& x, const Predicate& p)
  {
    for (std::size_t i=0; i< x.size(); ++i)
      if (!p(x[i]))
        return false;
    return true;
  }

  template<typename T, class Predicate>
  bool any(const M_Matrix<T>& x, const Predicate& p)
  {
    for (std::size_t i=0; i< x.size(); ++i)
      if (p(x[i]))
        return true;
    return false;
  }

  using std::isnan;

  template <typename T>
  bool isnan(const M_Matrix<T>& x)
  {
    return any(x,[](const T& e){return isnan(e);});
  }
}


namespace Matrix_Unary_Size_Functions
{

  template<typename T>

  std::size_t nrows(const M_Matrix<T>& x){return x.nrows();}

  template<typename T>
  std::size_t ncols(const M_Matrix<T>& x){return x.ncols();}

  template<typename T>
  std::size_t size(const M_Matrix<T>& x){return x.size();}

}
using namespace Matrix_Unary_Size_Functions;
namespace Matrix_Unary_Transformations
{
  template<class E,class F, typename T>
  M_Matrix<T>
  accumulate_by_Rows(const M_Matrix<E>& me,const F& f, const T& init);
  template<class E,class F,typename T>
  M_Matrix<T>
  accumulate_by_Cols(const M_Matrix<E>& me,const F& f, const T& start);
  template<typename T>
  M_Matrix<T> diag(const M_Matrix<T>& x);

  template<typename T>
  M_Matrix<T> Transpose(const M_Matrix<T>& x);


}
using namespace Matrix_Unary_Transformations;

namespace Matrix_Unary_Functions
{


  template<class E, class F, typename T>
  T
  accumulate(const M_Matrix<E> & x, const F& f, const T& start)
  {
    T out(start);
    for (auto it=x.begin(); it!=x.end(); ++it)
      out=f(out,*it);
    return out;
  }

  inline
  double
  fullSum(const double& x)
  {
    return x;
  }
  template <class T>
  double
  fullSum(const M_Matrix<T>& x)
  {
    return accumulate
        (x,[](double a, const T& b)
    {return a+fullSum(b);},0.0);
  }

  inline
  double
  logProduct(double x)
  {
    return log(x);
  }
  template <class T>
  double
  logProduct(const M_Matrix<T>& x)
  {
    return accumulate
        (x,[](double a, const T& b)
    {return a+logProduct(b);},0.0);
  }

  inline
  double maxAbs(double x)
  {
    return std::abs(x);
  }

  template<typename T>
  double maxAbs(const M_Matrix<T>& x)
  {
    return accumulate(x,[](double a, const T& b)
    {return std::max(a,maxAbs(b));},0.0);
  }

  template<typename T>
  double logDiagProduct(const M_Matrix<T>& x)
  {
    auto d=diag(x);
    return logProduct(d);
  }


  /**
       Product of the Diagonal of a Matrix

      */
  template<typename T>
  double  diagProduct(const M_Matrix<T>& x)
  {
    return std::exp(logDiagProduct(x));
  }


  template<typename T>
  double Trace(const M_Matrix<T>& x)
  {
    auto d=diag(x);
    return fullSum(d);
  }


  /**
    Maximum Value of the Sum of the absolute values in each row
   */
   template<typename T>
  T norm_inf(const M_Matrix<T>& x)
  {
  T n(0);
  for (size_t i=0; i<nrows(x); ++i)
  {
      T sum(0);
      for (size_t j=0; j<ncols(x); ++j)
      if (x(i,j)>0)
          sum+=x(i,j);
      else
      sum-=x(i,j);

      n=std::max(n,sum);
  }
  return n;
  }



  /**
  Returns the exponential of a Matrix
  @remarks  uses the Pade approximation
  @pre the matrix has to be square, ncols()==nrows
  @post assert the precondition

  */

  template<typename T>
  M_Matrix<T>  expm(const M_Matrix<T>& x)
  {
      assert(ncols(x)==nrows(x));
      assert(size(x)>0);
      using namespace Matrix_Generators;

      // Scale A by power of 2 so that its norm is < 1/2 .
      double e = std::ceil(log(norm_inf(x))/log(2));
      std::size_t s = std::max(0.0,e+1);

      M_Matrix<T> A = x/std::pow(2.0,int(s));



      // Pade approximation for exp(A)
      M_Matrix<T> X = A;
      double c = 0.5;
      M_Matrix<T> E = eye<T>(nrows(A)) + c*A;
      M_Matrix<T> D = eye<T>(nrows(A)) - c*A;


      std::size_t q = 6;
      bool p = true;
      for (std::size_t k = 2; k<q+1;++k)
      {
      c = c * (1.0*q-1.0*k+1) / (1.0*k*(2.0*q-k+1.0));
      X = A*X;
      M_Matrix<T> cX = c*X;
      E += cX;
      if (p)
          D += cX;
      else
          D -= cX;
      p = !p;
      }
      E = invSafe(D)*E;

      // Undo scaling by repeated squaring
      for (std::size_t k=0; k<s;k++)
       E = E*E;
      return E;

  }




}



namespace Matrix_Unary_Transformations
{
  using Matrix_Binary_Transformations::
  quadraticForm_BT_A_B;

  using Matrix_Binary_Transformations::
  quadraticForm_B_A_BT;

  template<class E,class F, typename T>
  M_Matrix<T>
  accumulate_by_Rows(const M_Matrix<E>& me,const F& f, const T& start)
  {
    M_Matrix<T> out(me.nrows(),1);
    switch(me.type())
      {
      case M_Matrix<T>::FULL:
      case M_Matrix<T>::SYMMETRIC:
      case M_Matrix<T>::SCALAR_FULL:
        {
          for (std::size_t i=0; i<me.nrows(); ++i)
            {
              T s=start;
              for (std::size_t j=0; j<me.ncols(); ++j)
                s=f(s,me(i,j));
              out(i,1)=s;
            }
          return out;
        }

      case M_Matrix<T>::DIAGONAL:
      case M_Matrix<T>::SCALAR_DIAGONAL:
        {
          for (std::size_t i=0; i<me.nrows(); ++i)
            out(i,1)=f(start,me(i,i));
          return out;
        }
      case M_Matrix<T>::ZERO:
        {
          return M_Matrix<T>(me.nrows(),1,f(start,me[0]));
        }
      };
  }

  template<class E,class F,typename T>
  M_Matrix<T>
  accumulate_by_Cols(const M_Matrix<E> &me, const F& f, const T& start)
  {
    M_Matrix<T> out(1, me.ncols());
    switch(me.type())
      {
      case M_Matrix<E>::FULL:
      case M_Matrix<E>::SYMMETRIC:
      case M_Matrix<E>::SCALAR_FULL:
        {
          for (std::size_t j=0; j<me.ncols(); ++j)
            {
              T s(start);
              for (std::size_t i=0; i<me.nrows(); ++i)
                s=f(s,me(i,j));
              out(1,j)=s;
            }
          return out;
        }

      case M_Matrix<E>::DIAGONAL:
      case M_Matrix<E>::SCALAR_DIAGONAL:
        {
          for (std::size_t i=0; i<me.ncols(); ++i)
            out(1,i)=f(start,me(i,i));
          return out;
        }
      case M_Matrix<E>::ZERO:
        {
          return M_Matrix<T>(1, me.ncols(), f(start,me[0]));
        }
      };
  }

  template<class T>
  M_Matrix<T>  Transpose(const M_Matrix<T>& x)
  {
    switch(x.type())
      {
      case M_Matrix<T>::FULL:
        {
          M_Matrix<T> out(x.ncols(),x.nrows());
          for (std::size_t i=0; i<x.nrows(); i++)
            for (std::size_t j=0; j<x.ncols(); ++j)
              out(j,i)=x(i,j);
          return out;
        }
      case M_Matrix<T>::ZERO:
      case M_Matrix<T>::SYMMETRIC:
      case M_Matrix<T>::DIAGONAL:
      case M_Matrix<T>::SCALAR_DIAGONAL:
      case M_Matrix<T>::SCALAR_FULL:
      default:
        {
          M_Matrix<T> out(x);
          return out;
        }
      }
  }




  template<typename T>
  M_Matrix<T>
  quadraticForm_XXT(const M_Matrix<T>& x)
  {
    M_Matrix<T> out(x.nrows(),x.nrows(),M_Matrix<T>::SYMMETRIC, T{});
    for (std::size_t i=0; i<x.nrows(); ++i)
      for (std::size_t j=0; j<i+1; ++j)
        for (std::size_t k=0; k<x.ncols(); ++k)
          out(i,j)+=x(i,k)*x(j,k);
    return out;
  }

  template<typename T>
  M_Matrix<T>
  quadraticForm_XTX(const M_Matrix<T>& x)
  {
    M_Matrix<T> out(x.ncols(),x.ncols(),M_Matrix<T>::SYMMETRIC, T{});
    for (std::size_t i=0; i<x.ncols(); ++i)
      for (std::size_t j=0; j<i+1; ++j)
        for (std::size_t k=0; k<x.nrows(); ++k)
          out(i,j)+=x(k,i)*x(k,j);
    return out;
  }


  template <class T>
  std::pair<M_Matrix<T>, std::string>
  inv(const M_Matrix<T>& x);

  template <class T>
  std::pair<M_Matrix<T>, std::string>
  chol(const M_Matrix<T>& x, const std::string& kind);



  namespace partition
  {
    template<class E>
    class FullPartition
    {
    public:
      FullPartition(const M_Matrix<E>& x, std::size_t n):
        A_(n,n),
        B_(n,x.ncols()-n),
        C_(x.nrows()-n,n),
        D_{x.nrows()-n,x.ncols()-n}
      {
        for (std::size_t i=0; i<n; ++i)
          {
            for (std::size_t j=0; j<n; ++j)
              A_(i,j)=x(i,j);
            for (std::size_t j=n; j<x.ncols(); ++j)
              B_(i,j-n)=x(i,j);
          }
        for (std::size_t i=n; i<x.nrows(); ++i)
          {
            for (std::size_t j=0; j<n; ++j)
              C_(i-n,j)=x(i,j);
            for (std::size_t j=n; j<x.nrows(); ++j)
              D_(i-n,j-n)=x(i,j);
          }

      }

      FullPartition
      (const M_Matrix<E>& A,
       const M_Matrix<E>& B,
       const M_Matrix<E>& C,
       const M_Matrix<E>& D)
        :A_(A),B_(B),C_(C),D_(D){}

      FullPartition( M_Matrix<E>&& A,
                     M_Matrix<E>&& B,
                     M_Matrix<E>&& C,
                     M_Matrix<E>&& D)
        :A_(A),B_(B),C_(C),D_(D){}

      FullPartition()=default;

      M_Matrix<E> full()const

      {
        std::size_t N=A_.nrows()+D_.nrows();
        M_Matrix<E> out(N,N);
        std::size_t n=A_.nrows();
        for (std::size_t i=0; i<n; ++i)
          {
            for (std::size_t j=0; j<n; ++j)
              out(i,j)=A_(i,j);
            for (std::size_t j=n; j<out.ncols(); ++j)
              out(i,j)=B_(i,j-n);
          }
        for (std::size_t i=n; i<out.nrows(); ++i)
          {
            for (std::size_t j=0; j<n; ++j)
              out(i,j)=C_(i-n,j);
            for (std::size_t j=n; j<out.nrows(); ++j)
              out(i,j)=D_(i-n,j-n);
          }
        return out;

      }

      std::pair<FullPartition, std::string>
      inverse_by_A()const
      {
        auto invA=inv(A_);
        if (!invA.second.empty())
          {
            return {FullPartition(), invA.second};
          }
        auto AinvB=invA.first*B_;
        auto Dinv=inv(D_-C_*AinvB);
        if (!Dinv.second.empty())
          {
            return {FullPartition{}, Dinv.second};
          }

        auto CAinv=C_*invA.first;
        auto Binv=-AinvB*Dinv.first;
        auto Cinv=-Dinv.first*CAinv;
        auto Ainv=invA.first+AinvB*Dinv.first*CAinv;
        return {FullPartition(Ainv,Binv,Cinv,Dinv.first),""};
      }
      std::pair<FullPartition, std::string>
      inverse_by_D()const
      {
        auto invD=inv(D_);
        if (!invD.second.empty())
          {
            return {FullPartition(), invD.second};
          }
        auto DinvC=invD.first*C_;
        auto Ainv=inv(A_-B_*DinvC);
        if (!Ainv.second.empty())
          {
            return {FullPartition(), Ainv.second};
          }


        auto BDinv=B_*invD.first;
        auto Cinv=-DinvC*Ainv.first;
        auto Binv=-Ainv.first*BDinv;
        auto Dinv=invD.first+DinvC*Ainv.first*BDinv;
        return {FullPartition(Ainv.first,Binv,Cinv,Dinv),""};
      }


    private:
      M_Matrix<E> A_;
      M_Matrix<E> B_;
      M_Matrix<E> C_;
      M_Matrix<E> D_;
    };

    template<class E>
    class SymmetricPartition
    {
    public:
      SymmetricPartition(const M_Matrix<E>& x, std::size_t n):
        A_{n,n,M_Matrix<E>::SYMMETRIC},
        B_{n,x.ncols()-n},
        D_{x.nrows()-n,x.ncols()-n,M_Matrix<E>::SYMMETRIC}
      {
        for (std::size_t j=0; j<n; ++j)
          {
            for (std::size_t i=j; i<n; ++i)
              A_(i,j)=x(i,j);
          }
        for (std::size_t j=n; j<x.nrows(); ++j)
          {
            for (std::size_t i=0; i<n; ++i)
              B_(i,j-n)=x(i,j);
            for (std::size_t i=j; i<x.nrows(); ++i)
              D_(i-n,j-n)=x(i,j);
          }

      }

      SymmetricPartition(const M_Matrix<E>& A,const M_Matrix<E>& B,const M_Matrix<E>& D)
        :A_(A),B_(B),D_(D){}
      SymmetricPartition( M_Matrix<E>&& A, M_Matrix<E>&& B,  M_Matrix<E>&& D)
        :A_(A),B_(B),D_(D){}

      SymmetricPartition()=default;

      M_Matrix<E> full()const
      {
        std::size_t N=A_.nrows()+D_.nrows();
        M_Matrix<E> x(N,N, M_Matrix<E>::SYMMETRIC);
        std::size_t n=A_.nrows();
        for (std::size_t j=0; j<n; ++j)
          {
            for (std::size_t i=j; i<n; ++i)
              x(i,j)=A_(i,j);
          }
        for (std::size_t j=n; j<x.nrows(); ++j)
          {
            for (std::size_t i=0; i<n; ++i)
              x(i,j)=B_(i,j-n);
            for (std::size_t i=j; i<x.nrows(); ++i)
              x(i,j)=D_(i-n,j-n);
          }
        return x;

      }

      std::pair<SymmetricPartition, std::string>
      inverse_by_A()const
      {
        auto invA=inv(A_);
        if (!invA.second.empty())
          {
            return{SymmetricPartition(), invA.second};
          }
        auto BTAinvB=quadraticForm_BT_A_B(invA.first,B_);
        auto Dinv=inv(D_-BTAinvB);
        if (!Dinv.second.empty())
          {
            return{SymmetricPartition(), Dinv.second};
          }

        auto AinvB=invA.first*B_;
        auto Binv=-AinvB*Dinv.first;
        auto Ainv=invA.first
            +quadraticForm_BT_A_B(Dinv.first,Transpose(AinvB));
        return {SymmetricPartition(Ainv,Binv,Dinv.first),""};
      }

      std::pair<SymmetricPartition, std::string>
      inverse_by_D()const
      {
        auto invD=inv(D_);
        if (!invD.second.empty())
          {
            return{SymmetricPartition(), invD.second};
          }

        auto BDinvBT=
            quadraticForm_BT_A_B(invD.first,Transpose(B_));
        auto Ainv=inv(A_-BDinvBT);
        if (!Ainv.second.empty())
          {
            return{SymmetricPartition(), Ainv.second};
          }
        auto BDinv=B_*invD.first;
        auto Binv=-Ainv.first*BDinv;
        auto Dinv=invD.first+
            +quadraticForm_BT_A_B(Ainv.first,BDinv);
        return {SymmetricPartition(Ainv.first,Binv,Dinv),""};
      }


      std::pair<FullPartition<E>, std::string>
      chol_by_A(const std::string& kind)const
      {
        if (kind=="upper")
          {
            auto Ua=chol(A_,kind);
            if (!Ua.second.empty()) return {{},"Cholesky block La: "+Ua.second};
            auto B_chol=Forward_Sustitution_Ly_b(Transpose(Ua.first),B_);
            auto Dschur=D_-quadraticForm_XTX(B_chol);
            auto D_chol=chol(Dschur,kind);
            if (!D_chol.second.empty())
              return {{},"Cholesky block Dschur: "+D_chol.second};
            auto Czero=M_Matrix<E>
                (B_chol.ncols(), B_chol.nrows(), M_Matrix<E>::ZERO);
            return {
                FullPartition<E>(Ua.first,B_chol,Czero,D_chol.first),""};
          }
        else
          {
            auto La=chol(A_,kind);
            if (!La.second.empty()) return {{},"Cholesky block La: "+La.second};
            auto C_chol=Transpose(Forward_Sustitution_Ly_b(La.first,B_));
            auto Dschur=D_-quadraticForm_XXT(C_chol);
            auto D_chol=chol(Dschur,kind);
            if (!D_chol.second.empty())
              return {{},"Cholesky block Dschur: "+D_chol.second};
            auto B_zero=M_Matrix<E>
                (C_chol.ncols(), C_chol.nrows(), M_Matrix<E>::ZERO);
            return {
                FullPartition<E>(La.first,B_zero,C_chol,D_chol.first),""};
          }
      }



    private:
      M_Matrix<E> A_;
      M_Matrix<E> B_;
      M_Matrix<E> D_;
    };
  }
  namespace matrix_inverse
  {
    using namespace partition;
    inline
    std::pair<M_Matrix<double>, std::string>
    full_inv(const M_Matrix<double>& a)

    {
      using Matrix_Unary_Transformations::Transpose;
      using lapack::dgetrf_;
      using lapack::dgetri_;
      if (a.size()==0)
        return {a,"EMPTY MATRIX"};
      else
        {
          assert(a.ncols()==a.nrows());
          int INFO=0;
          //  char msg[101];
          int N =a.ncols();
          auto IPIV=std::make_unique<int[]>(N);
          int LWORK;
          int M=N;
          M_Matrix<double> B=Transpose(a);
          int LDA=N;
          //A=new double[n*n];
          double *A= &B[0]; //more efficient code

          dgetrf_(&N, &M, A, &LDA,IPIV.get(),&INFO);

          LWORK= N*N;
          M_Matrix<double> W(N,N);
          double *WORK = &W[0];

          dgetri_(&N,A,&LDA,IPIV.get(),WORK,&LWORK,&INFO);

          if (INFO==0)
            return {B,""};
          else
            return {B,"Singular Matrix on i="+std::to_string(INFO)};
          ;
        }
    }


    template<typename T>
    std::pair<M_Matrix<T>, std::string>
    full_inv(const M_Matrix<T>& x)
    {
      auto p=FullPartition<T>(x,x.nrows()/2);
      auto pinv=p.inverse_by_A();
      if (!pinv.second.empty())
        return {M_Matrix<T>(),pinv.second};
      else
        return {pinv.first.full(),""};

    }

    inline
    std::pair<M_Matrix<double>
    ,std::string>
    symmetric_inv(const M_Matrix<double>& a,const std::string& kind="lower")

    {
      using lapack::dsytrf_;
      using lapack::dsytri_;

      if (a.size()==0)
        return {a,"EMPTY MATRIX"};
      else
        {
          assert(a.nrows()==a.ncols());

          /**
  Purpose:

       DSYTRF computes the factorization of a real symmetric matrix A using
       the Bunch-Kaufman diagonal pivoting method.  The form of the
       factorization is

          A = U*D*U**T  or  A = L*D*L**T

       where U (or L) is a product of permutation and unit upper (lower)
       triangular matrices, and D is symmetric and block diagonal with
       1-by-1 and 2-by-2 diagonal blocks.

       This is the blocked version of the algorithm, calling Level 3 BLAS.

  Parameters
      [in]	UPLO

                UPLO is CHARACTER*1
                = 'U':  Upper triangle of A is stored;
                = 'L':  Lower triangle of A is stored.

      [in]	N

                N is INTEGER
                The order of the matrix A.  N >= 0.

      [in,out]	A

                A is DOUBLE PRECISION array, dimension (LDA,N)
                On entry, the symmetric matrix A.  If UPLO = 'U', the leading
                N-by-N upper triangular part of A contains the upper
                triangular part of the matrix A, and the strictly lower
                triangular part of A is not referenced.  If UPLO = 'L', the
                leading N-by-N lower triangular part of A contains the lower
                triangular part of the matrix A, and the strictly upper
                triangular part of A is not referenced.

                On exit, the block diagonal matrix D and the multipliers used
                to obtain the factor U or L (see below for further details).

      [in]	LDA

                LDA is INTEGER
                The leading dimension of the array A.  LDA >= max(1,N).

      [out]	IPIV

                IPIV is INTEGER array, dimension (N)
                Details of the interchanges and the block structure of D.
                If IPIV(k) > 0, then rows and columns k and IPIV(k) were
                interchanged and D(k,k) is a 1-by-1 diagonal block.
                If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
                columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
                is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
                IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
                interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.

      [out]	WORK

                WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
                On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

      [in]	LWORK

                LWORK is INTEGER
                The length of WORK.  LWORK >=1.  For best performance
                LWORK >= N*NB, where NB is the block size returned by ILAENV.

                If LWORK = -1, then a workspace query is assumed; the routine
                only calculates the optimal size of the WORK array, returns
                this value as the first entry of the WORK array, and no error
                message related to LWORK is issued by XERBLA.

      [out]	INFO

                INFO is INTEGER
                = 0:  successful exit
                < 0:  if INFO = -i, the i-th argument had an illegal value
                > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
                      has been completed, but the block diagonal matrix D is
                      exactly singular, and division by zero will occur if it
                      is used to solve a system of equations.

        */

          char UPLO;
          if (kind=="lower")
            UPLO='U';
          else
            UPLO='L';

          int INFO=0;
          int N =a.ncols();
          auto IPIV=std::make_unique<int[]>(N);
          int LWORK=N*N; //
          M_Matrix<double> B(a.nrows(),a.ncols());
          if (kind!="lower")
            {
              for (std::size_t i=0; i<a.nrows(); ++i)
                for (std::size_t j=i; j<a.ncols(); ++j)
                  B(i,j)=a(i,j);
            } else
            {  for (std::size_t i=0; i<a.nrows(); ++i)
                for (std::size_t j=0; j<i+1; ++j)
                  B(i,j)=a(i,j);
            }

          int LDA=N;
          double *A= &B[0]; //more efficient code
          M_Matrix<double> W(N,N);
          double *WORK = &W[0];
          dsytrf_(&UPLO,&N,A,&LDA,IPIV.get(),WORK,&LWORK,&INFO);

          if (INFO<0)
            {
              return {{},"INVALID ARGUMENT"};
            }
          else if (INFO>0)
            {
              return {{},"SINGULAR MATRIX ON "+std::to_string(INFO)};
            }
          else
            {
              /**
             * dsytri()
  subroutine dsytri 	( 	character  	UPLO,
                  integer  	N,
                  double precision, dimension( lda, * )  	A,
                  integer  	LDA,
                  integer, dimension( * )  	IPIV,
                  double precision, dimension( * )  	WORK,
                  integer  	INFO
          )

  DSYTRI

  Download DSYTRI + dependencies [TGZ] [ZIP] [TXT]

  Purpose:

       DSYTRI computes the inverse of a real symmetric indefinite matrix
       A using the factorization A = U*D*U**T or A = L*D*L**T computed by
       DSYTRF.

  Parameters
      [in]	UPLO

                UPLO is CHARACTER*1
                Specifies whether the details of the factorization are stored
                as an upper or lower triangular matrix.
                = 'U':  Upper triangular, form is A = U*D*U**T;
                = 'L':  Lower triangular, form is A = L*D*L**T.

      [in]	N

                N is INTEGER
                The order of the matrix A.  N >= 0.

      [in,out]	A

                A is DOUBLE PRECISION array, dimension (LDA,N)
                On entry, the block diagonal matrix D and the multipliers
                used to obtain the factor U or L as computed by DSYTRF.

                On exit, if INFO = 0, the (symmetric) inverse of the original
                matrix.  If UPLO = 'U', the upper triangular part of the
                inverse is formed and the part of A below the diagonal is not
                referenced; if UPLO = 'L' the lower triangular part of the
                inverse is formed and the part of A above the diagonal is
                not referenced.

      [in]	LDA

                LDA is INTEGER
                The leading dimension of the array A.  LDA >= max(1,N).

      [in]	IPIV

                IPIV is INTEGER array, dimension (N)
                Details of the interchanges and the block structure of D
                as determined by DSYTRF.

      [out]	WORK

                WORK is DOUBLE PRECISION array, dimension (N)

      [out]	INFO

                INFO is INTEGER
                = 0: successful exit
                < 0: if INFO = -i, the i-th argument had an illegal value
                > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
                     inverse could not be computed.
             */
              dsytri_(&UPLO,&N,A,&LDA,IPIV.get(),WORK,&INFO);
              if (INFO!=0)
                {
                  return {B,"cannot invert a singular matrix"};
                }
              else
                {
                  M_Matrix<double> out
                      (B.nrows(),B.ncols(),M_Matrix<double>::SYMMETRIC);
                  if (kind=="lower")
                    {
                      for (std::size_t i=0; i<B.nrows(); ++i)
                        for(std::size_t j=0; j<=i; ++j)
                          out(i,j)=B(i,j);
                    }
                  else
                    {
                      for (std::size_t i=0; i<B.nrows(); ++i)
                        for(std::size_t j=0; j<=i; ++j)
                          out(i,j)=B(j,i);
                    }

                  /*    auto aa=a;
                 M_Matrix<double> test(a.nrows(),a.ncols(),M_Matrix<double>::FULL,a);

                 auto invtest=inv(test).first;

                 auto test_it=test*out;
                 auto test_a=aa*out;
                 auto invv=invtest*test;
             */
                  return {out,""};
                }
            }
        }
    }


    template<typename T>
    std::pair<M_Matrix<T>, std::string>
    symmetric_inv(const M_Matrix<T>& x)
    {
      auto p=SymmetricPartition<T>(x,x.nrows()/2);
      auto pinv=p.inverse_by_A();
      if (!pinv.second.empty())
        return {M_Matrix<T>(),pinv.second};
      else
        return {pinv.first.full(),""};

    }

    inline
    std::pair<M_Matrix<double>
    ,std::string>
    diagonal_inv(const M_Matrix<double>& a)
    {
      M_Matrix<double> out(a.nrows(),a.ncols(),M_Matrix<double>::DIAGONAL);
      for (std::size_t i=0; i<a.nrows(); ++i)
        if (a(i,i)==0)
          return {{},"Singular Matrix, Diagonal zero at "+std::to_string(i)};
        else
          out(i,i)=1.0/a(i,i);
      return {out,""};

    }


    template<typename T>
    std::pair<M_Matrix<T>, std::string>
    diagonal_inv(const M_Matrix<T>& a)
    {
      M_Matrix<T> out(a.nrows(),a.ncols(),M_Matrix<T>::DIAGONAL);
      for (std::size_t i=0; i<a.nrows(); ++i)
        {
          auto e=inv(a(i,i));
          if (!e.second.empty())
            {
              return {{},"Singular block="+std::to_string(i)+e.second};
            }
          else
            out(i,i)=std::move(e.first);
        }
      return {out,""};

    }

    inline
    std::pair<M_Matrix<double>
    ,std::string>
    scalar_diagonal_inv(const M_Matrix<double>& a)
    {
      if (a(0,0)==0)
        return {{},"Singular Matrix, Diagonal is zero"};
      else
        return {
            M_Matrix<double>
                (a.nrows(),a.ncols(),
                 M_Matrix<double>::SCALAR_DIAGONAL,
                 1.0/a(0,0))
                ,""
          };

    }


    template<typename T>
    std::pair<M_Matrix<T>, std::string>
    scalar_diagonal_inv(const M_Matrix<T>& a)
    {
      auto e=inv(a(0,0));
      if (!e.second.empty())
        {
          return {{},"Singular SCALAR DIAGONAL block: "+e.second};
        }
      else
        return {
            M_Matrix<T>
                (a.nrows(),a.ncols(),
                 M_Matrix<T>::SCALAR_DIAGONAL,
                 std::move(e.first))
                ,""
          };
    }

    template<typename E>
    std::pair<M_Matrix<E>, std::string>
    Matrix_inverse(const M_Matrix<E>& x)
    {
      assert(x.nrows()==x.ncols());
      switch(x.type())
        {
        case M_Matrix<E>::FULL:
          {
            return full_inv(x);
          }
        case M_Matrix<E>::SYMMETRIC:
          {
            return symmetric_inv(x);
          }
        case M_Matrix<E>::DIAGONAL:
          {
            return diagonal_inv(x);
          }
        case M_Matrix<E>::SCALAR_DIAGONAL:
          {
            return scalar_diagonal_inv(x);
          }
        case M_Matrix<E>::SCALAR_FULL:
          {
            return {{},"SCALAR FULL is Singular"};
          }
        case M_Matrix<E>::ZERO:
        default:
          {
            return {{},"ZERO Matrix is Singular"};
          }


        }
    }
  }

  namespace cholesky
  {
    using  Matrix_Unary_Transformations::Transpose;
    using namespace partition;
    inline
    std::pair<M_Matrix<double>, std::string>
    symmetric_chol(const M_Matrix<double>& x,const std::string& kind)
    {
      assert(x.nrows()==x.ncols());

      if (x.size()==0)
        return {{},"cholesky of ZERO MATRIX"};
      char UPLO='L';
      M_Matrix<double> res;
      if (kind=="lower")
        {
          UPLO='U';
          res=lapack::LT(x);
        }
      else
        {
          res=lapack::UT(x);
        }
      int N=x.nrows();
      int LDA=N;
      int INFO;


      lapack::dpotrf_(&UPLO,&N,&res[0],&LDA,&INFO);

      if (INFO!=0)
        {
          return {{},"Cholesky fails, zero diagonal at"+std::to_string(INFO)};
        }
      else
        {
          if (kind=="lower")
            {
              auto test=res*Transpose(res)-x;
            }
          else
            {
              auto test=Transpose(res)*res-x;

            }
          return {res,""};
        }

    }


    template<typename T>
    std::pair<M_Matrix<T>, std::string>
    symmetric_chol(const M_Matrix<T>& x,const std::string& kind)
    {
      auto p=SymmetricPartition<T>(x,x.nrows()/2);
      auto pchol=p.chol_by_A(kind);
      if (!pchol.second.empty())
        return {M_Matrix<T>(),"Block Cholesky fails "+pchol.second};
      else
        return {pchol.first.full(),""};
    }

    inline
    std::pair<M_Matrix<double>
    ,std::string>
    diagonal_chol(const M_Matrix<double>& a, const std::string& /*kind*/)
    {
      M_Matrix<double> out(a.nrows(),a.ncols(),M_Matrix<double>::DIAGONAL);
      for (std::size_t i=0; i<a.nrows(); ++i)
        if (a(i,i)<=0)
          return {{},"Negative Diagonal zero at "+std::to_string(i)+"value "+std::to_string(a(i,i))};
        else
          out(i,i)=std::sqrt(a(i,i));
      return {out,""};

    }


    template<typename T>
    std::pair<M_Matrix<T>, std::string>
    diagonal_chol(const M_Matrix<T>& a, const std::string& kind)
    {
      M_Matrix<T> out(a.nrows(),a.ncols(),M_Matrix<T>::DIAGONAL);
      for (std::size_t i=0; i<a.nrows(); ++i)
        {
          auto e=chol(a(i,i), kind);
          if (!e.second.empty())
            {
              return
              {{},"Non positive definite  block="
                  +std::to_string(i)+"  "+e.second};
            }
          else
            out(i,i)=std::move(e.first);
        }
      return {out,""};

    }


    inline
    std::pair<M_Matrix<double>
    ,std::string>
    scalar_diagonal_chol(const M_Matrix<double>& a, const std::string& /*kind*/)
    {
      if (a(0,0)<=0)
        return {{},"Non Definite positive, Diagonal is non positive"};
      else
        return {
            M_Matrix<double>
                (a.nrows(),a.ncols(),
                 M_Matrix<double>::SCALAR_DIAGONAL,
                 std::sqrt(a(0,0)))
                ,""
          };

    }


    template<typename T>
    std::pair<M_Matrix<T>, std::string>
    scalar_diagonal_chol(const M_Matrix<T>& a, const std::string& kind)
    {
      auto e=chol(a(0,0), kind);
      if (!e.second.empty())
        {
          return {{},"Not positive definite SCALAR DIAGONAL block: "+e.second};
        }
      else
        return {
            M_Matrix<T>
                (a.nrows(),a.ncols(),
                 M_Matrix<T>::SCALAR_DIAGONAL,
                 std::move(e.first))
                ,""
          };
    }



    template<typename E>
    std::pair<M_Matrix<E>, std::string>
    Matrix_cholesky(const M_Matrix<E>& x, const std::string& kind)

    {
      assert(x.nrows()==x.ncols());
      switch(x.type())
        {
        case M_Matrix<E>::FULL:
          {
            return {{},"Error, Cholesky on a FULL MATRIX"};
          }
        case M_Matrix<E>::SYMMETRIC:
          {
            return symmetric_chol(x, kind);
          }
        case M_Matrix<E>::DIAGONAL:
          {
            return diagonal_chol(x, kind);
          }
        case M_Matrix<E>::SCALAR_DIAGONAL:
          {
            return scalar_diagonal_chol(x, kind);
          }
        case M_Matrix<E>::SCALAR_FULL:
          {
            return {{},"SCALAR FULL is Not definite positive"};
          }
        case M_Matrix<E>::ZERO:
        default:
          {
            return {{},"ZERO Matrix is Not definite positive"};
          }


        }
    }




  }

  template <class T>
  std::pair<M_Matrix<T>, std::string>
  inv(const M_Matrix<T>& x)
  {
    return matrix_inverse::Matrix_inverse(x);
  }
  template <class T>
  std::pair<M_Matrix<T>, std::string>
  chol(const M_Matrix<T>& x, const std::string& kind)
  {
    return cholesky::Matrix_cholesky(x, kind);
  }






  template<typename T, class Compare>
  M_Matrix<T> sort(const M_Matrix<T>& x, Compare comp)
  {
    std::vector<T> o=x.toVector();
    std::sort(o.begin(), o.end(), comp);
    return M_Matrix<T>(x.nrows(),x.ncols(),o);

  }


  template<typename T>
  M_Matrix<T> sort(const M_Matrix<T>& x)
  {
    std::vector<T> o=x.toVector();
    std::sort(o.begin(), o.end());
    return M_Matrix<T>(x.nrows(),x.ncols(),o);

  }

  /**
       Diagonal of Matrix or Diagonal Matrix
       It has two behaviors:
       - If the input is a single column or a single row, it builds a diagonal
       Matrix with it
       - If the input is a Matrix, it returns the values of its diagonal

      */
  template<typename T>
  M_Matrix<T> diag(const M_Matrix<T>& x)
  {
    size_t nr=x.nrows();
    size_t nc=x.ncols();
    if ((nr>1)&(nc>1))
      {
        std::size_t n=std::min(nr,nc);
        M_Matrix<T> diagV(1,n);
        for (size_t i=0; i<n; ++i)
          diagV(0,i)=x(i,i);
        return diagV;
      }
    else
      {
        nr=std::max(nr,nc);
        M_Matrix<T> diagM(nr,nr,M_Matrix<T>::DIAGONAL);
        for (size_t i=0; i<nr; ++i)
          diagM(i,i)=x[i];
        return diagM;
      }

  }


  template<typename T>
  M_Matrix<T> col_vector(const M_Matrix<T>& x)
  {
    M_Matrix<T> colvec(x.size(),1);
    for (std::size_t i=0; i<x.size(); ++i)
      colvec[i]=x[i];
    return colvec;
  }
  template<typename T>
  M_Matrix<T> row_vector(const M_Matrix<T>& x)
  {
    M_Matrix<T> rowvec(1,x.size());
    for (std::size_t i=0; i<x.size(); ++i)
      rowvec[i]=x[i];
    return rowvec;
  }

}

/**
  Blas

SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  Purpose
*  =======
*
*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Arguments
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*/


namespace Matrix_Binary_Predicates
{
  using namespace Matrix_Unary_Transformations;

  /*
Matrix Equality Operator
*/
  template<typename T>
  bool operator==(const M_Matrix<T>& x,
                  const M_Matrix<T>& y)
  {
    if (x.type()!=y.type()) return false;
    if (x.size()!=y.size()) return false;
    else if (x.ncols()!=y.ncols()) return false;
    else for (std::size_t i=0; i< x.size(); ++i)
      if (x[i]!=y[i]) return false;
    return true;

  }



  /*
 Minor Operator based on a Lexicographic Comparison.
*/
  template<typename T>
  bool operator<(const M_Matrix<T>& x, const M_Matrix<T>& y)
  {
    if (x.size()<y.size()) return true;
    else if (y.size()<x.size()) return false;
    else if (x.nrows()<y.nrows()) return true;
    else if (y.nrows()<x.nrows()) return false;
    else for (std::size_t i=0; i< x.size(); ++x)
      {
        if (x[i]<y[i]) return true;
        else if (y[i]<x[i]) return false;
      }
    return false;

  }

}


namespace Matrix_Binary_Transformations
{
  using namespace Matrix_Unary_Transformations;

  template<typename T, typename S>
  typename M_Matrix<T>::TYPE
  result_of_sum
  (const  M_Matrix<T>& x,const  M_Matrix<S>&  y)
  {
    if ((x.type()==M_Matrix<T>::FULL)||
        (y.type()==M_Matrix<S>::FULL))
      return M_Matrix<T>::FULL;
    else if ((x.type()==M_Matrix<T>::SYMMETRIC)||
             (y.type()==M_Matrix<S>::SYMMETRIC))
      return M_Matrix<T>::SYMMETRIC;
    else if ((x.type()==M_Matrix<T>::DIAGONAL)||
             (y.type()==M_Matrix<S>::DIAGONAL))
      {
        if ((x.type()==M_Matrix<T>::SCALAR_FULL)||
            (y.type()==M_Matrix<S>::SCALAR_FULL))
          return M_Matrix<T>::SYMMETRIC;
        else
          return M_Matrix<T>::DIAGONAL;
      }
    else if ((x.type()==M_Matrix<T>::SCALAR_DIAGONAL)||
             (y.type()==M_Matrix<S>::SCALAR_DIAGONAL))
      {
        if ((x.type()==M_Matrix<T>::SCALAR_FULL)||
            (y.type()==M_Matrix<S>::SCALAR_FULL))
          return M_Matrix<T>::SYMMETRIC;
        else
          return M_Matrix<T>::SCALAR_DIAGONAL;
      }
    else if ((x.type()==M_Matrix<T>::SCALAR_FULL)||
             (y.type()==M_Matrix<S>::SCALAR_FULL))
      return M_Matrix<T>::SCALAR_FULL;
    else
      return M_Matrix<T>::ZERO;

  }

  /**
   *  use only when x.size()>=y.size()
   */
  template<typename T, typename S>
  typename M_Matrix<S>::ITER_TYPE
  iterator_of_sum
  (const  M_Matrix<T>& x,const  M_Matrix<S>&  y)
  {
    switch(x.type())
      {
      case M_Matrix<T>::FULL:
        {
          switch (y.type())
            {
            case M_Matrix<S>::FULL:
            case M_Matrix<S>::SYMMETRIC:
            case M_Matrix<S>::SCALAR_FULL:
              return M_Matrix<S>::FULL_IT;
            case M_Matrix<S>::DIAGONAL:
            case M_Matrix<S>::SCALAR_DIAGONAL:
              return M_Matrix<S>::DIAG_IT;
            case M_Matrix<S>::ZERO:
            default:
              return M_Matrix<S>::ZERO_IT;
            }
        }
      case M_Matrix<T>::SYMMETRIC:
        {
          switch (y.type())
            {
            case M_Matrix<S>::FULL:
              assert(false);
            case M_Matrix<S>::SYMMETRIC:
            case M_Matrix<S>::SCALAR_FULL:
              return M_Matrix<S>::LT_IT;
            case M_Matrix<S>::DIAGONAL:
            case M_Matrix<S>::SCALAR_DIAGONAL:
              return M_Matrix<S>::DIAG_IT;
            case M_Matrix<S>::ZERO:
            default:
              return M_Matrix<S>::ZERO_IT;
            }
        }
      case M_Matrix<T>::DIAGONAL:
        {
          switch (y.type())
            {
            case M_Matrix<S>::FULL:
            case M_Matrix<S>::SYMMETRIC:
            case M_Matrix<S>::SCALAR_FULL:
              assert(false);
            case M_Matrix<S>::DIAGONAL:
            case M_Matrix<S>::SCALAR_DIAGONAL:
              return M_Matrix<S>::DIAG_IT;
            case M_Matrix<S>::ZERO:
            default:
              return M_Matrix<S>::ZERO_IT;
            }
        }
      case M_Matrix<T>::SCALAR_DIAGONAL:
        {
          switch (y.type())
            {
            case M_Matrix<S>::FULL:
            case M_Matrix<S>::SYMMETRIC:
            case M_Matrix<S>::SCALAR_FULL:
            case M_Matrix<S>::DIAGONAL:
              assert(false);
            case M_Matrix<S>::SCALAR_DIAGONAL:
              return M_Matrix<S>::SCALAR_IT;
            case M_Matrix<S>::ZERO:
            default:
              return M_Matrix<S>::ZERO_IT;
            }
        }
      case M_Matrix<T>::SCALAR_FULL:
        {
          switch (y.type())
            {
            case M_Matrix<S>::FULL:
            case M_Matrix<S>::SYMMETRIC:
            case M_Matrix<S>::DIAGONAL:
            case M_Matrix<S>::SCALAR_DIAGONAL:
              assert(false);
            case M_Matrix<S>::SCALAR_FULL:
              return M_Matrix<S>::SCALAR_IT;
            case M_Matrix<S>::ZERO:
            default:
              return M_Matrix<S>::ZERO_IT;
            }
        }
      case M_Matrix<T>::ZERO:
      default:
        {
          switch (y.type())
            {
            case M_Matrix<S>::FULL:
            case M_Matrix<S>::SYMMETRIC:
            case M_Matrix<S>::DIAGONAL:
            case M_Matrix<S>::SCALAR_DIAGONAL:
            case M_Matrix<S>::SCALAR_FULL:
              assert(false);
            case M_Matrix<S>::ZERO:
            default:
              return M_Matrix<S>::ZERO_IT;
            }
        }

      }

  }





  template<typename T, typename S>
  std::pair<std::size_t, std::size_t>
  Size_of_Product(const M_Matrix<T>& one, const M_Matrix<S>& other, bool transposeOne,bool transposeOther)
  {
    assert(
          (transposeOne ? one.nrows():one.ncols())
          ==
          (transposeOther ? other.ncols():other.nrows())
          );
    return {
        transposeOne? one.ncols():one.nrows(),
            transposeOther? other.nrows():other.ncols()};
      }

    template<typename T, typename S>
    auto
    result_of_Product
    (const M_Matrix<T>& x, const M_Matrix<S>& y,bool transpose_x, bool transpose_y)
    -> M_Matrix<decltype(std::declval<T>()*std::declval<S>())>
    {
    typedef  decltype(std::declval<T>()*std::declval<S>()) R;
    auto s=Size_of_Product(x,y,transpose_x,transpose_y);

    if ((x.type()==M_Matrix<T>::ZERO)||
        (y.type()==M_Matrix<S>::ZERO))
      return M_Matrix<R>(s.first,s.second,M_Matrix<R>::ZERO);
    if ((x.type()==M_Matrix<T>::FULL)||
        (y.type()==M_Matrix<S>::FULL))
      return M_Matrix<R>(s.first,s.second,M_Matrix<R>::FULL);

    else if ((x.type()==M_Matrix<T>::SYMMETRIC)||
             (y.type()==M_Matrix<S>::SYMMETRIC))
      {
        if ((x.type()==M_Matrix<T>::SCALAR_DIAGONAL)||
            (y.type()==M_Matrix<S>::SCALAR_DIAGONAL))
          return M_Matrix<R>(s.first,s.second,M_Matrix<R>::SYMMETRIC);
        else
          return M_Matrix<R>(s.first,s.second,M_Matrix<R>::FULL);

      }
    else if ((x.type()==M_Matrix<T>::DIAGONAL)||
             (y.type()==M_Matrix<S>::DIAGONAL))
      {
        if ((x.type()==M_Matrix<T>::SCALAR_FULL)||
            (y.type()==M_Matrix<S>::SCALAR_FULL))
          return M_Matrix<R>(s.first,s.second,M_Matrix<R>::FULL);
        else
          return M_Matrix<R>(s.first,s.second,M_Matrix<R>::DIAGONAL);
      }
    else if ((x.type()==M_Matrix<T>::SCALAR_DIAGONAL)
             &&(y.type()==M_Matrix<S>::SCALAR_DIAGONAL))

      return M_Matrix<R>(s.first,s.second,M_Matrix<R>::SCALAR_DIAGONAL);
    else
      return M_Matrix<R>(s.first,s.second,M_Matrix<R>::SCALAR_FULL);
  }





  namespace product
  {
    inline M_Matrix<double> Lapack_Full_Product(const M_Matrix<double>& x,const M_Matrix<double>& y, bool transpose_x, bool transpose_y)
    {
      using lapack::dgemm_;

      std::size_t cols_i, rows_i, rows_e,cols_e;
      if (transpose_x)
        {
          rows_e=x.ncols();
          cols_i=x.nrows();
        }
      else
        {
          rows_e=x.nrows();
          cols_i=x.ncols();
        }

      if (transpose_y)
        {
          rows_i=y.ncols();
          cols_e=y.nrows();
        }
      else
        {
          rows_i=y.nrows();
          cols_e=y.ncols();
        }


      assert(rows_i==cols_i);
      // First it has to find out if the last dimension of x matches the
      //first of y
      // now we build the M_Matrix result
      M_Matrix<double> z(rows_e,cols_e,0.0);

      /***  as fortran uses the reverse order for matrices and we want to
            avoid a copying operation, we calculate
                Transpose(Z)=Transpose(y)*Transpose(x)

                Transpose(matrix)=just plain matrix in C++ format


            */
      char TRANSA;
      char TRANSB;

      if (transpose_y)
        TRANSA='T';
      else
        TRANSA='N';

      if (transpose_x)
        TRANSB='T';
      else
        TRANSB='N';

      int  	M=cols_e;
      int  	N=rows_e;
      int  	K=cols_i;

      double  ALPHA=1.0;
      double*  A=const_cast<double*> (&y[0]);
      int  	LDA;
      if (transpose_y)
        LDA=K;
      else LDA=M;


      double*  B=const_cast<double*> (&x[0]);

      int  	LDB;
      if (transpose_x)
        LDB=N;
      else
        LDB=K;

      double BETA=0.0;

      double * C=&z[0];

      int  	LDC=M;


      try
      {
        dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);
      }
      catch (...)
      {
        assert(false);
      }
      return z;
    }


    inline M_Matrix<double> Lapack_Symmetric_Regular_Product
    (const M_Matrix<double>& Sym,
     const M_Matrix<double>& Reg,
     bool SymReg=true,
     const std::string& kind="upper")
    {
      using lapack::dsymm_;

      /**
     *
     * SSYMM

  Purpose:

       SSYMM  performs one of the matrix-matrix operations

          C := alpha*A*B + beta*C,

       or

          C := alpha*B*A + beta*C,

       where alpha and beta are scalars,  A is a symmetric matrix and  B and
       C are  m by n matrices.

  Parameters
      [in]	SIDE

                SIDE is CHARACTER*1
                 On entry,  SIDE  specifies whether  the  symmetric matrix  A
                 appears on the  left or right  in the  operation as follows:

                    SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,

                    SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,

      [in]	UPLO

                UPLO is CHARACTER*1
                 On  entry,   UPLO  specifies  whether  the  upper  or  lower
                 triangular  part  of  the  symmetric  matrix   A  is  to  be
                 referenced as follows:

                    UPLO = 'U' or 'u'   Only the upper triangular part of the
                                        symmetric matrix is to be referenced.

                    UPLO = 'L' or 'l'   Only the lower triangular part of the
                                        symmetric matrix is to be referenced.


      [in]	LDA

                LDA is INTEGER
                 On entry, LDA specifies the first dimension of A as declared
                 in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
                 LDA must be at least  max( 1, m ), otherwise  LDA must be at
                 least  max( 1, n ).

      [in]	B

                B is REAL array, dimension ( LDB, N )
                 Before entry, the leading  m by n part of the array  B  must
                 contain the matrix B.

      [in]	LDB

                LDB is INTEGER
                 On entry, LDB specifies the first dimension of B as declared
                 in  the  calling  (sub)  program.   LDB  must  be  at  least
                 max( 1, m ).

      [in]	BETA

                BETA is REAL
                 On entry,  BETA  specifies the scalar  beta.  When  BETA  is
                 supplied as zero then C need not be set on input.

      [in,out]	C

                C is REAL array, dimension ( LDC, N )
                 Before entry, the leading  m by n  part of the array  C must
                 contain the matrix  C,  except when  beta  is zero, in which
                 case C need not be set on entry.
                 On exit, the array  C  is overwritten by the  m by n updated
                 matrix.

      [in]	LDC

                LDC is INTEGER
                 On entry, LDC specifies the first dimension of C as declared
                 in  the  calling  (sub)  program.   LDC  must  be  at  least
                 max( 1, m ).


     *
     * */
      /**
       * @brief UPLO
        [in]	UPLO

                  UPLO is CHARACTER*1
                   On  entry,   UPLO  specifies  whether  the  upper  or  lower
                   triangular  part  of  the  symmetric  matrix   A  is  to  be
                   referenced as follows:

                      UPLO = 'U' or 'u'   Only the upper triangular part of the
                                          symmetric matrix is to be referenced.

                      UPLO = 'L' or 'l'   Only the lower triangular part of the
                                          symmetric matrix is to be referenced.
       */
      char UPLO;
      if (kind=="lower")
        UPLO='U';
      else
        UPLO='L';

      M_Matrix<double>  S(Sym.nrows(),Sym.ncols());
      if (kind=="lower")
        {
          for (std::size_t i=0; i<Sym.nrows(); i++)
            for (std::size_t j=0; j<=i; ++j)
              S(i,j)=Sym(i,j);
        }
      else
        {
          for (std::size_t i=0; i<Sym.nrows(); i++)
            for (std::size_t j=i; j<Sym.ncols(); ++j)
              S(i,j)=Sym(i,j);
        }


      std::size_t rows_e;
      std::size_t cols_i;
      std::size_t rows_i;
      std::size_t cols_e;

      if (SymReg)  // calculo Sym * Reg
        {
          rows_e=Sym.nrows();
          cols_i=Sym.ncols();
          rows_i=Reg.nrows();
          cols_e=Reg.ncols();

        }
      else   // calculo  Reg^T * Sym
        {
          rows_e=Reg.nrows();
          cols_i=Reg.ncols();
          rows_i=Sym.nrows();
          cols_e=Sym.ncols();
        }
      assert(rows_i==cols_i);
      // First it has to find out if the last dimension of x matches the
      //first of y
      // now we build the M_Matrix result
      M_Matrix<double> z(rows_e,cols_e,0.0);

      /***  as fortran uses the reverse order for matrices and we want to
            avoid a copying operation, we calculate
                Transpose(Z)=Transpose(y)*Transpose(x)

                Transpose(matrix)=just plain matrix in C++ format



           */
      /**
     * @brief SIDE
     *   [in]	SIDE

                SIDE is CHARACTER*1
                 On entry,  SIDE  specifies whether  the  symmetric matrix  A
                 appears on the  left or right  in the  operation as follows:

                    SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,

                    SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,


     */
      char SIDE;
      if (SymReg)
        SIDE='R';
      else
        SIDE='L';

      /**
     * @brief M
      [in]	M

                M is INTEGER
                 On entry,  M  specifies the number of rows of the matrix  C.
                 M  must be at least zero.
     */

      int  	M=cols_e;

      /**
     * @brief N
      [in]	N

                N is INTEGER
                 On entry, N specifies the number of columns of the matrix C.
                 N  must be at least zero.

                 como uso la transpuesta de C, es rows_e

     */

      int  	N=rows_e;

      /**
     * @brief ALPHA
      [in]	ALPHA

                ALPHA is REAL
                 On entry, ALPHA specifies the scalar alpha.
     */
      double  	ALPHA=1.0;

      /**
     * @brief A
      [in]	A

                A is REAL array, dimension ( LDA, ka ), where ka is
                 m  when  SIDE = 'L' or 'l'  and is  n otherwise.
                 Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
                 the array  A  must contain the  symmetric matrix,  such that
                 when  UPLO = 'U' or 'u', the leading m by m upper triangular
                 part of the array  A  must contain the upper triangular part
                 of the  symmetric matrix and the  strictly  lower triangular
                 part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
                 the leading  m by m  lower triangular part  of the  array  A
                 must  contain  the  lower triangular part  of the  symmetric
                 matrix and the  strictly upper triangular part of  A  is not
                 referenced.
                 Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
                 the array  A  must contain the  symmetric matrix,  such that
                 when  UPLO = 'U' or 'u', the leading n by n upper triangular
                 part of the array  A  must contain the upper triangular part
                 of the  symmetric matrix and the  strictly  lower triangular
                 part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
                 the leading  n by n  lower triangular part  of the  array  A
                 must  contain  the  lower triangular part  of the  symmetric
                 matrix and the  strictly upper triangular part of  A  is not
                 referenced.

     */

      double * /*, dimension(lda,*)*/ A=const_cast<double*> (&S[0]);
      int  	LDA=S.ncols();
      double*  B=const_cast<double*> (&Reg[0]);

      int  	LDB=Reg.ncols();
      double   	BETA=0;
      double * /*, dimension(ldc,*) */ 	C;
      C=&z[0];
      int  	LDC=M;


      dsymm_(&SIDE,&UPLO,&M,&N,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);
      return z;
    }


    inline M_Matrix<double> Lapack_Transpose_Symmetric_Product(const M_Matrix<double>& Reg,const M_Matrix<double>& Sym, const std::string kind="upper")
    {
      return Lapack_Symmetric_Regular_Product(Sym,Reg,false, kind);
    }

    template<typename T, typename S>
    auto All_Product(const M_Matrix<T>& x, const M_Matrix<S>& y,
                     bool transpose_x, bool transpose_y)
    ->
    M_Matrix<decltype(std::declval<T>()*std::declval<S>())>
    {

      auto out=result_of_Product(x,y,transpose_x,transpose_y);
      if (x.isDiagonal())
        {
          if (! transpose_y)
            {
              for (auto it=out.begin();  it!=out.end(); ++it)
                {
                  auto i=it.iRow(); auto j=it.jCol();
                  *it=x(i,i)*y(i,j);
                }
            }
          else
            {
              for (auto it=out.begin();  it!=out.end(); ++it)
                {
                  auto i=it.iRow(); auto j=it.jCol();
                  *it=x(i,i)*y(j,i);
                }

            }
        }
      else if (y.isDiagonal())
        {
          if (! transpose_x)
            {
              for (auto it=out.begin();  it!=out.end(); ++it)
                {
                  auto i=it.iRow(); auto j=it.jCol();
                  *it=x(i,j)*y(j,j);
                }
            }
          else
            {
              for (auto it=out.begin();  it!=out.end(); ++it)
                {
                  auto i=it.iRow(); auto j=it.jCol();
                  *it=x(j,i)*y(j,j);
                }

            }
        }
      else  if(!transpose_x && ! transpose_y)
        {
          for (auto it=out.begin();  it!=out.end(); ++it)
            {
              auto i=it.iRow(); auto j=it.jCol();
              *it=x(i,0)*y(0,j);
              for (std::size_t k=1; k<x.ncols(); ++k)
                *it+=x(i,k)*y(k,j);
            }
        }
      else  if(transpose_x && ! transpose_y)
        {
          for (auto it=out.begin();  it!=out.end(); ++it)
            {
              auto i=it.iRow(); auto j=it.jCol();
              *it=x(0,i)*y(0,j);
              for (std::size_t k=1; k<x.nrows(); ++k)
                *it+=x(k,i)*y(k,j);
            }
        }
      else  if(!transpose_x &&  transpose_y)
        {
          for (auto it=out.begin();  it!=out.end(); ++it)
            {
              auto i=it.iRow(); auto j=it.jCol();
              *it=x(i,0)*y(j,0);
              for (std::size_t k=1; k<x.ncols(); ++k)
                *it+=x(i,k)*y(j,k);
            }
        }
      else  if(transpose_x &&  transpose_y)
        {
          for (auto it=out.begin();  it!=out.end(); ++it)
            {
              auto i=it.iRow(); auto j=it.jCol();
              *it=x(0,i)*y(j,0);
              for (std::size_t k=1; k<x.nrows(); ++k)
                *it+=x(k,i)*y(j,k);
            }
        }
      return out;
    }






    template<typename T, typename S>
    auto
    quadraticForm_Symmetric (const M_Matrix<T>& A, const M_Matrix<S>& B)
    ->M_Matrix<decltype(std::declval<T>()*std::declval<S>())>
    {
      typedef   decltype(std::declval<T>()*std::declval<S>()) R;

      auto AB=A*B;
      M_Matrix<R> out(B.ncols(),B.ncols(),M_Matrix<R>::SYMMETRIC, R{});
      for (std::size_t i=0; i<out.nrows(); ++i)
        for (std::size_t j=i; j<out.ncols(); ++j)
          for (std::size_t k=0; k<B.nrows(); ++k)
            out(i,j)+=B(k,i)*AB(k,j);
      return out;

    }

    template<typename T, typename S>
    auto
    quadraticForm_Symmetric_T (const M_Matrix<T>& A, const M_Matrix<S>& B)
    ->M_Matrix<decltype(std::declval<T>()*std::declval<S>())>
    {
      typedef   decltype(std::declval<T>()*std::declval<S>()) R;

      auto BA=B*A;
      M_Matrix<R> out(B.nrows(),B.nrows(),M_Matrix<R>::SYMMETRIC, R{});
      for (std::size_t i=0; i<out.nrows(); ++i)
        for (std::size_t j=i; j<out.ncols(); ++j)
          for (std::size_t k=0; k<B.ncols(); ++k)
            out(i,j)+=BA(i,k)*B(j,k);
      return out;

    }


    /*

    inline
    M_Matrix<double>
    Matrix_Mult_double
    (const M_Matrix<double>& one,
     const M_Matrix<double>& other,
     bool transposeOne,
     bool transposeOther)
    {
      std::size_t nrowsFirst, ncolsFirst, nrowsSecond,ncolsSecond;

      if (transposeOne)
        {
          nrowsFirst=one.ncols();
          ncolsFirst=one.nrows();
        }
      else
        {
          nrowsFirst=one.nrows();
          ncolsFirst=one.ncols();
        }
      if (transposeOther)
        {
          nrowsSecond=other.ncols();
          ncolsSecond=other.nrows();
        }
      else
        {
          nrowsSecond=other.nrows();
          ncolsSecond=other.ncols();
        }

      assert(ncolsFirst==nrowsSecond);
      if (one.type()==M_Matrix<double>::FULL)
        {
          if (other.type()==M_Matrix<double>::FULL)
            return Lapack_Full_Product
                (one,other,transposeOne,transposeOther);
          else if (other.type()==M_Matrix<double>::SYMMETRIC)
            {
              if (transposeOne)
                return Lapack_Transpose_Symmetric_Product(one,other);
              else {

                  auto tr=Transpose(one);
                  return Lapack_Symmetric_Transpose_Product(tr,other);
                }
            }
          else if ((other.type()==M_Matrix<double>::DIAGONAL)
                   || (other.type()==M_Matrix<double>::SCALAR_DIAGONAL))
            {
              M_Matrix<double> out(nrowsFirst,ncolsSecond);
              if (transposeOne)
                for (std::size_t i=0; i<out.ncols(); ++i)
                  {
                    for (std::size_t j=0; j<one.nrows(); ++j)
                      out(i,j)=one(j,i)*other(j,j);
                    for (std::size_t j=one.nrows(); j<out.ncols();++j)
                      out(i,j)=0;
                  }
              else
                for (std::size_t i=0; i<out.nrows(); ++i)
                  {
                    for (std::size_t j=0; j<one.ncols(); ++j)
                      out(i,j)=one(i,j)*other(j,j);
                    for (std::size_t j=one.ncols(); j<out.ncols();++j)
                      out(i,j)=0;
                  }
              return out;
            }
          else if (other.type()==M_Matrix<double>::SCALAR_FULL)
            {
              M_Matrix<double> out(nrowsFirst,ncolsSecond);
              if (transposeOne)
                for (std::size_t i=0; i<out.nrows(); ++i)
                  {
                    double sum=0;
                    for (std::size_t j=0; j<one.nrows(); ++j)
                      sum+=one(j,i);
                    sum*=other[0];
                    for (std::size_t k=0; k<other.ncols(); ++k)
                      out(i,k)=sum;
                  }
              else for (std::size_t i=0; i<out.nrows(); ++i)
                {
                  double sum=0;
                  for (std::size_t j=0; j<one.ncols(); ++j)
                    sum+=one(i,j);
                  sum*=other[0];
                  for (std::size_t k=0; k<other.ncols(); ++k)
                    out(i,k)=sum;
                }
              return out;
            }
          else // ZERO!
            {
              return M_Matrix<double>
                  (one.nrows(),other.ncols(),M_Matrix<double>::ZERO);

            }
        }
      else   if (one.type()==M_Matrix<double>::SYMMETRIC)
        {
          if (other.type()==M_Matrix<double>::FULL)
            {
              if (transposeOther)
                {
                  return Lapack_Symmetric_Transpose_Product(one,other,true);
                }
              else
                {
                  auto tr=Transpose(other);
                  return Lapack_Symmetric_Transpose_Product(one,tr,true);
                }
            }
          else if (other.type()==M_Matrix<double>::SYMMETRIC)
            {
              auto full=other.full();
              return Lapack_Symmetric_Transpose_Product(one,full,true);
            }
          else if ((other.type()==M_Matrix<double>::DIAGONAL)
                   || (other.type()==M_Matrix<double>::SCALAR_DIAGONAL))
            {
              M_Matrix<double> out
                  (nrowsFirst,ncolsSecond);
              for (std::size_t i=0; i<nrowsSecond; ++i)
                {
                  for (std::size_t j=0; j<ncolsSecond; ++j)
                    out(i,j)=one(i,j)*other(j,j);
                }
              return out;
            }
          else if (other.type()==M_Matrix<double>::SCALAR_FULL)
            {
              M_Matrix<double> out(one.nrows(),other.ncols());
              for (std::size_t i=0; i<out.nrows(); ++i)
                {
                  double sum=one(i,i);
                  for (std::size_t j=0; j<one.ncols(); ++j)
                    sum+=one(i,j);
                  sum*=other[0];
                  for (std::size_t k=0; k<other.ncols(); ++k)
                    out(i,k)=sum;
                }
              return out;
            }
          else // ZERO!
            {
              return M_Matrix<double>
                  (one.nrows(),other.ncols(),M_Matrix<double>::ZERO);

            }

        }
      else  if (one.type()==M_Matrix<double>::DIAGONAL)
        {
          if (other.type()==M_Matrix<double>::FULL)
            {
              if (transposeOther)
                {
                  M_Matrix<double> out(one.nrows(),other.nrows());
                  for (std::size_t i=0; i<other.ncols(); ++i)
                    {
                      for (std::size_t j=0; j<other.nrows(); ++j)
                        out(i,j)=one(i,i)*other(j,i);
                    }
                  for (std::size_t i=one.ncols(); i<one.ncols();++i)
                    for (std::size_t j=0; j<other.nrows(); ++j)
                      out(i,j)=0;

                  return out;

                }
              else
                {
                  M_Matrix<double> out(one.nrows(),other.ncols());
                  for (std::size_t i=0; i<other.nrows(); ++i)
                    {
                      for (std::size_t j=0; j<other.ncols(); ++j)
                        out(i,j)=one(i,i)*other(i,j);
                    }
                  for (std::size_t i=one.ncols(); i<one.ncols();++i)
                    for (std::size_t j=0; j<other.ncols(); ++j)
                      out(i,j)=0;

                  return out;
                }
            }
          else if (other.type()==M_Matrix<double>::SYMMETRIC)
            {
              M_Matrix<double> out
                  (one.nrows(),other.ncols());
              for (std::size_t i=0; i<other.nrows(); ++i)
                {
                  for (std::size_t j=0; j<=i; ++j)
                    out(i,j)=one(i,i)*other(i,j);
                }
              return out;
            }
          else if ((other.type()==M_Matrix<double>::DIAGONAL)
                   || (other.type()==M_Matrix<double>::SCALAR_DIAGONAL))
            {
              M_Matrix<double> out
                  (M_Matrix<double>::DIAGONAL,one.nrows(),other.ncols());
              for (std::size_t i=0; i<out.nrows(); ++i)
                {
                  out(i,i)=one(i,i)*other(i,i);
                }
              return out;
            }
          else if (other.type()==M_Matrix<double>::SCALAR_FULL)
            { M_Matrix<double> out(one.nrows(),other.ncols());
              for (std::size_t i=0; i<out.nrows(); ++i)
                {
                  double d=one(i,i)*other(i,i);
                  for (std::size_t j=0; j<other.ncols(); ++j)
                    out(i,j)=d;
                }
              return out;
            }
          else // ZERO!
            {
              return M_Matrix<double>
                  (one.nrows(),other.ncols(),M_Matrix<double>::ZERO);

            }

        }

      else  if (one.type()==M_Matrix<double>::SCALAR_DIAGONAL)
        {
          if (other.type()==M_Matrix<double>::FULL)
            {
              if (transposeOther)
                {
                  M_Matrix<double> out(one.nrows(),other.nrows());
                  for (std::size_t i=0; i<other.ncols(); ++i)
                    {
                      for (std::size_t j=0; j<other.nrows(); ++j)
                        out(i,j)=one(i,i)*other(j,i);
                    }
                  for (std::size_t i=one.ncols(); i<one.ncols();++i)
                    for (std::size_t j=0; j<other.nrows(); ++j)
                      out(i,j)=0;

                  return out;

                }
              else
                {
                  M_Matrix<double> out(one.nrows(),other.ncols());
                  for (std::size_t i=0; i<other.nrows(); ++i)
                    {
                      for (std::size_t j=0; j<other.ncols(); ++j)
                        out(i,j)=one(i,i)*other(i,j);
                    }
                  for (std::size_t i=one.ncols(); i<one.ncols();++i)
                    for (std::size_t j=0; j<other.ncols(); ++j)
                      out(i,j)=0;

                  return out;
                }
            }
          else if (other.type()==M_Matrix<double>::SYMMETRIC)
            {
              M_Matrix<double> out
                  (one.nrows(),other.ncols(),M_Matrix<double>::SYMMETRIC);
              for (std::size_t i=0; i<other.nrows(); ++i)
                {
                  for (std::size_t j=0; j<=i; ++j)
                    out(i,j)=one(i,i)*other(i,j);
                }
              return out;
            }
          else if (other.type()==M_Matrix<double>::DIAGONAL)

            {
              M_Matrix<double> out
                  (one.nrows(),other.ncols(),M_Matrix<double>::DIAGONAL);
              for (std::size_t i=0; i<out.nrows(); ++i)
                {
                  out(i,i)=one(i,i)*other(i,i);
                }
              return out;
            }
          else if ((other.type()==M_Matrix<double>::SCALAR_DIAGONAL))
            {
              return M_Matrix<double>
                  (one.nrows(),other.ncols(),M_Matrix<double>::SCALAR_DIAGONAL,
                   one[0]*other[0]);
            }
          else if (other.type()==M_Matrix<double>::SCALAR_FULL)
            {
              return  M_Matrix<double>
                  (one.nrows(),other.ncols(),
                   M_Matrix<double>::SCALAR_FULL,one[0]*other[0]);
            }
          else // ZERO!
            {
              return M_Matrix<double>
                  (one.nrows(),other.ncols(),M_Matrix<double>::ZERO);

            }

        }

      else  if (one.type()==M_Matrix<double>::SCALAR_FULL)
        {
          if (other.type()==M_Matrix<double>::FULL)
            {
              if (transposeOther)
                {
                  M_Matrix<double> out(one.nrows(),other.nrows());
                  for (std::size_t k=0; k<other.nrows(); ++k)
                    {
                      double sum=0;
                      for (std::size_t j=0; j<one.ncols(); ++j)
                        sum+=other(k,j);
                      sum*=one[0];
                      for (std::size_t i=0; i<out.nrows(); ++i)
                        out(i,k)=sum;
                    }
                  return out;

                }
              else
                {
                  M_Matrix<double> out(one.nrows(),other.ncols());
                  for (std::size_t k=0; k<other.ncols(); ++k)
                    {
                      double sum=0;
                      for (std::size_t j=0; j<one.ncols(); ++j)
                        sum+=other(j,k);
                      sum*=one[0];
                      for (std::size_t i=0; i<out.nrows(); ++i)
                        out(i,k)=sum;
                    }
                  return out;
                }
            }
          else if (other.type()==M_Matrix<double>::SYMMETRIC)
            {
              M_Matrix<double> out(one.nrows(),other.ncols());
              for (std::size_t i=0; i<other.nrows(); ++i)
                {
                  double sum=other(i,i);
                  for (std::size_t j=0; j<i; ++j)
                    sum+=other(i,j);
                  sum*=one[0];
                  for (std::size_t k=0; k<one.nrows(); ++k)
                    out(k,i)=sum;
                }
              return out;
            }
          else if (other.type()==M_Matrix<double>::DIAGONAL)
            {
              M_Matrix<double> out(one.nrows(),other.ncols());
              for (std::size_t i=0; i<out.nrows(); ++i)
                {
                  double d=one(i,i)*other(i,i);
                  for (std::size_t j=0; j<other.ncols(); ++j)
                    out(i,j)=d;
                }
              return out;
            }
          else if ((other.type()==M_Matrix<double>::SCALAR_DIAGONAL))
            {
              return M_Matrix<double>
                  (one.nrows(),other.ncols(),M_Matrix<double>::SCALAR_FULL,
                   one[0]*other[0]);
            }
          else if (other.type()==M_Matrix<double>::SCALAR_FULL)
            {
              return  M_Matrix<double>
                  (one.nrows(),other.ncols(),
                   M_Matrix<double>::SCALAR_FULL,one[0]*other[0]*one.ncols());
            }
          else // ZERO!
            {
              return M_Matrix<double>
                  (one.nrows(),other.ncols(),M_Matrix<double>::ZERO);

            }

        }

      else //ZERO
        {
          return M_Matrix<double>
              (one.nrows(),other.ncols(),M_Matrix<double>::ZERO);
        }
    }

*/

    inline
    M_Matrix<double>
    Matrix_Multiplication
    (const M_Matrix<double>&one,
     const M_Matrix<double>& other,
     bool transpose_x,
     bool transpose_y)

    {
      if ((one.type()==M_Matrix<double>::FULL)&&
          (other.type()==M_Matrix<double>::FULL))
        return Lapack_Full_Product(one,other,transpose_x,transpose_y);

      else if ((one.type()==M_Matrix<double>::FULL)&&
               (other.type()==M_Matrix<double>::SYMMETRIC))
        {
          if (transpose_x)
            {
              auto tr=Transpose(one);
              return Lapack_Symmetric_Regular_Product(other,tr,false);
            }
          else
            {
              return Lapack_Symmetric_Regular_Product(other,one,false);
            }
        }
      else if ((one.type()==M_Matrix<double>::SYMMETRIC)&&
               (other.type()==M_Matrix<double>::FULL))
        {
          if (transpose_y)
            {
              auto tr=Transpose(other);

              return Lapack_Symmetric_Regular_Product(one,tr,true);
            }
          else
            {
              return Lapack_Symmetric_Regular_Product(one,other,true);
            }
        }
      else if ((one.type()==M_Matrix<double>::SYMMETRIC)&&
               (other.type()==M_Matrix<double>::SYMMETRIC))
        {
          auto copy_other=M_Matrix<double>(other.nrows(),other.ncols(), M_Matrix<double>::FULL,other);
          return Lapack_Symmetric_Regular_Product(one,copy_other,true);
        }
      else
        return All_Product(one, other, transpose_x,transpose_y);
    }


    template<typename T, typename S>
    auto
    Matrix_Multiplication
    (const M_Matrix<T>& one,
     const M_Matrix<S>& other,
     bool transposeOne
     , bool transposeOther)
    ->M_Matrix<decltype(std::declval<T>()*std::declval<S>())>
    {
      return All_Product(one,other,transposeOne,transposeOther);
    }





  }

  template<typename T, typename S>
  auto
  operator *
  (const M_Matrix<T>& one,
   const M_Matrix<S>& other)
  ->M_Matrix<decltype(std::declval<T>()*std::declval<S>())>
  {
    return product::Matrix_Multiplication(one,other,false,false);
  }

  template<typename T, typename S>
  auto
  multTransp
  (const M_Matrix<T>& one,
   const M_Matrix<S>& other)
  ->M_Matrix<decltype(std::declval<T>()*std::declval<S>())>
  {
    return product::Matrix_Multiplication(one,other,false,true);
  }

  /**
     Transpose the first and multiply by the second
     @post transpMult(x,y)==Transpose(x)*y
     @remarks It is faster, since we save copying matrices
    */
  template<typename T, typename S>
  auto
  TranspMult
  (const M_Matrix<T>& one,
   const M_Matrix<S>& other)
  ->M_Matrix<decltype(std::declval<T>()*std::declval<S>())>
  {
    return product::Matrix_Multiplication(one,other,true,false);
  }


  inline
  M_Matrix<double>
  Forward_Sustitution_Ly_b(const M_Matrix<double>& L,
                           const M_Matrix<double>& b)
  {
    assert(L.nrows()==b.nrows());
    M_Matrix<double> y(L.ncols(),b.ncols());
    for (std::size_t k=0; k<b.ncols(); ++k)
      {
        y(0,k)=b(0,k)/L(0,0);
      }
    for (std::size_t i=0; i<L.nrows(); i++)
      {
        for (std::size_t k=0; k<b.ncols(); ++k)
          {
            auto bb=L(i,0)*y(0,k);
            for (std::size_t j=1; j<i; ++j)
              bb+=L(i,j)*y(j,k);
            bb=b(i,k)-bb;
            y(i,k)=bb/L(i,i);
          }
      }
    return y;
  }



  template<typename T>
  M_Matrix<T>
  Forward_Sustitution_Ly_b(const M_Matrix<T>& L, const M_Matrix<T>& b)
  {
    assert(L.nrows()==b.nrows());
    M_Matrix<T> y(L.ncols(),b.ncols());
    for (std::size_t k=0; k<b.ncols(); ++k)
      {
        y(0,k)=Forward_Sustitution_Ly_b(L(0,0),b(0,k));
      }
    for (std::size_t i=0; i<L.nrows(); i++)
      {
        for (std::size_t k=0; k<b.ncols(); ++k)
          {
            auto bb=L(i,0)*y(0,k);
            for (std::size_t j=1; j<i; ++j)
              bb+=L(i,j)*y(j,k);
            bb=b(i,k)-bb;
            y(i,k)=Forward_Sustitution_Ly_b(L(i,i),bb);
          }
      }
    return y;
  }


  inline
  M_Matrix<double>
  Back_Sustitution_Ux_b(const M_Matrix<double>& U, const M_Matrix<double>& b)
  {
    assert(U.nrows()==b.nrows());
    M_Matrix<double> x(U.ncols(),b.ncols());
    std::size_t n=U.nrows();
    for (std::size_t k=0; k<b.ncols(); ++k)
      {
        x(n-1,k)=b(n-1,k)/U(n-1,n-1);
      }
    for (std::size_t i=1; i<U.nrows(); i++)
      {
        for (std::size_t k=0; k<b.ncols(); ++k)
          {
            auto bb=U(n-1-i,n-1)*x(n-1,k);
            for (std::size_t j=1; j<i; ++j)
              bb+=U(n-1-i,n-1-j)*x(n-1-j,k);
            bb=b(n-1-i,k)-bb;
            x(n-1-i,k)=bb/U(n-1-i,n-1-i);
          }
      }
    return x;
  }



  template<typename T>
  M_Matrix<T>
  Back_Sustitution_Ux_b(const M_Matrix<T>& U, const M_Matrix<T>& b)
  {
    assert(U.nrows()==b.nrows());
    M_Matrix<T> x(U.ncols(),b.ncols());
    std::size_t n=U.nrows();
    for (std::size_t k=0; k<b.ncols(); ++k)
      {
        x(n-1,k)=Back_Sustitution_Ux_b (U(n-1,n-1),b(n-1,k));
      }
    for (std::size_t i=1; i<U.nrows(); i++)
      {
        for (std::size_t k=0; k<b.ncols(); ++k)
          {
            auto bb=U(n-1-i,n-1)*x(n-1,k);
            for (std::size_t j=1; j<i; ++j)
              bb+=U(n-1-i,n-1-j)*x(n-1-j,k);
            bb=b(n-1-i,k)-bb;
            x(n-1-i,k)=Back_Sustitution_Ux_b(U(n-1-i,n-1-i),bb);
          }
      }
    return x;
  }





  template<typename T, typename S>
  auto
  quadraticForm_BT_A_B (const M_Matrix<T>& A, const M_Matrix<S>& B)
  ->M_Matrix<decltype(std::declval<T>()*std::declval<S>())>
  {
    assert(A.ncols()==B.nrows());
    assert(B.nrows()==A.nrows());
    typedef   decltype(std::declval<T>()*std::declval<S>()) R;
    if ((A.type()==M_Matrix<T>::ZERO)||
        (B.type()==M_Matrix<S>::ZERO))
      return M_Matrix<R>(B.ncols(),B.ncols(),M_Matrix<R>::ZERO);
    else if (A.isSymmetric())
      return product::quadraticForm_Symmetric(A,B);
    else
      return TranspMult(B,A*B);
  }


  template<typename T, typename S>
  auto
  quadraticForm_B_A_BT (const M_Matrix<T>& A, const M_Matrix<S>& B)
  ->M_Matrix<decltype(std::declval<T>()*std::declval<S>())>
  {
    assert(A.ncols()==B.ncols());
    assert(B.ncols()==A.nrows());
    typedef   decltype(std::declval<T>()*std::declval<S>()) R;
    if ((A.type()==M_Matrix<T>::ZERO)||
        (B.type()==M_Matrix<S>::ZERO))
      return M_Matrix<R>(B.nrows(),B.nrows(),M_Matrix<R>::ZERO);
    else if (A.isSymmetric())
      return product::quadraticForm_Symmetric_T(A,B);
    else
      return multTransp(B*A,B);
  }



  /**
     Scalar Multiplication assignment.
     @returns a reference to itself
     @post all the values of the matrix are multiplied by the value x
     */
  template<typename E, typename T>
  M_Matrix<E>& operator*=(M_Matrix<E>& itself, T x)
  {
    if (itself.type()==M_Matrix<E>::ZERO)
      return itself;
    else
      {
        for (size_t i=0; i<itself.size(); i++)
          itself[i]*=x;
        return itself;
      }
  }


  /**
     Scalar Division assignment.
     @returns a reference to itself
     @post all the values of the matrix are divided by the value x
     */
  template<typename E,typename T>
  M_Matrix<E>& operator/=(M_Matrix<E>& itself, T x)
  {
    if (itself.type()!=M_Matrix<E>::ZERO)
      for (size_t i=0; i<itself.size(); i++)
        itself[i]/=x;

    return itself;
  }



  //@}
  /** @name Aritmetic Assigment Operations between two Matrices
      */
  //@{

  namespace additive
  {


    using namespace Matrix_Unary_Transformations;

    template<typename T, typename S>
    M_Matrix<T>& sum_assigment(M_Matrix<T>& itself,
                               const M_Matrix<S>& other)
    {
      auto resType=result_of_sum(itself,other);
      if (resType!=itself.type())
        {
          itself=M_Matrix<T>(itself.nrows(),itself.ncols(),resType,itself);
        }
      auto itType=iterator_of_sum(itself,other);
      for (auto it=other.begin(itType); it!=other.end(); ++it)
        itself(it)+=*it;
      return itself;

    }


    template<typename T, typename S>
    M_Matrix<T>& substraction_assigment(M_Matrix<T>& itself,
                                        const M_Matrix<S>& other)
    {
      auto resType=result_of_sum(itself,other);
      if (resType!=itself.type())
        {
          itself=M_Matrix<T>(itself.nrows(),itself.ncols(),resType,itself);
        }
      auto itType=iterator_of_sum(itself,other);
      for (auto it=other.begin(itType); it!=other.end(); ++it)
        itself(it)-=*it;
      return itself;

    }




    template<typename T, typename S>
    auto sum_Operator( const M_Matrix<T>& itself,
                       const M_Matrix<S>& other)
    ->M_Matrix<decltype(std::declval<T>()+std::declval<S>())>
    {
      typedef decltype(std::declval<T>()+std::declval<S>()) R;
      if (itself.size()<other.size())
        {
          M_Matrix<R> out(other);
          sum_assigment(out,itself);
          return out;
        }
      else
        {
          M_Matrix<R> out(itself);
          sum_assigment(out,other);
          return out;
        }
    }


    template<typename T, typename S>
    auto sum_Operator( M_Matrix<T>&& itself,
                       M_Matrix<S>&& other)
    ->M_Matrix<decltype(std::declval<T>()+std::declval<S>())>
    {
      typedef decltype(std::declval<T>()+std::declval<S>()) R;
      if (itself.size()<other.size())
        {
          M_Matrix<R> out(other);
          sum_assigment(out,itself);
          return out;
        }
      else
        {
          M_Matrix<R> out(itself);
          sum_assigment(out,other);
          return out;
        }
    }


    template<typename T, typename S>
    auto substraction_Operator( const M_Matrix<T>& itself,
                                const M_Matrix<S>& other)
    ->M_Matrix<decltype(std::declval<T>()+std::declval<S>())>
    {
      typedef decltype(std::declval<T>()+std::declval<S>()) R;
      M_Matrix<R> out(itself);
      substraction_assigment(out,other);
      return out;
    }
    template<typename T, typename S>
    auto substraction_Operator( M_Matrix<T>&& itself,
                                M_Matrix<S>&& other)
    ->M_Matrix<decltype(std::declval<T>()+std::declval<S>())>
    {
      typedef decltype(std::declval<T>()+std::declval<S>()) R;
      M_Matrix<R> out(itself);
      substraction_assigment(out,other);
      return out;
    }




  }







  template<typename T, typename S, class AsOp>
  M_Matrix<T>& Element_Wise_Multiplicative_Assigment_Operator
  (AsOp op,M_Matrix<T>& itself,
   const M_Matrix<S>& other)
  {
    if (itself.type()==M_Matrix<T>::ZERO)
      {
        return itself;
      }
    else if (other.type()==M_Matrix<S>::ZERO)
      {
        itself=other;
        return itself;
      }
    else if (itself.type()==other.type())
      {
        assert(itself.size()==other.size());
        for (size_t i=0; i<itself.size(); i++)
          op(itself[i],other[i]);
        return itself;
      }
    else if ((itself.type()==M_Matrix<T>::DIAGONAL)
             ||(itself.type()==M_Matrix<T>::SCALAR_DIAGONAL))
      {
        assert(itself.nrows()==other.nrows());
        assert(itself.ncols()==other.ncols());
        switch(other.type())
          {
          case M_Matrix<S>::ZERO: // not possible
          case M_Matrix<S>::FULL:
          case M_Matrix<S>::SYMMETRIC:
          case M_Matrix<S>::SCALAR_FULL:
          case M_Matrix<S>::SCALAR_DIAGONAL:
          case M_Matrix<S>::DIAGONAL:
            {
              for (std::size_t i=0; i<itself.nrows(); ++i)
                op(itself(i,i),other(i,i));
              return itself;
            }

          }
      }

    else if ((other.type()==M_Matrix<T>::DIAGONAL)
             ||(other.type()==M_Matrix<T>::SCALAR_DIAGONAL))
      {
        assert(itself.nrows()==other.nrows());
        assert(itself.ncols()==other.ncols());
        switch(itself.type())
          {
          M_Matrix<T> out(other);
          for (std::size_t i=0; i<std::min(out.nrows(), out.ncols()); ++i)
            op(out(i,i),itself(i,i));
          itself=std::move(other);
          return itself;

          }
      }
    else if(itself.type()==M_Matrix<T>::FULL)
      {
        assert(itself.nrows()==other.nrows());
        assert(itself.ncols()==other.ncols());
        switch(other.type())
          {
          case M_Matrix<S>::ZERO: // not
          case M_Matrix<S>::SCALAR_DIAGONAL:  //done
          case M_Matrix<S>::DIAGONAL: // done
          case M_Matrix<S>::FULL:
          case M_Matrix<S>::SYMMETRIC:
          case M_Matrix<S>::SCALAR_FULL:
            {
              for (std::size_t i=0; i<itself.nrows(); ++i)
                for (std::size_t j=0; j<itself.ncols(); ++j)
                  op(itself(i,j),other(i,j));
              return itself;
            }

          }
      }
    else if(other.type()==M_Matrix<T>::FULL)
      {
        assert(itself.nrows()==other.nrows());
        assert(itself.ncols()==other.ncols());
        M_Matrix<T> out(other);
        for (std::size_t i=0; i<itself.nrows(); ++i)
          for (std::size_t j=0; j<itself.ncols(); ++j)
            op(out(i,j),itself(i,j));
        itself=std::move(out);
        return itself;
      }
    else if(itself.type()==M_Matrix<T>::SYMMETRIC)
      {
        assert(itself.nrows()==other.nrows());
        assert(itself.ncols()==other.ncols());
        switch(other.type())
          {
          case M_Matrix<S>::ZERO:  //not possible
          case M_Matrix<S>::SCALAR_DIAGONAL: //not possible
          case M_Matrix<S>::DIAGONAL: //not possible
          case M_Matrix<S>::FULL:  //should not land here
          case M_Matrix<S>::SYMMETRIC: // should not land here
          case M_Matrix<S>::SCALAR_FULL:
            {
              for (std::size_t i=0; i<itself.nrows(); ++i)
                for (std::size_t j=0; j<=i; ++j)
                  op(itself(i,j),other(i,j));
              return itself;
            }
          }
      }
    else  // scalar_full and the other is symmetric
      {
        assert(itself.nrows()==other.nrows());
        assert(itself.ncols()==other.ncols());
        M_Matrix<T> out(other);
        for (std::size_t i=0; i<itself.nrows(); ++i)
          for (std::size_t j=0; j<=i; ++j)
            op(out(i,j),itself(i,j));
        itself=std::move(out);
        return itself;
      }
  }

  //@}


  /**
 Matrix sum, element wise.
 @pre \p x.nrows()==ncols(y) x.ncols()==ncols(y)
 @return z, where z.nrows()=rows(x), z.ncols()=ncols(y) and z(i,j)= sum
 on k of x(i,k)*y(k,j)
 @warning it \c assert the preconditions
 */
  template<typename T>
  M_Matrix<T> operator+(const M_Matrix<T>& x,const M_Matrix<T>& y)
  {
    auto out=additive::sum_Operator(x,y);
    return out;
  }

  template<typename T, typename S>
  M_Matrix<T>& operator+=(M_Matrix<T>& x,const M_Matrix<S>& y)
  {
    return additive::sum_assigment(x,y);
  }

  template<typename T, typename S>
  M_Matrix<T>& operator-=(M_Matrix<T>& x,const M_Matrix<S>& y)
  {
    return additive::substraction_assigment(x,y);
  }

  /**
 Matrix sustraction, element wise.
 @pre \p x.nrows()==ncols(y) x.ncols()==ncols(y)
 @return z, where z.nrows()=rows(x), z.ncols()=ncols(y) and
 z(i,j)= sum on k of x(i,k)*y(k,j)
 @warning it \c assert the preconditions
 */

  template<typename T>
  M_Matrix<T> operator-(const M_Matrix<T>& x,const M_Matrix<T>& y)
  {
    return additive::substraction_Operator(x,y);
  }




  /** @name Aritmetic operation between a Matrix and a scalar (single element)
      */
  //@{

  /**
     Scalar Addition.
     @returns a copy of the matrix with its values summed by x
     */
  template<typename T>
  M_Matrix<T> operator+(const M_Matrix<T>& x,T t)
  {    // we build the M_Matrix result
    M_Matrix<T> z(x);
    z+=t;
    return z;
  }

  /**
     Scalar Addition reverse order.
     */
  template<typename T>
  M_Matrix<T> operator+(T t,const M_Matrix<T>& x)
  {
    return x+t;
  }

  /**
     Scalar Subtraction.
     @returns a copy of the matrix with its values substracted by x
     */
  template<typename T>
  M_Matrix<T> operator-(const M_Matrix<T>& x,T t)
  {    // we build the M_Matrix result
    M_Matrix<T> z(x);
    z-=t;
    return z;
  }

  /**
     Scalar Subtraction reverse order.
     */
  template<typename T>
  M_Matrix<T> operator-(T t,const M_Matrix<T>& x)
  {
    return x-t;
  }

  /**
     Scalar Multiplication.
     @returns a copy of the matrix with its values multiplied by the value x
     */
  template<typename E, typename T>
  auto operator*(const M_Matrix<E> & x,T t)
  ->M_Matrix<decltype(std::declval<E>()*std::declval<T>())>
  {    // we build the M_Matrix result
    typedef decltype(std::declval<E>()*std::declval<T>()) R;
    M_Matrix<R> out(x);
    for (std::size_t i=0; i<x.size(); ++i)
      out[i]*=t;
    return out;
  }

  /**
     Scalar Multiplication reverse order.
     */

  template<typename T>
  M_Matrix<T> operator*(T t,const M_Matrix<T>& x)
  {
    return x*t;
  }


  /**
     Scalar Division.
     @returns a copy of the matrix with its values divided by x
     @returns a matrix of real numbers
 */
  template<typename T>
  M_Matrix<double> operator/(const M_Matrix<T>& x,T t)
  {    // we build the M_Matrix result
    M_Matrix<double> z(x);
    z/=double(t);
    return z;
  }

  /**
     Division by inhomogeneus types

     */

  template<typename T,typename S>
  M_Matrix<double> operator/(const M_Matrix<T>& x,S t)
  {    // we build the M_Matrix result
    M_Matrix<double> z(x);
    z/=double(t);
    return z;
  }



  /**
     Scalar Division reverse order.
     */
  template<typename T>
  M_Matrix<double> operator/(T t,const M_Matrix<T>& x)
  {
    M_Matrix<double> out(x.nrows(),x.ncols());
    for (std::size_t i=0;i<x.size();i++)
      out[i]=double(t)/double(x[i]);

    return out;
  }


  //@}
  /**
 @name Aritmetic operations applied between two Matrices
  */
  //@{





  /**
 Multiplication of the elements of two matrices.
  @pre \p x.nrows()==ncols(y) x.ncols()==ncols(y)
  @return z, where z.nrows()=rows(x), z.ncols()=ncols(y)
 and z(i,j)= sum on k of x(i,k)*y(k,j)
  @warning it \c assert the preconditions
 */
  template<typename T, typename S>
  M_Matrix<T>& elemMult(M_Matrix<T>& x,const M_Matrix<S>& y)
  {
    return Element_Wise_Multiplicative_Assigment_Operator([](T& a, const S& b){a*=b; return a;},x,y);
  }

  /**
 Division of the elements of two matrices.
  @pre \p x.nrows()==ncols(y) x.ncols()==ncols(y)
  @return z, where z.nrows()=rows(x), z.ncols()=ncols(y)
 and z(i,j)= sum on k of x(i,k)*y(k,j)
  @warning it \c assert the preconditions
 */
  template<typename T,typename S>
  M_Matrix<double> elemDiv(const M_Matrix<T>& x,const M_Matrix<S>& y)
  {
    M_Matrix<double> z(x);
    for (size_t i=0; i<z.nrows(); i++)
      for (size_t j=0; j<z.ncols(); j++)
        z(i,j)/=double(y(i,j));
    return z;
  }

  /**
 Safe Division of the elements of two matrices.
  @pre \p x.nrows()==ncols(y) x.ncols()==ncols(y)
  @return z, where z.nrows()=rows(x), z.ncols()=ncols(y)
 and z(i,j)= sum on k of x(i,k)*y(k,j)
  @warning it \c assert the preconditions
 */
  template<typename T,typename S>
  M_Matrix<double> elemDivSafe(const M_Matrix<T>& x,const M_Matrix<S>& y)
  {
    M_Matrix<double> z(x);
    for (size_t i=0; i<z.nrows(); i++)
      for (size_t j=0; j<z.ncols(); j++)
        if (y(i,j)!=0)
          z(i,j)/=double(y(i,j));
        else
          z(i,j)=0;
    return z;
  }
  // @}


  template<typename T>
  M_Matrix<double> operator<<(const M_Matrix<double>& A, const M_Matrix<T>& B)
  {
    M_Matrix<double> out
        (std::max(A.nrows(), B.nrows()), A.ncols()+B.ncols(), std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i=0; i<A.nrows();++i)
      {
        for (std::size_t j=0; j<A.ncols(); ++j)
          out(i,j)=A(i,j);
      }
    for (std::size_t i=0; i<B.nrows();++i)
      {
        for (std::size_t j=0; j<B.ncols(); ++j)
          out(i,j)=B(i,A.ncols()+j);
      }

    return out;

  }





  inline M_Matrix<double> xdiagXT(const M_Matrix<double>& x, const M_Matrix<double> Cdiag)
  {
    M_Matrix<double> o(x.nrows(), x.nrows(),0.0);
    for (std::size_t i=0;  i<x.nrows(); ++i)
      for (std::size_t j=0; j<x.nrows(); ++j)
        for (std::size_t k=0; k<x.ncols(); ++k)
          o(i,j)+=Cdiag[k]*x(i,k)*x(j,k);
    return o;
  }



  inline M_Matrix<double> MultDiag(const M_Matrix<double> &x, const M_Matrix<double> d)
  {
    M_Matrix<double> o(x.nrows(), x.ncols());
    for (std::size_t i=0;  i<x.nrows(); ++i)
      for (std::size_t j=0; j<x.ncols(); ++j)
        o(i,j)=x(i,j)*d[j];
    return o;
  }


  inline M_Matrix<double> DiagMult( const M_Matrix<double> d,const M_Matrix<double> &x)
  {
    M_Matrix<double> o(x.nrows(), x.ncols());
    for (std::size_t i=0;  i<x.nrows(); ++i)
      for (std::size_t j=0; j<x.ncols(); ++j)
        o(i,j)=x(i,j)*d[i];
    return o;
  }

  template<typename T>
  M_Matrix<T> diag_landa(const M_Matrix<T>& x,double landa)
  {
    double landa1=landa+1;
    M_Matrix<T> diagM(x);
    for (size_t i=0; i<x.nrows(); ++i)
      diagM(i,i)*=landa1;
    return diagM;

  }


  template<typename V>
  M_Matrix<double> JTd2J(const M_Matrix<double>& J, const V& D2)
  {
    assert(J.nrows()==D2.size());
    M_Matrix<double> out(J.ncols(),J.ncols());
    M_Matrix<double> Jc=Transpose(J);
    for (std::size_t i=0; i<Jc.nrows(); ++i)
      for (std::size_t j=0; j<Jc.ncols(); ++j)
        Jc(i,j)*=D2[j];
    out=Jc*J;
    return out;

  }

}


using namespace Vector_Binary_Transformations;
using namespace Container_Unary_Functions;
using namespace Container_Binary_Transformations;
using namespace Matrix_Generators;
using namespace Matrix_Stream_Operators;
using namespace Matrix_Unary_Predicates;
using namespace Matrix_Unary_Size_Functions;
using namespace Matrix_Unary_Functions;
using namespace Matrix_Unary_Transformations;
using namespace Matrix_Binary_Predicates;
using namespace Matrix_Binary_Transformations;


//@}



template<class T>
T
sum(const std::vector<T>& x)
{
    T o{};
    for (unsigned i=0; i<x.size(); ++i)
        o+=x[i];

  return o;
}










#endif // MATRIX_H
