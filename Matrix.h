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

#include <iterator> // std::iterator, std::input_iterator_tag

#include <iostream>
#include <algorithm>
#include <memory>

#include "mySerializer.h"
#include "mytests.h"
#include "myoptional.h"
#include "myderivatives.h"

template <typename T1, typename T2>
std::pair<T1, T2> &operator+=(std::pair<T1, T2> &x,
                              const std::pair<T1, T2> &other) {
  x.first += other.first;
  x.second += other.second;
  return x;
}

inline double average(double x, double y) { return 0.5 * (x + y); }

#ifdef TEST_MODE
#include <gtest/gtest.h>
#include <gmock/gmock-matchers.h>


TEST(teste, teste)
{
    EXPECT_EQ(1, 11);
    ASSERT_THAT(0, testing::Eq(0));
}
#endif

inline double sqr(double x) { return x * x; }

template <typename T> inline T sqr(const T &x);
namespace Vector_Unary_Index_Function {

template <typename T> std::size_t i_max(const std::vector<T> &x) {
  std::size_t j = 0;
  for (std::size_t i = 1; i < x.size(); ++i)
    if (x[i] > x[j])
      j = i;
  return j;
}
} // namespace Vector_Unary_Index_Function

namespace Container_Unary_Functions {

template <class V> double mean(const V &x) {
  double sum = 0;
  for (std::size_t i = 0; i < x.size(); ++i)
    sum += x[i];
  return sum / x.size();
}

template <class V> double stddev(const V &x) {
  double sum = 0;
  double sum2 = 0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    sum += x[i];
    sum2 += sqr(x[i]);
  }
  return std::sqrt(sum2 / x.size() - sqr(sum / x.size()));
}

} // namespace Container_Unary_Functions
template <template <typename...> class V, class... Ts>
bool all(const V<Ts...> &x) {
  for (std::size_t i = 0; i < x.size(); ++i)
    if (!x[i])
      return false;
  return true;
}




template <template <typename...> class V, class... Ts>
bool any(const V<bool, Ts...> &x) {
  for (std::size_t i = 0; i < x.size(); ++i)
    if (x[i])
      return true;
  return false;
}

namespace Container_Binary_Transformations {

namespace additive {

template <template <typename...> class M, typename... Ts>
M<Ts...> &operator_additive_assigment(M<Ts...> &itself, const M<Ts...> &x) {
  for (size_t i = 0; i < itself.size(); i++)
    itself[i] += x[i];
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

template <template <typename...> class M, typename T, typename... Ts>
M<Ts...> &operator_additive_assigment(M<Ts...> &itself, T x) {
  for (size_t i = 0; i < itself.size(); i++)
    itself[i] += x;
  return itself;
}
} // namespace additive
template <typename T>
std::vector<T> &operator+=(std::vector<T> &itself,
                           const std::vector<T> &other) {
  return additive::operator_additive_assigment(itself, other);
}


template <typename T, typename Number>
std::vector<decltype(std::declval<T>()*std::declval<Number>())> &operator*(const std::vector<T> &itself, Number x) {
    std::vector<decltype(std::declval<T>()*std::declval<Number>())> out(itself.size());
    for(std::size_t i=0; i<out.size(); ++i)
        out[i]=itself[i]*x;
    return out;
  }
  template <typename T, typename Number>
  std::vector<decltype(std::declval<T>()*std::declval<Number>())> &operator/(const std::vector<T> &itself, Number x) {
      std::vector<decltype(std::declval<T>()/std::declval<Number>())> out(itself.size());
      for(std::size_t i=0; i<out.size(); ++i)
          out[i]=itself[i]/x;
      return out;
  }




template <typename T>
std::vector<T> &operator+=(std::vector<T> &itself, const T &other) {
  return additive::operator_additive_assigment(itself, other);
}

} // namespace Container_Binary_Transformations

namespace lapack {
extern "C" void dgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV,
                        int *INFO);

extern "C" double dlange_(char *NORM, int *M, int *N, double *A, int *LDA,
                          double *WORK);

extern "C" void dgecon_(char *NORM, int *N, double *A, int *LDA, double *ANORM,
                        double *RCOND, double *WORK, int *IWORK, int *INFO);

extern "C" void dsycon_(char *UPLO, int *N,
                        double * /* precision, dimension( lda, * ) */ A,
                        int *LDA, int * /*, dimension( * ) */ IPIV,
                        double *ANORM, double *RCOND,
                        double * /* dimension( * )  */ WORK,
                        int * /* dimension( * ) */ IWORK, int *INFO);

extern "C" void dsytrf_(char *UPLO, int *N, double *A, int *LDA, int *IPIV,
                        double *WORK, int *LWORK, int *INFO);

extern "C" void dgetri_(int *n, double *B, int *dla, int *ipiv, double *work1,
                        int *lwork, int *info);

extern "C" void dsytri_(char *UPLO, int *N, double * /*dimension( lda, * )*/ A,
                        int *LDA, int * /* dimension( * ) */ IPIV,
                        double * /*dimension( * )*/ WORK, int *INFO);

extern "C" void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K,
                       double *ALPHA, double *A, int *LDA, double *B, int *LDB,
                       double *BETA, double *C, int *LDC);

extern "C" void dsymm_(char *SIDE, char *UPLO, int *M, int *N, double *ALPHA,
                       double * /*, dimension(lda,*)*/ A, int *LDA,
                       double * /*, dimension(ldb,*)*/ B, int *LDB,
                       double *BETA, double * /*, dimension(ldc,*) */ C,
                       int *LDC);

extern "C" void dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda,
                       double *wr, double *wi, double *vl, int *ldvl,
                       double *vr, int *ldvr, double *work, int *lwork,
                       int *info);

extern "C" void
dgeevx_(char *BALANC, char *JOBVL, char *JOBVR, char *SENSE, int *N,
        double * /* precision, dimension( lda, * ) */ A, int *LDA,
        double * /* precision, dimension( * ) */ WR,
        double * /* precision, dimension( * ) */ WI,
        double * /* precision, dimension( ldvl, * ) */ VL, int *LDVL,
        double * /* precision, dimension( ldvr, * ) */ VR, int *LDVR, int *ILO,
        int *IHI, double * /* precision, dimension( * ) */ SCALE, double *ABNRM,
        double * /* precision, dimension( * ) */ RCONDE,
        double * /* precision, dimension( * ) */ RCONDV,
        double * /* precision, dimension( * ) */ WORK, int *LWORK,
        int * /*dimension( * ) */ IWORK, int *INFO);

extern "C" void dsyevx_(char *JOBZ, char *RANGE, char *UPLO, int *N,
                        double * /*precision, dimension(lda, *) */ A, int *LDA,
                        double *VL, double *VU, int *IL, int *IU,
                        double *ABSTOL, int *M,
                        double * /*precision, dimension(*) */ W,
                        double * /*
                         dimension(ldz, *) */
                            Z,
                        int *LDZ, double * /*precision, dimension(*) */ WORK,


                       int *LWORK, int * /*, dimension(*) */ IWORK,
                        int * /*, dimension(*) */ IFAIL, int *INFO);


extern "C" void ddisna_(char *JOB, int * M, int *N,double * D, double *SEP,int * INFO );


extern "C" void
dgesdd_(char *JOBZ, int *M, int *N,
        double * /* precision, dimension( lda, * ) */ A, int *LDA,
        double * /*precision, dimension( * ) */ S,
        double * /*precision, dimension( ldu, * ) */ U, int *LDU,
        double * /* precision, dimension( ldvt, * ) */ VT, int *LDVT,
        double * /* precision, dimension( * ) */ WORK, int *LWORK,
        int * /*, dimension( * ) */ IWORK, int *INFO);

} // namespace lapack

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

enum class Matrix_TYPE {
  ZERO,
  FULL,
  SYMMETRIC,
  DIAGONAL,
  SCALAR_FULL,
  SCALAR_DIAGONAL
};
template <typename T> class M_Matrix ;
template <class> struct is_Matrix : public std::false_type {};

template <typename T> struct is_Matrix<M_Matrix<T>> : public std::true_type {};

template <class M> constexpr static bool is_Matrix_v = is_Matrix<M>::value;



template <typename T> class M_Matrix {
public:
  constexpr static auto className =
      my_static_string("matrix_") + my_trait<T>::className;

  static std::pair<std::size_t, std::size_t>
  pos_to_ij_Symmetric(std::size_t k) {
    std::size_t i = std::floor((std::sqrt(8 * k + 1) - 1) / 2);
    std::size_t j = k - (i * (i + 1)) / 2;
    assert(k == i * (i + 1) / 2 + j);
    return std::pair(i, j);
  }
  static std::pair<std::size_t, std::size_t> pos_to_ij_Full(std::size_t k,
                                                            std::size_t n) {
    std::size_t i = k / n;
    std::size_t j = k - i * n;
    assert(k == i * n + j);
    return std::pair(i, j);
  }

  static std::pair<std::size_t, std::size_t>
  pos_to_ij_Symmetric(std::size_t k, std::size_t n) {
    std::size_t i, j;
    std::tie(i, j) = pos_to_ij_Symmetric(n * (n + 1) / 2 - 1 - k);

    return std::pair(n - 1 - i, n - 1 - j);
  }

  class MyConstIterator;
  enum ITER_TYPE { FULL_IT, LT_IT, UT_IT, DIAG_IT, SCALAR_IT, ZERO_IT };

  // enum TYPE{ZERO,FULL,SYMMETRIC,DIAGONAL,SCALAR_FULL,SCALAR_DIAGONAL};
  typedef Matrix_TYPE TYPE;

  static ITER_TYPE iter_type(TYPE type) {
    switch (type) {
    case Matrix_TYPE::ZERO:
      return ZERO_IT;
    case Matrix_TYPE::FULL:
      return FULL_IT;
    case Matrix_TYPE::SYMMETRIC:
      return LT_IT;
    case Matrix_TYPE::DIAGONAL:
      return DIAG_IT;
    case Matrix_TYPE::SCALAR_DIAGONAL:
    case Matrix_TYPE::SCALAR_FULL:
      return SCALAR_IT;
    default:
      return {};
    }
  }
  static ITER_TYPE const_iter_type(TYPE type) {
    switch (type) {
    case Matrix_TYPE::FULL:
    case Matrix_TYPE::SYMMETRIC:
    case Matrix_TYPE::SCALAR_FULL:
      return FULL_IT;
    case Matrix_TYPE::DIAGONAL:
    case Matrix_TYPE::SCALAR_DIAGONAL:
      return DIAG_IT;
    case Matrix_TYPE::ZERO:
    default:
      return ZERO_IT;
    }
  }

  class MyIterator : public std::iterator<std::input_iterator_tag, T> {
    M_Matrix<T> &m;
    std::size_t i_;
    std::size_t j_;
    ITER_TYPE type_;

  public:
    friend class MyConstIterator;
    ITER_TYPE type() const { return type_; }
    std::size_t iRow() const { return i_; }
    std::size_t jCol() const { return j_; }

    MyIterator(M_Matrix<T> &x)
        : m(x), i_(0), j_(0), type_(iter_type(x.type())) {}

    MyIterator(M_Matrix<T> &x, ITER_TYPE t) : m(x), i_(0), j_(0), type_(t) {}

    MyIterator(M_Matrix<T> &x, std::size_t i, std::size_t j)
        : m(x), i_(i), j_(j), type_(iter_type(x.type())) {}

    MyIterator(const MyIterator &mit)
        : m(mit.m), i_(mit.i_), j_(mit.j_), type_(mit.type_) {}
    MyIterator &operator++() {
      ++j_;
      switch (type()) {
      case FULL_IT:
        if (j_ >= m.ncols()) {
          j_ = 0;
          ++i_;
        }
        break;
      case LT_IT:
        if (j_ > i_) {
          j_ = 0;
          ++i_;
        }
        break;
      case UT_IT:
        if (j_ >= m.ncols()) {
          ++i_;
          if (i_ >= m.nrows())
            j_ = 0;
          else
            j_ = i_;
        }
        break;

      case DIAG_IT:
        if (j_ >= m.ncols()) {
          j_ = 0;
          i_ = m.nrows();
        } else {
          i_ = j_;
        }
        break;
      case SCALAR_IT:
      case ZERO_IT:
        if (j_ > 0) {
          j_ = 0;
          i_ = m.nrows();
        }
        break;
      }
      return *this;
    }
    MyIterator operator++(int) {
      MyIterator tmp(*this);
      operator++();
      return tmp;
    }
    bool operator==(const MyIterator &rhs) {
      if (i_ != rhs.i_)
        return false;
      else if (j_ != rhs.j_)
        return false;
      else
        return true;
    }
    bool operator!=(const MyIterator &rhs) { return !(*this == rhs); }
    T &operator*() { return m(i_, j_); }
    const T &operator*() const { return m(i_, j_); }
  };

  class MyConstIterator : public std::iterator<std::input_iterator_tag, T> {
    const M_Matrix<T> &m;
    std::size_t i_;
    std::size_t j_;
    ITER_TYPE type_;

  public:
    std::size_t iRow() const { return i_; }
    std::size_t jCol() const { return j_; }
    ITER_TYPE type() const { return type_; }

    MyConstIterator(const M_Matrix<T> &x)
        : m(x), i_(0), j_(0), type_(const_iter_type(x.type())) {}

    MyConstIterator(const M_Matrix<T> &x, ITER_TYPE t)
        : m(x), i_(0), j_(0), type_(t) {}

    MyConstIterator(const M_Matrix<T> &x, std::size_t i, std::size_t j)
        : m(x), i_(i), j_(j), type_(const_iter_type(x.type())) {}

    MyConstIterator(const MyConstIterator &mit)
        : m(mit.m), i_(mit.i_), j_(mit.j_), type_(mit.type_) {}

    MyConstIterator(const MyIterator &mit)
        : m(mit.m), i_(mit.i_), j_(mit.j_), type_(mit.type_) {}

    MyConstIterator &operator++() {
      ++j_;
      switch (type()) {
      case FULL_IT:
        if (j_ >= m.ncols()) {
          j_ = 0;
          ++i_;
        }
        break;
      case LT_IT:
        if (j_ > i_) {
          j_ = 0;
          ++i_;
        }
        break;
      case UT_IT:
        if (j_ >= m.ncols()) {
          ++i_;
          if (i_ >= m.nrows())
            j_ = 0;
          else
            j_ = i_;
        }
        break;

      case DIAG_IT:
        if (j_ >= m.ncols()) {
          j_ = 0;
          i_ = m.nrows();
        } else {
          i_ = j_;
        }
        break;
      case SCALAR_IT:
      case ZERO_IT:
        if (j_ > 0) {
          j_ = 0;
          i_ = m.nrows();
        }
        break;
      }
      return *this;
    }
    MyConstIterator operator++(int) {
      MyConstIterator tmp(*this);
      operator++();
      return tmp;
    }
    bool operator==(const MyConstIterator &rhs) {
      if (i_ != rhs.i_)
        return false;
      else if (j_ != rhs.j_)
        return false;
      else
        return true;
    }
    bool operator!=(const MyConstIterator &rhs) { return !(*this == rhs); }
    const T &operator*() const { return m(i_, j_); }
  };

  typedef MyIterator iterator;
  typedef MyConstIterator const_iterator;

  iterator begin() {
    if (type() != Matrix_TYPE::ZERO) {
      MyIterator out(*this);
      return out;
    } else
      return end();
  }

  iterator begin(ITER_TYPE t) {
    if ((type() != Matrix_TYPE::ZERO) && (t != ZERO_IT)) {
      MyIterator out(*this, t);
      return out;
    } else
      return end();
  }

  iterator end() {
    MyIterator out(*this, nrows(), 0);
    return out;
  }

  const_iterator begin() const {
    if (type() != Matrix_TYPE::ZERO) {
      MyConstIterator out(*this);
      return out;
    } else
      return end();
  }

  const_iterator begin(ITER_TYPE t) const {
    if ((type() != Matrix_TYPE::ZERO) && (t != ZERO_IT)) {
      const_iterator out(*this, t);
      return out;
    } else
      return end();
  }
  const_iterator end() const {
    MyConstIterator out(*this, nrows(), 0);
    return out;
  }
  double getvalue() const {
    assert(size() == 1);
    return (*this)[0];
  }

  M_Matrix()=default;
 //     :type_{Matrix_TYPE::ZERO},_nrows{0ul},_ncols{0ul},_data{},zero_{}{} -> type_=ZERO

  M_Matrix(const M_Matrix<T> &sample)
      :type_{sample.type_},_nrows{sample._nrows},_ncols{sample._ncols},_data{sample._data},zero_{sample.zero_}{}

  M_Matrix(M_Matrix<T> &&sample)
      :type_{sample.type_},_nrows{sample._nrows},_ncols{sample._ncols},_data{std::move(sample._data)},zero_{std::move(sample.zero_)}{
      sample.type_=Matrix_TYPE::ZERO;
      sample._nrows=0; sample._ncols=0;
  }

  M_Matrix(std::size_t nrows, std::size_t ncols)
      : type_(Matrix_TYPE::FULL), _nrows(nrows), _ncols(ncols),
        _data(nrows * ncols) {}

  template <typename S, std::enable_if_t<std::is_convertible_v<S, T>, int> = 0>
  M_Matrix(const M_Matrix<S> &sample)
      : type_(TYPE(sample.type())), _nrows(sample.nrows()),
        _ncols(sample.ncols()), _data(sample.size())

  {
    for (auto it = begin(); it != end(); ++it)
      *it = sample(it);
  }
  M_Matrix(const std::vector<std::vector<T>> &sample)
      : type_(Matrix_TYPE::FULL), _nrows(sample.size()),
        _ncols(sample.size()>0?sample[0].size():0), _data(sample.size()>0?sample.size()*sample[0].size():0) {
    for (std::size_t i = 0; i < _nrows; ++i)
      for (std::size_t j = 0; j < _ncols; ++j)
        (*this)(i, j) = sample[i][j];
  }

  M_Matrix(TYPE t, const std::vector<std::vector<T>> &sample)
      : type_(t), _nrows(sample.size()), _ncols(sample[0].size()),
        _data(getSize(t, sample.size(), sample[0].size())) {
    for (auto it = begin(); it != end(); ++it)
      (*it) = sample[it.iRow()][it.jCol()];
  }

  M_Matrix(TYPE t, const M_Matrix<T> &sample)
      : type_(t), _nrows(sample.nrows()), _ncols(sample.ncols()),
        _data(getSize(t, sample.nrows(), sample.ncols())) {
    for (auto it = begin(); it != end(); ++it)
      (*it) = sample(it);
    init_zero();
  }

  template <typename S>
  M_Matrix(const std::vector<std::vector<S>> &sample)
      : type_(Matrix_TYPE::FULL), _nrows(sample.size()),
        _ncols(sample[0].size()), _data(sample.size() * sample[0].size()) {
    for (std::size_t i = 0; i < sample.size(); ++i)
      (*this)[i] = sample[i];
  }

  M_Matrix(std::size_t nrows_, std::size_t ncols_, TYPE t, T data)
      : type_(t), _nrows(nrows_), _ncols(ncols_),
        _data(getSize(t, nrows_, ncols_), data) {init_zero();}


  M_Matrix(std::size_t nrows_, std::size_t ncols_, TYPE t)
      : type_{t}, _nrows(nrows_), _ncols(ncols_),
        _data(getSize(t, nrows_, ncols_)) {init_zero();}

  template <typename Enum, std::enable_if_t<std::is_enum_v<Enum>, int> = 0>
  M_Matrix(std::size_t nrows_, std::size_t ncols_, Enum t, T x)
      : type_{t}, _nrows(nrows_), _ncols(ncols_),
        _data(getSize(type_, nrows_, ncols_), x) {init_zero();}

  M_Matrix(std::size_t nrows_, std::size_t ncols_, std::vector<T>&& data)
      : type_(Matrix_TYPE::FULL), _nrows(nrows_), _ncols(ncols_), _data(std::move(data)) {}

  M_Matrix(std::size_t nrows_, std::size_t ncols_, const std::vector<T>& data)
      : type_(Matrix_TYPE::FULL), _nrows(nrows_), _ncols(ncols_), _data(data) {}


  M_Matrix(std::size_t nrows_, std::size_t ncols_, const M_Matrix<T> &data)
      : type_(Matrix_TYPE::FULL), _nrows(nrows_), _ncols(ncols_),
        _data(nrows_ * ncols_) {
    assert(data.nrows() == nrows());
    assert(data.ncols() == ncols());
    for (size_t i = 0; i < nrows(); ++i)
      for (std::size_t j = 0; j < ncols(); ++j)
        (*this)(i, j) = data(i, j);
  }


  template <typename S>
  M_Matrix(std::size_t nrows_, std::size_t ncols_, const M_Matrix<S> &data)
      : type_(Matrix_TYPE::FULL), _nrows(nrows_), _ncols(ncols_),
        _data(nrows_ * ncols_)
  //   _data(new T[_ncells])
  {

    assert(data.nrows() == nrows());
    assert(data.ncols() == ncols());
    for (size_t i = 0; i < nrows(); ++i)
      for (std::size_t j = 0; j < ncols(); ++j)
        (*this)(i, j) = data(i, j);
  }

  template <typename S>
  M_Matrix(std::size_t nrows_, std::size_t ncols_, TYPE type,
           const M_Matrix<S> &data)
      : type_(type), _nrows(nrows_), _ncols(ncols_), _data(getSize(type_, nrows_, ncols_))
  //   _data(new T[_ncells])
  {

    assert(data.nrows() == nrows());
    assert(data.ncols() == ncols());
    for (auto it = begin(); it != end(); ++it)
      (*it) = data(it);
init_zero();
  }

  M_Matrix(std::size_t nrows_, std::size_t ncols_, T data)
      : type_(Matrix_TYPE::FULL), _nrows(nrows_), _ncols(ncols_),
        //   _data(new T[_ncells])
        _data(nrows_ * ncols_, data) {}

  M_Matrix &operator=(const M_Matrix<T> &sample) = default;
  M_Matrix &operator=(M_Matrix<T> &&sample) = default;

  ~M_Matrix() {
    //  if (_data>0)
    //	delete [] _data;
  }

  template <class F, typename...Ts> M_Matrix<T> apply(const F &f, const M_Matrix<Ts>&...xs) const {
    M_Matrix<T> out(nrows(), ncols(), type());
    for (std::size_t i = 0; i < size(); ++i)
      out[i] = std::invoke(f, (*this)[i],xs[i]...);
    return out;
  }

//  template <class F, typename...Ts> M_Matrix<T> apply(const F &f, const M_Matrix<Ts>&...xs) && {
//      for (std::size_t i = 0; i < size(); ++i)
//          (*this)[i] = std::invoke(f, (*this)[i],xs[i]...);
//      return *this;
//  }


  template <class F> auto auto_apply(const F &f) const {
    typedef std::invoke_result_t<F, T &> R;
    M_Matrix<R> out(nrows(), ncols(), type());
    for (std::size_t i = 0; i < size(); ++i)
      out[i] = std::invoke(f, (*this)[i]);
    return out;
  }

  template <class F, typename...Ts> M_Matrix &transform(const F &f, const M_Matrix<Ts>...xs) {
    for (std::size_t i = 0; i < size(); ++i)
      std::invoke(f, (*this)[i], xs[i]...);
    return *this;
  }

  M_Matrix<T> &operator=(T X) {
    for (std::size_t i = 0; i < size(); ++i)
      _data[i] = X;
    return *this;
  }

  size_t size() const { return _data.size(); }

  size_t nrows() const { return _nrows; }

  size_t ncols() const { return _ncols; }

  T &operator[](std::size_t n) { return _data[n]; }
  bool empty() const { return _data.empty(); }

  const T &operator[](std::size_t n) const { return _data[n]; }

  std::pair<std::size_t, std::size_t> pos_to_ij(std::size_t i) const {
    switch (type_) {
    case Matrix_TYPE::FULL: {
      return pos_to_ij_Full(i, _ncols);
    }
    case Matrix_TYPE::SYMMETRIC:
      return pos_to_ij_Symmetric(i);
    case Matrix_TYPE::DIAGONAL:
    case Matrix_TYPE::SCALAR_DIAGONAL:
      return {i, i};
    case Matrix_TYPE::SCALAR_FULL:
      return {0, 0};
    case Matrix_TYPE::ZERO:
    default:
      assert(false);
      return {0, 0};
    }
  }

  void set(std::size_t i, std::size_t j, const T& val){ (*this)(i,j)=val;}
  void set(std::size_t i,  const T& val){ (*this)[i]=val;}
  void add(std::size_t i, std::size_t j, const T& val){ (*this)(i,j)+=val;}

  T &operator()(std::size_t i, std::size_t j) {
    assert(i < nrows());
    assert(j < ncols());
    switch (type_) {
    case Matrix_TYPE::FULL:
      return (*this)[i * _ncols + j];
    case Matrix_TYPE::SYMMETRIC:
      if (i >= j)
        return (
            *this)[(i * (i + 1)) / 2 + j]; // stores at lower triangular portion
      else
        return (*this)[(j * (j + 1)) / 2 + i];
    case Matrix_TYPE::DIAGONAL:
      assert(i == j);
      return (*this)[i];
    case Matrix_TYPE::SCALAR_DIAGONAL:
      assert(i == j);
      return (*this)[0];
    case Matrix_TYPE::SCALAR_FULL:
      return (*this)[0];
    case Matrix_TYPE::ZERO:
      assert(false);
      return (*this)[0];
    default:
      assert(false);
      return (*this)[0];
    }
  }

  T const &operator()(std::size_t i, std::size_t j) const {
    if (i >= nrows())
      assert(false);
    assert(j < ncols());
    switch (type_) {
    case Matrix_TYPE::FULL:
      return (*this)[i * _ncols + j];
    case Matrix_TYPE::SYMMETRIC:
      if (i >= j)
        return (
            *this)[i * (i + 1) / 2 + j]; // stores at lower triangular portion
      else
        return (*this)[j * (j + 1) / 2 + i];
    case Matrix_TYPE::DIAGONAL:
      if (i == j)
        return (*this)[i];
      else
          return zero(i,j);
    case Matrix_TYPE::SCALAR_DIAGONAL:
      if (i == j)
        return (*this)[0];
      else
          return zero(0,0);
    case Matrix_TYPE::SCALAR_FULL:
      return (*this)[0];
    case Matrix_TYPE::ZERO:
        return zero(0,0);
    default:
        return zero(0,0);
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

  M_Matrix<T> &operator()(std::size_t iRow, std::string /*dummy*/,
                          const M_Matrix<T> &newValues) {
    assert(type_ == Matrix_TYPE::FULL);
    assert(iRow < nrows());              // number of rows
    assert(newValues.size() == ncols()); // number of columns
    for (std::size_t j = 0; j < std::min(ncols(), size()); j++)
      this->operator()(iRow, j) = newValues[j];
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
  M_Matrix<T> &operator()(const std::string /*dummy*/, std::size_t jColumn,
                          const V &newValues) {
    assert(type_ == Matrix_TYPE::FULL);

    for (std::size_t i = 0; i < std::min(nrows(), newValues.size()); i++)
      this->operator()(i, jColumn) = newValues[i];
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

  M_Matrix<T> operator()(std::size_t iRow, const std::string /*dummy*/
                         ) const {
    assert(type_ == Matrix_TYPE::FULL);
    M_Matrix<T> out(1, ncols());
    for (std::size_t j = 0; j < ncols(); j++)
      out[j] = this->operator()(iRow, j);
    return out;
  }

  template <class Iter> T &operator()(const Iter &it) {
    return (*this)(it.iRow(), it.jCol());
  }
  template <class Iter> const T &operator()(const Iter &it) const {
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

  M_Matrix<T> operator()(std::string /*dummy*/, std::size_t jColumn) const {
    assert(type_ == Matrix_TYPE::FULL);
    M_Matrix<T> out(nrows(), 1);
    for (std::size_t i = 0; i < nrows(); i++)
      out[i] = (*this)(i, jColumn);
    return out;
  }

  void clear() {
    type_ = Matrix_TYPE::ZERO;
    _nrows = 0;
    _ncols = 0;
    _data.clear();
  }

  std::vector<T> toVector() const { return _data; }

  M_Matrix<T> toVector_of_Rows() const { return M_Matrix<T>(size(), 1, _data); }
  M_Matrix<T> toVector_of_Cols() const { return M_Matrix<T>(1, size(), _data); }

  std::vector<std::vector<T>> toMatrix() const {
    std::vector<std::vector<T>> out(nrows(), std::vector<T>(ncols()));
    for (std::size_t i = 0; i < nrows(); ++i)
      for (std::size_t j = 0; j < ncols(); ++j)
        out[i][j] = (*this)(i, j);
    return out;
  }

  /**
   Returns a custom sized Matrix filled with ones
  @post (ones(n,m))(i,j)==T(1)
   */

  static M_Matrix unpackForLapack(const M_Matrix<double> &packedMatrix,
                                  char UPLO = 'L') {
    switch (packedMatrix.type()) {
    case Matrix_TYPE::SYMMETRIC: {
      M_Matrix<double> out(packedMatrix.nrows(), packedMatrix.ncols());
      if (UPLO == 'L') {
        for (std::size_t i = 0; i < out.nrows(); ++i)
          for (std::size_t j = 0; j <= i; ++j)
            out(i, j) = packedMatrix(i, j);
        return out;
      } else {
        for (std::size_t i = 0; i < out.nrows(); ++i)
          for (std::size_t j = i; j <= out.ncols(); ++j)
            out(i, j) = packedMatrix(i, j);
        return out;
      }
    }
    default:
      return M_Matrix<double>(packedMatrix);
    }
  }

  M_Matrix full() const {
    M_Matrix out(nrows(), ncols());
    for (std::size_t i = 0; i < nrows(); i++)
      for (std::size_t j = 0; j < ncols(); ++j)
        out(i, j) = (*this)(i, j);
    return out;
  }

  static std::size_t getSize(TYPE t, std::size_t nrows, std::size_t ncols) {
    switch (t) {
    case Matrix_TYPE::ZERO:
      return 0;
    case Matrix_TYPE::FULL:
      return nrows * ncols;
    case Matrix_TYPE::SYMMETRIC: {
      assert(nrows == ncols);
      return (nrows * (ncols + 1)) / 2;
    }
    case Matrix_TYPE::DIAGONAL: {
      return std::min(nrows, ncols);
      break;
    }
    case Matrix_TYPE::SCALAR_FULL:
    case Matrix_TYPE::SCALAR_DIAGONAL:
    default:
      return 1;
    }
  }
  std::size_t getSize(TYPE t, std::size_t nrows) {
    switch (t) {
    case Matrix_TYPE::ZERO:
      return 0;
    case Matrix_TYPE::FULL:
      return nrows * nrows;
    case Matrix_TYPE::SYMMETRIC: {
      return (nrows * (nrows + 1)) / 2;
    }
    case Matrix_TYPE::DIAGONAL: {
      return nrows;
      break;
    }
    case Matrix_TYPE::SCALAR_FULL:
    case Matrix_TYPE::SCALAR_DIAGONAL:
    default:
      return 1;
    }
  }

  TYPE type() const { return type_; }

  bool isSymmetric() const {
    switch (type_) {
    case Matrix_TYPE::ZERO:
    case Matrix_TYPE::SYMMETRIC:
      return true;
    case Matrix_TYPE::SCALAR_DIAGONAL:
    case Matrix_TYPE::DIAGONAL:
    case Matrix_TYPE::SCALAR_FULL:
      return ncols() == nrows();
    case Matrix_TYPE::FULL:
    default:
      return false;
    }
  }
  bool isZero() const {
      switch (type_) {
      case Matrix_TYPE::ZERO:
          return true;
      case Matrix_TYPE::SYMMETRIC:
      case Matrix_TYPE::SCALAR_DIAGONAL:
      case Matrix_TYPE::DIAGONAL:
      case Matrix_TYPE::SCALAR_FULL:
      case Matrix_TYPE::FULL:
      default:
          return false;
      }
  }


  bool isDiagonal() const {
    switch (type()) {
    case Matrix_TYPE::FULL:
    case Matrix_TYPE::SYMMETRIC:
    case Matrix_TYPE::SCALAR_FULL:
      return false;
    case Matrix_TYPE::DIAGONAL:
    case Matrix_TYPE::SCALAR_DIAGONAL:
    case Matrix_TYPE::ZERO:
    default:
      return true;
    }
  }

  /**
  Matrix of zeros with the shape of the provided Matrix
   @post nrows(ones(x))==x.nrows()   ncols(ones(x))=x.ncols()

  */

  M_Matrix(TYPE t, std::size_t nrows, std::size_t ncols, std::vector<T> data)
      : type_(t), _nrows(nrows), _ncols(ncols), _data(data) {}

  template <typename... Ts, template <typename...> class V>
  M_Matrix(std::size_t nrows, std::size_t ncols, TYPE t, V<Ts...> data)
      : type_(t), _nrows(nrows), _ncols(ncols), _data(data.size()) {
    for (std::size_t i = 0; i < data.size(); ++i)
      _data[i] = data[i];
    init_zero();
  }

  /*!
     Matrix Addition assignment.
     @returns a reference to itself
     @pre  same number of rows and columns
     @post all the values of the matrix are summed up by the corresponing
           values of the other matrix\n\n
           assert(nrows(*this)==nrows(other))&&
     assert(ncols(*this)==ncols(other))

     */


private:
    T const & zero(std::size_t i [[maybe_unused]], std::size_t j [[maybe_unused]])const
    {
        if constexpr (is_Matrix_v<T>)
            return zero_[i*ncols()+j];
        else
            return zero_;
    }
    void init_zero()
    {
        if constexpr (is_Matrix_v<T>)
        {
            if (type()==Matrix_TYPE::DIAGONAL)
            {
                zero_=std::vector<T>(nrows()*ncols());
                for (std::size_t i=0; i<nrows(); ++i)
                    for (std::size_t j=0; j<ncols(); ++j)
                        zero_[i*ncols()+j]=T((*this)[i].nrows(),(*this)[j].ncols(), Matrix_TYPE::ZERO);
          }
            else if (type()==Matrix_TYPE::SCALAR_DIAGONAL) {
              zero_=std::vector<T>(1,T((*this)[0].nrows(),(*this)[0].ncols(), Matrix_TYPE::ZERO));
            }
        }

    }

  TYPE type_ = Matrix_TYPE::ZERO;
  std::size_t _nrows; /**< number of rows */
  std::size_t _ncols; /**< number of columns */

  /** internal data
        @remarks can be a pointer or a vector
        @remarks now is a vector for debugging purposes
        */
  //	T                      *_data; /**< pointer to the data  */
  std::vector<T> _data;
  std::conditional_t<is_Matrix_v<T>,std::vector<T>,T> zero_=std::conditional_t<is_Matrix_v<T>,std::vector<T>,T>{};
};

template<typename T>
M_Matrix<M_Matrix<T>> to_Matrix(const std::vector<std::vector<M_Matrix<T>>>& x)
{
    std::size_t nrows=x.size();
    std::size_t ncols=0;
    std::size_t i_row_max_col=0;
    for (std::size_t i=0; i<nrows; ++i) if (x[i].size()>ncols)
        {
            ncols=x[i].size();
            i_row_max_col=i;
        }

    std::vector<std::size_t> n_subrows(nrows);
    for (std::size_t i=0; i<nrows; ++i) n_subrows[i]=x[i][0].nrows();
    std::vector<std::size_t> n_subcols(ncols);
    for (std::size_t j=0; j<ncols; ++j) n_subcols[j]=x[i_row_max_col][j].ncols();
    M_Matrix<M_Matrix<T>> out(nrows,ncols);
    for (std::size_t i=0; i<nrows; ++i)
    {

        for (std::size_t j=0; j<x[i].size(); ++j)
            out(i,j)=x[i][j];
        for (std::size_t j=x[i].size(); j<ncols; ++j)
            out(i,j)=M_Matrix<double>(n_subrows[i],n_subcols[j],Matrix_TYPE::ZERO);
    }
    return out;
}



template <typename T>
M_Matrix<T> full(const M_Matrix<M_Matrix<T>>& x, std::pair<std::size_t,std::size_t> i_range={0ul, std::string::npos}, std::pair<std::size_t,std::size_t> j_range={0ul, std::string::npos})
{
    if (i_range.second==std::string::npos) i_range.second=x.nrows();
    if (j_range.second==std::string::npos) j_range.second=x.ncols();

    std::size_t nrows=0;
    std::size_t ncols=0;
    if (!x.isDiagonal())
    {
        for (std::size_t i=i_range.first; i<i_range.second; ++i) nrows+=x(i,0).nrows();
        for (std::size_t i=j_range.first; i<j_range.second; ++i) ncols+=x(0,i).ncols();
    }
    else
    {
        for (std::size_t i=i_range.first; i<i_range.second; ++i) nrows+=x(i,i).nrows();
        for (std::size_t i=j_range.first; i<j_range.second; ++i) ncols+=x(i,i).ncols();
    }
    M_Matrix<T> out(nrows,ncols);
    std::size_t i_row=0;
    std::size_t i_col=0;
    for (std::size_t i=i_range.first; i<i_range.second; ++i)
    {

        i_col=0;
        for (std::size_t j=j_range.first; j<j_range.second; ++j)
        {
            assert(x(i,0).nrows()==x(i,j).nrows());
            for (std::size_t i2=0; i2<x(i,j).nrows(); ++i2)
                    for (std::size_t j2=0; j2<x(i,j).ncols(); ++j2)
                        out(i_row+i2,i_col+j2)=x(i,j)(i2,j2);
            i_col+=x(i,j).ncols();
        }
        i_row+=x(i,0).nrows();
     }
     return out;
}




template <bool output, template <bool, typename...> class test, typename T,
          class M, class ostream>
bool apply_test(const test<output, T> &t, const M_Matrix<T> &one,
                const M &other, ostream &os) {
  bool out = true;
  if constexpr (std::is_same_v<M, M_Matrix<T>>) {
    if ((one.ncols() != other.ncols()) || (one.nrows() != other.nrows())) {
      if constexpr (output) {
        os << "\n different size!!!"
           << " one.nrows=" << one.nrows() << " other.nrows=" << other.nrows();
        os << " one.ncols=" << one.ncols() << " other.ncols=" << other.ncols();
      }
      return false;
    }
  }
  for (std::size_t i = 0; i < one.nrows(); ++i) {
    for (std::size_t j = 0; j < one.ncols(); ++j) {
      bool mytest;
      if constexpr (std::is_same_v<M, M_Matrix<T>>)
        mytest = t.test(one(i, j), other(i, j), os);
      else
        mytest = t.test(one(i, j), other, os);

      if (!mytest) {
        if constexpr (output) {
          os << "  on i=" << i << " j=" << j;

          out = false;
        } else
          return false;
      }
    }
  }
  if (out)
    return true;
  else {
    if constexpr (output) {
      os << "\n one=\n" << one;
      os << "\n other=\n" << other;
    }
    return false;
  }
}

template <bool output, template <bool, typename> class test, typename T,
          class ostream>
bool apply_test(const test<output, T> &t, const M_Matrix<T> &one, ostream &os) {
  bool out = true;
  for (std::size_t i = 0; i < one.nrows(); ++i) {
    for (std::size_t j = 0; j < one.ncols(); ++j)
      if (!t.test(one(i, j), os)) {
        if constexpr (output) {
          os << "  on i=" << i << ", j=" << j;
          out = false;
        } else
          return false;
      }
  }
  if (out)
    return true;
  else {
    if constexpr (output) {
      os << "\n one=\n" << one;
    }
    return false;
  }
}

namespace Matrix_Unary_Functions {
template <typename T> double norm_1(const M_Matrix<T> &x);
template <typename T> double maxAbs(const M_Matrix<T> &x);
template <typename T> double minAbs(const M_Matrix<T> &x);


template <class E, class F, typename T>
T accumulate(const M_Matrix<E> &x, const F &f, const T &start) ;

inline double fullSum(const double &x) ;
template <class T> double fullSum(const M_Matrix<T> &x) ;

inline double colSum(const double &x) ;
template <class T> double colSum(const M_Matrix<T> &x) ;

inline double rowSum(const double &x) ;
template <class T> auto rowSum(const M_Matrix<T> &x) ;

inline double logProduct(double x);
template <class T> double logProduct(const M_Matrix<T> &x);





inline double maxAbs(double x) ;

template <typename T> double maxAbs(const M_Matrix<T> &x) ;

inline double minAbs(double x);

template <typename T> double minAbs(const M_Matrix<T> &x) ;


inline double max(const M_Matrix<double> &x) ;
inline double min(const M_Matrix<double> &x) ;

template <typename T> double min(const M_Matrix<T> &x) ;

template <typename T> double logDiagProduct(const M_Matrix<T> &x);

/**
       Product of the Diagonal of a Matrix

      */
template <typename T> double diagProduct(const M_Matrix<T> &x) ;

template <typename T> double Trace(const M_Matrix<T> &x) ;

inline double norm_1(const double x) ;
/**
    Maximum Value of the Sum of the absolute values in each row
   */
template <typename T> double norm_1(const M_Matrix<T> &x);

/**
    Maximum Value of the Sum of the absolute values in each row
   */
template <typename T> T norm_inf(const M_Matrix<T> &x) ;
template <typename T> std::size_t log2_norm(const M_Matrix<T> &x);


}

template <bool output, typename T> struct are_Equal<output, M_Matrix<T>> {
public:
  template <class ostream>
  bool test_sum(const M_Matrix<T> &one, const M_Matrix<T> &two,
                ostream &os = std::cerr) const {
    return apply_test(are_Equal<output, T>(absolute_, relative_), one, two, os);
  }

  template <class ostream>
  bool test(const M_Matrix<T> &one, const M_Matrix<T> &two,
            ostream &os = std::cerr) const {
    return test_sum(one, two, os);
  }
  template <class ostream>
  bool test_prod(const M_Matrix<T> &one, const M_Matrix<T> &two,
                 ostream &os = std::cerr) const {
    double N = std::max(Matrix_Unary_Functions::norm_1(one), 1.0);
    auto n = one.size();
    return apply_test(are_Equal<output, T>(std::sqrt(absolute_) * n * N,
                                           std::sqrt(relative_) * n),
                      one, two, os);
  }

  are_Equal(double absoluteError, double relativeError)
      : absolute_{absoluteError}, relative_{relativeError} {}
  are_Equal()
      : absolute_{std::numeric_limits<T>::epsilon() * 100},
        relative_{std::numeric_limits<T>::epsilon() * 100} {}

private:
  double absolute_;
  double relative_;
};

template <bool output, typename T> struct are_Equal<output, M_Matrix<M_Matrix<T>>> {
public:
    template <class ostream>
    bool test_sum(const M_Matrix<M_Matrix<T>> &one, const M_Matrix<M_Matrix<T>> &two,
                  ostream &os = std::cerr) const {
        return apply_test(are_Equal<output, M_Matrix<T>>(absolute_, relative_), one, two, os);
    }

    template <class ostream>
    bool test(const M_Matrix<M_Matrix<T>> &one, const M_Matrix<M_Matrix<T>> &two,
              ostream &os = std::cerr) const {
        return test_sum(one, two, os);
    }
    template <class ostream>
    bool test_prod(const M_Matrix<M_Matrix<T>> &one, const M_Matrix<M_Matrix<T>> &two,
                   ostream &os = std::cerr) const {
        double N = std::max(Matrix_Unary_Functions::norm_1(one), 1.0);
        auto n = one.size();
        return apply_test(are_Equal<output, M_Matrix<T>>(std::sqrt(absolute_) * n * N,
                                               std::sqrt(relative_) * n),
                          one, two, os);
    }

    are_Equal(double absoluteError, double relativeError)
        : absolute_{absoluteError}, relative_{relativeError} {}
    are_Equal()
        : absolute_{std::numeric_limits<T>::epsilon() * 100},
          relative_{std::numeric_limits<T>::epsilon() * 100} {}

private:
    double absolute_;
    double relative_;
};




template <bool output, typename T> struct are_non_negative<output, M_Matrix<T>> {
public:
  template <class ostream>
  bool test_sum(const M_Matrix<T> &one, ostream &os = std::cerr) const {
    return apply_test(are_non_negative<output, T>(absolute_), one, os);
  }
  template <class ostream>
  bool test_prod(const M_Matrix<T> &one, ostream &os = std::cerr) const {
    double N = norm_1(one);
    auto n = one.size();
    return apply_test(are_non_negative<output, T>(std::sqrt(absolute_) * n * N),
                      one, os);
  }

  are_non_negative(double absoluteError) : absolute_{absoluteError} {}
  are_non_negative() : absolute_{std::numeric_limits<T>::epsilon()} {}

private:
  double absolute_;
};

template <bool output, typename T> struct are_finite<output, M_Matrix<T>> {
public:
  template <class ostream>
  bool test(const M_Matrix<T> &one, ostream &os = std::cerr) const {
    return apply_test(are_finite<output, T>(), one, os);
  }
};

template <bool output, typename T> struct are_non_positive<output, M_Matrix<T>> {
public:
  template <class ostream>
  bool test_sum(const M_Matrix<T> &one, ostream &os = std::cerr) const {
    return apply_test(are_non_positive<output, T>(absolute_), one, os);
  }
  template <class ostream>
  bool test_prod(const M_Matrix<T> &one, ostream &os = std::cerr) const {
    double N = norm_1(one);
    auto n = one.size();
    return apply_test(are_non_positive<output, T>(std::sqrt(absolute_) * n * N),
                      one, os);
  }

  are_non_positive(double absoluteError) : absolute_{absoluteError} {}
  are_non_positive() : absolute_{std::numeric_limits<T>::epsilon()} {}

private:
  double absolute_;
};

template <bool output, typename T> struct are_not_less<output, M_Matrix<T>> {
public:
  template <class ostream>
  bool test(const M_Matrix<T> &one, ostream &os = std::cerr) const {
    //bool out = true;
    if (m_min_.size() == one.size())
      return apply_test(are_not_less<output, T>(absolute_error()), one, m_min_,
                        os);
    else
      return apply_test(are_not_less<output, T>(absolute_error()), one, min_,
                        os);
  }

  double min() const { return min_; }
  M_Matrix<double> const &Min() const { return m_min_; }
  double absolute_error() const { return absolute_; }
  are_not_less(bool missing, double min, double absoluteError)
      : missing_{missing}, min_{min}, absolute_{absoluteError}, m_min_{} {}
  are_not_less(bool missing, M_Matrix<double> &&min, double absoluteError)
      : missing_{missing}, min_{std::numeric_limits<double>::quiet_NaN()},
        absolute_{absoluteError}, m_min_{std::move(min)} {}

private:
  bool missing_;
  double min_;
  double absolute_ = std::numeric_limits<double>::epsilon();
  M_Matrix<double> m_min_;
};

template <bool output, typename T> struct are_not_more<output, M_Matrix<T>> {
public:
  template <class ostream>
  bool test(const M_Matrix<T> &one, ostream &os = std::cerr) const {
    if (m_max_.size() == one.size())
      return apply_test(are_not_more<output, T>(absolute_error()), one, m_max_,
                        os);
    else
      return apply_test(are_not_more<output, T>(absolute_error()), one, max_,
                        os);
  }

  double max() const { return max_; }
  M_Matrix<double> const &Max() const { return m_max_; }
  double absolute_error() const { return absolute_; }
  are_not_more(bool missing, double max, double absoluteError)
      : missing_{missing}, max_{max}, absolute_{absoluteError}, m_max_{} {}
  are_not_more(bool missing, M_Matrix<double> &&max, double absoluteError)
      : missing_{missing}, max_{std::numeric_limits<double>::quiet_NaN()},
        absolute_{absoluteError}, m_max_{std::move(max)} {}

private:
  bool missing_;
  double max_;
  double absolute_ = std::numeric_limits<double>::epsilon();
  M_Matrix<double> m_max_;
};

template <bool output, typename T> struct are_in_range<output, M_Matrix<T>> {
public:
  template <class ostream>
  bool test(const M_Matrix<T> &one, ostream &os = std::cerr) const {
    bool out = true;
    if (m_min_.size() == one.size()) {
      if (!apply_test(are_not_less<output, T>(absolute_error()), one, m_min_,
                      os))
        out = false;

    } else {
      if (!apply_test(are_not_less<output, T>(absolute_error()), one, min_, os))
        out = false;
    }
    if (m_max_.size() == one.size()) {
      if (!apply_test(are_not_more<output, T>(absolute_error()), one, m_max_,
                      os))
        out = false;

    } else {
      if (!apply_test(are_not_more<output, T>(absolute_error()), one, max_, os))
        out = false;
    }
    return out;
  }

  double min() const { return min_; }
  double max() const { return max_; }
  M_Matrix<double> const &Min() const { return m_min_; }
  M_Matrix<double> const &Max() const { return m_max_; }

  double absolute_error() const { return absolute_; }
  are_in_range(bool missing, double min, double max, double absoluteError)
      : missing_{missing}, min_{min}, max_{max}, absolute_{absoluteError},
        m_min_{}, m_max_{} {}
  are_in_range(bool missing, M_Matrix<double> &&min, M_Matrix<double> &&max,
               double absoluteError)
      : missing_{missing}, min_{std::numeric_limits<double>::quiet_NaN()},
        max_{std::numeric_limits<double>::quiet_NaN()},
        absolute_{absoluteError}, m_min_{std::move(min)}, m_max_{
                                                              std::move(max)} {}

private:
  bool missing_;
  double min_;
  double max_;
  double absolute_ = std::numeric_limits<double>::epsilon();
  M_Matrix<double> m_min_;
  M_Matrix<double> m_max_;
};

template <typename T> M_Matrix<T> operator-(const M_Matrix<T> &x) {
  M_Matrix<T> out(x.nrows(), x.ncols(), x.type());
  for (auto it = out.begin(); it != out.end(); ++it)
    *it = -x(it);
  return out;
}

template <typename T> M_Matrix<T> operator-(M_Matrix<T> &&x) {
  for (auto it = x.begin(); it != x.end(); ++it)
    *it = -*it;
  return std::move(x);
}

namespace Vector_Binary_Transformations {

template <typename T, typename K>
std::vector<T> operator-(const std::vector<T> &x, const std::vector<K> &y) {
  std::vector<T> out(x.size());
  for (std::size_t i = 0; i < x.size(); ++i)
    out[i] = x[i] - y[i];
  return out;
}

template <typename T>
std::vector<T> operator+(const std::vector<T> &x, const std::vector<T> &y) {
  std::vector<T> out(x.size());
  for (std::size_t i = 0; i < x.size(); ++i)
    out[i] = x[i] + y[i];
  return out;
}

template <typename T>
std::vector<T> elemMult(const std::vector<T> &x, const std::vector<T> &y) {
  std::vector<T> out(x.size());
  for (std::size_t i = 0; i < x.size(); ++i)
    out[i] = x[i] * y[i];
  return out;
}

template <typename T>
std::vector<T> elemDiv(const std::vector<T> &x, const std::vector<T> &y) {
  std::vector<T> out(x.size());
  for (std::size_t i = 0; i < x.size(); ++i)
    out[i] = x[i] / y[i];
  return out;
}

} // namespace Vector_Binary_Transformations
namespace Vector_Binary_Transformations {

inline std::vector<double> operator*(const std::vector<double> &x,
                                     const M_Matrix<double> &y) {
  std::vector<double> out(y.ncols(), 0);
  for (std::size_t i = 0; i < x.size(); ++i)
    for (std::size_t j = 0; j < y.ncols(); ++j)
      out[j] += x[i] * y(i, j);
  return out;
}

inline double xTSigmaX(const M_Matrix<double> &vector,
                       const M_Matrix<double> &matrix) {
  double sum = 0;
  if (matrix.isDiagonal()) {
    for (std::size_t i = 0; i < matrix.nrows(); ++i)
      sum += vector[i] * matrix(i, i) * vector[i];

  } else {
    for (std::size_t i = 0; i < matrix.nrows(); ++i) {
      sum += vector[i] * matrix(i, i) * vector[i];
      for (std::size_t j = i + 1; j < matrix.ncols(); ++j)
        sum += 2 * vector[i] * matrix(i, j) * vector[j];
    }
  }
  return sum;
}

template <class T>
double xTSigmaX(const M_Matrix<T> &vector, const M_Matrix<T> &matrix) {
  double sum = 0;
  if (matrix.isDiagonal()) {
    for (std::size_t i = 0; i < matrix.nrows(); ++i)
      sum += xTSigmaX(vector[i], matrix(i, i));
  } else {
    for (std::size_t i = 0; i < matrix.nrows(); ++i) {
      sum += xTSigmaX(vector[i], matrix(i, i));
      for (std::size_t j = i + 1; j < matrix.ncols(); ++j)
        sum += 2 * xTSigmaX(vector[i], matrix(i, j));
    }
  }
  return sum;
}

inline double xTSigmaX(const std::vector<double> &v,
                       const M_Matrix<double> &matrix) {
  double sum = 0;
  for (std::size_t i = 0; i < matrix.nrows(); ++i) {
    sum += v[i] * matrix(i, i) * v[i];
    for (std::size_t j = i + 1; j < matrix.ncols(); ++j)
      sum += 2 * v[i] * matrix(i, j) * v[j];
  }
  return sum;
}
inline double xTSigmaX(const std::vector<std::size_t> &v,
                       const M_Matrix<double> &matrix) {
    double sum = 0;
    for (std::size_t i = 0; i < matrix.nrows(); ++i) {
        sum += v[i] * matrix(i, i) * v[i];
        for (std::size_t j = i + 1; j < matrix.ncols(); ++j)
            sum += 2 * v[i] * matrix(i, j) * v[j];
    }
    return sum;
}
} // namespace Vector_Binary_Transformations

namespace Matrix_Unary_Transformations {

template <class T> M_Matrix<T> Symmetric_to_Full(const M_Matrix<T> &x) {
  assert(x.type() == Matrix_TYPE::SYMMETRIC);
  M_Matrix<T> out(x.nrows(), x.ncols(), Matrix_TYPE::FULL);
  for (std ::size_t i = 0; i < x.nrows(); ++i)
    for (std ::size_t j = 0; j < x.ncols(); ++j)
      out(i, j) = x(i, j);
  return out;
}

template <class T> myOptional_t<std::pair<M_Matrix<T>,double>> inv(const M_Matrix<T> &x);

template <class T>

myOptional_t<M_Matrix<T>> pinv(const M_Matrix<T> &x);

template <class T>
myOptional_t<M_Matrix<T>> chol(const M_Matrix<T> &x, const std::string &kind);

} // namespace Matrix_Unary_Transformations

using namespace Matrix_Unary_Transformations;

namespace lapack {

extern "C" void dpotrf_(char *UPLO, int *N, double *A, int *LDA, int *INFO);

inline M_Matrix<double> UT(const M_Matrix<double> &x) {
  M_Matrix<double> y(x.nrows(), x.ncols());
  for (std::size_t i = 0; i < x.nrows(); i++) {
    for (std::size_t j = 0; j < i; j++)
      y(i, j) = 0;
    for (std::size_t j = i; j < x.ncols(); j++)
      y(i, j) = x(i, j);
  }
  return y;
}

inline M_Matrix<double> LT(const M_Matrix<double> &x) {
  M_Matrix<double> y(x.nrows(), x.ncols());
  for (std::size_t i = 0; i < x.nrows(); i++) {
    for (std::size_t j = 0; j < i + 1; j++)
      y(i, j) = x(i, j);
    for (std::size_t j = i + 1; j < x.ncols(); j++)
      y(i, j) = 0;
  }
  return y;
}

} // namespace lapack

namespace Matrix_Binary_Transformations {
template <typename T>
M_Matrix<T> operator+(const M_Matrix<T> &x, const M_Matrix<T> &y);

template <typename T>
M_Matrix<T> operator-(const M_Matrix<T> &x, const M_Matrix<T> &y);

template <typename T> M_Matrix<T> operator-(M_Matrix<T> &&x, M_Matrix<T> &&y);

template <typename T, typename S>
M_Matrix<T> &operator+=(M_Matrix<T> &x, const M_Matrix<S> &y);

template <typename T, typename S>
M_Matrix<T> &operator-=(M_Matrix<T> &x, const M_Matrix<S> &y);

template <typename T, typename S>
auto operator*(const M_Matrix<T> &one, const M_Matrix<S> &other) -> M_Matrix<
    std::enable_if_t<is_Matrix_v<T> == is_Matrix_v<S>,
                     decltype(std::declval<T>() * std::declval<S>())>>;

template <typename T, typename S>
auto multTransp(const M_Matrix<T> &one, const M_Matrix<S> &other)
    -> M_Matrix<decltype(std::declval<T>() * std::declval<S>())>;

/**
     Transpose the first and multiply by the second
     @post transpMult(x,y)==Transpose(x)*y
     @remarks It is faster, since we save copying matrices
    */
template <typename T, typename S>
auto TranspMult(const M_Matrix<T> &one, const M_Matrix<S> &other)
    -> M_Matrix<decltype(std::declval<T>() * std::declval<S>())>;

template <typename T, typename S>
auto quadraticForm_BT_A_B(const M_Matrix<T> &A, const M_Matrix<S> &B)
    -> M_Matrix<decltype(std::declval<T>() * std::declval<S>())>;

template <typename T, typename S>
auto quadraticForm_B_A_BT(const M_Matrix<T> &A, const M_Matrix<S> &B)
    -> M_Matrix<decltype(std::declval<T>() * std::declval<S>())>;

template <typename T, typename S>
auto mean(const M_Matrix<T> &A, const M_Matrix<S> &B)
    -> M_Matrix<decltype(std::declval<T>() * std::declval<S>())>;





inline M_Matrix<double> Forward_Sustitution_Ly_b(const M_Matrix<double> &L,
                                                 const M_Matrix<double> &b);

template <typename T>
M_Matrix<T> Forward_Sustitution_Ly_b(const M_Matrix<T> &L,
                                     const M_Matrix<T> &b);

/**
     Scalar Multiplication assignment.
     @returns a reference to itself
     @post all the values of the matrix are multiplied by the value x
     */
template <typename E, typename T>
auto operator*=(M_Matrix<E> &itself, T x) -> M_Matrix<
    std::enable_if_t<!(is_Matrix_v<T> && !is_Matrix_v<E>),
                     decltype(std::declval<E>() * std::declval<T>())>> &;

/**
     Scalar Division assignment.
     @returns a reference to itself
     @post all the values of the matrix are divided by the value x
     */
template <typename E, typename T>
auto operator/=(M_Matrix<E> &itself, T x) -> M_Matrix<
    std::enable_if_t<!(is_Matrix_v<T> && !is_Matrix_v<E>),
                     decltype(std::declval<E>() * std::declval<T>())>> &;

/**
     Scalar Addition.
     @returns a copy of the matrix with its values summed by x
     */
template <typename T> M_Matrix<T> operator+(const M_Matrix<T> &x, T t);

/**
     Scalar Addition reverse order.
     */
template <typename T> M_Matrix<T> operator+(T t, const M_Matrix<T> &x);
/**
     Scalar Subtraction.
     @returns a copy of the matrix with its values substracted by x
     */
template <typename T> M_Matrix<T> operator-(const M_Matrix<T> &x, T t);

/**
     Scalar Subtraction reverse order.
     */
template <typename T> M_Matrix<T> operator-(T t, const M_Matrix<T> &x);
/**
     Scalar Multiplication.
     @returns a copy of the matrix with its values multiplied by the value x
     */
template <typename E, typename T>
auto operator*(const M_Matrix<E> &x, T t)
    -> std::enable_if_t<!(is_Matrix_v<T> && !is_Matrix_v<E>)&&std::is_same_v<
                            E, decltype(std::declval<E>() * std::declval<T>())>,
                        M_Matrix<E>>;
/**
     Scalar Multiplication reverse order.
     */

template <typename E, typename T>
auto operator*(T t, const M_Matrix<E> &x)
    -> std::enable_if_t<!(is_Matrix_v<T> && !is_Matrix_v<E>)&&std::is_same_v<
                            E, decltype(std::declval<E>() * std::declval<T>())>,
                        M_Matrix<E>>;

/**
     Scalar Division.
     @returns a copy of the matrix with its values divided by x
     @returns a matrix of real numbers
 */
template <typename T> M_Matrix<double> operator/(const M_Matrix<T> &x, T t);

/**
     Division by inhomogeneus types

     */

template <typename T, typename S>
M_Matrix<decltype(std::declval<T>() * std::declval<S>)>
operator/(const M_Matrix<T> &x, S t);

/**
     Scalar Division reverse order.
     */
template <typename T> M_Matrix<T> operator/(T t, const M_Matrix<T> &x);

} // namespace Matrix_Binary_Transformations

using namespace Matrix_Binary_Transformations;

namespace Matrix_Generators {
template <typename T> M_Matrix<T> ones(size_t nrows, size_t ncols) {
  return M_Matrix<T>(nrows, ncols, Matrix_TYPE::SCALAR_FULL, T(1));
}

template <typename T> M_Matrix<T> zeros(size_t nrows, size_t ncols) {
  return M_Matrix<T>(nrows, ncols, Matrix_TYPE::ZERO);
}

/**
  Identity Matrix of the specified size
 */
template <typename T> M_Matrix<T> eye(std::size_t n) {
  return M_Matrix<T>(n, n, Matrix_TYPE::SCALAR_DIAGONAL, T(1));
}

template <typename T> M_Matrix<T> eye(std::size_t nrows,std::size_t  ncols) {
    return M_Matrix<T>(nrows, ncols, Matrix_TYPE::SCALAR_DIAGONAL, T(1));
}


template <typename T> M_Matrix<T> Rand(const M_Matrix<T> &x) {
  std::normal_distribution<> normal;
  std::random_device rd;
  std::mt19937_64 sto(rd());
  auto out = M_Matrix<T>(x);
  for (std::size_t i = 0; i < out.size(); ++i)
    out[i] = normal(sto);
  return out;
}

template <class D>
M_Matrix<double> Rand(const M_Matrix<double> &x, D d, std::mt19937_64 &sto) {
  auto out = M_Matrix<double>(x.nrows(), x.ncols(), x.type());
  for (auto it = out.begin(); it != out.end(); ++it)
    *it = d(sto);
  return out;
}

template <typename T, class D>
M_Matrix<T> Rand(const M_Matrix<T> &x, D d, std::mt19937_64 &sto) {
  auto out = M_Matrix<T>(x.nrows(), x.ncols(), x.type());
  ;
  for (auto it = out.begin(); it != out.end(); ++it)
    *it = Rand(*it, d, sto);
  return out;
}
} // namespace Matrix_Generators

template <typename T>
std::ostream &operator<<(std::ostream &os, const M_Matrix<T> &x) {
  os << "[";
  for (std::size_t i = 0; i < x.nrows(); ++i) {
    for (std::size_t j = 0; j < x.ncols(); ++j)
      os << x(i, j) << " ";
    os << ";";
  }
  os << "]";
  return os;
}

template <typename T>
std::istream &operator>>(std::istream &is, M_Matrix<T> &x) {
  std::vector<T> o;
  std::size_t nrows = 0;
  char ch;
  while ((is >> ch) && (ch != '[')) {
  }
  if (ch != '[')
    return is;
  else
    while (ch != ']') {
      std::string s;
      while ((is.get(ch)) && ((ch != ']') && ch != ';')) {
        s.push_back(ch);
      }
      std::stringstream ss(s);
      T e;
      std::size_t i = o.size();
      while (ss >> e)
        o.push_back(e);
      if (o.size() > i)
        ++nrows;
    }
  std::size_t ncols = o.size() / nrows;
  x = M_Matrix<T>(nrows, ncols, o);
  return is;
}

namespace Matrix_Unary_Predicates {
template <typename T, class Predicate>
bool all(const M_Matrix<T> &x, const Predicate &p) {
  for (std::size_t i = 0; i < x.size(); ++i)
    if (!p(x[i]))
      return false;
  return true;
}

template <typename T, class Predicate>
bool any(const M_Matrix<T> &x, const Predicate &p) {
  for (std::size_t i = 0; i < x.size(); ++i)
    if (p(x[i]))
      return true;
  return false;
}

template <typename T, class Predicate>
std::set<std::size_t> find(const M_Matrix<T> &x, const Predicate &p) {
  std::set<std::size_t> out;
  for (std::size_t i = 0; i < x.size(); ++i)
    if (p(x[i]))
      out.insert(i);
  return out;
}

using std::isnan;
using std::isfinite;

template <typename T> bool isnan(const M_Matrix<T> &x) {
  return any(x, [](const T &e) { return isnan(e); });
}
template <typename T> bool isfinite(const M_Matrix<T> &x) {
    return all(x, [](const T &e) { return isfinite(e); });
}


} // namespace Matrix_Unary_Predicates

namespace Matrix_Unary_Size_Functions {

template <typename T>

std::size_t nrows(const M_Matrix<T> &x) {
  return x.nrows();
}

template <typename T> std::size_t ncols(const M_Matrix<T> &x) {
  return x.ncols();
}

template <typename T> std::size_t size(const M_Matrix<T> &x) {
  return x.size();
}

template <typename T> std::size_t rank_diag(const M_Matrix<T> &x, double tol=std::sqrt(std::numeric_limits<double>::epsilon())) {
  assert(x.isDiagonal());
  double e=Matrix_Unary_Functions::maxAbs(x)*tol;
  std::size_t r=0;
  for (std::size_t i=0; i<x.size(); ++i)
    if (std::abs(x[i])>e) ++r;
  return r;
}

} // namespace Matrix_Unary_Size_Functions
using namespace Matrix_Unary_Size_Functions;
namespace Matrix_Unary_Transformations {
template <class E, class F, typename T>
M_Matrix<T> accumulate_by_Rows(const M_Matrix<E> &me, const F &f,
                               const T &init);
template <class E, class F, typename T>
M_Matrix<T> accumulate_by_Cols(const M_Matrix<E> &me, const F &f,
                               const T &start);
template <typename T> M_Matrix<T> diag(const M_Matrix<T> &x);

template <typename T> M_Matrix<T> Transpose(const M_Matrix<T> &x);

template <typename T> M_Matrix<T> TransposeSum(const M_Matrix<T> &x);

template <typename T, std::size_t I0, std::size_t I1, std::size_t I2,
          std::size_t I3>
M_Matrix<M_Matrix<T>> Permute(const M_Matrix<M_Matrix<T>> &x);

} // namespace Matrix_Unary_Transformations
using namespace Matrix_Unary_Transformations;

namespace Matrix_Decompositions {

inline auto EigenSystem_full_real_eigenvalue_dgeev(const M_Matrix<double> &x) {
  typedef myOptional_t<
      std::tuple<M_Matrix<double>, M_Matrix<double>, M_Matrix<double>>>
      Op;

  using lapack::dgeev_;
  M_Matrix<double> VL(x.ncols(), x.nrows());
  M_Matrix<double> VR(x.ncols(), x.nrows());

  M_Matrix<double> L(x.nrows(), x.ncols(), Matrix_TYPE::DIAGONAL);
  M_Matrix<double> L_imag(x.nrows(), x.ncols(), Matrix_TYPE::DIAGONAL);

  M_Matrix<double> A = x;

  std::size_t N = x.ncols();

  char jobvl = 'N';
  char jobvr = 'V';

  int lda = N;
  int ldvl = N;
  int ldvr = N;
  auto work = std::make_unique<double[]>(N * N * 4);
  int lwork = N * N * 4;
  int n = N;
  int info;
  dgeev_(&jobvl, &jobvr, &n, &A(0, 0), &lda, &L(0, 0), &L_imag(0, 0), &VL(0, 0),
         &ldvl, &VR(0, 0), &ldvr, work.get(), &lwork, &info);

  auto invVR = inv(VR);
  if (info != 0)
    assert(false);

  are_zero<true, double> not_zero(
      std::sqrt(std::numeric_limits<double>::epsilon()));

  assert((apply_test(not_zero, L_imag, std::cerr)));

  return Op(std::make_tuple(
      VR, L, invVR.value().first)); // in reality VL, L, VR because of transposition
}

typedef std::tuple<M_Matrix<double>, M_Matrix<double>, M_Matrix<double>, M_Matrix<double>,M_Matrix<double>>
    eigensystem_type;

typedef std::tuple<M_Matrix<double>, M_Matrix<double>, M_Matrix<double>>
    SVD_type;

inline eigensystem_type sort(const eigensystem_type &e) {
  auto &[VR, L, VL,CL,CV] = e;
  std::vector<std::pair<double, std::size_t>> la(L.size());
  for (std::size_t i = 0; i < L.size(); ++i)
    la[i] = std::pair(L[i], i);
  std::sort(la.begin(), la.end());
  M_Matrix<double> Ls(L.nrows(), L.ncols(), L.type());
  M_Matrix<double> CLs(L.nrows(), L.ncols(), L.type());
  M_Matrix<double> CVs(L.nrows(), L.ncols(), L.type());

  M_Matrix<double> VRs(VR.nrows(), VR.ncols(), VR.type());
  M_Matrix<double> VLs(VL.nrows(), VL.ncols(), VL.type());
  for (std::size_t i = 0; i < la.size(); ++i) {
      std::size_t is = la[i].second;
    Ls[i] = L[is];
    CLs[i]=CL[is];
    CVs[i]=CV[is];

    for (std::size_t j = 0; j < la.size(); ++j) {
      VRs(j, i) = VR(j, is);
      VLs(i, j) = VL(is, j);
    }
  }
  return eigensystem_type(VRs, Ls, VLs,CLs,CVs);
}

inline M_Matrix<double> &Right_Eigenvector_Nelson_Normalization(M_Matrix<double> &X) {
  assert(X.nrows() == X.ncols());
  auto n = X.nrows();
  for (std::size_t k = 0; k < n; ++k) {
    double max = 0;
    for (std::size_t i = 0; i < n; ++i)
      if (std::abs(X(i, k)) >= std::abs(max)) {
        max = X(i, k);
      }
    for (std::size_t i = 0; i < n; ++i)
      X(i, k) = X(i, k) / max;
  }
  return X;
}
inline void Nelson_Normalization(M_Matrix<double> &VR, M_Matrix<double> &VL) {
  Right_Eigenvector_Nelson_Normalization(VR);
  for (std::size_t i = 0; i < VR.nrows(); ++i) {
    double sum = 0;
    for (std::size_t j = 0; j < VR.ncols(); ++j)
      sum += VL(i, j) * VR(j, i);
    for (std::size_t j = 0; j < VR.ncols(); ++j)
      VL(i, j) = VL(i, j) / sum;
  }
  //assert((are_Equal<true, M_Matrix<double>>().test_prod(
  //    VR * VL, Matrix_Generators::eye<double>(VR.ncols()), std::cerr)));
}

inline M_Matrix<double> &Left_Eigenvector_Nelson_Normalization(M_Matrix<double> &X) {
  assert(X.nrows() == X.ncols());
  auto n = X.nrows();
  for (std::size_t k = 0; k < n; ++k) {
    double max = 0;
    for (std::size_t i = 0; i < n; ++i)
      if (std::abs(X(k, i)) >= std::abs(max)) {
        max = X(k, i);
      }
    for (std::size_t i = 0; i < n; ++i)
      X(k, i) /= max;
  }
  return X;
}

inline auto EigenSystem_full_real_eigenvalue(const M_Matrix<double> &x, double min_inv_ConditionNumber=std::sqrt(std::numeric_limits<double>::epsilon())) {
  typedef myOptional_t<eigensystem_type> Op;

  using lapack::dgeevx_;

  /**
( 	char*  	BALANC,
                              );

*/

  /**
dgeevx()
subroutine dgeevx 	(
  character  	BALANC,
      character  	JOBVL,
      character  	JOBVR,
      character  	SENSE,
      integer  	N,
      double precision, dimension( lda, * )  	A,
      integer  	LDA,
      double precision, dimension( * )  	WR,
      double precision, dimension( * )  	WI,
      double precision, dimension( ldvl, * )  	VL,
      integer  	LDVL,
      double precision, dimension( ldvr, * )  	VR,
      integer  	LDVR,
      integer  	ILO,
      integer  	IHI,
      double precision, dimension( * )  	SCALE,
      double precision  	ABNRM,
      double precision, dimension( * )  	RCONDE,
      double precision, dimension( * )  	RCONDV,
      double precision, dimension( * )  	WORK,
      integer  	LWORK,
      integer, dimension( * )  	IWORK,
      integer  	INFO
  )

DGEEVX computes the eigenvalues and, optionally, the left and/or right
eigenvectors for GE matrices

Download DGEEVX + dependencies [TGZ] [ZIP] [TXT]

Purpose:

   DGEEVX computes for an N-by-N real nonsymmetric matrix A, the
   eigenvalues and, optionally, the left and/or right eigenvectors.

   Optionally also, it computes a balancing transformation to improve
   the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
   SCALE, and ABNRM), reciprocal condition numbers for the eigenvalues
   (RCONDE), and reciprocal condition numbers for the right
   eigenvectors (RCONDV).

   The right eigenvector v(j) of A satisfies
                    A * v(j) = lambda(j) * v(j)
   where lambda(j) is its eigenvalue.
   The left eigenvector u(j) of A satisfies
                 u(j)**H * A = lambda(j) * u(j)**H
   where u(j)**H denotes the conjugate-transpose of u(j).

   The computed eigenvectors are normalized to have Euclidean norm
   equal to 1 and largest component real.

   Balancing a matrix means permuting the rows and columns to make it
   more nearly upper triangular, and applying a diagonal similarity
   transformation D * A * D**(-1), where D is a diagonal matrix, to
   make its rows and columns closer in norm and the condition numbers
   of its eigenvalues and eigenvectors smaller.  The computed
   reciprocal condition numbers correspond to the balanced matrix.
   Permuting rows and columns will not change the condition numbers
   (in exact arithmetic) but diagonal scaling will.  For further
   explanation of balancing, see section 4.10.2 of the LAPACK
   Users' Guide.

Parameters

*/
  /**
  [in]	BALANC

            BALANC is CHARACTER*1
            Indicates how the input matrix should be diagonally scaled
            and/or permuted to improve the conditioning of its
            eigenvalues.
            = 'N': Do not diagonally scale or permute;
            = 'P': Perform permutations to make the matrix more nearly
                   upper triangular. Do not diagonally scale;
            = 'S': Diagonally scale the matrix, i.e. replace A by
                   D*A*D**(-1), where D is a diagonal matrix chosen
                   to make the rows and columns of A more equal in
                   norm. Do not permute;
            = 'B': Both diagonally scale and permute A.

            Computed reciprocal condition numbers will be for the matrix
            after balancing and/or permuting. Permuting does not change
            condition numbers (in exact arithmetic), but balancing does.
*/
  char BALANC = 'B';

  /**
  [in]	JOBVL

            JOBVL is CHARACTER*1
            = 'N': left eigenvectors of A are not computed;
            = 'V': left eigenvectors of A are computed.
            If SENSE = 'E' or 'B', JOBVL must = 'V'.
*/
  char JOBVL = 'V';
  /**
  [in]	JOBVR

            JOBVR is CHARACTER*1
            = 'N': right eigenvectors of A are not computed;
            = 'V': right eigenvectors of A are computed.
            If SENSE = 'E' or 'B', JOBVR must = 'V'.
 */

  char JOBVR = 'V';

  /**

  [in]	SENSE

            SENSE is CHARACTER*1
            Determines which reciprocal condition numbers are computed.
            = 'N': None are computed;
            = 'E': Computed for eigenvalues only;
            = 'V': Computed for right eigenvectors only;
            = 'B': Computed for eigenvalues and right eigenvectors.

            If SENSE = 'E' or 'B', both left and right eigenvectors
            must also be computed (JOBVL = 'V' and JOBVR = 'V').

  */
  char SENSE = 'B';

  /**

  [in]	N

            N is INTEGER
            The order of the matrix A. N >= 0.
  */
  int N = x.nrows();

  /**

  [in,out]	A

            A is DOUBLE PRECISION array, dimension (LDA,N)
            On entry, the N-by-N matrix A.
            On exit, A has been overwritten.  If JOBVL = 'V' or
            JOBVR = 'V', A contains the real Schur form of the balanced
            version of the input matrix A.

  */
  M_Matrix<double> /* precision; dimension( lda; * ) */ A(x);

  /**
  [in]	LDA

            LDA is INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

  */
  int LDA = N;

  /**
  [out]	WR

            WR is DOUBLE PRECISION array, dimension (N)

  */
  M_Matrix<double> /* precision; dimension( * ) */ WR(N, N,
                                                      Matrix_TYPE::DIAGONAL);

  /**
  [out]	WI

            WI is DOUBLE PRECISION array, dimension (N)
            WR and WI contain the real and imaginary parts,
            respectively, of the computed eigenvalues.  Complex
            conjugate pairs of eigenvalues will appear consecutively
            with the eigenvalue having the positive imaginary part
            first.

  */
  M_Matrix<double> /* precision; dimension( * ) */ WI(N, N,
                                                      Matrix_TYPE::DIAGONAL);

  /**
  [out]	VL

            VL is DOUBLE PRECISION array, dimension (LDVL,N)
            If JOBVL = 'V', the left eigenvectors u(j) are stored one
            after another in the columns of VL, in the same order
            as their eigenvalues.
            If JOBVL = 'N', VL is not referenced.
            If the j-th eigenvalue is real, then u(j) = VL(:,j),
            the j-th column of VL.
            If the j-th and (j+1)-st eigenvalues form a complex
            conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
            u(j+1) = VL(:,j) - i*VL(:,j+1).

  */
  M_Matrix<double> /* precision; dimension( ldvl; * ) */ VL_lapack(N, N);

  /**
  [in]	LDVL

            LDVL is INTEGER
            The leading dimension of the array VL.  LDVL >= 1; if
            JOBVL = 'V', LDVL >= N.

  */
  int LDVL = N;

  /**
  [out]	VR

            VR is DOUBLE PRECISION array, dimension (LDVR,N)
            If JOBVR = 'V', the right eigenvectors v(j) are stored one
            after another in the columns of VR, in the same order
            as their eigenvalues.
            If JOBVR = 'N', VR is not referenced.
            If the j-th eigenvalue is real, then v(j) = VR(:,j),
            the j-th column of VR.
            If the j-th and (j+1)-st eigenvalues form a complex
            conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
            v(j+1) = VR(:,j) - i*VR(:,j+1).

  */
  M_Matrix<double> /* precision; dimension( ldvr; * ) */ VR_lapack(N, N);

  /**
  [in]	LDVR

            LDVR is INTEGER
            The leading dimension of the array VR.  LDVR >= 1, and if
            JOBVR = 'V', LDVR >= N.
  */
  int LDVR = N;

  /**

  [out]	ILO

            ILO is INTEGER

  */
  int ILO;

  /**
  [out]	IHI

            IHI is INTEGER
            ILO and IHI are integer values determined when A was
            balanced.  The balanced A(i,j) = 0 if I > J and
            J = 1,...,ILO-1 or I = IHI+1,...,N.

  */
  int IHI;

  /**
  [out]	SCALE

            SCALE is DOUBLE PRECISION array, dimension (N)
            Details of the permutations and scaling factors applied
            when balancing A.  If P(j) is the index of the row and column
            interchanged with row and column j, and D(j) is the scaling
            factor applied to row and column j, then
            SCALE(J) = P(J),    for J = 1,...,ILO-1
                     = D(J),    for J = ILO,...,IHI
                     = P(J)     for J = IHI+1,...,N.
            The order in which the interchanges are made is N to IHI+1,
            then 1 to ILO-1.
  */
  M_Matrix<double> /* precision; dimension( * ) */ SCALE(N, N,
                                                         Matrix_TYPE::DIAGONAL);

  /**

  [out]	ABNRM

            ABNRM is DOUBLE PRECISION
            The one-norm of the balanced matrix (the maximum
            of the sum of absolute values of elements of any column).

  */
  double ABNRM;

  /**
  [out]	RCONDE

            RCONDE is DOUBLE PRECISION array, dimension (N)
            RCONDE(j) is the reciprocal condition number of the j-th
            eigenvalue.

  */
  M_Matrix<double> /* precision; dimension( * ) */ RCONDE(
      N, N, Matrix_TYPE::DIAGONAL);

  /**
  [out]	RCONDV

            RCONDV is DOUBLE PRECISION array, dimension (N)
            RCONDV(j) is the reciprocal condition number of the j-th
            right eigenvector.
  */
  M_Matrix<double> /* precision; dimension( * ) */ RCONDV(
      N, N, Matrix_TYPE::DIAGONAL);

  /**

  [out]	WORK

            WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

  */
  M_Matrix<double> /* precision; dimension( * ) */ WORK(1, 1, 0.0);

  /**
  [in]	LWORK

            LWORK is INTEGER
            The dimension of the array WORK.   If SENSE = 'N' or 'E',
            LWORK >= max(1,2*N), and if JOBVL = 'V' or JOBVR = 'V',
            LWORK >= 3*N.  If SENSE = 'V' or 'B', LWORK >= N*(N+6).
            For good performance, LWORK must generally be larger.

            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the WORK array, returns
            this value as the first entry of the WORK array, and no error
            message related to LWORK is issued by XERBLA.

  */
  int LWORK = -1;

  /**
  [out]	IWORK

            IWORK is INTEGER array, dimension (2*N-2)
            If SENSE = 'N' or 'E', not referenced.

  */
  std::vector<int> /*dimension( * ) */ IWORK(std::max(2*N - 2,0));

  /**
  [out]	INFO

            INFO is INTEGER
            = 0:  successful exit
            < 0:  if INFO = -i, the i-th argument had an illegal value.
            > 0:  if INFO = i, the QR algorithm failed to compute all the
                  eigenvalues, and no eigenvectors or condition numbers
                  have been computed; elements 1:ILO-1 and i+1:N of WR
                  and WI contain eigenvalues which have converged.

       */

  int INFO;
  dgeevx_(&BALANC, &JOBVL, &JOBVR, &SENSE, &N, &A(0, 0), &LDA, &WR(0, 0),
          &WI(0, 0), &VL_lapack[0], &LDVL, &VR_lapack[0], &LDVR, &ILO, &IHI,
          &SCALE[0], &ABNRM, &RCONDE[0], &RCONDV[0], &WORK[0], &LWORK,
          &IWORK[0], &INFO);

  M_Matrix<double> WORK_OPT(1, WORK[0]);
  LWORK = WORK[0];
  dgeevx_(&BALANC, &JOBVL, &JOBVR, &SENSE, &N, &A(0, 0), &LDA, &WR(0, 0),
          &WI(0, 0), &VL_lapack[0], &LDVL, &VR_lapack[0], &LDVR, &ILO, &IHI,
          &SCALE[0], &ABNRM, &RCONDE[0], &RCONDV[0], &WORK_OPT[0], &LWORK,
          &IWORK[0], &INFO);

  if (INFO != 0) {
    if (INFO > 0)
      return Op(
          false,
          std::string("the QR algorithm failed to compute all the  "
                      "eigenvalues, and no eigenvectors or condition numbers "
                      "have been computed; elements 1:ILO-1 and i+1:N of WR "
                      "and WI contain eigenvalues which have converged ILO=") +
              std::to_string(ILO) + " i=" + std::to_string(INFO));
    else {
      auto args = std::make_tuple(&BALANC, &JOBVL, &JOBVR, &SENSE, &N, &A(0, 0),
                                  &LDA, &WR(0, 0), &WI(0, 0), &VL_lapack[0],
                                  &LDVL, &VR_lapack[0], &LDVR, &ILO, &IHI,
                                  &SCALE[0], &ABNRM, &RCONDE[0], &RCONDV[0],
                                  &WORK_OPT[0], &LWORK, &IWORK[0], &INFO);
      std::string argumentsNames[] = {
          "BALANC", "JOBVL", "JOBVR", "SENSE", "N",     "A",
          "LDA",    "WR",    "WI",    "VL",    "LDVL",  "VR",
          "LDVR",   "ILO",   "IHI",   "SCALE", "ABNRM", "RCONDE",
          "RCONDV", "WORK",  "LWORK", "IWORK", "INFO"};

      std::stringstream ss;
      std::size_t i = -INFO;
      io::write_tuple_i(ss, args, i);
      return Op(false, "the" + std::to_string(i) + "-th argument " +
                           argumentsNames[i] +
                           " had the illegal value =" + ss.str());
    }

  } else {

      double invCondiNumber_la=Matrix_Unary_Functions::min(RCONDE)/ABNRM;
      double invCondiNumber_v=Matrix_Unary_Functions::min(RCONDV)/ABNRM;
      if ((invCondiNumber_la<min_inv_ConditionNumber )|| (invCondiNumber_v<min_inv_ConditionNumber)){
          return Op(false,"Eigendecomposition condition value is too big :"+ToString(1.0/invCondiNumber_la));
      }
      else {

    auto &VL_cpp = VR_lapack;
    auto VR_cpp = Transpose(VL_lapack);
    Nelson_Normalization(VR_cpp, VL_cpp);
                                  //assert((are_Equal<true, M_Matrix<double>>().test_prod(VR_cpp * WR * VL_cpp,
    //                                                      x, std::cerr)));
    return Op(sort(std::make_tuple(
        VR_cpp, WR,
        VL_cpp,ABNRM*RCONDE.apply([](auto x){return 1.0/x;}),ABNRM*RCONDV.apply([](auto x){return 1.0/x;}))));
      }
  }
}





inline auto EigenSystem_symm_real_eigenvalue(const M_Matrix<double> &x, std::string kind="lower") {
  typedef myOptional_t<eigensystem_type> Op;

  assert(x.isSymmetric());
  using lapack::dsyevx_;

  /** DSYEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices

Download DSYEVX + dependencies [TGZ] [ZIP] [TXT]

Purpose:

     DSYEVX computes selected eigenvalues and, optionally, eigenvectors
     of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
     selected by specifying either a range of values or a range of indices
     for the desired eigenvalues.
*/
  /**
Parameters
    [in]	JOBZ

              JOBZ is CHARACTER*1
              = 'N':  Compute eigenvalues only;
              = 'V':  Compute eigenvalues and eigenvectors.
*/

  char JOBZ = 'V';

  /**
    [in]	RANGE

              RANGE is CHARACTER*1
              = 'A': all eigenvalues will be found.
              = 'V': all eigenvalues in the half-open interval (VL,VU]
                     will be found.
              = 'I': the IL-th through IU-th eigenvalues will be found.


*/

  char RANGE = 'A';

  /**
    [in]	UPLO

              UPLO is CHARACTER*1
              = 'U':  Upper triangle of A is stored;
              = 'L':  Lower triangle of A is stored.


*/
  char UPLO;
  if (kind == "lower")
    UPLO = 'U';
  else
    UPLO = 'L';

  /*
    [in]	N

              N is INTEGER
              The order of the matrix A.  N >= 0.

*/

  int N = x.nrows();

/**



    [in,out]	A

              A is DOUBLE PRECISION array, dimension (LDA, N)
              On entry, the symmetric matrix A.  If UPLO = 'U', the
              leading N-by-N upper triangular part of A contains the
              upper triangular part of the matrix A.  If UPLO = 'L',
              the leading N-by-N lower triangular part of A contains
              the lower triangular part of the matrix A.
              On exit, the lower triangle (if UPLO='L') or the upper
              triangle (if UPLO='U') of A, including the diagonal, is
              destroyed.

*/

  auto A = Symmetric_to_Full(x);
 /*   M_Matrix<double> A(x.nrows(), x.ncols());
  if (kind != "lower") {
    for (std::size_t i = 0; i < x.nrows(); ++i)
      for (std::size_t j = i; j < x.ncols(); ++j)
        A(i, j) = x(i, j);
  } else {
    for (std::size_t i = 0; i < x.nrows(); ++i)
      for (std::size_t j = 0; j < i + 1; ++j)
        A(i, j) = x(i, j);
  }
*/


/**

    [in]	LDA

              LDA is INTEGER
              The leading dimension of the array A.  LDA >= max(1,N).

*/

  int LDA = N;


/**


    [in]	VL

              VL is DOUBLE PRECISION
              If RANGE='V', the lower bound of the interval to
              be searched for eigenvalues. VL < VU.
              Not referenced if RANGE = 'A' or 'I'.
*/

  double VL;


  /**


    [in]	VU

              VU is DOUBLE PRECISION
              If RANGE='V', the upper bound of the interval to
              be searched for eigenvalues. VL < VU.
              Not referenced if RANGE = 'A' or 'I'.
*/

  double VU;

  /**




    [in]	IL

              IL is INTEGER
              If RANGE='I', the index of the
              smallest eigenvalue to be returned.
              1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
              Not referenced if RANGE = 'A' or 'V'.


*/

  int IL;


  /**
    [in]	IU

              IU is INTEGER
              If RANGE='I', the index of the
              largest eigenvalue to be returned.
              1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
              Not referenced if RANGE = 'A' or 'V'.
*/
  int IU;


      /**

    [in]	ABSTOL

              ABSTOL is DOUBLE PRECISION
              The absolute error tolerance for the eigenvalues.
              An approximate eigenvalue is accepted as converged
              when it is determined to lie in an interval [a,b]
              of width less than or equal to

                      ABSTOL + EPS *   max( |a|,|b| ) ,

              where EPS is the machine precision.  If ABSTOL is less than
              or equal to zero, then  EPS*|T|  will be used in its place,
              where |T| is the 1-norm of the tridiagonal matrix obtained
              by reducing A to tridiagonal form.

              Eigenvalues will be computed most accurately when ABSTOL is
              set to twice the underflow threshold 2*DLAMCH('S'), not zero.
              If this routine returns with INFO>0, indicating that some
              eigenvectors did not converge, try setting ABSTOL to
              2*DLAMCH('S').

              See "Computing Small Singular Values of Bidiagonal Matrices
              with Guaranteed High Relative Accuracy," by Demmel and
              Kahan, LAPACK Working Note #3.
*/

  double ABSTOL = std ::numeric_limits<double>::epsilon()*100;

  /**
    [out]	M

              M is INTEGER
              The total number of eigenvalues found.  0 <= M <= N.
              If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
*/

  int M;

/**
    [out]	W

              W is DOUBLE PRECISION array, dimension (N)
              On normal exit, the first M elements contain the selected
              eigenvalues in ascending order.
*/

  M_Matrix<double> W(x.nrows(), x.ncols(), Matrix_TYPE::DIAGONAL);


  /**
    [out]	Z

              Z is DOUBLE PRECISION array, dimension (LDZ, max(1,M))
              If JOBZ = 'V', then if INFO = 0, the first M columns of Z
              contain the orthonormal eigenvectors of the matrix A
              corresponding to the selected eigenvalues, with the i-th
              column of Z holding the eigenvector associated with W(i).
              If an eigenvector fails to converge, then that column of Z
              contains the latest approximation to the eigenvector, and the
              index of the eigenvector is returned in IFAIL.
              If JOBZ = 'N', then Z is not referenced.
              Note: the user must ensure that at least max(1,M) columns are
              supplied in the array Z; if RANGE = 'V', the exact value of M
              is not known in advance and an upper bound must be used.



    [in]	LDZ

              LDZ is INTEGER
              The leading dimension of the array Z.  LDZ >= 1, and if
              JOBZ = 'V', LDZ >= max(1,N).



*/

  int LDZ = N;
  M_Matrix<double> Z(x.nrows(), x.nrows());
  /*

    [out]	WORK

              WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

*/

  std ::vector<double> WORK(1);

  /**
    [in]	LWORK

              LWORK is INTEGER
              The length of the array WORK.  LWORK >= 1, when N <= 1;
              otherwise 8*N.
              For optimal efficiency, LWORK >= (NB+3)*N,
              where NB is the max of the blocksize for DSYTRD and DORMTR
              returned by ILAENV.

              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA.
*/

  int LWORK = -1;

  /**
    [out]	IWORK

              IWORK is INTEGER array, dimension (5*N)

*/

  std ::vector<int> IWORK(5 * x.nrows());


  /**
    [out]	IFAIL

              IFAIL is INTEGER array, dimension (N)
              If JOBZ = 'V', then if INFO = 0, the first M elements of
              IFAIL are zero.  If INFO > 0, then IFAIL contains the
              indices of the eigenvectors that failed to converge.
              If JOBZ = 'N', then IFAIL is not referenced.
*/
  std::vector<int> IFAIL(x.nrows());

  /**
    [out]	INFO

              INFO is INTEGER
              = 0:  successful exit
              < 0:  if INFO = -i, the i-th argument had an illegal value
              > 0:  if INFO = i, then i eigenvectors failed to converge.
                    Their indices are stored in array IFAIL.
*/

  int INFO;


  dsyevx_(&JOBZ, &RANGE, &UPLO, &N, &A[0], &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, &W[0], &Z[0], &LDZ, &WORK[0], &LWORK, &IWORK[0], &IFAIL[0],&INFO );

  LWORK = WORK[0];
  WORK.resize(LWORK);
  dsyevx_(&JOBZ, &RANGE, &UPLO, &N, &A[0], &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, &W[0], &Z[0], &LDZ, &WORK[0], &LWORK, &IWORK[0], &IFAIL[0], &INFO);




  if (INFO != 0) {
    if (INFO > 0)
      return Op(
          false,

std::to_string(INFO)+
          std::string(" eigenvectors failed to converge. Their indices are: ") +
              io::ToString(IFAIL));
    else {
      auto args = std::make_tuple(&JOBZ, &RANGE, &UPLO, &N, &A[0], &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, &W, &Z[0], &LDZ, &WORK[0], &LWORK, &IWORK[0], &IFAIL[0],&INFO);
      std::string argumentsNames[] = {
          "JOBZ", "RANGE", "UPLO", "N", "A[0]", "LDA", "VL", "VU", "IL", "IU", "ABSTOL", "M", "W", "Z[0]", "LDZ", "WORK[0]", "LWORK",  "IWORK[0]","IFAIL[0]", "INFO"    };

      std::stringstream ss;
      std::size_t i = -INFO;
      io::write_tuple_i(ss, args, i);
      return Op(false, "the" + std::to_string(i) + "-th argument " +
                           argumentsNames[i] +
                           " had the illegal value =" + ss.str());
    }

  } else {
      /**
 * extern "C" void ddisna_(char *JOB, int * M, int *N,double * D, double *SEP,int * INFO );
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            INFO, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), SEP( * )
*     ..
*
*  Purpose
*  =======
*
*  DDISNA computes the reciprocal condition numbers for the eigenvectors
*  of a real symmetric or complex Hermitian matrix or for the left or
*  right singular vectors of a general m-by-n matrix. The reciprocal
*  condition number is the 'gap' between the corresponding eigenvalue or
*  singular value and the nearest other one.
*
*  The bound on the error, measured by angle in radians, in the I-th
*  computed vector is given by
*
*         DLAMCH( 'E' ) * ( ANORM / SEP( I ) )
*
*  where ANORM = 2-norm(A) = max( abs( D(j) ) ).  SEP(I) is not allowed
*  to be smaller than DLAMCH( 'E' )*ANORM in order to limit the size of
*  the error bound.
*
*  DDISNA may also be used to compute error bounds for eigenvectors of
*  the generalized symmetric definite eigenproblem.
*/

      /**
*  JOB     (input) CHARACTER*1
*          Specifies for which problem the reciprocal condition numbers
*          should be computed:
*          = 'E':  the eigenvectors of a symmetric/Hermitian matrix;
*          = 'L':  the left singular vectors of a general matrix;
*          = 'R':  the right singular vectors of a general matrix.
*/
      char JOB='E';
      /**  M       (input) INTEGER
*          The number of rows of the matrix. M >= 0.
*
*/

      /**
*  N       (input) INTEGER
*          If JOB = 'L' or 'R', the number of columns of the matrix,
*          in which case N >= 0. Ignored if JOB = 'E'.
*/

      /**  D       (input) DOUBLE PRECISION array, dimension (M) if JOB = 'E'
*                              dimension (min(M,N)) if JOB = 'L' or 'R'
*          The eigenvalues (if JOB = 'E') or singular values (if JOB =
*          'L' or 'R') of the matrix, in either increasing or decreasing
*          order. If singular values, they must be non-negative.
*
*/

      /**  SEP     (output) DOUBLE PRECISION array, dimension (M) if JOB = 'E'
*                               dimension (min(M,N)) if JOB = 'L' or 'R'
*          The reciprocal condition numbers of the vectors.
*/
          M_Matrix<double> RCOND(x.nrows(), x.ncols(), Matrix_TYPE::DIAGONAL);

      /**  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*/
          lapack::ddisna_(&JOB,&M,&N,&W[0],&RCOND[0] ,&INFO );
    auto VL_cpp = Z;
    auto VR_cpp = Transpose(Z);
    Nelson_Normalization(VR_cpp, VL_cpp);

    assert((are_Equal<true, M_Matrix<double>>().test_prod(Matrix_Generators
                                                          ::eye<double>(x.nrows()),
                                                          VL_cpp * VR_cpp, std::cerr)));
    assert((are_Equal<true, M_Matrix<double>>().test_prod(VR_cpp * W * VL_cpp,
                                                          x, std::cerr)));
    assert((are_Equal<true, M_Matrix<double>>().test_prod(VR_cpp * W * VL_cpp,
                                                          x, std::cerr)));
    return Op(std::make_tuple(
        VR_cpp, W,
        VL_cpp, Matrix_Generators::eye<double>(x.ncols()),RCOND));
  }
}

inline auto Singular_Value_decomposition(const M_Matrix<double> &A) {

  /***
   *
   * dgesdd()
subroutine dgesdd 	( 	character  	JOBZ,
      integer  	M,
      integer  	N,
      double precision, dimension( lda, * )  	A,
      integer  	LDA,
      double precision, dimension( * )  	S,
      double precision, dimension( ldu, * )  	U,
      integer  	LDU,
      double precision, dimension( ldvt, * )  	VT,
      integer  	LDVT,
      double precision, dimension( * )  	WORK,
      integer  	LWORK,
      integer, dimension( * )  	IWORK,
      integer  	INFO
  )

DGESDD

Download DGESDD + dependencies [TGZ] [ZIP] [TXT]

Purpose:

   DGESDD computes the singular value decomposition (SVD) of a real
   M-by-N matrix A, optionally computing the left and right singular
   vectors.  If singular vectors are desired, it uses a
   divide-and-conquer algorithm.

   The SVD is written

        A = U * SIGMA * transpose(V)

   where SIGMA is an M-by-N matrix which is zero except for its
   min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
   V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
   are the singular values of A; they are real and non-negative, and
   are returned in descending order.  The first min(m,n) columns of
   U and V are the left and right singular vectors of A.

   Note that the routine returns VT = V**T, not V.

   The divide and conquer algorithm makes very mild assumptions about
   floating point arithmetic. It will work on machines with a guard
   digit in add/subtract, or on those binary machines without guard
   digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
   Cray-2. It could conceivably fail on hexadecimal or decimal machines
   without guard digits, but we know of none.

Parameters

  [in]	JOBZ

            JOBZ is CHARACTER*1
            Specifies options for computing all or part of the matrix U:
            = 'A':  all M columns of U and all N rows of V**T are
                    returned in the arrays U and VT;
            = 'S':  the first min(M,N) columns of U and the first
                    min(M,N) rows of V**T are returned in the arrays U
                    and VT;
            = 'O':  If M >= N, the first N columns of U are overwritten
                    on the array A and all rows of V**T are returned in
                    the array VT;
                    otherwise, all columns of U are returned in the
                    array U and the first M rows of V**T are overwritten
                    in the array A;
            = 'N':  no columns of U or rows of V**T are computed.
*/

  char JOBZ = 'A';

  /*
  [in]	M

            M is INTEGER
            The number of rows of the input matrix A.  M >= 0.

   */

  int M = A.nrows();

  /*


  [in]	N

            N is INTEGER
            The number of columns of the input matrix A.  N >= 0.

            */

  int N = A.ncols();

  /*

  [in,out]	A

            A is DOUBLE PRECISION array, dimension (LDA,N)
            On entry, the M-by-N matrix A.
            On exit,
            if JOBZ = 'O',  A is overwritten with the first N columns
                            of U (the left singular vectors, stored
                            columnwise) if M >= N;
                            A is overwritten with the first M rows
                            of V**T (the right singular vectors, stored
                            rowwise) otherwise.
            if JOBZ .ne. 'O', the contents of A are destroyed.

 */

  M_Matrix<double> Ao;

  if (A.type() == Matrix_TYPE::SYMMETRIC)
    Ao = Symmetric_to_Full(A);
  else {
    Ao = Transpose(A);
  }

  /*
  [in]	LDA

            LDA is INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).


            */

  int LDA = M;

  /*
  [out]	S

            S is DOUBLE PRECISION array, dimension (min(M,N))
            The singular values of A, sorted so that S(i) >= S(i+1).
*/

  std::size_t minMN = std::min(M, N);

  M_Matrix<double> S(M, N, Matrix_TYPE::DIAGONAL);

  /*


  [out]	U

            U is DOUBLE PRECISION array, dimension (LDU,UCOL)
            UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N;
            UCOL = min(M,N) if JOBZ = 'S'.
            If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M
            orthogonal matrix U;
            if JOBZ = 'S', U contains the first min(M,N) columns of U
            (the left singular vectors, stored columnwise);
            if JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced.

  [in]	LDU

            LDU is INTEGER
            The leading dimension of the array U.  LDU >= 1; if
            JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M.

            */

  int LDU = M;

  M_Matrix<double> U(LDU, M);

  /*

  [out]	VT

            VT is DOUBLE PRECISION array, dimension (LDVT,N)
            If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the
            N-by-N orthogonal matrix V**T;
            if JOBZ = 'S', VT contains the first min(M,N) rows of
            V**T (the right singular vectors, stored rowwise);
            if JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced.

  [in]	LDVT

            LDVT is INTEGER
            The leading dimension of the array VT.  LDVT >= 1;
            if JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N;
            if JOBZ = 'S', LDVT >= min(M,N).
*/

  int LDVT = N;

  M_Matrix<double> VT(N, LDVT);

  /*


  [out]	WORK

            WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
*/

  std::vector<double> WORK(1);

  /*


  [in]	LWORK

            LWORK is INTEGER
            The dimension of the array WORK. LWORK >= 1.
            If LWORK = -1, a workspace query is assumed.  The optimal
            size for the WORK array is calculated and stored in WORK(1),
            and no other work except argument checking is performed.

            Let mx = max(M,N) and mn = min(M,N).
            If JOBZ = 'N', LWORK >= 3*mn + max( mx, 7*mn ).
            If JOBZ = 'O', LWORK >= 3*mn + max( mx, 5*mn*mn + 4*mn ).
            If JOBZ = 'S', LWORK >= 4*mn*mn + 7*mn.
            If JOBZ = 'A', LWORK >= 4*mn*mn + 6*mn + mx.
            These are not tight minimums in all cases; see comments inside code.
            For good performance, LWORK should generally be larger;
            a query is recommended.

            */

  int LWORK = -1;

  /*

  [out]	IWORK

            IWORK is INTEGER array, dimension (8*min(M,N))

  */

  std::vector<int> IWORK(8 * minMN);

  /*


  [out]	INFO

            INFO is INTEGER
            = 0:  successful exit.
            < 0:  if INFO = -i, the i-th argument had an illegal value.
            > 0:  DBDSDC did not converge, updating process failed.



*/

  int INFO;
  /*
//    ( 	char*  	JOBZ,
//                                  int*  	M,
//                                  int*  	N,
//                                  double  *  precision, dimension( lda, * )
A,
//                                  int *  	LDA,
//                                  double * precision, dimension( * )
//S,
  //                                  double * precision, dimension( ldu, * ) U,
  //                                  int*  	LDU,
  //                                  double*  precision, dimension( ldvt, * )
VT,
  //                                  int*  	LDVT,
  //                                  double * precision, dimension( * )
WORK,
  //                                  int*  	LWORK,
  //                                  int* , dimension( * )  	IWORK,
  //                                  int*  	INFO
  //                              );
  //    */

  lapack::dgesdd_(&JOBZ, &M, &N, &Ao[0], &LDA, &S[0], &U[0], &LDU, &VT[0],
                  &LDVT, &WORK[0], &LWORK, &IWORK[0], &INFO);

  LWORK = WORK[0];
  WORK.resize(LWORK);

  lapack::dgesdd_(&JOBZ, &M, &N, &Ao[0], &LDA, &S[0], &U[0], &LDU, &VT[0],
                  &LDVT, &WORK[0], &LWORK, &IWORK[0], &INFO);

  typedef myOptional_t<SVD_type> Op;

  if (INFO != 0) {
    if (INFO > 0)
      return Op(
          false,
          std::string("DBDSDC did not converge, updating process failed.") +
              std::to_string(INFO));
    else {
      auto args =
          std::make_tuple(&JOBZ, &M, &N, &Ao[0], &LDA, &S[0], &U[0], &LDU,
                          &VT[0], &LDVT, &WORK[0], &LWORK, &IWORK[0], &INFO);
      std::string argumentsNames[] = {"JOBZ", "M",     "N",     "A",   "LDA",
                                      "S",    "U",     "LDU",   "VT",  "LDVT",
                                      "WORK", "LWORK", "IWORK", "INFO"};

      std::stringstream ss;
      std::size_t i = -INFO;
      io::write_tuple_i(ss, args, i);
      return Op(false, "the" + std::to_string(i) + "-th argument " +
                           argumentsNames[i] +
                           " had the illegal value =" + ss.str());
    }

  } else {
    U = Transpose(U);
    VT = Transpose(VT);
    assert(are_Equal_v(A, U * S * VT, std::cerr,
                       "Matrix singular decomposition of ", A));
    return Op(std::tuple(std::move(U), std::move(S), std::move(VT)));
  }
}

} // namespace Matrix_Decompositions

namespace Matrix_Unary_Functions {

template <class E, class F, typename T>
T accumulate(const M_Matrix<E> &x, const F &f, const T &start) {
  T out(start);
  for (auto it = x.begin(); it != x.end(); ++it)
    out = f(out, *it);
  return out;
}

inline double fullSum(const double &x) { return x; }
template <class T> double fullSum(const M_Matrix<T> &x) {
  return accumulate(x, [](double a, const T &b) { return a + fullSum(b); },
                    0.0);
}

inline double colSum(const double &x) { return x; }
template <class T> double colSum(const M_Matrix<T> &x) {
  return accumulate_by_Cols(x, [](double a, const T &b) { return a + b; }, 0.0);
}

inline double rowSum(const double &x) { return x; }
template <class T> auto rowSum(const M_Matrix<T> &x) {
  return accumulate_by_Rows(x, [](double a, const T &b) { return a + b; }, 0.0);
}

inline double logProduct(double x) { return log(x); }
template <class T> double logProduct(const M_Matrix<T> &x) {
  return accumulate(x, [](double a, const T &b) { return a + logProduct(b); },
                    0.0);
}





inline double maxAbs(double x) { return std::abs(x); }

template <typename T> double maxAbs(const M_Matrix<T> &x) {
  return accumulate(
      x, [](double a, const T &b) { return std::max(a, maxAbs(b)); }, 0.0);
}

inline double minAbs(double x) { return std::abs(x); }

template <typename T> double minAbs(const M_Matrix<T> &x) {
    return accumulate(
        x, [](double a, const T &b) { return std::min(a, minAbs(b)); }, std::numeric_limits<double>::infinity());
}


inline double max(const M_Matrix<double> &x) {
  using std::max;
  return accumulate(x, [](double a, double b) { return max(a, b); },
                    -std::numeric_limits<double>::infinity());
}
inline double min(const M_Matrix<double> &x) {
  using std::min;
  return accumulate(x, [](double a, double b) { return min(a, b); },
                    std::numeric_limits<double>::infinity());
}

template <typename T> double min(const M_Matrix<T> &x) {
  using std::min;
  return accumulate(x, [](double a, const T &b) { return min(a, min(b)); },
                    std::numeric_limits<double>::infinity());
}

template <typename T> double logDiagProduct(const M_Matrix<T> &x) {
  auto d = diag(x);
  return logProduct(d);
}

/**
       Product of the Diagonal of a Matrix

      */
template <typename T> double diagProduct(const M_Matrix<T> &x) {
  return std::exp(logDiagProduct(x));
}

template <typename T> double Trace(const M_Matrix<T> &x) {
  auto d = diag(x);
  return fullSum(d);
}

inline double norm_1(const double x) { return std::abs(x); }

/**
    Maximum Value of the Sum of the absolute values in each row
   */
template <typename T> double norm_1(const M_Matrix<T> &x) {
  double n = 0;
  for (size_t i = 0; i < ncols(x); ++i) {
    double sum = 0;
    for (size_t j = 0; j < nrows(x); ++j)
      sum += norm_1(x(j, i));
    n = std::max(n, sum);
  }
  return n;
}
inline double norm_1_lapack(const M_Matrix<double>& x)
{
char NORM = '1';
  int N = x.ncols();
  auto a=x;
  int M = N;

  [[maybe_unused]] int INFO = 0;
  auto WORK_lange = std::make_unique<double[]>(N);
  int LDA = N;
  using lapack::dlange_;
  double ANORM = dlange_(&NORM, &M, &N, &a[0], &LDA, WORK_lange.get());
  return ANORM;
}

/**
    Maximum Value of the Sum of the absolute values in each row
   */

inline double norm_inf_lapack(const M_Matrix<double>& x)
{
    assert(x.type()==Matrix_TYPE::FULL);
    char NORM = 'I';
    int N = x.ncols();
    M_Matrix<double> a=x;
    int M = N;

    int INFO = 0;
    auto WORK_lange = std::make_unique<double[]>(N);
    int LDA = N;
    using lapack::dlange_;
    double ANORM = dlange_(&NORM, &M, &N, &a[0], &LDA, WORK_lange.get());
    return ANORM;
}
template <typename T> T norm_inf(const M_Matrix<T> &x) {
  T n(0);
  for (size_t i = 0; i < nrows(x); ++i) {
    T sum(0);
    for (size_t j = 0; j < ncols(x); ++j)
      if (x(i, j) > 0)
        sum += x(i, j);
      else
        sum -= x(i, j);

    n = std::max(n, sum);
  }
  return n;

}

template <typename T> std::size_t log2_norm(const M_Matrix<T> &x) {
  // Scale A by power of 2 so that its norm is < 1/2 .
  double e = std::ceil(log(norm_inf(x)) / log(2));
  std::size_t s = std::max(0.0, e + 1);
  return s;
}

} // namespace Matrix_Unary_Functions

namespace Matrix_Unary_Transformations {
using Matrix_Binary_Transformations::quadraticForm_BT_A_B;

using Matrix_Binary_Transformations::quadraticForm_B_A_BT;

template <class E, class F, typename T>
M_Matrix<T> accumulate_by_Rows(const M_Matrix<E> &me, const F &f,
                               const T &start) {
  M_Matrix<T> out(me.nrows(), 1);
  switch (
      me.type()) { // ZERO,FULL,SYMMETRIC,DIAGONAL,SCALAR_FULL,SCALAR_DIAGONAL
  case Matrix_TYPE::FULL:
  case Matrix_TYPE::SYMMETRIC:
  case Matrix_TYPE::SCALAR_FULL: {
    for (std::size_t i = 0; i < me.nrows(); ++i) {
      T s = start;
      for (std::size_t j = 0; j < me.ncols(); ++j)
        s = f(s, me(i, j));
      out(i, 0) = s;
    }
    return out;
  }

  case Matrix_TYPE::DIAGONAL:
  case Matrix_TYPE::SCALAR_DIAGONAL: {
    for (std::size_t i = 0; i < me.nrows(); ++i)
      out(i, 1) = f(start, me(i, i));
    return out;
  }
  case Matrix_TYPE::ZERO:
  default: { return M_Matrix<T>(me.nrows(), 1, f(start, me[0])); }
  };
}

template <class E, class F, typename T>
M_Matrix<T> accumulate_by_Cols(const M_Matrix<E> &me, const F &f,
                               const T &start) {
  M_Matrix<T> out(1, me.ncols());
  switch (me.type()) {
  case Matrix_TYPE::FULL:
  case Matrix_TYPE::SYMMETRIC:
  case Matrix_TYPE::SCALAR_FULL: {
    for (std::size_t j = 0; j < me.ncols(); ++j) {
      T s(start);
      for (std::size_t i = 0; i < me.nrows(); ++i)
        s = f(s, me(i, j));
      out(0, j) = s;
    }
    return out;
  }

  case Matrix_TYPE::DIAGONAL:
  case Matrix_TYPE::SCALAR_DIAGONAL: {
    for (std::size_t i = 0; i < me.ncols(); ++i)
      out(1, i) = f(start, me(i, i));
    return out;
  }
  case Matrix_TYPE::ZERO: {
    return M_Matrix<T>(1, me.ncols(), f(start, me[0]));
  }
  };
}

template <class T> M_Matrix<T> Transpose(const M_Matrix<T> &x) {
  switch (x.type()) {
  case Matrix_TYPE::FULL: {
    M_Matrix<T> out(x.ncols(), x.nrows());
    for (std::size_t i = 0; i < x.nrows(); i++)
      for (std::size_t j = 0; j < x.ncols(); ++j)
        out(j, i) = x(i, j);
    return out;
  }
    // beware: does not transpose in those cases!!
  case Matrix_TYPE::ZERO:
  case Matrix_TYPE::SYMMETRIC:
  case Matrix_TYPE::DIAGONAL:
  case Matrix_TYPE::SCALAR_DIAGONAL:
  case Matrix_TYPE::SCALAR_FULL:
  default: {
    M_Matrix<T> out(x);
    return out;
  }
  }
}

template <typename T, std::size_t I0, std::size_t I1, std::size_t I2,
          std::size_t I3>
M_Matrix<M_Matrix<T>> Permute(const M_Matrix<M_Matrix<T>> &x) {
  std::tuple<std::size_t, std::size_t, std::size_t, std::size_t> sizes;

  std::get<I0>(sizes) = x.nrows();
  std::get<I1>(sizes) = x[0].nrows();
  std::get<I2>(sizes) = x.ncols();
  std::get<I3>(sizes) = x[0].ncols();

  M_Matrix<M_Matrix<T>> out(std::get<0>(sizes), std::get<2>(sizes),
                            M_Matrix<M_Matrix<T>>::FULL,
                            M_Matrix<T>(std::get<1>(sizes), std::get<3>(sizes),
                                        M_Matrix<M_Matrix<T>>::FULL));

  for (std::size_t i = 0; i < x.nrows(); ++i)
    for (std::size_t j = 0; j < x.ncols(); ++j)
      for (std::size_t i1 = 0; i1 < x(i, j).nrows(); ++i1)
        for (std::size_t j1 = 0; j1 < x(i, j).ncols(); ++j1) {
          auto index = std::tuple(i, i1, j, j1);
          auto _i = std::get<I0>(index);
          auto _i1 = std::get<I1>(index);
          auto _j = std::get<I2>(index);
          auto _j1 = std::get<I3>(index);

          out(_i, _j)(_i1, _j1) = x(i, j)(i1, j1);
        }
  return out;
}

template <class T> M_Matrix<T> TransposeSum(const M_Matrix<T> &x) {
  assert(x.ncols() == x.nrows());
  switch (x.type()) {
  case Matrix_TYPE::FULL: {
    M_Matrix<T> out(x.ncols(), x.nrows(), Matrix_TYPE::SYMMETRIC);
    for (std::size_t i = 0; i < x.nrows(); i++)
      for (std::size_t j = 0; j < i + 1; ++j)
        out(j, i) = x(i, j) + x(j, i);
    return out;
  }
  case Matrix_TYPE::ZERO:
  case Matrix_TYPE::SYMMETRIC:
  case Matrix_TYPE::DIAGONAL:
  case Matrix_TYPE::SCALAR_DIAGONAL:
  case Matrix_TYPE::SCALAR_FULL:
  default: { return x * 2.0; }
  }
}

template <typename T> M_Matrix<T> exp(const M_Matrix<T> &M) {
  using std::exp;
  return M.apply([](auto &x) { return exp(x); });
}

template <typename T>
myOptional_t<M_Matrix<T>> expm_pade(const M_Matrix<T> &M) {
  typedef myOptional_t<M_Matrix<T>> Op;

  /// http://www2.humusoft.cz/www/papers/tcp08/017_brancik.pdf
  using namespace Matrix_Generators;

  M_Matrix<T> X = M;
  double c = 0.5;
  M_Matrix<T> F = eye<T>(nrows(M)) + c * M;
  M_Matrix<T> D = eye<T>(nrows(M)) - c * M;

  std::size_t q = 6;
  bool p = true;
  for (std::size_t k = 2; k < q + 1; ++k) {
    c = c * (1.0 * q - 1.0 * k + 1.0) / (1.0 * k * (2.0 * q - k + 1.0));
    X = M * X;

    M_Matrix<T> cX = c * X;
    F += cX;
    if (p)
      D += cX;
    else
      D -= cX;
    p = !p;
  }

  auto invD = inv(D);
  if (!invD)
    return Op(false, "cannot invert D: " + invD.error());
  else {
    F = invD.value().first * F;

    return Op(F);
  }
  /*    % Pade approximation of exp(M) and diff[exp(M)]
  X=M; Y=dM;
  c=1/2;
  F=eye(size(M))+c*M; dF=c*dM;
  D=eye(size(M))-c*M; dD=-c*dM;
  q=6;
  p=1;
  for
   k=2:q
     c=c*(q-k+1)/(k*(2*q-k+1));
     Y=dM*X+M*Y;
     X=M*X;
     cX=c*X; cY=c*Y;
     F=F+cX; dF=dF+cY;
  if
   p
       D=D+cX; dD=dD+cY;
  else
       D=D-cX; dD=dD-cY;
  end
     p=~p;
  end
  F=D\F;
  dF=D\(dF-dD*F);
  % Undo scaling by repeated squaring
  for
   k=1:r
      dF=dF*F+F*dF;
      F=F*F;
  en
  */
}

inline
    M_Matrix<double> expm_taylor(const M_Matrix<double> & x, std::size_t order=6)
{
    M_Matrix<double> out=x+Matrix_Generators::eye<double>(x.ncols());
    M_Matrix<double> xr=x;
    double a=1.0;
    for (std::size_t n=2; n+1<order; ++n)
    {
        a/=n;
        xr=xr*x;
        out=out+xr*a;
    }
    return out;

}


inline
    M_Matrix<double> expm_taylor_scaling_squaring(const M_Matrix<double> & x, std::size_t order=6)
{
    double max=Matrix_Unary_Functions::maxAbs(x);
    double desired=0.125;
    int k=std::ceil(std::log2(max/desired));
    std::size_t n=std::max(0,k);
    double scale=std::pow(2,-n);
    auto dx=x*scale;
    auto expm_dx=expm_taylor(dx,order);
    auto expm_run=expm_dx;
    for (std::size_t i=0; i<n; ++i)
    {
        expm_run=expm_run*expm_run;
    }
    return expm_run;
}



/**
  Returns the exponential of a Matrix
  @remarks  uses the Pade approximation
  @pre the matrix has to be square, ncols()==nrows
  @post assert the precondition

  */

template <typename T>
myOptional_t<M_Matrix<T>> full_expm(const M_Matrix<T> &x) {
  typedef myOptional_t<M_Matrix<T>> Op;
  assert(ncols(x) == nrows(x));
  assert(size(x) > 0);
  using namespace Matrix_Generators;

  // Scale A by power of 2 so that its norm is < 1/2 .
  std::size_t s = Matrix_Unary_Functions::log2_norm(x) + 1;

  M_Matrix<T> A = x * (1.0 / std::pow(2.0, int(s)));

  // Pade approximation for exp(A)
  auto eE = expm_pade(A);

  if (!eE)
    return Op(false, eE.error());
  else {

    auto E = eE.value();
    // Undo scaling by repeated squaring
    for (std::size_t k = 0; k < s; k++)
      E = E * E;
    return Op(E);
  }
}

template <typename T>
myOptional_t<M_Matrix<T>> diagonal_expm(const M_Matrix<T> &a) {
  typedef myOptional_t<M_Matrix<T>> Op;
  M_Matrix<double> out(a.nrows(), a.ncols(), Matrix_TYPE::DIAGONAL);
  if constexpr (std::is_arithmetic_v<T>) {
    for (std::size_t i = 0; i < a.nrows(); ++i)
      out(i, i) = exp(a(i, i));
  } else
    for (std::size_t i = 0; i < a.nrows(); ++i) {

      auto o = exp(a(i, i));
      if (!o)
        return Op(false,
                  " invalid exponential of diagonal i" + std::to_string(i));
      else
        out(i, i) = std::move(o.value());
    }
  return out;
}

template <typename T> M_Matrix<T> scalar_diagonal_expm(const M_Matrix<T> &a) {
  using std::exp;
  return M_Matrix<T>(a.nrows(), a.ncols(), Matrix_TYPE::SCALAR_DIAGONAL,
                     exp(a(0, 0)));
}

template <typename E> myOptional_t<M_Matrix<E>> expm(const M_Matrix<E> &x) {
  assert(x.nrows() == x.ncols());
  switch (x.type()) {
  case Matrix_TYPE::FULL: {
    return full_expm(x);
  }
  case Matrix_TYPE::SYMMETRIC: {
    return full_expm(x);
  }
  case Matrix_TYPE::DIAGONAL: {
    return diagonal_expm(x);
  }
  case Matrix_TYPE::SCALAR_DIAGONAL: {
    return scalar_diagonal_expm(x);
  }
  case Matrix_TYPE::SCALAR_FULL: {
    return full_expm(x);
  }
  case Matrix_TYPE::ZERO:
  default: { return Matrix_Generators::eye<E>(x.nrows()); }
  }
}

template <typename T> M_Matrix<T> quadraticForm_XXT(const M_Matrix<T> &x) {
  M_Matrix<T> out(x.nrows(), x.nrows(), Matrix_TYPE::SYMMETRIC, T{});
  for (std::size_t i = 0; i < x.nrows(); ++i)
    for (std::size_t j = 0; j < i + 1; ++j)
      for (std::size_t k = 0; k < x.ncols(); ++k)
        out(i, j) += x(i, k) * x(j, k);
  return out;
}

template <typename T> M_Matrix<T> quadraticForm_XTX(const M_Matrix<T> &x) {
  M_Matrix<T> out(x.ncols(), x.ncols(), Matrix_TYPE::SYMMETRIC, T{});
  for (std::size_t i = 0; i < x.ncols(); ++i)
    for (std::size_t j = 0; j < i + 1; ++j)
      for (std::size_t k = 0; k < x.nrows(); ++k)
        out(i, j) += x(k, i) * x(k, j);
  return out;
}

template <typename T> M_Matrix<double> quadraticForm_XTX_on_vector(const std::vector<T> &x) {
    M_Matrix<double> out(x.size(), x.size(), Matrix_TYPE::SYMMETRIC);
    for (std::size_t i = 0; i < x.size(); ++i)
        for (std::size_t j = 0; j < i + 1; ++j)
                out(i, j) = x[i] * x[j];
    return out;
}



namespace partition {
template <class E> class FullPartition {
public:
  FullPartition(const M_Matrix<E> &x, std::size_t n)
      : A_(n, n), B_(n, x.ncols() - n),
        C_(x.nrows() - n, n), D_{x.nrows() - n, x.ncols() - n} {
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j)
        A_(i, j) = x(i, j);
      for (std::size_t j = n; j < x.ncols(); ++j)
        B_(i, j - n) = x(i, j);
    }
    for (std::size_t i = n; i < x.nrows(); ++i) {
      for (std::size_t j = 0; j < n; ++j)
        C_(i - n, j) = x(i, j);
      for (std::size_t j = n; j < x.nrows(); ++j)
        D_(i - n, j - n) = x(i, j);
    }
  }

  FullPartition(const M_Matrix<E> &A, const M_Matrix<E> &B,
                const M_Matrix<E> &C, const M_Matrix<E> &D)
      : A_(A), B_(B), C_(C), D_(D) {}

  FullPartition(M_Matrix<E> &&A, M_Matrix<E> &&B, M_Matrix<E> &&C,
                M_Matrix<E> &&D)
      : A_(A), B_(B), C_(C), D_(D) {}

  FullPartition() = default;

  M_Matrix<E> full() const

  {
    std::size_t N = A_.nrows() + D_.nrows();
    M_Matrix<E> out(N, N);
    std::size_t n = A_.nrows();
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j)
        out(i, j) = A_(i, j);
      for (std::size_t j = n; j < out.ncols(); ++j)
        out(i, j) = B_(i, j - n);
    }
    for (std::size_t i = n; i < out.nrows(); ++i) {
      for (std::size_t j = 0; j < n; ++j)
        out(i, j) = C_(i - n, j);
      for (std::size_t j = n; j < out.nrows(); ++j)
        out(i, j) = D_(i - n, j - n);
    }
    return out;
  }

  auto inverse_by_A() const {
    typedef myOptional_t<FullPartition> Op;
    auto invA = inv(A_);
    if (!invA) {
      return Op(false, "cannot invert A :" + invA.error());
    }
    auto AinvB = invA.value() * B_;
    auto Dinv = inv(D_ - C_ * AinvB);
    if (!Dinv) {
      return Op(false, " cannot invert D" + Dinv.error());
    }

    auto CAinv = C_ * invA.first;
    auto Binv = -AinvB * Dinv.first;
    auto Cinv = -Dinv.first * CAinv;
    auto Ainv = invA.first + AinvB * Dinv.first * CAinv;
    return Op(FullPartition(Ainv, Binv, Cinv, Dinv.value()));
  }
  std::pair<FullPartition, std::string> inverse_by_D() const {
    typedef myOptional_t<FullPartition> Op;
    auto invD = inv(D_);
    if (!invD) {
      return Op(false, " cannot invert D" + invD.second);
    }
    auto DinvC = invD.first * C_;
    auto Ainv = inv(A_ - B_ * DinvC);
    if (!Ainv) {
      return Op(false, "cannot invert A :" + Ainv.error());
    }

    auto BDinv = B_ * invD.first;
    auto Cinv = -DinvC * Ainv.first;
    auto Binv = -Ainv.first * BDinv;
    auto Dinv = invD.first + DinvC * Ainv.first * BDinv;
    return Op(FullPartition(Ainv.value(), Binv, Cinv, Dinv));
  }

private:
  M_Matrix<E> A_;
  M_Matrix<E> B_;
  M_Matrix<E> C_;
  M_Matrix<E> D_;
};

template <class E> class SymmetricPartition {
public:
  SymmetricPartition(const M_Matrix<E> &x, std::size_t n)
      : A_{n, n, Matrix_TYPE::SYMMETRIC}, B_{n, x.ncols() - n},
        D_{x.nrows() - n, x.ncols() - n, Matrix_TYPE::SYMMETRIC} {
    for (std::size_t j = 0; j < n; ++j) {
      for (std::size_t i = j; i < n; ++i)
        A_(i, j) = x(i, j);
    }
    for (std::size_t j = n; j < x.nrows(); ++j) {
      for (std::size_t i = 0; i < n; ++i)
        B_(i, j - n) = x(i, j);
      for (std::size_t i = j; i < x.nrows(); ++i)
        D_(i - n, j - n) = x(i, j);
    }
  }

  SymmetricPartition(const M_Matrix<E> &A, const M_Matrix<E> &B,
                     const M_Matrix<E> &D)
      : A_(A), B_(B), D_(D) {}
  SymmetricPartition(M_Matrix<E> &&A, M_Matrix<E> &&B, M_Matrix<E> &&D)
      : A_(A), B_(B), D_(D) {}

  SymmetricPartition() = default;

  M_Matrix<E> full() const {
    std::size_t N = A_.nrows() + D_.nrows();
    M_Matrix<E> x(N, N, Matrix_TYPE::SYMMETRIC);
    std::size_t n = A_.nrows();
    for (std::size_t j = 0; j < n; ++j) {
      for (std::size_t i = j; i < n; ++i)
        x(i, j) = A_(i, j);
    }
    for (std::size_t j = n; j < x.nrows(); ++j) {
      for (std::size_t i = 0; i < n; ++i)
        x(i, j) = B_(i, j - n);
      for (std::size_t i = j; i < x.nrows(); ++i)
        x(i, j) = D_(i - n, j - n);
    }
    return x;
  }

  myOptional_t<SymmetricPartition> inverse_by_A() const {
    typedef myOptional_t<SymmetricPartition> Op;
    auto invA = inv(A_);
    if (!invA) {
      return Op(false, "cannot invert A :" + invA.error());
    }
    auto BTAinvB = quadraticForm_BT_A_B(invA.first, B_);
    auto Dinv = inv(D_ - BTAinvB);
    if (!Dinv) {
      return Op(false, "cannot invert D :" + Dinv.error());
    }

    auto AinvB = invA.first * B_;
    auto Binv = -AinvB * Dinv.first;
    auto Ainv =
        invA.value() + quadraticForm_BT_A_B(Dinv.value(), Transpose(AinvB));
    return Op(SymmetricPartition(Ainv, Binv, Dinv.value()));
  }

  myOptional_t<SymmetricPartition> inverse_by_D() const {
    typedef myOptional_t<SymmetricPartition> Op;
    auto invD = inv(D_);
    if (!invD) {
      return Op(false, "cannot invert D :" + invD.error());
    }

    auto BDinvBT = quadraticForm_BT_A_B(invD.first, Transpose(B_));
    auto Ainv = inv(A_ - BDinvBT);
    if (!Ainv.second.empty()) {
      return {SymmetricPartition(), Ainv.second};
    }
    auto BDinv = B_ * invD.first;
    auto Binv = -Ainv.first * BDinv;
    auto Dinv = invD.first + +quadraticForm_BT_A_B(Ainv.first, BDinv);
    return {SymmetricPartition(Ainv.first, Binv, Dinv), ""};
  }

  myOptional_t<SymmetricPartition> chol_by_A(const std::string &kind) const {
    typedef myOptional_t<SymmetricPartition> Op;

    if (kind == "upper") {
      auto Ua = chol(A_, kind);
      if (!Ua)
        return Op(false, "Cholesky block La: " + Ua.error());
      auto B_chol = Forward_Sustitution_Ly_b(Transpose(Ua.first), B_);
      auto Dschur = D_ - quadraticForm_XTX(B_chol);
      auto D_chol = chol(Dschur, kind);
      if (!D_chol.second.empty())
        return Op(false, "Cholesky block Dschur: " + D_chol.error());
      auto Czero =
          M_Matrix<E>(B_chol.ncols(), B_chol.nrows(), Matrix_TYPE::ZERO);
      return Op(FullPartition<E>(Ua.value(), B_chol, Czero, D_chol.value()));
    } else {
      auto La = chol(A_, kind);
      if (!La.second.empty())
        return Op(false, "Cholesky block La: " + La.error());
      auto C_chol = Transpose(Forward_Sustitution_Ly_b(La.first, B_));
      auto Dschur = D_ - quadraticForm_XXT(C_chol);
      auto D_chol = chol(Dschur, kind);
      if (!D_chol.second.empty())
        return Op(false, "Cholesky block Dschur: " + D_chol.error());
      auto B_zero =
          M_Matrix<E>(C_chol.ncols(), C_chol.nrows(), Matrix_TYPE::ZERO);
      return Op(FullPartition<E>(La.value(), B_zero, C_chol, D_chol.value()));
    }
  }

private:
  M_Matrix<E> A_;
  M_Matrix<E> B_;
  M_Matrix<E> D_;
};
} // namespace partition
namespace matrix_inverse {
using namespace partition;
inline auto full_inv(const M_Matrix<double> &a)

{
    typedef myOptional_t<std::pair<M_Matrix<double>,double>> Op;

  const double min_inv_condition_number = 1e-12;
  // using Matrix_Unary_Transformations::Transpose;
  using lapack::dgecon_;
  using lapack::dgetrf_;
  using lapack::dgetri_;
  using lapack::dlange_;
  if (a.size() == 0)
    return Op(false, "EMPTY MATRIX");
  else {
    assert(a.ncols() == a.nrows());

    char NORM = '1';
    int N = a.ncols();
    int M = N;

    int INFO = 0;
    //  char msg[101];
    auto IPIV = std::make_unique<int[]>(N);
    int LWORK;
    //        M_Matrix<double> B=Transpose(a);
    M_Matrix<double> B = a;
    int LDA = N;
    // A=new double[n*n];
    double *A = &B[0]; // more efficient code
    auto WORK_lange = std::make_unique<double[]>(N);

    double ANORM = dlange_(&NORM, &M, &N, A, &LDA, WORK_lange.get());

    dgetrf_(&N, &M, A, &LDA, IPIV.get(), &INFO);

    double RCOND;
    auto WORK_cond = std::make_unique<double[]>(N * 4);
    auto IWORK = std::make_unique<int[]>(N);
    int INFO_con;

    dgecon_(&NORM, &N, A, &LDA, &ANORM, &RCOND, WORK_cond.get(), IWORK.get(),
            &INFO_con);

    LWORK = N * N;
    M_Matrix<double> W(N, N);
    double *WORK = &W[0];

    dgetri_(&N, A, &LDA, IPIV.get(), WORK, &LWORK, &INFO);

    if (RCOND < min_inv_condition_number)
      return Op(false, "bad condition number RCOND=" + my_to_string(RCOND));
    if (INFO == 0)
        return Op(std::make_pair(B,RCOND));
    //  return Op({B,RCOND});
    else
      return Op(false, "Singular Matrix on i=" + std::to_string(INFO));
    ;
  }
}

template <typename T> myOptional_t<std::pair<M_Matrix<T>,double>> full_inv(const M_Matrix<T> &x) {
  typedef myOptional_t<std::pair<M_Matrix<T>,double>> Op;

  auto p = FullPartition<T>(x, x.nrows() / 2);
  auto pinv = p.inverse_by_A();
  if (!pinv)
    return Op(false, pinv.error());
  else
    return Op(pinv.value().full());
}

inline myOptional_t<std::pair<M_Matrix<double>,double>>
symmetric_inv(const M_Matrix<double> &a, const std::string &kind = "lower")

{
  typedef myOptional_t<std::pair<M_Matrix<double>,double>> Op;
  using lapack::dlange_;
  using lapack::dsycon_;
  using lapack::dsytrf_;
  using lapack::dsytri_;

  if (a.size() == 0)
    return Op(false, "EMPTY MATRIX");
  else {
    assert(a.nrows() == a.ncols());

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
    if (kind == "lower")
      UPLO = 'U';
    else
      UPLO = 'L';

    int INFO = 0;
    int N = a.ncols();
    auto IPIV = std::make_unique<int[]>(N);
    int LWORK = N * N; //
    M_Matrix<double> B(a.nrows(), a.ncols());
    if (kind != "lower") {
      for (std::size_t i = 0; i < a.nrows(); ++i)
        for (std::size_t j = i; j < a.ncols(); ++j)
          B(i, j) = a(i, j);
    } else {
      for (std::size_t i = 0; i < a.nrows(); ++i)
        for (std::size_t j = 0; j < i + 1; ++j)
          B(i, j) = a(i, j);
    }

    int LDA = N;
    double *A = &B[0]; // more efficient code
    M_Matrix<double> W(N, N);
    double *WORK = &W[0];
    dsytrf_(&UPLO, &N, A, &LDA, IPIV.get(), WORK, &LWORK, &INFO);
    double RCOND;
    auto WORK_cond = std::make_unique<double[]>(N * 4);
    auto IWORK = std::make_unique<int[]>(N);
    int INFO_con;

    /**
DSYCON

Download DSYCON + dependencies [TGZ] [ZIP] [TXT]

Purpose:

DSYCON estimates the reciprocal of the condition number (in the
1-norm) of a real symmetric matrix A using the factorization
A = U*D*U**T or A = L*D*L**T computed by DSYTRF.

An estimate is obtained for norm(inv(A)), and the reciprocal of the
condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

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

[in]	A

      A is DOUBLE PRECISION array, dimension (LDA,N)
      The block diagonal matrix D and the multipliers used to
      obtain the factor U or L as computed by DSYTRF.

[in]	LDA

      LDA is INTEGER
      The leading dimension of the array A.  LDA >= max(1,N).

[in]	IPIV

      IPIV is INTEGER array, dimension (N)
      Details of the interchanges and the block structure of D
      as determined by DSYTRF.

[in]	ANORM

      ANORM is DOUBLE PRECISION
      The 1-norm of the original matrix A.

[out]	RCOND

      RCOND is DOUBLE PRECISION
      The reciprocal of the condition number of the matrix A,
      computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
      estimate of the 1-norm of inv(A) computed in this routine.

[out]	WORK

      WORK is DOUBLE PRECISION array, dimension (2*N)

[out]	IWORK

      IWORK is INTEGER array, dimension (N)

[out]	INFO

      INFO is INTEGER
      = 0:  successful exit
      < 0:  if INFO = -i, the i-th argument had an illegal value*/

    char NORM = '1';
    int M = N;

    auto WORK_lange = std::make_unique<double[]>(N);

    double ANORM = dlange_(&NORM, &M, &N, A, &LDA, WORK_lange.get());

    //        dsycon_( 	char*  	UPLO, int *  	N,double */* precision,
    //        dimension( lda, * ) */  	A,
    //                                    int *  	LDA,double * ANORM,
    //                                    double *   	RCOND, double */*
    //                                    dimension( * )  */	WORK, int */*
    //                                    dimension( * ) */ 	IWORK, int *
    //                                    INFO);

    dsycon_(&UPLO, &N, A, &LDA, IPIV.get(), &ANORM, &RCOND, WORK_cond.get(),
            IWORK.get(), &INFO_con);

    if (INFO < 0) {
      std::string argNames[] = {"UPLO",  "N",     "A",         "LDA",   "IPIV",
                                "ANORM", "RCOND", "WORK_cond", "IWORK", "INFO"};
      return Op(false,
                "INVALID ARGUMENT " + std::to_string(INFO) + argNames[-INFO]);
    } else if (INFO > 0) {
      return Op(false, "SINGULAR MATRIX ON " + std::to_string(INFO));
    } else {
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
      dsytri_(&UPLO, &N, A, &LDA, IPIV.get(), WORK, &INFO);

      if (INFO != 0) {
        return Op(false,
                  "cannot invert a singular matrix " + std::to_string(INFO));
      } else {
        M_Matrix<double> out(B.nrows(), B.ncols(), Matrix_TYPE::SYMMETRIC);
        if (kind == "lower") {
          for (std::size_t i = 0; i < B.nrows(); ++i)
            for (std::size_t j = 0; j <= i; ++j)
              out(i, j) = B(i, j);
        } else {
          for (std::size_t i = 0; i < B.nrows(); ++i)
            for (std::size_t j = 0; j <= i; ++j)
              out(i, j) = B(j, i);
        }

        /*    auto aa=a;
         M_Matrix<double> test(a.nrows(),a.ncols(),Matrix_TYPE::FULL,a);

         auto invtest=inv(test).first;

         auto test_it=test*out;
         auto test_a=aa*out;
         auto invv=invtest*test;
     */
        return Op(std::make_pair(out,RCOND));
      }
    }
  }
}

template <typename T>
myOptional_t<M_Matrix<T>> symmetric_inv(const M_Matrix<T> &x) {
  typedef myOptional_t<M_Matrix<T>> Op;
  auto p = SymmetricPartition<T>(x, x.nrows() / 2);
  auto pinv = p.inverse_by_A();
  if (!pinv)
    return Op(false, pinv.error());
  else
    return Op(pinv.value().full());
}

inline auto diagonal_inv(const M_Matrix<double> &a) {
    typedef myOptional_t<std::pair<M_Matrix<double>,double>> Op;
  M_Matrix<double> out(a.nrows(), a.ncols(), Matrix_TYPE::DIAGONAL);
  double amax = 0;
  double imax = 0;
  for (std::size_t i = 0; i < a.nrows(); ++i)
    if (a(i, i) == 0)
      return Op(false,
                "Singular Matrix, Diagonal zero at " + std::to_string(i));
    else {
      out(i, i) = 1.0 / a(i, i);
      if (std::abs(a(i, i)) > amax)
        amax = std::abs(a(i, i));
      if (std::abs(out(i, i)) > imax)
        imax = std::abs(out(i, i));
    }
  return Op(std::make_pair(out,imax/amax));
}

template <typename T>
auto diagonal_pinv(const M_Matrix<T> &a, double min_value = 0) {
  M_Matrix<double> out(a.ncols(), a.nrows(), Matrix_TYPE::DIAGONAL);
  for (std::size_t i = 0; i < std::min(a.nrows(), a.ncols()); ++i)
    if (std::abs(a(i, i)) <= min_value)
      out(i, i) = 0;
    else {
      out(i, i) = 1.0 / a(i, i);
    }
  return out;
}

inline myOptional_t<M_Matrix<double>> full_pinv(const M_Matrix<double> &A) {

  typedef myOptional_t<M_Matrix<double>> Op;
  auto res = Matrix_Decompositions::Singular_Value_decomposition(A);
  if (!res.has_value())
    return Op(false, "SVD fails: " + res.error());
  else {
    auto [U, S, VT] = std::move(res).value();
    auto Spinv = diagonal_pinv(
        S, std::sqrt(std::numeric_limits<double>::epsilon()) * S[0]);
    auto VTS = TranspMult(VT, Spinv);
    auto out = multTransp(VTS, U);
    assert(are_Equal_v(A, A * out * A, std::cerr, "A =", A));
    assert(are_Equal_v(out, out * A * out, std::cerr, "A =", A));
    assert(are_Equal_v(Transpose(A * out), A * out, std::cerr));
    assert(are_Equal_v(Transpose(out * A), out * A, std::cerr));

    return Op(out);
  }
}

inline myOptional_t<M_Matrix<double>> symm_pinv(const M_Matrix<double> &A) {

  typedef myOptional_t<M_Matrix<double>> Op;
  auto res = Matrix_Decompositions::EigenSystem_symm_real_eigenvalue(A);
  if (!res.has_value())
    return Op(false, "Eigendecomposition fails: " + res.error());
  else {
    auto [VR, L, VL,RCOND_L,RCOND_V] = std::move(res).value();
    assert(are_Equal_v(Matrix_Generators::eye<double>(A.ncols()), VL * VR, std::cerr," VL *VR =", VL * VR));
    assert(are_Equal_v(A, VR * L * VL, std::cerr));
    auto Lpinv = diagonal_pinv(
        L,  std::sqrt(std::numeric_limits<double>::epsilon()) * std::abs(L[L.size() - 1]));
    auto Apinv = VR*Lpinv*VL;
    assert((are_Equal_v(A, A * Apinv * A, 0.1,0.1,std::cerr, "\nA =", A, " \nL=", L, " \nLpinv=", Lpinv, " \nVL=", VL, " \nVR=", VR)));
    assert(are_Equal_v(Apinv, Apinv * A * Apinv, 0.1,0.1,std::cerr, "A =", A));
    assert(are_Equal_v(Transpose(A * Apinv), A * Apinv, 0.1,0.1,std::cerr));
    assert(are_Equal_v(Transpose(Apinv * A), Apinv * A, 0.1,0.1,std::cerr));

    return Op(Apinv);
  }
}

template <typename T> auto diagonal_inv(const M_Matrix<T> &a) {
  typedef myOptional_t<M_Matrix<T>> Op;
  double amax = 0;
  double imax = 0;
  M_Matrix<T> out(a.nrows(), a.ncols(), Matrix_TYPE::DIAGONAL);
  for (std::size_t i = 0; i < a.nrows(); ++i) {
    double anorm = norm_1(a(i, i));
    auto e = inv(a(i, i));
    if (!e) {
      return Op(false, "Singular block=" + std::to_string(i) + e.error());
    } else {
      out(i, i) = std::move(e.first);
      double inorm = 1.0 / (anorm * e.second.first);
      if (anorm > amax)
        amax = anorm;
      if (inorm > imax)
        imax = inorm;
    }
  }
  //   return {out,{1.0/(amax*imax),""}};
  return Op(out);
}

inline auto scalar_diagonal_inv(const M_Matrix<double> &a) {
  typedef myOptional_t<std::pair<M_Matrix<double>,double>> Op;

  if (a(0, 0) == 0)
    return Op(false, "Singular Matrix, Diagonal is zero");
  else
      return Op(std::make_pair(M_Matrix<double>(a.nrows(), a.ncols(),
                                                Matrix_TYPE::SCALAR_DIAGONAL, 1.0 / a(0, 0)),1.0));
}

template <typename T> auto scalar_diagonal_inv(const M_Matrix<T> &a) {
  typedef myOptional_t<M_Matrix<T>> Op;
  auto e = inv(a(0, 0));
  if (!e.second.empty()) {
    return Op(false, "Singular SCALAR DIAGONAL block: " + e.error());
  } else
    return Op(M_Matrix<T>(a.nrows(), a.ncols(), Matrix_TYPE::SCALAR_DIAGONAL,
                          std::move(e.first)));
}

template <typename E> auto Matrix_inverse(const M_Matrix<E> &x) {
    typedef myOptional_t<std::pair<M_Matrix<E>,double>> Op;
  assert(x.nrows() == x.ncols());
  switch (x.type()) {
  case Matrix_TYPE::FULL: {
    return full_inv(x);
  }
  case Matrix_TYPE::SYMMETRIC: {
    return symmetric_inv(x);
  }
  case Matrix_TYPE::DIAGONAL: {
    return diagonal_inv(x);
  }
  case Matrix_TYPE::SCALAR_DIAGONAL: {
    return scalar_diagonal_inv(x);
  }
  case Matrix_TYPE::SCALAR_FULL: {
    return Op(false, "SCALAR FULL is Singular");
  }
  case Matrix_TYPE::ZERO:
  default: { return Op(false, "ZERO Matrix is is Singular"); }
  }
}

template <typename E> auto Matrix_pseudo_inverse(const M_Matrix<E> &x) {
  typedef myOptional_t<M_Matrix<E>> Op;
  switch (x.type()) {
  case Matrix_TYPE::FULL: {
    return full_pinv(x);
  }
  case Matrix_TYPE::SYMMETRIC: {
    return symm_pinv(x);
  }
  case Matrix_TYPE::DIAGONAL: {
    return Op(diagonal_pinv(x));
  }
  case Matrix_TYPE::SCALAR_DIAGONAL: {
    return Op(diagonal_pinv(x));
  }
  case Matrix_TYPE::SCALAR_FULL: {
    return Op(false, "SCALAR FULL is Singular");
  }
  case Matrix_TYPE::ZERO:
  default: { return Op(Transpose(x)); }
  }
}

} // namespace matrix_inverse

namespace cholesky {
using Matrix_Unary_Transformations::Transpose;
using namespace partition;
inline myOptional_t<M_Matrix<double>> symmetric_chol(const M_Matrix<double> &x,
                                                     const std::string &kind) {
  typedef myOptional_t<M_Matrix<double>> Op;
  assert(x.nrows() == x.ncols());

  if (x.size() == 0)
    return {{}, "cholesky of ZERO MATRIX"};
  char UPLO = 'L';
  M_Matrix<double> res;
  if (kind == "lower") {
    UPLO = 'U';
    res = lapack::LT(x);
  } else {
    res = lapack::UT(x);
  }
  int N = x.nrows();
  int LDA = N;
  int INFO;

  lapack::dpotrf_(&UPLO, &N, &res[0], &LDA, &INFO);

  if (INFO != 0) {
    return Op(false, "Cholesky fails, zero diagonal at" + std::to_string(INFO));
  } else {

    assert((kind == "lower"
                ? (are_Equal<true, M_Matrix<double>>(
                       std::numeric_limits<double>::epsilon() * 1000,
                       std::numeric_limits<double>::epsilon() * 1000)
                       .test_sum(res * Transpose(res), x, std::cerr))
                : (are_Equal<true, M_Matrix<double>>(
                       std::numeric_limits<double>::epsilon() * 1000,
                       std::numeric_limits<double>::epsilon() * 1000)
                       .test_sum(Transpose(res) * res, x, std::cerr))));
    return Op(res);
  }
}

template <typename T>
myOptional_t<M_Matrix<T>> symmetric_chol(const M_Matrix<T> &x,
                                         const std::string &kind) {
  typedef myOptional_t<M_Matrix<T>> Op;

  auto p = SymmetricPartition<T>(x, x.nrows() / 2);
  auto pchol = p.chol_by_A(kind);
  if (!pchol.second.empty())
    return Op(false, "Block Cholesky fails " + pchol.second);
  else
    return Op(pchol.value().full());
}

inline myOptional_t<M_Matrix<double>>
diagonal_chol(const M_Matrix<double> &a, const std::string & /*kind*/) {
  typedef myOptional_t<M_Matrix<double>> Op;

  M_Matrix<double> out(a.nrows(), a.ncols(), Matrix_TYPE::DIAGONAL);
  for (std::size_t i = 0; i < a.nrows(); ++i)
    if (a(i, i) <= 0)
      return Op(false, "Negative Diagonal zero at " + std::to_string(i) +
                           "value " + std::to_string(a(i, i)));
    else
      out(i, i) = std::sqrt(a(i, i));
  return Op(out);
}

template <typename T>
myOptional_t<M_Matrix<T>> diagonal_chol(const M_Matrix<T> &a,
                                        const std::string &kind) {
  typedef myOptional_t<M_Matrix<T>> Op;

  M_Matrix<T> out(a.nrows(), a.ncols(), Matrix_TYPE::DIAGONAL);
  for (std::size_t i = 0; i < a.nrows(); ++i) {
    auto e = chol(a(i, i), kind);
    if (!e.second.empty()) {
      return Op(false, "Non positive definite  block=" + std::to_string(i) +
                           "  " + e.error());
    } else
      out(i, i) = std::move(e.first);
  }
  return Op(out);
}

inline myOptional_t<M_Matrix<double>>
scalar_diagonal_chol(const M_Matrix<double> &a, const std::string & /*kind*/) {
  typedef myOptional_t<M_Matrix<double>> Op;

  if (a(0, 0) <= 0)
    return Op(false, "Non Definite positive, Diagonal is non positive");
  else
    return Op(M_Matrix<double>(a.nrows(), a.ncols(),
                               Matrix_TYPE::SCALAR_DIAGONAL,
                               std::sqrt(a(0, 0))));
}

template <typename T>
myOptional_t<M_Matrix<T>> scalar_diagonal_chol(const M_Matrix<T> &a,
                                               const std::string &kind) {
  typedef myOptional_t<M_Matrix<T>> Op;

  auto e = chol(a(0, 0), kind);
  if (!e.second.empty()) {
    return Op(false,
              "Not positive definite SCALAR DIAGONAL block: " + e.second);
  } else
    return Op(M_Matrix<T>(a.nrows(), a.ncols(), Matrix_TYPE::SCALAR_DIAGONAL,
                          std::move(e.first)));
  ;
}

template <typename E>
myOptional_t<M_Matrix<E>> Matrix_cholesky(const M_Matrix<E> &x,
                                          const std::string &kind)

{
  typedef myOptional_t<M_Matrix<E>> Op;
  assert(x.nrows() == x.ncols());
  switch (x.type()) {
  case Matrix_TYPE::FULL: {
    return Op(false, "Error, Cholesky on a FULL MATRIX");
  }
  case Matrix_TYPE::SYMMETRIC: {
    return symmetric_chol(x, kind);
  }
  case Matrix_TYPE::DIAGONAL: {
    return diagonal_chol(x, kind);
  }
  case Matrix_TYPE::SCALAR_DIAGONAL: {
    return scalar_diagonal_chol(x, kind);
  }
  case Matrix_TYPE::SCALAR_FULL: {
    return {{}, "SCALAR FULL is Not definite positive"};
  }
  case Matrix_TYPE::ZERO:
  default: { return {{}, "ZERO Matrix is Not definite positive"}; }
  }
}

} // namespace cholesky

template <class T> myOptional_t<std::pair<M_Matrix<T>,double>> inv(const M_Matrix<T> &x) {
  return matrix_inverse::Matrix_inverse(x);
}

template <class T> myOptional_t<M_Matrix<T>> pinv(const M_Matrix<T> &x) {
  return matrix_inverse::Matrix_pseudo_inverse(x);
}

template <class T>
myOptional_t<M_Matrix<T>> chol(const M_Matrix<T> &x, const std::string &kind) {
  return cholesky::Matrix_cholesky(x, kind);
}

template <typename T, class Compare>
M_Matrix<T> sort(const M_Matrix<T> &x, Compare comp) {
  std::vector<T> o = x.toVector();
  std::sort(o.begin(), o.end(), comp);
  return M_Matrix<T>(x.nrows(), x.ncols(), o);
}

template <typename T> M_Matrix<T> sort(const M_Matrix<T> &x) {
  std::vector<T> o = x.toVector();
  std::sort(o.begin(), o.end());
  return M_Matrix<T>(x.nrows(), x.ncols(), o);
}

/**
   Diagonal of Matrix or Diagonal Matrix
   It has two behaviors:
   - If the input is a single column or a single row, it builds a diagonal
   Matrix with it
   - If the input is a Matrix, it returns the values of its diagonal

  */
template <typename T> M_Matrix<T> diag(const M_Matrix<T> &x) {
  size_t nr = x.nrows();
  size_t nc = x.ncols();
  if ((nr > 1) && (nc > 1)) {
    std::size_t n = std::min(nr, nc);
    M_Matrix<T> diagV(nr, nc, Matrix_TYPE::DIAGONAL);
    for (size_t i = 0; i < n; ++i)
      diagV(i, i) = x(i, i);
    return diagV;
  } else {
    nr = std::max(nr, nc);
    M_Matrix<T> diagM(nr, nr, Matrix_TYPE::DIAGONAL);
    for (size_t i = 0; i < nr; ++i)
      diagM(i, i) = x[i];
    return diagM;
  }
}

template <typename T> M_Matrix<T> col_vector(const M_Matrix<T> &x) {
  M_Matrix<T> colvec(x.size(), 1);
  for (std::size_t i = 0; i < x.size(); ++i)
    colvec[i] = x[i];
  return colvec;
}
template <typename T> M_Matrix<T> row_vector(const M_Matrix<T> &x) {
  M_Matrix<T> rowvec(1, x.size());
  for (std::size_t i = 0; i < x.size(); ++i)
    rowvec[i] = x[i];
  return rowvec;
}

} // namespace Matrix_Unary_Transformations

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

namespace Matrix_Binary_Predicates {
using namespace Matrix_Unary_Transformations;

/*
Matrix Equality Operator
*/
template <typename T>
bool operator==(const M_Matrix<T> &x, const M_Matrix<T> &y) {
  if (x.type() != y.type())
    return false;
  if (x.size() != y.size())
    return false;
  else if (x.ncols() != y.ncols())
    return false;
  else
    for (std::size_t i = 0; i < x.size(); ++i)
      if (x[i] != y[i])
        return false;
  return true;
}

/*
 Minor Operator based on a Lexicographic Comparison.
*/
template <typename T>
bool operator<(const M_Matrix<T> &x, const M_Matrix<T> &y) {
  if (x.size() < y.size())
    return true;
  else if (y.size() < x.size())
    return false;
  else if (x.nrows() < y.nrows())
    return true;
  else if (y.nrows() < x.nrows())
    return false;
  else
    for (std::size_t i = 0; i < x.size(); ++i) {
      if (x[i] < y[i])
        return true;
      else if (y[i] < x[i])
        return false;
    }
  return false;
}

template <typename T> bool operator<(const M_Matrix<T> &x, const T &y) {
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (x[i] >= y)
      return false;
  }
  return true;
}
template <typename T> bool operator<=(const M_Matrix<T> &x, const T &y) {
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (x[i] > y)
      return false;
  }
  return true;
}

template <typename T> bool operator>=(const M_Matrix<T> &x, const T &y) {
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (x[i] < y)
      return false;
  }
  return true;
}
template <typename T> bool operator>(const M_Matrix<T> &x, const T &y) {
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (x[i] <= y)
      return false;
  }
  return true;
}

} // namespace Matrix_Binary_Predicates

namespace Matrix_Binary_Transformations {
using namespace Matrix_Unary_Transformations;

template <typename T, typename S>
typename M_Matrix<T>::TYPE result_of_sum(const M_Matrix<T> &x,
                                         const M_Matrix<S> &y) {
  if ((x.type() == Matrix_TYPE::FULL) || (y.type() == Matrix_TYPE::FULL))
    return Matrix_TYPE::FULL;
  else if ((x.type() == Matrix_TYPE::SYMMETRIC) ||
           (y.type() == Matrix_TYPE::SYMMETRIC))
    return Matrix_TYPE::SYMMETRIC;
  else if ((x.type() == Matrix_TYPE::DIAGONAL) ||
           (y.type() == Matrix_TYPE::DIAGONAL)) {
    if ((x.type() == Matrix_TYPE::SCALAR_FULL) ||
        (y.type() == Matrix_TYPE::SCALAR_FULL))
      return Matrix_TYPE::SYMMETRIC;
    else
      return Matrix_TYPE::DIAGONAL;
  } else if ((x.type() == Matrix_TYPE::SCALAR_DIAGONAL) ||
             (y.type() == Matrix_TYPE::SCALAR_DIAGONAL)) {
    if ((x.type() == Matrix_TYPE::SCALAR_FULL) ||
        (y.type() == Matrix_TYPE::SCALAR_FULL))
      return Matrix_TYPE::SYMMETRIC;
    else
      return Matrix_TYPE::SCALAR_DIAGONAL;
  } else if ((x.type() == Matrix_TYPE::SCALAR_FULL) ||
             (y.type() == Matrix_TYPE::SCALAR_FULL))
    return Matrix_TYPE::SCALAR_FULL;
  else
    return Matrix_TYPE::ZERO;
}

/**
 *  use only when x.size()>=y.size()
 */
template <typename T, typename S>
typename M_Matrix<S>::ITER_TYPE iterator_of_sum(const M_Matrix<T> &x,
                                                const M_Matrix<S> &y) {
  switch (x.type()) {
  case Matrix_TYPE::FULL: {
    switch (y.type()) {
    case Matrix_TYPE::FULL:
    case Matrix_TYPE::SYMMETRIC:
    case Matrix_TYPE::SCALAR_FULL:
      return M_Matrix<S>::FULL_IT;
    case Matrix_TYPE::DIAGONAL:
    case Matrix_TYPE::SCALAR_DIAGONAL:
      return M_Matrix<S>::DIAG_IT;
    case Matrix_TYPE::ZERO:
    default:
      return M_Matrix<S>::ZERO_IT;
    }
  }
  case Matrix_TYPE::SYMMETRIC: {
    switch (y.type()) {
    case Matrix_TYPE::FULL:
      assert(false);
    case Matrix_TYPE::SYMMETRIC:
    case Matrix_TYPE::SCALAR_FULL:
      return M_Matrix<S>::LT_IT;
    case Matrix_TYPE::DIAGONAL:
    case Matrix_TYPE::SCALAR_DIAGONAL:
      return M_Matrix<S>::DIAG_IT;
    case Matrix_TYPE::ZERO:
    default:
      return M_Matrix<S>::ZERO_IT;
    }
  }
  case Matrix_TYPE::DIAGONAL: {
    switch (y.type()) {
    case Matrix_TYPE::FULL:
    case Matrix_TYPE::SYMMETRIC:
    case Matrix_TYPE::SCALAR_FULL:
      assert(false);
    case Matrix_TYPE::DIAGONAL:
    case Matrix_TYPE::SCALAR_DIAGONAL:
      return M_Matrix<S>::DIAG_IT;
    case Matrix_TYPE::ZERO:
    default:
      return M_Matrix<S>::ZERO_IT;
    }
  }
  case Matrix_TYPE::SCALAR_DIAGONAL: {
    switch (y.type()) {
    case Matrix_TYPE::FULL:
    case Matrix_TYPE::SYMMETRIC:
    case Matrix_TYPE::SCALAR_FULL:
    case Matrix_TYPE::DIAGONAL:
      assert(false);
    case Matrix_TYPE::SCALAR_DIAGONAL:
      return M_Matrix<S>::SCALAR_IT;
    case Matrix_TYPE::ZERO:
    default:
      return M_Matrix<S>::ZERO_IT;
    }
  }
  case Matrix_TYPE::SCALAR_FULL: {
    switch (y.type()) {
    case Matrix_TYPE::FULL:
    case Matrix_TYPE::SYMMETRIC:
    case Matrix_TYPE::DIAGONAL:
    case Matrix_TYPE::SCALAR_DIAGONAL:
      assert(false);
    case Matrix_TYPE::SCALAR_FULL:
      return M_Matrix<S>::SCALAR_IT;
    case Matrix_TYPE::ZERO:
    default:
      return M_Matrix<S>::ZERO_IT;
    }
  }
  case Matrix_TYPE::ZERO:
  default: {
    switch (y.type()) {
    case Matrix_TYPE::FULL:
    case Matrix_TYPE::SYMMETRIC:
    case Matrix_TYPE::DIAGONAL:
    case Matrix_TYPE::SCALAR_DIAGONAL:
    case Matrix_TYPE::SCALAR_FULL:
      assert(false);
    case Matrix_TYPE::ZERO:
    default:
      return M_Matrix<S>::ZERO_IT;
    }
  }
  }
}

template <typename T, typename S>
std::pair<std::size_t, std::size_t>
Size_of_Product(const M_Matrix<T> &one, const M_Matrix<S> &other,
                bool transposeOne, bool transposeOther) {
  /*  assert(
              (transposeOne ? one.nrows():one.ncols())
              ==
              (transposeOther ? other.ncols():other.nrows())
              );*/
  return {transposeOne ? one.ncols() : one.nrows(),
          transposeOther ? other.nrows() : other.ncols()};
}

template <typename T, typename S>
auto result_of_Product(const M_Matrix<T> &x, const M_Matrix<S> &y,
                       bool transpose_x, bool transpose_y)
    -> M_Matrix<decltype(std::declval<T>() * std::declval<S>())> {
  typedef decltype(std::declval<T>() * std::declval<S>()) R;
  auto s = Size_of_Product(x, y, transpose_x, transpose_y);

  if ((x.type() == Matrix_TYPE::ZERO) || (y.type() == Matrix_TYPE::ZERO))
    return M_Matrix<R>(s.first, s.second, Matrix_TYPE::ZERO);
  if ((x.type() == Matrix_TYPE::FULL) || (y.type() == Matrix_TYPE::FULL))
    return M_Matrix<R>(s.first, s.second, Matrix_TYPE::FULL);

  else if ((x.type() == Matrix_TYPE::SYMMETRIC) ||
           (y.type() == Matrix_TYPE::SYMMETRIC)) {
    if ((x.type() == Matrix_TYPE::SCALAR_DIAGONAL) ||
        (y.type() == Matrix_TYPE::SCALAR_DIAGONAL))
      return M_Matrix<R>(s.first, s.second, Matrix_TYPE::SYMMETRIC);
    else
      return M_Matrix<R>(s.first, s.second, Matrix_TYPE::FULL);

  } else if ((x.type() == Matrix_TYPE::DIAGONAL) ||
             (y.type() == Matrix_TYPE::DIAGONAL)) {
    if ((x.type() == Matrix_TYPE::SCALAR_FULL) ||
        (y.type() == Matrix_TYPE::SCALAR_FULL))
      return M_Matrix<R>(s.first, s.second, Matrix_TYPE::FULL);
    else
      return M_Matrix<R>(s.first, s.second, Matrix_TYPE::DIAGONAL);
  } else if ((x.type() == Matrix_TYPE::SCALAR_DIAGONAL) &&
             (y.type() == Matrix_TYPE::SCALAR_DIAGONAL))

    return M_Matrix<R>(s.first, s.second, Matrix_TYPE::SCALAR_DIAGONAL);
  else
    return M_Matrix<R>(s.first, s.second, Matrix_TYPE::SCALAR_FULL);
}

namespace product {
template <typename T, typename S>
auto UnformattedProduct(const M_Matrix<T> &x, const M_Matrix<S> &y)
    -> M_Matrix<decltype(std::declval<T>() * std::declval<S>())> {
  typedef decltype(std::declval<T>() * std::declval<S>()) R;
  assert(x.ncols() == y.nrows());
  M_Matrix<decltype(std::declval<T>() * std::declval<S>())> out(
      x.nrows(), y.ncols(), R());
  for (std::size_t i = 0; i < x.nrows(); ++i)
    for (std::size_t j = 0; j < y.ncols(); ++j)
      for (std::size_t k = 0; k < x.ncols(); ++k)
        out(i, j) += x(i, k) * y(k, j);
  return out;
}

inline M_Matrix<double> Lapack_Full_Product(const M_Matrix<double> &x,
                                            const M_Matrix<double> &y,
                                            bool transpose_x,
                                            bool transpose_y) {
  using lapack::dgemm_;

  std::size_t cols_i, rows_e, cols_e;
  if (transpose_x) {
    rows_e = x.ncols();
    cols_i = x.nrows();
  } else {
    rows_e = x.nrows();
    cols_i = x.ncols();
  }

  if (transpose_y) {
    //  rows_i=y.ncols();
    cols_e = y.nrows();
  } else {
    //  rows_i=y.nrows();
    cols_e = y.ncols();
  }

  // assert(rows_i==cols_i);
  assert(((transpose_y ? y.ncols() : y.nrows()) == cols_i));
  // First it has to find out if the last dimension of x matches the
  // first of y
  // now we build the M_Matrix result
  M_Matrix<double> z(rows_e, cols_e, 0.0);

  /***  as fortran uses the reverse order for matrices and we want to
          avoid a copying operation, we calculate
              Transpose(Z)=Transpose(y)*Transpose(x)

              Transpose(matrix)=just plain matrix in C++ format


          */
  char TRANSA;
  char TRANSB;

  if (transpose_y)
    TRANSA = 'T';
  else
    TRANSA = 'N';

  if (transpose_x)
    TRANSB = 'T';
  else
    TRANSB = 'N';

  int M = cols_e;
  int N = rows_e;
  int K = cols_i;

  double ALPHA = 1.0;
  double *A = const_cast<double *>(&y[0]);
  int LDA;
  if (transpose_y)
    LDA = K;
  else
    LDA = M;

  double *B = const_cast<double *>(&x[0]);

  int LDB;
  if (transpose_x)
    LDB = N;
  else
    LDB = K;

  double BETA = 0.0;

  double *C = &z[0];

  int LDC = M;

  try {
    dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C,
           &LDC);
  } catch (...) {
    assert(false);
  }
  return z;
}

inline M_Matrix<double> Lapack_Symmetric_Regular_Product(
    const M_Matrix<double> &Sym, const M_Matrix<double> &Reg,
    bool SymReg = true, const std::string &kind = "upper") {
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
  if (kind == "lower")
    UPLO = 'U';
  else
    UPLO = 'L';

  M_Matrix<double> S(Sym.nrows(), Sym.ncols());
  if (kind == "lower") {
    for (std::size_t i = 0; i < Sym.nrows(); i++)
      for (std::size_t j = 0; j <= i; ++j)
        S(i, j) = Sym(i, j);
  } else {
    for (std::size_t i = 0; i < Sym.nrows(); i++)
      for (std::size_t j = i; j < Sym.ncols(); ++j)
        S(i, j) = Sym(i, j);
  }

  std::size_t rows_e;
  // std::size_t cols_i;
  // std::size_t rows_i;
  std::size_t cols_e;

  if (SymReg) // calculo Sym * Reg
  {
    rows_e = Sym.nrows();
    //    cols_i=Sym.ncols();
    //    rows_i=Reg.nrows();
    cols_e = Reg.ncols();

  } else // calculo  Reg^T * Sym
  {
    rows_e = Reg.nrows();
    //    cols_i=Reg.ncols();
    //    rows_i=Sym.nrows();
    cols_e = Sym.ncols();
  }
  //  assert(rows_i==cols_i);
  assert((SymReg ? Reg.nrows() == Sym.ncols() : Sym.nrows() == Reg.ncols()));
  // First it has to find out if the last dimension of x matches the
  // first of y
  // now we build the M_Matrix result
  M_Matrix<double> z(rows_e, cols_e, 0.0);

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
    SIDE = 'R';
  else
    SIDE = 'L';

  /**
   * @brief M
    [in]	M

              M is INTEGER
               On entry,  M  specifies the number of rows of the matrix  C.
               M  must be at least zero.
   */

  int M = cols_e;

  /**
   * @brief N
    [in]	N

              N is INTEGER
               On entry, N specifies the number of columns of the matrix C.
               N  must be at least zero.

               como uso la transpuesta de C, es rows_e

   */

  int N = rows_e;

  /**
   * @brief ALPHA
    [in]	ALPHA

              ALPHA is REAL
               On entry, ALPHA specifies the scalar alpha.
   */
  double ALPHA = 1.0;

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

  double * /*, dimension(lda,*)*/ A = const_cast<double *>(&S[0]);
  int LDA = S.ncols();
  double *B = const_cast<double *>(&Reg[0]);

  int LDB = Reg.ncols();
  double BETA = 0;
  double * /*, dimension(ldc,*) */ C;
  C = &z[0];
  int LDC = M;

  dsymm_(&SIDE, &UPLO, &M, &N, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);
  return z;
}

inline M_Matrix<double>
Lapack_Transpose_Symmetric_Product(const M_Matrix<double> &Reg,
                                   const M_Matrix<double> &Sym,
                                   const std::string kind = "upper") {
  return Lapack_Symmetric_Regular_Product(Sym, Reg, false, kind);
}

template <typename T, typename S>
auto All_Product(const M_Matrix<T> &x, const M_Matrix<S> &y, bool transpose_x,
                 bool transpose_y)
    -> M_Matrix<decltype(std::declval<T>() * std::declval<S>())> {

  auto out = result_of_Product(x, y, transpose_x, transpose_y);
  if (out.type()!=Matrix_TYPE::ZERO)
  {
  if (x.isDiagonal()) {
    if (!transpose_y) {
      for (auto it = out.begin(); it != out.end(); ++it) {
        auto i = it.iRow();
        auto j = it.jCol();
        *it = x(i, i) * y(i, j);
      }
    } else {
      for (auto it = out.begin(); it != out.end(); ++it) {
        auto i = it.iRow();
        auto j = it.jCol();
        *it = x(i, i) * y(j, i);
      }
    }
  } else if (y.isDiagonal()) {
    if (!transpose_x) {
      for (auto it = out.begin(); it != out.end(); ++it) {
        auto i = it.iRow();
        auto j = it.jCol();
        *it = x(i, j) * y(j, j);
      }
    } else {
      for (auto it = out.begin(); it != out.end(); ++it) {
        auto i = it.iRow();
        auto j = it.jCol();
        *it = x(j, i) * y(j, j);
      }
    }
  } else if (!transpose_x && !transpose_y) {
    for (auto it = out.begin(); it != out.end(); ++it) {
      auto i = it.iRow();
      auto j = it.jCol();
      *it = x(i, 0) * y(0, j);
      for (std::size_t k = 1; k < x.ncols(); ++k)
        *it += x(i, k) * y(k, j);
    }
  } else if (transpose_x && !transpose_y) {
    for (auto it = out.begin(); it != out.end(); ++it) {
      auto i = it.iRow();
      auto j = it.jCol();
      *it = x(0, i) * y(0, j);
      for (std::size_t k = 1; k < x.nrows(); ++k)
        *it += x(k, i) * y(k, j);
    }
  } else if (!transpose_x && transpose_y) {
    for (auto it = out.begin(); it != out.end(); ++it) {
      auto i = it.iRow();
      auto j = it.jCol();
      *it = x(i, 0) * y(j, 0);
      for (std::size_t k = 1; k < x.ncols(); ++k)
        *it += x(i, k) * y(j, k);
    }
  } else if (transpose_x && transpose_y) {
    for (auto it = out.begin(); it != out.end(); ++it) {
      auto i = it.iRow();
      auto j = it.jCol();
      *it = x(0, i) * y(j, 0);
      for (std::size_t k = 1; k < x.nrows(); ++k)
        *it += x(k, i) * y(j, k);
    }
  }
  }
  return out;
}

template <typename T, typename S>
auto quadraticForm_Symmetric(const M_Matrix<T> &A, const M_Matrix<S> &B)
    -> M_Matrix<decltype(std::declval<T>() * std::declval<S>())> {
  typedef decltype(std::declval<T>() * std::declval<S>()) R;

  auto AB = A * B;
  M_Matrix<R> out(B.ncols(), B.ncols(), Matrix_TYPE::SYMMETRIC, R{});
  for (std::size_t i = 0; i < out.nrows(); ++i)
    for (std::size_t j = i; j < out.ncols(); ++j)
      for (std::size_t k = 0; k < B.nrows(); ++k)
        out(i, j) += B(k, i) * AB(k, j);
  return out;
}

template <typename T, typename S>
auto quadraticForm_Diagonal(const M_Matrix<T> &A, const M_Matrix<S> &B)
    -> M_Matrix<decltype(std::declval<T>() * std::declval<S>())> {
  typedef decltype(std::declval<T>() * std::declval<S>()) R;

  M_Matrix<R> out(B.ncols(), B.ncols(), Matrix_TYPE::SYMMETRIC, R{});
  for (std::size_t i = 0; i < out.nrows(); ++i)
    for (std::size_t j = i; j < out.ncols(); ++j)
      for (std::size_t k = 0; k < B.nrows(); ++k)
        out(i, j) += B(k, i) * A(k, k) * B(k, j);
  return out;
}

template <typename T, typename S>
auto quadraticForm_Diagonal_T(const M_Matrix<T> &A, const M_Matrix<S> &B)
    -> M_Matrix<decltype(std::declval<T>() * std::declval<S>())> {
  typedef decltype(std::declval<T>() * std::declval<S>()) R;

  M_Matrix<R> out(B.nrows(), B.nrows(), Matrix_TYPE::SYMMETRIC, R{});
  for (std::size_t i = 0; i < out.nrows(); ++i)
    for (std::size_t j = i; j < out.ncols(); ++j)
      for (std::size_t k = 0; k < B.ncols(); ++k)
        out(i, j) += B(i, k) * A(k, k) * B(j, k);
  return out;
}

template <typename T, typename S>
auto quadraticForm_Symmetric_T(const M_Matrix<T> &A, const M_Matrix<S> &B)
    -> M_Matrix<decltype(std::declval<T>() * std::declval<S>())> {
  typedef decltype(std::declval<T>() * std::declval<S>()) R;

  auto BA = B * A;
  M_Matrix<R> out(B.nrows(), B.nrows(), Matrix_TYPE::SYMMETRIC, R{});
  for (std::size_t i = 0; i < out.nrows(); ++i)
    for (std::size_t j = i; j < out.ncols(); ++j)
      for (std::size_t k = 0; k < B.ncols(); ++k)
        out(i, j) += BA(i, k) * B(j, k);
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
      if (one.type()==Matrix_TYPE::FULL)
        {
          if (other.type()==Matrix_TYPE::FULL)
            return Lapack_Full_Product
                (one,other,transposeOne,transposeOther);
          else if (other.type()==Matrix_TYPE::SYMMETRIC)
            {
              if (transposeOne)
                return Lapack_Transpose_Symmetric_Product(one,other);
              else {

                  auto tr=Transpose(one);
                  return Lapack_Symmetric_Transpose_Product(tr,other);
                }
            }
          else if ((other.type()==Matrix_TYPE::DIAGONAL)
                   || (other.type()==Matrix_TYPE::SCALAR_DIAGONAL))
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
      else   if (one.type()==Matrix_TYPE::SYMMETRIC)
        {
          if (other.type()==Matrix_TYPE::FULL)
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
          else if (other.type()==Matrix_TYPE::SYMMETRIC)
            {
              auto full=other.full();
              return Lapack_Symmetric_Transpose_Product(one,full,true);
            }
          else if ((other.type()==Matrix_TYPE::DIAGONAL)
                   || (other.type()==Matrix_TYPE::SCALAR_DIAGONAL))
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
      else  if (one.type()==Matrix_TYPE::DIAGONAL)
        {
          if (other.type()==Matrix_TYPE::FULL)
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
          else if (other.type()==Matrix_TYPE::SYMMETRIC)
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
          else if ((other.type()==Matrix_TYPE::DIAGONAL)
                   || (other.type()==Matrix_TYPE::SCALAR_DIAGONAL))
            {
              M_Matrix<double> out
                  (Matrix_TYPE::DIAGONAL,one.nrows(),other.ncols());
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

      else  if (one.type()==Matrix_TYPE::SCALAR_DIAGONAL)
        {
          if (other.type()==Matrix_TYPE::FULL)
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
          else if (other.type()==Matrix_TYPE::SYMMETRIC)
            {
              M_Matrix<double> out
                  (one.nrows(),other.ncols(),Matrix_TYPE::SYMMETRIC);
              for (std::size_t i=0; i<other.nrows(); ++i)
                {
                  for (std::size_t j=0; j<=i; ++j)
                    out(i,j)=one(i,i)*other(i,j);
                }
              return out;
            }
          else if (other.type()==Matrix_TYPE::DIAGONAL)

            {
              M_Matrix<double> out
                  (one.nrows(),other.ncols(),Matrix_TYPE::DIAGONAL);
              for (std::size_t i=0; i<out.nrows(); ++i)
                {
                  out(i,i)=one(i,i)*other(i,i);
                }
              return out;
            }
          else if ((other.type()==Matrix_TYPE::SCALAR_DIAGONAL))
            {
              return M_Matrix<double>
                  (one.nrows(),other.ncols(),Matrix_TYPE::SCALAR_DIAGONAL,
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
          if (other.type()==Matrix_TYPE::FULL)
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
          else if (other.type()==Matrix_TYPE::SYMMETRIC)
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
          else if (other.type()==Matrix_TYPE::DIAGONAL)
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
          else if ((other.type()==Matrix_TYPE::SCALAR_DIAGONAL))
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

inline M_Matrix<double> Matrix_Multiplication(const M_Matrix<double> &one,
                                              const M_Matrix<double> &other,
                                              bool transpose_x,
                                              bool transpose_y)

{
  if ((one.type() == Matrix_TYPE::FULL) && (other.type() == Matrix_TYPE::FULL))
    return Lapack_Full_Product(one, other, transpose_x, transpose_y);

  else if ((one.type() == Matrix_TYPE::FULL) &&
           (other.type() == Matrix_TYPE::SYMMETRIC)) {
    if (transpose_x) {
      auto tr = Transpose(one);
      return Lapack_Symmetric_Regular_Product(other, tr, false);
    } else {
      return Lapack_Symmetric_Regular_Product(other, one, false);
    }
  } else if ((one.type() == Matrix_TYPE::SYMMETRIC) &&
             (other.type() == Matrix_TYPE::FULL)) {
    if (transpose_y) {
      auto tr = Transpose(other);

      return Lapack_Symmetric_Regular_Product(one, tr, true);
    } else {
      return Lapack_Symmetric_Regular_Product(one, other, true);
    }
  } else if ((one.type() == Matrix_TYPE::SYMMETRIC) &&
             (other.type() == Matrix_TYPE::SYMMETRIC)) {
    auto copy_other = M_Matrix<double>(other.nrows(), other.ncols(),
                                       Matrix_TYPE::FULL, other);
    return Lapack_Symmetric_Regular_Product(one, copy_other, true);
  } else
    return All_Product(one, other, transpose_x, transpose_y);
}

template <typename T, typename S>
auto Matrix_Multiplication(const M_Matrix<T> &one, const M_Matrix<S> &other,
                           bool transposeOne, bool transposeOther)
    -> M_Matrix<decltype(std::declval<T>() * std::declval<S>())> {
  typedef decltype(std::declval<T>() * std::declval<S>()) R;
  M_Matrix<R> out = All_Product(one, other, transposeOne, transposeOther);
  assert((are_Equal<true, M_Matrix<R>>().test_sum(
      out, UnformattedProduct(one, other), std::cerr)));
  return out;
}

} // namespace product

template <typename T, typename S>
auto operator*(const M_Matrix<T> &one, const M_Matrix<S> &other) -> M_Matrix<
    std::enable_if_t<is_Matrix_v<T> == is_Matrix_v<S>,
                     decltype(std::declval<T>() * std::declval<S>())>> {

  typedef decltype(std::declval<T>() * std::declval<S>()) R;

  M_Matrix<R> out = product::Matrix_Multiplication(one, other, false, false);
  // assert(apply_test(are_Equal<true,R>(std::numeric_limits<R>::epsilon(),std::sqrt(std::numeric_limits<R>::epsilon())
  // ),out,product::UnformattedProduct(one,other),std::cerr));
  return out;
}

template <typename T, typename S>
auto multTransp(const M_Matrix<T> &one, const M_Matrix<S> &other)
    -> M_Matrix<decltype(std::declval<T>() * std::declval<S>())> {
  return product::Matrix_Multiplication(one, other, false, true);
}

/**
     Transpose the first and multiply by the second
     @post transpMult(x,y)==Transpose(x)*y
     @remarks It is faster, since we save copying matrices
    */
template <typename T, typename S>
auto TranspMult(const M_Matrix<T> &one, const M_Matrix<S> &other)
    -> M_Matrix<decltype(std::declval<T>() * std::declval<S>())> {
  return product::Matrix_Multiplication(one, other, true, false);
}

inline M_Matrix<double> Forward_Sustitution_Ly_b(const M_Matrix<double> &L,
                                                 const M_Matrix<double> &b) {
  assert(L.nrows() == b.nrows());
  M_Matrix<double> y(L.ncols(), b.ncols());
  for (std::size_t k = 0; k < b.ncols(); ++k) {
    y(0, k) = b(0, k) / L(0, 0);
  }
  for (std::size_t i = 0; i < L.nrows(); i++) {
    for (std::size_t k = 0; k < b.ncols(); ++k) {
      auto bb = L(i, 0) * y(0, k);
      for (std::size_t j = 1; j < i; ++j)
        bb += L(i, j) * y(j, k);
      bb = b(i, k) - bb;
      y(i, k) = bb / L(i, i);
    }
  }
  return y;
}

template <typename T>
M_Matrix<T> Forward_Sustitution_Ly_b(const M_Matrix<T> &L,
                                     const M_Matrix<T> &b) {
  assert(L.nrows() == b.nrows());
  M_Matrix<T> y(L.ncols(), b.ncols());
  for (std::size_t k = 0; k < b.ncols(); ++k) {
    y(0, k) = Forward_Sustitution_Ly_b(L(0, 0), b(0, k));
  }
  for (std::size_t i = 0; i < L.nrows(); i++) {
    for (std::size_t k = 0; k < b.ncols(); ++k) {
      auto bb = L(i, 0) * y(0, k);
      for (std::size_t j = 1; j < i; ++j)
        bb += L(i, j) * y(j, k);
      bb = b(i, k) - bb;
      y(i, k) = Forward_Sustitution_Ly_b(L(i, i), bb);
    }
  }
  return y;
}

inline M_Matrix<double> Back_Sustitution_Ux_b(const M_Matrix<double> &U,
                                              const M_Matrix<double> &b) {
  assert(U.nrows() == b.nrows());
  M_Matrix<double> x(U.ncols(), b.ncols());
  std::size_t n = U.nrows();
  for (std::size_t k = 0; k < b.ncols(); ++k) {
    x(n - 1, k) = b(n - 1, k) / U(n - 1, n - 1);
  }
  for (std::size_t i = 1; i < U.nrows(); i++) {
    for (std::size_t k = 0; k < b.ncols(); ++k) {
      auto bb = U(n - 1 - i, n - 1) * x(n - 1, k);
      for (std::size_t j = 1; j < i; ++j)
        bb += U(n - 1 - i, n - 1 - j) * x(n - 1 - j, k);
      bb = b(n - 1 - i, k) - bb;
      x(n - 1 - i, k) = bb / U(n - 1 - i, n - 1 - i);
    }
  }
  return x;
}

template <typename T>
M_Matrix<T> Back_Sustitution_Ux_b(const M_Matrix<T> &U, const M_Matrix<T> &b) {
  assert(U.nrows() == b.nrows());
  M_Matrix<T> x(U.ncols(), b.ncols());
  std::size_t n = U.nrows();
  for (std::size_t k = 0; k < b.ncols(); ++k) {
    x(n - 1, k) = Back_Sustitution_Ux_b(U(n - 1, n - 1), b(n - 1, k));
  }
  for (std::size_t i = 1; i < U.nrows(); i++) {
    for (std::size_t k = 0; k < b.ncols(); ++k) {
      auto bb = U(n - 1 - i, n - 1) * x(n - 1, k);
      for (std::size_t j = 1; j < i; ++j)
        bb += U(n - 1 - i, n - 1 - j) * x(n - 1 - j, k);
      bb = b(n - 1 - i, k) - bb;
      x(n - 1 - i, k) = Back_Sustitution_Ux_b(U(n - 1 - i, n - 1 - i), bb);
    }
  }
  return x;
}

template <typename T, typename S>
auto quadraticForm_BT_A_B(const M_Matrix<T> &A, const M_Matrix<S> &B)
    -> M_Matrix<decltype(std::declval<T>() * std::declval<S>())> {
  assert(A.ncols() == B.nrows());
  assert(B.nrows() == A.nrows());
  typedef decltype(std::declval<T>() * std::declval<S>()) R;
  if ((A.type() == Matrix_TYPE::ZERO) || (B.type() == Matrix_TYPE::ZERO))
    return M_Matrix<R>(B.ncols(), B.ncols(), Matrix_TYPE::ZERO);
  else if (A.isDiagonal())
    return product::quadraticForm_Diagonal(A, B);
  else if (A.isSymmetric())
    return product::quadraticForm_Symmetric(A, B);
  else
    return TranspMult(B, A * B);
}




template <typename T, typename S>
auto quadraticForm_B_A_BT(const M_Matrix<T> &A, const M_Matrix<S> &B)
    -> M_Matrix<decltype(std::declval<T>() * std::declval<S>())> {
  assert(A.ncols() == B.ncols());
  assert(B.ncols() == A.nrows());
  typedef decltype(std::declval<T>() * std::declval<S>()) R;
  if ((A.type() == Matrix_TYPE::ZERO) || (B.type() == Matrix_TYPE::ZERO))
    return M_Matrix<R>(B.nrows(), B.nrows(), Matrix_TYPE::ZERO);
  else if (A.isDiagonal())
    return product::quadraticForm_Diagonal_T(A, B);
  else if (A.isSymmetric())
    return product::quadraticForm_Symmetric_T(A, B);
  else
    return multTransp(B * A, B);
}

template <typename T, typename S>
auto mean(const M_Matrix<T> &A, const M_Matrix<S> &B)
    -> M_Matrix<decltype(std::declval<T>() + std::declval<S>())> {
    assert(A.nrows() == B.nrows());
    assert(B.ncols() == A.ncols());
    return (A+B)*0.5;

 }
 template <class Norm,typename T, typename S>
 auto variance(const M_Matrix<T> &A, const M_Matrix<S> &B)
     -> M_Matrix<decltype(std::declval<T>() + std::declval<S>())> {
     assert(A.nrows() == B.nrows());
     assert(B.ncols() == A.ncols());
     return (A-B).apply([](auto x){return sqr(x)/2.0;});

 }




/**
     Scalar Multiplication assignment.
     @returns a reference to itself
     @post all the values of the matrix are multiplied by the value x
     */
template <typename E, typename T>
auto operator*=(M_Matrix<E> &itself, T x) -> M_Matrix<
    std::enable_if_t<!(is_Matrix_v<T> && !is_Matrix_v<E>),
                     decltype(std::declval<E>() * std::declval<T>())>> & {
  if (itself.type() == Matrix_TYPE::ZERO)
    return itself;
  else {
    for (size_t i = 0; i < itself.size(); i++)
      itself[i] *= x;
    return itself;
  }
}

/**
     Scalar Division assignment.
     @returns a reference to itself
     @post all the values of the matrix are divided by the value x
     */
template <typename E, typename T>
auto operator/=(M_Matrix<E> &itself, T x) -> M_Matrix<
    std::enable_if_t<!(is_Matrix_v<T> && !is_Matrix_v<E>),
                     decltype(std::declval<E>() * std::declval<T>())>> & {
  if (itself.type() != Matrix_TYPE::ZERO)
    for (size_t i = 0; i < itself.size(); i++)
      itself[i] /= x;

  return itself;
}

//@}
/** @name Aritmetic Assigment Operations between two Matrices
 */
//@{

namespace additive {

using namespace Matrix_Unary_Transformations;

template <typename T, typename S>
M_Matrix<T> &sum_assigment(M_Matrix<T> &itself, const M_Matrix<S> &other) {
  auto resType = result_of_sum(itself, other);
  if (resType != itself.type()) {
    itself = M_Matrix<T>(itself.nrows(), itself.ncols(), resType, itself);
  }
  auto itType = iterator_of_sum(itself, other);
  for (auto it = other.begin(itType); it != other.end(); ++it)
    itself(it) += *it;
  return itself;
}

template <typename T, typename S>
M_Matrix<T> &substraction_assigment(M_Matrix<T> &itself,
                                    const M_Matrix<S> &other) {
  auto resType = result_of_sum(itself, other);
  if (resType != itself.type()) {
    itself = M_Matrix<T>(itself.nrows(), itself.ncols(), resType, itself);
  }
  auto itType = iterator_of_sum(itself, other);
  for (auto it = other.begin(itType); it != other.end(); ++it)
    itself(it) -= *it;
  return itself;
}

template <typename T, typename S>
auto sum_Operator(const M_Matrix<T> &itself, const M_Matrix<S> &other)
    -> M_Matrix<decltype(std::declval<T>() + std::declval<S>())> {
  typedef decltype(std::declval<T>() + std::declval<S>()) R;
  if (itself.size() < other.size()) {
    M_Matrix<R> out(other);
    sum_assigment(out, itself);
    return out;
  } else {
    M_Matrix<R> out(itself);
    sum_assigment(out, other);
    return out;
  }
}

template <typename T, typename S>
auto sum_Operator(M_Matrix<T> &&itself, M_Matrix<S> &&other)
    -> M_Matrix<decltype(std::declval<T>() + std::declval<S>())> {
  typedef decltype(std::declval<T>() + std::declval<S>()) R;
  if (itself.size() < other.size()) {
    M_Matrix<R> out(other);
    sum_assigment(out, itself);
    return out;
  } else {
    M_Matrix<R> out(itself);
    sum_assigment(out, other);
    return out;
  }
}

template <typename T, typename S>
auto substraction_Operator(const M_Matrix<T> &itself, const M_Matrix<S> &other)
    -> M_Matrix<decltype(std::declval<T>() + std::declval<S>())> {
  typedef decltype(std::declval<T>() + std::declval<S>()) R;
  M_Matrix<R> out(itself);
  substraction_assigment(out, other);
  return out;
}
template <typename T, typename S>
auto substraction_Operator(M_Matrix<T> &&itself, M_Matrix<S> &&other)
    -> M_Matrix<decltype(std::declval<T>() + std::declval<S>())> {
  typedef decltype(std::declval<T>() + std::declval<S>()) R;
  M_Matrix<R> out(std::move(itself));
  substraction_assigment(out, other);
  return out;
}

} // namespace additive

template <typename T, typename S, class AsOp>
M_Matrix<T> &
Element_Wise_Multiplicative_Assigment_Operator(AsOp op, M_Matrix<T> &itself,
                                               const M_Matrix<S> &other) {
  if (itself.type() == Matrix_TYPE::ZERO) {
    return itself;
  } else if (other.type() == Matrix_TYPE::ZERO) {
    itself = other;
    return itself;
  } else if (itself.type() == other.type()) {
    assert(itself.size() == other.size());
    for (size_t i = 0; i < itself.size(); i++)
      op(itself[i], other[i]);
    return itself;
  } else if ((itself.type() == Matrix_TYPE::DIAGONAL) ||
             (itself.type() == Matrix_TYPE::SCALAR_DIAGONAL)) {
    assert(itself.nrows() == other.nrows());
    assert(itself.ncols() == other.ncols());
    switch (other.type()) {
    case Matrix_TYPE::ZERO: // not possible
    case Matrix_TYPE::FULL:
    case Matrix_TYPE::SYMMETRIC:
    case Matrix_TYPE::SCALAR_FULL:
    case Matrix_TYPE::SCALAR_DIAGONAL:
    case Matrix_TYPE::DIAGONAL: {
      for (std::size_t i = 0; i < itself.nrows(); ++i)
        op(itself(i, i), other(i, i));
      return itself;
    }
    }
  }

  else if ((other.type() == Matrix_TYPE::DIAGONAL) ||
           (other.type() == Matrix_TYPE::SCALAR_DIAGONAL)) {
    assert(itself.nrows() == other.nrows());
    assert(itself.ncols() == other.ncols());
    M_Matrix<T> out(other);
    for (std::size_t i = 0; i < std::min(out.nrows(), out.ncols()); ++i)
      op(out(i, i), itself(i, i));
    itself = std::move(other);
    return itself;

  } else if (itself.type() == Matrix_TYPE::FULL) {
    assert(itself.nrows() == other.nrows());
    assert(itself.ncols() == other.ncols());
    switch (other.type()) {
    case Matrix_TYPE::ZERO:            // not
    case Matrix_TYPE::SCALAR_DIAGONAL: // done
    case Matrix_TYPE::DIAGONAL:        // done
    case Matrix_TYPE::FULL:
    case Matrix_TYPE::SYMMETRIC:
    case Matrix_TYPE::SCALAR_FULL: {
      for (std::size_t i = 0; i < itself.nrows(); ++i)
        for (std::size_t j = 0; j < itself.ncols(); ++j)
          op(itself(i, j), other(i, j));
      return itself;
    }
    }
  } else if (other.type() == Matrix_TYPE::FULL) {
    assert(itself.nrows() == other.nrows());
    assert(itself.ncols() == other.ncols());
    M_Matrix<T> out(other);
    for (std::size_t i = 0; i < itself.nrows(); ++i)
      for (std::size_t j = 0; j < itself.ncols(); ++j)
        op(out(i, j), itself(i, j));
    itself = std::move(out);
    return itself;
  } else if (itself.type() == Matrix_TYPE::SYMMETRIC) {
    assert(itself.nrows() == other.nrows());
    assert(itself.ncols() == other.ncols());
    switch (other.type()) {
    case Matrix_TYPE::ZERO:            // not possible
    case Matrix_TYPE::SCALAR_DIAGONAL: // not possible
    case Matrix_TYPE::DIAGONAL:        // not possible
    case Matrix_TYPE::FULL:            // should not land here
    case Matrix_TYPE::SYMMETRIC:       // should not land here
    case Matrix_TYPE::SCALAR_FULL: {
      for (std::size_t i = 0; i < itself.nrows(); ++i)
        for (std::size_t j = 0; j <= i; ++j)
          op(itself(i, j), other(i, j));
      return itself;
    }
    }
  } else // scalar_full and the other is symmetric
  {
    assert(itself.nrows() == other.nrows());
    assert(itself.ncols() == other.ncols());
    M_Matrix<T> out(other);
    for (std::size_t i = 0; i < itself.nrows(); ++i)
      for (std::size_t j = 0; j <= i; ++j)
        op(out(i, j), itself(i, j));
    itself = std::move(out);
    return itself;
  }
  return itself; // lo pongo para que compile, no se porque piensa que llega
                 // aqui...
}

//@}

/**
 Matrix sum, element wise.
 @pre \p x.nrows()==ncols(y) x.ncols()==ncols(y)
 @return z, where z.nrows()=rows(x), z.ncols()=ncols(y) and z(i,j)= sum
 on k of x(i,k)*y(k,j)
 @warning it \c assert the preconditions
 */
template <typename T>
M_Matrix<T> operator+(const M_Matrix<T> &x, const M_Matrix<T> &y) {
  if (x.type() == Matrix_TYPE::ZERO)
    return y;
  else if (y.type() == Matrix_TYPE::ZERO)
    return x;
  else {
    auto out = additive::sum_Operator(x, y);
    return out;
  }
}

template <typename T, typename S>
M_Matrix<T> &operator+=(M_Matrix<T> &x, const M_Matrix<S> &y) {
  if (x.type() == Matrix_TYPE::ZERO) {
    x = y;
    return x;
  } else if (y.type() == Matrix_TYPE::ZERO)
    return x;
  else
    return additive::sum_assigment(x, y);
}

template <typename T, typename S>
M_Matrix<T> &operator-=(M_Matrix<T> &x, const M_Matrix<S> &y) {
  if (x.type() == Matrix_TYPE::ZERO) {
    x = -y;
    return x;
  } else if (y.type() == Matrix_TYPE::ZERO)
    return x;
  else
    return additive::substraction_assigment(x, y);
}

/**
 Matrix sustraction, element wise.
 @pre \p x.nrows()==ncols(y) x.ncols()==ncols(y)
 @return z, where z.nrows()=rows(x), z.ncols()=ncols(y) and
 z(i,j)= sum on k of x(i,k)*y(k,j)
 @warning it \c assert the preconditions
 */

template <typename T>
M_Matrix<T> operator-(const M_Matrix<T> &x, const M_Matrix<T> &y) {
  if (x.type() == Matrix_TYPE::ZERO) {
    return y;
  } else if (y.type() == Matrix_TYPE::ZERO)
    return x;
  else

    return additive::substraction_Operator(x, y);
}

template <typename T> M_Matrix<T> operator-(M_Matrix<T> &&x, M_Matrix<T> &&y) {
  if (x.type() == Matrix_TYPE::ZERO) {
    return -y;
  } else if (y.type() == Matrix_TYPE::ZERO)
      return std::move(x);
  else

    return additive::substraction_Operator(std::move(x), std::move(y));
}

/** @name Aritmetic operation between a Matrix and a scalar (single element)
 */
//@{

/**
     Scalar Addition.
     @returns a copy of the matrix with its values summed by x
     */
template <typename T>
M_Matrix<T>& operator+=(M_Matrix<T> &x,
                      T t) { // we build the M_Matrix result
    for (std::size_t i=0; i<x.size(); ++i)
    x[i] += t;
    return x;
}

template <typename T>
M_Matrix<T>& operator-=(M_Matrix<T> &x,
                        T t) { // we build the M_Matrix result
    for (std::size_t i=0; i<x.size(); ++i)
        x[i] -= t;
    return x;
}


/**
     Scalar Addition.
     @returns a copy of the matrix with its values summed by x
     */
template <typename T>
M_Matrix<T> operator+(const M_Matrix<T> &x,
                      T t) { // we build the M_Matrix result
  M_Matrix<T> z(x);
  z += t;
  return z;
}

/**
     Scalar Addition reverse order.
     */
template <typename T> M_Matrix<T> operator+(T t, const M_Matrix<T> &x) {
  return x + t;
}

/**
     Scalar Subtraction.
     @returns a copy of the matrix with its values substracted by x
     */
template <typename T>
M_Matrix<T> operator-(const M_Matrix<T> &x,
                      T t) { // we build the M_Matrix result
  M_Matrix<T> z(x);
  z -= t;
  return z;
}

/**
     Scalar Subtraction reverse order.
     */
template <typename T> M_Matrix<T> operator-(T t, const M_Matrix<T> &x) {
  return x - t;
}

/**
     Scalar Multiplication.
     @returns a copy of the matrix with its values multiplied by the value x
     */
template <typename E, typename T>
auto operator*(const M_Matrix<E> &x, T t)
    -> std::enable_if_t<!(is_Matrix_v<T> && !is_Matrix_v<E>)&&std::is_same_v<
                            E, decltype(std::declval<E>() * std::declval<T>())>,
                        M_Matrix<E>> { // we build the M_Matrix result
  typedef decltype(std::declval<E>() * std::declval<T>()) R;
  M_Matrix<R> out(x);
  for (std::size_t i = 0; i < x.size(); ++i)
    out[i] = out[i] * t;
  return out;
}

template <typename E, typename T>
auto operator*(T t, const M_Matrix<E> &x)
    -> std::enable_if_t<!(is_Matrix_v<T> && !is_Matrix_v<E>)&&std::is_same_v<
                            E, decltype(std::declval<E>() * std::declval<T>())>,
                        M_Matrix<E>> { // we build the M_Matrix result
  typedef decltype(std::declval<E>() * std::declval<T>()) R;
  M_Matrix<R> out(x);
  for (std::size_t i = 0; i < x.size(); ++i)
    out[i] = t * out[i];
  return out;
}

/**
     Scalar Division.
     @returns a copy of the matrix with its values divided by x
     @returns a matrix of real numbers
 */
template <typename T>
M_Matrix<double> operator/(const M_Matrix<T> &x,
                           T t) { // we build the M_Matrix result
  M_Matrix<double> z(x);
  z /= double(t);
  return z;
}

/**
     Division by inhomogeneus types

     */

template <typename T, typename S>
M_Matrix<double> operator/(const M_Matrix<T> &x,
                           S t) { // we build the M_Matrix result
  M_Matrix<double> z(x);
  z /= double(t);
  return z;
}

/**
     Scalar Division reverse order.
     */
template <typename T> M_Matrix<double> operator/(T t, const M_Matrix<T> &x) {
  M_Matrix<double> out(x.nrows(), x.ncols());
  for (std::size_t i = 0; i < x.size(); i++)
    out[i] = double(t) / double(x[i]);

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
template <typename T, typename S,
          std::enable_if_t<std::is_convertible_v<S, T>, int> = 0>
M_Matrix<T> elemMult(const M_Matrix<T> &x, const M_Matrix<S> &y) {
  auto out(x);
  return Element_Wise_Multiplicative_Assigment_Operator(
      [](T &a, const S &b) {
        a *= b;
        return a;
      },
      out, y);
}

template <typename T>
M_Matrix<M_Matrix<T>> elemMult_a(const M_Matrix<M_Matrix<T>> &x,
                                 const M_Matrix<T> &y) {
  return x.apply([&y](auto &e) { return elemMult(e, y); });
}
template <typename T>
M_Matrix<M_Matrix<T>> elemMult_a(const M_Matrix<T> &x,
                                 const M_Matrix<M_Matrix<T>> &y) {
  return y.apply([&x](auto &e) { return elemMult(x, e); });
}

template <typename T, typename S,
          std::enable_if_t<std::is_convertible_v<S, T>, int> = 0>
M_Matrix<T> elemMult(M_Matrix<T> &&x, M_Matrix<S> &&y) {
  return Element_Wise_Multiplicative_Assigment_Operator(
      [](T &a, const S &b) {
        a *= b;
        return a;
      },
      x, y);
}

template <typename T, typename S>
M_Matrix<T> elemDiv(const M_Matrix<T> &x, const M_Matrix<S> &y) {
  auto out(x);
  return Element_Wise_Multiplicative_Assigment_Operator(
      [](T &a, const S &b) {
        a /= b;
        return a;
      },
      out, y);
}

template <typename T, typename S>
M_Matrix<T> elemDiv(M_Matrix<T> &&x, M_Matrix<S> &&y) {
  return Element_Wise_Multiplicative_Assigment_Operator(
      [](T &a, const S &b) {
        a /= b;
        return a;
      },
      x, y);
}

inline double elemDivSafe(double x, double y,
                   double eps = std::numeric_limits<double>::epsilon()) {
  if (std::abs(y) > eps)
    return x / y;
  else
    return 0;
}

template <typename T, typename S>
auto UnformatedelemDivSafe(const M_Matrix<T> &x, const M_Matrix<S> &y,
                           double eps = std::numeric_limits<double>::epsilon())
    -> M_Matrix<decltype(std::declval<T>() * std::declval<S>())> {
  typedef decltype(std::declval<T>() * std::declval<S>()) R;
  assert(x.nrows() == y.nrows());
  assert(x.ncols() == y.ncols());
  M_Matrix<R> out(x.nrows(), x.ncols(), 0);
  for (std::size_t i = 0; i < x.nrows(); ++i)
    for (std::size_t j = 0; j < x.ncols(); ++j)
      out(i, j) = elemDivSafe(x(i, j), y(i, j), eps);
  return out;
}

template <typename T, typename S>
auto elemDivSafe(const M_Matrix<T> &x, const M_Matrix<S> &y,
                 double eps = std::numeric_limits<double>::epsilon()) {
  typedef decltype(std::declval<T>() * std::declval<S>()) R;
  //double norm = Matrix_Unary_Functions::norm_1(y);
  M_Matrix<R> out(x);
  out = Element_Wise_Multiplicative_Assigment_Operator(
      [eps](T &a, const S &b) {
        return a = elemDivSafe(a, b, eps);
        return a;
      },
      out, y);
  // assert((are_Equal<true,
  // M_Matrix<R>>().test_sum(out,UnformatedelemDivSafe(x,y,eps), std::cerr)));
  return out;
}

template <typename T, typename S>
M_Matrix<T> elemDivSafe(M_Matrix<T> &&x, M_Matrix<S> &&y,
                        double eps = std::numeric_limits<double>::epsilon()) {
  double norm = Matrix_Unary_Functions::norm_1(y);
  return Element_Wise_Multiplicative_Assigment_Operator(
      [eps, norm](T &a, const S &b) { return elemDivSafe(a, b, eps * norm); },
      x, y);
}

template <typename T>
auto elemDivSafe_a(const M_Matrix<M_Matrix<T>> &x, const M_Matrix<T> &y,
                   double eps = std::numeric_limits<double>::epsilon()) {
  return x.apply([&y, &eps](auto &e) { return elemDivSafe(e, y, eps); });
}

template <typename T>
auto elemDiv_a(const M_Matrix<M_Matrix<T>> &x, const M_Matrix<T> &y) {
  return x.apply([&y](auto &e) { return elemDiv(e, y); });
}

/**
 Division of the elements of two matrices.
  @pre \p x.nrows()==ncols(y) x.ncols()==ncols(y)
  @return z, where z.nrows()=rows(x), z.ncols()=ncols(y)
 and z(i,j)= sum on k of x(i,k)*y(k,j)
  @warning it \c assert the preconditions
 */
/*template<typename T,typename S>
M_Matrix<double> elemDiv(const M_Matrix<double>& x,const M_Matrix<S>& y)
{
    M_Matrix<double> z(x);
    for (size_t i=0; i<z.nrows(); i++)
        for (size_t j=0; j<z.ncols(); j++)
            z(i,j)/=double(y(i,j));
    return z;
}
*/
template <typename T>
M_Matrix<double> operator<<(const M_Matrix<double> &A, const M_Matrix<T> &B) {
  M_Matrix<double> out(std::max(A.nrows(), B.nrows()), A.ncols() + B.ncols(),
                       std::numeric_limits<double>::quiet_NaN());
  for (std::size_t i = 0; i < A.nrows(); ++i) {
    for (std::size_t j = 0; j < A.ncols(); ++j)
      out(i, j) = A(i, j);
  }
  for (std::size_t i = 0; i < B.nrows(); ++i) {
    for (std::size_t j = 0; j < B.ncols(); ++j)
      out(i, j) = B(i, A.ncols() + j);
  }

  return out;
}

inline M_Matrix<double> xdiagXT(const M_Matrix<double> &x,
                                const M_Matrix<double>& Cdiag) {
  M_Matrix<double> o(x.nrows(), x.nrows(), 0.0);
  for (std::size_t i = 0; i < x.nrows(); ++i)
    for (std::size_t j = 0; j < x.nrows(); ++j)
      for (std::size_t k = 0; k < x.ncols(); ++k)
        o(i, j) += Cdiag[k] * x(i, k) * x(j, k);
  return o;
}

inline M_Matrix<double> MultDiag(const M_Matrix<double> &x,
                                 const M_Matrix<double>& d) {
  M_Matrix<double> o(x.nrows(), x.ncols());
  for (std::size_t i = 0; i < x.nrows(); ++i)
    for (std::size_t j = 0; j < x.ncols(); ++j)
      o(i, j) = x(i, j) * d[j];
  return o;
}

inline M_Matrix<double> DiagMult(const M_Matrix<double>& d,
                                 const M_Matrix<double> &x) {
  M_Matrix<double> o(x.nrows(), x.ncols());
  for (std::size_t i = 0; i < x.nrows(); ++i)
    for (std::size_t j = 0; j < x.ncols(); ++j)
      o(i, j) = x(i, j) * d[i];
  return o;
}

template <typename T>
M_Matrix<T> diag_landa(const M_Matrix<T> &x, double landa) {
  double landa1 = landa + 1;
  M_Matrix<T> diagM(x);
  for (size_t i = 0; i < x.nrows(); ++i)
    diagM(i, i) *= landa1;
  return diagM;
}

template <typename V>
M_Matrix<double> JTd2J(const M_Matrix<double> &J, const V &D2) {
  assert(J.nrows() == D2.size());
  M_Matrix<double> out(J.ncols(), J.ncols());
  M_Matrix<double> Jc = Transpose(J);
  for (std::size_t i = 0; i < Jc.nrows(); ++i)
    for (std::size_t j = 0; j < Jc.ncols(); ++j)
      Jc(i, j) *= D2[j];
  out = Jc * J;
  return out;
}

} // namespace Matrix_Binary_Transformations

using namespace Vector_Binary_Transformations;
using namespace Container_Unary_Functions;
using namespace Container_Binary_Transformations;
using namespace Matrix_Generators;
// using namespace Matrix_Stream_Operators;
using namespace Matrix_Unary_Predicates;
using namespace Matrix_Unary_Size_Functions;
using namespace Matrix_Unary_Functions;
using namespace Matrix_Unary_Transformations;
using namespace Matrix_Binary_Predicates;
using namespace Matrix_Binary_Transformations;

//@}

template <class C> double sum(const C &x) {
  double o = 0;
  for (unsigned i = 0; i < x.size(); ++i)
    o += x[i];

  return o;
}

template <typename T> double Frobenius_norm(const M_Matrix<T> &x) {
  return std::sqrt(fullSum(elemMult(x, x)));
}

struct Frobenius_test : public invariant {
  template <typename T>
  static bool test(const M_Matrix<T> &one, const M_Matrix<T> &other,
                   double conditionNumber,
                   double eps = std::numeric_limits<double>::epsilon()) {
    double norm = Frobenius_norm(one);
    double normtest = Frobenius_norm(one - other) / norm * conditionNumber;
    if (normtest < eps)
      return true;
    else
      return false;
  }
};

inline
    double expm_taylor_error(const M_Matrix<double>& x, std::size_t order)
{
    double r=Matrix_Unary_Functions::norm_inf(x);
    double error=std::exp(r)*std::pow(r,order)/std::tgamma(order+1);
    return error;

}


template <typename T>
double calc_Commutator_Taylor_Series_error(const M_Matrix<T>  & Qx,const M_Matrix<T>  & GP, std::int8_t order)
{
    double r=Matrix_Unary_Functions::norm_inf(center(Qx))*2.0;
    double g=Matrix_Unary_Functions::norm_inf(center(GP));
    double error=std::exp(r)*std::pow(r,order)/std::tgamma(order+1)*g;
    return error;

}

inline M_Matrix<double> calc_Commutator_Taylor_Series(const M_Matrix<double>& GP, const M_Matrix<double> & Qx,std::int8_t order)
{
    M_Matrix<double> dGn = GP;
    auto Gtot = dGn;
    double a = 1.0;
    for (std::size_t n = 2; n + 1 < order; ++n) {
        a /= n;
        dGn = Qx * dGn - dGn * Qx;
        Gtot += dGn * a;
    }
    return Gtot;

}

inline std::vector<double> &pascal_triangle(std::vector<double> &coeff) {

    coeff.push_back(1.0);
    double a = 1.0;
    for (std::size_t i = 1; i < coeff.size(); ++i) {
        std::swap(coeff[i - 1], a);
        a = a + coeff[i];
    }
    return coeff;
}

template <class E>
void next_derivative(std::vector<M_Matrix<E>> &dGvar_n, M_Matrix<E> &Tau_n,
                     const M_Matrix<E> &Q, const M_Matrix<E> &GP) {
    M_Matrix<E> neg_2Q(Q * (-2.0));
    for (auto &e : dGvar_n)
        e = e * neg_2Q;
    Tau_n = Q * Tau_n + Tau_n * Q;
    dGvar_n.push_back(Tau_n * GP);
}

template <typename T>
double calc_square_Commutator_Taylor_Series_error(const M_Matrix<T>  & Qx,const M_Matrix<T>  & G, std::int8_t order)
{
    double q=Matrix_Unary_Functions::norm_inf(center(Qx));
    double g=Matrix_Unary_Functions::norm_inf(center(G));
    double r=4*q;
    double error=std::exp(r)*std::pow(r,order)/std::tgamma(order+1)*sqr(g)/sqr(q)/2;
    return error;

}

template <typename E>
M_Matrix<E> calc_square_Commutator_Taylor_Series(const M_Matrix<E> &Qx, const M_Matrix<E> &P,
                                               const M_Matrix<E> &G, std::size_t order) {
    M_Matrix<E> Tau_n = G;
    auto GP = G * P;
    std::vector<double> coeff(1, 1.0);
    std::vector<M_Matrix<E>> dGvar_n(1, G * GP);
    auto Gvar_tot = dGvar_n[0];
    double a = 2.0;
    for (std::size_t n = 3; n + 1 < order; ++n) {
        a /= n;
        pascal_triangle(coeff);
        next_derivative(dGvar_n, Tau_n, Qx, GP);
        for (std::size_t i = 0; i < coeff.size(); ++i)
            Gvar_tot += dGvar_n[i] * (coeff[i] * a);
    }
    return Gvar_tot;
}


template <typename T> inline T sqr(const T &x) { return x * x; }

#endif // MATRIX_H
