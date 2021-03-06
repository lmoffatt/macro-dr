#ifndef MATRIXDERIVATIVE_H
#define MATRIXDERIVATIVE_H
#include "Matrix.h"
#include "myDistributions.h" //for indexes
#include "myfields.h"

struct s_f {
  constexpr static auto const title = my_static_string("_f");
};

struct s_Derivative {  constexpr static auto const title = my_static_string("_der");};

struct transformable_tag{};
struct regular_tag{};
struct derivative_tag{};

template <typename...> struct Der;

template <typename T, typename = void>
struct is_transformable : public std::false_type {};

template<template<template<class...>class,typename...> class myClass,template<class...>class Tr,typename...Ts>
struct is_transformable<
    myClass<Tr,Ts...>, std::void_t<myClass<Tr,Ts...>>>
    : public std::true_type {};

template<template<template<class...>class,auto,auto...> class myClass,template<class...>class Tr,auto x,auto...xs>
struct is_transformable<
    myClass<Tr,x,xs...>, std::void_t<myClass<Tr,x,xs...>>>
    : public std::true_type {};


template <typename T>
inline constexpr bool is_transformable_v = is_transformable<T>::value;



static_assert (!is_transformable_v<double> );


template <typename T> struct Der<T,regular_tag> {

    typedef Derivative<T> type;

};

template <typename K, typename T, typename ...X>
struct Der<std::map<K,T,X...>,map_tag>
{
    typedef std::map<K, Derivative_t<T>> type;
};


template<template<template<class...>class,auto,auto...> class myClass,auto x,auto ...Ts>
struct Der<myClass<C,x,Ts...>,transformable_tag>{
    typedef myClass<Der,x,Ts...> type;
};

template<template<template<class...>class,typename...> class myClass,typename ...Ts>
struct Der<myClass<C,Ts...>,transformable_tag>{
    typedef myClass<Der,Ts...> type;
};


template <class X>
struct Der<X, derivative_tag> {
    typedef typename X::Derivative type;
};


template<class T>
struct my_transform_tag
{
    typedef  std::conditional_t<has_derivative_type_v<T>,derivative_tag,
        std::conditional_t<is_map<T>::value,map_tag,
        std::conditional_t<is_transformable_v<T>,transformable_tag,
                                                                     regular_tag>>> type;
};
template <typename...Ts>
using Der_t=typename Der<Ts...>::type;

template <typename T> struct Der<T> {
    typedef Der_t<T,typename my_transform_tag<T>::type> type;

};




template<>
struct my_template_trait<Der>
{
    constexpr static auto className = my_static_string("_derivative");

};





template <> class Derivative<double> {
  double f_;
  M_Matrix<double> const *x_;
  M_Matrix<double> dfdx_;

public:
  typedef double primitive_type;
  constexpr static auto const className =
      my_static_string("Derivative_") + my_trait<primitive_type>::className;
  typedef Derivative<primitive_type> self_type;
  typedef const Derivative<primitive_type> const_self_type;

  static auto get_constructor_fields() {
    double const &(self_type ::*myf)() const = &self_type::f;
    M_Matrix<double> const &(self_type ::*mydf)() const = &self_type::dfdx;
    myf = &self_type::f;
    return std::make_tuple(grammar::field(C<self_type>{}, "f", myf),
                           grammar::field(C<self_type>{}, "x", &self_type::x),
                           grammar::field(C<self_type>{}, "dfdx", mydf));
  }

  template <class Dindex>
  static auto get_data_index_static(Dindex) {
      double const &(self_type ::*myf)() const = &self_type::f;
    return Concatenate_tuple_static(
          std::make_tuple(make_data_static(
              std::tuple<>(), std::make_tuple(F_s(s_f{}, std::move(myf))))),
        make_data_static(
            std::make_tuple(I_s( Dindex{},
                                [](const self_type &s) { return s.x().size(); },
                                0)),
            std::make_tuple(
                F_s(s_Derivative{}, [](const self_type &s, std::size_t i) {
                  return s.dfdx()[i];
                }))));
  }

  double const &f() const { return f_; }
  double &f() { return f_; }
  M_Matrix<double> const &x() const { return *x_; }
  M_Matrix<double> const * x_ptr()const {return x_;}
  M_Matrix<double> const &dfdx() const { return dfdx_; }
  M_Matrix<double> &dfdx() { return dfdx_; }

  Derivative(double fx, const M_Matrix<double> &myx, M_Matrix<double> &&der)
      : f_{fx}, x_{&myx}, dfdx_{std::move(der)} {}
  Derivative(double fx, const M_Matrix<double> &myx,
             const M_Matrix<double> &der)
      : f_{fx}, x_{&myx}, dfdx_{der} {}

  Derivative(double fx, const M_Matrix<double> &myx)
      : f_{fx}, x_{&myx}, dfdx_{myx.nrows(), myx.ncols(), myx.type(),
                                std::isnan(fx) ? fx : 0} {}
  Derivative() = default;
  Derivative(double f): f_{f},x_{nullptr},dfdx_{}{}
  Derivative(const M_Matrix<double> &myx)
      : f_{0}, x_{&myx}, dfdx_{myx.nrows(), myx.ncols(), myx.type(), 0.0} {}
  Derivative &operator+=(const Derivative &other) {
      if (x_==nullptr)
      {
          x_=&other.x();
          dfdx_=M_Matrix<double>(x().nrows(),x().ncols(),x().type(),0.0);
      }
    assert(x() == other.x());
    f_ += other.f();
    dfdx_ += other.dfdx();
    return *this;
  }

  Derivative &operator*=(const Derivative &other) {
      if (x_==nullptr)
      {
          x_=&other.x();
          dfdx_=M_Matrix<double>(x().nrows(),x().ncols(),x().type(),0.0);
      }
      assert(x() == other.x());
      auto df = f() * other.dfdx();
      for (std::size_t i = 0; i < df.size(); ++i)
          df[i] += dfdx()[i] * other.f();
      f_ *= other.f();
      dfdx_ =df;
      return *this;
  }


  Derivative operator-()const &
  {
      Derivative out(*this);
      return -std::move(out);
  }
  Derivative operator-()&&
  {
      f_=-f_;
      dfdx_=-dfdx_;
      return *this;
  }



  Derivative &operator-=(const Derivative &other) {
      if (x_==nullptr)
      {
          x_=&other.x();
          dfdx_=M_Matrix<double>(x().nrows(),x().ncols(),x().type(),0.0);
      }
    assert(x() == other.x());
    f_ -= other.f();
    dfdx_ -= other.dfdx();
    return *this;
  }
};

inline auto& center(const Derivative<double>& x){return x.f();}




inline auto max(const Derivative<double> &one,double two) {
    if (center(one)>=two)
    {
        return one;
    }
    else
    {
        auto out=one;
        out.f()=two;
        return out;
    }
}

inline auto max(const Derivative<double> &one,const Derivative<double> & two) {
    if (center(one)>=center(two))
    {
        return one;
    }
    else
    {
        return two;
    }
}


template <class T>
std::ostream &operator<<(std::ostream &os, const Derivative<T> &x) {
  return io::output_operator_on_Object(os, x) << io::end_of_Object{};
}

// typedef typename
// decltype(std::declval<std::ostream&>()<<std::declval<Derivative<double>>())::kk
// kk;

static_assert(is_field_Object<Derivative<double>>::value, "is derivative");

template <> class Derivative<M_Matrix<double>> {
  M_Matrix<double> f_;
  M_Matrix<double> const *x_;
  M_Matrix<M_Matrix<double>> dfdx_;

public:
  typedef M_Matrix<double> primitive_type;
  constexpr static auto const className =
      my_static_string("Derivative_") + my_trait<primitive_type>::className;
  typedef Derivative<primitive_type> self_type;
  typedef const Derivative<primitive_type> const_self_type;

  static auto get_constructor_fields() {
    M_Matrix<double> const &(self_type ::*myf)() const = &self_type::f;
    M_Matrix<M_Matrix<double>> const &(self_type ::*mydf)() const =
        &self_type::dfdx;
    return std::make_tuple(grammar::field(C<self_type>{}, "f", myf),
                           grammar::field(C<self_type>{}, "x", &self_type::x),
                           grammar::field(C<self_type>{}, "dfdx", mydf));
  }
  template <class Findex, class Dindex>
  static auto get_data_index_static(Findex, Dindex) {
      return std::make_tuple(
        make_data_static(
            std::make_tuple(I_s(
                Findex{}, [](const self_type &s) { return s.f().size(); }, 0)),
            std::make_tuple(
                F_s(s_f{}, [](const self_type &s,
                              std::size_t i) { return s.f()[i]; }))),
        make_data_static(
            std::make_tuple(I_s(Findex{},
                                [](const self_type &s) { return s.f().size(); },
                                0),
                            I_s(Dindex{},
                                [](const self_type &s, std::size_t) {
                                  return s.x().size();
                                },
                                0)),
            std::make_tuple(
                F_s(s_Derivative{},
                    [](const self_type &s, std::size_t i, std::size_t j) {
                      return s.dfdx()[i][j];
                    }))));
  }

  template <class F_index_rows, class F_index_cols, class Dindex>
  static auto get_data_index_static(F_index_rows, F_index_cols,Dindex) {
      return std::make_tuple(
        make_data_static(
            std::make_tuple(
                I_s(F_index_rows{},
                    [](const self_type &s) { return s.f().nrows(); }, 0),
                I_s(F_index_cols{},
                    [](const self_type &s, std::size_t) {
                      return s.f().ncols();
                    },
                    0)),
            std::make_tuple(
                F_s(s_f{}, [](const self_type &s, std::size_t i,
                              std::size_t j) { return s.f()(i, j); }))),
        make_data_static(
            std::make_tuple(
                I_s(F_index_rows{},
                    [](const self_type &s) { return s.f().nrows(); }, 0),
                I_s(F_index_cols{},
                    [](const self_type &s, std::size_t) {
                      return s.f().ncols();
                    },
                    0),
                I_s(Dindex{},
                    [](const self_type &s, std::size_t, std::size_t) {
                      return s.x().size();
                    },
                    0)),
            std::make_tuple(
                F_s(s_Derivative{},
                    [](const self_type &s, std::size_t i, std::size_t j,
                       std::size_t ider) { return s.dfdx()[ider](i, j); }))));
  }

  Derivative &set_by_f_index(std::size_t i, std::size_t j,
                             const Derivative<double> &dfij) {
    assert(x() == dfij.x());
    f()(i, j) = dfij.f();
    for (std::size_t i = 0; i < dfdx().size(); ++i)
      dfdx()[i](i, j) = dfij.dfdx()[i];
    return *this;
  }

  M_Matrix<double> const &f() const { return f_; }
  M_Matrix<double> &f() { return f_; }
  M_Matrix<double> const &x() const { return *x_; }
  M_Matrix<double> const * x_ptr()const {return x_;}
  M_Matrix<M_Matrix<double>> const &dfdx() const { return dfdx_; }
  M_Matrix<M_Matrix<double>> &dfdx() { return dfdx_; }

  auto nrows()const {return f().nrows();}
  auto ncols()const {return f().ncols();}
  auto size()const {return f().size();}
  auto type()const {return f().type();}

  Derivative<double> operator()(std::size_t i, std::size_t j)const {
      M_Matrix<double> df(x().nrows(),x().ncols(), x().type());
      for (std::size_t n=0; n<dfdx().size(); ++n)
          df[n]=dfdx()[n](i,j);
      return Derivative<double>(f()(i,j),x(),std::move(df));
  }


  Derivative<double> operator[](std::size_t i)const {
      M_Matrix<double> df(x().nrows(),x().ncols(), x().type());
      for (std::size_t n=0; n<dfdx().size(); ++n)
          df[n]=dfdx()[n][i];
      return Derivative<double>(f()[i],x(),std::move(df));
  }


  void set(std::size_t i, std::size_t j, double val){
      f()(i,j)=val;
      if (x_!=nullptr)
      {
          for (std::size_t n = 0; n < dfdx().size(); ++n)
              dfdx()[n](i, j) = 0;
      }
  }
  void set(std::size_t i,  double val)
  {
      f()[i]=val;
      if (x_!=nullptr)
      {
          for (std::size_t n = 0; n < dfdx().size(); ++n)
              dfdx()[n][i] = 0;
      }
  }
  void set(std::size_t i,  const Derivative<double>& val){
      if (x_==nullptr)
      {
          x_=&val.x();
          dfdx_=M_Matrix<M_Matrix<double>>(
              val.x().nrows(), val.x().ncols(), val.x().type(),
              M_Matrix<double>(f().nrows(), f().ncols(),
                               f().type(),0.0));

      }
      assert(&x()==&val.x());
      f()[i]=val.f();
      for (std::size_t n = 0; n < dfdx().size(); ++n)
          dfdx()[n][i]=val.dfdx()[n];
  }
  void set(std::size_t i, std::size_t j, const Derivative<double>& val){
      assert((&val.x())!=nullptr);
      if (x_==nullptr)
      {
          x_=&val.x();
          dfdx_=M_Matrix<M_Matrix<double>>(
              val.x().nrows(), val.x().ncols(), val.x().type(),
              M_Matrix<double>(f().nrows(), f().ncols(),
                               f().type(),0.0));
      }
      assert(dfdx().nrows()==val.dfdx().nrows());
      assert(dfdx().ncols()==val.dfdx().ncols());
      f()(i,j)=val.f();
      for (std::size_t n = 0; n < dfdx().size(); ++n)
          dfdx()[n](i,j)=val.dfdx()[n];
  }
  void add(std::size_t i, std::size_t j, const Derivative<double>& val){
      assert((&val.x())!=nullptr);
      if (x_==nullptr)
      {
          x_=&val.x();
          dfdx_=M_Matrix<M_Matrix<double>>(
              val.x().nrows(), val.x().ncols(), val.x().type(),
              M_Matrix<double>(f().nrows(), f().ncols(),
                               f().type(),0.0));
      }
     assert(dfdx().nrows()==val.dfdx().nrows());
     assert(dfdx().ncols()==val.dfdx().ncols());
      f()(i,j)+=val.f();
      for (std::size_t n = 0; n < dfdx().size(); ++n)
          dfdx()[n](i,j)+=val.dfdx()[n];
  }

  bool check() {
    bool out = true;
    for (std ::size_t i = 0; i < dfdx_.size(); ++i)
      if (f_.type() != dfdx()[i].type())
        out = false;
    return out;
  }
  Derivative(std::size_t nrows, std::size_t ncols, Matrix_TYPE atype=Matrix_TYPE::FULL)
      :f_{nrows,ncols,atype},x_{nullptr},dfdx_{}{}

  Derivative(std::size_t nrows, std::size_t ncols, Matrix_TYPE atype, double val)
      :f_{nrows,ncols,atype,val},x_{nullptr},dfdx_{}{}


  Derivative(const M_Matrix<double> &fx, const M_Matrix<double> &myx,
             const M_Matrix<M_Matrix<double>> &der)
      : f_{fx}, x_{&myx}, dfdx_{der} {
    assert(check());
  }

  Derivative(const M_Matrix<double> &fx, const M_Matrix<double> &myx)
      : f_{fx}, x_{&myx}, dfdx_{M_Matrix<M_Matrix<double>>(
                              myx.nrows(), myx.ncols(), myx.type(),
                              M_Matrix<double>(fx.nrows(), fx.ncols(),
                                               fx.type()))} {
    assert(check());
  }

  template <class M_f, class M_d>
  Derivative(M_f &&fx, const M_Matrix<double> &myx, M_d &&der)
      : f_{std::forward<M_f>(fx)}, x_{&myx}, dfdx_{std::forward<M_d>(der)} {
    assert(check());
  }

  Derivative(std::size_t nrows, std::size_t ncols, M_Matrix<double>::TYPE t,
             const M_Matrix<double> &x)
      : f_{nrows, ncols, t}, x_{&x}, dfdx_{x.nrows(), x.ncols(), x.type(),
                                           M_Matrix<double>(nrows, ncols, t)} {
    assert(check());
  }
  Derivative() = default;

  Derivative(const M_Matrix<Derivative<double>> &x)
      : f_{x.nrows(), x.ncols(), x.type()}, x_{&x[0].x()},
        dfdx_{x[0].x().nrows(), x[0].x().ncols(), x[0].x().type(),
              M_Matrix<double>(x.nrows(), x.ncols(), x.type())} {
    for (std::size_t i = 0; i < x.size(); ++i) {
      f_[i] = x[i].f();
      for (std::size_t j = 0; j < x[0].x().size(); ++j)
        dfdx_[j][i] = x[i].dfdx()[j];
      assert(check());
    }
  }
  M_Matrix<Derivative<double>> to_Matrix() const {
    M_Matrix<Derivative<double>> out(f().nrows(), f().ncols(), f().type());
    for (std::size_t i = 0; i < f().size(); ++i) {
      M_Matrix<double> dfdxi(x());
      for (std::size_t j = 0; j < x().size(); ++j)
        dfdxi[j] = dfdx()[j][i];
      out[i] = Derivative<double>(f()[i], x(), std::move(dfdxi));
    }
    return out;
  }

  template <class F> Derivative apply(const F &f) {
    auto M = to_Matrix();
    return Derivative(M.apply(f));
  }
  Derivative &operator+=(const Derivative &other) {
    assert(x() == other.x());
    f_ += other.f();
    dfdx_ += other.dfdx();
    return *this;
  }


  Derivative<double> getvalue() && {
    {
      assert(f().size() == 1);
      M_Matrix<double> df(x().nrows(), x().ncols(), x().type());
      for (std::size_t i = 0; i < x().size(); ++i)
        df[i] = dfdx()[i][0];
      return Derivative<double>(f().getvalue(), x(), std::move(df));
    }
  }
};


inline auto& center(const Derivative<M_Matrix<double>>& x){return x.f();}
inline auto& center(Derivative<M_Matrix<double>>& x){return x.f();}


template <bool output> struct are_Equal<output, Derivative<double>> {
private:
  double absolute_;
  double relative_;

public:
  are_Equal(double abs, double rel) : absolute_{abs}, relative_{rel} {}

  bool test(const Derivative<double> &x, const Derivative<double> &y,
            std::ostream &os) {
    bool res = true;
    if (!are_Equal_v(x.f(), y.f(), absolute_, relative_, os)) {
      os << " failed in f()";
      res = false;
    }
    if (!are_Equal_v(x.x(), y.x(), absolute_, relative_, os)) {
      os << " failed in x()";
      res = false;
    }
    if (!are_Equal_v(x.dfdx(), y.dfdx(), absolute_, relative_, os)) {
      os << " failed in x()";
      res = false;
    }
    return res;
  }
};

template <bool output, typename T>
struct are_Equal<output, Derivative<M_Matrix<T>>> {
public:
  bool test(const Derivative<M_Matrix<T>> &one,
            const Derivative<M_Matrix<T>> &two, std::ostream &os) {
    return test_prod(one, two, os);
  }

  template <class ostream>
  bool test_sum(const Derivative<M_Matrix<T>> &one,
                const Derivative<M_Matrix<T>> &two,
                ostream &os = std::cerr) const {
    std::stringstream ss_f;
    bool equal_f = are_Equal<output, M_Matrix<T>>(absolute_, relative_)
                       .test_sum(one.f(), two.f(), ss_f);
    if (!equal_f)
      os << " error in function!!!\n " << ss_f.str()
         << " \n error in function  END!!!\n";
    bool equal_dfdx = true;
    for (std::size_t i = 0; i < one.x(); ++i)

    {
      std::stringstream ss_df;
      bool equal_dfdx_i = are_Equal<output, M_Matrix<T>>(absolute_, relative_)
                              .test_sum(one.dfdx()[i], two.dfdx()[i], ss_df);
      if (!equal_dfdx_i)
        os << " error in the" << i << "th derivative function!!!\n "
           << ss_df.str() << " error in the" << i
           << "th derivative function END!!!\n";
      if (!equal_dfdx_i)
        equal_dfdx = false;
    }

    return equal_f && equal_dfdx;
  }
  template <class ostream>
  bool test_prod(const Derivative<M_Matrix<T>> &one,
                 const Derivative<M_Matrix<T>> &two,
                 ostream &os = std::cerr) const {
    std::stringstream ss_f;
    bool equal_f = are_Equal<output, M_Matrix<T>>(absolute_, relative_)
                       .test_prod(one.f(), two.f(), ss_f);
    if (!equal_f)
      os << " error in function!!!\n " << ss_f.str()
         << " \n error in function  "
            "END----------------------!!!\n---------------------\n";
    bool equal_dfdx = true;
    for (std::size_t i = 0; i < one.x().size(); ++i)

    {
      std::stringstream ss_df;
      bool equal_dfdx_i = are_Equal<output, M_Matrix<T>>(absolute_, relative_)
                              .test_prod(one.dfdx()[i], two.dfdx()[i], ss_df);
      if (!equal_dfdx_i)
        os << " \n error in the" << i << "th derivative function!!!\n "
           << ss_df.str() << " error in the" << i
           << "\n th derivative function "
              "END-------------------------------!!!\n";
      if (!equal_dfdx_i)
        equal_dfdx = false;
    }

    return equal_f && equal_dfdx;
  }

  template <class ostream, class... Ts, std::size_t... Is>
  bool test_prod_impl(const std::tuple<Derivative<M_Matrix<Ts>>...> &one,
                      const std::tuple<Derivative<M_Matrix<Ts>>...> &two,
                      ostream &os, std::index_sequence<Is...>) const {
    return (test_prod(std::get<Is>(one), std::get<Is>(two), os) && ...);
  }

  template <class ostream, class... Ts>
  bool test_prod(const std::tuple<Derivative<M_Matrix<Ts>>...> &one,
                 const std::tuple<Derivative<M_Matrix<Ts>>...> &two,
                 ostream &os = std::cerr) const {
    return test_prod_impl(one, two, os, std::index_sequence_for<Ts...>());
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

inline double Taylor_first(const Derivative<double> &dx, std::size_t i,
                           double eps) {
  return dx.f() + dx.dfdx()[i] * eps;
}

inline bool isnan(const Derivative<double> &dx) { return std::isnan(dx.f()); }

template <class F>
Derivative<double> Incremental_ratio(double eps, const F &fun,
                                     const Derivative<double> &y) {
  double f = fun(y.f());
  double fpos = fun(y.f() + eps);
  double fneg = fun(y.f() - eps);
  double dfdy = (fpos - fneg) / (2 * eps);
  M_Matrix<double> dfdx = dfdy * y.dfdx();
  return Derivative<double>(f, y.x(), dfdx);
}

inline M_Matrix<double> Taylor_first(const Derivative<M_Matrix<double>> &dx,
                                     std::size_t i, double eps) {
  return dx.f() + dx.dfdx()[i] * eps;
}

inline bool isnan(const Derivative<M_Matrix<double>> &dx) {
  return isnan(dx.f());
}

inline double Directional_Derivative(const double &dx, const M_Matrix<double> &,
                                     std::size_t, double) {
  return dx;
}

inline Derivative<double> Directional_Derivative(const Derivative<double> &dx,
                                                 const M_Matrix<double> &new_x,
                                                 std::size_t i, double eps) {
  return Derivative<double>(Taylor_first(dx, i, eps), new_x,
                            M_Matrix<double>(1, 1, dx.dfdx()[i]));
}

inline Derivative<M_Matrix<double>>
Directional_Derivative(const Derivative<M_Matrix<double>> &dx,
                       const M_Matrix<double> &newx, std::size_t i,
                       double eps) {
  return Derivative<M_Matrix<double>>(
      Taylor_first(dx, i, eps), newx,
      M_Matrix<M_Matrix<double>>(1, 1, dx.dfdx()[i]));
}

template <class F>
Derivative<M_Matrix<double>>
Incremental_ratio(double eps, const F &fun,
                  const Derivative<M_Matrix<double>> &y) {
  M_Matrix<double> f = fun(y.f());
  M_Matrix<M_Matrix<double>> dfdx(y.x().nrows(), y.x().ncols(), y.x().type());
  for (std::size_t i = 0; i < y.x().size(); ++i) {
    double e = eps;

    if (std::abs(y.x()[i]) > 0)
      e *= std::abs(y.x()[i]);
    auto ypos = Taylor_first(y, i, e);
    auto yneg = Taylor_first(y, i, -e);
    auto fpos = fun(ypos);
    auto fneg = fun(yneg);
    dfdx[i] = (fpos - fneg) / (2 * e);
  }
  return Derivative<M_Matrix<double>>(f, y.x(), dfdx);
}

template <class T> M_Matrix<T> flatten(const M_Matrix<M_Matrix<T>> &dfdx) {
  auto ncols = dfdx.size();
  auto nrows = dfdx[0].size();
  M_Matrix<T> out(nrows, ncols, Matrix_TYPE::FULL);
  for (std::size_t i = 0; i < nrows; ++i)
    for (std::size_t j = 0; j < ncols; ++j)
      out(i, j) = dfdx[j][i];
  return out;
}

template <class T, class K, class L>
M_Matrix<M_Matrix<T>> unflatten(const M_Matrix<T> &dfdx, const M_Matrix<K> &f,
                                const M_Matrix<L> &x) {
  assert(dfdx.ncols() == x.size());
  assert(dfdx.nrows() == f.size());
  auto ncols = dfdx.ncols();
  auto nrows = dfdx.nrows();
  M_Matrix<M_Matrix<T>> out(x.nrows(), x.ncols(), x.type(),
                            M_Matrix<T>(f.nrows(), f.ncols(), f.type()));
  for (std::size_t i = 0; i < nrows; ++i)
    for (std::size_t j = 0; j < ncols; ++j)
      out[j][i] = dfdx(i, j);
  return out;
}

template <class T>
M_Matrix<M_Matrix<T>> compose(const M_Matrix<M_Matrix<T>> &dfdy,
                              const M_Matrix<M_Matrix<T>> &dydx) {
  return unflatten(flatten(dfdy) * flatten(dydx), dfdy[0], dydx);
}

template <class F>
Derivative<M_Matrix<double>>
Incremental_ratio_compose(double eps, const F &fun,
                          const Derivative<M_Matrix<double>> &y) {
  M_Matrix<double> f = fun(y.f());
  M_Matrix<M_Matrix<double>> dfdy(y.f().nrows(), y.f().ncols(), y.f().type());
  for (std::size_t i = 0; i < y.f().size(); ++i) {
    auto ypos = y.f();
    ypos[i] += eps;
    auto yneg = y.f();
    yneg[i] -= eps;

    auto fpos = fun(ypos);
    auto fneg = fun(yneg);
    dfdy[i] = (fpos - fneg) / (2 * eps);
  }
  return Derivative<M_Matrix<double>>(f, y.x(), compose(dfdy, y.dfdx()));
}

template <class F>
Derivative<M_Matrix<double>>
Incremental_ratio(double eps, const F &fun,
                  const Derivative<M_Matrix<double>> &one,
                  const Derivative<M_Matrix<double>> &other) {
  assert(one.x() == other.x());
  M_Matrix<double> f = fun(one.f(), other.f());
  M_Matrix<M_Matrix<double>> dfdx(one.x().nrows(), one.x().ncols(),
                                  one.x().type());
  for (std::size_t i = 0; i < one.x().size(); ++i) {
    auto one_pos = Taylor_first(one, i, eps);
    auto other_pos = Taylor_first(other, i, eps);
    auto one_neg = Taylor_first(one, i, -eps);
    auto other_neg = Taylor_first(other, i, -eps);

    auto fpos = fun(one_pos, other_pos);
    auto fneg = fun(one_neg, other_neg);
    dfdx[i] = (fpos - fneg) / (2 * eps);
  }
  return Derivative<M_Matrix<double>>(f, one.x(), dfdx);
}

template <class F, typename T, typename... Ds>
auto Incremental_ratio(double eps, const F &fun, const Derivative<T> &y0,
                       const Ds &... y) {
  assert(((y0.x() == y.x()) && ...));
  auto f = fun(y0.f(), Primitive(y)...);
  typedef decltype(f) R;
  auto &one = y0;
  M_Matrix<R> dfdx(one.x().nrows(), one.x().ncols(), one.x().type());
  for (std::size_t i = 0; i < one.x().size(); ++i) {
    double e = eps;
    if (one.x()[i] != 0)
      e *= std::abs(one.x()[i]);
    auto fpos = fun(Taylor_first(y0, i, e), Taylor_first(y, i, e)...);
    auto fneg = fun(Taylor_first(y0, i, -e), Taylor_first(y, i, -e)...);
    dfdx[i] = (fpos - fneg) / (2.0 * e);
  }
  return Derivative<R>(f, one.x(), dfdx);
}

template <class F>
std::tuple<Derivative<M_Matrix<double>>, Derivative<M_Matrix<double>>,
           Derivative<M_Matrix<double>>>
Incremental_ratio_tuple_3_compose(double eps, const F &fun,
                                  const Derivative<M_Matrix<double>> &y) {
  auto [f1, f2, f3] = std::invoke(fun, y.f());
  M_Matrix<M_Matrix<double>> df1dy(y.f().nrows(), y.f().ncols(), y.f().type());
  M_Matrix<M_Matrix<double>> df2dy(y.f().nrows(), y.f().ncols(), y.f().type());
  M_Matrix<M_Matrix<double>> df3dy(y.f().nrows(), y.f().ncols(), y.f().type());
  for (std::size_t i = 0; i < y.f().size(); ++i) {
    auto ypos = y.f();
    ypos[i] += eps;
    auto yneg = y.f();
    yneg[i] -= eps;

    auto [fpos1, fpos2, fpos3] = fun(ypos);
    auto [fneg1, fneg2, fneg3] = fun(yneg);
    df1dy[i] = (fpos1 - fneg1) / (2.0 * eps);
    df2dy[i] = (fpos2 - fneg2) / (2.0 * eps);
    df3dy[i] = (fpos3 - fneg3) / (2.0 * eps);
  }
  return std::tuple(
      Derivative<M_Matrix<double>>(f1, y.x(), compose(df1dy, y.dfdx())),
      Derivative<M_Matrix<double>>(f2, y.x(), compose(df2dy, y.dfdx())),
      Derivative<M_Matrix<double>>(f3, y.x(), compose(df3dy, y.dfdx())));
}

template <class F>
std::tuple<Derivative<M_Matrix<double>>, Derivative<M_Matrix<double>>,
           Derivative<M_Matrix<double>>>
Incremental_ratio_tuple_3(double eps, const F &fun,
                          const Derivative<M_Matrix<double>> &y) {
  auto [f1, f2, f3] = fun(y.f());
  M_Matrix<M_Matrix<double>> df1dx(y.x().nrows(), y.x().ncols(), y.x().type());
  M_Matrix<M_Matrix<double>> df2dx(y.x().nrows(), y.x().ncols(), y.x().type());
  M_Matrix<M_Matrix<double>> df3dx(y.x().nrows(), y.x().ncols(), y.x().type());
  for (std::size_t i = 0; i < y.x().size(); ++i) {
    double e = eps;
    if (y.x()[i] != 0)
      e *= std::abs(y.x()[i]);
    auto ypos = Taylor_first(y, i, e);
    auto yneg = Taylor_first(y, i, -e);

    auto [fpos1, fpos2, fpos3] = fun(ypos);
    auto [fneg1, fneg2, fneg3] = fun(yneg);
    df1dx[i] = (fpos1 - fneg1) / (2.0 * e);
    df2dx[i] = (fpos2 - fneg2) / (2.0 * e);
    df3dx[i] = (fpos3 - fneg3) / (2.0 * e);
  }
  return std::tuple(Derivative<M_Matrix<double>>(f1, y.x(), df1dx),
                    Derivative<M_Matrix<double>>(f2, y.x(), df2dx),
                    Derivative<M_Matrix<double>>(f3, y.x(), df3dx));
}

template <class C>
auto Primitive(const Derivative<C> &dy) -> std::enable_if_t<
    is_field_Object<C>::value && is_field_Object<Derivative<C>>::value, C> {
  auto fields = dy.get_constructor_fields();
  return std::apply(
      [&dy](auto &... m) {
        return C(Primitive(std::invoke(m.access_method, dy))...);
      },
      fields);
}

template <class C>
C Taylor_first(const Derivative<C> &dy, std::size_t i, double eps) {
  auto fields = dy.get_constructor_fields();
  return std::apply(
      [&i, &eps, &dy](auto &... m) {
        return C(Taylor_first(std::invoke(m.access_method, dy), i, eps)...);
      },
      fields);
}

template <class C>
Derivative<C> Directional_Derivative(const Derivative<C> &dy,
                                     const M_Matrix<double> &new_x,
                                     std::size_t i, double eps) {
  auto fields = dy.get_constructor_fields();
  return std::apply(
      [&i, &dy, &eps, &new_x](auto &... m) {
        return Derivative<C>(Directional_Derivative(
            std::invoke(m.access_method, dy), new_x, i, eps)...);
      },
      fields);
}

template <class C>
auto Incremental_ratio_calc(const std::vector<double> &eps, const C &f,
                            const M_Matrix<double> &x,
                            const std::vector<C> &fpos,
                            const std::vector<C> &fneg)
    -> std::enable_if_t<is_field_Object<C>::value &&
                            is_field_Object<Derivative<C>>::value,
                        Derivative<C>> {
  auto fields = f.get_constructor_fields();
  auto dfields = Derivative<C>::get_constructor_fields();
  return myApply(
      [&f, &x, &fpos, &fneg, &eps](auto &&... mdm) {
        return Derivative<C>(Incremental_ratio_method_calc(
            eps, mdm.first, mdm.second, f, x, fpos, fneg)...);
      },
      std::move(fields), std::move(dfields));
}

inline Derivative<double> Incremental_ratio_calc(
    const std::vector<double> &eps, double f, const M_Matrix<double> &x,
    const std::vector<double> &fpos, const std::vector<double> &fneg) {
  Derivative<double> out(f, x);
  for (std::size_t i = 0; i < x.size(); ++i)
    out.dfdx()[i] = (fpos[i] - fneg[i]) * (1.0 / (2. * eps[i]));
  return out;
}

inline Derivative<M_Matrix<double>>
Incremental_ratio_calc(const std::vector<double> &eps,
                       const M_Matrix<double> &f, const M_Matrix<double> &x,
                       const std::vector<M_Matrix<double>> &fpos,
                       const std::vector<M_Matrix<double>> &fneg) {
  Derivative<M_Matrix<double>> out(f, x);
  for (std::size_t i = 0; i < x.size(); ++i)
    out.dfdx()[i] = (fpos[i] - fneg[i]) * (1.0 / (2. * eps[i]));
  return out;
}

template <class C, class method, class dmethod>
auto Incremental_ratio_method_calc(
    const std::vector<double> &eps, const grammar::field<C, method> m,
    const grammar::field<Derivative<C>, dmethod> /*dm*/, const C &f,
    const M_Matrix<double> &x, const std::vector<C> &fpos,
    const std::vector<C> &fneg) {
  typedef std::invoke_result_t<method, C> R;
  typedef std::invoke_result_t<dmethod, Derivative<C>> dR;

  if constexpr (std::is_same_v<R, dR>) {
    return std::invoke(m.access_method, f);
  } else {
    auto rf = std::invoke(m.access_method, f);
    typedef decltype(rf) RR;
    std::vector<RR> rfpos(x.size());
    std::vector<RR> rfneg(x.size());

    for (std::size_t i = 0; i < x.size(); ++i) {
      rfpos[i] = std::invoke(m.access_method, fpos[i]);
      rfneg[i] = std::invoke(m.access_method, fneg[i]);
    }
    return Incremental_ratio_calc(eps, rf, x, rfpos, rfneg);
  }
}

template <class F, class Object>
auto Incremental_ratio_object(double eps, const F &fun,
                              const Derivative<Object> &y)
    -> Derivative<std::invoke_result_t<F, Object>> {

  typedef std::invoke_result_t<F, Object> R;
  auto f = std::invoke(fun, y.f());
  std::vector<R> fpos(y.x().size());
  std::vector<R> fneg(y.x().size());
  for (std::size_t i = 0; i < y.x().size(); ++i) {
    double e = eps;
    if (y.x()[i] != 0)
      e *= std::abs(y.x()[i]);
    fpos[i] = std::invoke(fun, Taylor_first(y, i, e));
    fneg[i] = std::invoke(fun, Taylor_first(y, i, -e));
  }
  return Incremental_ratio_calc(eps, f, y.x(), fpos, fneg);
}

inline bool mean_value_test(double eps, double tol, double f, double feps,
                            double df, double dfeps, std::ostream &os) {
  if (std::isnan(f + feps + df + dfeps) ||
      (std::abs((feps - f) - (df + dfeps) * eps / 2.0) >
       (200.0 * std::abs((dfeps - df) * eps) +
        std::sqrt(std::numeric_limits<double>::epsilon())) *
           tol)) {
    double deltaf = (feps - f) / eps;
    double dfm = (df + dfeps) / 2;
    double dferr = std::abs(dfeps - df);
    os << "\n eps= " << eps << " f= " << f << " feps= " << feps << " df= " << df
       << " dfeps= " << dfeps << "\n";
    os << "df_analytic_mean =\t " << dfm << "\t df_numeric =\t " << deltaf
       << "\tdifference * eps\t=" << std ::abs(dfm - deltaf) * eps;
    os << "\t df_analytic_error * eps=\t" << dferr * eps;
    return false;
  } else
    return true;
}
inline bool mean_value_test_calc(const std::vector<double> &eps, double tol,
                                 const Derivative<double> &y,
                                 const std::vector<Derivative<double>> &yeps,
                                 std::ostream &os) {
  bool out = true;
  for (std::size_t j = 0; j < yeps.size(); ++j) {
    std::stringstream ss;
    if (!mean_value_test(eps[j], tol, y.f(), yeps[j].f(), y.dfdx()[j],
                         yeps[j].dfdx()[0], ss)) {
      os << "\n Fails for the " << j << "th derivative: " << ss.str();
      out = false;
    }
  }
  return out;
}

inline bool mean_value_test_calc(
    std::vector<double> eps, double tol, const Derivative<M_Matrix<double>> &y,
    const std::vector<Derivative<M_Matrix<double>>> &yeps, std::ostream &os) {
  bool out = true;
  for (std::size_t j = 0; j < yeps.size(); ++j) {
    bool jder = true;
    std::stringstream ssder;
    for (std::size_t i = 0; i < y.f().size(); ++i) {
      std::stringstream ss;
      if (!mean_value_test(eps[j], tol, y.f()[i], yeps[j].f()[i],
                           y.dfdx()[j][i], yeps[j].dfdx()[0][i], ss)) {
        ssder << "\nat " << y.f().pos_to_ij(i) << ss.str();
        jder = false;
      }
    }
    if (!jder) {
      os << "\n Fails for the " << j << "th derivative: " << ssder.str();
      os << " \n y.f()\n"
         << y.f() << "\n yeps.f()\n"
         << yeps[j].f() << "\n y.dfdx[j]\n"
         << y.dfdx()[j] << "\n yeps[j].dfdx()\n"
         << yeps[j].dfdx()[0];
      os << "\n----------ends the " << j << "th derivative--\n";
      out = false;
    }
  }
  return out;
}

template <class Object, class Method, class DerivativeObject>
bool mean_value_test_method_calc(const std::vector<double> &eps, double tol [[maybe_unused]],
                                 const grammar::field<Object, Method> &m,
                                 const DerivativeObject &y,
                                 const std::vector<DerivativeObject> &yeps,
                                 std::ostream &os) {
  std::stringstream ss;
  typedef std::decay_t<std::invoke_result_t<Method, Object>> FieldObject;
  if constexpr (!is_Derivative_v<FieldObject>) {
    // typedef typename FieldObject::NO_DER gege;
    return true;

  } else {
    // typedef FieldObject DerFieldObject;
    // typedef typename DerFieldObject::IS_DER gege;

    auto my = std::invoke(m.access_method, y);
    std::vector<FieldObject> myeps(yeps.size());
    for (std::size_t i = 0; i < myeps.size(); ++i)
      myeps[i] = (std::invoke(m.access_method, yeps[i]));
    if (!mean_value_test_calc(eps, tol, my, myeps, ss)) {
      os << "\nMean value test fails for Field " << m.idField << "\n"
         << ss.str();
      return false;
    } else
      return true;
  }
}
template <class DerivativeObject>
auto mean_value_test_calc(const std::vector<double> &eps, double tol,
                          const DerivativeObject &y,
                          const std::vector<DerivativeObject> &yeps,
                          std::ostream &os)
    -> std::enable_if_t<is_field_Object<DerivativeObject>::value, bool> {
  std::stringstream ss;
  auto dfields = y.get_constructor_fields();
  bool out = std::apply(
      [&eps, &y, &yeps, &ss, &tol](auto &... m) {
        return (mean_value_test_method_calc(eps, tol, m, y, yeps, ss) * ...);
      },
      dfields);
  if (!out) {
    os << "\nMean value test fails for class "
       << my_trait<DerivativeObject>::className.str() << "\n"
       << ss.str();
    return false;
  } else
    return true;
}
template <class Diff, class DerX, class DerY, typename... Ts>
bool Derivative_correctness_mean_value_test(double eps, double tol,
                                            const Diff &Dfunction,
                                            const DerX &x0, const DerY &y,
                                            std::ostream &os, Ts... context) {

  std::vector<DerY> yeps(y.x().size());
  std::vector<double> e(y.x().size());

  for (std::size_t i = 0; i < y.x().size(); ++i) {
    e[i] = eps;
    if (std::abs(y.x()[i]) > 0)
      e[i] *= std ::abs(y.x()[i]);

    M_Matrix<double> new_x(1, 1, y.x()[i] + e[i]);

    auto xpos = Directional_Derivative(x0, new_x, i, e[i]);
    yeps[i] = std::invoke(Dfunction, xpos);
  }
  std::stringstream ss;
  bool out = mean_value_test_calc(e, tol, y, yeps, os);
  if (!out) {
    os << "\ncontext :";
    (os << ... << context);
  }
  return out;
}

template <class Diff, class... Ders, class Y, typename... Ts>
bool Derivative_correctness_mean_value_test(double eps, double tol,
                                            const Diff &Dfunction,
                                            const std::tuple<Ders...> &x0,
                                            const Derivative<Y> &y,
                                            std::ostream &os, Ts... context) {

  std::vector<Derivative<Y>> yeps(y.x().size());
  std::vector<double> e(y.x().size());
  for (std::size_t i = 0; i < y.x().size(); ++i) {
    e[i] = eps;
    if (std::abs(y.x()[i]) > 0)
      e[i] *= std ::abs(y.x()[i]);

    M_Matrix<double> new_x(1, 1, y.x()[i] + e[i]);
    yeps[i] = std::apply(
        [&Dfunction, &i, &new_x, &e](auto const &... t) {
          return std::invoke(Dfunction,
                             Directional_Derivative(t, new_x, i, e[i])...);
        },
        x0);
    if constexpr (std::is_same_v<Y, double> ||
                  std::is_same_v<Y, M_Matrix<double>>) {
      while (!isnan(y) && isnan(yeps[i])) {
        e[i] *= 0.125;
        M_Matrix<double> new_x(1, 1, y.x()[i] + e[i]);
        yeps[i] = std::apply(
            [&Dfunction, &i, &new_x, &e](auto const &... t) {
              return std::invoke(Dfunction,
                                 Directional_Derivative(t, new_x, i, e[i])...);
            },
            x0);
      }
    }
  }
  std::stringstream ss;
  bool out = mean_value_test_calc(e, tol, y, yeps, os);
  if (!out) {
    os << "\n context: ";
    (os << ... << context);
    os << "\n";
  }
  return out;
}


inline Derivative<double> operator+(const Derivative<double> &x,  double c) {
    return Derivative<double>(x.f() + c, x.x(), x.dfdx());
}

inline Derivative<double> operator+(double c, const Derivative<double> &x) {
    return Derivative<double>(c+x.f(), x.x(), x.dfdx());
}

inline Derivative<double> operator-(const Derivative<double> &x,  double c) {
    return Derivative<double>(x.f() - c, x.x(), x.dfdx());
}

inline Derivative<double> operator-(double c, const Derivative<double> &x) {
    return Derivative<double>(c-x.f(), x.x(), x.dfdx());
}




template <class T>
Derivative<T> operator+(const Derivative<T> &x, Constant<T> &&c) {
  return Derivative<T>(x.f() + c.value, x.x(), x.dfdx());
}

template <class T>
Derivative<T> operator+(Constant<T> &&c, const Derivative<T> &x) {
  return Derivative<T>(c.value + x.f(), x.x(), x.dfdx());
}
template <class T>
Derivative<T> operator-(Constant<T> &&c, const Derivative<T> &x) {
  return Derivative<T>(c.value - x.f(), x.x(), -x.dfdx());
}

template <class T>
Derivative<T> operator*(const Derivative<T> &x, const Constant<T> &c) {
  return Derivative<T>(x.f() * c.value, x.x(),
                       x.dfdx().apply([&c](auto &m) { return m * c.value; }));
}

template <class T>
Derivative<T> operator*(const Constant<T> &c, const Derivative<T> &x) {
  return Derivative<T>(c.value * x.f(), x.x(),
                       x.dfdx().apply([&c](auto &m) { return c.value * m; }));
}

inline auto exp(D, double x) { return std::exp(x); }

inline auto exp(const Derivative<double> &x) {
  return Derivative<double>(std::exp(x.f()), x.x(), exp(D(), x.f()) * x.dfdx());
}

inline auto log(const Derivative<double> &x) {
  return Derivative<double>(std::log(x.f()), x.x(), x.dfdx() / x.f());
}

inline auto sqrt(const Derivative<double> &x) {
  return Derivative<double>(std::sqrt(x.f()), x.x(),
                            (0.5 / std::sqrt(x.f())) * x.dfdx());
}

inline auto exp(const Derivative<M_Matrix<double>> &x) {
  return Derivative<M_Matrix<double>>(exp(x.f()), x.x(),
                                      elemMult_a(x.dfdx(), exp(x.f())));
}






inline auto expm1(const Derivative<double> &x) {
    return Derivative<double>(std::expm1(x.f()), x.x(),
                              exp(x.f())*x.dfdx());
}



inline auto abs(D, double x) {
  if (x >= 0)
    return 1.0;
  else
    return -1.0;
}

inline auto abs(const Derivative<double> &x) {
  return Derivative<double>(std::abs(x.f()), x.x(), abs(D(), x.f()) * x.dfdx());
}

inline auto pow(const Derivative<double> &x, const Derivative<double> &y) {
  assert(x.x() == y.x());
  auto f = std::pow(x.f(), y.f());

  auto dx = x.dfdx() * (y.f() / x.f() * f);
  auto dy = y.dfdx() * (std::log(x.f()) * f);
  return Derivative<double>(std::move(f), x.x(), dx + dy);
}

inline auto operator/(const Derivative<double> &x, const Derivative<double> &y) {
  assert(x.x() == y.x());
  auto f = x.f() / y.f();

  auto dx = x.dfdx() / y.f();
  auto dy = -y.dfdx() * x.f() / sqr(y.f());
  return Derivative<double>(std::move(f), x.x(), dx + dy);
}

inline auto operator/(const Derivative<double> &x, double y) {
    return Derivative<double>(x.f() / y, x.x(), x.dfdx() / y);
}

inline auto operator/(double x, const Derivative<double> &y) {
  auto f = x / y.f();

  auto dy = -y.dfdx() * x / sqr(y.f());
  return Derivative<double>(std::move(f), y.x(), std::move(dy));
}

template <typename T>
auto TranspMult(const Derivative<M_Matrix<T>> &x,
                const Derivative<M_Matrix<T>> &y) {
  assert(x.x() == y.x());
  auto f = TranspMult(x.f(), y.f());
  auto dfdx = x.dfdx().apply([&y](auto &m) { return TranspMult(m, y.f()); });
  auto dfdy = y.dfdx().apply([&x](auto &m) { return TranspMult(x.f(), m); });
  return Derivative<M_Matrix<T>>(std::move(f), x.x(), dfdx + dfdy);
}

template <typename T>
auto multTransp(const Derivative<M_Matrix<T>> &x,
                const Derivative<M_Matrix<T>> &y) {
  assert(x.x() == y.x());
  auto f = multTransp(x.f(), y.f());
  auto dfdx = x.dfdx().apply([&y](auto &m) { return multTransp(m, y.f()); });
  auto dfdy = y.dfdx().apply([&x](auto &m) { return multTransp(x.f(), m); });
  return Derivative<M_Matrix<T>>(std::move(f), x.x(), dfdx + dfdy);
}

template <typename T>
auto quadraticForm_BT_A_B(const Derivative<M_Matrix<T>> &A,
                          const Derivative<M_Matrix<T>> &B) {
  assert(A.x() == B.x());
  auto f = quadraticForm_BT_A_B(A.f(), B.f());
  auto dfdA =
      A.dfdx().apply([&B](auto &m) { return quadraticForm_BT_A_B(m, B.f()); });
  if (A.f().isSymmetric()) {
    auto dfdB = B.dfdx().apply([&A, &B](auto &m) {
      return TransposeSum(TranspMult(m, A.f() * B.f()));
    });
    return Derivative<M_Matrix<T>>(std::move(f), A.x(), dfdA + dfdB);
  } else {
    auto dfdB = B.dfdx().apply([&A, &B](auto &m) {
      return TranspMult(m, A.f() * B.f()) + TranspMult(B.f(), A.f() * m);
    });
    return Derivative<M_Matrix<T>>(std::move(f), A.x(), dfdA + dfdB);
  }
}

template <typename T>
auto elemMult(const Derivative<M_Matrix<T>> &x,
              const Derivative<M_Matrix<T>> &y) {
  assert(x.x() == y.x());
  return Derivative<M_Matrix<T>>(elemMult(x.f(), y.f()), x.x(),
                                 elemMult_a(x.dfdx(), y.f()) +
                                     elemMult_a(x.f(), y.dfdx()));
}

template <typename T>
auto elemDiv(const Derivative<M_Matrix<T>> &x,
             const Derivative<M_Matrix<T>> &y) {
  assert(x.x() == y.x());
  auto f = elemDiv(x.f(), y.f());
  auto dfdx = x.dfdx().apply([&y](auto &xi) { return elemDiv(xi, y.f()); });
  auto dfdys = elemDiv(f, y.f());
  auto dfdy =
      y.dfdx().apply([&dfdys](auto &yi) { return elemMult(yi, dfdys); });
  return Derivative<M_Matrix<T>>(std::move(f), x.x(), dfdx - dfdy);
}

template <typename T>
auto elemDivSafe(const Derivative<M_Matrix<T>> &x,
                 const Derivative<M_Matrix<T>> &y,
                 double eps = std::numeric_limits<double>::epsilon()) {
  assert(x.x() == y.x());
  auto f = elemDivSafe(x.f(), y.f(), eps);
  auto dfdx = elemDivSafe_a(x.dfdx(), y.f(), eps);
  auto dfdy = elemMult_a(elemDivSafe_a(y.dfdx(), y.f(), eps), f);

  return Derivative<M_Matrix<T>>(std::move(f), x.x(), dfdx - dfdy);
}

/*
template< class T>
auto compose (const Derivative<T>& dfdx, const Derivative<M_Matrix<double>>&
dxdy)
{
    assert(dfdx.x()==dxdy.f());
    return Derivative<T>(dfdx.f(), dxdy.x(), dfdx.dfdx()*dxdy.dfdx());
}
*/

template <typename T, typename S>
auto operator*(const Derivative<T> &x, S t)
    -> Derivative<std::enable_if_t<!is_Derivative_v<S>, T>> {
  return Derivative<T>(x.f() * t, x.x(), x.dfdx() * t);
}

template <typename T, typename S>
auto operator*(S t, const Derivative<T> &x)
    -> Derivative<std::enable_if_t<!is_Derivative_v<S>, T>> {
  return Derivative<T>(t * x.f(), x.x(), t * x.dfdx());
}

template <typename T>
Derivative<T> operator+(const Derivative<T> &x, const Derivative<T> &y)

{
  assert(x.x() == y.x());
  return Derivative<T>(x.f() + y.f(), x.x(), x.dfdx() + y.dfdx());
}

template <class T>
Derivative<M_Matrix<T>> TransposeSum(const Derivative<M_Matrix<T>> &x) {
  assert(x.f().ncols() == x.f().nrows());
  auto f = TransposeSum(x.f());
  auto df = x.dfdx().apply([](auto &x) { return TransposeSum(x); });
  return Derivative<M_Matrix<T>>(std::move(f), x.x(), std::move(df));
}

template <class T>
Derivative<M_Matrix<T>> quadraticForm_XTX(const Derivative<M_Matrix<T>> &x) {
  auto f = quadraticForm_XTX(x.f());
  auto df = x.dfdx().apply(
      [&x](auto &m) { return TransposeSum(TranspMult(m, x.f())); });
  return Derivative<M_Matrix<T>>(std::move(f), x.x(), std::move(df));
}

template <class T>
Derivative<M_Matrix<T>> Transpose(const Derivative<M_Matrix<T>> &x) {
  assert(x.f().ncols() == x.f().nrows());
  auto f = Transpose(x.f());
  auto df = x.dfdx().apply([](auto &x) { return Transpose(x); });
  return Derivative<M_Matrix<T>>(std::move(f), x.x(), std::move(df));
}

template <typename T>
Derivative<T> operator-(const Derivative<T> &x, const Derivative<T> &y) {
  assert(x.x() == y.x());
  return Derivative<T>(x.f() - y.f(), x.x(), x.dfdx() - y.dfdx());
}

template <typename T>
Derivative<T> operator*(const Derivative<T> &one, const Derivative<T> &other) {
  assert(one.x() == other.x());
  auto df1 = one.dfdx().apply([&other](auto &e) { return e * other.f(); });
  auto df2 = other.dfdx().apply([&one](auto &e) { return one.f() * e; });
  return Derivative<T>(one.f() * other.f(), one.x(), df1 + df2);
}

inline Derivative<M_Matrix<double>> operator*(const Derivative<double> &x,
                                       const Derivative<M_Matrix<double>> &M) {
  assert(x.x() == M.x());
  auto df = x.f() * M.dfdx();
  for (std::size_t i = 0; i < df.size(); ++i)
    df[i] += x.dfdx()[i] * M.f();
  return Derivative<M_Matrix<double>>(x.f() * M.f(), x.x(), std::move(df));
}

inline Derivative<M_Matrix<double>> operator*(const Derivative<M_Matrix<double>> &M,
                                       const Derivative<double> &x) {
  assert(x.x() == M.x());
  auto df = M.dfdx() * x.f();
  for (std::size_t i = 0; i < df.size(); ++i)
    df[i] += M.f() * x.dfdx()[i];
  return Derivative<M_Matrix<double>>(M.f() * x.f(), x.x(), std::move(df));
}

template <class T>
myOptional_t<Derivative<M_Matrix<T>>> inv(const Derivative<M_Matrix<T>> &x) {
  typedef myOptional_t<Derivative<M_Matrix<T>>> Op;

  auto invx = matrix_inverse::Matrix_inverse(x.f());
  if (!invx.has_value())
    return Op(false, "cannot invert x " + invx.error());
  else {
    //        auto fxxf=Permute<T,1,3,0,2>(x.dfdx());
    auto df = -invx.value().first * x.dfdx() * invx.value().first;
    //      df=Permute<T,2,0,3,1>(df);
    return Op(
        Derivative<M_Matrix<T>>(std::move(invx).value().first, x.x(), std::move(df)));
  }
}

inline
    Derivative<M_Matrix<double>> expm_taylor(const Derivative<M_Matrix<double>> & x, std::size_t order=6)
{
    Derivative<M_Matrix<double>> out=Constant(Matrix_Generators::eye<double>(x.f().ncols()))+x;
    Derivative<M_Matrix<double>> xr=x;
    double a=1.0;
    for (std::size_t n=2; n+1<order; ++n)
    {
        a/=n;
        xr=xr*x;
        out=out+xr*a;
    }
    return out;

}


namespace Matrix_Decompositions {


inline auto EigenSystem_full_real_eigenvalue(const Derivative<M_Matrix<double>> &Dx) {
  typedef myOptional_t<
      std::tuple<Derivative<M_Matrix<double>>, Derivative<M_Matrix<double>>,
                   Derivative<M_Matrix<double>>, M_Matrix<double>,M_Matrix<double>>>
      Op;
  auto res = Matrix_Decompositions::EigenSystem_full_real_eigenvalue(Dx.f());
  if (!res)
    return Op(false,
              "Derivative error cannot calculete function value" + res.error());
  else {
    auto [VR, landa, VL, CL,CV] = std::move(res).value();
    // auto fxxf=Permute<double,1,3,0,2>(Dx.dfdx());

    M_Matrix<Derivative<double>> dlanda(
        landa.nrows(), landa.ncols(), landa.type(), Derivative<double>(Dx.x()));

    // auto n = landa.size();
    for (std::size_t i = 0; i < landa.size(); ++i) {
      auto vT = VL(i, ":");
      auto u = VR(":", i);
      dlanda[i].f() = landa[i];
      dlanda[i].dfdx() = Dx.dfdx().auto_apply(
          [&vT, &u](auto &m) { return (vT * m * u).getvalue(); });
    }

    M_Matrix<M_Matrix<double>> C(
        Dx.x().nrows(), Dx.x().ncols(), Dx.x().type(),
        M_Matrix<double>(VR.nrows(), VR.ncols(), VR.type()));
    for (std::size_t k = 0; k < VR.nrows(); ++k) {
      auto uk = VR(":", k);
      std::size_t m = 0;
      for (std::size_t j = 0; j < VR.ncols(); ++j) {
        if (are_Equal<false, double>().test(uk[j], 1))
          m = j;
        if (k != j) {
          auto vTj = VL(j, ":");
          double dl = landa[k] - landa[j];
          for (std::size_t is = 0; is < Dx.x().size(); ++is)
            C[is](k, j) = (vTj * Dx.dfdx()[is] * uk).getvalue() / dl;
        }
      }
      for (std::size_t is = 0; is < Dx.x().size(); ++is) {
        C[is](k, k) = 0;
        for (std::size_t j = 0; j < VR.ncols(); ++j) {
          if (k != j)
            C[is](k, k) -= VR(m, j) * C[is](k, j);
        }
      }
    }
    auto VRR = VR;
    auto dVR = C.apply([&VRR](auto &m) { return multTransp(VRR, m); });
    Derivative<M_Matrix<double>> DVR(std::move(VR), Dx.x(), std::move(dVR));
    auto DVL = inv(DVR);
    if (!DVL)
      return Op(false, " fails to invert the left eigenvector");

    return Op(std::tuple(std::move(DVR), Derivative<M_Matrix<double>>(dlanda),
                         std::move(DVL).value(), std::move(CL), std::move(CV)));
  }
}
} // namespace Matrix_Decompositions

template <> struct myDerivative<Matrix_Decompositions::eigensystem_type> {
  typedef std::tuple<Derivative<M_Matrix<double>>, Derivative<M_Matrix<double>>,
                       Derivative<M_Matrix<double>>,M_Matrix<double>,M_Matrix<double>>
      type;
};
template <> struct Der<Matrix_Decompositions::eigensystem_type> {
    typedef std::tuple<Derivative<M_Matrix<double>>, Derivative<M_Matrix<double>>,
                       Derivative<M_Matrix<double>>,M_Matrix<double>,M_Matrix<double>>
        type;
};

template <bool output>
struct are_Equal<output,
                 Derivative_t<Matrix_Decompositions::eigensystem_type>> {
public:
  bool test(const Derivative_t<Matrix_Decompositions::eigensystem_type> &one,
            const Derivative_t<Matrix_Decompositions::eigensystem_type> &other,
            std::ostream &os) {
    are_Equal<output, Derivative<M_Matrix<double>>> t;
    std::stringstream ss_VR;
    bool test_VR = t.test_prod(std::get<0>(one), std::get<0>(other), ss_VR);
    if (!test_VR)
      os << "\n\n\n errror in VR-------!!!!\n\n\n"
         << ss_VR.str() << "\n\n\n errror in VR--END---!!!!\n\n\n";

    std::stringstream ss_L;
    bool test_L = t.test_prod(std::get<1>(one), std::get<1>(other), ss_L);
    if (!test_L)
      os << "\n\n\n errror in Landa-------!!!!\n\n\n"
         << ss_L.str() << "\n\n\n errror in Landa--END---!!!!\n\n\n";

    std::stringstream ss_VL;
    bool test_VL = t.test_prod(std::get<2>(one), std::get<2>(other), ss_VL);
    if (!test_VL)
      os << "\n\n\n errror in VLeft-------!!!!\n\n\n"
         << ss_VL.str() << "\n\n\n errror in VLeft--END---!!!!\n\n\n";

    return test_VR && test_L && test_VL;
  }
};

template <typename T>
Derivative<M_Matrix<T>> diag(const Derivative<M_Matrix<T>> &x) {
  return Derivative<M_Matrix<T>>(
      Matrix_Unary_Transformations::diag(x.f()), x.x(),
      x.dfdx().apply([](auto &df) { return diag(df); }));
}

#endif // MATRIXDERIVATIVE_H
