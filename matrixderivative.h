#ifndef MATRIXDERIVATIVE_H
#define MATRIXDERIVATIVE_H
#include "Matrix.h"
#include "myfields.h"

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
  double const &f() const { return f_; }
  double &f() { return f_; }
  M_Matrix<double> const &x() const { return *x_; }
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
  Derivative(const M_Matrix<double> &myx)
      : f_{0}, x_{&myx}, dfdx_{myx.nrows(), myx.ncols(), myx.type(), 0.0} {}
  Derivative &operator+=(const Derivative &other) {
    assert(x() == other.x());
    f_ += other.f();
    dfdx_ += other.dfdx();
    return *this;
  }
  Derivative &operator-=(const Derivative &other) {
    assert(x() == other.x());
    f_ -= other.f();
    dfdx_ -= other.dfdx();
    return *this;
  }
};

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
  M_Matrix<M_Matrix<double>> const &dfdx() const { return dfdx_; }
  M_Matrix<M_Matrix<double>> &dfdx() { return dfdx_; }

  Derivative(const M_Matrix<double> &fx, const M_Matrix<double> &myx,
             const M_Matrix<M_Matrix<double>> &der)
      : f_{fx}, x_{&myx}, dfdx_{der} {}

  Derivative(const M_Matrix<double> &fx, const M_Matrix<double> &myx)
      : f_{fx}, x_{&myx}, dfdx_{M_Matrix<M_Matrix<double>>(
                              myx.nrows(), myx.ncols(), myx.type(),
                              M_Matrix<double>(fx.nrows(), fx.ncols(),
                                               fx.type()))} {}

  template <class M_f, class M_d>
  Derivative(M_f &&fx, const M_Matrix<double> &myx, M_d &&der)
      : f_{std::forward<M_f>(fx)}, x_{&myx}, dfdx_{std::forward<M_d>(der)} {}

  Derivative(std::size_t nrows, std::size_t ncols, M_Matrix<double>::TYPE t,
             const M_Matrix<double> &x)
      : f_{nrows, ncols, t}, x_{&x}, dfdx_{x.nrows(), x.ncols(), x.type(),
                                           M_Matrix<double>(nrows, ncols, t)} {}
  Derivative() = default;

  Derivative(const M_Matrix<Derivative<double>> &x)
      : f_{x.nrows(), x.ncols(), x.type()}, x_{&x[0].x()},
        dfdx_{x[0].x().nrows(), x[0].x().ncols(), x[0].x().type(),
              M_Matrix<double>(x.nrows(), x.ncols(), x.type())} {
    for (std::size_t i = 0; i < x.size(); ++i) {
      f_[i] = x[i].f();
      for (std::size_t j = 0; j < x[0].x().size(); ++j)
        dfdx_[j][i] = x[i].dfdx()[j];
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
      return Derivative<double>(f().getvalue(), x(), std::move(dfdx()[0]));
    }
  }
};

template <bool output, typename T>
class are_Equal<output, Derivative<M_Matrix<T>>> {
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

template <class F, typename... Ts>
Derivative<M_Matrix<double>>
Incremental_ratio(double eps, const F &fun,
                  const Derivative<M_Matrix<double>> &y0,
                  const Derivative<M_Matrix<Ts>> &... y) {
  assert(((y0.x() == y.x()) && ...));
  M_Matrix<double> f = fun(y0.f(), y.f()...);
  auto &one = y0;
  M_Matrix<M_Matrix<double>> dfdx(one.x().nrows(), one.x().ncols(),
                                  one.x().type());
  for (std::size_t i = 0; i < one.x().size(); ++i) {
    double e = eps;
    if (one.x()[i] != 0)
      e *= std::abs(one.x()[i]);
    auto fpos = fun(Taylor_first(y0, i, e), Taylor_first(y, i, e)...);
    auto fneg = fun(Taylor_first(y0, i, -e), Taylor_first(y, i, -e)...);
    dfdx[i] = (fpos - fneg) / (2.0 * e);
  }
  return Derivative<M_Matrix<double>>(f, one.x(), dfdx);
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
auto Primitive(const Derivative<C> &dy)
    ->std::enable_if_t<is_field_Object<C>::value&&is_field_Object<Derivative<C>>::value,C>
{
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
auto Incremental_ratio_calc(const std::vector<double> &eps, const C &f,
                            const M_Matrix<double> &x,
                            const std::vector<C> &fpos,
                            const std::vector<C> &fneg)
    -> Derivative<std::decay_t<decltype(C::get_constructor_fields(),
                                        std::declval<C>())>> {
  auto fields = f.get_constructor_fields();
  auto dfields = Derivative<std::decay_t<C>>::get_constructor_fields();
  return myApply(
      [&f, &x, &fpos, &fneg, &eps](auto &&... mdm) {
        return Derivative<std::decay_t<C>>(Incremental_ratio_method_calc(
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
Incremental_ratio_calc(const std::vector<double> &eps, const M_Matrix<double> &f,
                       const M_Matrix<double> &x,
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
    const grammar::field<Derivative<C>, dmethod> dm, const C &f,
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

auto exp(D, double x) { return std::exp(x); }

auto exp(const Derivative<double> &x) {
  return Derivative<double>(std::exp(x.f()), x.x(), exp(D(), x.f()) * x.dfdx());
}

auto log(const Derivative<double> &x) {
  return Derivative<double>(std::log(x.f()), x.x(), x.dfdx() / x.f());
}

auto sqrt(const Derivative<double> &x) {
  return Derivative<double>(std::sqrt(x.f()), x.x(),
                            x.dfdx() / std::sqrt(x.f()));
}

auto exp(const Derivative<M_Matrix<double>> &x) {
  return Derivative<M_Matrix<double>>(exp(x.f()), x.x(),
                                      elemMult_a(x.dfdx(), exp(x.f())));
}

auto abs(D, double x) {
  if (x >= 0)
    return 1.0;
  else
    return -1.0;
}

auto abs(const Derivative<double> &x) {
  return Derivative<double>(std::abs(x.f()), x.x(), abs(D(), x.f()) * x.dfdx());
}

auto pow(const Derivative<double> &x, const Derivative<double> &y) {
  assert(x.x() == y.x());
  auto f = std::pow(x.f(), y.f());

  auto dx = x.dfdx() * (y.f() / x.f() * f);
  auto dy = y.dfdx() * (std::log(x.f()) * f);
  return Derivative<double>(std::move(f), x.x(), dx + dy);
}

auto operator/(const Derivative<double> &x, const Derivative<double> &y) {
  assert(x.x() == y.x());
  auto f = x.f() / y.f();

  auto dx = x.dfdx() / y.f();
  auto dy = -y.dfdx() * x.f() / sqr(y.f());
  return Derivative<double>(std::move(f), x.x(), dx + dy);
}

auto operator/(double x, const Derivative<double> &y) {
  auto f = x / y.f();

  auto dy = -y.dfdx() * x / sqr(y.f());
  return Derivative<double>(std::move(f), y.x(), std::move(dy));
}

template <class> struct is_Derivative : public std::false_type {};

template <typename T>
struct is_Derivative<Derivative<T>> : public std::true_type {};

template <class D>
constexpr static bool is_Derivative_v = is_Derivative<D>::value;

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
  auto dfdB = B.dfdx().apply([&A, &B](auto &m) {
    return TranspMult(m, A.f() * B.f()) + TranspMult(B.f(), A.f() * m);
  });
  return Derivative<M_Matrix<T>>(std::move(f), A.x(), dfdA + dfdB);
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
  auto dfdx = elemDiv_a(x.dfdx(), y.f());
  auto dfdy = elemMult_a(f, elemDiv_a(y.dfdx(), y.f()));

  return Derivative<M_Matrix<T>>(std::move(f), x.x(), dfdx + dfdy);
}

template <typename T>
auto elemDivSafe(const Derivative<M_Matrix<T>> &x,
                 const Derivative<M_Matrix<T>> &y,
                 double eps = std::numeric_limits<double>::epsilon()) {
  assert(x.x() == y.x());
  auto f = elemDivSafe(x.f(), y.f(), eps);
  auto dfdx = elemDivSafe_a(x.dfdx(), y.f(), eps);
  auto dfdy = elemMult_a(elemDivSafe_a(y.dfdx(), y.f(), eps), f);

  return Derivative<M_Matrix<T>>(std::move(f), x.x(), dfdx + dfdy);
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
      [&x](auto &m) { return TranspMult(m, x.f()) + TranspMult(x.f(), m); });
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

Derivative<M_Matrix<double>> operator*(const Derivative<double> &x,
                                       const Derivative<M_Matrix<double>> &M) {
  assert(x.x() == M.x());
  auto df = x.f() * M.dfdx();
  for (std::size_t i = 0; i < df.size(); ++i)
    df[i] += x.dfdx()[i] * M.f();
  return Derivative<M_Matrix<double>>(x.f() * M.f(), x.x(), std::move(df));
}

Derivative<M_Matrix<double>> operator*(const Derivative<M_Matrix<double>> &M,
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
    auto df = -invx.value() * x.dfdx() * invx.value();
    //      df=Permute<T,2,0,3,1>(df);
    return Op(
        Derivative<M_Matrix<T>>(std::move(invx).value(), x.x(), std::move(df)));
  }
}
namespace Matrix_Decompositions {

auto EigenSystem_full_real_eigenvalue(const Derivative<M_Matrix<double>> &Dx) {
  typedef myOptional_t<
      std::tuple<Derivative<M_Matrix<double>>, Derivative<M_Matrix<double>>,
                 Derivative<M_Matrix<double>>>>
      Op;
  auto res = Matrix_Decompositions::EigenSystem_full_real_eigenvalue(Dx.f());
  if (!res)
    return Op(false,
              "Derivative error cannot calculete function value" + res.error());
  else {
    auto [VR, landa, VL] = std::move(res).value();
    // auto fxxf=Permute<double,1,3,0,2>(Dx.dfdx());

    M_Matrix<Derivative<double>> dlanda(
        landa.nrows(), landa.ncols(), landa.type(), Derivative<double>(Dx.x()));

    auto n = landa.size();
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
                         std::move(DVL).value()));
  }
}
} // namespace Matrix_Decompositions

template <> struct myDerivative<Matrix_Decompositions::eigensystem_type> {
  typedef std::tuple<Derivative<M_Matrix<double>>, Derivative<M_Matrix<double>>,
                     Derivative<M_Matrix<double>>>
      type;
};

template <bool output>
class are_Equal<output, Derivative_t<Matrix_Decompositions::eigensystem_type>> {
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
