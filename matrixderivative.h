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
      : f_{fx}, x_{&myx}, dfdx_{myx.nrows(), myx.ncols(), myx.type(), 0.0} {}
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

  template <class M_f, class M_d>
  Derivative(M_f &&fx, const M_Matrix<double> &myx, M_d &&der)
      : f_{std::forward<M_f>(fx)}, x_{&myx}, dfdx_{std::forward<M_d>(der)} {}

  Derivative(std::size_t nrows, std::size_t ncols, M_Matrix<double>::TYPE t,
             const M_Matrix<double> &x)
      : f_{nrows, ncols, t}, x_{&x}, dfdx_{x.nrows(), x.ncols(), x.type(),
                                           M_Matrix<double>(nrows, ncols, t)} {}
  Derivative() = default;

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
  return Derivative<T>(x.f() * c.value, x.x(), x.dfdx());
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
    auto [VL, landa, VR] = std::move(res).value();
    // auto fxxf=Permute<double,1,3,0,2>(Dx.dfdx());
    auto fxxf = Dx.dfdx();

    M_Matrix<Derivative<double>> dlanda(
        landa.nrows(), landa.ncols(), landa.type(), Derivative<double>(Dx.x()));

    auto n = landa.size();
    for (std::size_t i = 0; i < landa.size(); ++i) {
      auto vT = VL(i, ":");
      auto u = VR(":", i);
      dlanda[i].f() = landa[i];
      dlanda[i].dfdx() =
          Dx.dfdx().auto_apply([&vT, &u](auto &m) { return (vT * m * u).getvalue(); });
    }

    M_Matrix<M_Matrix<double>> C(
        Dx.x().nrows(), Dx.x().ncols(), Dx.x().type(),
        M_Matrix<double>(VR.nrows(), VR.ncols(), VR.type()));
    for (std::size_t k = 0; k < VR.nrows(); ++k)
    {
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
          C[is](k, k) -= VR(m, j) * C[is]( k,j);
        }
      }
    }
    auto VRR = VR;
    auto dVR = C.apply([&VRR](auto &m) { return m*VRR; });
    Derivative<M_Matrix<double>> DVR(std::move(VR), Dx.x(), std::move(dVR));
    auto DVL = inv(DVR);
    if (!DVL)
      return Op(false, " fails to invert the left eigenvector");

    return Op(std::tuple(std::move(DVL).value(),
                         Derivative<M_Matrix<double>>(dlanda), std::move(DVR)));
  }
}
} // namespace Matrix_Decompositions

template <> struct myDerivative<Matrix_Decompositions::eigensystem_type> {
  typedef std::tuple<Derivative<M_Matrix<double>>, Derivative<M_Matrix<double>>,
                     Derivative<M_Matrix<double>>>
      type;
};

template <typename T>
Derivative<M_Matrix<T>> diag(const Derivative<M_Matrix<T>> &x) {
  return Derivative<M_Matrix<T>>(
      Matrix_Unary_Transformations::diag(x.f()), x.x(),
      x.dfdx().apply([](auto &df) { return diag(df); }));
}

#endif // MATRIXDERIVATIVE_H
