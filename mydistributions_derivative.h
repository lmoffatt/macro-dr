#ifndef MYDISTRIBUTIONS_DERIVATIVE_H
#define MYDISTRIBUTIONS_DERIVATIVE_H

#include "myDistributions.h"
#include "matrixderivative.h"

template <>
struct Derivative<Probability_distribution> : public Probability_distribution {
  typedef Probability_distribution base_type;

  template <bool output, class dP>
  static myOptional_t<void> testDerivative(const dP &p, double tolerance) {
    for (std::size_t j = 0; j < p.dfdx().size(); ++j) {
      double sum = 0;
      std::stringstream ss;
      for (std::size_t i = 0; i < p.f().size(); ++i)
        sum += p.dfdx()[j][i];
      if (std::abs(sum) >= tolerance)
      {
        if constexpr (output)
          return Op_void(false,
                         " p=\n" + ToString(p) +
                             " sum derivative on variable j= " + ToString(j) +
                             " sum=" + ToString(sum) + "\n");
        else
          return Op_void(false,"");

      }
      else
      return Op_void(true, "");
    }
    return Op_void(true, "");
  }

  template <bool output, class dP>
  static Op_void test(const dP &p, double tolerance) {
    auto test_f = base_type::test<output>(p.f(), tolerance);
    auto test_dfdx = testDerivative<output>(p, tolerance);

    return std::move(test_f) + std::move(test_dfdx);
  }

  template <class dP> dP static normalize_derivative(dP &&dp) {
    double dpos = 0;
    double dneg = 0;
    for (std::size_t i = 0; i < dp.size(); ++i) {
      if (dp[i] > 0)
        dpos += dp[i];
      else
        dneg -= dp[i];
    }
    if (dpos + dneg > 0) {
      double fpos = 2.0*dneg / (dpos + dneg);
      double fneg = 2.0*dpos / (dpos + dneg);
      for (std::size_t i = 0; i < dp.size(); ++i) {
        if (dp[i] > 0)
          dp[i] *= fpos;
        else
          dp[i] *= fneg;
      }
    }
    return dp;
  }

  template <class P>
  Derivative<P> static normalize(Derivative<P> &&p, double min_p) {
    double sum = 0;
    for (std::size_t i = 0; i < p.f().size(); ++i) {
      if (p.f()[i] < min_p) {
        p.f()[i] = 0;
        p.dfdx().transform([&i](auto &m) { m[i] = 0; });
      } else {
        if (p.f()[i] + min_p > 1) {
          p.f()[i] = 1;
          p.dfdx().transform([&i](auto &m) { m[i] = 0; });
        }
        sum += p.f()[i];
      }
    }
    for (std::size_t i = 0; i < p.f().size(); ++i)
      p.f()[i] = p.f()[i] / sum;
    p.dfdx().transform([](auto &m) { m = normalize_derivative(std::move(m)); });
    return p;
  }
};

template <>
struct Derivative<Probability_distribution_covariance>
    : public Probability_distribution_covariance {

  typedef Probability_distribution_covariance base_type;

  static M_Matrix<double>
  normalize_derivative(const M_Matrix<double> &dp0,
                       const std::set<std::size_t> &non_zero_i0) {
    M_Matrix<double> dp=dp0;
    auto non_zero_i=non_zero_i0;
    for (auto i : non_zero_i) {

      double dpos = 0;
      double dneg = 0;
      for (auto j : non_zero_i) {
        if (dp(i, j) > 0)
          dpos += dp(i, j);
        else
          dneg -= dp(i, j);
      }
      if ((dpos + dneg) > 0) {
        double fpos = 2.0*dneg / (dpos + dneg);
        double fneg = 2.0*dpos / (dpos + dneg);
        for (auto j : non_zero_i) {
          if (dp(i, j) > 0)
            dp(i, j) *= fpos;
          else
            dp(i, j) *= fneg;
        }
      }
    }
    return dp;
  }

  static Derivative<M_Matrix<double>>
  normalize(Derivative<M_Matrix<double>> &&p, double min_p) {
    std::set<std::size_t> non_zero_i;

    for (std::size_t i = 0; i < p.f().nrows(); ++i) {
      if (p.f()(i, i) < min_p) {
        for (std::size_t j = 0; j < p.f().ncols(); ++j) {
          p.f()(i, j) = 0;
          p.dfdx().transform([&i, &j](auto &m) { m(i, j) = 0; });
        }
      } else if (p.f()(i, i) + min_p > 1.0) {
        for (std::size_t n = 0; n < p.f().size(); ++n) {
          p.f()[n] = 0;
          p.dfdx().transform([&n](auto &m) { m[n] = 0; });
        }
        p.f()(i, i) = 1;

        return p;
      } else
        non_zero_i.insert(i);
    }
    for (auto i : non_zero_i) {

      double sum = 0;
      for (auto j : non_zero_i)
        if (j != i)
          sum += p.f()(i, j);

      if (-sum != p.f()(i, i)) {
        auto sum_new = (sum - p.f()(i, i)) * 0.5;
        double f = sum_new / sum;
        p.f()(i, i) = -sum_new;
        for (auto j : non_zero_i)
          if (i != j)
            p.f()(i, j) *= f;
      }
    }
    p.dfdx().transform([&non_zero_i](auto &m) {
      m = normalize_derivative(std::move(m), non_zero_i);
    });

    return p;
  }

  template <bool output>
  static Op_void test_derivative(const M_Matrix<M_Matrix<double>> &dp,
                                 double tolerance) {
    for (std::size_t npar = 0; npar < dp.size(); ++npar) {
      auto &p = dp[npar];
      for (std::size_t i = 0; i < p.nrows(); ++i) {
        double sum = 0;
        for (std::size_t j = 0; j < p.ncols(); ++j) {

          if (i != j)
            sum += p(i, j);
        }
        if (std::abs(p(i, i) + sum) > tolerance) {
          if constexpr (output) {
            std::stringstream ss;
            ss << "tolerance=" << tolerance << "\n";
            ss << " dp[" + ToString(npar) + "]=\n"
               << p << "\n i=" << i << " dp(i,j)=" << p(i, i) << " sum=" << sum
               << "\n";
            return Op_void(false, ss.str());
          }
          return Op_void(false, "");
        }
      }
    }
    return Op_void(true, "");
  }

  template <bool output, class dP>
  static Op_void test(const dP &p, double tolerance) {
    auto test_f = base_type::test<output>(p.f(), tolerance);
    auto test_df = test_derivative<output>(p.dfdx(), tolerance);
    return std::move(test_f) + std::move(test_df);
  }
};

template <> struct Derivative<variance_value> : public variable_value {

  static Derivative<double> adjust(Derivative<double> &&value,
                                   double min_variance) {
    value.f() = std::max(value.f(), min_variance);
    return value;
  }
};

template <typename T> class Derivative<Base_Distribution<T>> {
public:
  virtual Derivative<Base_Distribution<T>> *clone() const = 0;

  constexpr static auto const className =
      my_static_string("Base_Distribution_derivative_") +
      my_trait<T>::className;
  virtual std::string myClass() const = 0;

  virtual T sample(std::mt19937_64 &mt) const = 0;

  virtual Derivative<double> p(const T &x) const = 0;

  virtual Derivative<double> logP(const T &x) const = 0;

  virtual M_Matrix<Derivative<T>> const &param() const = 0;

  // virtual M_Matrix<T> Score(const T& x)const =0;    // V(x)

  virtual M_Matrix<T> Fisher_Information() const = 0; // I()

  virtual Derivative<T> dlogL_dx(const T &x) const = 0;

  virtual Derivative<T> dlogL_dx2(const T &x) const = 0;

  virtual Derivative<T> mean() const = 0;

  virtual Derivative<T> stddev() const = 0;

  M_Matrix<double> const &x() const { return mean().x(); }

  virtual Derivative<double> expected_logP() const { return {0, mean().x()}; }
  virtual Derivative<double> variance_logP() const { return {0.5, mean().x()}; }
  virtual M_Matrix<double> FIM() const = 0;

  virtual ~Derivative<Base_Distribution<T>>() {}
};

template <>
class Derivative<Normal_Distribution<double>>
    : public Derivative<Base_Distribution<double>> {
public:
  typedef Derivative<Base_Distribution<double>> base_type;
  typedef Normal_Distribution<double> primitive_type;

  constexpr static auto const className =
      my_static_string("Normal_Distribution_derivative");
  std::string myClass() const override { return className.str(); }

  typedef Normal_Distribution<double> self_type;
  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "mean", &self_type::mean),
        grammar::field(C<self_type>{}, "variance", &self_type::variance));
  }

  virtual Derivative<Normal_Distribution<double>> *clone() const override {
    return new Derivative(*this);
  };

  virtual double sample(std::mt19937_64 &mt) const override {
    return std::normal_distribution<double>{mean().f(),
                                            std::sqrt(variance().f())}(mt);
  }

  virtual Derivative<double> p(const double &x) const override {
    return 1.0 / sqrt(2.0 * PI * variance()) *
           exp(-0.5 * sqr(Constant(x) - mean()) / variance());
  }

  virtual Derivative<double> logP(const double &x) const override {
    return Constant(-0.5 * std::log(2 * PI)) - 0.5 * log(variance()) -
           0.5 * sqr(Constant(x) - mean()) / variance();
  }

  /*
  virtual M_Matrix<double> Score(const double& x) const override
  { return
  M_Matrix<double>(1,2,std::vector<double>{dlogL_dmean(x),dlogL_dstddev(x)});};
*/
  virtual M_Matrix<double> Fisher_Information() const override {
    return M_Matrix<double>(
        2, 2, Matrix_TYPE::DIAGONAL,
        std::vector<double>{d2logL_dmean2(), d2logL_dvariance2()});
  };

  virtual Derivative<double> mean() const override { return param_[0]; }

  virtual Derivative<double> stddev() const override {
    return sqrt(variance());
  };

  virtual Derivative<double> variance() const { return param_[1]; }

  virtual M_Matrix<Derivative<double>> const &param() const override {
    return param_;
  }

  Derivative(Derivative<double> mean, Derivative<double> variance)
      : param_(1, 2, std::vector<Derivative<double>>{mean, variance}) {}

  Derivative() = default;
  virtual ~Derivative() {}
  virtual Derivative<double> dlogL_dx(const double &x) const override {
    return (mean() + Constant(-x)) / variance();
  }

  virtual Derivative<double> dlogL_dx2(const double &) const override {
    return -1.0 / variance();
  }

  virtual Derivative<double> expected_logP() const override {
    return Constant(-0.5) - log((2.0 * PI) * variance()) * 0.5;
  }

  virtual Derivative<double> variance_logP() const override {
    return {0.5, mean().x()};
  }

  Derivative(const Derivative &) = default;
  Derivative(Derivative &&) = default;
  Derivative &operator=(const Derivative &) = default;
  Derivative &operator=(Derivative &&) = default;

  M_Matrix<double> const &x() { return mean().x(); }

  virtual M_Matrix<double> FIM() const override {
    if (mean().dfdx().nrows() > mean().dfdx().ncols())
      return quadraticForm_XXT(mean().dfdx()) * d2logL_dmean2() +
             quadraticForm_XXT(variance().dfdx()) * d2logL_dvariance2();
    else
      return quadraticForm_XTX(mean().dfdx()) * d2logL_dmean2() +
             quadraticForm_XTX(variance().dfdx()) * d2logL_dvariance2();
  }

protected:
  M_Matrix<Derivative<double>> param_;

  double d2logL_dmean2() const { return 1.0 / variance().f(); }

  double d2logL_dvariance2() const { return 0.5 / sqr(variance().f()); }
};

#endif // MYDISTRIBUTIONS_DERIVATIVE_H
