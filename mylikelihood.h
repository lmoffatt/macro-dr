#ifndef MYLIKELIHOOD_H
#define MYLIKELIHOOD_H
#include "Matrix.h"
#include "myDistributions.h"
#include "myparameters.h"
#include "mydataframe.h"
#include <iomanip>

namespace evidence {

/*
 *
 *  supongamos en el medio que tenemos un Modelo que tira
 *
 *
 *   Modelo + vector de Parametros  +  vector de Datos -> vector de distribucion
 * de los parametros  y de los datos.
 *
 *
 *
 *
 * Como obtengo lo de arriba a partir de esto?
 * Posibilidad 1, a partir de derivada por diferencia
 *
 * Posibilidad 2 a partir de derivada explicita. (esta posibilidad es una
 * pesadilla para macrodr no tiene sentido en esta etapa)
 *
 *
 *
 * */

template <class E, class D>
std::tuple<double, double, double, double>
calculate_Likelihood(const Base_Distribution<E> &p, const D &data) {
  if (std::isfinite(data))
    return {p.logP(data), p.expected_logP(),
            sqr(p.logP(data) - p.expected_logP()), p.variance_logP()};
  else
    return {0.0, 0.0, 0.0, 0.0};
}

template <template <class...> class V, template <class> class Distribution,
          class T, class D>
std::tuple<double, double, double, double>
calculate_Likelihood(const V<Distribution<T>> &P, const D &data) {

  assert(P.size() == data.num_measurements());
  double logL = 0;
  double elogL = 0;
  double vlogL = 0;
  double evlogL = 0;
  for (std::size_t i = 0; i < data.num_measurements(); ++i) {
    auto [lik, elik, vlik, evlik] = calculate_Likelihood(P[i], data[i]);
    logL += lik;
    elogL += elik;
    vlogL += vlik;
    evlogL += evlik;
  }
  return {logL, elogL, vlogL, evlogL};
}

template <class E, class D>
double calculate_Gradient(const Base_Distribution<E> &dn,
                          const Base_Distribution<E> &dp, const D &data) {
  if (std::isfinite(data))
    return dp.logP(data) - dn.logP(data);
  else
    return 0;
}

template <class E, class D>
double calculate_vGradient(const Base_Distribution<E> &dn,
                           const Base_Distribution<E> &dp,
                           const Base_Distribution<E> &dn2,
                           const Base_Distribution<E> &dp2, const D &data) {
  if (std::isfinite(data))
    return (dp.logP(data) - dn.logP(data)) *
           ((dp2.logP(data) - dn2.logP(data)));
  else
    return 0;
}

template <template <class...> class V, template <class...> class Distributions,
          typename T, class D>
double calculate_Gradient(const V<Distributions<T>> &dn,
                          const V<Distributions<T>> &dp, const D &data) {
  assert(dn.size() == dp.size());
  double out = 0;
  for (std::size_t i = 0; i < dn.size(); ++i) {
    out += calculate_Gradient(dn[i], dp[i], data[i]);
  }
  return out;
}

template <template <class...> class V, template <class...> class Distributions,
          typename T, class D>
double calculate_vGradient(const V<Distributions<T>> &dn,
                           const V<Distributions<T>> &dp,
                           const V<Distributions<T>> &dn2,
                           const V<Distributions<T>> &dp2, const D &data) {
  assert(dn.size() == dp.size());
  assert(dn.size() == dp2.size());

  double out = 0;
  for (std::size_t i = 0; i < dn.size(); ++i) {
    out += calculate_vGradient(dn[i], dp[i], dn2[i], dp2[i], data[i]);
  }
  return out;
}

template <template <class> class Distributions, typename T, class D>
double calculate_Gradient(const Distributions<T> &dn,
                          const Distributions<T> &dp, const D &data,
                          double eps) {
  return calculate_Gradient(dn, dp, data) / eps;
}
template <template <class> class Distributions, typename T, class D>
double
calculate_vGradient(const Distributions<T> &dn, const Distributions<T> &dp,
                    const Distributions<T> &dn2, const Distributions<T> &dp2,
                    const D &data, double eps, double eps2) {
  return calculate_vGradient(dn, dp, dn2, dp2, data) / (eps * eps2);
}

template <template <class...> class Distributions, typename T, class D>
M_Matrix<double> calculate_Gradient(const Distributions<T> &d0,
                                    const std::vector<Distributions<T>> &dp,
                                    const D &data,
                                    const std::vector<double> &eps) {
  // typedef myOptional_t<M_Matrix<double>> Op;
  M_Matrix<double> out(1, dp.size());
  for (std::size_t i = 0; i < dp.size(); ++i) {
    out[i] = calculate_Gradient(d0, dp[i], data) / eps[i];
  }
  return out;
}

template <template <class...> class Distributions, typename T, class D>
M_Matrix<double> calculate_vGradient(const Distributions<T> &d0,
                                     const std::vector<Distributions<T>> &dp,
                                     const D &data,
                                     const std::vector<double> &eps) {
  // typedef myOptional_t<M_Matrix<double>> Op;
  M_Matrix<double> out(dp.size(), dp.size(), Matrix_TYPE::SYMMETRIC);
  for (std::size_t i = 0; i < dp.size(); ++i)
    for (std::size_t j = 0; j <= i; ++j) {
      out(i, j) =
          calculate_vGradient(d0, dp[i], d0, dp[j], data) / (eps[i] * eps[j]);
    }
  return out;
}

template <template <class...> class Distributions, typename T, class D>
M_Matrix<double>
calculate_Gradient_center(const std::vector<Distributions<T>> &dn,
                          const std::vector<Distributions<T>> &dp,
                          const D &data, const std::vector<double> &eps) {
  // typedef myOptional_t<M_Matrix<double>> Op;
  M_Matrix<double> out(1, dp.size());
  for (std::size_t i = 0; i < dp.size(); ++i) {
    out[i] = calculate_Gradient(dn[i], dp[i], data) / (eps[i]);
  }
  return out;
}

template <template <class...> class Distributions, typename T, class D>
M_Matrix<double>
calculate_vGradient_center(const std::vector<Distributions<T>> &dn,
                           const std::vector<Distributions<T>> &dp,
                           const D &data, const std::vector<double> &eps) {
  // typedef myOptional_t<M_Matrix<double>> Op;
  M_Matrix<double> out(dp.size(), dp.size(), Matrix_TYPE::SYMMETRIC);
  for (std::size_t i = 0; i < dp.size(); ++i)
    for (std::size_t j = 0; j <= i; ++j) {
      out(i, j) = calculate_vGradient(dn[i], dp[i], dn[j], dp[j], data) /
                  (eps[i] * eps[j]);
    }
  return out;
}

template <class E> auto getParameter(const Base_Distribution<E> &d) {
  return d.param();
}

template <class E>
auto getParameter(const std::unique_ptr<Base_Distribution<E>> &d) {
  return d->param();
}

template <class E> auto getParameter(const Base_Distribution<E> *d) {
  return d->param();
}

template <class E> auto getFIM(const Base_Distribution<E> &d) {
  return d.Fisher_Information();
}

template <class E> auto getFIM(const std::unique_ptr<Base_Distribution<E>> &d) {
  return d->Fisher_Information();
}

template <class E> auto getFIM(const Base_Distribution<E> *d) {
  return d->Fisher_Information();
}

template <class E> auto get_elogL(const Base_Distribution<E> &d) {
  return d.expected_logP();
}

template <class E>
auto get_elogL(const std::unique_ptr<Base_Distribution<E>> &d) {
  return d->expected_logP();
}

template <class E> auto get_elogL(const Base_Distribution<E> *d) {
  return d->expected_logP();
}

template <class Distributions>
auto calculate_Hessian(std::size_t i, const std::vector<Distributions> &d0,
                       const std::vector<std::vector<Distributions>> &d,
                       const std::vector<double> &eps) {
  auto k = d.size();
  // auto n=d0.size();
  //  assert(n==data.size());

  auto &d_i = d0[i];
  auto param = getParameter(d_i);
  //   os<<"\n -------------i="<<i<<"------------\n";
  //   os<<"\nparam "<<param;
  auto FIM = getFIM(d_i);
  //    os<<"\nFIM\n"<<FIM;
  //    os<<"eps="<<eps;
  auto npar = param.size();
  M_Matrix<double> J(npar, k);
  for (std::size_t j = 0; j < k; ++j) {
    auto dpar = (getParameter(d.at(j).at(i)) - param) * 1.0 / eps[j];
    for (std::size_t jj = 0; jj < npar; ++jj)
      J(jj, j) = dpar[jj];
  }
  auto H = quadraticForm_BT_A_B(FIM, J);
  return -H;
}

template <class Distributions>
auto calculate_Hessian_center(std::size_t i,
                              const std::vector<Distributions> &d0,
                              const std::vector<std::vector<Distributions>> &dn,
                              const std::vector<std::vector<Distributions>> &dp,
                              std::vector<double> eps) {
  auto k = dp.size();
  // auto n=d0.size();
  //  assert(n==data.size());

  auto &d_i = d0[i];
  auto param = getParameter(d_i);
  //   os<<"\n -------------i="<<i<<"------------\n";
  //   os<<"\nparam "<<param;
  auto FIM = getFIM(d_i);
  //    os<<"\nFIM\n"<<FIM;
  //    os<<"eps="<<eps;
  auto npar = param.size();
  M_Matrix<double> J(npar, k);
  for (std::size_t j = 0; j < k; ++j) {
    auto dpar = (getParameter(dp.at(j).at(i)) - getParameter(dn.at(j).at(i))) *
                1.0 / eps[j];
    for (std::size_t jj = 0; jj < npar; ++jj)
      J(jj, j) = dpar[jj];
  }
  auto H = quadraticForm_BT_A_B(FIM, J);
  return -H;
}

template <class Distributions, class D>
auto calculate_Hessian(const std::vector<Distributions> &d0,
                       const std::vector<std::vector<Distributions>> &d,
                       const D &, const std::vector<double> &eps,
                       std::ostream &) {
  auto k = d.size();
  auto n = d0.size();

  M_Matrix<double> out(k, k, Matrix_TYPE::SYMMETRIC, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    out += calculate_Hessian(i, d0, d, eps);
  }
  return out;
}

template <class Distributions, class D>
auto calculate_Hessian_center(const std::vector<Distributions> &d0,
                              const std::vector<std::vector<Distributions>> &dn,
                              const std::vector<std::vector<Distributions>> &dp,
                              const D &, std::vector<double> eps,
                              std::ostream &) {
  auto k = dp.size();
  auto n = d0.size();

  M_Matrix<double> out(k, k, Matrix_TYPE::SYMMETRIC, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    out += calculate_Hessian_center(i, d0, dn, dp, eps);
  }
  return out;
}

class logLikelihood {
public:
  // std::string myClass()const  { return className.str();}

  constexpr static auto const className = my_static_string("logLikelihood");

  typedef logLikelihood self_type;
  static auto get_constructor_fields() {
    double (logLikelihood::*myLogL)() const = &self_type::logL;
    return std::make_tuple(
        grammar::field(C<self_type>{}, "logL", myLogL),
        grammar::field(C<self_type>{}, "elogL", &self_type::elogL),
        grammar::field(C<self_type>{}, "vlogL", &self_type::vlogL));
  }

  double logL() const { return logL_; }
  double elogL() const { return elogL_; }
  double vlogL() const { return vlogL_; }
  double evlogL() const { return evlogL_; }

  double chilogL() const { return (logL() - elogL()) / std::sqrt(vlogL()); }

  double chi2logL() const { return sqr(logL() - elogL()) / vlogL(); }

  logLikelihood(double value, double elogL, double vlogL, double evlogL)
      : logL_{value}, elogL_{elogL}, vlogL_{vlogL}, evlogL_{evlogL} {}
  logLikelihood(std::tuple<double, double, double, double> logL)
      : logL_{std::get<0>(logL)}, elogL_{std::get<1>(logL)},
        vlogL_{std::get<2>(logL)}, evlogL_{std::get<3>(logL)} {}

  logLikelihood() : logL_{std::numeric_limits<double>::quiet_NaN()} {}
  operator bool() const { return std::isfinite(logL_); }
  void set_logL(double logLik, double elogL, double vlogL, double evlogL) {
    logL_ = logLik;
    elogL_ = elogL;
    vlogL_ = vlogL;
    evlogL_ = evlogL;
  }

  template <class DataFrame>
  static void insert_col(DataFrame &d, const std::string &pre) {
    d.insert_column(pre + "logL", C<double>{});
    d.insert_column(pre + "elogL", C<double>{});
    d.insert_column(pre + "vlogL", C<double>{});
    d.insert_column(pre + "evlogL", C<double>{});
    d.insert_column(pre + "chilogL", C<double>{});
  }
  auto data_row(std::size_t) const {
    return std::tuple(logL(), elogL(), vlogL(), evlogL(), chilogL());
  }
  std::size_t data_size() const { return 1; }

  logLikelihood &operator+=(const logLikelihood &other) {
    logL_ += other.logL_;
    elogL_ += other.elogL_;
    vlogL_ += other.vlogL_;
    evlogL_ += other.evlogL_;
    return *this;
  }

private:
  double logL_;
  double elogL_;
  double vlogL_;
  double evlogL_;
};
} // namespace evidence
template <> class moments<evidence::logLikelihood> {
public:
  // std::string myClass()const  { return className.str();}

  constexpr static auto const className =
      my_static_string("logLikelihood_moments");

  typedef moments<evidence::logLikelihood> self_type;
  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "logL", &self_type::logL),
        grammar::field(C<self_type>{}, "elogL", &self_type::elogL),
        grammar::field(C<self_type>{}, "vlogL", &self_type::vlogL),
        grammar::field(C<self_type>{}, "evlogL", &self_type::evlogL));
  }

  moments<double> logL() const { return logL_; }
  moments<double> elogL() const { return elogL_; }
  moments<double> vlogL() const { return vlogL_; }
  moments<double> evlogL() const { return evlogL_; }
  moments<double> chilogL() const { return chilogL_; }

  moments(const std::vector<evidence::logLikelihood> &l)
      : logL_{l, &evidence::logLikelihood::logL},
        elogL_{l, &evidence::logLikelihood::elogL},
        vlogL_{l, &evidence::logLikelihood::vlogL},
        evlogL_{l, &evidence::logLikelihood::evlogL},
        chilogL_{l, &evidence::logLikelihood::chilogL} {}

  template <class logL>
  moments(const std::vector<logL> &l)
      : logL_{l, &logL::logL}, elogL_{l, &logL::elogL}, vlogL_{l, &logL::vlogL},
        evlogL_{l, &logL::evlogL}, chilogL_{l, &logL::chilogL} {}

  template <class logL, class Op, typename... Ts>
  moments(const std::vector<logL> &l, const Op &u, Ts... t)
      : logL_{l, [&u, &t...](const logL &d) { return u(d, t...).logL(); }},
        elogL_{l, [&u, &t...](const logL &d) { return u(d, t...).elogL(); }},
        vlogL_{l, [&u, &t...](const logL &d) { return u(d, t...).vlogL(); }},
        evlogL_{l, [&u, &t...](const logL &d) { return u(d, t...).evlogL(); }},
        chilogL_{l,
                 [&u, &t...](const logL &d) { return u(d, t...).chilogL(); }} {}

  moments(moments<double> value, moments<double> elogL, moments<double> vlogL,
          moments<double> evlogL, moments<double> chilogL)
      : logL_{value}, elogL_{elogL}, vlogL_{vlogL}, evlogL_{evlogL},
        chilogL_{chilogL} {}

  moments() = default;
  operator bool() const { return std::isfinite(logL_.mean()); }
  void set_logL(moments<double> logLik, moments<double> elogL,
                moments<double> vlogL, moments<double> evlogL,
                moments<double> chilogL) {
    logL_ = logLik;
    elogL_ = elogL;
    vlogL_ = vlogL;
    evlogL_ = evlogL;
    chilogL_ = chilogL;
  }

  auto data_row(MOMENTS m) const {
    return std::tuple_cat(logL().data_row(m), elogL().data_row(m),
                          vlogL().data_row(m), evlogL().data_row(m),
                          chilogL().data_row(m));
  }

private:
  moments<double> logL_;
  moments<double> elogL_;
  moments<double> vlogL_;
  moments<double> evlogL_;

  moments<double> chilogL_;
};

namespace evidence {
class DlogLikelihood : public logLikelihood {
public:
  typedef logLikelihood base_type;

  constexpr static auto const className =
      my_static_string("DlogLikelihood_") + base_type::className;
  // std::string myClass()const  { return className.str();}

  typedef DlogLikelihood self_type;
  static auto get_constructor_fields() {
    double (logLikelihood::*myLogL)() const = &base_type::logL;
    M_Matrix<double> const &(self_type::*myG)() const = &self_type::G;
    const M_Matrix<double> &(self_type::*myH)() const = &self_type::H;

    return std::make_tuple(
        grammar::field(C<self_type>{}, "logL", myLogL),
        grammar::field(C<self_type>{}, "elogL", &self_type::elogL),
        grammar::field(C<self_type>{}, "vlogL", &self_type::vlogL),
        grammar::field(C<self_type>{}, "evlogL", &self_type::evlogL),
        grammar::field(C<self_type>{}, "Gradient", myG),
        grammar::field(C<self_type>{}, "Hessian", myH));
  }

  const M_Matrix<double> &G() const { return G_; }
  const M_Matrix<double> &H() const { return H_; }

  operator bool() const {
    return (G_.size() > 0) && (H_.size() > 0) && logLikelihood::operator bool();
  }
  void set_G(M_Matrix<double> &&gradient) { G_ = std::move(gradient); }
  void set_H(M_Matrix<double> &&hessian) {
    assert(hessian.isSymmetric());
    H_ = std::move(hessian);
  }

  DlogLikelihood(double logL, double elogL, double vlogL, double evlogL,
                 const M_Matrix<double> &Gradient,
                 const M_Matrix<double> &Hessian)
      : logLikelihood(logL, elogL, vlogL, evlogL), G_{Gradient}, H_{Hessian} {
    assert(H_.isSymmetric());
  }
  DlogLikelihood(std::tuple<double, double, double, double, M_Matrix<double>,
                            M_Matrix<double>> &&logL)
      : DlogLikelihood(std::get<0>(logL), std::get<1>(logL), std::get<2>(logL),
                       std::get<3>(logL), std::move(std::get<4>(logL)),
                       std::move(std::get<5>(logL))) {}

  DlogLikelihood(std::tuple<double, double, double, double> logL,
                 M_Matrix<double> &&Gradient, M_Matrix<double> &&Hessian)
      : DlogLikelihood(std::get<0>(logL), std::get<1>(logL), std::get<2>(logL),
                       std::get<3>(logL), std::move(Gradient),
                       std::move(Hessian)) {}
  DlogLikelihood(double logL, double elogL, double vlogL, double evlog,
                 M_Matrix<double> &&Gradient, M_Matrix<double> &&Hessian)
      : logLikelihood(logL, elogL, vlogL, evlog), G_{std::move(Gradient)},
        H_{std::move(Hessian)} {
    assert(H_.isSymmetric());
  }
  DlogLikelihood(double) : logLikelihood(0.0, 0.0, 0.0, 0.0), G_{}, H_{} {}

  DlogLikelihood() = default;

  DlogLikelihood &operator+=(const DlogLikelihood &other) {
    logLikelihood::operator+=(other);
    G_ += other.G_;
    H_ += other.H_;
    return *this;
  }

  template <class DataFrame>
  static void insert_col(DataFrame &d, const std::string pre) {
    base_type::insert_col(d, pre);
    d.insert_column(pre + "Gradient", C<double>{});
    d.insert_column(pre + "Hessian", C<double>{});
  }
  auto data_row(std::size_t i, std::size_t j) const {
    return std::tuple_cat(base_type::data_row(i), std::tuple(G_[i], H_(i, j)));
  }
  template <class Parameters_Distribution>
  auto data_row(std::size_t k, const Parameters_Distribution &prior) const {
    std::size_t n = G_.size();
    auto [i, j] = M_Matrix<double>::pos_to_ij_Symmetric(k, n);
    return std::tuple_cat(
        base_type::data_row(k),
        std::tuple(prior.name(j), prior.name(i), G_[j], H_(j - i, i + j)));
  }
  std::size_t data_size() const { return G_.size(); }
  std::size_t data_big_size() const { return H_.size(); }

private:
  M_Matrix<double> G_;
  M_Matrix<double> H_;
};

class DVlogLikelihood : public DlogLikelihood {
public:
  typedef DlogLikelihood base_type;

  constexpr static auto const className = my_static_string("DVlogLikelihood");
  // std::string myClass()const  { return className.str();}

  typedef DVlogLikelihood self_type;
  static auto get_constructor_fields() {
    double (logLikelihood::*myLogL)() const = &base_type::logL;
    M_Matrix<double> const &(self_type::*myG)() const = &self_type::G;
    const M_Matrix<double> &(self_type::*myH)() const = &self_type::H;

    return std::make_tuple(
        grammar::field(C<self_type>{}, "logL", myLogL),
        grammar::field(C<self_type>{}, "elogL", &self_type::elogL),
        grammar::field(C<self_type>{}, "vlogL", &self_type::vlogL),
        grammar::field(C<self_type>{}, "evlogL", &self_type::evlogL),
        grammar::field(C<self_type>{}, "Gradient", myG),
        grammar::field(C<self_type>{}, "Hessian", myH),
        grammar::field(C<self_type>{}, "varG", &self_type::vG),
        grammar::field(C<self_type>{}, "scalingMatrix", &self_type::scaleM),
        grammar::field(C<self_type>{}, "gradient_chi_value",
                       &self_type::gradient_chi_value),
        grammar::field(C<self_type>{}, "gradient_fim_value",
                       &self_type::gradient_fim_value));
  }
  static double calc_chi2_value(const M_Matrix<double> &myG,
                                const M_Matrix<double> &myH) {
    auto Hinv = pinv(myH);
    if (!Hinv)
      return std::numeric_limits<double>::quiet_NaN();
    else
      return xTSigmaX(myG, Hinv.value());
  }
  static myOptional_t<M_Matrix<double>> calc_Scale_Matrix(
      M_Matrix<double> const &vG, M_Matrix<double> const &H,
      double tol = std::sqrt(std::numeric_limits<double>::epsilon())) {
    typedef myOptional_t<M_Matrix<double>> Op;
    auto fim = Matrix_Decompositions::EigenSystem_symm_real_eigenvalue(-H);
    auto cov = Matrix_Decompositions::EigenSystem_symm_real_eigenvalue(vG);
    if (fim.has_value() && cov.has_value()) {
      auto [fVR, fL, fVL] = std::move(fim).value();
      auto [cVR, cL, cVL] = std::move(cov).value();
      auto sqrt_cL = cL.apply([](auto const &x) { return std::sqrt(x); });
      auto sqrtinv_fL = fL.apply([&tol](auto const &x) {
        if (x < tol)
          return 0.0;
        else
          return 1.0 / std::sqrt(x);
      });
      auto sM = cVR * sqrt_cL * sqrtinv_fL * fVL;

      assert((rank_diag(sqrtinv_fL) == fL.size()
                 ? (are_Equal_v(quadraticForm_B_A_BT(-H, sM), vG, std::cerr,
                                 "\ncL\n", cL, "\n sqrt_cL\n ", sqrt_cL,
                                 "\n fL\n", fL, "\n sqrtinv_fL \n", sqrtinv_fL))
                  : true));
      return Op(sM);
    } else
      return Op(false, "error calculating the scale Matrix: fim:" +
                           fim.error() + "  cov:" + cov.error());
  }

  bool has_Scale_Matrix() const { return hasScaleMatrix_; }
  void calc(double tol = std::sqrt(std::numeric_limits<double>::epsilon())) {
    auto sM = calc_Scale_Matrix(vG(), H(), tol);
    if (sM) {
      sM_ = std::move(sM).value();
      gradient_chi_value_ = calc_chi2_value(G(), vG());
      gradient_fim_value_ = calc_chi2_value(G(), -H());

      hasScaleMatrix_ = true;
    } else {
      hasScaleMatrix_ = false;
    }
  }

  const M_Matrix<double> &vG() const { return vG_; }
  const M_Matrix<double> &scaleM() const { return sM_; }
  double gradient_chi_value() const { return gradient_chi_value_; }
  double gradient_fim_value() const { return gradient_fim_value_; }

  void set_vG(M_Matrix<double> &&variance_of_gradient) {
    assert(variance_of_gradient.isSymmetric());
    vG_ = std::move(variance_of_gradient);
    hasScaleMatrix_ = false;
  }

  DVlogLikelihood(double logL, double elogL, double vlogL, double evlogL,
                  const M_Matrix<double> &Gradient,
                  const M_Matrix<double> &Hessian,
                  const M_Matrix<double> &var_Gradient,
                  const M_Matrix<double> &scale_Matrix,
                  double gradient_chi_value, double gradient_fim_value)
      : DlogLikelihood(logL, elogL, vlogL, evlogL, Gradient, Hessian),
        vG_{var_Gradient}, hasScaleMatrix_(true), sM_{scale_Matrix},
        gradient_chi_value_{gradient_chi_value}, gradient_fim_value_{
                                                     gradient_fim_value} {
    assert(vG_.isSymmetric());
  }
  DVlogLikelihood(bool do_calc, double logL, double elogL, double vlogL,
                  double evlogL, const M_Matrix<double> &Gradient,
                  const M_Matrix<double> &Hessian,
                  const M_Matrix<double> &var_Gradient)
      : DlogLikelihood(logL, elogL, vlogL, evlogL, Gradient, Hessian),
        vG_{var_Gradient}, hasScaleMatrix_(false), sM_{} {
    assert(vG_.isSymmetric());
    if (do_calc)
      calc();
  }
  DVlogLikelihood(bool do_calc, double logL, double elogL, double vlogL,
                  double evlogL, const M_Matrix<double> &&Gradient,
                  const M_Matrix<double> &&Hessian,
                  const M_Matrix<double> &&var_Gradient)
      : DlogLikelihood(logL, elogL, vlogL, evlogL, std::move(Gradient),
                       std::move(Hessian)),
        vG_{std::move(var_Gradient)}, hasScaleMatrix_(false), sM_{} {
    assert(vG_.isSymmetric());
    if (do_calc)
      calc();
  }
  DVlogLikelihood(double)
      : DlogLikelihood(0.0), vG_{}, hasScaleMatrix_(false), sM_{},
        gradient_chi_value_{}, gradient_fim_value_{} {}

  DVlogLikelihood(bool do_calc,
                  const std::tuple<double, double, double, double> &logL,
                  const M_Matrix<double> &&Gradient,
                  const M_Matrix<double> &&Hessian,
                  const M_Matrix<double> &&var_Gradient)
      : DVlogLikelihood(do_calc, std::get<0>(logL), std::get<1>(logL),
                        std::get<2>(logL), std::get<3>(logL),
                        std::move(Gradient), std::move(Hessian),
                        std::move(var_Gradient)) {}

  DVlogLikelihood() = default;

  DVlogLikelihood &operator+=(const DVlogLikelihood &other) {
    if (other.evlogL()>0)
    {
      DlogLikelihood::operator+=(other);
      vG_ += other.vG_;
      hasScaleMatrix_ = false;
    }
    return *this;

  }

  template <class DataFrame> static void insert_col(DataFrame &d) {
    base_type::insert_col(d, "");
    d.insert_column("Gradient_Variance", C<double>{});
    d.insert_column("Scaling_Matrix", C<double>{});
    d.insert_column("Gradient_chi2_value", C<double>());
    d.insert_column("Gradient_fim_value", C<double>());
  }
  auto data_row(std::size_t i, std::size_t j) const {
    return std::tuple_cat(base_type::data_row(i, j),
                          std::tuple(vG()(i, j), scaleM()(i, j),
                                     gradient_chi_value(),
                                     gradient_fim_value()));
  }

private:
  M_Matrix<double> vG_;
  bool hasScaleMatrix_;
  M_Matrix<double> sM_;
  double gradient_chi_value_;
  double gradient_fim_value_;
};

} // namespace evidence

template <>
class moments<evidence::DlogLikelihood>
    : public moments<evidence::logLikelihood> {
public:
  typedef moments<evidence::logLikelihood> base_type;

  constexpr static auto const className =
      my_static_string("moments_DlogLikelihood_") + base_type::className;
  // std::string myClass()const  { return className.str();}

  typedef moments self_type;
  static auto get_constructor_fields() {

    return std::tuple_cat(
        base_type::get_constructor_fields(),
        std::tuple(grammar::field(C<self_type>{}, "Gradient", &self_type::G),
                   grammar::field(C<self_type>{}, "Hessian", &self_type::H)));
  }

  const moments<M_Matrix<double>> &G() const { return G_; }
  const moments_matrix<double> &H() const { return H_; }
  operator bool() const {
    return (G().mean().size() > 0) && (H_.size() > 0) &&
           base_type::operator bool();
  }
  void set_G(moments<M_Matrix<double>> &&gradient) { G_ = std::move(gradient); }
  void set_H(moments_matrix<double> &&hessian) { H_ = std::move(hessian); }

  moments(const moments<double> &logL, moments<double> const &elogL,
          moments<double> const &vlogL, moments<double> const &evlogL,
          moments<double> const &chilogL,
          const moments<M_Matrix<double>> &Gradient,
          const moments_matrix<double> &Hessian)
      : base_type(logL, elogL, vlogL, evlogL, chilogL), G_{Gradient},
        H_{Hessian} {}

  moments(const std::vector<evidence::DlogLikelihood> &l)
      : base_type(l), G_{l, &evidence::DlogLikelihood::G},
        H_{l, &evidence::DlogLikelihood::H} {}

  template <class DlogL>
  moments(const std::vector<DlogL> &l)
      : base_type(l), G_{l, &DlogL::G}, H_{l, &DlogL::H} {}

  template <class logL, class Op, typename... Ts>
  moments(const std::vector<logL> &l, const Op &u, Ts... t)
      : base_type(l, u, t...), G_{l,
                                  [&u, &t...](const logL &d) {
                                    return u(d, t...).G();
                                  }},
        H_{l, [&u, &t...](const logL &d) { return u(d, t...).H(); }} {}

  moments() = default;

  auto data_row(MOMENTS m, std::size_t i, std::size_t j) const {
    return std::tuple_cat(base_type::data_row(m), G().data_row(m, i, j),
                          H().data_row(m, i, j));
  }

  std::size_t data_dim(MOMENTS m) const {
    return std::max(G().data_dim(m), H().data_dim(m));
  }
  std::size_t data_size() const { return G().data_size(); }

private:
  moments<M_Matrix<double>> G_;
  moments_matrix<double> H_;
};

template <>
class moments<evidence::DVlogLikelihood>
    : public moments<evidence::DlogLikelihood> {
public:
  typedef moments<evidence::DlogLikelihood> base_type;

  constexpr static auto const className =
      my_static_string("DVlogLikelihood_moments");
  // std::string myClass()const  { return className.str();}

  typedef moments<evidence::DVlogLikelihood> self_type;

  auto &vG() const { return vG_; }
  auto &scaleM() const { return sM_; }
  auto &gradient_chi_value() const { return gradient_chi_value_; }
  auto &gradient_fim_value() const { return gradient_fim_value_; }

  static auto get_constructor_fields() {

    return std::tuple_cat(
        base_type::get_constructor_fields(),
        std::tuple(
            grammar::field(C<self_type>{}, "varG", &self_type::vG),
            grammar::field(C<self_type>{}, "scalingMatrix", &self_type::scaleM),
            grammar::field(C<self_type>{}, "gradient_chi_value",
                           &self_type::gradient_chi_value),
            grammar::field(C<self_type>{}, "gradient_fim_value",
                           &self_type::gradient_fim_value))

                );
  }

  moments(const moments<double> &logL, const moments<double> &elogL,
          const moments<double> &vlogL, const moments<double> &evlogL,
          const moments<double> &chilogL,
          const moments<M_Matrix<double>> &Gradient,
          const moments_matrix<double> &Hessian,
          const moments_matrix<double> &var_Gradient,
          const moments_matrix<double> &scale_Matrix)
      : base_type(logL, elogL, vlogL, evlogL, chilogL, Gradient, Hessian),
        vG_{var_Gradient}, sM_{scale_Matrix} {}

  moments(const std::vector<evidence::DVlogLikelihood> &l)
      : base_type(l), vG_{l, &evidence::DVlogLikelihood::vG},
        sM_{l, &evidence::DVlogLikelihood::scaleM},
        gradient_chi_value_{l, &evidence::DVlogLikelihood::gradient_chi_value},
        gradient_fim_value_{l, &evidence::DVlogLikelihood::gradient_fim_value} {
  }

  template <class DlogL>
  moments(const std::vector<DlogL> &l)
      : base_type(l), vG_{l, &DlogL::vG}, sM_{l, &DlogL::scaleM},
        gradient_chi_value_{l, &DlogL::gradient_chi_value},
        gradient_fim_value_{l, &DlogL::gradient_fim_value}

  {}

  template <class logL, class Op, typename... Ts>
  moments(const std::vector<logL> &l, const Op &u, Ts... t)
      : base_type(l, u, t...), vG_{l,
                                   [&u, &t...](const logL &d) {
                                     return u(d, t...).vG();
                                   }},
        sM_{l, [&u, &t...](const logL &d) { return u(d, t...).scaleM(); }},
        gradient_chi_value_{l,
                            [&u, &t...](const logL &d) {
                              return u(d, t...).gradient_chi_value();
                            }},
        gradient_fim_value_{l,
                            [&u, &t...](const logL &d) {
                              return u(d, t...).gradient_fim_value();
                            }}

  {}

  moments() = default;

  auto data_row(MOMENTS m, std::size_t i, std::size_t j) const {
    return std::tuple_cat(base_type::data_row(m, i, j), vG().data_row(m, i, j),
                          scaleM().data_row(m, i, j),
                          gradient_chi_value().data_row(m),
                          gradient_fim_value().data_row(m));
  }

private:
  moments_matrix<double> vG_;
  moments_matrix<double> sM_;
  moments<double> gradient_chi_value_;
  moments<double> gradient_fim_value_;
};

namespace evidence {

template <class... Auxs> class PartialDLogLikelihood : public DVlogLikelihood {
public:
  typedef DVlogLikelihood base_type;

  constexpr static bool hasAux = sizeof...(Auxs) > 0;

  typedef std::tuple_element_t<0, std::tuple<Auxs...>> Aux;

  constexpr static auto const className =
      my_static_string("PartialDLogLikelihood");
  // std::string myClass()const  { return className.str();}

  typedef PartialDLogLikelihood self_type;
  static auto get_constructor_fields() {
    double (logLikelihood::*myLogL)() const = &base_type::logL;
    M_Matrix<double> const &(self_type::*myG)() const = &self_type::G;
    const M_Matrix<double> &(self_type::*myH)() const = &self_type::H;

    if constexpr (hasAux)
      return std::make_tuple(
          grammar::field(C<self_type>{}, "logL", myLogL),
          grammar::field(C<self_type>{}, "elogL", &self_type::elogL),
          grammar::field(C<self_type>{}, "vlogL", &self_type::vlogL),
          grammar::field(C<self_type>{}, "evlogL", &self_type::evlogL),
          grammar::field(C<self_type>{}, "Gradient", myG),
          grammar::field(C<self_type>{}, "Hessian", myH),
          grammar::field(C<self_type>{}, "varG", &self_type::vG),
          grammar::field(C<self_type>{}, "scalingMatrix", &self_type::scaleM),
          grammar::field(C<self_type>{}, "gradient_chi_value",
                         &self_type::gradient_chi_value),
          grammar::field(C<self_type>{}, "gradient_fim_value",
                         &self_type::gradient_fim_value),
          grammar::field(C<self_type>{}, "Partial_DlogL",
                         &self_type::partial_DlogL),
          grammar::field(C<self_type>{}, my_trait<Aux>::className.c_str(),
                         &self_type::getAuxiliar));
    else
      return std::make_tuple(
          grammar::field(C<self_type>{}, "logL", myLogL),
          grammar::field(C<self_type>{}, "elogL", &self_type::elogL),
          grammar::field(C<self_type>{}, "vlogL", &self_type::vlogL),
          grammar::field(C<self_type>{}, "evlogL", &self_type::evlogL),
          grammar::field(C<self_type>{}, "Gradient", myG),
          grammar::field(C<self_type>{}, "Hessian", myH),
          grammar::field(C<self_type>{}, "varG", &self_type::vG),
          grammar::field(C<self_type>{}, "scalingMatrix", &self_type::scaleM),
          grammar::field(C<self_type>{}, "gradient_chi_value",
                         &self_type::gradient_chi_value),
          grammar::field(C<self_type>{}, "gradient_fim_value",
                         &self_type::gradient_fim_value),
          grammar::field(C<self_type>{}, "Partial_DlogL",
                         &self_type::partial_DlogL));
  }

  template <class DataFrame> static void insert_col(DataFrame &d) {

    my_trait<Aux>::insert_col(d, "");
    base_type::insert_col(d);
    DlogLikelihood::insert_col(d, "partial_");
  }
  auto data_row(std::size_t isample, std::size_t i, std::size_t j) const {
    return std::tuple_cat(my_trait<Aux>::data_row(getAuxiliar()[isample]),
                          base_type::data_row(i, j),
                          partial_DlogL()[isample].data_row(i, j));
  }

  std::enable_if_t<hasAux, std::vector<Aux> const &> getAuxiliar() const {
    return aux_;
  }

  std::vector<DlogLikelihood> const &partial_DlogL() const { return partial_; }

  PartialDLogLikelihood(DVlogLikelihood &&dlogL,
                        std::vector<DlogLikelihood> &&dist,
                        std::vector<Aux> &&aux)
      : DVlogLikelihood(std::move(dlogL)), partial_{std::move(dist)},
        aux_{std::move(aux)} {}

  PartialDLogLikelihood(double logL, double elogL, double vlogL, double evlogL,
                        const M_Matrix<double> &Gradient,
                        const M_Matrix<double> &Hessian,
                        const M_Matrix<double> &var_Gradient,
                        const M_Matrix<double> &scale_Matrix,
                        double gradient_chi_value, double gradient_fim_value,
                        const std::vector<DlogLikelihood> &dist,
                        const std::vector<Aux> &aux)
      : DVlogLikelihood(logL, elogL, vlogL, evlogL, Gradient, Hessian,
                        var_Gradient, scale_Matrix, gradient_chi_value,
                        gradient_fim_value),
        partial_{dist}, aux_{aux} {}
  PartialDLogLikelihood(double logL, double elogL, double vlogL, double evlogL,
                        M_Matrix<double> &&Gradient, M_Matrix<double> &&Hessian,
                        M_Matrix<double> &&var_Gradient,
                        M_Matrix<double> &&scale_Matrix,
                        double gradient_chi_value, double gradient_fim_value,

                        std::vector<DlogLikelihood> &&dist,
                        std::vector<Aux> &&aux)
      : DVlogLikelihood(logL, elogL, vlogL, evlogL, std::move(Gradient),
                        std::move(Hessian), std::move(var_Gradient),
                        std::move(scale_Matrix), gradient_chi_value,
                        gradient_fim_value),
        partial_{std::move(dist)}, aux_{std::move(aux)} {
    calc();
  }
  PartialDLogLikelihood() = default;
  PartialDLogLikelihood(std::size_t n)
      : DVlogLikelihood(0.0), partial_{std::vector<DlogLikelihood>(n)},
        aux_{std::vector<Aux>(n)} {}

  PartialDLogLikelihood &add_one_i(std::size_t i, const DVlogLikelihood &other,
                                   Aux a) {

    DVlogLikelihood::operator+=(other);
    partial_[i] = other;
    aux_[i] = a;
    return *this;
  }
  PartialDLogLikelihood &add_one_i(std::size_t i, const DlogLikelihood &other) {
    DlogLikelihood::operator+=(other);
    partial_[i] = other;
  }

private:
  std::vector<DlogLikelihood> partial_;
  std::vector<Aux> aux_;
};

} // namespace evidence

template <typename Aux>
class moments<evidence::PartialDLogLikelihood<Aux>>
    : public moments<evidence::DVlogLikelihood> {
public:
  typedef moments<evidence::DVlogLikelihood> base_type;
  typedef evidence::PartialDLogLikelihood<Aux> element_type;

  constexpr static auto const className =
      my_static_string("moments_DlogLikelihood_");
  // std::string myClass()const  { return className.str();}

  typedef moments self_type;
  static auto get_constructor_fields() {

    return std::tuple_cat(
        base_type::get_constructor_fields(),
        std::tuple(grammar::field(C<self_type>{}, "Partial_DlogL",
                                  &self_type::partial_DlogL),
                   grammar::field(C<self_type>{}, "Partial_Auxiliar",
                                  &self_type::partial_Auxiliar)));
  }

  std::vector<moments<evidence::DlogLikelihood>> const &partial_DlogL() const {
    return partial_;
  }
  std::vector<moments<Aux>> const &partial_Auxiliar() const { return aux_; }

  moments(const std::vector<evidence::PartialDLogLikelihood<Aux>> &data)
      : base_type(data), partial_(calc_moments(data)),
                                            aux_(calc_aux_moments(data)) {}

  moments(const base_type &base, const moments<double> &GHG,
          const std::vector<moments<evidence::DlogLikelihood>> &dlog,
          const std::vector<moments<Aux>> &newAux)
      : base_type(base), partial_{dlog}, aux_{newAux} {}
  moments() = default;
  auto data_row(MOMENTS m, std::size_t isample, std::size_t i,
                std::size_t j) const {
    return std::tuple_cat(partial_Auxiliar()[isample].data_row(m),
                          base_type::data_row(m, i, j),
                          partial_DlogL()[isample].data_row(m, i, j));
  }

private:
  std::vector<moments<evidence::DlogLikelihood>> partial_;
  std::vector<moments<Aux>> aux_;

  std::vector<moments<evidence::DlogLikelihood>>
  calc_moments(const std::vector<evidence::PartialDLogLikelihood<Aux>> &data) {
    std::size_t n = 0;
    if (data.size() > 0) {
      n = data[0].partial_DlogL().size();
    }
    std::vector<moments<evidence::DlogLikelihood>> out(n);
    for (std::size_t i = 0; i < n; ++i) {
      out[i] = moments<evidence::DlogLikelihood>(
          data,
          [](const evidence::PartialDLogLikelihood<Aux> &p, std::size_t ii) {
            return p.partial_DlogL()[ii];
          },
          i);
    }
    return out;
  }  std::vector<moments<Aux>> calc_aux_moments(
      const std::vector<evidence::PartialDLogLikelihood<Aux>> &data) {
    auto n = data[0].partial_DlogL().size();
    std::vector<moments<Aux>> out(n);
    for (std::size_t i = 0; i < n; ++i) {
      out[i] = moments<Aux>(data,
                            [](const evidence::PartialDLogLikelihood<Aux> &p,
                               std::size_t ii) { return p.getAuxiliar()[ii]; },
                            i);
    }
    return out;
  }
};

namespace evidence {

class Distribution_Function_Model {

  template <class Data, class Parameters>
  auto compute_Distribution(const Data &d, const Parameters &p) {}
  template <class Data, class Parameters>
  auto compute_Distribution_aux(const Data &data, const Parameters &p,
                                std::ostream &os) {}

  template <class Data, class Parameters, class Aux>
  auto compute_Distribution_aux(const Data &data, const Parameters &p,
                                std::ostream &os, const Aux &aux) {}
};

/// Aproximacion por Fisher Information Matrix al Hessiano
///
///

template <class Distribution_Model, class Data> class Likelihood_Model {
public:
  typedef M_Matrix<double> Parameters;

  typedef Likelihood_Model self_type;
  typedef Cs<Distribution_Model, Data> template_types;
  constexpr static auto const className = my_static_string("Likelihood_Model") +
                                          my_trait<template_types>::className;

  myOptional_t<logLikelihood> compute_Likelihood(const Parameters &p) const {
    typedef myOptional_t<logLikelihood> Op;
    auto D0 = l_.compute_Distribution(d_, p);
    if (D0.has_value()) {
      auto logL = calculate_Likelihood(D0, d_);
      std::stringstream ss;
      if (are_finite<true, double>().test(logL, ss))
        return Op(logLikelihood(logL));
      else
        return Op(false, ss.str());
    } else
      return Op(false, "fails to compute model data" + D0.error());
  }
  Likelihood_Model(const Distribution_Model &l, const Data &d) : l_{l}, d_{d} {}

  const Distribution_Model &model() const { return l_; }
  const Data &data() const { return d_; }
  void set_Data(Data &&d) { d_ = std::move(d); }

protected:
  Distribution_Model l_;
  Data d_;
};

template <class Distribution_Model, class Data>
class FIM_Model : public Likelihood_Model<Distribution_Model, Data> {
public:
  typedef FIM_Model self_type;
  typedef Likelihood_Model<Distribution_Model, Data> base_type;
  constexpr static auto const className =
      my_static_string("FIM_Model") + my_trait<base_type>::className;
  typedef M_Matrix<double> Parameters;
  typedef typename Distribution_Model::aux_type Aux;

  typedef Likelihood_Model<Distribution_Model, Data> L;

  auto compute_DLikelihood_init(const Parameters &p, std::ostream &os) const {
    return compute_DLikelihood_init(p, os, L::data());
  }

  template <class mySample>
  auto compute_DLikelihood(const Parameters &p, const mySample &current,
                           std::ostream &os) const {
    return compute_DLikelihood(p, current, os, L::data());
  }

  auto compute_Distributions(const Parameters &p, std::ostream &os,
                             const Data &data,
                             const std::vector<double> &eps) const {
    std::string error;
    auto D0res = L::model().compute_Distribution_aux(data, p, os);
    typedef decltype(D0res.value().first) myDistr;
    typedef decltype(D0res.value().second) myAux;

    typedef myOptional_t<std::tuple<myDistr, std::vector<myDistr>, myAux>> Op;

    if (!D0res.has_value())
      return Op(false, "fails to compute the model :" + D0res.error());
    else {

      auto [D0, aux] = std::move(D0res).value();

      std::vector<std::decay_t<decltype(D0)>> D(p.size());
      for (std::size_t i = 0; i < p.size(); ++i) {
        Parameters x(p);
        x[i] = x[i] + eps[i];
        auto OpD = L::model().compute_Distribution_aux(data, x, os, aux);
        if (OpD.has_value())
          D[i] = std::move(OpD).value();
        else
          return Op(false, " getDLikelihood error at i=" + ToString(i) +
                               "th parameter  :" + OpD.error());
      }
      return Op(std::tuple(D0, D, aux));
    }
  }

  auto compute_Distributions_center(const Parameters &p, std::ostream &os,
                                    const Data &data,
                                    const std::vector<double> &eps) const {
    std::string error;
    auto D0res = L::model().compute_Distribution_aux(data, p, os);
    typedef decltype(D0res.value().first) myDistr;
    typedef decltype(D0res.value().second) myAux;

    typedef myOptional_t<
        std::tuple<myDistr, std::vector<myDistr>, std::vector<myDistr>, myAux>>
        Op;

    if (!D0res.has_value())
      return Op(false, "fails to compute the model :" + D0res.error());
    else {

      auto [D0, aux] = std::move(D0res).value();

      std::vector<std::decay_t<decltype(D0)>> Dp(p.size());
      std::vector<std::decay_t<decltype(D0)>> Dn(p.size());
      for (std::size_t i = 0; i < p.size(); ++i) {
        Parameters xp(p);
        xp[i] = xp[i] + 0.5*eps[i];
        auto OpDp = L::model().compute_Distribution_aux(data, xp, os, aux);
        if (OpDp.has_value())
          Dp[i] = std::move(OpDp).value();
        else
          return Op(false, " getDLikelihood error at i=" + ToString(i) +
                               "th parameter  :" + OpDp.error());
        Parameters xn(p);
        xn[i] = xn[i] - 0.5*eps[i];
        auto OpDn = L::model().compute_Distribution_aux(data, xn, os, aux);
        if (OpDn.has_value())
          Dn[i] = std::move(OpDn).value();
        else
          return Op(false, " getDLikelihood error at i=" + ToString(i) +
                               "th parameter  :" + OpDn.error());
      }
      return Op(std::tuple(D0, Dn, Dp, aux));
    }
  }

  template <class mySample>
  auto compute_DLikelihood(const Parameters &p, const mySample &current,
                           std::ostream &os, const Data &data) const

  {
    std::vector<double> eps(p.size());
    for (std::size_t i = 0; i < p.size(); ++i)
      eps[i] = 2.0 * std::sqrt(eps_f() * data.size() *
                               std::numeric_limits<double>::epsilon() *
                               std::abs(current.logL() / current.H()(i, i)));

    return compute_DLikelihood(p, os, data, eps);
  }

  auto compute_DLikelihood(const Parameters &p, std::ostream &os,
                           const Data &data,
                           const std::vector<double> &epsG) const

  {
    typedef myOptional_t<DlogLikelihood> Op;

    auto res = compute_Distributions(p, os, data, epsG);
    if (!res)
      return Op(false, res.error());
    else {

      auto [D0, D, aux] = std::move(res).value();
      return compute_DLikelihood(D0, D, os, data, epsG);
    }
  }

  auto compute_DLikelihood_center(const Parameters &p, std::ostream &os,
                                  const Data &data,
                                  const std::vector<double> &epsG) const

  {
    typedef myOptional_t<DlogLikelihood> Op;

    auto res = compute_Distributions_center(p, os, data, epsG);
    if (!res)
      return Op(false, res.error());
    else {

      auto [D0, Dp, Dn, aux] = std::move(res).value();
      return compute_DLikelihood_center(D0, Dp, Dn, os, data, epsG);
    }
  }

  auto compute_DLikelihood_init(const Parameters &p, std::ostream &os,
                                const Data &data) const

  {
    std::vector<double> epsG(p.size(), eps_);
    return compute_DLikelihood(p, os, data, epsG);
  }

  template <class Distribution>
  auto compute_DLikelihood(const Distribution &D0,
                           const std::vector<Distribution> &D, std::ostream &os,
                           const Data &data,
                           const std::vector<double> &eps) const {
    typedef myOptional_t<DlogLikelihood> Op;

    auto logL = calculate_Likelihood(D0, data);
    std::stringstream ss;
    if (!are_finite<true, double>().test(std::get<0>(logL), ss))
      return Op(false, " getDLikelihood error for getLikelihood " + ss.str());
    auto G = calculate_Gradient(D0, D, data, eps);
    if (!are_finite<true, M_Matrix<double>>().test(G, ss))
      return Op(false, " getDLikelihood error for gradient " + ss.str());

    auto H = calculate_Hessian(D0, D, data, eps, os);
    if (!are_finite<true, M_Matrix<double>>().test(H, ss))
      return Op(false, " getDLikelihood error for Hessian " + ss.str());
    return Op(DlogLikelihood(logL, std::move(G), std::move(H)));
  }

  template <class Distribution>
  auto compute_DVLikelihood(const Distribution &D0,
                            const std::vector<Distribution> &D,
                            std::ostream &os, const Data &data,
                            const std::vector<double> &eps) const {
    typedef myOptional_t<DVlogLikelihood> Op;

    auto logL = calculate_Likelihood(D0, data);
    std::stringstream ss;
    if (!are_finite<true, double>().test(std::get<0>(logL), ss))
      return Op(false, " getDLikelihood error for getLikelihood " + ss.str());
    auto G = calculate_Gradient(D0, D, data, eps);
    if (!are_finite<true, M_Matrix<double>>().test(G, ss))
      return Op(false, " getDLikelihood error for gradient " + ss.str());
    auto vG = calculate_vGradient(D0, D, data, eps);
    if (!are_finite<true, M_Matrix<double>>().test(vG, ss))
      return Op(false,
                " getDLikelihood error for gradient variance" + ss.str());

    auto H = calculate_Hessian(D0, D, data, eps, os);
    if (!are_finite<true, M_Matrix<double>>().test(H, ss))
      return Op(false, " getDLikelihood error for Hessian " + ss.str());
    auto DV =
        DVlogLikelihood(true, logL, std::move(G), std::move(H), std::move(vG));
    DV.calc();
    return Op(std::move(DV));
  }

  template <class Distribution>
  auto compute_DLikelihood_center(const Distribution &D0,
                                  const std::vector<Distribution> &Dn,
                                  const std::vector<Distribution> &Dp,
                                  std::ostream &os, const Data &data,
                                  const std::vector<double> &epsG) const {
    typedef myOptional_t<DlogLikelihood> Op;

    auto logL = calculate_Likelihood(D0, data);
    std::stringstream ss;
    if (!are_finite<true, double>().test(std::get<0>(logL), ss))
      return Op(false, " getDLikelihood error for getLikelihood " + ss.str());
    auto G = calculate_Gradient_center(Dn, Dp, data, epsG);
    if (!are_finite<true, M_Matrix<double>>().test(G, ss))
      return Op(false, " getDLikelihood error for gradient " + ss.str());

    auto H = calculate_Hessian_center(D0, Dn, Dp, data, epsG, os);
    if (!are_finite<true, M_Matrix<double>>().test(H, ss))
      return Op(false, " getDLikelihood error for Hessian " + ss.str());
    return Op(DlogLikelihood(logL, std::move(G), std::move(H)));
  }

  template <class Distribution>
  auto compute_DVLikelihood_center(const Distribution &D0,
                                   const std::vector<Distribution> &Dn,
                                   const std::vector<Distribution> &Dp,
                                   std::ostream &os, const Data &data,
                                   const std::vector<double> &epsG) const {
    typedef myOptional_t<DVlogLikelihood> Op;

    auto logL = calculate_Likelihood(D0, data);
    std::stringstream ss;
    if (!are_finite<true, double>().test(std::get<0>(logL), ss))
      return Op(false, " getDLikelihood error for getLikelihood " + ss.str());
    auto G = calculate_Gradient_center(Dn, Dp, data, epsG);
    if (!are_finite<true, M_Matrix<double>>().test(G, ss))
      return Op(false, " getDLikelihood error for gradient " + ss.str());
    auto vG = calculate_vGradient_center(Dn, Dp, data, epsG);

    auto H = calculate_Hessian_center(D0, Dn, Dp, data, epsG, os);
    if (!are_finite<true, M_Matrix<double>>().test(H, ss))
      return Op(false, " getDLikelihood error for Hessian " + ss.str());
    return Op(
        DVlogLikelihood(true, logL, std::move(G), std::move(H), std::move(vG)));
  }

  auto compute_PartialDLikelihood_init(const Parameters &p, std::ostream &os,
                                       const Data &data) const {
    std::vector<double> eps(p.size(), eps_);
    return compute_PartialDLikelihood(p, os, data, eps);
  }

  auto compute_PartialDLikelihood(const Parameters &p, std::ostream &os,
                                  const Data &data,
                                  const std::vector<double> &epsG) const {
    auto res = compute_Distributions(p, os, data, epsG);
    typedef std::decay_t<decltype(res.value())> resType;
    typedef typename std::tuple_element_t<2, resType>::value_type Aux;
    typedef myOptional_t<PartialDLogLikelihood<Aux>> Op;
    if (!res)
      return Op(false, res.error());
    else {
      auto [D0, D, aux] = std::move(res).value();

      auto Dlik = compute_DVLikelihood(D0, D, os, data, epsG);
      if (!Dlik)
        return Op(false, Dlik.error());

      std::vector<DlogLikelihood> partials(D0.size());

      for (std::size_t n = 0; n < D0.size(); ++n) {

        auto logL = calculate_Likelihood(D0[n], data[n]);
        std::stringstream ss;
        if (!are_finite<true, double>().test(std::get<0>(logL), ss))
          return Op(false,
                    " getDLikelihood error for getLikelihood " + ss.str());
        M_Matrix<double> G(1, p.size());
        for (std::size_t i = 0; i < p.size(); ++i) {
          G[i] = calculate_Gradient(D0[n], D[i][n], data[n], epsG[i]);
        }
        if (!are_finite<true, M_Matrix<double>>().test(G, ss))
          return Op(false, " getDLikelihood error for gradient " + ss.str());

        auto H = calculate_Hessian(n, D0, D, epsG);
        if (!are_finite<true, M_Matrix<double>>().test(H, ss))
          return Op(false, " getDLikelihood error for Hessian " + ss.str());
        partials[n] = DlogLikelihood(logL, std::move(G), std::move(H));
      }

      return Op(PartialDLogLikelihood<typename decltype(aux)::value_type>(
          std::move(Dlik).value(), std::move(partials), std::move(aux)));
    }
  }

  auto
  compute_PartialDLikelihood_center(const Parameters &p, std::ostream &os,
                                    const Data &data,
                                    const std::vector<double> &epsG) const {
    auto res = compute_Distributions_center(p, os, data, epsG);
    typedef std::decay_t<decltype(res.value())> resType;
    typedef typename std::tuple_element_t<3, resType>::value_type Aux;
    typedef myOptional_t<PartialDLogLikelihood<Aux>> Op;
    if (!res)
      return Op(false, res.error());
    else {
      auto [D0, Dn, Dp, aux] = std::move(res).value();

      auto Dlik = compute_DVLikelihood_center(D0, Dn, Dp, os, data, epsG);
      if (!Dlik)
        return Op(false, Dlik.error());

      std::vector<DlogLikelihood> partials(D0.size());

      for (std::size_t n = 0; n < D0.size(); ++n) {

        auto logL = calculate_Likelihood(D0[n], data[n]);
        std::stringstream ss;
        if (!are_finite<true, double>().test(std::get<0>(logL), ss))
          return Op(false,
                    " getDLikelihood error for getLikelihood " + ss.str());
        M_Matrix<double> G(1, p.size());
        for (std::size_t i = 0; i < p.size(); ++i) {
          G[i] = calculate_Gradient(Dn[i][n], Dp[i][n], data[n], epsG[i]);
        }
        if (!are_finite<true, M_Matrix<double>>().test(G, ss))
          return Op(false, " getDLikelihood error for gradient " + ss.str());

        auto H = calculate_Hessian_center(n, D0, Dn, Dp, epsG);
        if (!are_finite<true, M_Matrix<double>>().test(H, ss))
          return Op(false, " getDLikelihood error for Hessian " + ss.str());
        partials[n] = DlogLikelihood(logL, std::move(G), std::move(H));
      }

      return Op(PartialDLogLikelihood<Aux>(
          std::move(Dlik).value(), std::move(partials), std::move(aux)));
    }
  }

  FIM_Model(const Distribution_Model &l, const Data &d, double eps, double epsf)
      : Likelihood_Model<Distribution_Model, Data>(l, d), eps_{eps},
        eps_f_{epsf} {}

  void set_Data(Data &&d) { base_type::set_Data(std::move(d)); }

  double eps_f() const { return eps_f_; }

private:
  double eps_;
  double eps_f_;
  };

  template <class Parameters_Distribution, class Parameters_Values,
            class ExperimentData, class logLikelihood>
  class Likelihood_Analisis {

  public:
    typedef Likelihood_Analisis self_type;
    typedef Cs<Parameters_Distribution, Parameters_Values, ExperimentData,
               logLikelihood>
        template_types;
    constexpr static auto const className =
        my_static_string("Likelihood_Analisis") +
        my_trait<template_types>::className;

    static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "prior",
                       &self_type::get_Parameters_Distribution),
        grammar::field(C<self_type>{}, "simulations",
                       &self_type::get_Parameters),
        grammar::field(C<self_type>{}, "simulations",
                         &self_type::get_Simulations),
          grammar::field(C<self_type>{}, "likelihoods",
                         &self_type::get_Likelihoods));
    }

    auto size() const { return s_.size(); }

    double mean_logL() const {
      double out = s_[0].logL();
      for (std::size_t i = 1; i < s_.size(); ++i)
        out += s_[i].logL();
      return out / s_.size();
    }

    M_Matrix<double> mean_Gradient() const {
      M_Matrix<double> out = s_[0].G();
      for (std::size_t i = 1; i < s_.size(); ++i)
        out += s_[i].G();
      return out / s_.size();
    }
    std::vector<M_Matrix<double>> partial_mean_Gradient() const {
      auto ns = s_[0].partial_DlogL().size();
      std::vector<M_Matrix<double>> out(ns);
      for (std::size_t n = 0; n < ns; ++n) {
        M_Matrix<double> m = s_[0].partial_DlogL()[n].G();
        for (std::size_t i = 1; i < s_.size(); ++i)
          m += s_[i].partial_DlogL()[n].G();

        out[n] = m / s_.size();
      }
      return out;
    }
    M_Matrix<double> mean_sqr_Gradient() const {
      M_Matrix<double> out = quadraticForm_XTX(s_[0].G());
      for (std::size_t i = 1; i < s_.size(); ++i)
        out += quadraticForm_XTX(s_[i].G());
      return out / s_.size();
    }
    std::vector<M_Matrix<double>> partial_sqr_Gradient() const {
      auto ns = s_[0].partial_DlogL().size();
      std::vector<M_Matrix<double>> out(ns);
      for (std::size_t n = 0; n < ns; ++n) {
        M_Matrix<double> m = quadraticForm_XTX(s_[0].partial_DlogL()[n].G());
        for (std::size_t i = 1; i < s_.size(); ++i)
          m += quadraticForm_XTX(s_[i].partial_DlogL()[n].G());
        out[n] = m / s_.size();
      }
      return out;
    }

    M_Matrix<double> mean_Hessian() const {
      M_Matrix<double> out = s_[0].H();
      for (std::size_t i = 1; i < s_.size(); ++i)
        out += s_[i].H();
      return out / s_.size();
    }

    std::vector<M_Matrix<double>> partial_mean_Hessian() const {
      auto ns = s_[0].partial_DlogL().size();
      std::vector<M_Matrix<double>> out(ns);
      for (std::size_t n = 0; n < ns; ++n) {
        M_Matrix<double> m = s_[0].partial_DlogL()[n].H();
        for (std::size_t i = 1; i < s_.size(); ++i)
          m += s_[i].partial_DlogL()[n].H();

        out[n] = m / s_.size();
      }
      return out;
    }

    auto &get_Parameters_Distribution() const { return prior_; }
    auto &get_Parameters() const { return p_; }

    std::vector<ExperimentData> const &get_Simulations() const { return e_; }

    std::vector<logLikelihood> const &get_Likelihoods() const { return s_; }

    moments<logLikelihood> get_Likelihood_Moments() const {
      return moments<logLikelihood>(get_Likelihoods());
    }

    moments<ExperimentData> get_Experiment_Moments() const {
      return moments<ExperimentData>(get_Simulations());
    }

    moments<Parameters_Values> get_Parameters_Moments() const {
      return moments<Parameters_Values>(get_Parameters());
    }

  Likelihood_Analisis(Parameters_Distribution &&prior,
                        std::vector<Parameters_Values> &&par,
                        std::vector<ExperimentData> &&e,
                        std::vector<logLikelihood> &&s)
        : prior_{std::move(prior)}, p_{std::move(par)}, e_{std::move(e)},
          s_{std::move(s)} {}

  Likelihood_Analisis(const Parameters_Distribution &prior,
                      const std::vector<Parameters_Values> &par,
                      const std::vector<ExperimentData> &e,
                      const std::vector<logLikelihood> &s)
      : prior_{prior}, p_{par}, e_{e}, s_{s} {}

  Likelihood_Analisis() = default;

  io::myDataFrame<double, std::size_t, std::string> make_Data_Frame() const {
    io::myDataFrame<double, std::size_t, std::string> d;
    d.insert_column("statistic", C<std::string>());
    d.insert_column("n_sample", C<std::size_t>());

    ExperimentData::insert_col(d);

    d.insert_column("ParName", C<std::string>());
    d.insert_column("Par_trValue", C<double>());
    d.insert_column("Par_Value", C<double>());
    d.insert_column("ParName_T", C<std::string>());
    d.insert_column("Par_trValue_T", C<double>());
    d.insert_column("Par_Value_T", C<double>());

    logLikelihood::insert_col(d);
    auto nsamples = s_.size();
    for (std::size_t i = 0; i < s_.size(); ++i) {
      auto data_sample = std::tuple(std::string("value"), i);
      auto &e = e_[e_.size() == nsamples ? i : 0];
      auto &p = p_[p_.size() == nsamples ? i : 0];
      auto &l = s_[i];
      assert(e.steps().size() == l.partial_DlogL().size());
      auto nsteps = e.steps().size();
      for (std::size_t i_step = 0; i_step < nsteps; ++i_step) {
        auto data_exp = e.data_row(i_step);
        auto npar = p.size();
        assert(npar == l.G().size());
        for (std::size_t i_par = 0; i_par < npar; ++i_par) {
          if constexpr (false)
            for (std::size_t i_parT = 0; i_parT <= i_par; ++i_parT) {
              auto data_lik = l.data_row(i_step, i_par, i_parT);
              auto data_par = std::tuple(
                  prior_.name(i_par), p[i_par],
                  prior_.tr_to_Parameter(p[i_par], i_par), prior_.name(i_parT),
                  p[i_parT], prior_.tr_to_Parameter(p[i_parT], i_parT));
              auto data_t =
                  std::tuple_cat(data_sample, data_exp, data_par, data_lik);
              d.push_back_t(data_t);
            }
          else {
            auto data_lik = l.data_row(i_step, i_par, i_par);
            auto data_par = std::tuple(prior_.name(i_par), p[i_par],
                                       prior_.tr_to_Parameter(p[i_par], i_par),
                                       prior_.name(i_par), p[i_par],
                                       prior_.tr_to_Parameter(p[i_par], i_par));
            auto data_t =
                std::tuple_cat(data_sample, data_exp, data_par, data_lik);
            d.push_back_t(data_t);
          }
        }
      }
    }

    auto likmom = get_Likelihood_Moments();
    auto expmom = get_Experiment_Moments();
    auto parmom = get_Parameters_Moments();

    for (auto m : MOMENTS_VALUES) {
      auto data_sample = std::tuple(MOMENTS_string[m], 0ul);
      assert(expmom.steps().size() == likmom.partial_DlogL().size());
      auto nsteps = expmom.steps().size();
      for (std::size_t i_step = 0; i_step < nsteps; ++i_step) {
        auto data_exp = expmom.data_row(m, i_step);
        auto npar = prior_.size();
        assert(npar == likmom.G().mean().size());
        for (std::size_t i_par = 0; i_par < npar; ++i_par) {
          for (std::size_t i_parT = 0; i_parT <= i_par; ++i_parT) {
            auto data_lik = likmom.data_row(m, i_step, i_par, i_parT);
            auto data_par = std::tuple_cat(
                std::tuple(prior_.name(i_par)),
                parmom.data_row(m, i_par, i_parT),
                std::tuple(
                    prior_.tr_to_Parameter(
                        std::get<0>(parmom.data_row(m, i_par, i_par)), i_par),
                    prior_.name(i_parT)),
                parmom.data_row(m, i_parT, i_par),
                std::tuple(prior_.tr_to_Parameter(
                    std::get<0>(parmom.data_row(m, i_parT, i_parT)), i_parT)));
            auto data_t =
                std::tuple_cat(data_sample, data_exp, data_par, data_lik);
            d.push_back_t(data_t);
          }
        }
      }
    }

    return d;
  }

private:
  Parameters_Distribution prior_;
  std::vector<Parameters_Values> p_;
  std::vector<ExperimentData> e_;
  std::vector<logLikelihood> s_;
  };

  template <class ExperimentSimulation, class PartialDLogLikelihood>
  class Likelihood_sample {
  public:
    constexpr static auto const className =
        my_static_string("Likelihood_Test_sample");

    typedef Likelihood_sample self_type;
    static auto get_constructor_fields() {
    return std::make_tuple(grammar::field(C<self_type>{}, "simulations",
                                          &self_type::getSimulations),
                             grammar::field(C<self_type>{}, "likelihoods",
                                            &self_type::getLikelihoods));
    }

    auto size() const { return s_.size(); }
    M_Matrix<double> mean_Gradient() const {
      M_Matrix<double> out = s_[0].G();
      for (std::size_t i = 1; i < s_.size(); ++i)
        out += s_[i].G();
      return out / s_.size();
    }
    std::vector<M_Matrix<double>> partial_mean_Gradient() const {
      auto ns = s_[0].partial_DlogL().size();
      std::vector<M_Matrix<double>> out(ns);
      for (std::size_t n = 0; n < ns; ++n) {
        M_Matrix<double> m = s_[0].partial_DlogL()[n].G();
        for (std::size_t i = 1; i < s_.size(); ++i)
          m += s_[i].partial_DlogL()[n].G();

        out[n] = m / s_.size();
      }
      return out;
    }
    M_Matrix<double> mean_sqr_Gradient() const {
      M_Matrix<double> out = quadraticForm_XTX(s_[0].G());
      for (std::size_t i = 1; i < s_.size(); ++i)
        out += quadraticForm_XTX(s_[i].G());
      return out / s_.size();
    }
    std::vector<M_Matrix<double>> partial_sqr_Gradient() const {
      auto ns = s_[0].partial_DlogL().size();
      std::vector<M_Matrix<double>> out(ns);
      for (std::size_t n = 0; n < ns; ++n) {
        M_Matrix<double> m = quadraticForm_XTX(s_[0].partial_DlogL()[n].G());
        for (std::size_t i = 1; i < s_.size(); ++i)
          m += quadraticForm_XTX(s_[i].partial_DlogL()[n].G());
        out[n] = m / s_.size();
      }
      return out;
    }

    M_Matrix<double> mean_Hessian() const {
      M_Matrix<double> out = s_[0].H();
      for (std::size_t i = 1; i < s_.size(); ++i)
        out += s_[i].H();
      return out / s_.size();
    }

    std::vector<M_Matrix<double>> partial_mean_Hessian() const {
      auto ns = s_[0].partial_DlogL().size();
      std::vector<M_Matrix<double>> out(ns);
      for (std::size_t n = 0; n < ns; ++n) {
        M_Matrix<double> m = s_[0].partial_DlogL()[n].H();
        for (std::size_t i = 1; i < s_.size(); ++i)
          m += s_[i].partial_DlogL()[n].H();

        out[n] = m / s_.size();
      }
      return out;
    }

    std::vector<ExperimentSimulation> const &getSimulations() const { return e_; }

    std::vector<PartialDLogLikelihood> const &getLikelihoods() const {
      return s_;
    }

    Likelihood_sample(std::vector<ExperimentSimulation> &&e,
                      std::vector<PartialDLogLikelihood> &&s)
        : e_{std::move(e)}, s_{std::move(s)} {}

    Likelihood_sample(const std::vector<ExperimentSimulation> &e,
                      const std::vector<PartialDLogLikelihood> &s)
        : e_{e}, s_{s} {}

    Likelihood_sample() = default;

  private:
    std::vector<ExperimentSimulation> e_;
    std::vector<PartialDLogLikelihood> s_;
  };

  class Likelihood_Test {
  public:
    static myOptional_t<std::pair<double, std::size_t>>
    chitest(const M_Matrix<double> &mean, const M_Matrix<double> &Cov,
            std::size_t nsamples, std::ostream &os) {
      typedef myOptional_t<std::pair<double, std::size_t>> Op;
      auto Cov_inv = inv(Cov);
      if (!Cov_inv) {
        return Op(false, "Covariance matrix is not inversible");
      } else {
        double chi2 = xTSigmaX(mean, Cov_inv.value()) * nsamples;
        os << "\nchi2 =" << chi2 << " df=" << mean.size();

        return Op({chi2, nsamples * mean.size()});
      }
    }

    static std::pair<std::vector<double>, std::size_t>
    t_test(const M_Matrix<double> &mean, const M_Matrix<double> &Cov,
           std::size_t nsamples, std::ostream &os) {

      std::vector<double> out(mean.size());
      for (std::size_t i = 0; i < mean.size(); ++i)
        if (mean[i] == 0)
          out[i] = 0;
        else
          out[i] = mean[i] / std::sqrt(Cov(i, i) / nsamples);
      os << "\nt test =" << out << " df=" << nsamples;
      return {out, nsamples};
    }

    static myOptional_t<std::pair<double, std::size_t>>
    Cov_test_1(const M_Matrix<double> &Cov, const M_Matrix<double> &CovExp,
               std::size_t n, std::ostream &os) {
      typedef myOptional_t<std::pair<double, std::size_t>> Op;
      auto p = CovExp.ncols();
      auto Cov_inv = inv(CovExp);
      if (!Cov_inv.has_value())
        return Op(false, "singular target covariance");
      auto S = Cov * Cov_inv.value();
      auto chol_S = chol(S, "lower");
      double logdet = logDiagProduct(chol_S.value());

      double TLR = n * (p * std::log(1) - logdet + Trace(S) - p);
      os << "\n Trace(S)=" << Trace(S) << " logdet=" << logdet;
      os << " TLR=" << TLR;
      TLR = (1.0 - 1.0 / (6 * n - 1) * (2 * p + 1 - 2.0 / (p + 1))) * TLR;
      os << " TLRc=" << TLR << " df=" << (p * (p + 1)) / 2;
      return {{TLR, (p * (p + 1)) / 2}};
    }

    static myOptional_t<std::pair<double, std::size_t>>
    Cov_test_2(const M_Matrix<double> &Cov, const M_Matrix<double> &CovExp,
               std::size_t, std::ostream &os) {
      typedef myOptional_t<std::pair<double, std::size_t>> Op;
      auto p = CovExp.ncols();
      auto Cov_inv = inv(CovExp);
      if (!Cov_inv.has_value())
        return Op(false, "singular target covariance");
      auto S = Cov * Cov_inv.value();
      auto TJn = 1.0 / p *
                 Trace(sqr(S / (1.0 / p * Trace(S)) -
                           Matrix_Generators::eye<double>(p)));
      os << "\nTJn= " << TJn - 1 << " df=" << (p * (p + 1)) / 2;

      return {{TJn, (p * (p + 1)) / 2 - 1}};
    }

    static auto sample_Cov(std::mt19937_64 &mt, const M_Matrix<double> &Cov,
                           std::size_t nsamples) {
      typedef myOptional_t<std::pair<M_Matrix<double>, M_Matrix<double>>> Op;
      auto k = Cov.ncols();
      M_Matrix<double> mean(1, k, 0.0);
      auto normal = Normal_Distribution<M_Matrix<double>>::make(mean, Cov);
      if (!normal.has_value()) {
        return Op(false, "\nsingular Covariance!!\n");
      } else {
        M_Matrix<double> S(k, k, Matrix_TYPE::SYMMETRIC, 0.0);
        M_Matrix<double> m(1, k, 0.0);
        for (std::size_t i = 0; i < nsamples; ++i) {
          auto x = normal.value().sample(mt);
          S += quadraticForm_XTX(x);
          m += x;
        }
        return Op(std::pair(m / nsamples, S / nsamples));
      }
    }

  template <class Simulation_Model, class FIM_Model, class Data,
            class ParametersDistribution, class Parameters>
  static auto compute_Sample(std::ostream &os, const Simulation_Model &sim,
                             const FIM_Model &fim, const Data &e,
                             const ParametersDistribution &prior,
                               const Parameters &p, std::vector<double> &epsG,
                               bool centered, std::mt19937_64 &mt) {
    auto data = sim.compute_simulation(e, p, mt);
    typedef std::decay_t<decltype(data.value())> sim_type;
    typedef myOptional_t<
        std::pair<sim_type, PartialDLogLikelihood<typename FIM_Model::Aux>>>
        Op;
    if (!data.has_value())
      return Op(false, "Simulation failed " + data.error());
    else {
      auto Par = prior.Parameter_to_tr(p);
      if (!centered) {
        auto Dlik = fim.compute_PartialDLikelihood(Par, os, data.value(), epsG);
        if (!Dlik)
          return Op(false, "Likelihood failed!!: " + Dlik.error());

        return Op(std::pair(std::move(data).value(), std::move(Dlik).value()));
      } else {
        auto Dlik =
            fim.compute_PartialDLikelihood_center(Par, os, data.value(), epsG);
        if (!Dlik)
          return Op(false, "Likelihood failed!!: " + Dlik.error());

        return Op(std::pair(std::move(data).value(), std::move(Dlik).value()));
      }
    }
  }

  template <class Simulation_Model, class der_Model, class Data,
            class ParametersDistributionDerivative, class Parameters>
  static auto
  compute_Sample_derivative(std::ostream &os, const Simulation_Model &sim,
                            const der_Model &dm, const Data &e,
                            const ParametersDistributionDerivative &dprior,
                            const Parameters &p, std::mt19937_64 &mt) {
    auto data = sim.compute_simulation(e, p, mt);
    typedef std::decay_t<decltype(data.value())> sim_type;
    typedef std::invoke_result_t<
        decltype(&ParametersDistributionDerivative::Parameter_to_tr),
        ParametersDistributionDerivative, Parameters>
        partype;
    typedef typename std::invoke_result_t<
        decltype(&der_Model::template compute_PartialDLikelihood<sim_type>),
        der_Model, partype, std::ostream &, sim_type>::value_type lik_type;
    typedef myOptional_t<std::pair<sim_type, lik_type>> Op;
    if (!data.has_value())
      return Op(false, "Simulation failed " + data.error());
    else {
      auto Par = dprior.Parameter_to_tr(p);
      auto Dlik = dm.compute_PartialDLikelihood(Par, os, data.value());
      if (!Dlik) {
        return Op(false, "Likelihood failed!!: " + Dlik.error());
      }

      return Op(std::pair(std::move(data).value(), std::move(Dlik).value()));
    }
  }

  template <class Simulation_Model, class FIM_Model, class Data,
            class ParametersDistribution, class Parameters>
  static auto compute_Sample_init(std::ostream &os, const Simulation_Model &sim,
                                  const FIM_Model &fim, const Data &e,
                                  const ParametersDistribution &prior,
                                  const Parameters &p, std::mt19937_64 &mt) {
    auto data = sim.compute_simulation(e, p, mt);
    typedef std::decay_t<decltype(data.value())> sim_type;
    typedef myOptional_t<
        std::pair<sim_type, PartialDLogLikelihood<typename FIM_Model::Aux>>>
        Op;
    if (!data.has_value())
      return Op(false, "Simulation failed " + data.error());
    else {
      auto Par = prior.Parameter_to_tr(p);
      auto Dlik = fim.compute_PartialDLikelihood_init(Par, os, data.value());
      if (!Dlik)
        return Op(false, "Likelihood failed!!: " + Dlik.error());

      return Op(std::pair(std::move(data).value(), std::move(Dlik).value()));
    }
  }
  static std::vector<std::mt19937_64> mts(std::mt19937_64 &mt, std::size_t n) {
    std::uniform_int_distribution<typename std::mt19937_64::result_type> useed;
    std::vector<std::mt19937_64> out;
    for (std::size_t i = 0; i < n; ++i)
    {
      auto s=useed(mt);
      std::cerr<<"seed "<<i<<" ="<<s<<"\n";
      out.emplace_back(s);
    }
    return out;
  }

  template <class Simulation_Model, class FIM_Model, class Data,
            class ParametersDistribution, class Parameters>
  static auto get_sample(std::ostream &os, const Simulation_Model &sim,
                         const FIM_Model &fim, const Data &e,
                         const ParametersDistribution &prior,
                         const Parameters &p, std::mt19937_64 &mt,
                         std::size_t nsamples, bool init, bool centered) {
    typedef std::decay_t<decltype(sim.compute_simulation(e, p, mt).value())>
        sim_type;

    typedef myOptional_t<Likelihood_sample<
        sim_type, PartialDLogLikelihood<typename FIM_Model::Aux>>>
        Op;
    std::stringstream ss;
    std::vector<sim_type> simuls(nsamples);
    std::vector<PartialDLogLikelihood<typename FIM_Model::Aux>> s(nsamples);

    std::vector<std::stringstream> oss(8);
    auto mtvec = mts(mt, 8);

    if (!init) {
      auto logL = compute_Sample_init(os, sim, fim, e, prior, p, mt);
      if (!logL.has_value())
        return Op(false, "erron in get_sample"+logL.error());
      else {
        std::vector<double> eps(p.size());
        for (std::size_t i = 0; i < 1; ++i) {
          double myLogL = logL.value().second.logL();
          auto myH = logL.value().second.H();
          for (std::size_t i = 0; i < p.size(); ++i)
            eps[i] = 2.0 * std::sqrt(fim.eps_f() * e.size() *
                                     std::numeric_limits<double>::epsilon() *
                                     std::abs(myLogL / myH(i, i)));
          os << "\n eps=" << eps;
          auto logL =
              compute_Sample(os, sim, fim, e, prior, p, eps, centered, mt);
          if (!logL.has_value())
            return Op(false, "erron in get_sample"+logL.error());
        }
        std::vector<bool> succed(8,true);
        std::vector<std::string> error(8,"");

        for (std::size_t nc = 0; nc < nsamples / 8; ++nc) {
#pragma omp parallel for
          for (std::size_t ni = 0; ni < 8; ++ni) {
            std::size_t i = nc * 8 + ni;
            auto logL = compute_Sample(oss[ni], sim, fim, e, prior, p, eps,
                                       centered, mtvec[ni]);
            if (!logL.has_value())
            {
              succed[ni] = false; error[ni]="in i=" +ToString(ni)+"  "+logL.error()+"\n";
            }
            else {
              auto pair = std::move(logL).value();
              s[i] = std::move(pair.second);
              simuls[i] = std::move(pair.first);
              //             os<<"\n in sample\n"<<s[i];
            }
          }
          for (auto &ess : oss) {
            os << ess.str();
            ess.str("");
          }
          if(!all(succed))
            return  Op(false, consolidate(error));
        }
      }
    } else {

      auto mtvec = mts(mt, 8);
      std::vector<bool> succed(8,true);
      std::vector<std::string> error(8,"");
      for (std::size_t nc = 0; nc < nsamples / 8; ++nc) {
#pragma omp parallel for
        for (std::size_t ni = 0; ni < 8; ++ni) {
          std::size_t i = nc * 8 + ni;
          auto logL =
              compute_Sample_init(oss[ni], sim, fim, e, prior, p, mtvec[ni]);
          if (!logL.has_value())
          {succed[ni] = false; error[ni]="in i=" +ToString(ni)+"  "+logL.error()+"\n";}
          else {
            auto pair = std::move(logL).value();
            s[i] = std::move(pair.second);
            simuls[i] = std::move(pair.first);
            //             os<<"\n in sample\n"<<s[i];
          }
        }
        for (auto &ess : oss) {
          os << ess.str();
          ess.str("");
        }
        if(!all(succed))
          return  Op(false, consolidate(error));

      }
    }
    os << "\n successfully obtained " + ToString(nsamples) +
              " samples for Gradient_Expectancy_Test \n";

    auto samples = Likelihood_sample(std::move(simuls), std::move(s));
    return Op(samples);
  }

  template <class Simulation_Model, class Der_Model, class Data,
            class ParametersDistributionDerivative, class Parameters>
  static auto get_sample_derivative(
      std::ostream &os, const Simulation_Model &sim, const Der_Model &dm,
      const Data &e, const ParametersDistributionDerivative &dprior,
      const Parameters &p, std::mt19937_64 &mt, std::size_t nsamples) {
    typedef std::decay_t<decltype(
        compute_Sample_derivative(os, sim, dm, e, dprior, p, mt))>
        elem_type;
    typedef typename elem_type::value_type::first_type sim_type;
    typedef typename elem_type::value_type::second_type logL_type;

    typedef myOptional_t<Likelihood_sample<sim_type, logL_type>> Op;
    std::stringstream ss;
    std::vector<sim_type> simuls(nsamples);
    std::vector<logL_type> s(nsamples);

    std::vector<std::stringstream> oss(8);

    bool succed = true;

    auto mtvec = mts(mt, 8);
    for (std::size_t nc = 0; nc < nsamples / 8; ++nc) {
#pragma omp parallel for
      for (std::size_t ni = 0; ni < 8; ++ni) {
        std::size_t i = nc * 8 + ni;
        auto logL = compute_Sample_derivative(oss[ni], sim, dm, e, dprior, p,
                                              mtvec[ni]);
        if (!logL.has_value()) {
          std::cerr << logL.error();
          assert(false);
          succed = false;
        } else {
          auto pair = std::move(logL).value();
          s[i] = std::move(pair.second);
          simuls[i] = std::move(pair.first);
          //             os<<"\n in sample\n"<<s[i];
        }
      }
      for (auto &ess : oss) {
        os << ess.str();
        ess.str("");
      }
    }
    if (!succed)

      return Op(false, "failed in one of the samples");
    os << "\n successfully obtained " + ToString(nsamples) +
              " samples for Gradient_Expectancy_Test \n";
    auto samples = Likelihood_sample(std::move(simuls), std::move(s));
    return Op(samples);
  }

  static Op_void
  Gradient_Expectancy_Test(const M_Matrix<double> &Gradient,
                           const M_Matrix<double> &GradientVariance,
                           std::size_t nsamples, double pvalue,
                           std::ostream &os)

  {

    auto t_stud = t_test(Gradient, GradientVariance, nsamples, os);
    auto chi2_FIM = chitest(Gradient, GradientVariance, nsamples, os);
    if (!chi2_FIM)
      return Op_void(false,
                     "Gradient_Expectancy_Test could not be performed: " +
                         chi2_FIM.error());
    else {
      auto chivalue =
          1.0 - chi2_cdf(chi2_FIM.value().first, chi2_FIM.value().second);
      if (chivalue < pvalue)
        return Op_void(false, "Gradient_Expectancy_Test indicates a low "
                              "probability of being true: chi2=" +
                                  std::to_string(chi2_FIM.value().first) +
                                  " : p(" + std::to_string(chivalue) + ")<" +
                                  ToString(pvalue) + " : " + chi2_FIM.error());
      else
        return Op_void(
            true,
            "Gradient_Expectancy_Test fall within the expected values: chi2=" +
                std::to_string(chi2_FIM.value().first) + " p(" +
                std::to_string(chivalue) + ")>" + ToString(pvalue) + " : " +
                chi2_FIM.error());
    }
  }

  static Op_void Gradient_Variance_Expectancy_Test(
      const M_Matrix<double> &GradientVariance,
      const M_Matrix<double> &ExpectedGradientVariance, std::size_t nsamples,
      double pvalue, std::ostream &os) {

    auto chi2 =
        Cov_test_1(GradientVariance, ExpectedGradientVariance, nsamples, os);
    auto chi22 =
        Cov_test_2(GradientVariance, ExpectedGradientVariance, nsamples, os);
    if (!chi2)
      return Op_void(
          false, "Gradient_Variance_Expectancy_Test could not be performed: " +
                     chi2.error());
    else {
      auto chivalue = chi2_cdf(chi2.value().first, chi2.value().second);
      if (chivalue < pvalue)
        return Op_void(false, "Gradient_Variance_Expectancy_Test indicates a "
                              "low probability of being true: p(" +
                                  std::to_string(chivalue) + ")<" +
                                  ToString(pvalue) + " : " + chi2.error());
      else
        return Op_void(true, "Gradient_Variance_Expectancy_Test fall within "
                             "the expected values: p(" +
                                 std::to_string(chivalue) + ")>" +
                                 ToString(pvalue) + " : " + chi2.error());
    }
  }

  template <class Simulation_Model, class FIM_Model, class Data,
            class ParametersDistribution, class Parameters>
  static auto compute_test(std::ostream &os, const Simulation_Model &sim,
                           const FIM_Model &fim, const Data &e,
                           const ParametersDistribution &prior,
                           const Parameters &p, std::mt19937_64 &mt,
                           std::size_t nsamples, double pvalue, bool init,

                           bool centergradient) {


    auto mysamples = get_sample(os, sim, fim, e, prior, p, mt, nsamples, init,
                                centergradient);



    if constexpr (false) {
      auto Gmean = mysamples.value().mean_Gradient();
      os << "\n Gmean \n" << Gmean;

      auto partialGmean = mysamples.value().partial_mean_Gradient();
      os << "\n partialGmean \n";
      for (std::size_t n = 0; n < partialGmean.size(); ++n)
        os << "\n" << partialGmean[n];

      auto partialGsqr = mysamples.value().partial_sqr_Gradient();
      auto partialH = mysamples.value().partial_mean_Hessian();

      auto Gsqr = mysamples.value().mean_sqr_Gradient();
      os << "\n Gsqr \n" << Gsqr;

      auto H = -mysamples.value().mean_Hessian();
      os << "\n H \n" << H;
      os << "\n elemDiv(Gsqr,H)\n" << elemDiv(Gsqr, H);

      // lets simulate results...
      auto resres = sample_Cov(mt, H, nsamples);
      M_Matrix<double> Smean = resres.value().first;
      M_Matrix<double> Ssqr = resres.value().second;

      os << "\n TEST gradient_information_test \n";
      auto gradient_information_test_test =
          Gradient_Expectancy_Test(Smean, H, nsamples, pvalue, os);

      os << gradient_information_test_test.error() << "\n\n";

      os << "\n gradient_information_test \n";
      auto gradient_information_test =
          Gradient_Expectancy_Test(Gmean, H, nsamples, pvalue, os);
      os << gradient_information_test.error() << "\n\n";
      os << "\n partial gradient_information_test \n";
      for (std::size_t n = 0; n < partialGmean.size(); ++n) {
        auto partial_gradient_information_test = Gradient_Expectancy_Test(
            partialGmean[n], partialH[n], nsamples, pvalue, os);
        os << partial_gradient_information_test.error() << "\n\n";
      }

      os << "\n TEST gradient_variance_test \n";
      auto gradient_variance_test_test =
          Gradient_Expectancy_Test(Smean, Ssqr, nsamples, pvalue, os);
      os << gradient_variance_test_test.error() << "\n\n";

      os << "\n gradient_variance_test \n";
      auto gradient_variance_test =
          Gradient_Expectancy_Test(Gmean, Gsqr, nsamples, pvalue, os);
      os << gradient_variance_test.error() << "\n\n";
      os << "\n partial gradient_variance_test \n";
      for (std::size_t n = 0; n < partialGmean.size(); ++n) {
        auto partial_gradient_variance_test = Gradient_Expectancy_Test(
            partialGmean[n], partialGsqr[n], nsamples, pvalue, os);
        os << partial_gradient_variance_test.error() << "\n\n";
      }

      os << "\n TEST gradient_variance_information_test \n";
      auto gradient_variance_information_test_test =
          Gradient_Variance_Expectancy_Test(Ssqr, H, nsamples, pvalue, os);
      os << gradient_variance_information_test_test.error() << "\n\n";

      os << "\n gradient_variance_information_test \n";
      auto gradient_variance_information_test =
          Gradient_Variance_Expectancy_Test(Gsqr, H, nsamples, pvalue, os);
      os << gradient_variance_information_test.error() << "\n\n";
      ;

      Op_void all_tests = gradient_information_test
                          << gradient_variance_information_test
                          << gradient_variance_test;
    }
    auto par = std::vector({prior.Parameter_to_tr(p)});


    auto out=Likelihood_Analisis(prior, par, mysamples.value().getSimulations(),
                                  mysamples.value().getLikelihoods());
    typedef myOptional_t<std::decay_t<decltype(out)>> Op;
    if (!mysamples)
      return Op(false, "Sampling failed!! \n"+mysamples.error());
    else return Op(out);

  }

  template <class Simulation_Model, class Der_Model, class Data,
            class ParametersDistributionDerivative, class Parameters>
  static auto
  compute_test_derivative(std::ostream &os, const Simulation_Model &sim,
                          const Der_Model &dm, const Data &e,
                          const ParametersDistributionDerivative &dprior,
                          const Parameters &p, std::mt19937_64 &mt,
                          std::size_t nsamples, double pvalue) {
    auto mysamples =
        get_sample_derivative(os, sim, dm, e, dprior, p, mt, nsamples);
    if constexpr (false) {
      auto Gmean = mysamples.value().mean_Gradient();
      os << "\n Gmean \n" << Gmean;

      auto partialGmean = mysamples.value().partial_mean_Gradient();
      os << "\n partialGmean \n";
      for (std::size_t n = 0; n < partialGmean.size(); ++n)
        os << "\n" << partialGmean[n];

      auto partialGsqr = mysamples.value().partial_sqr_Gradient();
      auto partialH = mysamples.value().partial_mean_Hessian();

      auto Gsqr = mysamples.value().mean_sqr_Gradient();
      os << "\n Gsqr \n" << Gsqr;

      auto H = -mysamples.value().mean_Hessian();
      os << "\n H \n" << H;
      os << "\n elemDiv(Gsqr,H)\n" << elemDiv(Gsqr, H);

      // lets simulate results...
      auto resres = sample_Cov(mt, H, nsamples);
      M_Matrix<double> Smean = resres.value().first;
      M_Matrix<double> Ssqr = resres.value().second;

      os << "\n TEST gradient_information_test \n";
      auto gradient_information_test_test =
          Gradient_Expectancy_Test(Smean, H, nsamples, pvalue, os);

      os << gradient_information_test_test.error() << "\n\n";

      os << "\n gradient_information_test \n";
      auto gradient_information_test =
          Gradient_Expectancy_Test(Gmean, H, nsamples, pvalue, os);
      os << gradient_information_test.error() << "\n\n";
      os << "\n partial gradient_information_test \n";
      for (std::size_t n = 0; n < partialGmean.size(); ++n) {
        auto partial_gradient_information_test = Gradient_Expectancy_Test(
            partialGmean[n], partialH[n], nsamples, pvalue, os);
        os << partial_gradient_information_test.error() << "\n\n";
      }

      os << "\n TEST gradient_variance_test \n";
      auto gradient_variance_test_test =
          Gradient_Expectancy_Test(Smean, Ssqr, nsamples, pvalue, os);
      os << gradient_variance_test_test.error() << "\n\n";

      os << "\n gradient_variance_test \n";
      auto gradient_variance_test =
          Gradient_Expectancy_Test(Gmean, Gsqr, nsamples, pvalue, os);
      os << gradient_variance_test.error() << "\n\n";
      os << "\n partial gradient_variance_test \n";
      for (std::size_t n = 0; n < partialGmean.size(); ++n) {
        auto partial_gradient_variance_test = Gradient_Expectancy_Test(
            partialGmean[n], partialGsqr[n], nsamples, pvalue, os);
        os << partial_gradient_variance_test.error() << "\n\n";
      }

      os << "\n TEST gradient_variance_information_test \n";
      auto gradient_variance_information_test_test =
          Gradient_Variance_Expectancy_Test(Ssqr, H, nsamples, pvalue, os);
      os << gradient_variance_information_test_test.error() << "\n\n";

      os << "\n gradient_variance_information_test \n";
      auto gradient_variance_information_test =
          Gradient_Variance_Expectancy_Test(Gsqr, H, nsamples, pvalue, os);
      os << gradient_variance_information_test.error() << "\n\n";
      ;

      Op_void all_tests = gradient_information_test
                          << gradient_variance_information_test
                          << gradient_variance_test;
    }
    auto par = std::vector({dprior.Parameter_to_tr(p)});
    return Likelihood_Analisis(dprior.f(), par,
                               mysamples.value().getSimulations(),
                               mysamples.value().getLikelihoods());
  }
  };

  } // namespace evidence

#endif // MYLIKELIHOOD_H
