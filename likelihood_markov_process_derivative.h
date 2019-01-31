#ifndef LIKELIHOOD_MARKOV_PROCESS_DERIVATIVE_H
#define LIKELIHOOD_MARKOV_PROCESS_DERIVATIVE_H

#include "likelihood_markov_process.h"
#include "matrixderivative.h"
#include "qmodel_derivative.h"
#include "mydistributions_derivative.h"

template <>
class Derivative<markov::mp_state_information> {
  Derivative<M_Matrix<double>> P_mean_;
  Derivative<M_Matrix<double>> P_cov_;

  Derivative<double> y_mean_;
  Derivative<double> y_var_;
  Derivative<double> plogL_;
  Derivative<double> eplogL_;
  double vplogL_;

public:
  Derivative(Derivative<M_Matrix<double>> &&P_mean__,
             Derivative<M_Matrix<double>> &&P_cov__,
             Derivative<double> &&y_mean__, Derivative<double> &&y_var__,
             Derivative<double> &&plogL__, Derivative<double> &&eplogL__, double vplogL__)
      : P_mean_{std::move(P_mean__)}, P_cov_{std::move(P_cov__)},
        y_mean_{std::move(y_mean__)}, y_var_{std::move(y_var__)},
        plogL_{std::move(plogL__)}, eplogL_{std::move(eplogL__)}, vplogL_{vplogL__} {}

  Derivative() = default;
  auto &P_mean() const { return P_mean_; }
  auto &P_cov() const { return P_cov_; }

  auto &y_mean() const { return y_mean_; }
  auto &y_var() const { return y_var_; }
  auto &plogL() const { return plogL_; }
  auto &eplogL() const { return eplogL_; }
  auto vplogL() const { return vplogL_; }

  auto &x() const { return y_mean().x(); }
  typedef Derivative self_type;
  constexpr static auto className =
      my_static_string("markov::mp_state_information_Derivative");
  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "P_mean", &self_type::P_mean),
        grammar::field(C<self_type>{}, "P_cov", &self_type::P_cov),
        grammar::field(C<self_type>{}, "y_mean", &self_type::y_mean),
        grammar::field(C<self_type>{}, "y_var", &self_type::y_var),
        grammar::field(C<self_type>{}, "plogL", &self_type::plogL),
        grammar::field(C<self_type>{}, "eplogL", &self_type::eplogL),
        grammar::field(C<self_type>{}, "vplogL", &self_type::vplogL));
  }

  static Derivative<markov::mp_state_information>
  adjust(Derivative<M_Matrix<double>> &&P_mean__,
         Derivative<M_Matrix<double>> &&P_cov__, Derivative<double> &&y_mean__,
         Derivative<double> &&y_var__, Derivative<double> &&plogL__,
         Derivative<double> &&eplogL__, double vplogL__,double min_p, double min_var) {
    return Derivative<markov::mp_state_information>(
        Derivative<Probability_distribution>::normalize(std::move(P_mean__),
                                                        min_p),
        Derivative<Probability_distribution_covariance>::normalize(
            std::move(P_cov__), min_p),
        std::move(y_mean__),
        Derivative<variance_value>::adjust(std::move(y_var__), min_var),
        std::move(plogL__), std::move(eplogL__), vplogL__);
  }

  static Op_void test(const Derivative<M_Matrix<double>> &P_mean__,
                      const Derivative<M_Matrix<double>> &P_cov__,
                      double tolerance) {
    auto ck_mean =
        Derivative<Probability_distribution>::test<true>(P_mean__, tolerance);
    auto ck_cov = Derivative<Probability_distribution_covariance>::test<true>(
        P_cov__, tolerance);
    if (ck_mean.has_value() && ck_cov.has_value())
      return Op_void(true, "");
    else {
      std::stringstream ss;
      ss << " Pmean test: " << ck_mean.error()
         << " Pcov test: " << ck_cov.error();
      return Op_void(false, ss.str());
    }
  }

  static Op_void test(const Derivative<M_Matrix<double>> &P_mean__,
                      const Derivative<M_Matrix<double>> &P_cov__,
                      const Derivative<double> &y_mean__,
                      const Derivative<double> &y_var__,
                      const Derivative<double> &plogL__,
                      const Derivative<double> &eplogL__, double minvariance,
                      double tolerance) {
    auto ck_mean =
        Derivative<Probability_distribution>::test<true>(P_mean__, tolerance);
    auto ck_cov = Derivative<Probability_distribution_covariance>::test<true>(
        P_cov__, tolerance);
    auto ck_ymean = variable_value::test<true>(y_mean__.f());
    auto ck_yvar =
        variance_value::test<true>(y_var__.f(), minvariance, tolerance);
    auto ck_plogL = logLikelihood_value::test<true>(plogL__.f());
    auto ck_eplogL = logLikelihood_value::test<true>(eplogL__.f());
    if (ck_mean.has_value() && ck_cov.has_value() && ck_ymean.has_value() &&
        ck_yvar.has_value() && ck_plogL.has_value() && ck_eplogL.has_value())
      return Op_void(true, "");
    else {
      std::stringstream ss;
      ss << " Pmean test: " << ck_mean.error()
         << " Pcov test: " << ck_cov.error()
         << " yvar test:" << ck_yvar.error();
      ss << " pLogL test: " << ck_plogL.error() << " eplogL test"
         << ck_eplogL.error();
      return Op_void(false, ss.str());
    }
  }
};

std::ostream& operator<<(std::ostream& os, const Derivative<markov::mp_state_information>& d){ return io::output_operator_on_Object(os, d);}


template <>
class Derivative<markov::hidden::MacroDMNR> : public markov::hidden::MacroDMNR {
public:
  typedef markov::hidden::MacroDMNR base_type;

  base_type const &f() const { return *this; }

  template <class Model, class Step>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior, Model &m,
      const Step &p, std::ostream &os) const {
    markov::MACROR alg;
    return run(prior, m, p, os, alg);
  }

  template <class Model, class Step, class MACROR>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior, Model &m,
      const Step &p, std::ostream &os, MACROR &alg) const {
    typedef myOptional_t<Derivative<markov::mp_state_information>> Op;
    auto Q_dto = m.get_P(p, 0);
    if (!Q_dto)
      return Op(false, "fails in auto Q_dt=m.get_P(p,0) :" + Q_dto.error());
    Derivative<Markov_Transition_step_double> Q_dt = std::move(Q_dto).value();
    return run(prior, Q_dt, m, p, os, alg);
  }

  template <class Model, class Step>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior,
      const Derivative<Markov_Transition_step_double> Q_dt, Model &m,
      const Step &p, std::ostream &os, const markov::MACROR &alg) const {
    typedef myOptional_t<Derivative<markov::mp_state_information>> Op;
    switch (alg) {
    case markov::MACRO_DMNR:
      return Derivative < markov::hidden::MacroDMNR > (tolerance()).run(prior, Q_dt, m, p, os);
    case markov::MACRO_DVNR:
      return Op(false, "Not a contemplated algorithm :" +
                           markov::MACROR_string[markov::MACRO_DMNR] +
                           " is not valid for " + className.str());
    case markov::MACRO_DMR:
      return Op(false, "Not a contemplated algorithm :" +
                           markov::MACROR_string[markov::MACRO_DVNR] +
                           " is not valid for " + className.str());
    case markov::MACRO_DVR:
    default:
      return Op(false, "Not a contemplated algorithm :" +
                           markov::MACROR_string[markov::MACRO_DVR] +
                           " is not valid for " + className.str());
    }
  }

  template <class Model, class Step>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior,
      const Derivative<Markov_Transition_step_double> &Q_dt, Model &m,
      const Step &p, std::ostream &os, markov::MACROR &alg) const {
    alg = markov::MACRO_DMNR;
    const markov::MACROR alg2 = alg;
    return run(prior, Q_dt, m, p, os, alg2);
  }

  template <class Model, class Step, class... Aux>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior,
      Derivative<Markov_Transition_step_double> Q_dt, Model &m, const Step &p,
      std::ostream &) const {
    typedef myOptional_t<Derivative<markov::mp_state_information>> Op;
    Derivative<double> e = m.noise_variance(p.nsamples());
    Derivative<double> N = m.AverageNumberOfChannels();
    M_Matrix<double> u(prior.P_mean().f().size(), 1, 1.0);

    auto SmD = prior.P_cov() - diag(prior.P_mean());

    assert((are_Equal_v(SmD,
                        Incremental_ratio(1e-6,
                                          [](auto &cov, auto &mean) {
      return cov - diag(mean);
                                          },
                                          prior.P_cov(), prior.P_mean()),
                        std::cerr)));
    Derivative<double> gSg =
        (TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.gmean_i()).getvalue() +
        (prior.P_mean() *
         (elemMult(Q_dt.gtotal_ij(), Q_dt.gmean_ij()) * Constant(u)))
            .getvalue();

    assert((are_Equal_v(
        gSg,
        Incremental_ratio(
            1e-6,
            [&u](auto const &gmean_i, auto const &SmD, auto const &P_mean,
                 auto const &gtotal_ij, auto const &gmean_ij) {
      return (TranspMult(gmean_i, SmD) * gmean_i).getvalue() +
             (P_mean * (elemMult(gtotal_ij, gmean_ij) * u)).getvalue();
            },
            Q_dt.gmean_i(), SmD, prior.P_mean(), Q_dt.gtotal_ij(),
            Q_dt.gmean_ij()),
        std::cerr)));

    Derivative<double> ms = (prior.P_mean() * Q_dt.gvar_i()).getvalue();

    auto e_mu = e + N * ms;
    auto y_mean = N * (prior.P_mean() * Q_dt.gmean_i()).getvalue();
    auto y_var = e_mu + N * gSg;
    auto y = p.y();
    if (std::isnan(y)) {
      Derivative<double> plogL(0.0,
                               prior.y_mean().x());
      Derivative<double> eplogL(0.0,
                                prior.y_mean().x());
      double vplogL=0;
      auto P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P());
      assert(are_Equal_v(P__cov,
                         Incremental_ratio(1e-4,
                                           [](auto const &P, auto const &SmD_) {
        return quadraticForm_BT_A_B(SmD_,
                                    P);
                                           },
                                           Q_dt.P(), SmD),
                         std::cerr));
      auto P_mean = prior.P_mean() * Q_dt.P();
      assert(are_Equal_v(
          P_mean,
          Incremental_ratio(
              1e-4, [](auto const &Pmean, auto const &P) { return Pmean * P; },
              prior.P_mean(), Q_dt.P()),
          std::cerr));

      P__cov += diag(P_mean);

    // std::cerr<<"\nPcov nana corr\n"<<P__cov<<"\nP_mean nana
      // corr\n"<<P_mean<<"\nQ.P \n"<<Q_dt.P();
      auto test = markov::mp_state_information::test(P_mean.f(), P__cov.f(),
                                                     tolerance());
      if (test.has_value())
        return Op(Derivative<markov::mp_state_information>::adjust(
            std::move(P_mean), std::move(P__cov), std::move(y_mean),
            std::move(y_var), std::move(plogL), std::move(eplogL), vplogL,Q_dt.min_P(),
            e.f()));
      else
        return Op(false, "fails at intertrace prediction!!: " + test.error());
    }
    auto dy = Constant(y) - y_mean;
    auto chi = dy / y_var;
    auto P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P());

    assert(are_Equal_v(P__cov,
                       Incremental_ratio(1e-4,
                                         [](auto const &P, auto const &SmD_) {
      return quadraticForm_BT_A_B(SmD_, P);
                                         },
                                         Q_dt.P(), SmD),
                       std::cerr));

    auto P_mean = prior.P_mean() * Q_dt.P();
    assert(are_Equal_v(P_mean,
                       Incremental_ratio(1e-4,
                                         [](auto const &Pmean, auto const &P) {
      return Pmean * P;
                                         },
                                         prior.P_mean(), Q_dt.P()),
                       std::cerr));
    P__cov += diag(P_mean);

    auto chi2 = dy * chi;

    Derivative<double> plogL(prior.y_mean().x());
    if (y_var.f() > 0)
      plogL = -0.5 * log(2 * PI * y_var) - 0.5 * chi2;
    else
      plogL.f() = std::numeric_limits<double>::infinity();

    Derivative<double> eplogL = -0.5 * log(2 * PI * y_var) +
                                Constant(-0.5); // e_mu+N*gSg"-N*zeta*sqr(sSg)"
    // double chilogL=(eplogL-plogL)/std::sqrt(0.5);

    double vplogL=0.5;
    auto test = markov::mp_state_information::test(
        P_mean.f(), P__cov.f(), y_mean.f(), y_var.f(), plogL.f(), eplogL.f(),
        e.f(), tolerance());
    if (!test) {
      std::stringstream ss;

      ss << "\nP_mean \n" << P_mean;
      ss << "\nPcov \n" << P__cov;
      // ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;

      return Op(false, "\nfails in trace!!!; error=" + test.error() + ss.str());
    } else
      return Op(Derivative<markov::mp_state_information>::adjust(
          std::move(P_mean), std::move(P__cov), std::move(y_mean),
          std::move(y_var), std::move(plogL), std::move(eplogL), vplogL,Q_dt.min_P(),
          e.f()));
  }

  template <class Model, class Step>
  myOptional_t<Derivative<markov::mp_state_information>>
  start(Model &m, const Step &p, double min_p) const {
    typedef myOptional_t<Derivative<markov::mp_state_information>> Op;

    auto P_mean = m.Peq(p.begin()->x());
    //  std::cerr<<"gslreijgsorjgps INIT!!!"<<P_mean.value();
    if (!P_mean)
      return Op(false, "fails to get Peq :" + P_mean.error());
    else {
      auto P_cov = diag(P_mean.value()) - quadraticForm_XTX(P_mean.value());
      //     std::cerr<<"gslreijgsorjgps INIT!!!"<<P_cov;
      auto nan = std::numeric_limits<double>::quiet_NaN();
      auto &x = P_mean.value().x();
      return Op(Derivative<markov::mp_state_information>::adjust(
          std::move(P_mean).value(), std::move(P_cov),
          Derivative<double>(nan, x), Derivative<double>(nan, x),
          Derivative<double>(nan, x), Derivative<double>(nan, x), min_p, 0));
    }
  }

  Derivative(double tolerance) : base_type{tolerance} {}
  Derivative() = default;
  };

  template <>
  class Derivative<markov::hidden::MacroDVNR> : public markov::hidden::MacroDVNR {
  public:
    inline constexpr static auto const className = my_static_string("MacroDVNR");
    typedef markov::hidden::MacroDVNR base_type;
    base_type const &f() const { return *this; }

  template <class Model, class Step>
    myOptional_t<Derivative<markov::mp_state_information>>
    run(const Derivative<markov::mp_state_information> &prior, Model &m,
        const Step &p, std::ostream &os) const {
    markov::MACROR alg;
    return run(prior, m, p, os, alg);
  }

  template <class Model, class Step, class MACROR>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior, Model &m,
      const Step &p, std::ostream &os, MACROR &alg) const {
    typedef myOptional_t<Derivative<markov::mp_state_information>> Op;
    auto Q_dto = m.get_P(p, 0);
    if (!Q_dto)
      return Op(false, "fails in auto Q_dt=m.get_P(p,0) :" + Q_dto.error());
    Derivative<Markov_Transition_step_double> Q_dt = std::move(Q_dto).value();
    return run(prior, Q_dt, m, p, os, alg);
  }

  template <class Model, class Step>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior,
      const Derivative<Markov_Transition_step_double> &Q_dt, Model &m,
      const Step &p, std::ostream &os, const markov::MACROR &alg) const {
    typedef myOptional_t<Derivative<markov::mp_state_information>> Op;
    switch (alg) {
    case markov::MACRO_DMNR:
      return Derivative<markov::hidden::MacroDMNR>(tolerance())
          .run(prior, Q_dt, m, p, os);
    case markov::MACRO_DVNR:
      return Derivative<markov::hidden::MacroDVNR>(tolerance(), Variance_magic_number())
          .run(prior, Q_dt, m, p, os);
    case markov::MACRO_DMR:
      return Op(false, "Not a contemplated algorithm :" +
                           markov::MACROR_string[markov::MACRO_DMR] +
                           " is not valid for " + className.str());
    case markov::MACRO_DVR:
    default:
      return Op(false, "Not a contemplated algorithm :" +
                           markov::MACROR_string[markov::MACRO_DVR] +
                           " is not valid for " + className.str());
    }
  }

  template <class Model, class Step>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior,
      const Derivative<Markov_Transition_step_double> &Q_dt, Model &m,
      const Step &p, std::ostream &os, markov::MACROR &alg) const {
    auto y = p.y();
    if (std::isnan(y))
      alg = markov::MACRO_DMNR;
    else {
      double mg = (prior.P_mean().f() * Q_dt.gmean_i().f()).getvalue();
      double g_max = max(Q_dt.g().f());
      double g_min = min(Q_dt.g().f());
      double g_range = g_max - g_min;
      double N = m.AverageNumberOfChannels().f();
      auto p_bi = (g_max - mg) / g_range;
      auto q_bi = (mg - g_min) / g_range;
      auto test_Binomial =
          markov::mp_state_information::is_Binomial_Approximation_valid(
              N, p_bi, q_bi, Variance_magic_number());
      if (!test_Binomial.has_value())
        alg = markov::MACRO_DMNR;
      else
        alg = markov::MACRO_DVNR;
    }
    const markov::MACROR alg2 = alg;
    return run(prior, Q_dt, m, p, os, alg2);
  }

  template <class Model, class Step>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior,
      const Derivative<Markov_Transition_step_double> &Q_dt, Model &m,
      const Step &p, std::ostream & /*os*/) const {
    typedef myOptional_t<Derivative<markov::mp_state_information>> Op;
    auto y = p.y();

    Derivative<double> N = m.AverageNumberOfChannels();

    Derivative<double> e = m.noise_variance(p.nsamples());
    M_Matrix<double> u(prior.P_mean().f().size(), 1, 1.0);

    auto SmD = prior.P_cov() - diag(prior.P_mean());

    Derivative<double> gSg =
        (TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.gmean_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_ij(), Q_dt.gmean_ij()) * u))
            .getvalue();

    Derivative<double> sSg =
        (TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.gmean_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_var_ij(), Q_dt.gmean_ij()) * u))
            .getvalue();
    Derivative<double> sSs =
        (TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.gvar_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_var_ij(), Q_dt.gvar_ij()) * u))
            .getvalue();
    auto sS = TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.P() +
              prior.P_mean() * Q_dt.gtotal_var_ij();
    auto gS = TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.P() +
              prior.P_mean() * Q_dt.gtotal_ij();

    Derivative<double> ms = (prior.P_mean() * Q_dt.gvar_i()).getvalue();
    Derivative<double> delta_emu = sqr(ms + e / N) - sSs / N * 2.0;
    delta_emu.f() = std::max(delta_emu.f(), 0.0);

    Derivative<double> ms0 = (ms - e / N) * 0.5 + sqrt(delta_emu) * 0.5;

    auto e_mu = e + N * ms0;
    auto y_mean =
        N * (prior.P_mean() * Q_dt.gmean_i()).getvalue() - N * 0.5 / e_mu * sSg;
    auto zeta = N / (2 * sqr(e_mu) + N * sSs);
    auto y_var = e_mu + N * gSg - N * zeta * sqr(sSg);
    y_var.f() = std::max(y_var.f(), e.f());

    auto dy = Constant(y) - y_mean;
    auto chi = dy / y_var;

    auto chi2 = dy * chi;

    Derivative<double> plogL(prior.y_mean().x());
    if (y_var.f() > 0)
      plogL = -0.5 * log(2 * PI * y_var) - 0.5 * chi2;
    else
      plogL.f() = std::numeric_limits<double>::infinity();

    Derivative<double> eplogL = -0.5 * log(2 * PI * y_var) +
                                Constant(-0.5); // e_mu+N*gSg"-N*zeta*sqr(sSg)"
    // double chilogL=(eplogL-plogL)/std::sqrt(0.5);

    double vplogL=0.5;
    //      std::cerr<<" \n\n----test----\n"<<test.error()<<"\n";
    auto P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P());
    auto P_mean = prior.P_mean() * Q_dt.P();
    P__cov += diag(P_mean);
    auto test = markov::mp_state_information::test(
        P_mean.f(), P__cov.f(), y_mean.f(), y_var.f(), plogL.f(), eplogL.f(),
        e.f(), tolerance());
    if (!test) {
      std::stringstream ss;
    ss << "\n step=" << p << " "
         << " y=" << y << " ymean=" << y_mean << " yvar=" << y_var << " e=" << e
         << " e_mu" << e_mu << " N*gSg" << N * gSg
         << " -N*zeta*sSg=" << N * zeta * sqr(sSg);
    ss << "\n ms0=(ms-e/N)/2+std::sqrt(delta_emu)/2;\nms0=" << ms0
       << "\n ms=" << ms << "\n e/N/2=" << e / N;

    ss << "\n std::sqrt(delta_emu)/2=" << sqrt(delta_emu) * 0.5;
    ss << "\ndelta_emu=sqr(ms+e/N)-2.0/N*sSs \nsqr(ms+e/N)="
       << sqr(ms + e / N) << " 2.0*sSs/N=" << 2.0 * sSs / N << " sSs=" << sSs;
    ss << "\nP_mean \n" << P_mean;
    ss << " N*zeta*sqr(sSg)=" << N * zeta * sqr(sSg) << " plogL=" << plogL
       << " eplogL=" << eplogL << "dif=" << plogL - eplogL << " sSs=" << sSs;
    ss << "\nPcov \n" << P__cov;
    // ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;

    return Op(false, "\nfails in trace!!!; error=" + test.error() + ss.str());
    } else
      return Op(Derivative<markov::mp_state_information>::adjust(
          std::move(P_mean), std::move(P__cov), std::move(y_mean),
          std::move(y_var), std::move(plogL), std::move(eplogL),vplogL, Q_dt.min_P(),
          e.f()));
  }

  template <class Model, class Step>
  myOptional_t<Derivative<markov::mp_state_information>>
  start(Model &m, const Step &p, double min_p) const {
    typedef myOptional_t<Derivative<markov::mp_state_information>> Op;

    auto P_mean = m.Peq(p.begin()->x());
    //  std::cerr<<"gslreijgsorjgps INIT!!!"<<P_mean.value();
    if (!P_mean)
      return Op(false, "fails to get Peq :" + P_mean.error());
    else {
      auto P_cov = diag(P_mean.value()) - quadraticForm_XTX(P_mean.value());
      //     std::cerr<<"gslreijgsorjgps INIT!!!"<<P_cov;
      auto nan = std::numeric_limits<double>::quiet_NaN();
      auto &x = P_mean.value().x();
      return Op(Derivative<markov::mp_state_information>::adjust(
          std::move(P_mean).value(), std::move(P_cov),
          Derivative<double>(nan, x), Derivative<double>(nan, x),
          Derivative<double>(nan, x), Derivative<double>(nan, x), min_p, 0));
    }
  }

  Derivative() = default;
  Derivative(double tolerance, double variance_magical)
      : base_type(tolerance, variance_magical) {}
  };

  template <> class Derivative<markov::hidden::MacroDMR> : public markov::hidden::MacroDMR {
  public:
    inline constexpr static auto const className = my_static_string("MacroDMR");

    typedef markov::hidden::MacroDMR base_type;
    base_type const &f() const { return *this; }

  template <class Model, class Step>
    myOptional_t<Derivative<markov::mp_state_information>>
    run(const Derivative<markov::mp_state_information> &prior, Model &m,
        const Step &p, std::ostream &os) const {
    markov::MACROR alg;
    return run(prior, m, p, os, alg);
  }

  template <class Model, class Step, class MACROR>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior, Model &m,
      const Step &p, std::ostream &os, MACROR &alg) const {
    typedef myOptional_t<Derivative<markov::mp_state_information>> Op;
    auto Q_dto = m.get_P(p, 0);
    if (!Q_dto)
      return Op(false, "fails in auto Q_dt=m.get_P(p,0) :" + Q_dto.error());

    assert((Derivative_correctness_mean_value_test(
        1e-2,1e4, [&p](auto &mo) { return std::move(mo.get_P(p, 0)).value(); }, m,
        Q_dto.value(), std::cerr, " step: ", p, "MACROR = ", alg)));
    std::vector<double> myps = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8};

    std::vector<std::decay_t<decltype(Q_dto.value())>> dQs;
    for (auto eps : myps)
      dQs.push_back(Incremental_ratio_model(
          eps, [&p](auto &mo) { return mo.get_P(p, 0).value(); }, m));

    auto fields = dQs[0].get_constructor_fields();
    auto f = [&Q_dto, &dQs, &myps](auto const &m) {
      std::cerr << m.idField << ": \n";
      std::cerr << "\n-------Derivative---------\n";
      std::cerr << std::invoke(m.access_method, Q_dto.value())
                << "-----------\n";

      for (std::size_t i = 0; i < myps.size(); ++i) {
        std::cerr << "\n-------" << myps[i] << "---------\n";
        std::cerr << std::invoke(m.access_method, dQs[i]) << "-----------\n";
      }
      return 0;
    };

    //  std::apply([&f](auto const &... m) { return (f(m) + ...); }, fields);

    // assert
    if constexpr (false)
      ((are_Equal_v(
          Q_dto.value(),
          Incremental_ratio_model(
              1e-3, [&p](auto &mo) { return mo.get_P(p, 0).value(); }, m),
          std::cerr)));
    Derivative<Markov_Transition_step_double> Q_dt = std::move(Q_dto).value();
    return run(prior, Q_dt, m, p, os, alg);
  }

  template <class Model, class Step>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior,
      const Derivative<Markov_Transition_step_double> Q_dt, Model &m,
      const Step &p, std::ostream &os, const markov::MACROR &alg) const {
    typedef myOptional_t<Derivative<markov::mp_state_information>> Op;
    switch (alg) {
    case markov::MACRO_DMNR:
      return Derivative < markov::hidden::MacroDMNR > (tolerance())
          .run(prior, Q_dt, m, p, os);
    case markov::MACRO_DVNR:
      return Op(false, "Not a contemplated algorithm :" +
                           markov::MACROR_string[markov::MACRO_DVNR] +
                           " is not valid for " + className.str());
    case markov::MACRO_DMR:
      return Derivative<markov::hidden::MacroDMR>(tolerance(), Binomial_magic_number())
          .run(prior, Q_dt, m, p, os);
    case markov::MACRO_DVR:
    default:
      return Op(false, "Not a contemplated algorithm :" +
                           markov::MACROR_string[markov::MACRO_DVR] +
                           " is not valid for " + className.str());
    }
  }

  template <class Model, class Step>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior,
      const Derivative<Markov_Transition_step_double> &Q_dt, Model &m,
      const Step &p, std::ostream &os, markov::MACROR &alg) const {
    auto y = p.y();
    if (std::isnan(y))
      alg = markov::MACRO_DMNR;
    else {
      double mg = (prior.P_mean().f() * Q_dt.gmean_i().f()).getvalue();
      double g_max = max(Q_dt.g().f());
      double g_min = min(Q_dt.g().f());
      double g_range = g_max - g_min;
      double N = m.AverageNumberOfChannels().f();
      auto p_bi = (g_max - mg) / g_range;
      auto q_bi = (mg - g_min) / g_range;
      auto test_Binomial =
          markov::mp_state_information::is_Binomial_Approximation_valid(
              N, p_bi, q_bi, Binomial_magic_number());
      if (!test_Binomial.has_value())
        alg = markov::MACRO_DMNR;
      else
        alg = markov::MACRO_DMR;
    }
    const markov::MACROR alg2 = alg;
    return run(prior, Q_dt, m, p, os, alg2);
  }

  template <class Model, class Step>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior,
      const Derivative<Markov_Transition_step_double> &Q_dt, Model &m,
      const Step &p, std::ostream & /*os*/) const {
    typedef myOptional_t<Derivative<markov::mp_state_information>> Op;
    auto y = p.y();

    // double mg=(prior.P_mean()*Q_dt.gmean_i()).getvalue();
    // double g_max=max(Q_dt.g());
    // double g_min=min(Q_dt.g());
    //     double g_range=g_max-g_min;
    Derivative<double> N = m.AverageNumberOfChannels();
    Derivative<double> e = m.noise_variance(p.nsamples());
    M_Matrix<double> u(prior.P_mean().f().size(), 1, 1.0);
    auto SmD = prior.P_cov() - diag(prior.P_mean());

    //    auto SmD = prior.P_cov() - diag(prior.P_mean());

    assert((are_Equal_v(SmD,
                        Incremental_ratio(1e-6,
                                          [](auto &cov, auto &mean) {
      return cov - diag(mean);
                                          },
                                          prior.P_cov(), prior.P_mean()),
                        std::cerr)));

    Derivative<double> gSg =
        (TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.gmean_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_ij(), Q_dt.gmean_ij()) * u))
            .getvalue();

    //    double gSg =
    //        (TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.gmean_i()).getvalue() +
    //        (prior.P_mean() * (elemMult(Q_dt.gtotal_ij(), Q_dt.gmean_ij()) *
    //        u))
    //            .getvalue();

    assert((are_Equal_v(
        gSg,
        Incremental_ratio(
            1e-6,
            [&u](auto const &gmean_i, auto const &SmD, auto const &P_mean,
                 auto const &gtotal_ij, auto const &gmean_ij) {
      return (TranspMult(gmean_i, SmD) * gmean_i).getvalue() +
             (P_mean * (elemMult(gtotal_ij, gmean_ij) * u)).getvalue();
            },
            Q_dt.gmean_i(), SmD, prior.P_mean(), Q_dt.gtotal_ij(),
            Q_dt.gmean_ij()),
        std::cerr)));

    //   auto mu_n=;
    auto gS = TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.P() +
              prior.P_mean() * Q_dt.gtotal_ij();

    Derivative<double> ms = (prior.P_mean() * Q_dt.gvar_i()).getvalue();

    auto e_mu = e + N * ms;
    auto y_mean = N * (prior.P_mean() * Q_dt.gmean_i()).getvalue();
    auto y_var = e_mu + N * gSg;
    auto dy = Constant(y) - y_mean;
    auto chi = dy / y_var;
    auto P_mean = prior.P_mean() * Q_dt.P() + chi * gS;
    assert(are_Equal_v(P_mean,
                       Incremental_ratio(1e-4,
                                         [](auto const &Pmean, auto const &P,
                                            auto const &chi_, auto const &gS_) {
      return Pmean * P + chi_ * gS_;
                                         },
                                         prior.P_mean(), Q_dt.P(), chi, gS),
                       std::cerr));

    auto P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P()) +
                  diag(prior.P_mean() * Q_dt.P()) -
                  (N / y_var) * quadraticForm_XTX(gS);

    assert(are_Equal_v(
        P__cov,
        Incremental_ratio(
            1e-4,
            [](auto const &Pmean, auto const &P, auto const &gS_,
               auto const &SmD_, auto const &N_, auto const &y_var_) {
      return quadraticForm_BT_A_B(SmD_, P) + diag(Pmean * P) -
                     (N_ / y_var_) * quadraticForm_XTX(gS_);
            },
            prior.P_mean(), Q_dt.P(), gS, SmD, N, y_var),
        std::cerr));

    auto chi2 = dy * chi;

    Derivative<double> plogL(y_mean.x());
    if (y_var.f() > 0)
      plogL = -0.5 * log(2.0 * PI * y_var) - 0.5 * chi2;
    else
      plogL.f() = std::numeric_limits<double>::infinity();

    /*

//
    //   auto mu_n=;
    auto gS = TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.P() +
              prior.P_mean() * Q_dt.gtotal_ij();

    double ms = (prior.P_mean() * Q_dt.gvar_i()).getvalue();

    auto e_mu = e + N * ms;
    auto y_mean = N * (prior.P_mean() * Q_dt.gmean_i()).getvalue();
    auto y_var = e_mu + N * gSg;
    auto dy = y - y_mean;
    auto chi = dy / y_var;
    auto P_mean = prior.P_mean() * Q_dt.P() + chi * gS;

    auto P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P()) +
                  diag(prior.P_mean() * Q_dt.P()) -
                  (N / y_var) * quadraticForm_XTX(gS);

    auto chi2 = dy * chi;

    double plogL;
    if (y_var > 0)
      plogL = -0.5 * log(2 * PI * y_var) - 0.5 * chi2;
    else
      plogL = std::numeric_limits<double>::infinity();
*/
    Derivative<double> eplogL = -0.5 * log(2.0 * PI * y_var) +
                                Constant(-0.5); // e_mu+N*gSg"-N*zeta*sqr(sSg)"
    // double chilogL=(eplogL-plogL)/std::sqrt(0.5);

    auto test = Derivative<markov::mp_state_information>::test(
        P_mean, P__cov, y_mean, y_var, plogL, eplogL, e.f(), tolerance());
    if (!test) {
      std::stringstream ss;

      ss << "\nP_mean \n" << P_mean;
      ss << "\nPcov \n" << P__cov;
      // ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;

      return Op(false, "\nfails in trace!!!; error=" + test.error() + ss.str());
    } else
      return Op(Derivative<markov::mp_state_information>::adjust(
          std::move(P_mean), std::move(P__cov), std::move(y_mean),
          std::move(y_var), std::move(plogL), std::move(eplogL), Q_dt.min_P(),
          e.f()));
  }

  template <class Model, class Step>
  myOptional_t<Derivative<markov::mp_state_information>>
  start(Model &m, const Step &p, double min_p) const {
    typedef myOptional_t<Derivative<markov::mp_state_information>> Op;

    auto P_mean = m.Peq(p.begin()->x());
    //  std::cerr<<"gslreijgsorjgps INIT!!!"<<P_mean.value();
    if (!P_mean)
      return Op(false, "fails to get Peq :" + P_mean.error());
    else {
      auto P_cov = diag(P_mean.value()) - quadraticForm_XTX(P_mean.value());
      //     std::cerr<<"gslreijgsorjgps INIT!!!"<<P_cov;
      auto nan = std::numeric_limits<double>::quiet_NaN();
      auto &x = P_mean.value().x();
      return Op(Derivative<markov::mp_state_information>::adjust(
          std::move(P_mean).value(), std::move(P_cov),
          Derivative<double>(nan, x), Derivative<double>(nan, x),
          Derivative<double>(nan, x), Derivative<double>(nan, x), min_p, 0));
    }
  }
  Derivative(double tolerance, double binomial_magical)
      : base_type(tolerance, binomial_magical) {}
  Derivative() = default;

  template <class Model, class Step>
  myOptional_t<Derivative<markov::mp_state_information>>
  start(Derivative<Model> &m, const Step &p, double min_p) const {
    typedef myOptional_t<Derivative<markov::mp_state_information>> Op;

    auto P_mean = m.Peq(p.begin()->x());
    //  std::cerr<<"gslreijgsorjgps INIT!!!"<<P_mean.value();
    if (!P_mean)
      return Op(false, "fails to get Peq :" + P_mean.error());
    else {
      auto P_cov = diag(P_mean.value()) - quadraticForm_XTX(P_mean.value());
      //     std::cerr<<"gslreijgsorjgps INIT!!!"<<P_cov;
      double nan = std::numeric_limits<double>::quiet_NaN();
      return Op(markov::mp_state_information::adjust(std::move(P_mean).value(),
                                                     std::move(P_cov), nan, nan,
                                                     nan, nan, min_p, 0));
    }
  }
  };

  template <>
  class Derivative<markov::hidden::MacroDVR> : public markov::hidden::MacroDVR

  {
  public:
    inline constexpr static auto const className = my_static_string("MacroDVR");
    typedef markov::hidden::MacroDVR base_type;
    base_type const &f() const { return *this; }

  template <class Model, class Step>
    myOptional_t<Derivative<markov::mp_state_information>>
    run(const Derivative<markov::mp_state_information> &prior, Model &m,
        const Step &p, std::ostream &os) const {
    markov::MACROR alg;
    return run(prior, m, p, os, alg);
  }

  template <class Model, class Step, class MACROR>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior, Model &m,
      const Step &p, std::ostream &os, MACROR &alg) const {
    typedef myOptional_t<Derivative<markov::mp_state_information>> Op;
    auto Q_dto = m.get_P(p, 0);
    if (!Q_dto)
      return Op(false, "fails in auto Q_dt=m.get_P(p,0) :" + Q_dto.error());
    Derivative<Markov_Transition_step_double> Q_dt = std::move(Q_dto).value();
    return run(prior, Q_dt, m, p, os, alg);
  }

  template <class Model, class Step>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior,
      const Derivative<Markov_Transition_step_double> Q_dt, Model &m,
      const Step &p, std::ostream &os, const markov::MACROR &alg) const {
    switch (alg) {
    case markov::MACRO_DMNR:
      return Derivative<markov::hidden::MacroDMNR>(tolerance())
          .run(prior, Q_dt, m, p, os);
    case markov::MACRO_DVNR:
      return Derivative<markov::hidden::MacroDVNR>(tolerance(), Variance_magic_number())
          .run(prior, Q_dt, m, p, os);
    case markov::MACRO_DMR:
      return Derivative<markov::hidden::MacroDMR>(tolerance(), Binomial_magic_number())
          .run(prior, Q_dt, m, p, os);
    case markov::MACRO_DVR:
    default:
      return Derivative<markov::hidden::MacroDVR>(tolerance(), Binomial_magic_number(),
                                                  Variance_magic_number())
          .run(prior, Q_dt, m, p, os);
    }
  }

  template <class Model, class Step>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior,
      const Derivative<Markov_Transition_step_double> &Q_dt, Model &m,
      const Step &p, std::ostream &os, markov::MACROR &alg) const {
    auto y = p.y();
    if (std::isnan(y))
      alg = markov::MACRO_DMNR;
    else {
      double mg = (prior.P_mean().f() * Q_dt.gmean_i().f()).getvalue();
      double g_max = max(Q_dt.g().f());
      double g_min = min(Q_dt.g().f());
      double g_range = g_max - g_min;
      double N = m.AverageNumberOfChannels().f();
      std::pair<Op_void, Op_void> test_Binomial;
      auto p_bi = (g_max - mg) / g_range;
      auto q_bi = (mg - g_min) / g_range;
      test_Binomial =
          markov::mp_state_information::is_Binomial_Approximation_valid(
              N, p_bi, q_bi, Binomial_magic_number(), 0);
      if (!test_Binomial.first.has_value()) {
        if (!test_Binomial.second.has_value())
          alg = markov::MACRO_DMNR;
        else
          alg = markov::MACRO_DVNR;
      } else if (!test_Binomial.second.has_value())
        alg = markov::MACRO_DMR;
      else
        alg = markov::MACRO_DVR;
    }
    const markov::MACROR alg2 = alg;
    return run(prior, Q_dt, m, p, os, alg2);
  }

  template <class Model, class Step>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior,
      const Derivative<Markov_Transition_step_double> &Q_dt, Model &m,
      const Step &p, std::ostream & /*os*/) const {
    auto y = p.y();
    typedef myOptional_t<Derivative<markov::mp_state_information>> Op;
    Derivative<double> N = m.AverageNumberOfChannels();

    Derivative<double> e = m.noise_variance(p.nsamples());
    //
    M_Matrix<double> u(prior.P_mean().f().size(), 1, 1.0);

    auto SmD = prior.P_cov() - diag(prior.P_mean());

    Derivative<double> gSg =
        (TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.gmean_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_ij(), Q_dt.gmean_ij()) * u))
            .getvalue();

    Derivative<double> sSg =
        (TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.gmean_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_var_ij(), Q_dt.gmean_ij()) * u))
            .getvalue();
    Derivative<double> sSs =
        (TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.gvar_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_var_ij(), Q_dt.gvar_ij()) * u))
            .getvalue();
    //   auto mu_n=;
    auto sS = TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.P() +
              prior.P_mean() * Q_dt.gtotal_var_ij();
    auto gS = TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.P() +
              prior.P_mean() * Q_dt.gtotal_ij();

    Derivative<double> ms = (prior.P_mean() * Q_dt.gvar_i()).getvalue();
    Derivative<double> delta_emu = sqr(ms + e / N) - sSs / N * 2.0;
    delta_emu.f() = std::max(delta_emu.f(), 0.0);

    Derivative<double> ms0 = (ms - e / N) * 0.5 + sqrt(delta_emu) * 2.0;

    auto e_mu = e + N * ms0;
    auto y_mean =
        N * (prior.P_mean() * Q_dt.gmean_i()).getvalue() - N * 0.5 / e_mu * sSg;
    auto zeta = N / (2 * sqr(e_mu) + N * sSs);
    auto y_var = e_mu + N * gSg - N * zeta * sqr(sSg);
    y_var.f() = std::max(y_var.f(), e.f());

    auto dy = Constant(y) - y_mean;
    auto chi = dy / y_var;
    auto P_mean = prior.P_mean() * Q_dt.P() + chi * gS -
                  (chi * zeta * sSg + 0.5 / e_mu) * sS;

    auto P__cov =
        quadraticForm_BT_A_B(SmD, Q_dt.P()) + diag(prior.P_mean() * Q_dt.P()) -
        (zeta + N / y_var * sqr(zeta * sSg)) * quadraticForm_XTX(sS) +
        (2.0 * N / y_var * zeta * sSg) * TransposeSum(TranspMult(sS, gS)) -
        (N / y_var) * quadraticForm_XTX(gS);

    auto chi2 = dy * chi;

    Derivative<double> plogL(prior.y_mean().x());
    if (y_var.f() > 0)
      plogL = -0.5 * log(2 * PI * y_var) - 0.5 * chi2;
    else
      plogL.f() = std::numeric_limits<double>::infinity();

    Derivative<double> eplogL = -0.5 * log(2 * PI * y_var) +
                                Constant(-0.5); // e_mu+N*gSg"-N*zeta*sqr(sSg)"
    // double chilogL=(eplogL-plogL)/std::sqrt(0.5);
    auto test = markov::mp_state_information::test(
        P_mean.f(), P__cov.f(), y_mean.f(), y_var.f(), plogL.f(), eplogL.f(),
        e.f(), tolerance());
    if (!test) {
      std::stringstream ss;
    ss << "\n step=" << p << " "
         << " y=" << y << " ymean=" << y_mean << " yvar=" << y_var << " e=" << e
         << " e_mu" << e_mu << " N*gSg" << N * gSg
         << " -N*zeta*sSg=" << N * zeta * sqr(sSg);
    ss << "\n ms0=(ms-e/N)/2+std::sqrt(delta_emu)/2;\nms0=" << ms0
       << "\n ms=" << ms << "\n e/N/2=" << e / N;

    ss << "\n std::sqrt(delta_emu)/2=" << sqrt(delta_emu) * 0.5;
    ss << "\ndelta_emu=sqr(ms+e/N)-2.0/N*sSs \nsqr(ms+e/N)="
       << sqr(ms + e / N) << " 2.0/N*sSs=" << 2.0 / N * sSs << " sSs=" << sSs;
    ss << "\nP_mean \n" << P_mean;
    ss << " -N*zeta*sqr(sSg)=" << N * zeta * sqr(sSg) << " plogL=" << plogL
       << " eplogL=" << eplogL << "dif=" << plogL - eplogL << " sSs=" << sSs;
    ss << "\nPcov \n" << P__cov;
    // ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;

    return Op(false, "\nfails in trace!!!; error=" + test.error() + ss.str());
    } else
      return Op(Derivative<markov::mp_state_information>::adjust(
          std::move(P_mean), std::move(P__cov), std::move(y_mean),
          std::move(y_var), std::move(plogL), std::move(eplogL), Q_dt.min_P(),
          e.f()));
  }

  template <class Model, class Step>
  myOptional_t<Derivative<markov::mp_state_information>>
  start(Model &m, const Step &p, double min_p) const {
    typedef myOptional_t<Derivative<markov::mp_state_information>> Op;

    auto P_mean = m.Peq(p.begin()->x());
    //  std::cerr<<"gslreijgsorjgps INIT!!!"<<P_mean.value();
    if (!P_mean)
      return Op(false, "fails to get Peq :" + P_mean.error());
    else {
      auto P_cov = diag(P_mean.value()) - quadraticForm_XTX(P_mean.value());
      //     std::cerr<<"gslreijgsorjgps INIT!!!"<<P_cov;
      auto nan = std::numeric_limits<double>::quiet_NaN();
      auto &x = P_mean.value().x();
      return Op(Derivative<markov::mp_state_information>::adjust(
          std::move(P_mean).value(), std::move(P_cov),
          Derivative<double>(nan, x), Derivative<double>(nan, x),
          Derivative<double>(nan, x), Derivative<double>(nan, x), min_p, 0));
    }
  }
  Derivative(double tolerance, double binomial_magical, double variance_magical)
      : base_type(tolerance, binomial_magical, variance_magical) {}
  Derivative() = default;
  };


  template <bool recursive, int averaging, bool variance>
  class Derivative<markov::Macro_R<recursive, averaging, variance>> : public markov::Macro_R<recursive, averaging, variance>

  {
  public:
    using markov::Macro_R<recursive, averaging, variance>::tolerance;
    using markov::Macro_R<recursive, averaging, variance>::Binomial_magic_number;
    using markov::Macro_R<recursive, averaging, variance>::Variance_magic_number;


    typedef markov::Macro_R<recursive, averaging, variance> base_type;
    inline constexpr static auto const className = my_trait<base_type>::className+my_static_string("_derivative");
    base_type const &f() const { return *this; }

  template <class Model, class Step>
    myOptional_t<Derivative<markov::mp_state_information>>
    run(const Derivative<markov::mp_state_information> &prior, Model &m,
        const Step &p, std::ostream &os) const {
    markov::MACROR alg;
    return run(prior, m, p, os, alg);
  }

  template <class Model, class Step, class MACROR>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior, Model &m,
      const Step &p, std::ostream &os, MACROR &alg) const {
    typedef myOptional_t<Derivative<markov::mp_state_information>> Op;
    auto Q_dto = m.get_P(p, 0);
    if (!Q_dto)
      return Op(false, "fails in auto Q_dt=m.get_P(p,0) :" + Q_dto.error());
    Derivative<Markov_Transition_step_double> Q_dt = std::move(Q_dto).value();
    return run(prior, Q_dt, m, p, os, alg);
  }

  template <class Model, class Step>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior,
      const Derivative<Markov_Transition_step_double> Q_dt, Model &m,
      const Step &p, std::ostream &os, const markov::MACROR &alg) const {
    switch (alg) {
    case markov::MACRO_DMNR:
      return Derivative<markov::MacroDMNR>(tolerance(), Binomial_magic_number(),
                                           Variance_magic_number())
          .run(prior, Q_dt, m, p, os);
    case markov::MACRO_DVNR:
      return Derivative<markov::MacroDVNR>(tolerance(), Binomial_magic_number(),
                                           Variance_magic_number())
          .run(prior, Q_dt, m, p, os);
    case markov::MACRO_DMR:
      return Derivative<markov::MacroDMR>(tolerance(), Binomial_magic_number(),
                                          Variance_magic_number())
          .run(prior, Q_dt, m, p, os);
    case markov::MACRO_DVR:
    default:
      return Derivative<markov::MacroDVR>(tolerance(), Binomial_magic_number(),
                                          Variance_magic_number())
          .run(prior, Q_dt, m, p, os);
    }
  }

  template <class Model, class Step>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior,
      const Derivative<Markov_Transition_step_double> &Q_dt, Model &m,
      const Step &p, std::ostream &os, markov::MACROR &alg) const {
    auto y = p.y();
    if (std::isnan(y))
      alg = markov::MACRO_DMNR;
    else {
      double mg = (prior.P_mean().f() * Q_dt.gmean_i().f()).getvalue();
      double g_max = max(Q_dt.g().f());
      double g_min = min(Q_dt.g().f());
      double g_range = g_max - g_min;
      double N = m.AverageNumberOfChannels().f();
      std::pair<Op_void, Op_void> test_Binomial;
      auto p_bi = (g_max - mg) / g_range;
      auto q_bi = (mg - g_min) / g_range;
      test_Binomial =
          markov::mp_state_information::is_Binomial_Approximation_valid(
              N, p_bi, q_bi, Binomial_magic_number(), Variance_magic_number());
      if (!recursive || !test_Binomial.first.has_value()) {
        if (!variance || !test_Binomial.second.has_value())
          alg = markov::MACRO_DMNR;
        else
          alg = markov::MACRO_DVNR;
      } else if (!variance || !test_Binomial.second.has_value())
        alg = markov::MACRO_DMR;
      else
        alg = markov::MACRO_DVR;


    }
    const markov::MACROR alg2 = alg;
    return run(prior, Q_dt, m, p, os, alg2);
  }

  template <class Model, class Step>
  myOptional_t<Derivative<markov::mp_state_information>>
  run(const Derivative<markov::mp_state_information> &prior,
      const Derivative<Markov_Transition_step_double> &Q_dt, Model &m,
      const Step &p, std::ostream & /*os*/) const {
    auto y = p.y();
    typedef myOptional_t<Derivative<markov::mp_state_information>> Op;
    Derivative<double> N = m.AverageNumberOfChannels();

    Derivative<double> e = m.noise_variance(p.nsamples());
    //
    M_Matrix<double> u(prior.P_mean().f().size(), 1, 1.0);

    auto SmD = prior.P_cov() - diag(prior.P_mean());

    assert((Derivative_correctness_mean_value_test(
        1e-2, 1e4,[&u](auto P_cov,auto P_mean) {
          return
              P_cov - diag(P_mean);
        },
        std::forward_as_tuple(prior.P_cov(), prior.P_mean()), SmD, std::cerr,
        " step=",p,"prior=",prior,"Q_dt= ", Q_dt, "  SmD=", SmD)));

    Derivative<double> gSg =
        (TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.gmean_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_ij(), Q_dt.gmean_ij()) * Constant(u)))
            .getvalue();

    assert((Derivative_correctness_mean_value_test(
        1e-2, 1e4,[&u](auto gmean_i,auto SmD_,auto P_mean,auto gtotal_ij,auto gmean_ij) {
      return
          (TranspMult(gmean_i, SmD_) *gmean_i).getvalue() +
          (P_mean * (elemMult(gtotal_ij, gmean_ij) * Constant(u)))
              .getvalue();
        },
        std::forward_as_tuple(Q_dt.gmean_i(), SmD,prior.P_mean(),Q_dt.gtotal_ij(), Q_dt.gmean_ij()), gSg, std::cerr,
        " step=",p,"prior=",prior,"Q_dt= ", Q_dt, "  SmD=", SmD)));



    Derivative<double> ms = (prior.P_mean() * Q_dt.gvar_i()).getvalue();

    assert((Derivative_correctness_mean_value_test(
        1e-2, 1e4,[&u](auto P_mean, auto gvar_i) {
          return
              (P_mean * gvar_i).getvalue();;
        },
        std::forward_as_tuple(prior.P_mean(),Q_dt.gvar_i()), ms, std::cerr,
        " step=",p,"prior=",prior,"Q_dt= ", Q_dt)));


    Derivative<double> e_mu;
    Derivative<double> y_mean;
    Derivative<double> y_var;

    Derivative<double> sSg;
    Derivative<double> sSs;
    Derivative<double> zeta;
    if constexpr ((!variance) && (!recursive)) {
      e_mu = e + N * ms;
      y_mean = N * (prior.P_mean() * Q_dt.gmean_i()).getvalue();
      y_var = e_mu + N * gSg;




    } else if constexpr (!variance && recursive) {
    auto gS = TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.P() +
                prior.P_mean() * Q_dt.gtotal_ij();

    Derivative<double> ms = (prior.P_mean() * Q_dt.gvar_i()).getvalue();

    e_mu = e + N * ms;
    y_mean = N * (prior.P_mean() * Q_dt.gmean_i()).getvalue();
    y_var = e_mu + N * gSg;
    }
    else // (variance && (recursive || !recursive))
    {
      sSg =
          (TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.gmean_i()).getvalue() +
          (prior.P_mean() * (elemMult(Q_dt.gtotal_var_ij(), Q_dt.gmean_ij()) * Constant(u)))
              .getvalue();

    assert((Derivative_correctness_mean_value_test(
        1e-2, 1e4,[&u](auto gvar_i,auto SmD_,auto gmean_i,auto P_mean,auto gtotal_var_ij,auto gmean_ij) {
        return
            (TranspMult(gvar_i, SmD_) *gmean_i).getvalue() +
            (P_mean * (elemMult(gtotal_var_ij, gmean_ij) * Constant(u)))
                .getvalue();
},
          std::forward_as_tuple(Q_dt.gvar_i(), SmD,Q_dt.gmean_i(),prior.P_mean(),Q_dt.gtotal_var_ij(), Q_dt.gmean_ij()), sSg, std::cerr,
          " step=",p,"prior=",prior,"Q_dt= ", Q_dt, "  SmD=", SmD)));

    sSs =
        (TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.gvar_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_var_ij(), Q_dt.gvar_ij()) * Constant(u)))
            .getvalue();


    assert((Derivative_correctness_mean_value_test(
        1e-2, 1e4,[&u](auto gvar_i,auto SmD_,auto P_mean,auto gtotal_var_ij,auto gvar_ij) {
        return
            (TranspMult(gvar_i, SmD_) *gvar_i).getvalue() +
            (P_mean * (elemMult(gtotal_var_ij, gvar_ij) * Constant(u)))
                .getvalue();
        },
        std::forward_as_tuple(Q_dt.gvar_i(), SmD,prior.P_mean(),Q_dt.gtotal_var_ij(), Q_dt.gvar_ij()), sSs, std::cerr,
        " step=",p,"prior=",prior,"Q_dt= ", Q_dt, "  SmD=", SmD)));


    Derivative<double> delta_emu = sqr(ms + e / N) - sSs / N * 2.0;

    assert((Derivative_correctness_mean_value_test(
        1e-2, 1e4,[&u](auto ms_, auto e_, auto N_, auto sSs_) {
      return
          sqr(ms_ + e_ / N_) - sSs_ / N_ * 2.0;        },
        std::forward_as_tuple(ms,e,N,sSs),
        delta_emu,
        std::cerr,
        " step=",p,"prior=",prior,"Q_dt= ", Q_dt, "  SmD=", SmD)));



     delta_emu.f() = std::max(delta_emu.f(), 0.0);



     Derivative<double> ms0;
     if (delta_emu.f()>0){
     ms0= (ms - e / N) * 0.5 + sqrt(delta_emu) * 0.5;
    assert((Derivative_correctness_mean_value_test(
        1e-2, 1e4,[](auto ms_, auto e_, auto N_, auto delta_emu_) {
        return
            (ms_ - e_ / N_) * 0.5 + sqrt(delta_emu_) * 0.5;        },
         std::forward_as_tuple(ms,e,N,delta_emu),
         ms0,
         std::cerr,
         " step=",p,"ms=",ms,"e=",e,"N=",N,"delta_emu=",delta_emu)));


     }else {
       ms0= (ms - e / N) * 0.5 ;

     }


    e_mu = e + N * ms0;
    y_mean =
        N * (prior.P_mean() * Q_dt.gmean_i()).getvalue() - N * 0.5 / e_mu * sSg;

    assert((Derivative_correctness_mean_value_test(
        1e-2, 1e4,[](auto N_, auto P_mean_, auto gmean_i_, auto e_mu_, auto sSg_) {
      return
              N_ * (P_mean_ * gmean_i_).getvalue() - N_ * 0.5 / e_mu_ * sSg_;        },
        std::forward_as_tuple(N,prior.P_mean(),Q_dt.gmean_i(),e_mu,sSg),
        y_mean,
        std::cerr,
        " step=",p,"emu=",e_mu, "sSg=",sSg)));


    zeta = N / (2 * sqr(e_mu) + N * sSs);


    assert((Derivative_correctness_mean_value_test(
        1e-2, 1e4,[](auto N_, auto e_mu_, auto sSs_) {
      return
          N_ / (2 * sqr(e_mu_) + N_ * sSs_);        },
        std::forward_as_tuple(N,e_mu,sSs),
        zeta,
        std::cerr,
        " step=",p,"zeta=",zeta)));




    y_var =e_mu + N * gSg - N * zeta * sqr(sSg);

    assert((Derivative_correctness_mean_value_test(
        1e-2, 1e4,[](auto N_, auto e_mu_, auto sSg_,auto gSg_,auto zeta_) {
          return
              e_mu_ + N_ * gSg_ - N_ * zeta_ * sqr(sSg_);        },
        std::forward_as_tuple(N,e_mu,sSg,gSg,zeta),
        y_var,
        std::cerr,
        " step=",p,"yvar=",y_var,"N=",N, " e_mu=", e_mu," sSg=",sSg," gSg=",gSg," zeta=",zeta)));





    //  y_var.f() = std::max(y_var.f(), e.f());
    }


    if (std::isnan(y)) {
      Derivative<double> plogL(0.0,
                               prior.y_mean().x());
      Derivative<double> eplogL(0.0,
                                prior.y_mean().x());
      double vplogL=0.0;
      auto P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P());
      assert(are_Equal_v(P__cov,
                         Incremental_ratio(1e-4,
                                           [](auto const &P, auto const &SmD_) {
      return quadraticForm_BT_A_B(SmD_,
                                                                         P);
                                           },
                                           Q_dt.P(), SmD),
                         std::cerr));
      auto P_mean = prior.P_mean() * Q_dt.P();
      assert(are_Equal_v(
          P_mean,
          Incremental_ratio(
              1e-4, [](auto const &Pmean, auto const &P) { return Pmean * P; },
              prior.P_mean(), Q_dt.P()),
          std::cerr));

      P__cov += diag(P_mean);

      // std::cerr<<"\nPcov nana corr\n"<<P__cov<<"\nP_mean nana
      // corr\n"<<P_mean<<"\nQ.P \n"<<Q_dt.P();
      auto test = markov::mp_state_information::test(P_mean.f(), P__cov.f(),
                                                     tolerance());
      if (test.has_value())
        return Op(Derivative<markov::mp_state_information>::adjust(
            std::move(P_mean), std::move(P__cov), std::move(y_mean),
            std::move(y_var), std::move(plogL), std::move(eplogL), vplogL,Q_dt.min_P(),
            e.f()));
      else
        return Op(false, "fails at intertrace prediction!!: " + test.error());
    }

    auto dy = Constant(y) - y_mean;
    auto chi = dy / y_var;

    assert((Derivative_correctness_mean_value_test(
        1e-2, 1e4,[&y](auto y_mean_, auto y_var_) {
          return
              (Constant(y) - y_mean_) / y_var_;
        },
        std::forward_as_tuple(y_mean,y_var),
        chi,
        std::cerr,
        " step=",p)));



    Derivative<M_Matrix<double>> P_mean;
    Derivative<M_Matrix<double>> P__cov;
    if constexpr (!recursive) {
      P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P());




      P_mean = prior.P_mean() * Q_dt.P();
      P__cov += diag(P_mean);



    assert((Derivative_correctness_mean_value_test(
        1e-2, 1e4,[&u](auto SmD_, auto P, auto P_mean) {
    return
        quadraticForm_BT_A_B(SmD_, P)+diag(P_mean);
    ;        },
        std::forward_as_tuple(SmD,Q_dt.P(),P_mean),
          P__cov,
          std::cerr,
          " step=",p
          // ,"prior=",prior,"Q_dt= ", Q_dt, "  SmD=", SmD
                                                                               )));



    } else if constexpr (!variance) {
      auto gS = TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.P() +
                prior.P_mean() * Q_dt.gtotal_ij();
      P_mean = prior.P_mean() * Q_dt.P() + chi * gS;

      P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P()) +
               diag(prior.P_mean() * Q_dt.P()) -
               (N / y_var) * quadraticForm_XTX(gS);

    } else {
      auto gS = TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.P() +
                prior.P_mean() * Q_dt.gtotal_ij();
      auto sS = TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.P() +
                prior.P_mean() * Q_dt.gtotal_var_ij();
      P_mean = prior.P_mean() * Q_dt.P() + chi * gS -
               (chi * zeta * sSg + 0.5 / e_mu) * sS;

    P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P()) +
             diag(prior.P_mean() * Q_dt.P()) -
               (zeta + N / y_var * sqr(zeta * sSg)) * quadraticForm_XTX(sS) +
               (2.0 * N / y_var * zeta * sSg) * TransposeSum(TranspMult(sS, gS)) -
               (N / y_var) * quadraticForm_XTX(gS);
    }

    auto chi2 = dy * chi;
    assert((Derivative_correctness_mean_value_test(
        1e-2, 1e4,[&y](auto dy_, auto chi_) {
          return
              dy_ * chi_;
        },
        std::forward_as_tuple(dy,chi),
        chi2,
        std::cerr,
        " step=",p)));


    Derivative<double> plogL(y_mean.x());
    if (y_var.f() > 0)
      plogL = -0.5 * log((2.0 * PI) * y_var) - 0.5 * chi2;
    else
      plogL.f() = std::numeric_limits<double>::infinity();

    double vplogL=0.5;
    assert((Derivative_correctness_mean_value_test(
      1e-2, 1e4,[](auto y_var_, auto chi2_) {
          return
              -0.5 * log((2.0 * PI) * y_var_) - 0.5 * chi2_;
          ;
        },
        std::forward_as_tuple(y_var,chi2),
        plogL,
        std::cerr,
        " step=",p)));




    Derivative<double> eplogL =
        -0.5 * log((2 * PI) * y_var) +Constant(-0.5);

    auto test = Derivative<markov::mp_state_information>::test(
        P_mean, P__cov, y_mean, y_var, plogL, eplogL, e.f(), tolerance());
    if (!test) {
      std::stringstream ss;

      ss << "\nP_mean \n" << P_mean;
      ss << "\nPcov \n" << P__cov;
      // ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;

      return Op(false, "\nfails in trace!!!; error=" + test.error() + ss.str());
    } else if constexpr (variance)
    return Op(Derivative<markov::mp_state_information>::adjust(
          std::move(P_mean), std::move(P__cov), std::move(y_mean),
          std::move(y_var), std::move(plogL), std::move(eplogL), vplogL,Q_dt.min_P(),
          e.f()));
    else {
      return Op(Derivative<markov::mp_state_information>(
          std::move(P_mean), std::move(P__cov), std::move(y_mean),
          std::move(y_var), std::move(plogL),std::move(eplogL),vplogL));

    }
    }
    template <class Model, class Step>
    myOptional_t<Derivative<markov::mp_state_information>>
    start(Model &m, const Step &p, double min_p) const {
      typedef myOptional_t<Derivative<markov::mp_state_information>> Op;

      auto P_mean = m.Peq(p.begin()->x());
      //  std::cerr<<"gslreijgsorjgps INIT!!!"<<P_mean.value();
      if (!P_mean)
        return Op(false, "fails to get Peq :" + P_mean.error());
      else {
        auto P_cov = diag(P_mean.value()) - quadraticForm_XTX(P_mean.value());
        //     std::cerr<<"gslreijgsorjgps INIT!!!"<<P_cov;
        auto nan = std::numeric_limits<double>::quiet_NaN();
        auto &x = P_mean.value().x();
      return Op(Derivative<markov::mp_state_information>::adjust(
          std::move(P_mean).value(), std::move(P_cov),
            Derivative<double>(nan, x), Derivative<double>(nan, x),
            Derivative<double>(nan, x), Derivative<double>(nan, x), 0.0,min_p, 0));
      }
    }
    Derivative(double tolerance, double binomial_magical, double variance_magical)
        : base_type(tolerance, binomial_magical, variance_magical) {}
    Derivative() = default;
    };




    template <> struct Derivative<markov::logLikelihood_function> {
      template <class Experiment>
      std::tuple<double, double, double, M_Matrix<double>, M_Matrix<double>>
      operator()(const Experiment &) const {
        return std::tuple(0, 0, 0, M_Matrix<double>(), M_Matrix<double>());
      }
      // template<class mp_state_information>
  void operator()(const Derivative<markov::mp_state_information> &mp,
                      std::tuple<double, double, double, M_Matrix<double>,
                                 M_Matrix<double>> &logLsum,
                      std::size_t &i) const {
    if (std::isfinite(mp.plogL().f())&& mp.plogL().f()!=0.0) {
      Derivative<Normal_Distribution<double>> n(mp.y_mean(), mp.y_var());
      std::get<0>(logLsum) += mp.plogL().f();
      std::get<1>(logLsum) += mp.eplogL().f();
      std::get<2>(logLsum) += 0.5;
      std::get<3>(logLsum) += mp.plogL().dfdx();
      std::get<4>(logLsum) += n.FIM();
      ++i;
    }
  }
    };

    template <> struct Derivative<markov::partialDistribution_function> {
      template <class Experiment> auto operator()(const Experiment &e) const {
        std::size_t n = e.num_measurements();

        return std::vector<Derivative<Normal_Distribution<double>>>(n);
      }

      // template<class mp_state_information>
      void operator()(const Derivative<markov::mp_state_information> &mp,
                      std::vector<Derivative<Normal_Distribution<double>>> &v,
                      std::size_t &i) const {
        v[i] = Derivative<Normal_Distribution<double>>(mp.y_mean(), mp.y_var());
        ++i;
      }
    };

    template <class aux, class PartialDlogLikelihood>
    struct Derivative_partialDistribution_function_aux {
      template <class Experiment>
      auto operator()(const Experiment &e, std::vector<aux> &v) const {
        std::size_t n = e.num_measurements();
        v.resize(n);
        return PartialDlogLikelihood(n);
      }

      void operator()(const Derivative<markov::mp_state_information> &mp,
                      PartialDlogLikelihood &v, std::size_t &i,
                      const aux &macror_alg) const {
        if constexpr (true) {
          Derivative<Normal_Distribution<double>> n(mp.y_mean(), mp.y_var());

      typename PartialDlogLikelihood::base_type logL(false,
                                                     mp.plogL().f(), mp.eplogL().f(),sqr(mp.plogL().f()- mp.eplogL().f()), mp.vplogL(), mp.plogL().dfdx(), -n.FIM(),quadraticForm_XTX(mp.plogL().dfdx()));
      v.add_one_i(i, logL, macror_alg);
      ++i;
        } else { // this fork does not make much sense for derivatives, or it does?
          if (std::isfinite(mp.plogL().f())) {
            Derivative<Normal_Distribution<double>> n(mp.y_mean(), mp.y_var());

            typename PartialDlogLikelihood::base_type logL(
                mp.plogL().f(), mp.eplogL().f(), 0.5, mp.plogL().dfdx(), n.FIM());
            v.add_one_i(i, logL, macror_alg);
            ++i;
          } else // hack to avoid including intervals in the Jacobian calculation
          {
            typename PartialDlogLikelihood::base_type logL(
                0, 0, 0, M_Matrix<double>(), M_Matrix<double>());
            v.add_one_i(i, logL, macror_alg);
            ++i;
          }
        }
      }


      void operator()(PartialDlogLikelihood &v)const {
        v.calc();
      }

    };

    namespace markov {

    template <class F, class MacroDR, class Model,
          template <class, class> class Experiment, class Point, class Measure,
              class... Aux>
    auto logLikelihood_experiment_calculation_derivative(
        const F &f, const MacroDR &a, Model &m, const Experiment<Point, Measure> &e,
        std::ostream &os, Aux &... aux) {
      auto out = f(e, aux...);
      typedef myOptional_t<std::decay_t<decltype(out)>> Op;

      auto first_step = *e.begin_begin();
      auto prior = a.start(m, first_step, m.min_P());

      if (!prior.has_value())
        return Op(false, "calculation interrupted at start :" + prior.error());
      assert(are_Equal_v(prior.value(), Incremental_ratio_model(
                                            1e-6,
                                            [&first_step, &a](auto &mo) {
    return std::move(a.f().start(mo, first_step, mo.min_P())).value();
                                            },
                                            m), std::cerr));

      std::size_t i = 0;
      for (auto it = e.begin_begin(); it != e.end_end(); ++it) {
        auto post = a.run(prior.value(), m, *it, os, aux[i]...);
        if (!post.has_value())
          return Op(false, "calculation interrupted at " + ToString(*it) + "  :" +
                               post.error());
        else {
          if constexpr (false) {
        auto dpost = Incremental_ratio_model(
            1e-6,
                [&a, &it, &os, &i, &aux...](auto &mo, auto const &p) {
                  return std::move(a.f().run(p, mo, *it, os, aux[i]...)).value();
                },
                m, prior.value());
        assert((are_Equal_v(post.value(), dpost, std::cerr, "i= ", i,
                            "  position: *it=", *it)));
          }
          if constexpr(false)(Derivative_correctness_mean_value_test(
          1e-2,1e4,
          [&a, &it, &os, &i, &aux...](auto mo, auto p) {
        return std::move(a.run(p, mo, *it, os, aux[i]...)).value();
          },
          std::forward_as_tuple(m, prior.value()), post.value(), std::cerr
          ,"i= ", i, "  position: *it=", *it," aux=",aux[i]..., "\n prior= ", prior.value(),"\n post= ", post.value()
                                                                                                                                                                                                                                  ));
      f(post.value(), out, i, aux[i]...);
      prior.value() = std::move(post).value();
        }
      }
      f(out);
      return Op(out);
    }

template <class MacroDR, class Model, class Experiment>
myOptional_t<
    std::tuple<double, double, double, M_Matrix<double>, M_Matrix<double>>>
    logLikelihood_derivative(const Derivative<MacroDR> &a, Model &m,
                             const Experiment e, std::ostream &os) {
  return logLikelihood_experiment_calculation_derivative(
      Derivative<logLikelihood_function>(), a, m, e, os);
}

template <class MacroDR, class Model, class Experiment>
myOptional_t<std::vector<Derivative<Normal_Distribution<double>>>>
partialDistribution_derivative(const Derivative<MacroDR> &a, Model &m,
                               const Experiment e) {
  return logLikelihood_experiment_calculation_derivative(
      Derivative<markov::partialDistribution_function>(), a, m, e);
}

template <class aux, class PartialDlogLikelihood, class DMacroDR, class Model,
          class Experiment>
auto partialDlikelihood_derivative(const DMacroDR &a, Model &m,
                                   const Experiment &e, std::ostream &os) {
  std::vector<aux> dummy_a;
  return logLikelihood_experiment_calculation_derivative(
      Derivative_partialDistribution_function_aux<aux, PartialDlogLikelihood>(),
      a, m, e, os, dummy_a);
}
} // namespace markov

#endif // LIKELIHOOD_MARKOV_PROCESS_DERIVATIVE_H
