#ifndef LIKELIHOOD_MARKOV_PROCESS_DERIVATIVE_H
#define LIKELIHOOD_MARKOV_PROCESS_DERIVATIVE_H

#include "likelihood_markov_process.h"
#include "matrixderivative.h"
#include "qmodel_derivative.h"
#include "mydistributions_derivative.h"

template <> class Derivative<markov::mp_state_information> {
  Derivative<M_Matrix<double>> P_mean_;
  Derivative<M_Matrix<double>> P_cov_;

  Derivative<double> y_mean_;
  Derivative<double> y_var_;
  Derivative<double> plogL_;
  Derivative<double> eplogL_;

public:
  Derivative(Derivative<M_Matrix<double>> &&P_mean__,
             Derivative<M_Matrix<double>> &&P_cov__,
             Derivative<double> &&y_mean__, Derivative<double> &&y_var__,
             Derivative<double> &&plogL__, Derivative<double> &&eplogL__)
      : P_mean_{std::move(P_mean__)}, P_cov_{std::move(P_cov__)},
        y_mean_{std::move(y_mean__)}, y_var_{std::move(y_var__)},
        plogL_{std::move(plogL__)}, eplogL_{std::move(eplogL__)} {}

  Derivative() = default;
  auto &P_mean() const { return P_mean_; }
  auto &P_cov() const { return P_cov_; }

  auto &y_mean() const { return y_mean_; }
  auto &y_var() const { return y_var_; }
  auto &plogL() const { return plogL_; }
  auto &eplogL() const { return eplogL_; }

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
        grammar::field(C<self_type>{}, "eplogL", &self_type::eplogL));
  }

  static Derivative<markov::mp_state_information>
  adjust(Derivative<M_Matrix<double>> &&P_mean__,
         Derivative<M_Matrix<double>> &&P_cov__, Derivative<double> &&y_mean__,
         Derivative<double> &&y_var__, Derivative<double> &&plogL__,
         Derivative<double> &&eplogL__, double min_p, double min_var) {
    return Derivative<markov::mp_state_information>(
        Derivative<Probability_distribution>::normalize(std::move(P_mean__),
                                                        min_p),
        Derivative<Probability_distribution_covariance>::normalize(
            std::move(P_cov__), min_p),
        std::move(y_mean__),
        Derivative<variance_value>::adjust(std::move(y_var__), min_var),
        std::move(plogL__), std::move(eplogL__));
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
                      const Derivative<double> &eplogL__,
                      double minvariance,
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








template <> class Derivative<markov::MacroDMNR> : public markov::MacroDMNR {
public:
  typedef markov::MacroDMNR base_type;

  base_type const & f()const { return *this;}


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
      return Derivative<MacroDMNR>(tolerance()).run(prior, Q_dt, m, p, os);
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
    Derivative<double> gSg =
        (TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.gmean_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_ij(), Q_dt.gmean_ij()) * u))
            .getvalue();

    Derivative<double> ms = (prior.P_mean() * Q_dt.gvar_i()).getvalue();

    auto e_mu = e + N * ms;
    auto y_mean = N * (prior.P_mean() * Q_dt.gmean_i()).getvalue();
    auto y_var = e_mu + N * gSg;
    auto y = p.y();
    if (std::isnan(y)) {
      Derivative<double> plogL(std::numeric_limits<double>::quiet_NaN(),
                               prior.y_mean().x());
      Derivative<double> eplogL(std::numeric_limits<double>::quiet_NaN(),
                                prior.y_mean().x());
      auto P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P());
      auto P_mean = prior.P_mean() * Q_dt.P();
      P__cov += diag(P_mean);
      // std::cerr<<"\nPcov nana corr\n"<<P__cov<<"\nP_mean nana
      // corr\n"<<P_mean<<"\nQ.P \n"<<Q_dt.P();
      auto test = markov::mp_state_information::test(P_mean.f(), P__cov.f(),
                                                     tolerance());
      if (test.has_value())
        return Op(Derivative<markov::mp_state_information>::adjust(
            std::move(P_mean), std::move(P__cov), std::move(y_mean),
            std::move(y_var), std::move(plogL), std::move(eplogL), Q_dt.min_P(),
            e.f()));
      else
        return Op(false, "fails at intertrace prediction!!: " + test.error());
    }
    auto dy = Constant(y) - y_mean;
    auto chi = dy / y_var;
    auto P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P());
    auto P_mean = prior.P_mean() * Q_dt.P();
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
      auto x = P_mean.value().x();
      return Op(Derivative<markov::mp_state_information>::adjust(
          std::move(P_mean).value(), std::move(P_cov),
          Derivative<double>(nan, x), Derivative<double>(nan, x),
          Derivative<double>(nan, x), Derivative<double>(nan, x), min_p, 0));
    }
  }

  Derivative(double tolerance) : base_type{tolerance} {}
  Derivative() = default;
};

template <> class Derivative<markov::MacroDVNR> : public markov::MacroDVNR {
public:
  inline constexpr static auto const className = my_static_string("MacroDVNR");
  typedef markov::MacroDVNR base_type;
  base_type const & f()const { return *this;}

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
      return Derivative<markov::MacroDMNR>(tolerance())
          .run(prior, Q_dt, m, p, os);
    case markov::MACRO_DVNR:
      return Derivative<markov::MacroDVNR>(tolerance(), Variance_magic_number())
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
      auto x = P_mean.value().x();
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

template <> class Derivative<markov::MacroDMR> : public markov::MacroDMR {
public:
  inline constexpr static auto const className = my_static_string("MacroDMR");

  typedef markov::MacroDMR base_type;
  base_type const & f()const { return *this;}

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
      return Derivative<markov::MacroDMNR>(tolerance())
          .run(prior, Q_dt, m, p, os);
    case markov::MACRO_DVNR:
      return Op(false, "Not a contemplated algorithm :" +
                           markov::MACROR_string[markov::MACRO_DVNR] +
                           " is not valid for " + className.str());
    case markov::MACRO_DMR:
      return Derivative<markov::MacroDMR>(tolerance(), Binomial_magic_number())
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

    Derivative<double> gSg =
        (TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.gmean_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_ij(), Q_dt.gmean_ij()) * u))
            .getvalue();

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

    auto P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P()) +
                  diag(prior.P_mean() * Q_dt.P()) -
                  (N / y_var) * quadraticForm_XTX(gS);

    auto chi2 = dy * chi;

    Derivative<double> plogL(y_mean.x());
    if (y_var.f() > 0)
      plogL = -0.5 * log(2.0 * PI * y_var) - 0.5 * chi2;
    else
      plogL.f() = std::numeric_limits<double>::infinity();

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
      auto x = P_mean.value().x();
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
class Derivative<markov::MacroDVR> : public markov::MacroDVR

{
public:
  inline constexpr static auto const className = my_static_string("MacroDVR");
  typedef markov::MacroDVR base_type;
  base_type const & f()const { return *this;}

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
      return Derivative<markov::MacroDMNR>(tolerance())
          .run(prior, Q_dt, m, p, os);
    case markov::MACRO_DVNR:
      return Derivative<markov::MacroDVNR>(tolerance(), Variance_magic_number())
          .run(prior, Q_dt, m, p, os);
    case markov::MACRO_DMR:
      return Derivative<markov::MacroDMR>(tolerance(), Binomial_magic_number())
          .run(prior, Q_dt, m, p, os);
    case markov::MACRO_DVR:
    default:
      return Derivative<MacroDVR>(tolerance(), Binomial_magic_number(),
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
    Derivative<double> delta_emu = sqr(ms + e / N) - sSs / N * 0.5;
    delta_emu.f() = std::max(delta_emu.f(), 0.0);

    Derivative<double> ms0 = (ms - e / N) * 0.5 + sqrt(delta_emu) * 0.5;

    auto e_mu = e + N * ms0;
    auto y_mean =
        N * (prior.P_mean() * Q_dt.gmean_i()).getvalue() - N * 0.5 / e_mu * sSg;
    auto zeta = N / (2 * sqr(e_mu) + N * sSs);
    auto y_var = e_mu + N * gSg - N * zeta * sqr(sSg);
    y_var.f() = std::max(y_var.f() , e.f());

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
      auto x = P_mean.value().x();
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
    if (std::isfinite(mp.plogL().f())) {
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
  auto operator()(const Experiment &e,  std::vector<aux> & v) const {
    std::size_t n = e.num_measurements();
    v.resize(n);
    return PartialDlogLikelihood(n);
  }

  void operator()(const Derivative<markov::mp_state_information> &mp,
                  PartialDlogLikelihood &v, std::size_t &i,
                  const aux &macror_alg) const {
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
};

namespace markov {

template <class F, class MacroDR, class Model,
          template <class, class> class Experiment, class Point, class Measure,
          class... Aux>
auto logLikelihood_experiment_calculation_derivative(const F &f, const MacroDR &a,
                                          Model &m,
                                          const Experiment<Point, Measure> &e,
                                          std::ostream &os, Aux &... aux) {
  auto out = f(e, aux...);
  typedef myOptional_t<std::decay_t<decltype(out)>> Op;

  auto first_step = *e.begin_begin();
  auto prior = a.start(m, first_step, m.min_P());
  assert(are_Equal_v(prior.value(),Incremental_ratio_model([&first_step,&m,&a](auto &&mo){return a.f().start(mo,first_step,m.min_P());},m,1e-6),os));

  if (!prior.has_value())
    return Op(false, "calculation interrupted at start :" + prior.error());
  std::size_t i = 0;
  for (auto it = e.begin_begin(); it != e.end_end(); ++it) {
    auto post = a.run(prior.value(), m, *it, os, aux[i]...);
    if (!post.has_value())
      return Op(false, "calculation interrupted at " + ToString(*it) + "  :" +
                           post.error());
    else {
      f(post.value(), out, i, aux[i]...);
      prior.value() = std::move(post).value();
    }
  }
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
