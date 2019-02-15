#ifndef MYEVIDENCE_H
#define MYEVIDENCE_H
#include "Matrix.h"
#include "myDistributions.h"
#include "mylikelihood.h"
#include "myparameters.h"
#include <algorithm>
#include <filesystem>
#include <iomanip>
namespace evidence {
inline std::vector<std::string>
concatenate_unique(std::vector<std::string> &&one,
                   std::vector<std::string> &&two) {
  // both one and two has all unique labels
  std::vector<std::string> out2;
  for (auto &&e : std::move(two)) {
    if (std::find(one.begin(), one.end(), e) == one.end())
      out2.push_back(std::move(e));
  }
  one.insert(one.end(), out2.begin(), out2.end());
  return std::move(one);
}

inline std::vector<std::string>
concatenate_add_prefix_if_same(std::vector<std::string> &&one,
                               std::vector<std::string> &&two,
                               const std::string &prefix) {
  for (auto &e : two) {
    if (std::find(one.begin(), one.end(), e) != one.end())
      e = prefix + e;
  }
  one.insert(one.end(), two.begin(), two.end());
  return std::move(one);
}

struct s_beta {
  constexpr static auto const title = my_static_string("beta");
};
struct s_thermo {
  constexpr static auto const title = my_static_string("thermo_");
};
struct s_prior {
  constexpr static auto const title = my_static_string("prior_");
};
struct s_lik {
  constexpr static auto const title = my_static_string("lik_");
};
struct s_LM_landa {
  constexpr static auto const title = my_static_string("LM_landa");
};
struct s_accepting {
  constexpr static auto const title = my_static_string("accepting");
};

struct s_stretch_alpha {
  constexpr static auto const title = my_static_string("stretch_alpha");
};
struct s_landa50 {
  constexpr static auto const title = my_static_string("landa50");
};
struct s_hill_coeff {
  constexpr static auto const title = my_static_string("hill_coeff");
};

struct i_beta {
  constexpr static auto const title = my_static_string("i_beta");
};
struct i_walker {
  constexpr static auto const title = my_static_string("i_walker");
};
struct s_scout_number {
  constexpr static auto const title = my_static_string("scout_number");
};
struct i_sample {
  constexpr static auto const title = my_static_string("i_sample");
};
struct s_sample {
  constexpr static auto const title = my_static_string("n_sample");
};

struct i_stretch_alpha {
  constexpr static auto const title = my_static_string("i_stretch_alpha");
};
struct i_LM_landa {
  constexpr static auto const title = my_static_string("i_LM_landa");
};

template <class Model, class Parameters_distribution> class Prior_Model {
public:
  typedef Prior_Model self_type;
  typedef Cs<Model> template_types;
  constexpr static auto const className =
      my_static_string("Prior_Model") + my_trait<template_types>::className;

  static auto get_constructor_fields() {
    return std::make_tuple(grammar::field(C<self_type>{},
                                          "Parameters_distribution",
                                          &self_type::Parameters_Distribution));
  }

  typedef M_Matrix<double> Parameters;
  auto const &Parameters_Distribution() const { return prior_; }

  auto name(std::size_t ipar) const {
    return Parameters_Distribution().name(ipar);
  }

  auto size() const { return Parameters_Distribution().size(); }
  Parameters sample(std::mt19937_64 &mt) const { return prior_.sample(mt); }
  myOptional_t<logLikelihood>
  compute_Likelihood(const M_Matrix<double> &x) const {
    typedef myOptional_t<logLikelihood> Op;
    double logL = prior_.logP(x);
    if (std::isfinite(logL)) {
      double logL = prior_.logP(x);
      double vlogL = prior_.vlogP(x);
      double elogL = prior_.expected_logP();
      double evlogL = prior_.variance_logP();
      return Op(logLikelihood(logL, elogL, vlogL, evlogL));
    } else
      return Op(false, "not a finite value=" + ToString(logL));
  }
  myOptional_t<DlogLikelihood>
  compute_DLikelihood(const M_Matrix<double> &x) const {
    typedef myOptional_t<DlogLikelihood> Op;
    double logL = prior_.logP(x);
    double vlogL = prior_.vlogP(x);
    double elogL = prior_.expected_logP();
    double evlogL = prior_.variance_logP();

    auto G = prior_.dlogL_dx(x);
    auto H = prior_.dlogL_dx2(x);
    std::stringstream ss;
    if (are_finite<true, double>().test(logL, ss) &&
        are_finite<true, M_Matrix<double>>().test(G, ss) &&
        are_finite<true, M_Matrix<double>>().test(H, ss))
      return Op(DlogLikelihood(logL, elogL, vlogL, evlogL, std::move(G),
                               std::move(H)));
    else
      return Op(false, ss.str());
  }

  Prior_Model(const Parameters_distribution &prior) : prior_{prior} {}
  Prior_Model() = default;

private:
  Parameters_distribution prior_;
};

class ThlogLikelihood : public logLikelihood {
public:
  typedef logLikelihood base_type;
  typedef ThlogLikelihood self_type;
  constexpr static auto const className = my_static_string("ThlogLikelihood");

  double beta() const { return beta_; }
  auto &prior() const { return prior_; }
  auto &likelihood() const { return lik_; }
  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "prior", &self_type::prior),
        grammar::field(C<self_type>{}, "likelihood", &self_type::likelihood),
        grammar::field(C<self_type>{}, "beta", &self_type::beta));
  }

  ThlogLikelihood(logLikelihood prior, logLikelihood lik, double beta)
      : logLikelihood(prior.logL() + beta * lik.logL(),
                      prior.elogL() + beta * lik.elogL(),
                      prior.vlogL() + beta * lik.vlogL(),
                      prior.evlogL() + beta * lik.evlogL()),
        prior_{prior}, lik_{lik}, beta_{beta} {}
  ThlogLikelihood() = default;
  void set_beta(double beta) {
    beta_ = beta;
    set_logL(prior().logL() + beta * likelihood().logL(),
             prior().elogL() + beta * likelihood().elogL(),
             prior().vlogL() + beta * likelihood().vlogL(),
             prior().evlogL() + beta * likelihood().evlogL());
  }
  static auto get_data_index_static() {
    using namespace std::literals::string_literals;
    return Concatenate_tuple_static(
        std::make_tuple(make_data_static(
            std::tuple<>(), std::make_tuple(F_s(s_beta{}, &self_type::beta)))),
        Compose_static(
            F_s(s_thermo{},
                [](const self_type &s) -> base_type const & { return s; }),
            base_type::get_data_index_static()),
        Compose_static(F_s(s_prior{}, [](const self_type &s) -> auto const & {
                         return s.prior();
                       }),
                       base_type::get_data_index_static()),
        Compose_static(F_s(s_lik{}, [](const self_type &s) -> auto const & {
                         return s.likelihood();
                       }),
                       base_type::get_data_index_static()));
  }

private:
  logLikelihood prior_;
  logLikelihood lik_;
  double beta_;
};

class ThDlogLikelihood : public DlogLikelihood {
public:
  typedef DlogLikelihood base_type;

  constexpr static auto const className =
      my_static_string("ThDlogLikelihood_") + base_type::className;
  // std::string myClass()const  { return className.str();}

  typedef ThDlogLikelihood self_type;
  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "prior", &self_type::prior),
        grammar::field(C<self_type>{}, "likelihood", &self_type::likelihood),
        grammar::field(C<self_type>{}, "beta", &self_type::beta));
  }

  DlogLikelihood const &prior() const { return prior_; }
  DlogLikelihood const &likelihood() const { return lik_; }
  double beta() const { return beta_; }
  ThDlogLikelihood(const DlogLikelihood &prior, const DlogLikelihood &lik,
                   double beta)
      : DlogLikelihood(prior.logL() + beta * lik.logL(),
                       prior.elogL() + beta * lik.elogL(),
                       prior.vlogL() + beta * lik.vlogL(),
                       prior.evlogL() + beta * lik.evlogL(),

                       prior.G() + lik.G() * beta, prior.H() + lik.H() * beta),
        prior_{prior}, lik_{lik}, beta_{beta} {}
  ThDlogLikelihood(DlogLikelihood &&prior, DlogLikelihood &&lik, double beta)
      : DlogLikelihood(prior.logL() + beta * lik.logL(),
                       prior.elogL() + beta * lik.elogL(),
                       prior.vlogL() + beta * lik.vlogL(),
                       prior.evlogL() + beta * lik.evlogL(),
                       prior.G() + lik.G() * beta, prior.H() + lik.H() * beta),
        prior_{std::move(prior)}, lik_{std::move(lik)}, beta_{beta} {}
  ThDlogLikelihood() = default;
  void set_beta(double beta) {
    beta_ = beta;
    set_logL(prior().logL() + beta * likelihood().logL(),
             prior().elogL() + beta * likelihood().elogL(),
             prior().vlogL() + beta * likelihood().vlogL(),
             prior().evlogL() + beta * likelihood().evlogL());
    set_G(prior().G() + likelihood().G() * beta);
    set_H(prior().H() + likelihood().H() * beta);
  }
  static auto get_data_index_static() {
    using namespace std::literals::string_literals;
    return Concatenate_tuple_static(
        std::make_tuple(make_data_static(
            std::tuple<>(), std::make_tuple(F_s(s_beta{}, &self_type::beta)))),
        Compose_static(
            F_s(s_thermo{},
                [](const self_type &s) -> base_type const & { return s; }),
            base_type::get_data_index_static()),
        Compose_static(F_s(s_prior{}, [](const self_type &s) -> auto const & {
                         return s.prior();
                       }),
                       base_type::get_data_index_static()),
        Compose_static(F_s(s_lik{}, [](const self_type &s) -> auto const & {
                         return s.likelihood();
                       }),
                       base_type::get_data_index_static()));
  }

private:
  DlogLikelihood prior_;
  DlogLikelihood lik_;
  double beta_;
};

template <class... Aux> class ThPartialDlogLikelihood : public DlogLikelihood {
public:
  typedef DlogLikelihood base_type;

  constexpr static auto const className =
      my_static_string("ThPartialDlogLikelihood_") + base_type::className;
  // std::string myClass()const  { return className.str();}

  typedef ThPartialDlogLikelihood self_type;
  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "prior", &self_type::prior),
        grammar::field(C<self_type>{}, "likelihood", &self_type::likelihood),
        grammar::field(C<self_type>{}, "beta", &self_type::beta));
  }

  auto &prior() const { return prior_; }
  auto &likelihood() const { return lik_; }
  double beta() const { return beta_; }
  ThPartialDlogLikelihood(const DlogLikelihood &prior,
                          const PartialDLogLikelihood<Aux...> &lik, double beta)
      : DlogLikelihood(prior.logL() + beta * lik.logL(),
                       prior.elogL() + beta * lik.elogL(),
                       prior.vlogL() + beta * lik.vlogL(),
                       prior.evlogL() + beta * lik.evlogL(),

                       prior.G() + lik.G() * beta, prior.H() + lik.H() * beta),
        prior_{prior}, lik_{lik}, beta_{beta} {}
  ThPartialDlogLikelihood(DlogLikelihood &&prior,
                          PartialDLogLikelihood<Aux...> &&lik, double beta)
      : DlogLikelihood(prior.logL() + beta * lik.logL(),
                       prior.elogL() + beta * lik.elogL(),
                       prior.vlogL() + beta * lik.vlogL(),
                       prior.evlogL() + beta * lik.evlogL(),
                       prior.G() + lik.G() * beta, prior.H() + lik.H() * beta),
        prior_{std::move(prior)}, lik_{std::move(lik)}, beta_{beta} {}
  ThPartialDlogLikelihood() = default;
  void set_beta(double beta) {
    beta_ = beta;
    set_logL(prior().logL() + beta * likelihood().logL(),
             prior().elogL() + beta * likelihood().elogL(),
             prior().vlogL() + beta * likelihood().vlogL(),
             prior().evlogL() + beta * likelihood().evlogL());
    set_G(prior().G() + likelihood().G() * beta);
    set_H(prior().H() + likelihood().H() * beta);
  }

  static auto get_data_index_static() {
    typedef std::decay_t<decltype(std::declval<self_type>().likelihood())>
        PartialD;
    return Concatenate_tuple_static(
        std::make_tuple(make_data_static(
            std::tuple<>(), std::make_tuple(F_s(s_beta{}, &self_type::beta)))),
        Compose_static(
            F_s(s_thermo{},
                [](const self_type &s) -> base_type const & { return s; }),
            base_type::get_data_index_static()),
        Compose_static(F_s(s_prior{}, [](const self_type &s) -> auto const & {
                         return s.prior();
                       }),
                       base_type::get_data_index_static()),
        Compose_static(F_s(s_lik{}, [](const self_type &s) -> auto const & {
                         return s.likelihood();
                       }),
                       PartialD::get_data_index_static()));
  }

private:
  DlogLikelihood prior_;
  PartialDLogLikelihood<Aux...> lik_;
  double beta_;
};

template <class PriorModel, class LikelihoodModel> class Thermodynamic_Model {
  static_assert(std::is_same_v<typename PriorModel::Parameters,
                               typename LikelihoodModel::Parameters>);

public:
  typedef Thermodynamic_Model self_type;
  typedef Cs<PriorModel, LikelihoodModel> template_types;
  constexpr static auto const className =
      my_static_string("Thermodynamic_Model_") +
      my_trait<template_types>::className;

  typedef typename PriorModel::Parameters Parameters;

  Parameters sample(std::mt19937_64 &mt) const { return prior().sample(mt); }

  typedef ThDlogLikelihood DLikelihoodResult;
  typedef ThlogLikelihood LikelihoodResult;

  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "prior", &self_type::prior),
        grammar::field(C<self_type>{}, "likelihood", &self_type::likelihood),
        grammar::field(C<self_type>{}, "beta", &self_type::beta));
  }

  Thermodynamic_Model(const PriorModel &prior, const LikelihoodModel &lik,
                      double beta)
      : prior_{prior}, lik_{lik}, beta_{beta} {}
  //  Thermodynamic_Model(const Thermodynamic_Model &other)
  //      : prior_{other.prior_}, lik_{other.lik_}, beta_{other.beta_} {};
  //  Thermodynamic_Model(Thermodynamic_Model &&other)
  //      : prior_{other.prior_}, lik_{other.lik_}, beta_{
  //                                                    std::move(other.beta_)}
  //                                                    {};
  //  Thermodynamic_Model &operator=(Thermodynamic_Model &&) = default;
  Thermodynamic_Model() = default;
  double beta() const { return beta_; }

  std::string name(std::size_t i) const { return prior().name(i); }

  auto size() const { return prior().size(); }

  void set_beta(double _beta) { beta_ = _beta; }

  template <class mySample>
  auto compute_DLikelihood(const Parameters &x, const mySample &current,
                           std::ostream &os) const {
    typedef myOptional_t<ThDlogLikelihood> Op;
    auto p = prior_.compute_DLikelihood(x);
    if (!p.has_value())
      return Op(false, "invalid prior DlogLik :" + p.error());
    auto l = lik_.compute_DLikelihood(x, current, os);
    if (!l.has_value())
      return Op(false, "fails in likelihood :" + l.error());
    else
      return Op(
          ThDlogLikelihood(std::move(p).value(), std::move(l).value(), beta()));
  }

  template <class Parameters>
  auto compute_DLikelihood_init(const Parameters &x, std::ostream &os) const {
    typedef myOptional_t<ThDlogLikelihood> Op;
    auto p = prior_.compute_DLikelihood(x);
    if (!p.has_value())
      return Op(false, "invalid prior DlogLik :" + p.error());
    auto l = lik_.compute_DLikelihood_init(x, os);
    if (!l.has_value())
      return Op(false, "fails in likelihood :" + l.error());
    else
      return Op(
          ThDlogLikelihood(std::move(p).value(), std::move(l).value(), beta()));
  }

  auto compute_Likelihood(const Parameters &x) const {
    typedef myOptional_t<ThlogLikelihood> Op;
    auto p = prior_.compute_Likelihood(x);
    if (!p.has_value())
      return Op(false, "invalid prior logLik :" + p.error());
    auto l = lik_.compute_Likelihood(x);
    if (!l.has_value())
      return Op(false, "fails in likelihood :" + l.error());
    else
      return Op(
          ThlogLikelihood(std::move(p).value(), std::move(l).value(), beta()));
  }

  template <class Sample>
  auto compute_PartialDLikelihood(const Parameters &x, const Sample &current,
                                  std::ostream &os) const {
    auto p = prior_.compute_DLikelihood(x);
    auto l = lik_.compute_PartialDLikelihood(x, current, os);
    typedef myOptional_t<decltype(ThPartialDlogLikelihood(
        std::move(p).value(), std::move(l).value(), beta()))>
        Op;
    if (!p.has_value())
      return Op(false, "invalid prior DlogLik :" + p.error());
    if (!l.has_value())
      return Op(false, "fails in likelihood :" + l.error());
    else
      return Op(ThPartialDlogLikelihood(std::move(p).value(),
                                        std::move(l).value(), beta()));
  }

  PriorModel const &prior() const { return prior_; };
  LikelihoodModel const &likelihood() const { return lik_; }

private:
  PriorModel prior_;
  LikelihoodModel lik_;
  double beta_;
};

template <class PriorModel, class LikelihoodModel>
class Thermodynamic_Model_Series {
  static_assert(std::is_same_v<typename PriorModel::Parameters,
                               typename LikelihoodModel::Parameters>);

public:
  typedef Thermodynamic_Model_Series self_type;
  typedef Cs<PriorModel, LikelihoodModel> template_types;
  constexpr static auto const className =
      my_static_string("Thermodynamic_Model_Series") +
      my_trait<template_types>::className;

public:
  typedef typename PriorModel::Parameters Parameters;

  typedef Thermodynamic_Model<PriorModel, LikelihoodModel> subModel;

  typedef ThDlogLikelihood DLikelihoodResult;
  typedef ThlogLikelihood LikelihoodResult;
  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "prior", &self_type::prior),
        grammar::field(C<self_type>{}, "likelihood", &self_type::likelihood),
        grammar::field(C<self_type>{}, "betas", &self_type::betas));
  }

  Thermodynamic_Model_Series(const PriorModel &prior,
                             const LikelihoodModel &lik,
                             std::vector<double> beta)
      : prior_{prior}, lik_{lik}, beta_{beta}, models_{getModels(prior_, lik_,
                                                                 beta)} {}

  double beta(std::size_t i) const { return models_[i].beta(); }

  auto &betas() const { return beta_; }

  auto name(std::size_t i) const { return model(0).name(i); }

  auto name_size() const { return model(0).size(); }

  std::size_t size() const { return models_.size(); }

  void set_beta(double abeta, std::size_t i) {
    beta_[i] = abeta;
    models_[i].set_beta(abeta);
  }

  Thermodynamic_Model<PriorModel, LikelihoodModel> &model(std::size_t i) {
    return models_[i];
  }
  Thermodynamic_Model<PriorModel, LikelihoodModel> const &
  model(std::size_t i) const {
    return models_[i];
  }

  Thermodynamic_Model_Series(const Thermodynamic_Model_Series &) = default;

  auto &prior() const { return prior_; }
  auto &likelihood() const { return lik_; }
  Thermodynamic_Model_Series() = default;

private:
  PriorModel prior_;
  LikelihoodModel lik_;
  std::vector<double> beta_;
  std::vector<Thermodynamic_Model<PriorModel, LikelihoodModel>> models_;

  static std::vector<Thermodynamic_Model<PriorModel, LikelihoodModel>>
  getModels(const PriorModel &prior, const LikelihoodModel &lik,
            std::vector<double> beta) {
    std::vector<Thermodynamic_Model<PriorModel, LikelihoodModel>> out;
    for (std::size_t i = 0; i < beta.size(); ++i)
      out.push_back(Thermodynamic_Model<PriorModel, LikelihoodModel>(prior, lik,
                                                                     beta[i]));
    return out;
  }
};

template <class Model> class LevenbergMarquardt {
public:
  typedef typename Model::Parameters Parameters;

  class myDistribution : public Normal_Distribution<Parameters> {
  public:
    typedef Normal_Distribution<Parameters> base_type;
    constexpr static auto const className =
        my_static_string("Normal_Distribution_landa") +
        my_trait<Parameters>::className;
    std::string myClass() const override { return className.str(); }
    virtual myDistribution *clone() const override {
      return new myDistribution(*this);
    };

    typedef myDistribution self_type;
    static auto get_constructor_fields() {
      return std::make_tuple(
          grammar::field(C<self_type>{}, "Normal",
                         &self_type::get_Normal_Distribution),
          grammar::field(C<self_type>{}, "landa", &self_type::landa)

      );
    }
    myDistribution(Normal_Distribution<Parameters> &&d, double landa)
        : base_type{std::move(d)}, landa_{landa} {}

    myDistribution(const Normal_Distribution<Parameters> &d, double landa)
        : base_type{d}, landa_{landa} {}

    double landa() const { return landa_; }
    myDistribution() = default;
    static Data_Index_scheme data_index() {
      Data_Index_scheme out;
      out.push_back("mean", {"i_Parameter"});
      out.push_back("Covariance", {"i_Parameter", "i_Parameter_T"});
      out.push_back("logDetCov", {});
      out.push_back("landa", {});
      return out;
    }

    static auto get_data_index_static()

    {
      using namespace std::literals::string_literals;
      return Concatenate_tuple_static(
          std::make_tuple(make_data_static(
              std::tuple<>(),
              std::make_tuple(
                  F_s(CT_s<s_propDist, s_LM_landa>{}, &self_type::landa)))),
          Compose_static(
              F_s(s_propDist{},
                  [](const self_type &s) -> base_type const & { return s; }),
              base_type::get_data_index_static(i_Parameter{})));
    }

    base_type const &get_Normal_Distribution() const { return *this; }

  private:
    double landa_;
  };

  typedef myDistribution Distribution;
  typedef typename Model::DLikelihoodResult LikelihoodResult;
  typedef LevenbergMarquardt self_type;
  typedef Cs<Model> template_types;
  typedef s_LM_landa label_type;
  typedef i_LM_landa index_type;

  constexpr static auto const className =
      my_static_string("LevenbergMarquardt") +
      my_trait<template_types>::className;

  // std::string myClass()const  { return className.str();}

  static std::string index() { return "landa_LM"; }

  struct myAcceptProb {
    constexpr static auto const className =
        my_trait<LevenbergMarquardt>::className +
        my_static_string("_myAcceptProb");
    static std::tuple<> get_constructor_fields() { return std::tuple<>(); }

    static Data_Index_scheme data_index() { return {}; }
    struct myParameters {
      double first;
      double second;
      double landa_50() const { return first; }
      double hill_coeff() const { return second; }

      typedef i_Fit_Parameter index_type;
      typedef myParameters self_type;

      myParameters(double l50, double myh) : first(l50), second(myh) {}
      myParameters() = default;

      constexpr static auto className = my_static_string("Accept_parameters");
      bool operator<(const myParameters &other) const {
        if (first < other.first)
          return true;
        else if (first > other.first)
          return false;
        else
          return second < other.second;
      }
      static Data_Index_scheme data_index() {
        Data_Index_scheme out;
        out.push_back("landa_50", {"i_hill", "i_landa50"});
        out.push_back("hill_coeff", {"i_hill", "i_landa50"});
        return out;
      }

      static auto get_data_index_static() {
        return std::make_tuple(make_data_static(
            std::tuple<>(),
            std::make_tuple(F_s(s_LM_landa{}, &self_type::landa_50),
                            F_s(s_hill_coeff{}, &self_type::hill_coeff))));
      }

      static auto get_constructor_fields() {
        return std::make_tuple(
            grammar::field(C<self_type>{}, "landa_50", &self_type::landa_50),
            grammar::field(C<self_type>{}, "hill_coeff",
                           &self_type::hill_coeff));
      }
    };

    typedef myParameters Parameters;
    double operator()(const LevenbergMarquardt<Model> &LM,
                      const Parameters &param) const {
      double landa50 = param.first;
      double h = param.second;
      return 1.0 / (1.0 + std::pow(landa50 / (LM.landa() + 1.0), h));
    }
    double operator()(const Parameters &param,
                      const LevenbergMarquardt<Model> &LM) const {
      return operator()(LM, param);
    }
    static std::map<Parameters, double>
    uniform_parameter_prior(const std::vector<std::vector<double>> &v,
                            double p = -1) {
      std::map<Parameters, double> out;
      if (p == -1)
        p = 1.0 / (v[0].size() * v[1].size());
      for (std::size_t i = 0; i < v[0].size(); ++i)
        for (std::size_t j = 0; j < v[1].size(); ++j)
          out[{v[0][i], v[1][j]}] += p;
      return out;
    }
  };
  typedef myAcceptProb AcceptanceProbability;

  struct myExpectVelocity {
    constexpr static auto const className =
        my_trait<LevenbergMarquardt>::className +
        my_static_string("_myExpectVelocity");
    static std::tuple<> get_constructor_fields() { return std::tuple<>(); }
    static Data_Index_scheme data_index() { return {}; }
    double operator()(const LevenbergMarquardt<Model> &LM) const {
      return (1.0 / (1.0 + LM.landa()));
    }

    std::string operator()() const { return className.str(); }
  };

  typedef myExpectVelocity ExpectedVelocity;

  double landa() const { return landa_; }

  double operator()() const { return landa(); }

  myOptional<Distribution, reg_tag>
  calculate_Distribution(const Model &, const Parameters x,
                         const LikelihoodResult &logL) const {
    typedef myOptional<Distribution, reg_tag> Op;
    auto &H = logL.H();
    auto &G = logL.G();
    auto inv_cov = -(H + diag(H) * landa_);
    auto cov = inv(inv_cov);
    if (cov) {
      Parameters d = (G * cov.value());
      Parameters newpoint = x + d;
      auto chol_cov = chol(cov.value(), "lower");
      if (chol_cov.has_value()) {
        return Op(
            Distribution(typename Distribution::base_type(
                             newpoint, cov.value(), inv_cov, chol_cov.value()),
                         landa_));
      } else
        return Op(false,
                  " invalid cholesky decomposition: " + chol_cov.error());

    } else
      return Op(false, "cannot invert Hessian landa matrix " + cov.error());
  }

  template <class mySample>
  auto compute_Likelihood(const Model &m, const Parameters &x,
                          const mySample &current, std::ostream &os) const {
    return m.compute_DLikelihood(x, current, os);
  }

  auto compute_Likelihood_init(const Model &m, const Parameters &x,
                               std::ostream &os) const {
    return m.compute_DLikelihood_init(x, os);
  }

  LevenbergMarquardt(double _landa) : landa_{_landa} {}
  LevenbergMarquardt() = default;

  bool operator<(const LevenbergMarquardt &other) const {
    return landa() < other.landa();
  }

  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "landa", &self_type::landa));
  }

private:
  double landa_;
};

template <class Model> class AdaptiveCovariance {
public:
  typedef typename Model::Parameters Parameters;

  typedef Cs<Model> template_types;
  typedef Normal_Distribution<Parameters> Distribution;
  typedef typename Model::LikelihoodResult LikelihoodResult;

  typedef AdaptiveCovariance self_type;
  constexpr static auto const className =
      my_static_string("AdaptiveCovariance") +
      my_trait<template_types>::className;

  Distribution getDistribution(const Model &, const Parameters &x,
                               const LikelihoodResult &) const {
    return Normal_Distribution<Parameters>(x, cov_);
  }

  LikelihoodResult compute_Likelihood(const Model &m, const Parameters &x) {
    return m.compute_Likelihood(x);
  }

  AdaptiveCovariance(double r, const M_Matrix<double> &c)
      : cov_{c * r}, ratio_(r) {}

  double ratio() const { return ratio_; }

private:
  M_Matrix<double> cov_;
  double ratio_;
};

template <class Model, class aLikelihoodResult, class aParameters,
          class aDistribution>
class mcmc_sample : public aLikelihoodResult {
public:
  typedef aParameters Parameters;
  typedef aDistribution Distribution;
  typedef aLikelihoodResult LikelihoodResult;

  typedef aLikelihoodResult base_type;

  typedef mcmc_sample self_type;
  typedef Cs<Model, aParameters, aDistribution> template_types;
  constexpr static auto const className = my_static_string("mcmc_sample") +
                                          my_trait<base_type>::className +
                                          my_trait<template_types>::className;

  std::string myClass() const { return className.str(); }

  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "Parameters", &self_type::x),
        grammar::field(C<self_type>{}, "likelihood_result",
                       &self_type::likelihood_result),
        grammar::field(C<self_type>{}, "Distribution", &self_type::myg));
  }

  static std::vector<std::string> data_index_titles() {
    std::vector<std::string> out;
    out.push_back("ParameterName");
    out.push_back("ParameterNameT");
    out = concatenate_unique(std::move(out),
                             LikelihoodResult::data_index_titles());
    out = concatenate_unique(std::move(out), Distribution::data_index_titles());
    return out;
  }

  static Data_Index_scheme data_index() {
    Data_Index_scheme out;
    out.push_back("ParameterName", {"i_Parameter"});
    out.push_back("ParameterName", {"i_Parameter", "i_Parameter_T"});
    out.push_back("ParameterName_T", {"i_Parameter", "i_Parameter_T"});
    out.push_back("ParameterValue", {"i_Parameter"});
    out.push_back("ParameterValue_T", {"i_Parameter", "i_Parameter_T"});

    out.push_back("Accepting", {});
    out = concatenate(std::move(out), LikelihoodResult::data_index(),
                      Distribution::data_index());
    return out;
  }

  static auto get_data_index_static() {

    // typedef typename
    // decltype(LikelihoodResult::get_data_index_static())::test error;
    return Concatenate_tuple_static(
        std::make_tuple(
            make_data_static(
                std::tuple<>(),
                std::make_tuple(F_s(s_accepting{}, &self_type::accept))),

            make_data_static(
                std::make_tuple(I_s(
                    i_Parameter{},
                    [](const self_type &self) { return self.x().size(); }, 0)),
                std::make_tuple(
                    F_s(s_ParameterValue{},
                        [](const self_type &self, std::size_t i_par) {
                          return self.x()[i_par];
                        })))),
        Compose_static([](const self_type &self) -> auto const & {
          return self.likelihood_result();
        },
                       LikelihoodResult::get_data_index_static()));
    //     Compose_static([](const self_type& self) {return self.g();},
    //                    Distribution::get_data_index_static()));
  }

  auto getIndexedData(std::set<std::string> &index) {
    Data_Index_point out;

    auto it = index.find("i_Parameter");
    if (it != index.end()) {
      out.push_back("i_Parameter", x_.size());
      index.erase(it);
    }
    out.push_back("i_Parameter_T", x_.size());
    LikelihoodResult &l = *this;
    auto ind = l.getIndexedData(index);
  }

  LikelihoodResult const &likelihood_result() const { return *this; }
  Parameters const &x() const { return x_; }
  Distribution const &myg() const { return g_; }
  Distribution const &g() const & { return g_; }
  Distribution g() && { return std::move(g_); }

  bool accept() const { return accepted_; }
  void push_rejected() { accepted_ = false; }
  void push_accepted() { accepted_ = true; }

  void set_g(Distribution &&new_g) { g_ = std::move(new_g); }
  mcmc_sample(const Parameters &par, LikelihoodResult &&lik,
              Distribution &&dist)
      : LikelihoodResult(std::move(lik)), x_{par}, g_{std::move(dist)},
        accepted_{true} {}

  mcmc_sample() = default;

private:
  Parameters x_;
  Distribution g_;
  bool accepted_;
};

template <class Model, class Distribution_Generator>
using mcmc_sample_t =
    mcmc_sample<Model, typename Distribution_Generator::LikelihoodResult,
                typename Model::Parameters,
                typename Distribution_Generator::Distribution>;

template <class Model, class Distribution_Generator, class mySample,
          class Parameters = typename Model::Parameters>
static myOptional_t<mcmc_sample_t<Model, Distribution_Generator>>
calculate(const Model &M, const Distribution_Generator &G, Parameters &&x,
          const mySample &current, std::ostream &os) {
  typedef myOptional_t<mcmc_sample_t<Model, Distribution_Generator>> Op;

  auto logL = G.compute_Likelihood(M, x, current, os);
  os << "\nlogL.likelihood().H()\n" << logL.value().likelihood().H();
  os << "\nlogL.H()\n" << logL.value().H();
  //   std::cerr<<" gets a logL!!----------="<<logL<<"\n";
  if (!logL.has_value())
    return Op(false, "cannot compute logLikelihood because " + logL.error());
  auto g = G.calculate_Distribution(M, x, logL.value());
  if (!g)
    return Op(false, "cannot compute distribution because " + g.error());
  return Op(mcmc_sample_t<Model, Distribution_Generator>(
      std::move(x), std::move(logL).value(), std::move(g).value()));
}

template <class Model, class Distribution_Generator,
          class Parameters = typename Model::Parameters>
static myOptional_t<mcmc_sample_t<Model, Distribution_Generator>>
calculate_init(const Model &M, const Distribution_Generator &G, Parameters &&x,
               std::ostream &os) {
  typedef myOptional_t<mcmc_sample_t<Model, Distribution_Generator>> Op;

  auto logL = G.compute_Likelihood_init(M, x, os);
  os << "\nlogL.likelihood().H()\n" << logL.value().likelihood().H();
  os << "\nlogL.H()\n" << logL.value().H();
  //   std::cerr<<" gets a logL!!----------="<<logL<<"\n";
  if (!logL.has_value())
    return Op(false, "cannot compute logLikelihood because " + logL.error());
  auto g = G.calculate_Distribution(M, x, logL.value());
  if (!g)
    return Op(false, "cannot compute distribution because " + g.error());
  return Op(mcmc_sample_t<Model, Distribution_Generator>(
      std::move(x), std::move(logL).value(), std::move(g).value()));
}

template <bool Hastings> constexpr static auto str_Hastings() {
  if constexpr (Hastings)
    return my_static_string("_Hastings");
  else
    return my_static_string("");
}

template <class RandomEngine>
static std::vector<RandomEngine> mts(RandomEngine &mt, std::size_t n) {
  std::uniform_int_distribution<typename RandomEngine::result_type> useed;
  std::vector<RandomEngine> out;
  for (std::size_t i = 0; i < n; ++i)
    out.emplace_back(useed(mt));
  return out;
}

template <bool Hastings = true> class Metropolis_H {
public:
  typedef Metropolis_H self_type;
  constexpr static auto const className =
      my_static_string("Metropolis") + str_Hastings<Hastings>();

  static auto get_constructor_fields() {
    return std::make_tuple(grammar::field(C<self_type>{},
                                          "max_number_of_trials",
                                          &self_type::max_number_of_trials));
  }

  template <class mcmcsample, class stream>
  static double Acceptance(const mcmcsample &current,
                           const mcmcsample &candidate, stream &s) {
    double logL_current = double(current.logL());
    double logL_candidate = double(candidate.logL());
    if constexpr (Hastings) {
      double logP_forward = current.g().logP(candidate.x());
      double logP_backward = candidate.g().logP(current.x());
      double logA =
          (logL_candidate + logP_backward) - (logL_current + logP_forward);
      double A = std::min(1.0, std::exp(logA));
      s << "\n "
           "logA=(logL_candidate+logP_backward)-(logL_current+logP_forward)\n "
           "logA="
        << logA << "\n";
      s << "\nlogL_candidate=" << logL_candidate << "\n";
      s << "\nlogL_current=" << logL_current << "\n";

      s << "\nlogP_foward=" << logP_forward << "\n";
      s << "\nlogP_backward=" << logP_backward << "\n";
      s << "\ncandidate.beta()=" << candidate.beta() << "\n";
      s << "\ncurrent.beta()=" << current.beta() << "\n";

      s << "\ncandidate.g().landa()=" << candidate.g().landa() << "\n";
      s << "\ncurrent.g().landa()=" << current.g().landa() << "\n";

      s << "\n logDetCov_forward =" << current.g().logDetCov() << "\n";
      s << "\n logDetCov_backward =" << candidate.g().logDetCov() << "\n";

      auto dcandidate = candidate.x() - candidate.g().mean();
      auto dforward = candidate.x() - current.g().mean();
      auto dback = current.x() - candidate.g().mean();
      auto dcurrent = current.x() - current.g().mean();

      s << "\ndcandidate\n" << dcandidate;
      s << "\ndforward\n" << dforward;
      s << "\n dback\n" << dback;
      s << "\n dcurrent\n" << dcurrent;

      s << "\nxTSigmaX(dforward, current.g().CovInv())=\n"
        << xTSigmaX(current.g().CovInv(), dforward) << "\n";

      s << "\nxTSigmaX(dforward, current.g().CovInv())=\n"
        << xTSigmaX(dforward, current.g().CovInv()) << "\n";
      s << "\nxTSigmaX(dcurrent, current.g().CovInv())=\n"
        << xTSigmaX(dcurrent, current.g().CovInv()) << "\n";
      s << "\n INV xTSigmaX(dcandidate, current.g().CovInv())=\n"
        << xTSigmaX(dcandidate, current.g().CovInv()) << "\n";
      s << "\n INV xTSigmaX(dback, current.g().CovInv())=\n"
        << xTSigmaX(dback, current.g().CovInv()) << "\n";

      s << "\nxTSigmaX(dforward, candidate.g().CovInv())=\n"
        << xTSigmaX(dforward, candidate.g().CovInv()) << "\n";
      s << "\nxTSigmaX(dcurrent, candidate.g().CovInv())=\n"
        << xTSigmaX(dcurrent, candidate.g().CovInv()) << "\n";
      s << "\n INV xTSigmaX(dcandidate, candidate.g().CovInv())=\n"
        << xTSigmaX(dcandidate, candidate.g().CovInv()) << "\n";
      s << "\n INV xTSigmaX(dback, candidate.g().CovInv())=\n"
        << xTSigmaX(dback, candidate.g().CovInv()) << "\n";

      s << "\nquadraticForm_B_A_BT( current.g().CovInv(),dforward)=\n"
        << quadraticForm_B_A_BT(current.g().CovInv(), dforward) << "\n";

      s << "\nquadraticForm_B_A_BT(candidate.g().CovInv(),dback)=\n"
        << quadraticForm_B_A_BT(candidate.g().CovInv(), dback) << "\n";

      s << "\n current.x() \n" << current.x();
      s << "\n current.g().mean() \n" << current.g().mean();

      s << "\n candidate.x() \n" << candidate.x();
      s << "\n candidate.g().mean() \n" << candidate.g().mean();

      s << "\n current.x()-candidate.x() \n" << current.x() - candidate.x();
      s << "\n ds=-current.g().mean()+candidate.x() \n"
        << -current.g().mean() + candidate.x();
      s << "\n current.g().mean()-candidate.g().mean() \n"
        << current.g().mean() - candidate.g().mean();

      s << "\ncurrent.g().Chol()\n" << current.g().Chol();

      s << "\ncandidate.g().Chol()\n" << candidate.g().Chol();

      s << "\ncurrent.g().Chol()-candidate.g().Chol()\n"
        << current.g().Chol() - candidate.g().Chol();

      s << "\ncurrent.g().CovInv()\n" << current.g().CovInv();
      s << "\ncandidate.g().CovInv()\n" << candidate.g().CovInv();

      // s << "\n current \n" << current;
      // s << "\n candidate \n" << candidate;
      s << "-------------------------------------------------------------------"
           "------\n";

      return A;
    } else {
      double logA = (logL_candidate) - (logL_current);
      double A = std::min(1.0, std::exp(logA));
      return A;
    }
  }

  template <class Model, class Distribution_Generator,
            class Parameters = typename Model::Parameters>
  myOptional_t<mcmc_sample_t<Model, Distribution_Generator>>
  start(const Model &M, std::mt19937_64 &mt, const Distribution_Generator &G,
        std::ostream &os) const {
    std::size_t n = 0;
    typedef myOptional_t<mcmc_sample_t<Model, Distribution_Generator>> Op;

    //  std::cerr<<"\nmax number of trials="<<max_number_of_trials();
    std::string error;
    while (n < max_number_of_trials()) {

      Parameters x = M.sample(mt);
      auto out = calculate_init(M, G, x, os);
      if (out.has_value()) {
        return out;
      } else {
        error += out.error();
        ++n;
      }
    }
    return Op(false, "fails to start after n=" + ToString(n) +
                         " trials. Error: " + error);
  }

  template <class Model, class Distribution_Generator, class stream,
            class mySample = mcmc_sample_t<Model, Distribution_Generator>>
  Op_void next(const Model &M, std::mt19937_64 &mt,
               const Distribution_Generator &G, mySample &current,
               stream &s) const {
    auto candidate_point = current.g().sample(mt);
    // current.g().autoTest(mt,900);
    myOptional_t<mySample> candidate =
        calculate(M, G, candidate_point, current, s);
    if (!candidate) {
      current.push_rejected();
    } else {
      // auto gg=G.calculate_Distribution(M,candidate.value().x(),
      // candidate.value());
      //            s<<"gg.mean()"<<gg.value().mean();
      //            s<<"\ncandidate.g().mean()\n"<<candidate.value().g().mean();
      //            s<<"gg.Cov()"<<gg.value();
      //            s<<"\ncandidate.g()\n"<<candidate.value().g();
      //            s<<"\ncandidate.g().Cov()-
      //            gg.value().Cov()\n"<<candidate.value().g().Cov()-
      //            gg.value().Cov();
      double A = Acceptance(current, candidate.value(), s);
      std::uniform_real_distribution<double> u(0, 1);
      double r = u(mt);
      bool accept = r < A;
      if (accept) {
        current = std::move(candidate.value());
        current.push_accepted();
      } else
        current.push_rejected();
    }
    return Op_void(true, "");
  }

  template <class Model, class mySample>
  auto compute_PartialDLikelihood(const Model &M, const mySample &current,
                                  std::ostream &os) const {
    return M.compute_PartialDLikelihood(current.x(), os);
  }

  std::size_t max_number_of_trials() const { return ntrials_; }
  Metropolis_H<Hastings>(std::size_t max_trials) : ntrials_{max_trials} {}
  Metropolis_H() = default;

private:
  std::size_t ntrials_;
};

typedef Metropolis_H<true> Metropolis_Hastings;

typedef Metropolis_H<false> Metropolis;

template <bool Hastings = true> class Adaptive_Metropolis_H {
public:
  typedef Adaptive_Metropolis_H self_type;
  constexpr static auto const className =
      my_static_string("Adaptive_Metropolis") + str_Hastings<Hastings>();

  static auto get_constructor_fields() {
    return std::make_tuple(grammar::field(C<self_type>{},
                                          "max_number_of_trials",
                                          &self_type::max_number_of_trials));
  }
  Adaptive_Metropolis_H() = default;
  auto max_number_of_trials() const { return mh_.max_number_of_trials(); }

  template <class Model, class Adapative_Distribution_Generator,
            class mySample = myOptional_t<
                mcmc_sample_t<Model, Adapative_Distribution_Generator>>>
  mySample start(const Model &M, std::mt19937_64 &mt,
                 Adapative_Distribution_Generator &adaptive,
                 std::ostream &os) const {
    auto G = adaptive.sample(mt);
    std::cerr << " start " << G;
    auto out = get_Metropolis().start(M, mt, G, os);
    return out;
  }

  template <
      class Model, class Adapative_Distribution_Generator, class stream,
      class mySample = mcmc_sample_t<Model, Adapative_Distribution_Generator>>
  myOptional_t<void> next(const Model &M, std::mt19937_64 &mt,
                          Adapative_Distribution_Generator &adaptive,
                          mySample &current, stream &s) const {
    auto G = adaptive.sample(mt);
    auto gg =
        G.calculate_Distribution(M, current.x(), current.likelihood_result());
    if (!gg) {
      return myOptional_t<void>(false,
                                "cannot calculate distribution " + gg.error());
    } else {
      //    std::cerr<<"\nold current.g()\n"<<current.g();

      //  auto gold=gg.value();
      current.set_g(std::move(gg).value());
      //   s<<"\n new
      //   current.g().mean()-ggold.mean()\n"<<current.g().mean()-gold.mean();
      //   s<<"\n new
      //   current.g().Cov()-ggold.COv()\n"<<current.g().Cov()-gold.Cov();

      get_Metropolis().next(M, mt, G, current, s);
      adaptive.push_outcome(G, current, current.accept());
      return myOptional_t<void>(true, "");
    }
  }

  Metropolis_H<Hastings> const &get_Metropolis() const { return mh_; }

  Adaptive_Metropolis_H(Metropolis_H<Hastings> mh) : mh_{std::move(mh)} {}
  Adaptive_Metropolis_H(std::size_t max_num_trials) : mh_{max_num_trials} {}

private:
  Metropolis_H<Hastings> mh_;
};

typedef Adaptive_Metropolis_H<true> Adaptive_Metropolis_Hastring;
typedef Adaptive_Metropolis_H<false> Adaptive_Metropolis;

template <class Model, class AdaptiveMover> class emcee_sample {
public:
  typedef emcee_sample self_type;
  typedef Cs<Model, AdaptiveMover> template_types;
  constexpr static auto const className =
      my_static_string("emcee_sample") + my_trait<template_types>::className;

  typedef mcmc_sample_t<Model, AdaptiveMover> mcmc_s;
  mcmc_s &Walker(std::size_t i) { return walkers_[i]; }
  mcmc_s const &Walker(std::size_t i) const { return walkers_[i]; }
  std::size_t numWalkers() const { return walkers_.size(); }
  bool accept(std::size_t i) const { return Walker(i).accept(); }

  template <class Random_Engine>
  static emcee_sample evaluate(const Model &m,
                               std::vector<Random_Engine> & /*mt*/,
                               AdaptiveMover &mov, std::size_t nwalkers) {
    auto w = getWalkers(m, mov, nwalkers);
  }

  auto likelihood() const {
    auto lik = Walker(0).likelihood();
    for (std::size_t i = 1; i < numWalkers(); ++i)
      lik += Walker(i).likelihood();
    lik /= numWalkers();
    return lik;
  }

  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "walkers", &self_type::walkers));
  }

  auto &walkers() const { return walkers_; }
  // auto& acceptance()const {return  accept_;}
  emcee_sample() = default;
  emcee_sample(std::vector<mcmc_s> &&walkers) : walkers_{std::move(walkers)} {}
  emcee_sample(const std::vector<mcmc_s> &walkers, const std::vector<bool> &)
      : walkers_{walkers} {}

  static Data_Index_scheme data_index() {
    Data_Index_scheme out;
    out.push_back("aceptance", {"i_walker"});
    auto vmc = Insert_Index(mcmc_s::data_index(), {"i_walker"});
    out = concatenate(std::move(out), std::move(vmc));
    return out;
  }

  static auto get_data_index_static() {
    return Compose_static(
        std::make_tuple(I_s(
            i_walker{},
            [](const self_type &self) { return self.walkers().size(); }, 0)),
        [](const self_type &self, std::size_t i_w) -> auto const & {
          return self.walkers()[i_w];
        },
        mcmc_s::get_data_index_static());
  }

private:
  std::vector<mcmc_s> walkers_;

  template <class Random_Generator>
  static myOptional_t<mcmc_s>
  getWalker(const Model &model, Random_Generator &mt,
            const AdaptiveMover &adaptive, std::size_t n_tries) {
    typedef myOptional_t<mcmc_s> Op;
    std::size_t n = 0;
    std::string error;
    while (true) {
      auto x = model.sample(mt);
      auto logL = model.compute_Likelihood(x);
      if (logL) {
        auto d = adaptive.sample(mt);
        return Op(mcmc_s(std::move(x), std::move(logL).value(), std::move(d)));
      } else if (n < n_tries) {
        ++n;
        error += logL.error();
      } else
        return Op(false, error);
    }
  }

public:
  template <class Random_Generator>
  static myOptional_t<std::vector<mcmc_s>>
  getWalkers(const Model &model, Random_Generator &mt,
             const AdaptiveMover &adaptive, std::size_t n,
             std::size_t n_tries) {
    typedef myOptional_t<std::vector<mcmc_s>> Op;
    std::vector<mcmc_s> out(n);
    std::vector<Op_void> op(n, Op_void(true, ""));
    //#pragma omp parallel for
    for (std::size_t i = 0; i < n; ++i) {
      auto s = getWalker(model, mt, adaptive, n_tries);
      if (s)
        out[i] = std::move(s).value();
      else
        op[i] = Op_void(false, s.error());
    }
    auto res = consolidate(std::move(op));
    if (res)
      return Op(out);
    else
      return Op(false, res.error());
  }
};

class Ensemble_Metropolis_Hastings {
private:
  std::size_t numWalkers_;
  std::size_t numTrials_;

public:
  typedef Ensemble_Metropolis_Hastings self_type;
  constexpr static auto const className =
      my_static_string("Ensemble_Metropolis_Hastings");

  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "numWalkers", &self_type::numWalkers),
        grammar::field(C<self_type>{}, "numTrials", &self_type::numTrials));
  }

  std::size_t numWalkers() const { return numWalkers_; }
  std::size_t numTrials() const { return numTrials_; }

  Ensemble_Metropolis_Hastings(std::size_t numberWalkers, std::size_t num_tries)
      : numWalkers_{numberWalkers}, numTrials_{num_tries} {}
  Ensemble_Metropolis_Hastings() = default;
  template <class Model> class stretch_move {
  public:
    typedef typename Model::LikelihoodResult LikelihoodResult;
    typedef stretch_move_Distribution Distribution;
    typedef typename Model::Parameters Parameters;

    static myOptional_t<LikelihoodResult> compute_Likelihood(const Model &m,
                                                             Parameters &x) {
      return m.compute_Likelihood(x);
    }

    constexpr static auto const className =
        my_static_string("stretch_move") + my_trait<Model>::className;

    typedef stretch_move self_type;
    static auto get_constructor_fields() { return std::tuple<>(); }

    static double Acceptance(double Z, const logLikelihood &candidateLogLik,
                             const logLikelihood &currentLogLik,
                             std::size_t numParam) {
      double logA = candidateLogLik.logL() - currentLogLik.logL() +
                    log(Z) * (numParam - 1);
      return std::min(0.0, logA);
    }

    template <class Adaptive_Strecth_Move>
    static void move(const Model &model, std::mt19937_64 &mt,
                     const stretch_move_Distribution &d,
                     const emcee_sample<Model, Adaptive_Strecth_Move> &s,
                     const std::vector<std::size_t> &index,
                     mcmc_sample_t<Model, Adaptive_Strecth_Move> &current) {
      std::uniform_int_distribution<std::size_t> u(1, index.size());
      std::size_t i = index[u(mt) - 1];
      double z = d.sample(mt);
      Parameters candidate =
          (current.x() - (s.Walker(i).x())) * z + s.Walker(i).x();
      auto lik = compute_Likelihood(model, candidate);
      if (!lik) {
        current.push_rejected();

      } else {
        auto logA = Acceptance(z, lik.value(), current, candidate.size());
        double A = std::exp(logA);
        std::uniform_real_distribution<double> real(0, 1);
        double r = real(mt);
        bool accept = r < A;
        if (accept) {
          current = mcmc_sample_t<Model, Adaptive_Strecth_Move>(
              std::move(candidate), std::move(lik).value(),
              std::move(current).g());
          current.push_accepted();
        } else
          current.push_rejected();
      }
    }
  };

  template <class Model> class adaptive_stretch_mover {
  public:
    typedef stretch_move<Model> Mover;
    typedef typename Model::LikelihoodResult LikelihoodResult;
    typedef typename Mover::Distribution Distribution;
    typedef typename Mover::Parameters Parameters;

    constexpr static auto const className =
        my_static_string("adaptive_stretch_mover_") +
        my_trait<Model>::className;

    typedef adaptive_stretch_mover self_type;
    static auto get_constructor_fields() {
      return std::make_tuple(
          grammar::field(C<self_type>{}, "target", &self_type::target),
          grammar::field(C<self_type>{}, "alfa_map", &self_type::alfa_map),
          grammar::field(C<self_type>{}, "alfa", &self_type::alfa));
    }
    stretch_move_Distribution sample(std::mt19937_64 &mt) const {
      return stretch_move_Distribution(alfa_.sample(mt));
    }

    virtual ~adaptive_stretch_mover() {}
    virtual std::ostream &put(std::ostream &os) const { return os << *this; }

    void push_outcome(const stretch_move_Distribution &d, bool accept) {
      if (accept)
        alfa_map_[d.alpha()].push_accept();
      else
        alfa_map_[d.alpha()].push_reject();
    }

    void actualize() { alfa_ = alfa_map_.Distribute_on_p(target_); }
    adaptive_stretch_mover(double target, const std::vector<double> &alfas)
        : target_{target}, alfa_map_{Beta_map<double>::UniformPrior(alfas)},
          alfa_{} {
      actualize();
    }
    adaptive_stretch_mover() = default;

    adaptive_stretch_mover(double target, const Beta_map<double> &alfa_map,
                           const Probability_map<double> &alfa)
        : target_{target}, alfa_map_{alfa_map}, alfa_{alfa} {}

    double target() const { return target_; }
    Beta_map<double> const &alfa_map() const { return alfa_map_; }
    Probability_map<double> const &alfa() const { return alfa_; }

    adaptive_stretch_mover(const adaptive_stretch_mover &other)
        : target_{other.target_}, alfa_map_{other.alfa_map_},
          alfa_{other.alfa_} {}
    adaptive_stretch_mover(adaptive_stretch_mover &&other)
        : target_{other.target_}, alfa_map_{std::move(other.alfa_map_)},
          alfa_{std::move(other.alfa_)} {}

    adaptive_stretch_mover &operator=(const adaptive_stretch_mover &other) {
      adaptive_stretch_mover tmp(other);
      *this = std::move(tmp);
      return *this;
    }
    adaptive_stretch_mover &operator=(adaptive_stretch_mover &&other) {
      target_ = other.target_;
      alfa_map_ = std::move(other.alfa_map_);
      alfa_ = std::move(other.alfa_);
    }

    static Data_Index_scheme data_index() {
      Data_Index_scheme out;
      out.push_back("target", {"const"});
      out = concatenate(
          std::move(out),
          Beta_map<double>::data_index("stretch_alfa", "stretch_alfa"),
          Probability_map<double>::data_index("stretch_alfa", "stretch_alfa"));

      return out;
    }

    static auto get_data_index_static() {
      return Concatenate_tuple_static(
          std::make_tuple(make_data_static(
              std::make_tuple(I_s(i_stretch_alpha(),
                                  [](const self_type &self) {
                                    return self.alfa_map().x().size();
                                  },
                                  0)),
              std::make_tuple(F_s(s_stretch_alpha{},
                                  [](const self_type &s, std::size_t i_a) {
                                    return s.alfa().x()[i_a];
                                  }),
                              F_s(CT_s<s_stretch_alpha, s_use, s_probability>{},
                                  [](const self_type &s, std::size_t i_a) {
                                    return s.alfa().p().at(s.alfa().x()[i_a]);
                                  })))),
          Compose_static(
              std::make_tuple(I_s(i_stretch_alpha(),
                                  [](const self_type &self) {
                                    return self.alfa_map().x().size();
                                  },
                                  0)),
              F_s(CT_s<s_stretch_alpha, s_obs, s_probability>{},
                  [](const self_type &self, std::size_t i_a) -> auto const & {
                    return self.alfa_map().get_Map().at(
                        self.alfa_map().x()[i_a]);
                  }),
              Beta_Distribution::get_data_index_static()));
    }

  private:
    double target_;
    Beta_map<double> alfa_map_;
    Probability_map<double> alfa_;
  };

  template <class Random_Engine>
  static std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
  random_split(Random_Engine &mt, std::size_t n) {
    std::vector<std::size_t> index(n);
    for (std::size_t i = 0; i < n; ++i)
      index[i] = i;

    std::shuffle(index.begin(), index.end(), mt);
    std::size_t half = n / 2;

    std::vector<std::size_t> one(index.begin(), index.begin() + half);
    std::vector<std::size_t> two(index.begin() + half, index.end());
    return {one, two};
  }

  template <class Model, class Randome_Engine, class Adaptive_Mover,
            class Mover = typename Adaptive_Mover::Mover>
  myOptional_t<emcee_sample<Model, Adaptive_Mover>>
  start(const Model &model, Randome_Engine &mt,
        Adaptive_Mover &adaptive) const {
    typedef myOptional_t<emcee_sample<Model, Adaptive_Mover>> Op;
    auto w = emcee_sample<Model, Adaptive_Mover>::getWalkers(
        model, mt, adaptive, numWalkers(), numTrials());
    if (w)
      return Op(emcee_sample<Model, Adaptive_Mover>(std::move(w).value()));
    else {
      return Op(false, w.error());
    }
  }

  template <class Model, class Random_Engine,
            template <class> class Adaptive_Mover,
            class Mover = typename Adaptive_Mover<Model>::Mover>
  bool next_para(const Model &model, std::vector<Random_Engine> &mt,
                 adaptive_stretch_mover<Model> &adaptive,
                 emcee_sample<Model, stretch_move<Model>> &current) {

    //   typedef typename Adaptive_Mover<Model>::Mover stretch_move;
    auto tu = random_split(mt[0], current.numWalkers());
    auto one = std::move(std::get<0>(tu));
    auto two = std::move(std::get<1>(tu));
    //    auto [one, two] =tu;
#pragma omp parallel for
    for (std::size_t ii = 0; ii < one.size(); ++ii) {
      auto i = one[ii];
      current.Walker(i).set_g(adaptive.sample(mt[ii]));
      Mover::move(model, current.mt(i), current.Walker(i).g(), current, one,
                  current.Walker(i));
    }
    for (std::size_t ii = 0; ii < one.size(); ++ii) {
      auto i = one[ii];
      adaptive.push_outcome(current.Walker(i).g(), current.accept(i));
    }

#pragma omp parallel for
    for (std::size_t ii = 0; ii < two.size(); ++ii) {
      auto i = two[ii];
      current.Walker(i).g() = adaptive.sample(mt);
      stretch_move<Model>::move(model, current.mt(i), current.Walker(i).g(),
                                current, two, current.Walker(i));
    }
    for (std::size_t ii = 0; ii < two.size(); ++ii) {
      auto i = two[ii];
      adaptive.push_outcome(current.Walker(i).g(), current.accept(i));
    }
    adaptive.actualize();
    return current;
  }

  template <class Model, class Random_Engine, class Adaptive_Mover,
            class Sample>
  Op_void next(const Model &model, Random_Engine &mt, Adaptive_Mover &adaptive,
               Sample &current) const {

    typedef typename Adaptive_Mover::Mover Mover;
    auto tu = random_split(mt, current.numWalkers());
    auto one = std::move(std::get<0>(tu));
    auto two = std::move(std::get<1>(tu));
    //    auto [one, two] =tu;
    //#pragma omp parallel for
    for (std::size_t ii = 0; ii < one.size(); ++ii) {
      auto i = one[ii];
      current.Walker(i).set_g(adaptive.sample(mt));
      Mover::move(model, mt, current.Walker(i).g(), current, one,
                  current.Walker(i));
    }
    for (std::size_t ii = 0; ii < one.size(); ++ii) {
      auto i = one[ii];
      adaptive.push_outcome(current.Walker(i).g(), current.accept(i));
    }

    //#pragma omp parallel for
    for (std::size_t ii = 0; ii < two.size(); ++ii) {
      auto i = two[ii];
      current.Walker(i).set_g(adaptive.sample(mt));
      stretch_move<Model>::move(model, mt, current.Walker(i).g(), current, two,
                                current.Walker(i));
    }
    for (std::size_t ii = 0; ii < two.size(); ++ii) {
      auto i = two[ii];
      adaptive.push_outcome(current.Walker(i).g(), current.accept(i));
    }
    adaptive.actualize();
    return Op_void(true, "");
  }

  template <class Model, class mySample>
  auto compute_PartialDLikelihood(const Model &M, const mySample &current,
                                  std::ostream &os) const {
    std::vector<
        std::decay_t<decltype(M.compute_PartialDLikelihood(
                                   current.Walker(0).x(), current.Walker(0), os)
                                  .value())>>
        out(current.numWalkers());

    for (std::size_t i = 0; i < current.numWalkers(); ++i) {
      out[i] = M.compute_PartialDLikelihood(current.Walker(i).x(),
                                            current.Walker(i), os)
                   .value();
    }
    return out;
  }

  static std::vector<std::mt19937_64> mts(std::mt19937_64 &mt, std::size_t n) {
    std::uniform_int_distribution<typename std::mt19937_64::result_type> useed;
    std::vector<std::mt19937_64> out(n);
    for (std::size_t i = 0; i < n; ++i)
      out[i].seed(useed(mt));
    return out;
  }
};

template <class Parameterized_Distribution /*=Landa*/,
          class Gain = typename Parameterized_Distribution::ExpectedVelocity,
          class Likelihood = typename Parameterized_Distribution::myAcceptProb>
class Adaptive_Parameterized_Distribution_Generator {
public:
  typedef
      typename Parameterized_Distribution::LikelihoodResult LikelihoodResult;
  typedef typename Parameterized_Distribution::Distribution Distribution;

  typedef Adaptive_Parameterized_Distribution_Generator self_type;
  typedef Cs<Parameterized_Distribution, Gain, Likelihood> template_types;
  constexpr static auto const className =
      my_static_string("Adaptive_Parameterized_Distribution_Generator") +
      my_trait<template_types>::className;

  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "likelihoodFunction",
                       &self_type::getLikelihood),
        grammar::field(C<self_type>{}, "gainFunction", &self_type::getGain),
        grammar::field(C<self_type>{}, "gainMoment", &self_type::getGainMoment),
        grammar::field(C<self_type>{}, "probabilityMap",
                       &self_type::get_ProbabilityMap),
        grammar::field(C<self_type>{}, "likelihoodProbabilityMap",
                       &self_type::get_Parameters_Map),
        grammar::field(C<self_type>{}, "logEvidence",
                       &self_type::getLogEvidence)

    );
  }

  Parameterized_Distribution sample(std::mt19937_64 &mt) const {
    return p_.sample(mt);
  }

  template <class mcmc>
  void push_outcome(const Parameterized_Distribution &landa, const mcmc &,
                    bool accept) {
    if (accept)
      push_acceptance(landa);
    else
      push_rejection(landa);
    actualize();
  }

  void push_acceptance(const Parameterized_Distribution &landa) {
    parDist_ = logBayes_rule(Log_of(lik_), landa, parDist_);
  }

  void push_rejection(const Parameterized_Distribution &landa) {
    parDist_ = logBayes_rule(Log_of(Complement_prob(lik_)), landa, parDist_);
  }

  void actualize(double nmax) {
    actualize();
    parDist_.reduce(nmax);
  }

  static std::pair<Probability_map<Parameterized_Distribution>, double>
  Distribute_on_gain(
      const Gain &g, const Likelihood &lik,
      const logLikelihood_map<typename Likelihood::Parameters> &par,
      const Probability_map<Parameterized_Distribution> &landas,
      double moment) {
    auto out = landas.p();
    auto p_par = par.p();
    for (auto &e : out) {
      double meangain = Expectance(g, lik, p_par, e.first);
      e.second = std::pow(meangain, moment);
    }
    auto o = Probability_map<Parameterized_Distribution>::normalize(
        out, landas.nsamples());
    double sum = 0;
    for (auto &e : o.first.p()) {
      sum += e.second * std::pow(out[e.first], 1.0 / moment);
    }
    return {o.first, sum};
  }

  void actualize() {
    auto pold = p_;
    auto o = Distribute_on_gain(g_, lik_, parDist_, p_, gainMoment_);
    p_ = o.first;
  }

  /*Adaptive_Parameterized_Distribution_Generator(const
std::map<Parameterized_Distribution,double>& prior_landa, const
std::map<typename Likelihood::Parameters, double>& prior_par, double
nsamples,std::size_t gainMoment): gainMoment_(gainMoment),
p_{prior_landa,nsamples},
parDist_(prior_par,nsamples){}
*/

  template <template <typename...> class V>
  Adaptive_Parameterized_Distribution_Generator(
      const V<Parameterized_Distribution> &landa,
      const std::vector<std::vector<double>> &par, double gainMoment)
      : lik_{}, g_{}, gainMoment_(gainMoment), p_{landa},
        parDist_{Likelihood::uniform_parameter_prior(par, 0.0), 0}, logEvidence{
                                                                        0} {}

  Adaptive_Parameterized_Distribution_Generator() {}

  Adaptive_Parameterized_Distribution_Generator(
      Likelihood lik, Gain g, std::size_t gainMoment,
      const Probability_map<Parameterized_Distribution> &pmap,
      const logLikelihood_map<typename Likelihood::Parameters> &pardist,
      double logEvide)
      : lik_{lik}, g_{g}, gainMoment_{gainMoment}, p_{pmap}, parDist_{pardist},
        logEvidence{logEvide} {}

  Likelihood const &getLikelihood() const { return lik_; }

  Gain const &getGain() const { return g_; }

  std::size_t getGainMoment() const { return gainMoment_; }
  Probability_map<Parameterized_Distribution> get_ProbabilityMap() const {
    return p_;
  }
  logLikelihood_map<typename Likelihood::Parameters> const &
  get_Parameters_Map() const {
    return parDist_;
  }
  double getLogEvidence() const { return logEvidence; }

  static Data_Index_scheme data_index() {
    Data_Index_scheme out =
        concatenate(Likelihood::data_index(), Gain::data_index());
    out.push_back("gainMoment", {"const"});
    out = concatenate(
        std::move(out),
        Probability_map<Parameterized_Distribution>::data_index(
            Parameterized_Distribution::index(),
            Parameterized_Distribution::index()),
        logLikelihood_map<typename Likelihood::Parameters>::data_index(
            "likelihood"));

    out.push_back("logEvidence", {});
    return out;
  }

  static auto get_data_index_static() {
    return Concatenate_tuple_static(
        std::make_tuple(
            make_data_static(
                std::tuple<>(),
                std::make_tuple(F_s(s_logEvidence{}, &self_type::getLogEvidence)

                                    )),
            make_data_static(
                std::make_tuple(
                    I_s(typename Parameterized_Distribution::index_type{},
                        [](const self_type &self) {
                          return self.get_ProbabilityMap().x().size();
                        },
                        0)),
                std::make_tuple(
                    F_s(typename Parameterized_Distribution::label_type{},
                        [](const self_type &self, std::size_t i_par) {
                          return self.get_ProbabilityMap().x()[i_par]();
                        }),
                    F_s(CT_s<typename Parameterized_Distribution::label_type,
                             s_use, s_probability>{},
                        [](const self_type &self, std::size_t i_par) {
                          return self.get_ProbabilityMap().p().at(
                              self.get_ProbabilityMap().x()[i_par]());
                        })))),
        Compose_static(F_s(s_obs{}, &self_type::get_Parameters_Map),
                       logLikelihood_map<typename Likelihood::Parameters>::
                           get_data_index_static(
                               typename Likelihood::Parameters::index_type{})));
  }

private:
  Likelihood lik_;
  Gain g_;
  std::size_t gainMoment_;
  Probability_map<Parameterized_Distribution> p_;
  logLikelihood_map<typename Likelihood::Parameters> parDist_;
  double logEvidence;
};

} // namespace evidence

template <class Adaptive_mover>
struct myData_Index<std::vector<Adaptive_mover>> {
  typedef std::vector<Adaptive_mover> self_type;
  static Data_Index_scheme data_index() {
    Data_Index_scheme out = Adaptive_mover::data_index();
    out.insert_index("i_beta");
    return out;
  }

  static auto get_data_index_static() {
    return Compose_static(
        std::make_tuple(I_s(evidence::i_beta{}, &self_type::size, 0)),
        [](const self_type &self, std::size_t i_beta) -> auto const & {
          return self[i_beta];
        },
        Adaptive_mover::get_data_index_static());
  }
};

namespace evidence {

template <class Model>
class AdaptiveCovarianceGenerator

{
  typedef AdaptiveCovarianceGenerator self_type;
  typedef Cs<Model> template_types;
  constexpr static auto const className =
      my_static_string("AdaptiveCovarianceGenerator") +
      my_trait<template_types>::className;

  typedef typename Model::Parameters Parameters;
  AdaptiveCovariance<Model> sample(std::mt19937_64 &mt) {
    double r = r_.sample(mt);
    return AdaptiveCovariance<Model>(r, cov_);
  }

  void actualize() {
    auto m = Sx_ * (1.0 / nsamples_);
    cov_ = Sxx_ * (1.0 / nsamples_) - TranspMult(m, m);
    r_ = rmap_.Distribute_on_p(target_);
  }

  void reset(std::size_t nmax) {
    double f = 1.0 * nmax / nsamples_;
    if (f < 1) {
      nsamples_ = nmax;
      Sx_ *= f;
      Sxx_ *= f;
      rmap_.reduce(nmax);
      r_ = rmap_.Distribute_on_p(target_);
    }
  }
  void push_outcome(const AdaptiveCovariance<Model> &cov, const Parameters &p,
                    bool accept) {
    Sxx_ += TranspMult(p, p);
    Sx_ += p;
    ++nsamples_;
    if (accept)
      rmap_[cov.ratio()].push_accept();
    else
      rmap_[cov.ratio()].push_reject();
    if (nsamples_ % nkip_ == 0) {
      actualize();
    }
  }

private:
  double target_;
  std::size_t nkip_;
  Beta_map<double> rmap_;
  std::size_t nsamples_;
  M_Matrix<double> Sxx_;
  M_Matrix<double> Sx_;
  M_Matrix<double> cov_;
  Probability_map<double> r_;
};

template <class DistGen /*=Landa*/,
          class EV = typename DistGen::ExpectedVelocity>
class Adaptive_Probability_Distribution_Generator {
public:
  typedef typename DistGen::LikelihoodResult LikelihoodResult;
  typedef typename DistGen::Distribution Distribution;

  typedef Adaptive_Probability_Distribution_Generator self_type;
  typedef Cs<EV, DistGen> template_types;
  constexpr static auto const className =
      my_static_string("Adaptive_Probability_Distribution_Generator") +
      my_trait<template_types>::className;

  DistGen sample(std::mt19937_64 &mt) const { return p_.sample(mt); }

  template <class mcmc>
  void push_outcome(const DistGen &landa, const mcmc &, bool accept) {
    if (accept)
      push_acceptance(landa);
    else
      push_rejection(landa);
    actualize();
  }

  void push_acceptance(const DistGen &landa) {
    landaDist_[landa].push_accept();
  }

  void push_rejection(const DistGen &landa) { landaDist_[landa].push_reject(); }

  void actualize(double nmax) {
    actualize();
    landaDist_.reduce(nmax);
  }

  void actualize() {
    if (target_prob() > 0)
      p_ = landaDist_.Distribute_on_p(target_prob());
    else
      p_ = landaDist_.Distribute_on_f(f_);
  }

  Adaptive_Probability_Distribution_Generator(
      const EV &f, const std::map<DistGen, double> &prior_landa,
      double target_prob)
      : f_(f), p_{prior_landa},
        landaDist_{Beta_map<DistGen>::UnInformativePrior(prior_landa)},
        target_prob_{target_prob} {}

  Adaptive_Probability_Distribution_Generator(
      const Probability_map<DistGen> &prior_landa,
      const Beta_map<DistGen> &landa_history, double target_prob)
      : f_(), p_{prior_landa}, landaDist_{landa_history}, target_prob_{
                                                              target_prob} {}

  template <template <typename...> class V, typename... As>
  Adaptive_Probability_Distribution_Generator(const EV &f,
                                              const V<DistGen, As...> &landa,
                                              double target_prob)
      : Adaptive_Probability_Distribution_Generator(f, uniform_prior(landa),
                                                    target_prob) {}

  Adaptive_Probability_Distribution_Generator(
      const std::map<DistGen, double> &prior_landa, double target_prob)
      : f_(), p_{prior_landa}, landaDist_{Beta_map<DistGen>::UnInformativePrior(
                                   prior_landa)},
        target_prob_{target_prob} {}

  template <template <typename...> class V, typename... As>
  Adaptive_Probability_Distribution_Generator(const V<DistGen, As...> &landa,
                                              double target_prob)
      : Adaptive_Probability_Distribution_Generator(uniform_prior(landa),
                                                    target_prob) {}

  Adaptive_Probability_Distribution_Generator() {}

  auto &Expected_Velocity() const { return f_; }

  static auto get_data_index_static() {
    return Concatenate_tuple_static(
        std::make_tuple(make_data_static(
            std::make_tuple(I_s(typename DistGen::index_type{},
                                [](const self_type &self) {
                                  return self.get_landa_Map().x().size();
                                },
                                0)),
            std::make_tuple(
                F_s(typename DistGen::label_type{},
                    [](const self_type &s, std::size_t i_a) {
                      return s.get_ProbabilityMap().x()[i_a];
                    }),
                F_s(CT_s<typename DistGen::label_type, s_use, s_probability>{},
                    [](const self_type &s, std::size_t i_a) {
                      return s.get_ProbabilityMap().p().at(
                          s.get_ProbabilityMap().x()[i_a]);
                    })))),
        Compose_static(
            std::make_tuple(I_s(typename DistGen::index_type(),
                                [](const self_type &self) {
                                  return self.get_landa_Map().x().size();
                                },
                                0)),
            F_s(CT_s<typename DistGen::label_type, s_obs, s_probability>{},
                [](const self_type &self, std::size_t i_a) -> auto const & {
                  return self.get_landa_Map().get_Map().at(
                      self.get_landa_Map().x()[i_a]);
                }),
            Beta_Distribution::get_data_index_static()));
  }

  auto &get_ProbabilityMap() const { return p_; }

  auto &get_landa_Map() const { return landaDist_; }

  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "probabilityMap",
                       &self_type::get_ProbabilityMap),
        grammar::field(C<self_type>{}, "historyMap", &self_type::get_landa_Map),
        grammar::field(C<self_type>{}, "target_probability",
                       &self_type::target_prob));
  }

  double target_prob() const { return target_prob_; }

private:
  EV f_;
  Probability_map<DistGen> p_;
  Beta_map<DistGen> landaDist_;
  double target_prob_;

  template <template <typename...> class V, class... As>
  static std::map<DistGen, double> uniform_prior(const V<DistGen, As...> &v) {
    std::map<DistGen, double> out;
    double p = 1.0 / v.size();
    for (auto it = v.begin(); it != v.end(); ++it)
      out[*it] += p;
    return out;
  }
};

template <bool Hastings> class Parallel_Tempering {
public:
  typedef Parallel_Tempering self_type;
  constexpr static auto const className =
      my_static_string("Parallel_Tempering") + str_Hastings<Hastings>();

  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "Metropolis",
                       &self_type::get_Metropolis),
        grammar::field(C<self_type>{}, "Jumping_prob", &self_type::pJump));
  }

  template <class priorModel, class LikelihoodModel,
            class Adapative_Distribution_Generator>
  class samples {
  public:
    samples() = default;
    typedef Thermodynamic_Model<priorModel, LikelihoodModel> Model;
    typedef Thermodynamic_Model_Series<priorModel, LikelihoodModel> ModelSeries;
    auto Scout_i(std::size_t i) const { return i_to_scout.at(i); }
    auto &Scout(std::size_t i) { return scouts_[i_to_scout[i]]; }
    auto &Scout(std::size_t i) const { return scouts_[i_to_scout.at(i)]; }
    void swap(std::size_t i, const ModelSeries &model) {
      std::swap(i_to_scout[i], i_to_scout[i + 1]);
      Scout(i).set_beta(model.beta(i));
      Scout(i + 1).set_beta(model.beta(i + 1));
    }

    double Accept(std::size_t i, const ModelSeries &model) const {
      double logA =
          (model.beta(i) - model.beta(i + 1)) *
          (Scout(i + 1).likelihood().logL() - Scout(i).likelihood().logL());
      std::cerr << "\n thermologA" << logA << "  " << model.beta(i) << " "
                << model.beta(i + 1) << "  " << Scout(i).likelihood().logL()
                << " " << Scout(i + 1).likelihood().logL() << "\n";
      return std::min(1.0, std::exp(logA));
    }
    std::size_t size() const { return scouts_.size(); }

    template <class Random_Engine>
    static myOptional<samples, reg_tag>
    evaluate(const Adaptive_Metropolis_H<Hastings> &adm,
             const Thermodynamic_Model_Series<priorModel, LikelihoodModel> &m,
             std::vector<Random_Engine> &mt,
             std::vector<Adapative_Distribution_Generator> &adaptive) {
      typedef myOptional<samples, reg_tag> Op;
      auto itoscout = init_map(m.size());
      std::vector<mcmc_sample_t<Model, Adapative_Distribution_Generator>>
          myscouts(m.size());
      std::vector<Op_void> res(m.size());
      std::cerr << "\n por empezar\n";
      std::vector<std::stringstream> os(m.size());

#pragma omp parallel for
      for (std::size_t i = 0; i < m.size(); ++i) {
        std::cerr << "\n" << i << " aqui\n";
        auto s = adm.start(m.model(i), mt[i], adaptive[i], os[i]);
        //     std::cerr<<"\n"<<i <<" next\n";
        //     std::cerr<<"\n"<<i <<s.value()<<" value\n";

        res[i] = s;
        //     std::cerr<<"\n"<<i <<" res\n";
        myscouts[i] = std::move(s).value();
        //     std::cerr<<"\n"<<i <<" scouts\n";
      }
      for (auto &ss : os)
        std::cerr << ss.str();
      auto ops = consolidate(std::move(res));
      if (ops.has_value())
        return Op(samples(std::move(itoscout), std::move(myscouts)));
      else
        return Op(false, "fails to build a scout :" + ops.error());
    }

    samples(std::map<std::size_t, std::size_t> &&i_to_scouts,
            std::vector<mcmc_sample_t<Model, Adapative_Distribution_Generator>>
                &&scouts)
        : i_to_scout{std::move(i_to_scouts)}, scouts_{std::move(scouts)} {}

    typedef samples self_type;
    typedef Parallel_Tempering enclosing_type;
    constexpr static auto className =
        enclosing_type::className + my_static_string("_samples");
    static auto get_constructor_fields() {
      return std::make_tuple(
          grammar::field(C<self_type>{}, "i_to_scout",
                         &self_type::get_i_to_scout),
          grammar::field(C<self_type>{}, "scouts", &self_type::get_scouts));
    }
    std::map<std::size_t, std::size_t> const &get_i_to_scout() const {
      return i_to_scout;
    }
    std::vector<mcmc_sample_t<Model, Adapative_Distribution_Generator>> const &
    get_scouts() const {
      return scouts_;
    }

    static std::vector<std::string> data_index_titles() {
      std::vector<std::string> out;
      out.push_back("i_beta");
      out.push_back("i_scout");
      out = concatenate_unique(
          std::move(out),
          mcmc_sample_t<Model,
                        Adapative_Distribution_Generator>::data_index_titles());
      return out;
    }

    static Data_Index_scheme data_index() {
      Data_Index_scheme out;
      out.push_back("beta", {"i_beta"});
      out.push_back("scout_number", {"i_beta"});
      out = concatenate(
          std::move(out),
          mcmc_sample_t<Model, Adapative_Distribution_Generator>::data_index());
      return out;
    }

    static auto get_data_index_static() {

      return Concatenate_tuple_static(
          std::make_tuple(make_data_static(
              std::make_tuple(I_s(i_beta{},
                                  [](const self_type &self) {
                                    return self.get_i_to_scout().size();
                                  },
                                  0)),
              std::make_tuple(F_s(s_scout_number{},
                                  [](const self_type &self, std::size_t i_b) {
                                    return self.Scout_i(i_b);
                                  })))),
          Compose_static(
              std::make_tuple(I_s(i_beta{}, &self_type::size, 0)),
              [](const self_type &self, std::size_t i_b) -> auto const & {
                return self.Scout(i_b);
              },
              mcmc_sample_t<Model, Adapative_Distribution_Generator>::
                  get_data_index_static()));
    }

    auto getIndexedData(const std::set<std::string> &,
                        const std::set<std::string> &) const {
      Data_Index_point out;
      out.push_back("i_beta", scouts_.size());
      auto ind = scouts_[0].getIndexedData;
    }

    static std::vector<std::string> data_titles() {
      std::vector<std::string> out;
      out.push_back("i_beta");
      out.push_back("i_scout");
      auto v =
          mcmc_sample_t<Model, Adapative_Distribution_Generator>::data_title();
      out.insert(out.end(), v.begin(), v.end());
      return out;
    }

    bool append_data(const std::string &idname, std::size_t, Model &) {
      std::ofstream f(idname);
      for (std::size_t is = 0; is < size(); ++is) {
      }
      return true;
    }

  private:
    std::map<std::size_t, std::size_t> i_to_scout;
    std::vector<mcmc_sample_t<Model, Adapative_Distribution_Generator>> scouts_;

    static std::map<std::size_t, std::size_t> init_map(std::size_t n) {
      std::map<std::size_t, std::size_t> out;
      for (std::size_t i = 0; i < n; ++i)
        out[i] = i;
      return out;
    }
  };

  template <class Model, class samples>
  static double Evidence(const Model &m, const samples &s) {
    auto n = m.size();
    double beta0 = m.beta(0);
    double sum = 0;
    double sumdb = 0;
    double logLik0 = s.Scout(0).likelihood().logL();
    for (std::size_t i = 1; i < m.size(); ++i) {
      double beta = m.beta(i);
      double db = beta0 - beta;
      double logLik = s.Scout(i).likelihood().logL();
      sum += db * (logLik0 + logLik) / 2;
      sumdb += db;
      if ((i == n - 1) && (beta > 0)) {
        double db0 = beta;
        sum +=
            (db0 * (1 + db0 / db / 2)) * logLik - sqr(db0) / db / 2 * logLik0;
        sumdb += db0;
      }
      beta0 = beta;
      logLik0 = logLik;
    }
    return sum;
  }

  template <class priorModel, class LikelihoodModel,
            class Adapative_Distribution_Generator, class Random_engine>
  myOptional<
      samples<priorModel, LikelihoodModel, Adapative_Distribution_Generator>,
      reg_tag>
  start(const Thermodynamic_Model_Series<priorModel, LikelihoodModel> &model,
        Random_engine &mt,
        std::vector<Adapative_Distribution_Generator> &g) const {
    return samples<priorModel, LikelihoodModel,
                   Adapative_Distribution_Generator>::evaluate(get_Metropolis(),
                                                               model, mt, g);
  }

  template <class priorModel, class LikelihoodModel, class RandomEngine>
  static std::vector<RandomEngine> start_Random(
      const Thermodynamic_Model_Series<priorModel, LikelihoodModel> &model,
      RandomEngine &mt) {
    return mts(mt, model.size());
  }

  template <class priorModel, class LikelihoodModel,
            class Adapative_Distribution_Generator>
  Op_void
  next(const Thermodynamic_Model_Series<priorModel, LikelihoodModel> &model,
       std::vector<std::mt19937_64> &mt,
       std::vector<Adapative_Distribution_Generator> &adaptives,
       samples<priorModel, LikelihoodModel, Adapative_Distribution_Generator>
           &current) const {
    std::vector<Op_void> res(model.size());
    std::vector<std::stringstream> ss(model.size());
#pragma omp parallel for
    for (std::size_t i = 0; i < model.size(); ++i) {
      res[i] = get_Metropolis().next(model.model(i), mt[i], adaptives[i],
                                     current.Scout(i), ss[i]);
    }

    for (auto &s : ss)
      std::cerr << s.str();
    auto ops = consolidate(std::move(res));
    if (!ops.has_value())
      return Op_void(false, "get_Metropolis fails: " + ops.error());
    std::uniform_real_distribution<> U;
    double r = U(mt[0]);

    std::cerr << " test thermo r=" << r << "P_jump_=" << P_jump_;
    if (r < P_jump_) {
      std::size_t j;
      if (U(mt[0]) < 0.5)
        j = 0;
      else
        j = 1;

      std::cerr << " j=" << j << "\n";
#pragma omp parallel for
      for (std::size_t i = j; i < model.size() - 1; i += 2) {
        double A = current.Accept(i, model);
        double r = std::uniform_real_distribution<>(0, 1)(mt[i]);
        if (r < A)
          current.swap(i, model);
      }
    }
    return Op_void(true, "");
  }

  template <class priorModel, class LikelihoodModel, class Samples>
  auto compute_PartialDLikelihood(
      const Thermodynamic_Model_Series<priorModel, LikelihoodModel> &model,
      const Samples &current, std::ostream &os) const {
    std::vector<std::decay_t<decltype(
        model.model(0)
            .compute_PartialDLikelihood(current.Scout(0).x(), current.Scout(0),
                                        std::declval<std::ostream &>())
            .value())>>
        out(model.size());

    std::vector<std::stringstream> ss(model.size());

#pragma omp parallel for
    for (std::size_t i = 0; i < model.size(); ++i) {
      out[i] = model.model(i)
                   .compute_PartialDLikelihood(current.Scout(i).x(),
                                               current.Scout(i), ss[i])
                   .value();
    }
    for (std::size_t i = 0; i < model.size(); ++i)
      os << ss[i].str();
    return out;
  }

  Parallel_Tempering(Adaptive_Metropolis_H<Hastings> amh, double P_jump)
      : amh_{std::move(amh)}, P_jump_{P_jump} {}

  Adaptive_Metropolis_H<Hastings> const &get_Metropolis() const { return amh_; }

  auto pJump() const { return P_jump_; }
  Parallel_Tempering() = default;

private:
  Adaptive_Metropolis_H<Hastings> amh_;
  double P_jump_;
};
template <class Model, class Adaptive_Mover> class Ensemble_samples {
public:
  typedef Ensemble_samples self_type;
  typedef typename Model::subModel subModel;

  typedef Cs<Model, Adaptive_Mover> template_types;
  constexpr static auto const className =
      my_static_string("Ensemble_Metropolis_Hastings_samples") +
      my_trait<template_types>::className;

  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "ij_history", &self_type::ij_history),
        grammar::field(C<self_type>{}, "scouts", &self_type::scouts));
  }
  static Data_Index_scheme data_index() {
    Data_Index_scheme out;
    out.push_back("n_scout", {"i_beta", "i_scout"});
    auto vemcee = Insert_Index(
        emcee_sample<subModel, Adaptive_Mover>::data_index(), {"i_beta"});
    out = concatenate(std::move(out), std::move(vemcee));
    return out;
  }

  static auto get_data_index_static() {
    return Concatenate_tuple_static(
        std::make_tuple(make_data_static(
            std::make_tuple(
                I_s(i_beta{}, [](const self_type &self) { return self.size(); },
                    0),
                I_s(i_walker{},
                    [](const self_type &self, std::size_t i_beta) {
                      return self.Scout(i_beta).numWalkers();
                    },
                    0)),
            std::make_tuple(
                F_s(s_scout_number{},
                    [](const self_type &self, std::size_t i_b,
                       std::size_t i_w) { return self.Scout_i(i_b, i_w); })))),
        Compose_static(
            std::make_tuple(
                I_s(i_beta{}, [](const self_type &self) { return self.size(); },
                    0)),
            [](const self_type &self, std::size_t i_b) -> auto const & {
              return self.Scout(i_b);
            },
            emcee_sample<subModel, Adaptive_Mover>::get_data_index_static()));
  }

  std::size_t Scout_i(std::size_t iscout, std::size_t iwalker) const {
    return ij_history()[iscout][iwalker].first * Scout(0).numWalkers() +
           ij_history()[iscout][iwalker].second;
  }

  emcee_sample<subModel, Adaptive_Mover> &Scout(std::size_t i) {
    return scouts_[i];
  }
  emcee_sample<subModel, Adaptive_Mover> const &Scout(std::size_t i) const {
    return scouts_[i];
  }

  void swap(std::size_t i, std::size_t j, std::size_t k, const Model &model) {
    std::cerr << "\n swap i=" << i << "  j=" << j << "  k=" << k << "  log_j="
              << Scout(i).Walker(j).likelihood_result().likelihood().logL()
              << "log_k="
              << Scout(i + 1).Walker(k).likelihood_result().likelihood().logL()
              << "\n";
    std::swap(Scout(i).Walker(j), Scout(i + 1).Walker(k));
    std::swap(ij_history_[i][j], ij_history_[i + 1][k]);
    Scout(i).Walker(j).set_beta(model.beta(i));
    Scout(i + 1).Walker(k).set_beta(model.beta(i + 1));
    std::cerr << "result i=" << i << "  j=" << j << "  k=" << k << "  log_j="
              << Scout(i).Walker(j).likelihood_result().likelihood().logL()
              << "log_k="
              << Scout(i + 1).Walker(k).likelihood_result().likelihood().logL()
              << "\n";
  }

  std::size_t size() const { return scouts_.size(); }

  double Accept(std::size_t i, std::size_t j, std::size_t k,
                const Model &model) {
    double logA =
        -(model.beta(i) - model.beta(i + 1)) *
        (Scout(i).Walker(j).likelihood_result().likelihood().logL() -
         Scout(i + 1).Walker(k).likelihood_result().likelihood().logL());
    return std::min(1.0, std::exp(logA));
  }

  auto &ij_history() const { return ij_history_; }
  auto &scouts() const { return scouts_; }

  Ensemble_samples(
      const std::vector<std::vector<std::pair<std::size_t, std::size_t>>>
          &ij_history,
      const std::vector<emcee_sample<subModel, Adaptive_Mover>> &scouts)
      : ij_history_{ij_history}, scouts_{scouts} {}

  Ensemble_samples(std::vector<std::vector<std::pair<std::size_t, std::size_t>>>
                       &&ij_history,
                   std::vector<emcee_sample<subModel, Adaptive_Mover>> &&scouts)
      : ij_history_{std::move(ij_history)}, scouts_{std::move(scouts)} {}

  Ensemble_samples() = default;

  template <class Random_Engine>
  static myOptional_t<Ensemble_samples>
  evaluate(const Ensemble_Metropolis_Hastings &emh, const Model &m,
           std::vector<Random_Engine> &_mt,
           std::vector<Adaptive_Mover> &adaptives) {
    typedef myOptional_t<Ensemble_samples> Op;
    auto ij_history = init_history(m.size(), emh.numWalkers());

    std::vector<emcee_sample<subModel, Adaptive_Mover>> scouts_(m.size());
    std::vector<Op_void> op(m.size(), Op_void(true, ""));
    for (std::size_t i = 0; i < m.size(); ++i) {
      auto s = emh.start(m.model(i), _mt[i], adaptives[i]);
      if (s)
        scouts_[i] = std::move(s).value();
      else
        op[i] = Op_void(false, s.error());
    }
    auto res = consolidate(std::move(op));
    if (res)
      return Op(Ensemble_samples(std::move(ij_history), std::move(scouts_)));
    else {
      return Op(false, res.error());
    }
  }

private:
  std::vector<std::vector<std::pair<std::size_t, std::size_t>>> ij_history_;
  std::vector<emcee_sample<subModel, Adaptive_Mover>> scouts_;

  std::vector<std::mt19937_64> build_mts(std::mt19937_64 &mt, std::size_t n) {
    std::uniform_int_distribution<typename std::mt19937_64::result_type> useed;
    std::vector<std::mt19937_64> out(n);
    for (std::size_t i = 0; i < n; ++i)
      out[i].seed(useed(mt));
    return out;
  }

  static std::vector<std::vector<std::pair<std::size_t, std::size_t>>>
  init_history(std::size_t nScouts, std::size_t nWalkers) {
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> out(
        nScouts, std::vector<std::pair<std::size_t, std::size_t>>(nWalkers));
    for (std::size_t i = 0; i < nScouts; ++i)
      for (std::size_t j = 0; j < nWalkers; ++j)
        out[i][j] = std::pair(i, j);
    return out;
  }
};

class Ensemble_Parallel_Tempering {

public:
  typedef Ensemble_Parallel_Tempering self_type;
  constexpr static auto const className =
      my_static_string("Ensemble_Parallel_Tempering");

  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "Ensemble_Metropolis",
                       &self_type::get_Ensamble_Metropolis),
        grammar::field(C<self_type>{}, "pJump", &self_type::pJump));
  }

  template <class priorModel, class LikelihoodModel, class RandomEngine>
  static std::vector<RandomEngine> start_Random(
      const Thermodynamic_Model_Series<priorModel, LikelihoodModel> &model,
      RandomEngine &mt) {
    return mts(mt, model.size());
  }

  template <class Model, class Random_Engine, class AdaptiveMover>
  myOptional<Ensemble_samples<Model, AdaptiveMover>, reg_tag>
  start(const Model &model, std::vector<Random_Engine> &mt,
        std::vector<AdaptiveMover> adaptives) const {
    return Ensemble_samples<Model, AdaptiveMover>::evaluate(
        get_Ensamble_Metropolis(), model, mt, adaptives);
  }

  template <class Model, class Random_Engine, class AdaptiveMover>
  Op_void next(const Model &model, std::vector<Random_Engine> &mt,
               std::vector<AdaptiveMover> &adaptives,
               Ensemble_samples<Model, AdaptiveMover> &current) const {

    std::vector<Op_void> op(model.size());
#pragma omp parallel for
    for (std::size_t i = 0; i < model.size(); ++i) {
      op[i] = get_Ensamble_Metropolis().next(model.model(i), mt[i],
                                             adaptives[i], current.Scout(i));
    }

    auto res = consolidate(std::move(op));
    if (res) {
      std::uniform_real_distribution<> U;
      double r = U(mt[0]);

      if (r < P_jump_) {
        std::size_t j;
        if (U(mt[0]) < 0.5)
          j = 0;
        else
          j = 1;

#pragma omp parallel for
        for (std::size_t i = j; i < model.size() - 1; i += 2) {
          std::size_t j = std::uniform_int_distribution<std::size_t>(
              0, current.Scout(i).numWalkers() - 1)(mt[i]);
          std::size_t k = std::uniform_int_distribution<std::size_t>(
              0, current.Scout(i + 1).numWalkers() - 1)(mt[i]);
          double A = current.Accept(i, j, k, model);
          double r = std::uniform_real_distribution<>(0, 1)(mt[i]);
          if (r < A)
            current.swap(i, j, k, model);
        }
      }
    }
    return res;
  }

  Ensemble_Metropolis_Hastings const &get_Ensamble_Metropolis() const {
    return emh_;
  }
  Ensemble_Parallel_Tempering(Ensemble_Metropolis_Hastings ensamble_mh,
                              double probability_jump)
      : emh_{ensamble_mh}, P_jump_{probability_jump} {}

  double pJump() const { return P_jump_; }
  Ensemble_Parallel_Tempering() = default;
  template <class Model, class mySample>
  auto compute_PartialDLikelihood(const Model &M, const mySample &current,
                                  std::ostream &os) const {
    std::vector<decltype(get_Ensamble_Metropolis().compute_PartialDLikelihood(
        M.model(0), current.Scout(0), os))>
        out(M.size());

#pragma omp parallel for
    for (std::size_t i = 0; i < M.size(); ++i) {
      out[i] = get_Ensamble_Metropolis().compute_PartialDLikelihood(
          M.model(i), current.Scout(i), os);
    }
    return out;
  }

private:
  Ensemble_Metropolis_Hastings emh_;
  double P_jump_;
};

template <class Random_Generator, class Metropolis_algorithm, class ModelSeries,
          class Distribution_Generator, class StateSamplingFunctor,
          class GenSamplingFunctor, class AnaSamplingFunctor>
struct OutputGenerator {
  typedef Cs<Metropolis_algorithm, ModelSeries, Distribution_Generator>
      template_types;
  constexpr static auto const className = my_static_string("OutputGenerator_") +
                                          my_trait<template_types>::className;

  typedef decltype(std::declval<Metropolis_algorithm &>().start(
      std::declval<ModelSeries &>(), std::declval<Random_Generator &>(),
      std::declval<Distribution_Generator &>())) OpState;
  typedef std::decay_t<decltype(std::declval<OpState>().value())> State;

  static constexpr const auto info_end = "--END of Info--";
  void save(const std::string &idname, const std::string &info,
            const Metropolis_algorithm &mcmc, const ModelSeries &M,
            std::size_t nsamples) {
    std::ofstream f(idname + run_ext);
    f << className.str() << std::endl;
    f << info << std::endl;
    f << info_end << std::endl;
    f << mcmc << std::endl;
    f << M << std::endl;
    f << "nsamples" << std::endl;
    f << nsamples << std::endl;
    f.close();
  }
  void load(const std::string &idname, std::string &info,
            Metropolis_algorithm &mcmc, ModelSeries &M, std::size_t &nsamples) {
    std::ifstream f(idname + run_ext);
    std::string line;
    std::getline(f, line);
    assert((line.find(className.str()) != line.npos));
    std::getline(f, line);
    info.clear();
    info = line;
    std::getline(f, line);
    while ((f) && (line.find(info_end) == line.npos)) {
      info += "\n" + line;
      std::getline(f, line);
    }
    assert(line.find(info_end) != line.npos);
    f >> mcmc;
    std::getline(f, line);
    f >> M;
    std::getline(f, line);
    std::getline(f, line);
    assert(line == "nsamples");
    f >> nsamples;
    f.close();
  }

  static constexpr const char *save_ext = "_save.txt";
  static constexpr const char *run_ext = "_run.txt";

  void save(const std::string &idname, std::size_t isample, const State &state,
            const Distribution_Generator &G) {
    std::string tmp = idname + save_ext + ".tmp";
    std::string old = idname + save_ext + ".old";
    std::string name = idname + save_ext;
    std::ofstream f(tmp);
    f << "isample"
      << "\n";
    f << isample << "\n";
    f << state << "\n";
    f << G << "\n";
    f.close();
    std::filesystem::remove(old);
    if (std::filesystem::exists(name))
      std::filesystem::rename(name, old);
    std::filesystem::rename(tmp, name);
  }
  bool load(const std::string &idname, std::size_t &isample, State &state,
            Distribution_Generator &G) {
    std::string name = idname + save_ext;
    std::ifstream f(name);
    std::string line;
    std::getline(f, line);
    assert(line == "isample");
    f >> isample;
    std::getline(f, line);
    f >> state;
    std::getline(f, line);
    f >> G;
    bool res = f.good();
    f.close();
    return res;
  }

  static constexpr auto state_ext = "state";
  static constexpr auto adap_dist_ext = "ADG";

  template <class State_Data_Index_tuple, class Gen_Data_Index_tuple,
            class Ana_Data_Index_tuple>
  void print_title(const State_Data_Index_tuple &statedata,
                   const Gen_Data_Index_tuple &gendata,
                   const Ana_Data_Index_tuple &anadata,
                   const std::string &idname) {

    //       typedef typename State_Data_Index_tuple::format state_format;
    std::apply(
        [&idname, this](auto &... d) {
          (d.print_title(idname + "_sta_", ".txt", sep) && ...);
        },
        statedata);
    std::apply(
        [&idname, this](auto &... d) {
          (d.print_title(idname + "_gen_", ".txt", sep) && ...);
        },
        gendata);
    std::apply(
        [&idname, this](auto &... d) {
          (d.print_title(idname + "_ana_", ".txt", sep) && ...);
        },
        anadata);
  }
  bool is_sampling(std::size_t isample) const {
    return max_index_state(isample) > 0;
  }

  std::size_t max_index_state(std::size_t i_sample) const {
    return std::invoke(f_sta_, i_sample);
  }

  std::size_t max_index_gen(std::size_t i_sample) const {
    return std::invoke(f_gen_, i_sample);
  }
  std::size_t max_index_ana(std::size_t i_sample) const {
    return std::invoke(f_ana_, i_sample);
  }

  template <class State_Data_Index_tuple, class Gen_Data_Index_tuple>
  void print_element(const State_Data_Index_tuple &statedata,
                     const Gen_Data_Index_tuple &gendata,
                     const std::string &idname, std::size_t isample,
                     const ModelSeries &M, const State &state,
                     const Distribution_Generator &G) {

    std::tuple<std::size_t, State const *, ModelSeries const *> s(isample,
                                                                  &state, &M);

    auto sta_max_index = max_index_state(isample);
    std::apply(
        [&s, &idname, &sta_max_index, this](auto &... t) {
          ((t.NumIndex() <= sta_max_index
                ? t.print_data(s, idname + "_sta_", ".txt", sep)
                : true) &&
           ...);
        },
        statedata);

    auto gen_max_index = max_index_gen(isample);
    if (gen_max_index > 0) {
      std::tuple<std::size_t, Distribution_Generator const *,
                 ModelSeries const *>
          g(isample, &G, &M);
      std::apply(
          [&g, &idname, &gen_max_index, this](auto &... t) {
            ((t.NumIndex() <= gen_max_index
                  ? t.print_data(g, idname + "_gen_", ".txt", sep)
                  : true) &&
             ...);
          },
          gendata);
    }
  }

  template <class State_Data_Index_tuple, class Gen_Data_Index_tuple,
            class Ana_Data_Index_tuple, class AnaSample>
  void print_element(const State_Data_Index_tuple &statedata,
                     const Gen_Data_Index_tuple &gendata,
                     const Ana_Data_Index_tuple &anadata,
                     const std::string &idname, std::size_t isample,
                     const ModelSeries &M, const State &state,
                     const Distribution_Generator &G, const AnaSample &ana) {
    print_element(statedata, gendata, idname, isample, M, state, G);

    auto ana_max_index = max_index_ana(isample);
    if (ana_max_index > 0) {
      std::tuple<std::size_t, AnaSample const *, ModelSeries const *> s(
          isample, &ana, &M);

      std::apply(
          [&s, &idname, &ana_max_index, this](auto &... t) {
            ((t.NumIndex() <= ana_max_index
                  ? t.print_data(s, idname + "_ana_", ".txt", sep)
                  : true) &&
             ...);
          },
          anadata);
    }
  }

  OutputGenerator(Cs<Random_Generator>, Cs<Metropolis_algorithm>,
                  Cs<ModelSeries>, Cs<Distribution_Generator>,
                  const StateSamplingFunctor &fsta,
                  const GenSamplingFunctor &fgen,
                  const AnaSamplingFunctor &fana, std::ostream &ost,
                  std::string separator = "\t")
      : f_sta_{fsta}, f_gen_{fgen}, f_ana_{fana}, os(ost), sep{separator} {}

private:
  StateSamplingFunctor f_sta_;
  GenSamplingFunctor f_gen_;
  AnaSamplingFunctor f_ana_;
  std::ostream &os;

  std::string sep;
};

template <class self_type, class ModelSeries>
auto compose_state_model_static() {
  typedef std::tuple<std::size_t, self_type const *, ModelSeries const *>
      self_type_model;
  //      typedef typename decltype(self_type::get_data_index_static())::dsgfs
  //      re;

  return Concatenate_tuple_static(
      std::make_tuple(make_data_static(
          std::make_tuple(
              I_s(i_sample{}, [](const self_type_model &) { return 1; },
                  [](const self_type_model &s, std::size_t) -> std::size_t {
                    return std::get<0>(s);
                  }),
              I_s(i_beta{},
                  [](const self_type_model &self, std::size_t) {
                    return std::get<2>(self)->betas().size();
                  },
                  [](const self_type_model &self, std::size_t,
                     std::size_t i_beta) {
                    return std::get<2>(self)->beta(i_beta);
                  }),
              I_s(i_Parameter{},
                  [](const self_type_model &self, std::size_t, std::size_t) {
                    return std::get<2>(self)->name_size();
                  },
                  [](const self_type_model &self, std::size_t, std::size_t,
                     std::size_t i_par) {
                    return std::get<2>(self)->name(i_par);
                  }),
              I_s(CT_s<i_Parameter, s_Transpose>{},
                  [](const self_type_model &self, std::size_t, std::size_t,
                     std::size_t) { return std::get<2>(self)->name_size(); },
                  [](const self_type_model &self, std::size_t, std::size_t,
                     std::size_t, std::size_t i_par) {
                    return std::get<2>(self)->name(i_par);
                  })),
          std::tuple<>())),
      Compose_static(
          std::make_tuple(
              I_s(i_sample{}, [](const self_type_model &) { return 1; },
                  [](const self_type_model &s, std::size_t) -> std::size_t {
                    return std::get<0>(s);
                  })),
          [](const self_type_model &self, std::size_t) -> auto const & {
            return *std::get<1>(self);
          },
          self_type::get_data_index_static()));
}

template <class self_type, class ModelSeries>
auto compose_ana_state_model_static(Cs<std::vector<self_type>>,
                                    Cs<ModelSeries>) {
  typedef std::tuple<std::size_t, std::vector<self_type> const *,
                     ModelSeries const *>
      self_type_model;
  return Concatenate_tuple_static(
      Compose_static(
          std::make_tuple(
              I_s(i_sample{}, [](const self_type_model &) { return 1; },
                  [](const self_type_model &s, std::size_t) -> std::size_t {
                    return std::get<0>(s);
                  }),

              I_s(i_beta{},
                  [](const self_type_model &self, std::size_t) {
                    return std::get<2>(self)->betas().size();
                  },
                  [](const self_type_model &self, std::size_t,
                     std::size_t i_beta) {
                    return std::get<2>(self)->beta(i_beta);
                  })),
          [](const self_type_model &self, std::size_t,
             std::size_t i) -> auto const & { return (*std::get<1>(self))[i]; },
          self_type::get_data_index_static()),
      std::make_tuple(make_data_static(

          std::make_tuple(
              I_s(i_sample{}, [](const self_type_model &) { return 1; },
                  [](const self_type_model &s, std::size_t) -> std::size_t {
                    return std::get<0>(s);
                  }),
              I_s(i_beta{},
                  [](const self_type_model &self, std::size_t) {
                    return std::get<2>(self)->betas().size();
                  },
                  [](const self_type_model &self, std::size_t,
                     std::size_t i_beta) {
                    return std::get<2>(self)->beta(i_beta);
                  }),
              I_s(i_measure{},
                  [](const self_type_model &self, std::size_t, std::size_t) {
                    return std::get<2>(self)->betas().size();
                  },
                  [](const self_type_model &self, std::size_t, std::size_t,
                     std::size_t i_beta) {
                    return std::get<2>(self)->beta(i_beta);
                  }),
              I_s(i_Parameter{},
                  [](const self_type_model &self, std::size_t, std::size_t,
                     std::size_t) { return std::get<2>(self)->name_size(); },
                  [](const self_type_model &self, std::size_t, std::size_t,
                     std::size_t, std::size_t i_par) {
                    return std::get<2>(self)->name(i_par);
                  }),
              I_s(CT_s<i_Parameter, s_Transpose>{},
                  [](const self_type_model &self, std::size_t, std::size_t,
                     std::size_t,
                     std::size_t) { return std::get<2>(self)->name_size(); },
                  [](const self_type_model &self, std::size_t, std::size_t,
                     std::size_t, std::size_t, std::size_t i_par) {
                    return std::get<2>(self)->name(i_par);
                  })),
          std::tuple<>())));
}

template <class self_type, class ModelSeries>
auto compose_ana_state_model_static(Cs<std::vector<std::vector<self_type>>>,
                                    Cs<ModelSeries>)

{
  typedef std::tuple<std::size_t, std::vector<std::vector<self_type>> const *,
                     ModelSeries const *>
      self_type_model;

  return Concatenate_tuple_static(
      Compose_static(
          std::make_tuple(
              I_s(i_sample{}, [](const self_type_model &) { return 1; },
                  [](const self_type_model &s, std::size_t) -> std::size_t {
                    return std::get<0>(s);
                  }),
              I_s(i_beta{},
                  [](const self_type_model &self, std::size_t) {
                    return std::get<2>(self)->betas().size();
                  },
                  [](const self_type_model &self, std::size_t,
                     std::size_t i_beta) {
                    return std::get<2>(self)->beta(i_beta);
                  }),
              I_s(i_walker{},
                  [](const self_type_model &self, std::size_t,
                     std::size_t i_beta) {
                    return std::get<1>(self)->operator[](i_beta).size();
                  },
                  0)),
          [](const self_type_model &self, std::size_t, std::size_t i_beta,
             std::size_t i_walker) -> auto const & {
            return std::get<1>(self)->operator[](i_beta)[i_walker];
          },
          self_type::get_data_index_static()),
      std::make_tuple(make_data_static(
          std::make_tuple(
              I_s(i_sample{}, [](const self_type_model &) { return 1; },
                  [](const self_type_model &s, std::size_t) -> std::size_t {
                    return std::get<0>(s);
                  }),
              I_s(i_beta{},
                  [](const self_type_model &self, std::size_t) {
                    return std::get<2>(self)->betas().size();
                  },
                  [](const self_type_model &self, std::size_t,
                     std::size_t i_beta) {
                    return std::get<2>(self)->beta(i_beta);
                  }),
              I_s(i_walker{},
                  [](const self_type_model &self, std::size_t,
                     std::size_t i_beta) {
                    return std::get<1>(self)->operator[](i_beta).size();
                  },
                  0),
              I_s(i_measure{},
                  [](const self_type_model &self, std::size_t, std::size_t,
                     std::size_t i_beta) {
                    return std::get<1>(self)->operator[](i_beta).size();
                  },
                  0),
              I_s(i_Parameter{},
                  [](const self_type_model &self, std::size_t, std::size_t,
                     std::size_t,
                     std::size_t) { return std::get<2>(self)->name_size(); },
                  [](const self_type_model &self, std::size_t, std::size_t,
                     std::size_t, std::size_t, std::size_t i_par) {
                    return std::get<2>(self)->name(i_par);
                  }),
              I_s(CT_s<i_Parameter, s_Transpose>{},
                  [](const self_type_model &self, std::size_t, std::size_t,
                     std::size_t, std::size_t,
                     std::size_t) { return std::get<2>(self)->name_size(); },
                  [](const self_type_model &self, std::size_t, std::size_t,
                     std::size_t, std::size_t, std::size_t, std::size_t i_par) {
                    return std::get<2>(self)->name(i_par);
                  })),
          std::tuple<>())));
}

template <class self_type, class ModelSeries> auto compose_adg_model_static() {
  typedef std::tuple<std::size_t, self_type const *, ModelSeries const *>
      self_type_model;
  return Concatenate_tuple_static(
      std::make_tuple(make_data_static(
          std::make_tuple(
              I_s(i_sample{}, [](const self_type_model &) { return 1; },
                  [](const self_type_model &s, std::size_t) -> std::size_t {
                    return std::get<0>(s);
                  }),
              I_s(i_beta{},
                  [](const self_type_model &self, std::size_t) {
                    return std::get<2>(self)->betas().size();
                  },
                  [](const self_type_model &self, std::size_t,
                     std::size_t i_beta) {
                    return std::get<2>(self)->beta(i_beta);
                  })),
          std::tuple<>())),
      Compose_static(
          std::make_tuple(
              I_s(i_sample{}, [](const self_type_model &) { return 1; },
                  [](const self_type_model &s, std::size_t) -> std::size_t {
                    return std::get<0>(s);
                  })),
          [](const self_type_model &self, std::size_t) -> self_type const & {
            return *std::get<1>(self);
          },
          myData_Index<self_type>::get_data_index_static()));
}

template <class State, class Random_generator, class Metropolis_algorithm,
          class ModelSeries, class Distribution_Generator, class GenData,
          class StateData, class AnaData, class OutputGenerator>
Op_void continue_Montecarlo_Markov_Chain(
    const std::string &idname, std::size_t isample, State &state,
    Random_generator &mts, const Metropolis_algorithm &mcmc,
    const ModelSeries &M, Distribution_Generator &G, std::size_t nsamples,
    const StateData &statedata, const GenData &gendata, const AnaData &anadata,
    OutputGenerator &output) {

  output.print_title(statedata, gendata, anadata, idname);
  for (; isample < nsamples; ++isample) {
    auto res = mcmc.next(M, mts, G, state);
    if (!res)
      return Op_void(false, " interrupted at sample i=" + ToString(isample) +
                                res.error());

    if (output.max_index_ana(isample)>0) {
      auto ana_sample = mcmc.compute_PartialDLikelihood(M, state, std::cerr);

      output.print_element(statedata, gendata, anadata, idname, isample, M,
                           state, G, ana_sample);

    } else {
      output.print_element(statedata, gendata, idname, isample, M, state, G);
    }
    output.save(idname, isample, state, G);
  }
  return Op_void(true, "run for " + ToString(nsamples) + " samples");
}

template <class Metropolis_algorithm, class ModelSeries,
          class Distribution_Generator, class OutputGenerator>
Op_void run_Montecarlo_Markov_Chain(const std::string idname,
                                    const std::string &info,
                                    const Metropolis_algorithm &mcmc,
                                    std::mt19937_64 &mt, const ModelSeries &M,
                                    Distribution_Generator &G,
                                    std::size_t nsamples, OutputGenerator &O) {

  O.save(idname, info, mcmc, M, nsamples);

  assert(([&O, &idname, &info, &M, &mcmc, &nsamples]() -> bool {
    std::string info_;
    Metropolis_algorithm mcmc_;
    ModelSeries M_;
    std::size_t nsamples_;
    double eps = std::numeric_limits<double>::epsilon() * 10;
    O.load(idname, info_, mcmc_, M_, nsamples_);
    bool out = true;
    if (!(info_ == info)) {
      std::cerr << " info \n original: " << info << "\n loaded: " << info_;
      out = false;
    }
    if (!are_Equal_v(M_, M, eps, eps, std::cerr))
      out = false;
    if (!are_Equal_v(mcmc_, mcmc, eps, eps, std::cerr))
      out = false;
    return out;
  }()));

  auto mts = mcmc.start_Random(M, mt);
  auto state = mcmc.start(M, mts, G);
  if (!state) {
    std::cerr << state.error();
    return Op_void(false, " does not start at all " + state.error());
  }

  std::size_t isample = 0;

  typedef std::decay_t<decltype(
      mcmc.compute_PartialDLikelihood(M, state.value(), std::cerr))>
      AnaState;

  auto state_data =
      compose_state_model_static<std::decay_t<decltype(state.value())>,
                                 ModelSeries>();
  auto gen_data =
      compose_adg_model_static<Distribution_Generator, ModelSeries>();
  auto ana_data =
      compose_ana_state_model_static(Cs<AnaState>(), Cs<ModelSeries>());

  return continue_Montecarlo_Markov_Chain(idname, isample, state.value(), mts,
                                          mcmc, M, G, nsamples, state_data,
                                          gen_data, ana_data, O);
}

template <class Random_generator, class Metropolis_algorithm, class ModelSeries,
          class Distribution_Generator,
          template <class, class, class, class> class OutputGenerator>
Op_void continue_Montecarlo_Markov_Chain(
    const std::string &idname, const std::string & /*new_info*/,
    OutputGenerator<Random_generator, Metropolis_algorithm, ModelSeries,
                    Distribution_Generator> &output) {
  ModelSeries M;
  std::size_t nsamples;
  Metropolis_algorithm mcmc;
  std::string info_old;
  output.load(idname, info_old, mcmc, M, nsamples);
  Random_generator mts;
  Distribution_Generator G;
  decltype(mcmc.start(M, mts, G).value()) state;
  auto data_output = std::decay_t<decltype(state)>::get_data_index_model();
  std::size_t isample;
  output.load(idname, isample, mts, G, state);
  return continue_Montecarlo_Markov_Chain(idname, isample, state, mts, mcmc, M,
                                          G, nsamples, data_output, output);
}

template <class PriorModel, class Likelihood_Model, class OutputGenerator>
Op_void run_Thermo_Levenberg_ProbVel(
    const std::string id, const std::string &info, const PriorModel &prior,
    const Likelihood_Model &lik, std::mt19937_64 &mt,
    const std::vector<double> &betas, const std::vector<double> &landa,
    const std::vector<std::vector<double>> &landa_50_hill, double pjump,
    double gain_moment, std::size_t nSamples, std::size_t ntrials,
    OutputGenerator &output) {
  typedef Thermodynamic_Model<PriorModel, Likelihood_Model> Model;

  //  typedef std::vector<std::mt19937> RG;
  // typedef  std::vector<
  //     Adaptive_Parameterized_Distribution_Generator<LevenbergMarquardt<Model>>>
  //     Adaptive;
  //  typedef  Thermodynamic_Model_Series<PriorModel, Likelihood_Model> Models;
  //  typedef Parallel_Tempering<true> MCMC;

  Thermodynamic_Model_Series<PriorModel, Likelihood_Model> model(prior, lik,
                                                                 betas);

  Parallel_Tempering<true> PT(Adaptive_Metropolis_H<true>(ntrials), pjump);
  std::vector<LevenbergMarquardt<Model>> las(landa.size());
  for (std::size_t i = 0; i < las.size(); ++i) {
    las[i] = LevenbergMarquardt<Model>(landa[i]);
  }
  std::vector<
      Adaptive_Parameterized_Distribution_Generator<LevenbergMarquardt<Model>>>
      ala(betas.size(),
          Adaptive_Parameterized_Distribution_Generator<
              LevenbergMarquardt<Model>>(las, landa_50_hill, gain_moment));

  return run_Montecarlo_Markov_Chain(id, info, PT, mt, model, ala, nSamples,
                                     output);
}

template <class PriorModel, class Likelihood_Model, class OutputGenerator>
Op_void run_Thermo_Levenberg_Prob(
    const std::string id, const std::string &info, const PriorModel &prior,
    const Likelihood_Model &lik, std::mt19937_64 &mt,
    const std::vector<double> &betas, const std::vector<double> &landa,
    double target_probability, double pjump, std::size_t nSamples,
    std::size_t ntrials, OutputGenerator &output) {
  typedef Thermodynamic_Model<PriorModel, Likelihood_Model> Model;

  //  typedef std::vector<std::mt19937> RG;
  // typedef  std::vector<
  //     Adaptive_Parameterized_Distribution_Generator<LevenbergMarquardt<Model>>>
  //     Adaptive;
  //  typedef  Thermodynamic_Model_Series<PriorModel, Likelihood_Model> Models;
  //  typedef Parallel_Tempering<true> MCMC;

  Thermodynamic_Model_Series<PriorModel, Likelihood_Model> model(prior, lik,
                                                                 betas);

  Parallel_Tempering<true> PT(Adaptive_Metropolis_H<true>(ntrials), pjump);
  std::vector<LevenbergMarquardt<Model>> las(landa.size());
  for (std::size_t i = 0; i < las.size(); ++i) {
    las[i] = LevenbergMarquardt<Model>(landa[i]);
  }
  std::vector<
      Adaptive_Probability_Distribution_Generator<LevenbergMarquardt<Model>>>
      ala(betas.size(),
          Adaptive_Probability_Distribution_Generator<
              LevenbergMarquardt<Model>>(las, target_probability));

  return run_Montecarlo_Markov_Chain(id, info, PT, mt, model, ala, nSamples,
                                     output);
}

template <class PriorModel, class Likelihood_Model, class OutputGenerator>
Op_void run_Thermo_emcee(const std::string id, const std::string info,
                         const PriorModel &prior, const Likelihood_Model &lik,
                         std::mt19937_64 &mt, const std::vector<double> &betas,
                         std::size_t numWalkers,
                         const std::vector<double> &alfas, double pjump,
                         double target_prob, std::size_t nSamples,
                         std::size_t ntrials, OutputGenerator &output) {

  typedef Thermodynamic_Model<PriorModel, Likelihood_Model> Model;
  typedef
      typename evidence::Ensemble_Metropolis_Hastings::adaptive_stretch_mover<
          Model>
          Adaptive;
  // typedef  Thermodynamic_Model_Series<PriorModel, Likelihood_Model> Models;
  // typedef Ensemble_Parallel_Tempering MCMC;
  Thermodynamic_Model_Series<PriorModel, Likelihood_Model> model(prior, lik,
                                                                 betas);

  Ensemble_Parallel_Tempering ET(
      Ensemble_Metropolis_Hastings(numWalkers, ntrials), pjump);
  std::vector<Adaptive> alfa_adap(betas.size(), Adaptive(target_prob, alfas));

  return run_Montecarlo_Markov_Chain(id, info, ET, mt, model, alfa_adap,
                                     nSamples, output);
}

} // namespace evidence

#endif // MYEVIDENCE_H
