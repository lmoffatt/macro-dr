#ifndef LIKELIHOOD_MARKOV_PROCESS_H
#define LIKELIHOOD_MARKOV_PROCESS_H
#include "Experiment.h"
#include "Matrix.h"
#include "myDistributions.h"
#include "mymath.h"
#include "myoptimization.h"
#include "qmodel.h"
#include "matrixcalculations.h"
#include <type_traits>

namespace markov {

struct s_Pmean {
  constexpr static auto const title = my_static_string("Pmean");
};
struct s_Pcov {
  constexpr static auto const title = my_static_string("Pcov");
};
struct s_MacroAlgorithm {
  constexpr static auto const title = my_static_string("MacroAlgorithm");
};
struct i_state {
  constexpr static auto const title = my_static_string("i_state");
};

template<template<class...> class Tr>
class mp_information_;
typedef mp_information_<C> mp_information;

inline bool is_Binomial_Approximation_valid(double N, double p, double q,
                                               double Np_min) {
    return (N * p>= Np_min) && (N * q>= Np_min);
 }


 template<template<class, class, bool>class Error, class Norm, bool diff>
 bool is_Binomial_Approximation_valid(double N, const Error<double,Norm,diff>& p, const Error<double,Norm,diff>&  q,
                                             double Np_min) {
     return (N * p.center()>= Np_min) && (N * q.center()>= Np_min);
 }


 inline std::pair<bool, bool>
is_Binomial_Approximation_valid(double N, double p, double q, double Np_min,
                                double Vp_min) {
    return {is_Binomial_Approximation_valid(N, p, q, Np_min),
            is_Binomial_Approximation_valid(N, p, q, Vp_min)};
}


 template<template<class, class, bool>class Error, class Norm, bool diff>
std::pair<bool,bool> is_Binomial_Approximation_valid(double N, const Error<double,Norm,diff>& p, const Error<double,Norm,diff>&  q,
                                     double Np_min,
                                     double Vp_min) {
    return {is_Binomial_Approximation_valid(N, p, q, Np_min),
            is_Binomial_Approximation_valid(N, p, q, Vp_min)};
}

template<template<class...> class Tr>
class mp_information_ {
public:
    template <class T> using Tr_t=typename Tr<T>::type;
private:

    Tr_t<double> y_mean_;
    Tr_t<double> y_var_;
    Tr_t<double> plogL_;
    Tr_t<double> eplogL_;
    Tr_t<double> vplogL_;

public:

    typedef mp_information_ self_type;
    constexpr static auto className = my_template_trait<Tr>::className+my_static_string("mp_state_information");
    static auto get_constructor_fields() {
        return std::make_tuple(
            grammar::field(C<self_type>{}, "y_mean", &self_type::y_mean),
            grammar::field(C<self_type>{}, "y_var", &self_type::y_var),
            grammar::field(C<self_type>{}, "plogL", &self_type::plogL),
            grammar::field(C<self_type>{}, "eplogL", &self_type::eplogL),
            grammar::field(C<self_type>{}, "vplogL", &self_type::vplogL));
    }




    static mp_information_
    adjust(Tr_t<double> y_mean__, Tr_t<double> y_var__, Tr_t<double> plogL__, Tr_t<double> eplogL__,
           Tr_t<double> vplogL__, [[maybe_unused]] double min_p, double min_var) {
        return mp_information(y_mean__, variance_value::adjust(y_var__, min_var), plogL__, eplogL__,
            vplogL__);
    }
    mp_information_(Tr_t<double> y_mean__, Tr_t<double> y_var__, Tr_t<double> plogL__,
                         Tr_t<double> eplogL__, Tr_t<double> vplogL__)
        : y_mean_{y_mean__}, y_var_{y_var__}, plogL_{plogL__}, eplogL_{eplogL__},
          vplogL_{vplogL__} {}

    mp_information_() = default;

    Tr_t<double> y_mean() const { return y_mean_; }
    Tr_t<double> y_var() const { return y_var_; }
    Tr_t<double> plogL() const { return plogL_; }
    Tr_t<double> eplogL() const { return eplogL_; }
    Tr_t<double> vplogL() const { return vplogL_; }

    template <class... Indexes> static auto get_data_index_static(Indexes...) {
        return std::make_tuple(
            make_data_static(
                std::tuple<>(),
                std::make_tuple(
                    F_s(CT_s<s_pred, s_mean>{}, &self_type::y_mean),
                    F_s(CT_s<s_pred, s_variance>{}, &self_type::y_var),
                    F_s(CT_s<s_obs, s_logL>{}, &self_type::plogL),
                    F_s(CT_s<s_pred, s_logL>{}, &self_type::eplogL),
                    F_s(CT_s<s_variance, s_logL>{}, &self_type::vplogL)))
           );
    }
};


template<template<class...> class Tr>
class mp_macro_information_;
typedef mp_macro_information_<C> mp_macro_information;

template<template<class...> class Tr>
class mp_macro_information_ : public mp_information_<Tr>
{
public:
    template <class T> using Tr_t=typename Tr<T>::type;
private:

    Tr_t<M_Matrix<double>> P_mean_;
    Tr_t<M_Matrix<double>> P_cov_;


public:
    typedef mp_macro_information_ self_type;
    typedef Tr_t<mp_information> base_type;
    constexpr static auto className = my_template_trait<Tr>::className+ my_static_string("mp_macro_information");
    static auto get_constructor_fields() {
        return std::tuple_cat(
            std::make_tuple(
                grammar::field(C<self_type>{}, "P_mean", &self_type::P_mean),
                grammar::field(C<self_type>{}, "P_cov", &self_type::P_cov)
                    ),
        base_type::get_constructor_fields());
    }


    static mp_macro_information_
    adjust(Tr_t<M_Matrix<double>> &&P_mean__, Tr_t<M_Matrix<double>> &&P_cov__,
           Tr_t<double> y_mean__, Tr_t<double> y_var__, Tr_t<double> plogL__, Tr_t<double> eplogL__,
           Tr_t<double> vplogL__, double min_p, double min_var) {
        return mp_macro_information_(
            Probability_distribution_new_<Tr>::normalize(std::move(P_mean__), min_p),
            Probability_distribution_covariance_new_<Tr>::normalize(std::move(P_cov__),
                                                           min_p),
            y_mean__, variance_value_new_<Tr>::adjust(std::move(y_var__), min_var), plogL__, eplogL__,
            vplogL__);
    }
    mp_macro_information_(Tr_t<M_Matrix<double>> &&P_mean__, Tr_t<M_Matrix<double>> &&P_cov__,
                         Tr_t<double> y_mean__, Tr_t<double> y_var__, Tr_t<double> plogL__,
                         Tr_t<double> eplogL__, Tr_t<double> vplogL__)
        : base_type(y_mean__,y_var__,plogL__,eplogL__,vplogL__),
          P_mean_{std::move(P_mean__)}, P_cov_{std::move(P_cov__)}{}

    mp_macro_information_(const Tr_t<M_Matrix<double>> &P_mean__,
                         const Tr_t<M_Matrix<double>> &P_cov__, Tr_t<double> y_mean__,
                         Tr_t<double> y_var__, Tr_t<double> plogL__, Tr_t<double> eplogL__,
                         Tr_t<double> vplogL__)
        :base_type(y_mean__,y_var__,plogL__,eplogL__,vplogL__),
          P_mean_{P_mean__}, P_cov_{P_cov__} {}

    mp_macro_information_() = default;
    const Tr_t<M_Matrix<double>> &P_mean() const { return P_mean_; }
    const Tr_t<M_Matrix<double>> &P_cov() const { return P_cov_; }

    template <class... Indexes> static auto get_data_index_static(Indexes... i) {
        return std::make_tuple(
            base_type::get_data_index_static(i...),
            make_data_static(
                std::make_tuple(I_s(
                    i_state{},
                    [](const self_type &self) { return self.P_mean().size(); }, 0)),
                std::make_tuple(F_s(s_Pmean{},
                                    [](const self_type &self, std::size_t i) {
                                        return self.P_mean()[i];
                                    }))),
            make_data_static(
                std::make_tuple(
                    I_s(i_state{},
                        [](const self_type &self) { return self.P_mean().size(); },
                        0),
                    I_s(CT_s<i_state, s_Transpose>{},
                        [](const self_type &self, std::size_t) {
                            return self.P_mean().size();
                        },
                        0)),
                std::make_tuple(F_s(
                    s_Pcov{}, [](const self_type &self, std::size_t i,
                                 std::size_t j) { return self.P_cov()(i, j); }))));
    }
};


class mp_state_information
{
  M_Matrix<double> P_mean_;
  M_Matrix<double> P_cov_;

  double y_mean_;
  double y_var_;
  double plogL_;
  double eplogL_;
  double vplogL_;

public:
  typedef mp_state_information self_type;
  constexpr static auto className = my_static_string("mp_state_information");
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

  static Op_void close_to_zero_test(const M_Matrix<double> &P_mean__,
                                    double tolerance) {
    std::stringstream ss;
    auto ck_mean = are_in_range<true, M_Matrix<double>>(false, 0, 1, tolerance)
                       .test(P_mean__, ss);
    if (ck_mean)
      return Op_void(true, "");
    else {
      return Op_void(false, "a probability value crossed zero: " + ss.str());
    }
  }

  static Op_void is_Binomial_Approximation_valid(double N, double p, double q,
                                                 double Np_min) {
    std::stringstream ss;
    if (!are_not_less<true, double>(std::numeric_limits<double>::epsilon())
             .test(N * p, Np_min, ss))
      return Op_void(false, "Np failed: " + ss.str());
    else if (!are_not_less<true, double>(std::numeric_limits<double>::epsilon())
                  .test(N * q, Np_min, ss))
      return Op_void(false, "Nq failed: " + ss.str());
    else {
      return Op_void(true, "");
    }
  }

  static std::pair<Op_void, Op_void>
  is_Binomial_Approximation_valid(double N, double p, double q, double Np_min,
                                  double Vp_min) {
    return {is_Binomial_Approximation_valid(N, p, q, Np_min),
            is_Binomial_Approximation_valid(N, p, q, Vp_min)};
  }

  static Op_void test(const M_Matrix<double> &P_mean__,
                      const M_Matrix<double> &P_cov__, double tolerance) {
    auto ck_mean = Probability_distribution::test<true>(P_mean__, tolerance);
    auto ck_cov =
        Probability_distribution_covariance::test<true>(P_cov__, tolerance);
    if (ck_mean.has_value() && ck_cov.has_value())
      return Op_void(true, "");
    else {
      std::stringstream ss;
      ss << " Pmean test: " << ck_mean.error()
         << " Pcov test: " << ck_cov.error();
      return Op_void(false, ss.str());
    }
  }

  static Op_void test(const M_Matrix<double> &P_mean__,
                      const M_Matrix<double> &P_cov__, double y_mean__,
                      double y_var__, double plogL__, double eplogL__,
                      double minvariance, double tolerance) {
    auto ck_mean = Probability_distribution::test<true>(P_mean__, tolerance);
    auto ck_cov =
        Probability_distribution_covariance::test<true>(P_cov__, tolerance);
    auto ck_ymean = variable_value::test<true>(y_mean__);
    auto ck_yvar = variance_value::test<true>(y_var__, minvariance, tolerance);
    auto ck_plogL = logLikelihood_value::test<true>(plogL__);
    auto ck_eplogL = logLikelihood_value::test<true>(eplogL__);
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

  static mp_state_information
  adjust(M_Matrix<double> &&P_mean__, M_Matrix<double> &&P_cov__,
         double y_mean__, double y_var__, double plogL__, double eplogL__,
         double vplogL__, double min_p, double min_var) {
    return mp_state_information(
        Probability_distribution::normalize(std::move(P_mean__), min_p),
        Probability_distribution_covariance::normalize(std::move(P_cov__),
                                                       min_p),
        y_mean__, variance_value::adjust(y_var__, min_var), plogL__, eplogL__,
        vplogL__);
  }
  mp_state_information(M_Matrix<double> &&P_mean__, M_Matrix<double> &&P_cov__,
                       double y_mean__, double y_var__, double plogL__,
                       double eplogL__, double vplogL__)
      : P_mean_{std::move(P_mean__)}, P_cov_{std::move(P_cov__)},
        y_mean_{y_mean__}, y_var_{y_var__}, plogL_{plogL__}, eplogL_{eplogL__},
        vplogL_{vplogL__} {}

  mp_state_information(const M_Matrix<double> &P_mean__,
                       const M_Matrix<double> &P_cov__, double y_mean__,
                       double y_var__, double plogL__, double eplogL__,
                       double vplogL__)
      : P_mean_{P_mean__}, P_cov_{P_cov__}, y_mean_{y_mean__}, y_var_{y_var__},
        plogL_{plogL__}, eplogL_{eplogL__}, vplogL_{vplogL__} {}

  mp_state_information() = default;
  const M_Matrix<double> &P_mean() const { return P_mean_; }
  const M_Matrix<double> &P_cov() const { return P_cov_; }

  double y_mean() const { return y_mean_; }
  double y_var() const { return y_var_; }
  double plogL() const { return plogL_; }
  double eplogL() const { return eplogL_; }
  double vplogL() const { return vplogL_; }

  template <class... Indexes> static auto get_data_index_static(Indexes...) {
    return std::make_tuple(
        make_data_static(
            std::tuple<>(),
            std::make_tuple(
                F_s(CT_s<s_pred, s_mean>{}, &self_type::y_mean),
                F_s(CT_s<s_pred, s_variance>{}, &self_type::y_var),
                F_s(CT_s<s_obs, s_logL>{}, &self_type::plogL),
                F_s(CT_s<s_pred, s_logL>{}, &self_type::eplogL),
                F_s(CT_s<s_variance, s_logL>{}, &self_type::vplogL))),
        make_data_static(
            std::make_tuple(I_s(
                i_state{},
                [](const self_type &self) { return self.P_mean().size(); }, 0)),
            std::make_tuple(F_s(s_Pmean{},
                                [](const self_type &self, std::size_t i) {
                                  return self.P_mean()[i];
                                }))),
        make_data_static(
            std::make_tuple(
                I_s(i_state{},
                    [](const self_type &self) { return self.P_mean().size(); },
                    0),
                I_s(CT_s<i_state, s_Transpose>{},
                    [](const self_type &self, std::size_t) {
                      return self.P_mean().size();
                    },
                    0)),
            std::make_tuple(F_s(
                s_Pcov{}, [](const self_type &self, std::size_t i,
                             std::size_t j) { return self.P_cov()(i, j); }))));
  }
};

class mp_single_state_information {
protected:
  M_Matrix<double> P_;
  double y_mean_;
  double y_var_;
  double plogL_;
  double eplogL_;
  double vplogL_;

public:
  typedef mp_single_state_information self_type;
  constexpr static auto className =
      my_static_string("mp_single_state_information");
  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "P", &self_type::P),
        grammar::field(C<self_type>{}, "y_mean", &self_type::y_mean),
        grammar::field(C<self_type>{}, "y_var", &self_type::y_var),
        grammar::field(C<self_type>{}, "plogL", &self_type::plogL),
        grammar::field(C<self_type>{}, "eplogL", &self_type::eplogL),
        grammar::field(C<self_type>{}, "vplogL", &self_type::vplogL));
  }

  static Op_void test(const M_Matrix<double> &P__, double y_mean__,
                      double y_var__, double plogL__, double eplogL__,
                      double minvariance, double tolerance) {
    auto ck_mean = Probability_distribution::test<true>(P__, tolerance);
    auto ck_ymean = variable_value::test<true>(y_mean__);
    auto ck_yvar = variance_value::test<true>(y_var__, minvariance, tolerance);
    auto ck_plogL = logLikelihood_value::test<true>(plogL__);
    auto ck_eplogL = logLikelihood_value::test<true>(eplogL__);
    if (ck_mean.has_value() && ck_ymean.has_value() && ck_yvar.has_value() &&
        ck_plogL.has_value() && ck_eplogL.has_value())
      return Op_void(true, "");
    else {
      std::stringstream ss;
      ss << " P test: " << ck_mean.error() << " yvar test:" << ck_yvar.error();
      ss << " pLogL test: " << ck_plogL.error() << " eplogL test"
         << ck_eplogL.error();
      return Op_void(false, ss.str());
    }
  }

  static mp_single_state_information adjust(M_Matrix<double> &&P__,
                                            double y_mean__, double y_var__,
                                            double plogL__, double eplogL__,
                                            double vplogL__, double min_p,
                                            double min_var) {
    return mp_single_state_information(
        Probability_distribution::normalize(std::move(P__), min_p), y_mean__,
        variance_value::adjust(y_var__, min_var), plogL__, eplogL__, vplogL__);
  }
  mp_single_state_information(M_Matrix<double> &&P__, double y_mean__,
                              double y_var__, double plogL__, double eplogL__,
                              double vplogL__)
      : P_{std::move(P__)}, y_mean_{y_mean__}, y_var_{y_var__}, plogL_{plogL__},
        eplogL_{eplogL__}, vplogL_{vplogL__} {}

  mp_single_state_information(const M_Matrix<double> &P__, double y_mean__,
                              double y_var__, double plogL__, double eplogL__,
                              double vplogL__)
      : P_{P__}, y_mean_{y_mean__}, y_var_{y_var__}, plogL_{plogL__},
        eplogL_{eplogL__}, vplogL_{vplogL__} {}

  mp_single_state_information() = default;
  const M_Matrix<double> &P() const { return P_; }
  M_Matrix<double> &release_P() { return P_; }

  double y_mean() const { return y_mean_; }
  double y_var() const { return y_var_; }
  double plogL() const { return plogL_; }
  double eplogL() const { return eplogL_; }
  double vplogL() const { return vplogL_; }

  template <class... Indexes> static auto get_data_index_static(Indexes...) {
    return std::make_tuple(
        make_data_static(
            std::tuple<>(),
            std::make_tuple(
                F_s(CT_s<s_pred, s_mean>{}, &self_type::y_mean),
                F_s(CT_s<s_pred, s_variance>{}, &self_type::y_var),
                F_s(CT_s<s_obs, s_logL>{}, &self_type::plogL),
                F_s(CT_s<s_pred, s_logL>{}, &self_type::eplogL),
                F_s(CT_s<s_variance, s_logL>{}, &self_type::vplogL))),
        make_data_static(
            std::make_tuple(
                I_s(i_state{},
                    [](const self_type &self) { return self.P().size(); }, 0)),
            std::make_tuple(
                F_s(s_Pmean{}, [](const self_type &self, std::size_t i) {
                  return self.P()[i];
                }))));
  }
};

template<template<class...> class Tr>
class mp_single_information_;
typedef mp_single_information_<C> mp_single_information;

template<template<class...> class Tr>
class mp_single_information_ : public mp_information_<Tr>
{
public:
    template <class T> using Tr_t=typename Tr<T>::type;
private:

protected:
    Tr_t<M_Matrix<double>> P_;

public:
    typedef mp_single_information_ self_type;
    typedef Tr_t<mp_information> base_type;
    constexpr static auto className = my_template_trait<Tr>::className+
        my_static_string("mp_single_information");
    static auto get_constructor_fields() {

        return std::tuple_cat(
std::make_tuple(
            grammar::field(C<self_type>{}, "P", &self_type::P)),
            base_type::get_constructor_fields()
            );
    }


    static mp_single_information_ adjust(Tr_t<M_Matrix<double>> &&P__,
                                              Tr_t<double> y_mean__, Tr_t<double> y_var__,
                                              Tr_t<double> plogL__, Tr_t<double> eplogL__,
                                              Tr_t<double> vplogL__, Tr_t<double> min_p,
                                              Tr_t<double> min_var) {
        return mp_single_information(
            Probability_distribution::normalize(std::move(P__), min_p), y_mean__,
            variance_value::adjust(y_var__, min_var), plogL__, eplogL__, vplogL__);
    }
    mp_single_information_(Tr_t<M_Matrix<double>> &&P__, Tr_t<double> y_mean__,
                                Tr_t<double> y_var__, Tr_t<double> plogL__, Tr_t<double> eplogL__,
                                Tr_t<double> vplogL__)
        : Tr_t<mp_information>(y_mean__,y_var__,plogL__,eplogL__,vplogL__),
          P_{std::move(P__)} {}

    mp_single_information_(const Tr_t<M_Matrix<double>> &P__, Tr_t<double> y_mean__,
                                Tr_t<double> y_var__, Tr_t<double> plogL__, Tr_t<double> eplogL__,
                                Tr_t<double> vplogL__)
        :  Tr_t<mp_information>(y_mean__,y_var__,plogL__,eplogL__,vplogL__),
            P_{P__} {}

    mp_single_information_() = default;
    const Tr_t<M_Matrix<double>> &P() const { return P_; }
    Tr_t<M_Matrix<double>> &release_P() { return P_; }


    template <class... Indexes> static auto get_data_index_static(Indexes...i) {
        return std::make_tuple(
            base_type::get_data_index_static(i...),
            make_data_static(
                std::make_tuple(
                    I_s(i_state{},
                        [](const self_type &self) { return self.P().size(); }, 0)),
                std::make_tuple(
                    F_s(s_Pmean{}, [](const self_type &self, std::size_t i) {
                        return self.P()[i];
                    }))));
    }
};

template<template<class...> class Tr>
class mp_micro_information_;
typedef mp_micro_information_<C> mp_micro_information;


template<template<class...> class Tr>
class mp_micro_information_ : public mp_single_information_<Tr> {
public:
    template <class T> using Tr_t=typename Tr<T>::type;
private:
    Tr_t<Microscopic_description_new> mi_;


public:
    typedef mp_single_information_<Tr> base_type;

    mp_micro_information_() = default;
    auto &microscopic_description() const { return mi_; }

    auto to_Macroscopic() const {
        return microscopic_description().convert_to_Macroscopic(P());
    }
    typedef mp_micro_information_ self_type;
    constexpr static auto className = my_template_trait<Tr>::className+
        my_static_string("mp_micro_information");
    static auto get_constructor_fields() {
        return std::tuple_cat(std::make_tuple(
            grammar::field(C<self_type>{}, "Microscopic_description",
                                                 &self_type::P)),
                              base_type::get_constructor_fields());}


    static mp_micro_information_
    adjust(Tr_t<Microscopic_description_new> &&mi__, Tr_t<M_Matrix<double>> &&P__,
           Tr_t<double> y_mean__, Tr_t<double> y_var__, Tr_t<double> plogL__, Tr_t<double> eplogL__,
           Tr_t<double> vplogL__, double min_p, double min_var, double reduce_by_p) {
        mi__.reduce(P__, reduce_by_p);
        //P__= mi__.add_states_to_fill_gaps_P(P__,connections_to_x,connections_from_x);
        return mp_micro_information_(
            std::move(mi__),
            Probability_distribution::normalize(std::move(P__), min_p), y_mean__,
            variance_value::adjust(y_var__, min_var), plogL__, eplogL__, vplogL__);
    }

    mp_micro_information_(Tr_t<Microscopic_description_new> &&mi__,
                          Tr_t<M_Matrix<double>> &&P__, Tr_t<double> y_mean__,
                                  Tr_t<double> y_var__, Tr_t<double> plogL__, Tr_t<double> eplogL__,
                                  Tr_t<double> vplogL__)
        : base_type(std::move(P__), y_mean__, y_var__, plogL__,
                                      eplogL__, vplogL__),
          mi_{std::move(mi__)} {}

    mp_micro_information_(const Tr_t<Microscopic_description_new> &mi__,
                          const Tr_t<M_Matrix<double>> &P__, Tr_t<double> y_mean__,
                                  Tr_t<double> y_var__, Tr_t<double> plogL__, Tr_t<double> eplogL__,
                                  Tr_t<double> vplogL__)
        : base_type(P__, y_mean__, y_var__, plogL__, eplogL__,
                                      vplogL__),
          mi_{mi__} {}

    const Tr_t<M_Matrix<double>> &P() const { return base_type::P_; }

    template <class... Indexes> static auto get_data_index_static(Indexes... i) {
        return std::make_tuple(
            base_type::get_data_index_static(i...),
            make_data_static(
                std::make_tuple(
                    I_s(i_state{},
                        [](const self_type &self) { return self.P().size(); }, 0)),
                std::make_tuple(
                    F_s(s_Pmean{}, [](const self_type &self, std::size_t i) {
                        return self.P()[i];
                    }))));
    }
};



class mp_ensamble_state_information : public mp_single_state_information {
    Microscopic_description mi_;

public:
    mp_ensamble_state_information() = default;
    auto &microscopic_description() const { return mi_; }

    auto to_Macroscopic() const {
        return microscopic_description().convert_to_Macroscopic(P());
    }
    typedef mp_ensamble_state_information self_type;
    constexpr static auto className =
        my_static_string("mp_single_state_information");
    static auto get_constructor_fields() {
        return std::make_tuple(
            grammar::field(C<self_type>{}, "Microscopic_description",
                           &self_type::P),
            grammar::field(C<self_type>{}, "P", &self_type::P),
            grammar::field(C<self_type>{}, "y_mean", &self_type::y_mean),
            grammar::field(C<self_type>{}, "y_var", &self_type::y_var),
            grammar::field(C<self_type>{}, "plogL", &self_type::plogL),
            grammar::field(C<self_type>{}, "eplogL", &self_type::eplogL),
            grammar::field(C<self_type>{}, "vplogL", &self_type::vplogL));
    }

    static Op_void test(const M_Matrix<double> &P__, double y_mean__,
                        double y_var__, double plogL__, double eplogL__,
                        double minvariance, double tolerance) {
        auto ck_mean = Probability_distribution::test<true>(P__, tolerance);
        auto ck_ymean = variable_value::test<true>(y_mean__);
        auto ck_yvar = variance_value::test<true>(y_var__, minvariance, tolerance);
        auto ck_plogL = logLikelihood_value::test<true>(plogL__);
        auto ck_eplogL = logLikelihood_value::test<true>(eplogL__);
        if (ck_mean.has_value() && ck_ymean.has_value() && ck_yvar.has_value() &&
            ck_plogL.has_value() && ck_eplogL.has_value())
            return Op_void(true, "");
        else {
            std::stringstream ss;
            ss << " P test: " << ck_mean.error() << " yvar test:" << ck_yvar.error();
            ss << " pLogL test: " << ck_plogL.error() << " eplogL test"
               << ck_eplogL.error();
            return Op_void(false, ss.str());
        }
    }

    static mp_ensamble_state_information
    adjust(Microscopic_description &&mi__, M_Matrix<double> &&P__,
           double y_mean__, double y_var__, double plogL__, double eplogL__,
           double vplogL__, double min_p, double min_var, double reduce_by_p,
           [[maybe_unused]] const std::vector<std::vector<std::size_t>> &connections_to_x ,
           [[maybe_unused]] const std::vector<std::vector<std::size_t>> &connections_from_x ) {
        mi__.reduce(P__, reduce_by_p);
        //P__= mi__.add_states_to_fill_gaps_P(P__,connections_to_x,connections_from_x);
        return mp_ensamble_state_information(
            std::move(mi__),
            Probability_distribution::normalize(std::move(P__), min_p), y_mean__,
            variance_value::adjust(y_var__, min_var), plogL__, eplogL__, vplogL__);
    }
    mp_ensamble_state_information(Microscopic_description &&mi__,
                                  M_Matrix<double> &&P__, double y_mean__,
                                  double y_var__, double plogL__, double eplogL__,
                                  double vplogL__)
        : mp_single_state_information(std::move(P__), y_mean__, y_var__, plogL__,
                                      eplogL__, vplogL__),
          mi_{std::move(mi__)} {}

    mp_ensamble_state_information(const Microscopic_description &mi__,
                                  const M_Matrix<double> &P__, double y_mean__,
                                  double y_var__, double plogL__, double eplogL__,
                                  double vplogL__)
        : mp_single_state_information(P__, y_mean__, y_var__, plogL__, eplogL__,
                                      vplogL__),
          mi_{mi__} {}

    const M_Matrix<double> &P() const { return P_; }
    double y_mean() const { return y_mean_; }
    double y_var() const { return y_var_; }
    double plogL() const { return plogL_; }
    double eplogL() const { return eplogL_; }
    double vplogL() const { return vplogL_; }

    template <class... Indexes> static auto get_data_index_static(Indexes...) {
        return std::make_tuple(
            make_data_static(
                std::tuple<>(),
                std::make_tuple(
                    F_s(CT_s<s_pred, s_mean>{}, &self_type::y_mean),
                    F_s(CT_s<s_pred, s_variance>{}, &self_type::y_var),
                    F_s(CT_s<s_obs, s_logL>{}, &self_type::plogL),
                    F_s(CT_s<s_pred, s_logL>{}, &self_type::eplogL),
                    F_s(CT_s<s_variance, s_logL>{}, &self_type::vplogL))),
            make_data_static(
                std::make_tuple(
                    I_s(i_state{},
                        [](const self_type &self) { return self.P().size(); }, 0)),
                std::make_tuple(
                    F_s(s_Pmean{}, [](const self_type &self, std::size_t i) {
                        return self.P()[i];
                    }))));
    }
};


class mp_ensamble_state_information_new : public mp_single_state_information {
    Microscopic_description_new mi_;

public:
    mp_ensamble_state_information_new() = default;
    auto &microscopic_description() const { return mi_; }

    auto to_Macroscopic() const {
        return microscopic_description().convert_to_Macroscopic(p());
    }
    typedef mp_ensamble_state_information_new self_type;
    constexpr static auto className =
        my_static_string("mp_single_state_information");
    static auto get_constructor_fields() {
        return std::make_tuple(
            grammar::field(C<self_type>{}, "Microscopic_description",
                           &self_type::p),
            grammar::field(C<self_type>{}, "P", &self_type::p),
            grammar::field(C<self_type>{}, "y_mean", &self_type::y_mean),
            grammar::field(C<self_type>{}, "y_var", &self_type::y_var),
            grammar::field(C<self_type>{}, "plogL", &self_type::plogL),
            grammar::field(C<self_type>{}, "eplogL", &self_type::eplogL),
            grammar::field(C<self_type>{}, "vplogL", &self_type::vplogL));
    }

    static Op_void test(const M_Matrix<double> &P__, double y_mean__,
                        double y_var__, double plogL__, double eplogL__,
                        double minvariance, double tolerance) {
        auto ck_mean = Probability_distribution::test<true>(P__, tolerance);
        auto ck_ymean = variable_value::test<true>(y_mean__);
        auto ck_yvar = variance_value::test<true>(y_var__, minvariance, tolerance);
        auto ck_plogL = logLikelihood_value::test<true>(plogL__);
        auto ck_eplogL = logLikelihood_value::test<true>(eplogL__);
        if (ck_mean.has_value() && ck_ymean.has_value() && ck_yvar.has_value() &&
            ck_plogL.has_value() && ck_eplogL.has_value())
            return Op_void(true, "");
        else {
            std::stringstream ss;
            ss << " P test: " << ck_mean.error() << " yvar test:" << ck_yvar.error();
            ss << " pLogL test: " << ck_plogL.error() << " eplogL test"
               << ck_eplogL.error();
            return Op_void(false, ss.str());
        }
    }

    static mp_ensamble_state_information_new
    adjust(Microscopic_description_new mi__, M_Matrix<double> &&P__,
           double y_mean__, double y_var__, double plogL__, double eplogL__,
           double vplogL__, double min_p, double min_var, double reduce_by_p) {
        reduce_P_mi(mi__,P__, reduce_by_p);
        return mp_ensamble_state_information_new(
            std::move(mi__),
            Probability_distribution::normalize(std::move(P__), min_p), y_mean__,
            variance_value::adjust(y_var__, min_var), plogL__, eplogL__, vplogL__);
    }
    mp_ensamble_state_information_new(Microscopic_description_new &&mi__,
                                      M_Matrix<double> &&P__, double y_mean__,
                                      double y_var__, double plogL__, double eplogL__,
                                      double vplogL__)
        : mp_single_state_information(std::move(P__), y_mean__, y_var__, plogL__,
                                      eplogL__, vplogL__),
          mi_{std::move(mi__)} {}

    mp_ensamble_state_information_new(const Microscopic_description_new &mi__,
                                      const M_Matrix<double> &P__, double y_mean__,
                                      double y_var__, double plogL__, double eplogL__,
                                      double vplogL__)
        : mp_single_state_information(P__, y_mean__, y_var__, plogL__, eplogL__,
                                      vplogL__),
          mi_{mi__} {}

    const M_Matrix<double> &p() const { return P_; }
    double y_mean() const { return y_mean_; }
    double y_var() const { return y_var_; }
    double plogL() const { return plogL_; }
    double eplogL() const { return eplogL_; }
    double vplogL() const { return vplogL_; }

    template <class... Indexes> static auto get_data_index_static(Indexes...) {
        return std::make_tuple(
            make_data_static(
                std::tuple<>(),
                std::make_tuple(
                    F_s(CT_s<s_pred, s_mean>{}, &self_type::y_mean),
                    F_s(CT_s<s_pred, s_variance>{}, &self_type::y_var),
                    F_s(CT_s<s_obs, s_logL>{}, &self_type::plogL),
                    F_s(CT_s<s_pred, s_logL>{}, &self_type::eplogL),
                    F_s(CT_s<s_variance, s_logL>{}, &self_type::vplogL))),
            make_data_static(
                std::make_tuple(
                    I_s(i_state{},
                        [](const self_type &self) { return self.p().size(); }, 0)),
                std::make_tuple(
                    F_s(s_Pmean{}, [](const self_type &self, std::size_t i) {
                        return self.p()[i];
                    }))));
    }
};







}
// namespace markov
namespace markov {

enum MACROR_new { MACRO_DMNR_new, MACRO_DVNR_new, MACRO_DMR_new, MACRO_DVR_new };

const std::string MACROR_string_new[] = {"MacroDMNR_new", "MacroDVNR_new", "MacroDMR_new",
                                     "MacroDVR_new"};

constexpr MACROR_new MACROR_VALUES_new[] = {MACRO_DMNR_new, MACRO_DVNR_new, MACRO_DMR_new,
                                    MACRO_DVR_new};

inline std::ostream &operator<<(std::ostream &s, const MACROR_new &m) {
    s << MACROR_string_new[m];
    return s;
}

inline std::istream &operator>>(std::istream &is, MACROR_new &m) {
    std::string s;
    is >> s;
    for (auto v : MACROR_VALUES_new) {
        if (s == MACROR_string_new[v])
            m = v;
        return is;
    }
    is.setstate(std::ios::failbit);
    return is;
}



enum MACROR { MACRO_DMNR, MACRO_DVNR, MACRO_DMR, MACRO_DVR };

const std::string MACROR_string[] = {"MacroDMNR", "MacroDVNR", "MacroDMR",
                                     "MacroDVR"};

constexpr MACROR MACROR_VALUES[] = {MACRO_DMNR, MACRO_DVNR, MACRO_DMR,
                                    MACRO_DVR};

inline std::ostream &operator<<(std::ostream &s, const MACROR &m) {
  s << MACROR_string[m];
  return s;
}

inline std::istream &operator>>(std::istream &is, MACROR &m) {
  std::string s;
  is >> s;
  for (auto v : MACROR_VALUES) {
    if (s == MACROR_string[v])
      m = v;
    return is;
  }
  is.setstate(std::ios::failbit);
  return is;
}

enum MICROR { MICRO_DMR, MICRO_DVR };

const std::string MICROR_string[] = {"MicroDMR", "MicroDVR"};

constexpr MICROR MICROR_VALUES[] = {MICRO_DMR, MICRO_DVR};

inline std::ostream &operator<<(std::ostream &s, const MICROR &m) {
  s << MICROR_string[m];
  return s;
}

inline std::istream &operator>>(std::istream &is, MICROR &m) {
  std::string s;
  is >> s;
  for (auto v : MICROR_VALUES) {
    if (s == MICROR_string[v])
      m = v;
    return is;
  }
  is.setstate(std::ios::failbit);
  return is;
}

//enum RMICROR { RMICRO_DMR, RMICRO_DVR };

//const std::string RMICROR_string[] = {"RMicroDMR", "RMicroDVR"};

//constexpr RMICROR RMICROR_VALUES[] = {RMICRO_DMR, RMICRO_DVR};

//inline std::ostream &operator<<(std::ostream &s, const RMICROR &m) {
//  s << RMICROR_string[m];
//  return s;
//}

//inline std::istream &operator>>(std::istream &is, RMICROR &m) {
//  std::string s;
//  is >> s;
//  for (auto v : RMICROR_VALUES) {
//    if (s == RMICROR_string[v])
//      m = v;
//    return is;
//  }
//  is.setstate(std::ios::failbit);
//  return is;
//}

enum RMICROR_new { RMICRO_DMR_new, RMICRO_DVR_new };

const std::string RMICROR_new_string[] = {"RMicroDMR_new", "RMicroDVR_new"};

constexpr RMICROR_new RMICROR_new_VALUES[] = {RMICRO_DMR_new, RMICRO_DVR_new};

inline std::ostream &operator<<(std::ostream &s, const RMICROR_new &m) {
    s << RMICROR_new_string[m];
    return s;
}

inline std::istream &operator>>(std::istream &is, RMICROR_new &m) {
    std::string s;
    is >> s;
    for (auto v : RMICROR_new_VALUES) {
        if (s == RMICROR_new_string[v])
            m = v;
        return is;
    }
    is.setstate(std::ios::failbit);
    return is;
}

} // namespace markov

// namespace markov
template <> struct my_trait<markov::MACROR> {
  constexpr static auto className = my_static_string("Macror_Algorithm");

  template <class DataFrame>
  static void insert_col(DataFrame &d, const std::string &pre) {
    for (auto e : markov::MACROR_VALUES)
      d.insert_column(pre + markov::MACROR_string[e].c_str(), C<double>{});
  }

  static auto data_row(markov::MACROR d) {
    switch (d) {
    case markov::MACRO_DMNR:
      return std::tuple(1.0, 0.0, 0.0, 0.0);
    case markov::MACRO_DVNR:
      return std::tuple(0.0, 1.0, 0.0, 0.0);
    case markov::MACRO_DMR:
      return std::tuple(0.0, 0.0, 1.0, 0.0);
    case markov::MACRO_DVR:
    default:
      return std::tuple(0.0, 0.0, 0.0, 1.0);
    }
  }
};

template <> struct my_trait<markov::MACROR_new> {
    constexpr static auto className = my_static_string("Macror_Algorithm_new");

    template <class DataFrame>
    static void insert_col(DataFrame &d, const std::string &pre) {
        for (auto e : markov::MACROR_VALUES_new)
            d.insert_column(pre + markov::MACROR_string_new[e].c_str(), C<double>{});
    }

    static auto data_row(markov::MACROR d) {
        switch (d) {
        case markov::MACRO_DMNR:
            return std::tuple(1.0, 0.0, 0.0, 0.0);
        case markov::MACRO_DVNR:
            return std::tuple(0.0, 1.0, 0.0, 0.0);
        case markov::MACRO_DMR:
            return std::tuple(0.0, 0.0, 1.0, 0.0);
        case markov::MACRO_DVR:
        default:
            return std::tuple(0.0, 0.0, 0.0, 1.0);
        }
    }
};


template <> struct my_trait<markov::MICROR> {
    constexpr static auto className = my_static_string("Macror_Algorithm");

    template <class DataFrame>
    static void insert_col(DataFrame &d, const std::string &pre) {
        for (auto e : markov::MICROR_VALUES)
            d.insert_column(pre + markov::MICROR_string[e].c_str(), C<double>{});
    }

    static auto data_row(markov::MICROR d) {
        switch (d) {
        case markov::MICRO_DMR:
            return std::tuple(1.0, 0.0);
        case markov::MICRO_DVR:
        default:
            return std::tuple( 0.0, 1.0);
        }
    }
};

//template <> struct my_trait<markov::RMICROR> {
//    constexpr static auto className = my_static_string("Macror_Algorithm");

//    template <class DataFrame>
//    static void insert_col(DataFrame &d, const std::string &pre) {
//        for (auto e : markov::RMICROR_VALUES)
//            d.insert_column(pre + markov::RMICROR_string[e].c_str(), C<double>{});
//    }

//    static auto data_row(markov::RMICROR d) {
//        switch (d) {
//        case markov::RMICRO_DMR:
//            return std::tuple(1.0, 0.0);
//        case markov::RMICRO_DVR:
//        default:
//            return std::tuple( 0.0, 1.0);
//        }
//    }
//};

template <> struct my_trait<markov::RMICROR_new> {
    constexpr static auto className = my_static_string("Macror_Algorithm");

    template <class DataFrame>
    static void insert_col(DataFrame &d, const std::string &pre) {
        for (auto e : markov::RMICROR_new_VALUES)
            d.insert_column(pre + markov::RMICROR_new_string[e].c_str(), C<double>{});
    }

    static auto data_row(markov::RMICROR_new d) {
        switch (d) {
        case markov::RMICRO_DMR_new:
            return std::tuple(1.0, 0.0);
        case markov::RMICRO_DVR_new:
        default:
            return std::tuple( 0.0, 1.0);
        }
    }
};


template <> struct myData_Index<markov::MACROR> {

  typedef markov::MACROR self_type;
  static Data_Index_scheme data_index() {
    Data_Index_scheme out;
    for (auto e : markov::MACROR_VALUES) {
      out.push_back(markov::MACROR_string[e].c_str(), {});
    }
    return out;
  }

  template <class... Indexes> static auto get_data_index_static(Indexes...) {
    return std::make_tuple(make_data_static(
        std::tuple<>(), std::make_tuple(F_s(markov::s_MacroAlgorithm{}, [
        ](self_type s) -> auto const & { return markov::MACROR_string[s]; }))));
  }
};

template <> class moments<markov::MACROR> {
  static auto calc(const std::vector<markov::MACROR> &d) {
    std::vector<moments<double>> data(size);
    for (auto e : d)
      data[e].push_back(1.0);
    return data;
  }

  template <class C, class F, class... Ts>
  static auto calc(const std::vector<C> &d, const F &f, Ts... t) {
    std::vector<moments<double>> data(size);
    for (auto e : d) {
      markov::MACROR m = std::invoke(f, e, t...);
      for (auto mi : markov::MACROR_VALUES)
        if (mi == m)
          data[mi].push_back(1.0);
        else
          data[mi].push_back(0.0);
    }
    return data;
  }

  template <std::size_t... Is>
  auto data_row(MOMENTS m, std::index_sequence<Is...>) const {

    return std::tuple_cat(data_[Is].data_row(m)...);
  }

public:
  typedef moments self_type;
  constexpr static auto className =
      my_static_string("moments_") + my_trait<markov::MACROR>::className;
  static auto get_constructor_fields() {
    return std::tuple(
        grammar::field(C<self_type>{}, "values", &self_type::getMoments));
  }

  template <class C, class F, class... Ts>
  moments(const std::vector<C> &data, const F &f, const Ts... t)
      : data_(calc(data, f, t...)) {}
  moments(const std::vector<markov::MACROR> &data) : data_{calc(data)} {}

  moments(const std::vector<moments<double>> &data) : data_{data} {}

  std::size_t data_size() const { return size; }

  auto data_row(MOMENTS m) const {
    return data_row(m, std::make_index_sequence<size>());
  }

  std::vector<moments<double>> const &getMoments() const { return data_; }

  moments() = default;

  moments &operator+=(const moments &other) {
    for (std::size_t i = 0; i < data_.size(); ++i)
      data_[i] += other.getMoments()[i];
    return *this;
  }

  constexpr static const std::size_t size =
      sizeof(markov::MACROR_VALUES) / sizeof(markov::MACROR_VALUES[0]);

private:
  std::vector<moments<double>> data_;
};

namespace markov {

namespace hidden {
class MacroDMNR {
public:
  static constexpr MACROR self_kind = MACRO_DMNR;

  inline constexpr static auto const className = my_static_string("MacroDMNR");
  template <class Model, class Step>
  myOptional_t<mp_state_information> run(const mp_state_information &prior,
                                         Model &m, const Step &p,
                                         std::ostream &os) const {
    MACROR alg;
    return run(prior, m, p, os, alg);
  }

  template <class Model, class Step, class MACROR>
  myOptional_t<mp_state_information> run(const mp_state_information &prior,
                                         Model &m, const Step &p,
                                         std::ostream &os, MACROR &alg) const {
    typedef myOptional_t<mp_state_information> Op;
    auto Q_dto = m.get_P(p, 0);
    if (!Q_dto)
      return Op(false, "fails in auto Q_dt=m.get_P(p,0) :" + Q_dto.error());
    Markov_Transition_step_double Q_dt = std::move(Q_dto).value();
    return run(prior, Q_dt, m, p, os, alg);
  }

  template <class Model, class Step>
  myOptional_t<mp_state_information>
  run(const mp_state_information &prior,
      const Markov_Transition_step_double& Q_dt,  Model &m, const Step &p,
      std::ostream &os, const MACROR &alg) const {
    typedef myOptional_t<mp_state_information> Op;
    switch (alg) {
    case MACRO_DMNR:
      return MacroDMNR(tolerance()).run(prior, Q_dt, m, p, os);
    case MACRO_DVNR:
      return Op(false,
                "Not a contemplated algorithm :" + MACROR_string[MACRO_DMNR] +
                    " is not valid for " + className.str());
    case MACRO_DMR:
      return Op(false,
                "Not a contemplated algorithm :" + MACROR_string[MACRO_DVNR] +
                    " is not valid for " + className.str());
    case MACRO_DVR:
    default:
      return Op(false,
                "Not a contemplated algorithm :" + MACROR_string[MACRO_DVR] +
                    " is not valid for " + className.str());
    }
  }

  template <class Model, class Step>
  myOptional_t<mp_state_information>
  run(const mp_state_information &prior,
      const Markov_Transition_step_double &Q_dt, Model &m, const Step &p,
      std::ostream &os, MACROR &alg) const {
    alg = MACRO_DMNR;
    const MACROR alg2 = alg;
    return run(prior, Q_dt, m, p, os, alg2);
  }

  template <class Model, class Step, class... Aux>
  myOptional_t<mp_state_information>
  run(const mp_state_information &prior, Markov_Transition_step_double Q_dt,
      Model &m, const Step &p, std::ostream &) const {
    typedef myOptional_t<mp_state_information> Op;
    auto y = p.y();
    double e = m.noise_variance(p.nsamples());
    double N = m.AverageNumberOfChannels();
    M_Matrix<double> u(prior.P_mean().size(), 1, 1.0);

    auto SmD = prior.P_cov() - diag(prior.P_mean());
    double gSg =
        (TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.gmean_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_ij(), Q_dt.gmean_ij()) * u))
            .getvalue();

    double ms = (prior.P_mean() * Q_dt.gvar_i()).getvalue();

    auto e_mu = e + N * ms;
    auto y_mean = N * (prior.P_mean() * Q_dt.gmean_i()).getvalue();
    auto y_var = e_mu + N * gSg;
    if (std::isnan(y)) {
      double vplogL = 0.0;
      double plogL = std::numeric_limits<double>::quiet_NaN();
      double eplogL = std::numeric_limits<double>::quiet_NaN();
      auto P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P());
      auto P_mean = prior.P_mean() * Q_dt.P();
      P__cov += diag(P_mean);
      // std::cerr<<"\nPcov nana corr\n"<<P__cov<<"\nP_mean nana
      // corr\n"<<P_mean<<"\nQ.P \n"<<Q_dt.P();
      auto test = mp_state_information::test(P_mean, P__cov, tolerance_);
      if (test.has_value())
        return Op(mp_state_information::adjust(
            std::move(P_mean), std::move(P__cov), y_mean, y_var, plogL, eplogL,
            vplogL, Q_dt.min_P(), e));
      else
        return Op(false, "fails at intertrace prediction!!: " + test.error());
    }
    auto dy = y - y_mean;
    auto chi = dy / y_var;
    auto P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P());
    auto P_mean = prior.P_mean() * Q_dt.P();
    P__cov += diag(P_mean);

    auto chi2 = dy * chi;

    double plogL;
    if (y_var > 0)
      plogL = -0.5 * log(2 * PI * y_var) - 0.5 * chi2;
    else
      plogL = std::numeric_limits<double>::infinity();

    double eplogL =
        -0.5 * log(2 * PI * y_var) - 0.5; // e_mu+N*gSg"-N*zeta*sqr(sSg)"
    // double chilogL=(eplogL-plogL)/std::sqrt(0.5);

    double vplogL = 0.5;

    auto test = mp_state_information::test(P_mean, P__cov, y_mean, y_var, plogL,
                                           eplogL, e, tolerance());
    if (!test) {
      std::stringstream ss;

      ss << "\nP_mean \n" << P_mean;
      ss << "\nPcov \n" << P__cov;
      // ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;

      return Op(false, "\nfails in trace!!!; error=" + test.error() + ss.str());
    } else
      return Op(mp_state_information::adjust(
          std::move(P_mean), std::move(P__cov), y_mean, y_var, plogL, eplogL,
          vplogL, Q_dt.min_P(), e));
  }

  template <class Model, class Step>
  myOptional_t<mp_state_information> start(Model &m, const Step &p,
                                           double min_p) const {
    typedef myOptional_t<mp_state_information> Op;

    auto P_mean = m.Peq(p.begin()->x());
    //  std::cerr<<"gslreijgsorjgps INIT!!!"<<P_mean.value();
    if (!P_mean)
      return Op(false, "fails to get Peq :" + P_mean.error());
    else {
      auto P_cov = diag(P_mean.value()) - quadraticForm_XTX(P_mean.value());
      double nan = std::numeric_limits<double>::quiet_NaN();
      return Op(mp_state_information::adjust(std::move(P_mean).value(),
                                             std::move(P_cov), nan, nan, nan,
                                             nan, min_p, 0));
    }
  }
  double tolerance() const { return tolerance_; }

  MacroDMNR(double tolerance) : tolerance_{tolerance} {}
  MacroDMNR() = default;

private:
  double tolerance_ = 1e-5;
};

class MacroDVNR {
public:
  inline constexpr static auto const className = my_static_string("MacroDVNR");

  template <class Model, class Step>
  myOptional_t<mp_state_information> run(const mp_state_information &prior,
                                         Model &m, const Step &p,
                                         std::ostream &os) const {
    MACROR alg;
    return run(prior, m, p, os, alg);
  }

  template <class Model, class Step, class MACROR>
  myOptional_t<mp_state_information> run(const mp_state_information &prior,
                                         Model &m, const Step &p,
                                         std::ostream &os, MACROR &alg) const {
    typedef myOptional_t<mp_state_information> Op;
    auto Q_dto = m.get_P(p, 0);
    if (!Q_dto)
      return Op(false, "fails in auto Q_dt=m.get_P(p,0) :" + Q_dto.error());
    Markov_Transition_step_double Q_dt = std::move(Q_dto).value();
    return run(prior, Q_dt, m, p, os, alg);
  }

  template <class Model, class Step>
  myOptional_t<mp_state_information>
  run(const mp_state_information &prior,
      const Markov_Transition_step_double& Q_dt, Model &m, const Step &p,
      std::ostream &os, const MACROR &alg) const {
    typedef myOptional_t<mp_state_information> Op;
    switch (alg) {
    case MACRO_DMNR:
      return MacroDMNR(tolerance()).run(prior, Q_dt, m, p, os);
    case MACRO_DVNR:
      return MacroDVNR(tolerance(), Variance_magic_number())
          .run(prior, Q_dt, m, p, os);
    case MACRO_DMR:
      return Op(false,
                "Not a contemplated algorithm :" + MACROR_string[MACRO_DMR] +
                    " is not valid for " + className.str());
    case MACRO_DVR:
    default:
      return Op(false,
                "Not a contemplated algorithm :" + MACROR_string[MACRO_DVR] +
                    " is not valid for " + className.str());
    }
  }

  template <class Model, class Step>
  myOptional_t<mp_state_information>
  run(const mp_state_information &prior,
      const Markov_Transition_step_double &Q_dt, Model &m, const Step &p,
      std::ostream &os, MACROR &alg) const {
    auto y = p.y();
    if (std::isnan(y))
      alg = MACRO_DMNR;
    else {
      double mg = (prior.P_mean() * Q_dt.gmean_i()).getvalue();
      double g_max = max(Q_dt.g());
      double g_min = min(Q_dt.g());
      double g_range = g_max - g_min;
      double N = m.AverageNumberOfChannels();
      auto p_bi = (g_max - mg) / g_range;
      auto q_bi = (mg - g_min) / g_range;
      auto test_Binomial =
          mp_state_information::is_Binomial_Approximation_valid(
              N, p_bi, q_bi, Variance_magic_number());
      if (!test_Binomial.has_value())
        alg = MACRO_DMNR;
      else
        alg = MACRO_DVNR;
    }
    const MACROR alg2 = alg;
    return run(prior, Q_dt, m, p, os, alg2);
  }

  template <class Model, class Step>
  myOptional_t<mp_state_information>
  run(const mp_state_information &prior,
      const Markov_Transition_step_double &Q_dt, Model &m, const Step &p,
      std::ostream & /*os*/) const {
    typedef myOptional_t<mp_state_information> Op;
    auto y = p.y();

    double N = m.AverageNumberOfChannels();

    double e = m.noise_variance(p.nsamples());
    M_Matrix<double> u(prior.P_mean().size(), 1, 1.0);

    auto SmD = prior.P_cov() - diag(prior.P_mean());
    double ms = (prior.P_mean() * Q_dt.gvar_i()).getvalue();

    double gSg =
        (TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.gmean_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_ij(), Q_dt.gmean_ij()) * u))
            .getvalue();

    double sSg =
        (TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.gmean_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_var_ij(), Q_dt.gmean_ij()) * u))
            .getvalue();
    double sSs =
        (TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.gvar_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_var_ij(), Q_dt.gvar_ij()) * u))
            .getvalue();
    auto sS = TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.P() +
              prior.P_mean() * Q_dt.gtotal_var_ij();
    auto gS = TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.P() +
              prior.P_mean() * Q_dt.gtotal_ij();

    double delta_emu = std::max(sqr(ms + e / N) - 2.0 / N * sSs, 0.0);
    double ms0 = (ms - e / N) / 2 + std::sqrt(delta_emu) / 2;

    auto e_mu = e + N * ms0;
    auto y_mean =
        N * (prior.P_mean() * Q_dt.gmean_i()).getvalue() - N * 0.5 / e_mu * sSg;
    auto zeta = N / (2 * sqr(e_mu) + N * sSs);
    auto y_var = std::max(e_mu + N * gSg - N * zeta * sqr(sSg), e);

    auto dy = y - y_mean;
    auto chi = dy / y_var;

    auto chi2 = dy * chi;

    double plogL;
    if (y_var > 0)
      plogL = -0.5 * log(2 * PI * y_var) - 0.5 * chi2;
    else
      plogL = std::numeric_limits<double>::infinity();

    double eplogL =
        -0.5 * log(2 * PI * y_var) - 0.5; // e_mu+N*gSg"-N*zeta*sqr(sSg)"
    // double chilogL=(eplogL-plogL)/std::sqrt(0.5);

    double vplogL = 0.5;
    //      std::cerr<<" \n\n----test----\n"<<test.error()<<"\n";
    auto P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P());
    auto P_mean = prior.P_mean() * Q_dt.P();
    P__cov += diag(P_mean);
    auto test = mp_state_information::test(P_mean, P__cov, y_mean, y_var, plogL,
                                           eplogL, e, tolerance());
    if (!test) {
      std::stringstream ss;
      ss << "\n step=" << p << " "
         << " y=" << y << " ymean=" << y_mean << " yvar=" << y_var << " e=" << e
         << " e_mu" << e_mu << " N*gSg" << N * gSg
         << " -N*zeta*sSg=" << -N * zeta * sqr(sSg);
      ss << "\n ms0=(ms-e/N)/2+std::sqrt(delta_emu)/2;\nms0=" << ms0
         << "\n ms=" << ms << "\n e/N/2=" << e / N;

      ss << "\n std::sqrt(delta_emu)/2=" << std::sqrt(delta_emu) / 2;
      ss << "\ndelta_emu=sqr(ms+e/N)-2.0/N*sSs \nsqr(ms+e/N)="
         << sqr(ms + e / N) << " 2.0/N*sSs=" << 2.0 / N * sSs << " sSs=" << sSs;
      ss << "\nP_mean \n" << P_mean;
      ss << " -N*zeta*sqr(sSg)=" << -N * zeta * sqr(sSg) << " plogL=" << plogL
         << " eplogL=" << eplogL << "dif=" << plogL - eplogL << " sSs=" << sSs;
      ss << "\nPcov \n" << P__cov;
      // ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;

      return Op(false, "\nfails in trace!!!; error=" + test.error() + ss.str());
    } else
      return Op(mp_state_information::adjust(
          std::move(P_mean), std::move(P__cov), y_mean, y_var, plogL, eplogL,
          vplogL, Q_dt.min_P(), e));
  }

  template <class Model, class Step>
  myOptional_t<mp_state_information> start(Model &m, const Step &p,
                                           double min_p) const {
    typedef myOptional_t<mp_state_information> Op;

    auto P_mean = m.Peq(p.begin()->x());
    //  std::cerr<<"gslreijgsorjgps INIT!!!"<<P_mean.value();
    if (!P_mean)
      return Op(false, "fails to get Peq :" + P_mean.error());
    else {
      auto P_cov = diag(P_mean.value()) - quadraticForm_XTX(P_mean.value());
      //     std::cerr<<"gslreijgsorjgps INIT!!!"<<P_cov;
      double nan = std::numeric_limits<double>::quiet_NaN();
      return Op(mp_state_information::adjust(std::move(P_mean).value(),
                                             std::move(P_cov), nan, nan, nan,
                                             nan, min_p, 0));
    }
  }
  double tolerance() const { return tolerance_; }

  MacroDVNR() = default;
  double Variance_magic_number() const { return Variance_magic; }
  MacroDVNR(double tolerance, double variance_magical)
      : tolerance_{tolerance}, Variance_magic(variance_magical) {}

private:
  double tolerance_ = 1e-5;
  double Variance_magic = 1;
};
class MacroDMR {
public:
  inline constexpr static auto const className = my_static_string("MacroDMR");

  template <class Model, class Step>
  myOptional_t<mp_state_information> run(const mp_state_information &prior,
                                         Model &m, const Step &p,
                                         std::ostream &os) const {
    MACROR alg;
    return run(prior, m, p, os, alg);
  }

  template <class Model, class Step, class MACROR>
  myOptional_t<mp_state_information> run(const mp_state_information &prior,
                                         Model &m, const Step &p,
                                         std::ostream &os, MACROR &alg) const {
    typedef myOptional_t<mp_state_information> Op;
    auto Q_dto = m.get_P(p, 0);
    if (!Q_dto)
      return Op(false, "fails in auto Q_dt=m.get_P(p,0) :" + Q_dto.error());
    Markov_Transition_step_double Q_dt = std::move(Q_dto).value();
    return run(prior, Q_dt, m, p, os, alg);
  }

  template <class Model, class Step>
  myOptional_t<mp_state_information>
  run(const mp_state_information &prior,
      const Markov_Transition_step_double& Q_dt, Model &m, const Step &p,
      std::ostream &os, const MACROR &alg) const {
    typedef myOptional_t<mp_state_information> Op;
    switch (alg) {
    case MACRO_DMNR:
      return MacroDMNR(tolerance()).run(prior, Q_dt, m, p, os);
    case MACRO_DVNR:
      return Op(false,
                "Not a contemplated algorithm :" + MACROR_string[MACRO_DVNR] +
                    " is not valid for " + className.str());
    case MACRO_DMR:
      return MacroDMR(tolerance(), Binomial_magic_number())
          .run(prior, Q_dt, m, p, os);
    case MACRO_DVR:
    default:
      return Op(false,
                "Not a contemplated algorithm :" + MACROR_string[MACRO_DVR] +
                    " is not valid for " + className.str());
    }
  }

  template <class Model, class Step>
  myOptional_t<mp_state_information>
  run(const mp_state_information &prior,
      const Markov_Transition_step_double &Q_dt, Model &m, const Step &p,
      std::ostream &os, MACROR &alg) const {
    auto y = p.y();
    if (std::isnan(y))
      alg = MACRO_DMNR;
    else {
      double mg = (prior.P_mean() * Q_dt.gmean_i()).getvalue();
      double g_max = max(Q_dt.g());
      double g_min = min(Q_dt.g());
      double g_range = g_max - g_min;
      double N = m.AverageNumberOfChannels();
      auto p_bi = (g_max - mg) / g_range;
      auto q_bi = (mg - g_min) / g_range;
      auto test_Binomial =
          mp_state_information::is_Binomial_Approximation_valid(
              N, p_bi, q_bi, Binomial_magic_number());
      if (!test_Binomial.has_value())
        alg = MACRO_DMNR;
      else
        alg = MACRO_DMR;
    }
    const MACROR alg2 = alg;
    return run(prior, Q_dt, m, p, os, alg2);
  }

  template <class Model, class Step>
  myOptional_t<mp_state_information>
  run(const mp_state_information &prior,
      const Markov_Transition_step_double &Q_dt, Model &m, const Step &p,
      std::ostream & /*os*/) const {
    typedef myOptional_t<mp_state_information> Op;
    auto y = p.y();

    // double mg=(prior.P_mean()*Q_dt.gmean_i()).getvalue();
    // double g_max=max(Q_dt.g());
    // double g_min=min(Q_dt.g());
    //     double g_range=g_max-g_min;
    double N = m.AverageNumberOfChannels();
    double e = m.noise_variance(p.nsamples());
    M_Matrix<double> u(prior.P_mean().size(), 1, 1.0);
    auto SmD = prior.P_cov() - diag(prior.P_mean());

    double gSg =
        (TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.gmean_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_ij(), Q_dt.gmean_ij()) * u))
            .getvalue();

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

    double eplogL =
        -0.5 * log(2 * PI * y_var) - 0.5; // e_mu+N*gSg"-N*zeta*sqr(sSg)"
    // double chilogL=(eplogL-plogL)/std::sqrt(0.5);

    auto test = mp_state_information::test(P_mean, P__cov, y_mean, y_var, plogL,
                                           eplogL, e, tolerance());
    if (!test) {
      std::stringstream ss;

      ss << "\nP_mean \n" << P_mean;
      ss << "\nPcov \n" << P__cov;
      // ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;

      return Op(false, "\nfails in trace!!!; error=" + test.error() + ss.str());
    } else
      return Op(mp_state_information::adjust(std::move(P_mean),
                                             std::move(P__cov), y_mean, y_var,
                                             plogL, eplogL, Q_dt.min_P(), e));
  }

  template <class Model, class Step>
  myOptional_t<mp_state_information> start(Model &m, const Step &p,
                                           double min_p) const {
    typedef myOptional_t<mp_state_information> Op;

    auto P_mean = m.Peq(p.begin()->x());
    //  std::cerr<<"gslreijgsorjgps INIT!!!"<<P_mean.value();
    if (!P_mean)
      return Op(false, "fails to get Peq :" + P_mean.error());
    else {
      auto P_cov = diag(P_mean.value()) - quadraticForm_XTX(P_mean.value());
      //     std::cerr<<"gslreijgsorjgps INIT!!!"<<P_cov;
      double nan = std::numeric_limits<double>::quiet_NaN();
      return Op(mp_state_information::adjust(std::move(P_mean).value(),
                                             std::move(P_cov), nan, nan, nan,
                                             nan, min_p, 0));
    }
  }
  double tolerance() const { return tolerance_; }
  double Binomial_magic_number() const { return binomial_magic; }
  MacroDMR(double tolerance, double binomial_magical)
      : tolerance_{tolerance}, binomial_magic(binomial_magical) {}
  MacroDMR() = default;

  template <class Model, class Step>
  myOptional_t<mp_state_information> start(Derivative<Model> &m, const Step &p,
                                           double min_p) const {
    typedef myOptional_t<mp_state_information> Op;

    auto P_mean = m.Peq(p.begin()->x());
    //  std::cerr<<"gslreijgsorjgps INIT!!!"<<P_mean.value();
    if (!P_mean)
      return Op(false, "fails to get Peq :" + P_mean.error());
    else {
      auto P_cov = diag(P_mean.value()) - quadraticForm_XTX(P_mean.value());
      //     std::cerr<<"gslreijgsorjgps INIT!!!"<<P_cov;
      double nan = std::numeric_limits<double>::quiet_NaN();
      return Op(mp_state_information::adjust(std::move(P_mean).value(),
                                             std::move(P_cov), nan, nan, nan,
                                             nan, min_p, 0));
    }
  }

private:
  double tolerance_ = 1e-5;
  double binomial_magic = 5;
};

class MacroDVR {
public:
  inline constexpr static auto const className = my_static_string("MacroDVR");

  template <class Model, class Step>
  myOptional_t<mp_state_information> run(const mp_state_information &prior,
                                         Model &m, const Step &p,
                                         std::ostream &os) const {
    MACROR alg;
    return run(prior, m, p, os, alg);
  }

  template <class Model, class Step, class MACROR>
  myOptional_t<mp_state_information> run(const mp_state_information &prior,
                                         Model &m, const Step &p,
                                         std::ostream &os, MACROR &alg) const {
    typedef myOptional_t<mp_state_information> Op;
    auto Q_dto = m.get_P(p, 0);
    if (!Q_dto)
      return Op(false, "fails in auto Q_dt=m.get_P(p,0) :" + Q_dto.error());
    Markov_Transition_step_double Q_dt = std::move(Q_dto).value();
    return run(prior, Q_dt, m, p, os, alg);
  }

  template <class Model, class Step>
  myOptional_t<mp_state_information>
  run(const mp_state_information &prior,
      const Markov_Transition_step_double& Q_dt, Model &m, const Step &p,
      std::ostream &os, const MACROR &alg) const {
    switch (alg) {
    case MACRO_DMNR:
      return MacroDMNR(tolerance()).run(prior, Q_dt, m, p, os);
    case MACRO_DVNR:
      return MacroDVNR(tolerance(), Variance_magic_number())
          .run(prior, Q_dt, m, p, os);
    case MACRO_DMR:
      return MacroDMR(tolerance(), Binomial_magic_number())
          .run(prior, Q_dt, m, p, os);
    case MACRO_DVR:
    default:
      return MacroDVR(tolerance(), Binomial_magic_number(),
                      Variance_magic_number())
          .run(prior, Q_dt, m, p, os);
    }
  }

  template <class Model, class Step>
  myOptional_t<mp_state_information>
  run(const mp_state_information &prior,
      const Markov_Transition_step_double &Q_dt, Model &m, const Step &p,
      std::ostream &os, MACROR &alg) const {
    auto y = p.y();
    if (std::isnan(y))
      alg = MACRO_DMNR;
    else {
      double mg = (prior.P_mean() * Q_dt.gmean_i()).getvalue();
      double g_max = max(Q_dt.g());
      double g_min = min(Q_dt.g());
      double g_range = g_max - g_min;
      double N = m.AverageNumberOfChannels();
      std::pair<Op_void, Op_void> test_Binomial;
      auto p_bi = (g_max - mg) / g_range;
      auto q_bi = (mg - g_min) / g_range;
      test_Binomial = mp_state_information::is_Binomial_Approximation_valid(
          N, p_bi, q_bi, Binomial_magic_number(), 0);
      if (!test_Binomial.first.has_value()) {
        if (!test_Binomial.second.has_value())
          alg = MACRO_DMNR;
        else
          alg = MACRO_DVNR;
      } else if (!test_Binomial.second.has_value())
        alg = MACRO_DMR;
      else
        alg = MACRO_DVR;
    }
    const MACROR alg2 = alg;
    return run(prior, Q_dt, m, p, os, alg2);
  }

  template <class Model, class Step>
  myOptional_t<mp_state_information>
  run(const mp_state_information &prior,
      const Markov_Transition_step_double &Q_dt, Model &m, const Step &p,
      std::ostream & /*os*/) const {
    auto y = p.y();
    typedef myOptional_t<mp_state_information> Op;
    double N = m.AverageNumberOfChannels();

    double e = m.noise_variance(p.nsamples());
    //
    M_Matrix<double> u(prior.P_mean().size(), 1, 1.0);

    auto SmD = prior.P_cov() - diag(prior.P_mean());

    double gSg =
        (TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.gmean_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_ij(), Q_dt.gmean_ij()) * u))
            .getvalue();

    double sSg =
        (TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.gmean_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_var_ij(), Q_dt.gmean_ij()) * u))
            .getvalue();
    double sSs =
        (TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.gvar_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_var_ij(), Q_dt.gvar_ij()) * u))
            .getvalue();
    //   auto mu_n=;
    auto sS = TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.P() +
              prior.P_mean() * Q_dt.gtotal_var_ij();
    auto gS = TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.P() +
              prior.P_mean() * Q_dt.gtotal_ij();

    double ms = (prior.P_mean() * Q_dt.gvar_i()).getvalue();
    double delta_emu = std::max(sqr(ms + e / N) - 2.0 / N * sSs, 0.0);
    double ms0 = (ms - e / N) / 2 + std::sqrt(delta_emu) / 2;

    auto e_mu = e + N * ms0;
    auto y_mean =
        N * (prior.P_mean() * Q_dt.gmean_i()).getvalue() - N * 0.5 / e_mu * sSg;
    auto zeta = N / (2 * sqr(e_mu) + N * sSs);
    auto y_var = std::max(e_mu + N * gSg - N * zeta * sqr(sSg), e);
    auto dy = y - y_mean;
    auto chi = dy / y_var;
    auto P_mean = prior.P_mean() * Q_dt.P() + chi * gS -
                  (chi * zeta * sSg + 0.5 / e_mu) * sS;

    auto P__cov =
        quadraticForm_BT_A_B(SmD, Q_dt.P()) + diag(prior.P_mean() * Q_dt.P()) -
        (zeta + N / y_var * sqr(zeta * sSg)) * quadraticForm_XTX(sS) +
        (2.0 * N / y_var * zeta * sSg) * TransposeSum(TranspMult(sS, gS)) -
        (N / y_var) * quadraticForm_XTX(gS);

    auto chi2 = dy * chi;

    double plogL;
    if (y_var > 0)
      plogL = -0.5 * log(2 * PI * y_var) - 0.5 * chi2;
    else
      plogL = std::numeric_limits<double>::infinity();

    double eplogL =
        -0.5 * log(2 * PI * y_var) - 0.5; // e_mu+N*gSg"-N*zeta*sqr(sSg)"
    // double chilogL=(eplogL-plogL)/std::sqrt(0.5);
    auto test = mp_state_information::test(P_mean, P__cov, y_mean, y_var, plogL,
                                           eplogL, e, tolerance());
    if (!test) {
      std::stringstream ss;
      ss << "\n step=" << p << " "
         << " y=" << y << " ymean=" << y_mean << " yvar=" << y_var << " e=" << e
         << " e_mu" << e_mu << " N*gSg" << N * gSg
         << " -N*zeta*sSg=" << -N * zeta * sqr(sSg);
      ss << "\n ms0=(ms-e/N)/2+std::sqrt(delta_emu)/2;\nms0=" << ms0
         << "\n ms=" << ms << "\n e/N/2=" << e / N;

      ss << "\n std::sqrt(delta_emu)/2=" << std::sqrt(delta_emu) / 2;
      ss << "\ndelta_emu=sqr(ms+e/N)-2.0/N*sSs \nsqr(ms+e/N)="
         << sqr(ms + e / N) << " 2.0/N*sSs=" << 2.0 / N * sSs << " sSs=" << sSs;
      ss << "\nP_mean \n" << P_mean;
      ss << " -N*zeta*sqr(sSg)=" << -N * zeta * sqr(sSg) << " plogL=" << plogL
         << " eplogL=" << eplogL << "dif=" << plogL - eplogL << " sSs=" << sSs;
      ss << "\nPcov \n" << P__cov;
      // ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;

      return Op(false, "\nfails in trace!!!; error=" + test.error() + ss.str());
    } else
      return Op(mp_state_information::adjust(std::move(P_mean),
                                             std::move(P__cov), y_mean, y_var,
                                             plogL, eplogL, Q_dt.min_P(), e));
  }

  template <class Model, class Step>
  myOptional_t<mp_state_information> start(Model &m, const Step &p,
                                           double min_p) const {
    typedef myOptional_t<mp_state_information> Op;

    auto P_mean = m.Peq(p.begin()->x());
    //  std::cerr<<"gslreijgsorjgps INIT!!!"<<P_mean.value();
    if (!P_mean)
      return Op(false, "fails to get Peq :" + P_mean.error());
    else {
      auto P_cov = diag(P_mean.value()) - quadraticForm_XTX(P_mean.value());
      //     std::cerr<<"gslreijgsorjgps INIT!!!"<<P_cov;
      double nan = std::numeric_limits<double>::quiet_NaN();
      return Op(mp_state_information::adjust(std::move(P_mean).value(),
                                             std::move(P_cov), nan, nan, nan,
                                             nan, min_p, 0));
    }
  }
  double tolerance() const { return tolerance_; }

  double Binomial_magic_number() const { return binomial_magic; }
  double Variance_magic_number() const { return variance_magic; }
  MacroDVR(double tolerance, double binomial_magical, double variance_magical)
      : tolerance_{tolerance}, binomial_magic(binomial_magical),
        variance_magic(variance_magical) {}
  MacroDVR() = default;

private:
  double tolerance_ = 1e-5;
  double binomial_magic = 5;
  double variance_magic = 1;
};
} // namespace hidden

template <bool recursive, int averaging, bool variance> class Macro_R;
template <template <class...> class Tr,bool recursive, int averaging, bool variance> class Macro_R_new_;
template <template <class...> class Tr,bool recursive, int averaging, bool variance> class Macro_Sel_;


template <bool recursive, int averaging, bool variance> using  Macro_R_new=Macro_R_new_<C,recursive,averaging,variance>;
template <bool recursive, int averaging, bool variance> using  Macro_Sel=Macro_Sel_<C,recursive,averaging,variance>;


template <int averaging, bool variance> class Micro_R;

template <int averaging, bool variance> class RMicro_R;
template <template <class...> class Tr,int averaging, bool variance> class RMicro_R_new_;

template <int averaging, bool variance> using  RMicro_R_new=RMicro_R_new_<C,averaging,variance>;

typedef Macro_R<false, 2, false> MacroDMNR;
typedef Macro_R<false, 2, true> MacroDVNR;
typedef Macro_R<true, 2, false> MacroDMR;
typedef Macro_R<true, 2, true> MacroDVR;

typedef Macro_R_new<false, 2, false> MacroDMNR_new;
typedef Macro_R_new<false, 2, true> MacroDVNR_new;
typedef Macro_R_new<true, 2, false> MacroDMR_new;
typedef Macro_R_new<true, 2, true> MacroDVR_new;

typedef Macro_Sel<false, 2, false> MacroSDMNR;
typedef Macro_Sel<false, 2, true> MacroSDVNR;
typedef Macro_Sel<true, 2, false> MacroSDMR;
typedef Macro_Sel<true, 2, true> MacroSDVR;




typedef Micro_R<2, false> MicroDMR;
typedef Micro_R<2, true> MicroDVR;

//typedef RMicro_R<2, false> RMicroDMR;
//typedef RMicro_R<2, true> RMicroDVR;


typedef RMicro_R_new<2, false> RMicroDMR_new;
typedef RMicro_R_new<2, true> RMicroDVR_new;

template<template <class...> class Tr,bool eigenvalue>
class Markov_Model_calculations_new_;

template<bool eigenvalue> using Markov_Model_calculations_new=Markov_Model_calculations_new_<C,eigenvalue>;


template<template <class...> class Tr,bool eigenvalue>
class Markov_Model_calculations_new_ {

    double fs_;
    double min_P_;

public:
    template<class T> using Tr_t=typename Tr<T>::type;
    Markov_Model_calculations_new_(double fs, double min_P): fs_{fs},min_P_{min_P}{}
    double fs()const {return fs_;}
    double min_P()const {return min_P_;}
    auto calc_Markov_Transition_rate(Tr_t<M_Matrix<double>> &&_Qrun, const Tr_t<M_Matrix<double>> &_g
                                     )const
    {
        typedef myOptional<Tr_t<Markov_Transition_rate_new>,reg_tag> Op;
        auto eig = calculate_eigenvalues()(_Qrun);
        if (eig) {
            //if (auto eigtest = test(eig.value(), tolerance); eigtest.has_value())
            return Op(Tr_t<Markov_Transition_rate_new>(std::move(_Qrun), _g,
                                             std::move(eig).value()));
            //            else
            //                return Op(false, "invalid eigenvalues " + eig.error());
        } else
            return Op(false, "eigenvalue decomposition fails:" + eig.error());
    }


    auto calc_step_variance(const Tr_t<Markov_Transition_rate_new>& Qx, std::size_t nsamples, double fs, double min_P)const
    {

        return Tr_t<Markov_Transition_step_calculations>().calc_step_variance(Qx,nsamples,fs,min_P);

    }



    template< template <template<class...> class,auto...> class Macro_R, auto...Ts, class Model, class X>
    auto get_Qrun(const Macro_R<Tr,Ts...> /*macro*/, const Model& m, const X& x)const
    {

        if constexpr (eigenvalue)
        {
            return calc_Markov_Transition_rate(m.Q(x),m.g(x));
        }
    }


    template< template <template<class...> class,auto...Ts> class Macro_R, auto...Ts, class Model, class Step>
    auto get_Peq(const Macro_R<Tr,Ts...>& macro,const Model& m, const Step& s)const
    {
        if constexpr (eigenvalue)
        {
            typedef myOptional_t<std::invoke_result_t<
                decltype(&Tr_t<Markov_Transition_rate_new>::calc_Peq), Tr_t<Markov_Transition_rate_new>>>
                Op;
            //  auto Qx=get_Qx(x);
            auto x=s.begin()->x();
            auto Qx = get_Qrun(macro,m,x);
            // std::cerr<<"calcQx"<<Qx.value().Qrun();
            if (!Qx)
                return Op(false, "cannot calculate Qx: " + Qx.error());
            else
                return Op(Qx.value().calc_Peq());

        }
    }

    template <class Model,class Step>
    auto operator()(const Tr_t<MacroDVR_new>& macro,const Model& m,
                    const mp_macro_information &/*prior*/,
                    const Step &p,
                    std::ostream &/*os*/) const {
        typedef myOptional_t<Tr_t<Markov_Transition_step_variance>> Op;
        if constexpr (eigenvalue)
        {
            auto Qrun=get_Qrun(macro,m,p.begin()->x());
            if (!Qrun)
                return Op(false, Qrun.error());
            else
                return minimum_to_step_variance(Qrun.value(),p.begin()->nsamples(),fs(), min_P());


        }

    }




    template <class Model,class Step>
    auto operator()(const MacroDMR_new& macro,const Model& m,
                    const mp_macro_information &/*prior*/,
                    const Step &p,
                    std::ostream &/*os*/) const {
        typedef myOptional_t<Markov_Transition_step_variance> Op;
        if constexpr (eigenvalue)
        {
            auto Qrun=get_Qrun(macro,m,p.begin()->x());
            if (!Qrun)
                return Op(false, Qrun.error());
            else
                return calculate_step_variance<Tr>::minimum_to_step_variance(std::move(Qrun).value(),p.begin()->nsamples(),fs(), min_P());


        }

    }

    template <class Model,class Step,bool recursive, int averaging, bool variance>
    auto operator()(const Macro_R_new_<Tr,recursive,averaging,variance>& macro,const Model& m,
                    const Tr_t<mp_macro_information> &/*prior*/,
                    const Step &p,
                    std::ostream &/*os*/) const {
        typedef myOptional<Tr_t<Markov_Transition_step_variance>, reg_tag> Op;
        if constexpr (eigenvalue)
        {
            auto Qrun=get_Qrun(macro,m,p.begin()->x());
            if (!Qrun)
                return Op(false, Qrun.error());
            else
                return calculate_step_variance<Tr>::minimum_to_step_variance(std::move(Qrun).value(),p.begin()->nsamples(),fs(), min_P());


        }

    }

        template <class Model,class Step,bool recursive, int averaging, bool variance>
    auto operator()(const Macro_Sel_<Tr,recursive,averaging,variance>& macro,const Model& m,
                    const Tr_t<mp_macro_information> &/*prior*/,
                    const Step &p,
                    std::ostream &/*os*/) const {
            typedef myOptional_t<Tr_t<Markov_Transition_step_variance>> Op;
        if constexpr (eigenvalue)
        {
            auto Qrun=get_Qrun(macro,m,p.begin()->x());
            if (!Qrun)
                return Op(false, Qrun.error());
            else
                return calculate_step_variance<Tr>::minimum_to_step_variance(std::move(Qrun).value(),p.begin()->nsamples(),fs(), min_P());


        }

    }


};

















template <template<class...> class Tr,bool recursive, int averaging, bool variance> class Macro_R_new_ {
public:
    template <typename T> using Tr_t=typename Tr<T>::type;
    inline constexpr static auto myClass() {
        if constexpr (!recursive) {
            if constexpr (!variance)
                return my_static_string("MacroDMNR_new");
            else
                return my_static_string("MacroDVNR_new");
        } else {
            if constexpr (!variance)
                return my_static_string("MacroDMR_new");
            else
                return my_static_string("MacroDVR_new");
        }
    }

    inline constexpr static auto const className = myClass();


    template <class Model,class Calculator, class Step>
    myOptional_t<Tr_t<mp_macro_information>> run(const Tr_t<mp_macro_information> &prior,
                                           Calculator& c,const Model &m, const Step &p,
                                           std::ostream &os)const  {
        typedef myOptional<Tr_t<mp_macro_information>,reg_tag> Op;
        auto Q_dto = c(*this,m,prior,p,os);
        if (!Q_dto)
            return Op(false, "fails in auto Q_dt=m.get_P(p,0) :" + Q_dto.error());
        auto  Q_dt = std::move(Q_dto).value();
        return run(prior, c,Q_dt, m, p, os);
    }


    template <class Model, class Calculator,class Step>
    myOptional_t<Tr_t<mp_macro_information>>
    run(const Tr_t<mp_macro_information> &prior,Calculator& c,
        const Tr_t<Markov_Transition_step_variance> &Q_dt, Model &m, const Step &p,
        std::ostream & /*os*/) const {
        typedef myOptional_t<Tr_t<mp_macro_information>> Op;
        auto y = p.y();
        auto e = m.noise_variance(p.nsamples(),c.fs());
        auto N = m.AverageNumberOfChannels();
        M_Matrix<double> u(prior.P_mean().size(), 1, 1.0);

        auto SmD = prior.P_cov() - diag(prior.P_mean());
        auto gSg =
            (TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.gmean_i()).getvalue() +
            (prior.P_mean() * (elemMult(Q_dt.gtotal_ij(), Q_dt.gmean_ij()) * u))
                .getvalue();

        auto ms = (prior.P_mean() * Q_dt.gvar_i()).getvalue();

        Tr_t<double> e_mu;
        Tr_t<double> y_mean;
        Tr_t<double> y_var;

        Tr_t<double> sSg;
        Tr_t<double> sSs;
        Tr_t<double> zeta;

        if constexpr ((!variance) && (!recursive)) {
            e_mu = e + N * ms;
            y_mean = N * (prior.P_mean() * Q_dt.gmean_i()).getvalue();
            y_var = e_mu + N * gSg;

        } else if constexpr (!variance && recursive) {

            Tr_t<double> ms = (prior.P_mean() * Q_dt.gvar_i()).getvalue();
            e_mu = e + N * ms;
            y_mean = N * (prior.P_mean() * Q_dt.gmean_i()).getvalue();
            y_var = e_mu + N * gSg;

        } else // (variance && (recursive || !recursive))
        {
            sSg = (TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.gmean_i()).getvalue() +
                (prior.P_mean() *
                 (elemMult(Q_dt.gtotal_var_ij(), Q_dt.gmean_ij()) * u))
                    .getvalue();
            sSs = (TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.gvar_i()).getvalue() +
                (prior.P_mean() *
                 (elemMult(Q_dt.gtotal_var_ij(), Q_dt.gvar_ij()) * u))
                    .getvalue();

            //typename decltype (sqr(ms + e / N))::ki k;
            using std::max;
            Tr_t<double> delta_emu = max(sqr(ms + e / N) - 2.0 / N * sSs, 0.0);
            Tr_t<double> ms0 = (ms - e / N) / 2 + sqrt(delta_emu) / 2;

            e_mu = e + N * ms0;
            y_mean = N * (prior.P_mean() * Q_dt.gmean_i()).getvalue() -
                N * 0.5 / e_mu * sSg;
            zeta = N / (2 * sqr(e_mu) + N * sSs);
            y_var = max(e_mu + N * gSg - N * zeta * sqr(sSg), e);
        }
        if (std::isnan(center(y))) {

            Tr_t<double> plogL = std::numeric_limits<double>::quiet_NaN();
            Tr_t<double> eplogL = std::numeric_limits<double>::quiet_NaN();
            Tr_t<double> vplogL = 0.0;
            auto P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P());
            auto P_mean = prior.P_mean() * Q_dt.P();
            P__cov += diag(P_mean);
            // std::cerr<<"\nPcov nana corr\n"<<P__cov<<"\nP_mean nana
            // corr\n"<<P_mean<<"\nQ.P \n"<<Q_dt.P();
//            auto test = mp_macro_information::test(P_mean, P__cov, tolerance_);
//            if (test.has_value())
            return Op(Tr_t<mp_macro_information>::adjust(
                    std::move(P_mean), std::move(P__cov), y_mean, y_var, plogL, eplogL,
                vplogL, Q_dt.min_P(), center(e)));
//            else
//                return Op(false, "fails at intertrace prediction!!: " + test.error());
        }

        Tr_t<double> dy = y - y_mean;
        Tr_t<double> chi = dy / y_var;
        Tr_t<M_Matrix<double>> P_mean;
        Tr_t<M_Matrix<double>> P__cov;
        if constexpr (!recursive) {
            P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P());

            P_mean = prior.P_mean() * Q_dt.P();
            P__cov += diag(P_mean);
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

            P__cov =
                quadraticForm_BT_A_B(SmD, Q_dt.P()) +
                diag(prior.P_mean() * Q_dt.P()) -
                (zeta + N / y_var * sqr(zeta * sSg)) * quadraticForm_XTX(sS) +
                (2.0 * N / y_var * zeta * sSg) * TransposeSum(TranspMult(sS, gS)) -
                (N / y_var) * quadraticForm_XTX(gS);
        }

        auto chi2 = dy * chi;

        Tr_t<double> plogL;
        if (center(y_var) > 0)
            plogL = -0.5 * log(2 * PI * y_var) - 0.5 * chi2;
        else
            plogL = std::numeric_limits<double>::infinity();

        auto eplogL =
            -0.5 * log(2 * PI * y_var) - 0.5; // e_mu+N*gSg"-N*zeta*sqr(sSg)"
        // double chilogL=(eplogL-plogL)/std::sqrt(0.5);

        double vplogL = 0.5;
//        auto test = mp_macro_information::test(P_mean, P__cov, y_mean, y_var, plogL,
//                                               eplogL, e, tolerance());
//        if constexpr (false) {
//            if (!test) {
//                std::stringstream ss;

//                ss << "\nP_mean \n" << P_mean;
//                ss << "\nPcov \n" << P__cov;
//                // ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;

//                return Op(false,
//                          "\nfails in trace!!!; error=" + test.error() + ss.str());
//            }
//        } else
        return Op(Tr_t<mp_macro_information>::adjust(
                std::move(P_mean), std::move(P__cov), y_mean, y_var, plogL, eplogL,
            vplogL, Q_dt.min_P(), center(e)));
    }




  };

  template <template <class...>class Tr,bool recursive, int averaging, bool variance> class Macro_Sel_ {
public:
    template <class T> using Tr_t=typename Tr<T>::type;
    inline constexpr static auto myClass() {
        if constexpr (!recursive) {
            if constexpr (!variance)
                return my_static_string("MacroSDMNR");
            else
                return my_static_string("MacroSDVNR");
        } else {
            if constexpr (!variance)
                return my_static_string("MacroSDMR");
            else
                return my_static_string("MacroSDVR");
        }
    }

    inline constexpr static auto const className = myClass();

    template <class Model,class Calculator, class Step>
    myOptional_t<Tr_t<mp_macro_information>> run(const Tr_t<mp_macro_information> &prior,
                                           Calculator& calc,
                                           const Model &m, const Step &p,
                                           std::ostream &os) const {
        MACROR_new alg;
        return run(prior, calc,m, p, os, alg);
    }



    template <class Model,class Calculator, class Step>
    myOptional_t<Tr_t<mp_macro_information>> run(const Tr_t<mp_macro_information> &prior,
                                           Calculator& c,const Model &m, const Step &p,
                                           std::ostream &os, MACROR_new &alg) const {
        typedef myOptional_t<Tr_t<mp_macro_information>> Op;
        auto Q_dto = c(*this,m,prior,p,os);
        if (!Q_dto)
            return Op(false, "fails in auto Q_dt=m.get_P(p,0) :" + Q_dto.error());
        auto Q_dt=std::move(Q_dto).value();
        if constexpr (!recursive && !variance)
                alg = MACRO_DMNR_new;

            else {
                auto y = p.y();

                if (std::isnan(y))
                    alg = MACRO_DMNR_new;
                else {
                    auto mg = (center(prior.P_mean()) * center(Q_dt.gmean_i())).getvalue();
                    auto g_max = max(center(Q_dt.g()));
                    auto g_min = min(center(Q_dt.g()));
                    auto g_range = g_max - g_min;
                    auto N = center(m.AverageNumberOfChannels());
                    auto p_bi = (g_max - mg) / g_range;
                    auto q_bi = (mg - g_min) / g_range;
                    auto test_Binomial = is_Binomial_Approximation_valid(
                        N, p_bi, q_bi, Binomial_magic_number(), Variance_magic_number());
                    if (!recursive || !test_Binomial.first) {
                        if (!variance || !test_Binomial.second)
                            alg = MACRO_DMNR_new;
                        else
                            alg = MACRO_DVNR_new;
                    } else if (!variance || !test_Binomial.second)
                        alg = MACRO_DMR_new;
                    else
                        alg = MACRO_DVR_new;
                }
            }
            const MACROR_new alg2 = alg;

            return run(prior, c, m, p, os, alg2);
        }



    template <class Model, class Calculator,class Step>
        myOptional_t<Tr_t<mp_macro_information>>
        run(const Tr_t<mp_macro_information> &prior,Calculator& calc,
        const Model &m, const Step &p,
        std::ostream &os, const MACROR_new &alg) const {
        switch (alg) {
        case MACRO_DMNR_new:
            return Tr_t<MacroDMNR_new>().run(prior, calc, m, p, os);
        case MACRO_DVNR_new:
            return Tr_t<MacroDVNR_new>().run(prior, calc, m, p, os);
        case MACRO_DMR_new:
            return Tr_t<MacroDMR_new>().run(prior, calc, m, p, os);
        case MACRO_DVR_new:
        default:
            return Tr_t<MacroDVR_new>().run(prior, calc, m, p, os);
        }
    }

    template <class Model, class Calculation,class Step>
    myOptional_t<Tr_t<mp_macro_information>> start(Calculation& calc,const Model &m, const Step &p,
                                                    double min_p)  const {
        typedef myOptional_t<Tr_t<mp_macro_information>> Op;

        auto P_mean = calc.get_Peq(*this,m,p);
        if (!P_mean)
            return Op(false, "fails to get Peq :" + P_mean.error());
        else {
            auto P_cov = diag(P_mean.value()) - quadraticForm_XTX(P_mean.value());
            Tr_t<double> nan = Tr_t<double>(std::numeric_limits<double>::quiet_NaN());
            Tr_t<double> zero = Tr_t<double>(0.0);


            return Op(Tr_t<mp_macro_information>::adjust(std::move(P_mean).value(),
                                                std::move(P_cov), nan, nan, nan,
                                                         nan, zero, min_p, 0));
        }
    }

    double tolerance() const { return tolerance_; }

    double Binomial_magic_number() const { return binomial_magic; }
    double Variance_magic_number() const { return variance_magic; }
    Macro_Sel_(double tolerance, double binomial_magical, double variance_magical)
        : tolerance_{tolerance}, binomial_magic(binomial_magical),
          variance_magic(variance_magical) {}
    Macro_Sel_() = default;

private:
    double tolerance_ = 1e-5;
    double binomial_magic = 5;
    double variance_magic = 1;
};







template <bool recursive, int averaging, bool variance> class Macro_R {
public:
  inline constexpr static auto myClass() {
    if constexpr (!recursive) {
      if constexpr (!variance)
        return my_static_string("MacroDMNR");
      else
        return my_static_string("MacroDVNR");
    } else {
      if constexpr (!variance)
        return my_static_string("MacroDMR");
      else
        return my_static_string("MacroDVR");
    }
  }

  inline constexpr static auto const className = myClass();

  template <class Model, class Step>
  myOptional_t<mp_state_information> run(const mp_state_information &prior,
                                         Model &m, const Step &p,
                                         std::ostream &os) const {
    MACROR alg;
    return run(prior, m, p, os, alg);
  }

  template <class Model, class Step, class MACROR>
  myOptional_t<mp_state_information> run(const mp_state_information &prior,
                                         Model &m, const Step &p,
                                         std::ostream &os, MACROR &alg) const {
    typedef myOptional_t<mp_state_information> Op;
    auto Q_dto = m.get_P(p, 0);
    if (!Q_dto)
      return Op(false, "fails in auto Q_dt=m.get_P(p,0) :" + Q_dto.error());
    Markov_Transition_step_double Q_dt = std::move(Q_dto).value();
    return run(prior, Q_dt, m, p, os, alg);
  }

  template <class Model, class Step>
  myOptional_t<mp_state_information>
  run(const mp_state_information &prior,
      const Markov_Transition_step_double& Q_dt, Model &m, const Step &p,
      std::ostream &os,
      const std::tuple<MACROR, mp_state_information> &alg) const {
    return run(prior, Q_dt, m, p, os, std::get<0>(alg));
  }
  template <class Model, class Step>
  myOptional_t<mp_state_information>
  run(const mp_state_information &prior,
      const Markov_Transition_step_double& Q_dt, Model &m, const Step &p,
      std::ostream &os, const MACROR &alg) const {
    switch (alg) {
    case MACRO_DMNR:
      return MacroDMNR(tolerance(), Binomial_magic_number(),
                       Variance_magic_number())
          .run(prior, Q_dt, m, p, os);
    case MACRO_DVNR:
      return MacroDVNR(tolerance(), Binomial_magic_number(),
                       Variance_magic_number())
          .run(prior, Q_dt, m, p, os);
    case MACRO_DMR:
      return MacroDMR(tolerance(), Binomial_magic_number(),
                      Variance_magic_number())
          .run(prior, Q_dt, m, p, os);
    case MACRO_DVR:
    default:
      return MacroDVR(tolerance(), Binomial_magic_number(),
                      Variance_magic_number())
          .run(prior, Q_dt, m, p, os);
    }
  }

  template <class Model, class Step>
  myOptional_t<mp_state_information>
  run(const mp_state_information &prior,
      const Markov_Transition_step_double &Q_dt, Model &m, const Step &p,
      std::ostream &os, MACROR &alg) const {
    if constexpr (!recursive && !variance)
      alg = MACRO_DMNR;

    else {
      auto y = p.y();

      if (std::isnan(y))
        alg = MACRO_DMNR;
      else {
        double mg = (prior.P_mean() * Q_dt.gmean_i()).getvalue();
        double g_max = max(Q_dt.g());
        double g_min = min(Q_dt.g());
        double g_range = g_max - g_min;
        double N = m.AverageNumberOfChannels();
        std::pair<Op_void, Op_void> test_Binomial;
        auto p_bi = (g_max - mg) / g_range;
        auto q_bi = (mg - g_min) / g_range;
        test_Binomial = mp_state_information::is_Binomial_Approximation_valid(
            N, p_bi, q_bi, Binomial_magic_number(), Variance_magic_number());
        if (!recursive || !test_Binomial.first.has_value()) {
          if (!variance || !test_Binomial.second.has_value())
            alg = MACRO_DMNR;
          else
            alg = MACRO_DVNR;
        } else if (!variance || !test_Binomial.second.has_value())
          alg = MACRO_DMR;
        else
          alg = MACRO_DVR;
      }
    }
    const MACROR alg2 = alg;

    return run(prior, Q_dt, m, p, os, alg2);
  }

  template <class Model, class Step>
  myOptional_t<mp_state_information>
  run(const mp_state_information &prior,
      const Markov_Transition_step_double &Q_dt, Model &m, const Step &p,
      std::ostream & /*os*/) const {
    typedef myOptional_t<mp_state_information> Op;
    auto y = p.y();
    double e = m.noise_variance(p.nsamples());
    double N = m.AverageNumberOfChannels();
    M_Matrix<double> u(prior.P_mean().size(), 1, 1.0);

    auto SmD = prior.P_cov() - diag(prior.P_mean());
    double gSg =
        (TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.gmean_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_ij(), Q_dt.gmean_ij()) * u))
            .getvalue();

    double ms = (prior.P_mean() * Q_dt.gvar_i()).getvalue();

    double e_mu;
    double y_mean;
    double y_var;

    double sSg;
    double sSs;
    double zeta;

    if constexpr ((!variance) && (!recursive)) {
      e_mu = e + N * ms;
      y_mean = N * (prior.P_mean() * Q_dt.gmean_i()).getvalue();
      y_var = e_mu + N * gSg;

    } else if constexpr (!variance && recursive) {
      auto gS = TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.P() +
                prior.P_mean() * Q_dt.gtotal_ij();

      double ms = (prior.P_mean() * Q_dt.gvar_i()).getvalue();

      e_mu = e + N * ms;
      y_mean = N * (prior.P_mean() * Q_dt.gmean_i()).getvalue();
      y_var = e_mu + N * gSg;

    } else // (variance && (recursive || !recursive))
    {
      sSg = (TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.gmean_i()).getvalue() +
            (prior.P_mean() *
             (elemMult(Q_dt.gtotal_var_ij(), Q_dt.gmean_ij()) * u))
                .getvalue();
      sSs = (TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.gvar_i()).getvalue() +
            (prior.P_mean() *
             (elemMult(Q_dt.gtotal_var_ij(), Q_dt.gvar_ij()) * u))
                .getvalue();

      double delta_emu = std::max(sqr(ms + e / N) - 2.0 / N * sSs, 0.0);
      double ms0 = (ms - e / N) / 2 + std::sqrt(delta_emu) / 2;

      e_mu = e + N * ms0;
      y_mean = N * (prior.P_mean() * Q_dt.gmean_i()).getvalue() -
               N * 0.5 / e_mu * sSg;
      zeta = N / (2 * sqr(e_mu) + N * sSs);
      y_var = std::max(e_mu + N * gSg - N * zeta * sqr(sSg), e);
    }
    if (std::isnan(y)) {

      double plogL = std::numeric_limits<double>::quiet_NaN();
      double eplogL = std::numeric_limits<double>::quiet_NaN();
      double vplogL = 0.0;
      auto P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P());
      auto P_mean = prior.P_mean() * Q_dt.P();
      P__cov += diag(P_mean);
      // std::cerr<<"\nPcov nana corr\n"<<P__cov<<"\nP_mean nana
      // corr\n"<<P_mean<<"\nQ.P \n"<<Q_dt.P();
      auto test = mp_state_information::test(P_mean, P__cov, tolerance_);
      if (test.has_value())
        return Op(mp_state_information::adjust(
            std::move(P_mean), std::move(P__cov), y_mean, y_var, plogL, eplogL,
            vplogL, Q_dt.min_P(), e));
      else
        return Op(false, "fails at intertrace prediction!!: " + test.error());
    }

    auto dy = y - y_mean;
    auto chi = dy / y_var;
    M_Matrix<double> P_mean;
    M_Matrix<double> P__cov;
    if constexpr (!recursive) {
      P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P());

      P_mean = prior.P_mean() * Q_dt.P();
      P__cov += diag(P_mean);
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

      P__cov =
          quadraticForm_BT_A_B(SmD, Q_dt.P()) +
          diag(prior.P_mean() * Q_dt.P()) -
          (zeta + N / y_var * sqr(zeta * sSg)) * quadraticForm_XTX(sS) +
          (2.0 * N / y_var * zeta * sSg) * TransposeSum(TranspMult(sS, gS)) -
          (N / y_var) * quadraticForm_XTX(gS);
    }

    auto chi2 = dy * chi;

    double plogL;
    if (y_var > 0)
      plogL = -0.5 * log(2 * PI * y_var) - 0.5 * chi2;
    else
      plogL = std::numeric_limits<double>::infinity();

    double eplogL =
        -0.5 * log(2 * PI * y_var) - 0.5; // e_mu+N*gSg"-N*zeta*sqr(sSg)"
    // double chilogL=(eplogL-plogL)/std::sqrt(0.5);

    double vplogL = 0.5;
    auto test = mp_state_information::test(P_mean, P__cov, y_mean, y_var, plogL,
                                           eplogL, e, tolerance());
    if constexpr (false) {
      if (!test) {
        std::stringstream ss;

        ss << "\nP_mean \n" << P_mean;
        ss << "\nPcov \n" << P__cov;
        // ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;

        return Op(false,
                  "\nfails in trace!!!; error=" + test.error() + ss.str());
      }
    } else
      return Op(mp_state_information::adjust(
          std::move(P_mean), std::move(P__cov), y_mean, y_var, plogL, eplogL,
          vplogL, Q_dt.min_P(), e));
  }




  template <class Model, class Step>
  myOptional_t<mp_state_information> start(Model &m, const Step &p,
                                           double min_p) const {
    typedef myOptional_t<mp_state_information> Op;

    auto P_mean = m.Peq(p.begin()->x());
    if (!P_mean)
      return Op(false, "fails to get Peq :" + P_mean.error());
    else {
      auto P_cov = diag(P_mean.value()) - quadraticForm_XTX(P_mean.value());
      double nan = std::numeric_limits<double>::quiet_NaN();

      return Op(mp_state_information::adjust(std::move(P_mean).value(),
                                             std::move(P_cov), nan, nan, nan,
                                             nan, 0.0, min_p, 0));
    }
  }
  double tolerance() const { return tolerance_; }

  double Binomial_magic_number() const { return binomial_magic; }
  double Variance_magic_number() const { return variance_magic; }
  Macro_R(double tolerance, double binomial_magical, double variance_magical)
      : tolerance_{tolerance}, binomial_magic(binomial_magical),
        variance_magic(variance_magical) {}
  Macro_R() = default;

private:
  double tolerance_ = 1e-5;
  double binomial_magic = 5;
  double variance_magic = 1;
};






template <int averaging, bool variance> class Micro_R {
public:
  inline constexpr static auto myClass() {
    if constexpr (!variance)
      return my_static_string("MicroDMR");
    else
      return my_static_string("MicroDVR");
  }

  inline constexpr static auto const className = myClass();

  template <class Model, class Step>
  myOptional_t<mp_single_state_information>
  run(const mp_single_state_information &prior, Model &m, const Step &p,
      std::ostream &os) const {
    MICROR alg;
    return run(prior, m, p, os, alg);
  }

  template <class Model, class Step, class MICROR>
  myOptional_t<mp_single_state_information>
  run(const mp_single_state_information &prior, Model &m, const Step &p,
      std::ostream &os, MICROR &alg) const {
    typedef myOptional_t<mp_single_state_information> Op;
    auto Q_dto = m.get_P(p, 0);
    if (!Q_dto)
      return Op(false, "fails in auto Q_dt=m.get_P(p,0) :" + Q_dto.error());
    Markov_Transition_step_double Q_dt = std::move(Q_dto).value();
    return run(prior, Q_dt, m, p, os, alg);
  }

  template <class Model, class Step>
  myOptional_t<mp_single_state_information>
  run(const mp_single_state_information &prior,
      const Markov_Transition_step_double& Q_dt, Model &m, const Step &p,
      std::ostream &os,
      const std::tuple<MICROR, mp_single_state_information> &alg) const {
    return run(prior, Q_dt, m, p, os, std::get<0>(alg));
  }
  template <class Model, class Step>
  myOptional_t<mp_single_state_information>
  run(const mp_single_state_information &prior,
      const Markov_Transition_step_double& Q_dt, Model &m, const Step &p,
      std::ostream &os, const MICROR &alg) const {
    switch (alg) {
    case MICRO_DMR:
      return MicroDMR().run(prior, Q_dt, m, p, os);
    case MICRO_DVR:
    default:
      return MicroDVR().run(prior, Q_dt, m, p, os);
    }
  }

  template <class Model, class Step>
  myOptional_t<mp_single_state_information>
  run(const mp_single_state_information &prior,
      const Markov_Transition_step_double &Q_dt, Model &m, const Step &p,
      std::ostream &os, MICROR &alg) const {
    if constexpr (!variance)
      alg = MICRO_DMR;

    else
      alg = MICRO_DVR;

    const MICROR alg2 = alg;

    return run(prior, Q_dt, m, p, os, alg2);
  }

  template <class Model, class Step>
  myOptional_t<mp_single_state_information>
  run(const mp_single_state_information &prior,
      const Markov_Transition_step_double &Q_dt, Model &m, const Step &p,
      std::ostream & /*os*/) const {
    typedef myOptional_t<mp_single_state_information> Op;
    auto y = p.y();
    double e = m.noise_variance(p.nsamples());
    auto k = prior.P().size();
    M_Matrix<double> Post(1, prior.P().size(), 0);
    M_Matrix<double> logLik(prior.P().size(), prior.P().size(), 0);
    double partialLik = 0;
    [[maybe_unused]] double epartialLik = 0;
    double y_mean = 0;
    double y_var = e;
    double eplogL;
    double vplogL;
    double plogL;
    if (std::isnan(y)) {
      Post = prior.P() * Q_dt.P();
      y_mean = std::numeric_limits<double>::quiet_NaN();
      plogL = std::numeric_limits<double>::quiet_NaN();
      y_var = std::numeric_limits<double>::quiet_NaN();
      eplogL = std::numeric_limits<double>::quiet_NaN();
      vplogL = std::numeric_limits<double>::quiet_NaN();
    } else {
        double maxlogLik=-std::numeric_limits<double>::infinity();
        for (std::size_t j = 0; j < k; ++j) {
            for (std::size_t i = 0; i < k; ++i) {
                if (prior.P()[i] * Q_dt.P()(i, j) > Q_dt.min_P()) {
                    logLik(i,j) = log_gaussian_function(y, Q_dt.gmean_ij()(i, j),
                                          e + std::abs(Q_dt.gvar_ij()(i, j)));
                    if (logLik(i,j)>maxlogLik) maxlogLik=logLik(i,j);
                }
            }
        }
      for (std::size_t j = 0; j < k; ++j) {
        for (std::size_t i = 0; i < k; ++i) {
          if (prior.P()[i] * Q_dt.P()(i, j) > Q_dt.min_P()) {
              Post[j] += prior.P()[i] * Q_dt.P()(i, j) *std::exp(logLik(i,j)-maxlogLik);
            //   y_mean +=prior.P()[i] * Q_dt.P()(i, j) *Q_dt.gmean_ij()(i, j);
            //   y_var +=prior.P()[i] * Q_dt.P()(i, j)
            //   *std::abs(Q_dt.gvar_ij()(i, j));
          }
        }
        partialLik += Post[j];
        y_mean += prior.P()[j] * Q_dt.gmean_i()[j];
        y_var += prior.P()[j] * Q_dt.gvar_i()[j];
      }
      for (std::size_t j = 0; j < k; ++j) {
        Post[j] /= partialLik;
      }
      assert(partialLik>0);
      plogL = std::log(partialLik)+maxlogLik;
      eplogL = -0.5 * log(2 * PI * y_var) - 0.5;
      vplogL = 0.5;
      assert(std::isfinite(eplogL));
    }
    return Op(mp_single_state_information::adjust(std::move(Post), y_mean,
                                                  y_var, plogL, eplogL, vplogL,
                                                  Q_dt.min_P(), e));
  }

  template <class Model, class Step>
  myOptional_t<mp_single_state_information> start(Model &m, const Step &p,
                                                  double min_p) const {
    typedef myOptional_t<mp_single_state_information> Op;

    auto Peq = m.Peq(p.begin()->x());
    if (!Peq)
      return Op(false, "fails to get Peq :" + Peq.error());
    else {
      double nan = std::numeric_limits<double>::quiet_NaN();

      return Op(mp_single_state_information::adjust(
          std::move(Peq).value(), nan, nan, nan, nan, 0.0, min_p, 0));
    }
  }
  Micro_R() = default;
};

//template <int averaging, bool variance> class RMicro_R {
//public:
//  inline constexpr static auto myClass() {
//    if constexpr (!variance)
//      return my_static_string("RMicroDMR");
//    else
//      return my_static_string("RMicroDVR");
//  }

//  inline constexpr static auto const className = myClass();

//  template <class Model, class Step>
//  myOptional_t<mp_ensamble_state_information>
//  run(const mp_ensamble_state_information &prior, Model &m, const Step &p,
//      std::ostream &os) const {
//    RMICROR alg;
//    return run(prior, m, p, os, alg);
//  }

//  template <class Model, class Step, class RMICROR>
//  myOptional_t<mp_ensamble_state_information>
//  run(const mp_ensamble_state_information &prior, Model &m, const Step &p,
//      std::ostream &os, RMICROR &alg) const {
//    typedef myOptional_t<mp_ensamble_state_information> Op;
//    auto Q_dto = m.get_P(p, 0);
//    if (!Q_dto)
//      return Op(false, "fails in auto Q_dt=m.get_P(p,0) :" + Q_dto.error());
//    Markov_Transition_step_double Q_dt = std::move(Q_dto).value();
//    return run(prior, Q_dt, m, p, os, alg);
//  }

//  template <class Model, class Step>
//  myOptional_t<mp_ensamble_state_information>
//  run(const mp_ensamble_state_information &prior,
//      const Markov_Transition_step_double& Q_dt, Model &m, const Step &p,
//      std::ostream &os,
//      const std::tuple<RMICROR, mp_ensamble_state_information> &alg) const {
//    return run(prior, Q_dt, m, p, os, std::get<0>(alg));
//  }
//  template <class Model, class Step>
//  myOptional_t<mp_ensamble_state_information>
//  run(const mp_ensamble_state_information &prior,
//      const Markov_Transition_step_double& Q_dt, Model &m, const Step &p,
//      std::ostream &os, const RMICROR &alg) const {
//    switch (alg) {
//    case RMICRO_DMR:
//      return RMicroDMR(maxChi2(), reduce_by_p(), min_connection_ratio())
//          .run(prior, Q_dt, m, p, os);
//    case RMICRO_DVR:
//    default:
//      return RMicroDVR(maxChi2(), reduce_by_p(), min_connection_ratio())
//          .run(prior, Q_dt, m, p, os);
//    }
//  }

//  template <class Model, class Step>
//  myOptional_t<mp_ensamble_state_information>
//  run(const mp_ensamble_state_information &prior,
//      const Markov_Transition_step_double &Q_dt, Model &m, const Step &p,
//      std::ostream &os, RMICROR &alg) const {
//    if constexpr (!variance)
//      alg = RMICRO_DMR;

//    else
//      alg = RMICRO_DVR;

//    const RMICROR alg2 = alg;

//    return run(prior, Q_dt, m, p, os, alg2);
//  }
//  Markov_Transition_step &
//  calc_Qdt(std::map<double, Markov_Transition_step> &Qdt_map,
//           const Markov_Transition_rate &Qx, std::size_t nsamples, double min_P,
//           double fs) const {
//    if (Qdt_map.find(fs) == Qdt_map.end()) {
//      auto Q_dtrun = Markov_Transition_step(Qx, nsamples, fs, min_P);
//      Qdt_map[fs] = std::move(Q_dtrun);
//    }
//    return Qdt_map[fs];
//  }

//  auto calc_new_Microscopic_Description(
//      const std::vector<std::vector<std::size_t>> &connections_to_x,
//      const std::vector<std::vector<std::size_t>> &connections_from_x,

//      std::size_t number_of_channels, const Markov_Transition_step &Q_dtrun,
//      const M_Matrix<double> &Pmean_0, const M_Matrix<double> &Pcov_0) const {
//    auto P_mean_1 = Pmean_0 * Q_dtrun.P();
//    auto SmD_1 = Pcov_0 - diag(Pmean_0);
//    auto P_cov_1 = quadraticForm_BT_A_B(SmD_1, Q_dtrun.P());
//    P_cov_1 += diag(P_mean_1);
//    auto mi_1 =
//        Microscopic_description(number_of_channels, P_mean_1, P_cov_1,
//                                connections_to_x,connections_from_x, maxChi2());
//    return std::make_tuple(P_mean_1, P_cov_1, mi_1);
//  }

//  M_Matrix<double> extend_description(
//      const Markov_Transition_rate &Qx,
//      const std::vector<std::vector<std::size_t>> &connections_to_x,
//      const std::vector<std::vector<std::size_t>> &connections_from_x,
//      const M_Matrix<double> &P, Microscopic_description &mi,
//      Microscopic_description &&mi_final, Microscopic_description &&mi_1,
//      M_Matrix<double> &&Pmean_0, M_Matrix<double> &&Pcov_0,
//      M_Matrix<double> &&Pmean_1, M_Matrix<double> &&Pcov_1,
//      std::map<double, Markov_Transition_step> &&Qdt_map, std::size_t nsamples,
//      double min_P, double fs,
//      std::vector<std::tuple<double, M_Matrix<double>, M_Matrix<double>,
//                             Microscopic_description>> &&history) const {
//    if (mi.connecting_ratio_closure(mi_1,connections_to_x, connections_from_x) < min_connection_ratio()) {
//      fs = fs * 2;
//      history.emplace_back(fs, std::move(Pmean_1), std::move(Pcov_1),
//                           std::move(mi_1));
//      auto &Q_dt = calc_Qdt(Qdt_map, Qx, nsamples, min_P, fs);
//      auto [Pmean_1, Pcov_1, mi_1] = calc_new_Microscopic_Description(
//          connections_to_x, connections_from_x, mi.number_of_channels(), Q_dt,
//          Pmean_0, Pcov_0);
//   //   std::cerr<<"fs*2 ="<<fs<<"Pmean_0"<<Pmean_0<<"Pmean_1"<<Pmean_1<<"\n";
//      return extend_description(
//          Qx, connections_to_x, connections_from_x, P, mi, std::move(mi_final),
//          std::move(mi_1), std::move(Pmean_0), std::move(Pcov_0),
//          std::move(Pmean_1), std::move(Pcov_1), std::move(Qdt_map), nsamples,
//          min_P, fs, std::move(history));

//    } else {
//        mi.merge(std::move(mi_1));
//        if (mi.connecting_ratio_closure(mi_final, connections_to_x, connections_from_x) < min_connection_ratio()) {
//          Pmean_0 = std::move(Pmean_1);
//          Pcov_0 = std::move(Pcov_1);
//          auto [fs, Pmean_1, Pcov_1, mi_1] = std::move(history.back());
//          history.pop_back();
//         // std::cerr<<"fs pop ="<<fs<<"Pmean_0"<<Pmean_0<<"Pmean_1"<<Pmean_1<<"\n";
//          return extend_description(
//              Qx, connections_to_x, connections_from_x ,P, mi, std::move(mi_final),
//              std::move(mi_1), std::move(Pmean_0), std::move(Pcov_0),
//              std::move(Pmean_1), std::move(Pcov_1), std::move(Qdt_map),
//              nsamples, min_P, fs, std::move(history));

//        } else
//          return mi.merge_P(P, std::move(mi_final),connections_to_x, connections_from_x);
//      }
//  }

//  template <class Model, class Step>
//  myOptional_t<mp_ensamble_state_information>
//  run(const mp_ensamble_state_information &prior,
//      const Markov_Transition_step_double &Q_dt, Model &m, const Step &p,
//      std::ostream &os) const {
//    typedef myOptional_t<mp_ensamble_state_information> Op;
//    /**
//     * Primero tengo que extender la descripcion microscopica a los estados
//     * plausibles luego de esta step
//     *
//     */
//    Microscopic_description mi = prior.microscopic_description();
//    auto [Pmean_0, Pcov_0] = prior.to_Macroscopic();

//    auto x=p.x().x();
//    auto [Pmean_1, Pcov_1, mi_1] = calc_new_Microscopic_Description(
//        m.connections_to_x(x), m.connections_from_x(x), mi.number_of_channels(), Q_dt,
//        Pmean_0, Pcov_0);
//    auto Qx = m.calc_Qx(x);

//    double fs = Q_dt.fs() * 2;
//    auto nsamples = Q_dt.nsamples();
//    auto P_ext = extend_description(
//        Qx.value(), m.connections_to_x(x),m.connections_from_x(x), prior.P(), mi,
//        std::move(mi_1), Microscopic_description(mi_1), std::move(Pmean_0), std::move(Pcov_0),
//        std::move(Pmean_1), std::move(Pcov_1),
//        std::map<double, Markov_Transition_step>(), nsamples, Q_dt.min_P(), fs,
//        std::vector<std::tuple<double, M_Matrix<double>, M_Matrix<double>,
//                               Microscopic_description>>());
//    assert(P_ext.size()==mi.size());
//    auto Qmi_dt = m.get_Pmi(mi, p, 0);
//    if (!Qmi_dt)
//      return Op(false, Qmi_dt.error());
//    mp_single_state_information prior_ext(std::move(P_ext), prior.y_mean(),
//                                          prior.y_var(), prior.plogL(),
//                                          prior.eplogL(), prior.vplogL());
//    auto post = MicroDVR().run(prior_ext, Qmi_dt.value(), m, p, os);
//    if (post.has_value()) {
//      mp_single_state_information mp = std::move(post).value();
//      auto e = m.noise_variance(p.nsamples());
//      return Op(mp_ensamble_state_information::adjust(
//          std::move(mi), std::move(mp.release_P()), mp.y_mean(), mp.y_var(),
//          mp.plogL(), mp.eplogL(), mp.vplogL(), Q_dt.min_P(), e,
//          reduce_by_p(),m.connections_to_x(x),m.connections_from_x(x)));
//    } else
//      return Op(false, post.error());
//  }

//  template <class Model, class Step>
//  myOptional_t<mp_ensamble_state_information_new>
//  run(const mp_ensamble_state_information_new &prior,
//      const Markov_Transition_step_double &Q_dt, Model &m, const Step &p,
//      std::ostream &os) const {
//      typedef myOptional_t<mp_ensamble_state_information_new> Op;
//      /**
//     * Primero tengo que extender la descripcion microscopica a los estados
//     * plausibles luego de esta step
//     *
//     */
//      Microscopic_description_new mi = prior.microscopic_description();
//      auto [Pmean_0, Pcov_0] = prior.to_Macroscopic();

//      auto x=p.x().x();
//      auto [Pmean_1, Pcov_1, mi_1] = calc_new_Microscopic_Description(
//          m.connections_to_x(x), m.connections_from_x(x), mi.number_of_channels(), Q_dt,
//          Pmean_0, Pcov_0);
//      auto Qx = m.calc_Qx(x);

//      double fs = Q_dt.fs() * 2;
//      auto nsamples = Q_dt.nsamples();
//      auto P_ext = extend_description(
//          Qx.value(), m.connections_to_x(x),m.connections_from_x(x), prior.p(), mi,
//          std::move(mi_1), Microscopic_description(mi_1), std::move(Pmean_0), std::move(Pcov_0),
//          std::move(Pmean_1), std::move(Pcov_1),
//          std::map<double, Markov_Transition_step>(), nsamples, Q_dt.min_P(), fs,
//          std::vector<std::tuple<double, M_Matrix<double>, M_Matrix<double>,
//                                 Microscopic_description>>());
//      assert(P_ext.size()==mi.size());
//      auto Qmi_dt = m.get_Pmi(mi, p, 0);
//      if (!Qmi_dt)
//          return Op(false, Qmi_dt.error());
//      mp_single_state_information prior_ext(std::move(P_ext), prior.y_mean(),
//                                            prior.y_var(), prior.plogL(),
//                                            prior.eplogL(), prior.vplogL());
//      auto post = MicroDVR().run(prior_ext, Qmi_dt.value(), m, p, os);
//      if (post.has_value()) {
//          mp_single_state_information mp = std::move(post).value();
//          auto e = m.noise_variance(p.nsamples());
//          return Op(mp_ensamble_state_information::adjust(
//              std::move(mi), std::move(mp.release_P()), mp.y_mean(), mp.y_var(),
//              mp.plogL(), mp.eplogL(), mp.vplogL(), Q_dt.min_P(), e,
//              reduce_by_p(),m.connections_to_x(x),m.connections_from_x(x)));
//      } else
//          return Op(false, post.error());
//  }


//  template <class Model, class Step>
//  myOptional_t<mp_ensamble_state_information> start(Model &m, const Step &p,
//                                                    double min_p) const {
//    typedef myOptional_t<mp_ensamble_state_information> Op;
//    auto x=p.begin()->x();
//    auto Peq = m.Peq(x);
//    if (!Peq)
//      return Op(false, "fails to get Peq :" + Peq.error());
//    else {
//      auto P_cov = diag(Peq.value()) - quadraticForm_XTX(Peq.value());

//      double nan = std::numeric_limits<double>::quiet_NaN();
//      std::size_t N = m.AverageNumberOfChannels();
//      Microscopic_description mi(m, x,N, Peq.value(), P_cov, maxChi2());

//      auto P = mi.convert_to_Microscopic_maxent(Peq.value(), P_cov);
//      if (!P)
//          return Op(false,"Maximum entropy failed: "+P.error());
//      else

//      return Op(mp_ensamble_state_information::adjust(
//              std::move(mi), std::move(P).value(), nan, nan, nan, nan, 0.0, min_p, 0,
//          reduce_by_p(),m.connections_to_x(x),m.connections_from_x(x)));
//    }
//  }


//  template <class Model, class Step>
//  myOptional_t<mp_ensamble_state_information_new> start(Model &m, const Step &p,
//                                                    double min_p,std::size_t order) const {
//      typedef myOptional_t<mp_ensamble_state_information_new> Op;
//      auto x=p.begin()->x();
//      auto Peq = m.Peq(x);
//      if (!Peq)
//          return Op(false, "fails to get Peq :" + Peq.error());
//      else {
//          auto P_cov = diag(Peq.value()) - quadraticForm_XTX(Peq.value());

//          double nan = std::numeric_limits<double>::quiet_NaN();
//          std::size_t N = m.AverageNumberOfChannels();
//          Microscopic_description mi(m, x,N, Peq.value(), P_cov, maxChi2());

//          auto P = mi.convert_to_Microscopic_maxent(Peq.value(), P_cov);
//          if (!P)
//              return Op(false,"Maximum entropy failed: "+P.error());
//          else

//          return Op(mp_ensamble_state_information::adjust(
//              std::move(mi), std::move(P).value(), nan, nan, nan, nan, 0.0, min_p, 0,
//              reduce_by_p(),m.connections_to_x(x),m.connections_from_x(x)));
//      }
//  }


//  RMicro_R(double maxChi2, double reduce_by_p, double min_connection_ratio__)
//      : maxChi2_{maxChi2}, reduce_by_p_{reduce_by_p},
//        min_connection_ratio_{min_connection_ratio__} {}
//  RMicro_R() = default;
//  double maxChi2() const { return maxChi2_; }
//  double reduce_by_p() const { return reduce_by_p_; }
//  double min_connection_ratio() const { return min_connection_ratio_; }

//private:
//  double maxChi2_;
//  double reduce_by_p_;
//  double min_connection_ratio_;
//};

template <template <class...> class Tr,int averaging, bool variance> class RMicro_R_new_ {
public:
    template<class T> using Tr_t=typename Tr<T>::type;
    inline constexpr static auto myClass() {
        if constexpr (!variance)
            return my_static_string("RMicroDMR_new");
        else
            return my_static_string("RMicroDVR_new");
    }

    inline constexpr static auto const className = myClass();

    template <class Model, class Step>
    myOptional_t<Tr_t<mp_ensamble_state_information_new>>
    run(const mp_ensamble_state_information_new &prior, Model &m, const Step &p,
        std::ostream &os) const {
        MICROR alg;
        return run(prior, m, p, os, alg);
    }

    template <class Model, class Step, class MICROR>
    myOptional_t<mp_ensamble_state_information_new>
    run(const mp_ensamble_state_information_new &prior, Model &m, const Step &s,
        std::ostream &os, MICROR &alg) const {
        typedef myOptional_t<mp_ensamble_state_information_new> Op;
        auto Q_dto = m.get_Pmi_new(prior.microscopic_description(),prior.p(),s, order_number(), reduce_by_p());
        if (!Q_dto)
            return Op(false, "fails in auto Q_dt=m.get_P(p,0) :" + Q_dto.error());
        Microscopic_Transition_step_double Q_dt = std::move(Q_dto).value();
        return run(prior, Q_dt, m, s, os, alg);
    }

    template <class Model, class Step>
    myOptional_t<mp_ensamble_state_information_new>
    run(const mp_ensamble_state_information_new &prior,
        const Microscopic_Transition_step_double& Q_dt, Model &m, const Step &p,
        std::ostream &os,
        const std::tuple<MICROR, mp_ensamble_state_information_new> &alg) const {
        return run(prior, Q_dt, m, p, os, std::get<0>(alg));
    }
    template <class Model, class Step>
    myOptional_t<mp_ensamble_state_information_new>
    run(const mp_ensamble_state_information_new &prior,
        const Microscopic_Transition_step_double& Q_dt, Model &m, const Step &p,
        std::ostream &os, const MICROR &alg) const {
        switch (alg) {
        case MICRO_DMR:
            return RMicroDMR_new().run(prior, Q_dt, m, p, os);
        case MICRO_DVR:
        default:
            return RMicroDVR_new().run(prior, Q_dt, m, p, os);
        }
    }

    template <class Model, class Step>
    myOptional_t<mp_ensamble_state_information_new>
    run(const mp_ensamble_state_information_new &prior,
        const Microscopic_Transition_step_double &Q_dt, Model &m, const Step &p,
        std::ostream &os, MICROR &alg) const {
        if constexpr (!variance)
            alg = MICRO_DMR;

        else
            alg = MICRO_DVR;

        const MICROR alg2 = alg;

        return run(prior, Q_dt, m, p, os, alg2);
    }

    template <class Model, class Step>
    myOptional_t<mp_ensamble_state_information_new>
    run(const mp_ensamble_state_information_new &prior,
        const Microscopic_Transition_step_double &Q_dt, Model &m, const Step &p,
        std::ostream & /*os*/) const {
        typedef myOptional_t<mp_ensamble_state_information_new> Op;
        auto y = p.y();
        double e = m.noise_variance(p.nsamples());
        auto k = prior.p().size();
        M_Matrix<double> Post(1, Q_dt.P().ncols(), 0);
        M_Matrix<double> logLik(Q_dt.P().nrows(), Q_dt.P().ncols(), 0);
        double partialLik = 0;
        [[maybe_unused]] double epartialLik = 0;
        double y_mean = 0;
        double y_var = e;
        double eplogL;
        double vplogL;
        double plogL;
        if (std::isnan(y)) {
            Post = prior.p() * Q_dt.P();
            y_mean = std::numeric_limits<double>::quiet_NaN();
            plogL = std::numeric_limits<double>::quiet_NaN();
            y_var = std::numeric_limits<double>::quiet_NaN();
            eplogL = std::numeric_limits<double>::quiet_NaN();
            vplogL = std::numeric_limits<double>::quiet_NaN();
        } else {
            double maxlogLik=-std::numeric_limits<double>::infinity();
            for (std::size_t j = 0; j < k; ++j) {
                for (std::size_t i = 0; i < k; ++i) {
                    if (prior.p()[i] * Q_dt.P()(i, j) > Q_dt.min_P()) {
                        logLik(i,j) = log_gaussian_function(y, Q_dt.gmean_ij()(i, j),
                                                             e + std::abs(Q_dt.gvar_ij()(i, j)));
                        if (logLik(i,j)>maxlogLik) maxlogLik=logLik(i,j);
                    }
                }
            }
            for (std::size_t j = 0; j < Q_dt.P().ncols(); ++j) {
                for (std::size_t i = 0; i < Q_dt.P().nrows(); ++i) {
                    if (prior.p()[i] * Q_dt.P()(i, j) > Q_dt.min_P()) {
                        Post[j] += prior.p()[i] * Q_dt.P()(i, j) *std::exp(logLik(i,j)-maxlogLik);
                        //   y_mean +=prior.P()[i] * Q_dt.P()(i, j) *Q_dt.gmean_ij()(i, j);
                        //   y_var +=prior.P()[i] * Q_dt.P()(i, j)
                        //   *std::abs(Q_dt.gvar_ij()(i, j));
                    }
                }
                partialLik += Post[j];
                y_mean += prior.p()[j] * Q_dt.gmean_i()[j];
                y_var += prior.p()[j] * Q_dt.gvar_i()[j];
            }
            for (std::size_t j = 0; j < Q_dt.P().ncols(); ++j) {
                Post[j] /= partialLik;
            }
            assert(partialLik>0);

            plogL = std::log(partialLik)+maxlogLik;
            eplogL = -0.5 * log(2 * PI * y_var) - 0.5;
            vplogL = 0.5;
            assert(std::isfinite(eplogL));
            std::cerr<<"y="<<y<<"\t"<<y_mean<<"\t"<<plogL<<"\n";
        }
        return Op(mp_ensamble_state_information_new::adjust(Q_dt.mi1(),std::move(Post), y_mean,
                                                      y_var, plogL, eplogL, vplogL,
                                                            Q_dt.min_P(), e,reduce_by_p()));
    }

    template <class Model, class Step>
    myOptional_t<mp_ensamble_state_information_new> start(Model &m, const Step &p,
                                                    double min_p) const {
        typedef myOptional_t<mp_ensamble_state_information_new> Op;
        auto x=p.begin()->x();

        auto Peq = m.Peq(x);
        if (!Peq)
            return Op(false, "fails to get Peq :" + Peq.error());
        else {
            auto P_cov = diag(Peq.value()) - quadraticForm_XTX(Peq.value());

                double nan = std::numeric_limits<double>::quiet_NaN();
            std::size_t N = m.AverageNumberOfChannels();
            double e = m.noise_variance(p.nsamples());

            Microscopic_description_new mi(N, Peq.value(), P_cov, max_Chi2());

                auto P = mi.convert_to_Microscopic_maxent(Peq.value(), P_cov);
            if (!P)
                return Op(false,"Maximum entropy failed: "+P.error());
            else

                    return Op(mp_ensamble_state_information_new::adjust(
                        std::move(mi), std::move(P).value(), nan, nan, nan, nan, 0.0, min_p, e,
                        reduce_by_p()));
            }
    }
    RMicro_R_new_() = default;
    RMicro_R_new_(double min_p__, double reduce_by_p__, double max_Chi2__, std::size_t order_number__):min_p_{min_p__},reduce_by_p_{reduce_by_p__},
                                                                                                        max_Chi2_{max_Chi2__}, order_number_{order_number__}{}
    double reduce_by_p()const {return reduce_by_p_;}
    double max_Chi2()const {return max_Chi2_;}
    std::size_t order_number()const {return order_number_;}
private:
    double min_p_;
    double reduce_by_p_;
    double max_Chi2_;
    std::size_t order_number_;
};


class MacroDR {
public:
  inline constexpr static auto const className = my_static_string("MacroDR");

  template <class Model, class Step>
  myOptional_t<mp_state_information> run(const mp_state_information &prior,
                                         Model &m, const Step &p,
                                         std::ostream &os) const {
    MACROR alg;
    return run(prior, m, p, os, alg);
  }

  template <class Model, class Step, class... Aux>
  myOptional_t<mp_state_information> run(const mp_state_information &prior,
                                         Model &m, const Step &p,
                                         std::ostream &os, Aux &... aux) const {
    typedef myOptional_t<mp_state_information> Op;
    auto Q_dto = m.get_P(p, 0);
    if (!Q_dto)
      return Op(false, "fails in auto Q_dt=m.get_P(p,0) :" + Q_dto.error());

    Markov_Transition_step_double Q_dt = std::move(Q_dto).value();
    //  Markov_Transition_step_double Q_dt2=m.get_P(p,20,eps);
    //  assert(class_Invariants<Markov_Transition_step_double>::test<true>(Q_dt,
    //  prior.tolerance()));
    //  assert(class_Invariants<Markov_Transition_step_double>::test<true>(Q_dt2,
    //  prior.tolerance())); areEqual(Q_dt,Q_dt2,1e-6);

    double e = m.noise_variance(p.nsamples());
    //
    double N = m.AverageNumberOfChannels();
    M_Matrix<double> u(prior.P_mean().size(), 1, 1.0);

    auto SmD = prior.P_cov() - diag(prior.P_mean());

    double gSg =
        (TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.gmean_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_ij(), Q_dt.gmean_ij()) * u))
            .getvalue();

    double sSg =
        (TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.gmean_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_var_ij(), Q_dt.gmean_ij()) * u))
            .getvalue();
    double sSs =
        (TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.gvar_i()).getvalue() +
        (prior.P_mean() * (elemMult(Q_dt.gtotal_var_ij(), Q_dt.gvar_ij()) * u))
            .getvalue();
    //   auto mu_n=;
    auto sS = TranspMult(Q_dt.gvar_i(), SmD) * Q_dt.P() +
              prior.P_mean() * Q_dt.gtotal_var_ij();
    auto gS = TranspMult(Q_dt.gmean_i(), SmD) * Q_dt.P() +
              prior.P_mean() * Q_dt.gtotal_ij();

    double ms = (prior.P_mean() * Q_dt.gvar_i()).getvalue();
    double delta_emu = std::max(sqr(ms + e / N) - 2.0 / N * sSs, 0.0);
    double ms0 = (ms - e / N) / 2 + std::sqrt(delta_emu) / 2;
    // std::cerr<<"(ms0==ms-sSs/(2*(e+N*ms0)))="<<(ms0==ms-sSs/(2*(e+N*ms0)))<<
    // "ms0="<<ms0;
    // assert((are_Equal<true,double>(std::numeric_limits<double>::epsilon()*1000).test(ms0,ms-sSs/(2*(e+N*ms0)))));

    auto e_mu = e + N * ms0;
    auto y_mean =
        N * (prior.P_mean() * Q_dt.gmean_i()).getvalue() - N * 0.5 / e_mu * sSg;
    auto zeta = N / (2 * sqr(e_mu) + N * sSs);
    auto y_var = std::max(e_mu + N * gSg - N * zeta * sqr(sSg), e);
    auto y = p.y();
    if (std::isnan(y)) {
      double plogL = std::numeric_limits<double>::quiet_NaN();
      double eplogL = std::numeric_limits<double>::quiet_NaN();
      double vplogL = 0.0;
      auto P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P());
      auto P_mean = prior.P_mean() * Q_dt.P();
      P__cov += diag(P_mean);
      // std::cerr<<"\nPcov nana corr\n"<<P__cov<<"\nP_mean nana
      // corr\n"<<P_mean<<"\nQ.P \n"<<Q_dt.P();
      auto test = mp_state_information::test(P_mean, P__cov, tolerance_);
      if (test.has_value())
        return Op(mp_state_information::adjust(
            std::move(P_mean), std::move(P__cov), y_mean, y_var, plogL, eplogL,
            vplogL, Q_dt.min_P(), e));
      else
        return Op(false, "fails at intertrace prediction!!: " + test.error());
    } else {
      auto dy = y - y_mean;
      auto chi = dy / y_var;
      auto P_mean = prior.P_mean() * Q_dt.P() + chi * gS -
                    (chi * zeta * sSg + 0.5 / sqr(e_mu)) * sS;

      auto P__cov =
          quadraticForm_BT_A_B(SmD, Q_dt.P()) +
          diag(prior.P_mean() * Q_dt.P()) -
          (zeta + N / y_var * sqr(zeta * sSg)) * quadraticForm_XTX(sS) +
          (2.0 * N / y_var * zeta * sSg) * TransposeSum(TranspMult(sS, gS)) -
          (N / y_var) * quadraticForm_XTX(gS);

      auto chi2 = dy * chi;

      double plogL;
      if (y_var > 0)
        plogL = -0.5 * log(2 * PI * y_var) - 0.5 * chi2;
      else
        plogL = std::numeric_limits<double>::infinity();

      double eplogL =
          -0.5 * log(2 * PI * y_var) - 0.5; // e_mu+N*gSg"-N*zeta*sqr(sSg)"
      double chilogL = (eplogL - plogL) / std::sqrt(0.5);
      if (sqr(chilogL) > 10) {

        double dt = p.nsamples() / m.fs();
        os << "\n--------------------------------------------------------------"
              "-----------------\n";
        os << "\n---------------------------chi_logL "
              "big----------------------------------------------------\n";
        os << "\n--------------------------------------------------------------"
              "-----------------\n";

        os << "\n t= " << p.x().t() << " nsamples=" << p.nsamples()
           << " dt=" << dt << " x=" << p.x().x() << " y=" << p.y()
           << " y_mean=" << y_mean;
        os << " chi2=" << chi2 << " y_var" << y_var;
        os << "N*(prior.P_mean()*Q_dt.gmean_i()).getvalue()-N*0.5/e_mu*sSg";
        os << "N*(prior.P_mean()*Q_dt.gmean_i()).getvalue()="
           << N * (prior.P_mean() * Q_dt.gmean_i()).getvalue() << "\n";
        os << "N*0.5/e_mu*sSg=" << N * 0.5 / e_mu * sSg << "\n";
        os << "0.5/e_mu=" << 0.5 / e_mu << "\n";
        os << "sSg=" << sSg << "\n";

        os << "\nprior.P_mean()\n" << prior.P_mean();

        //         os<<"\n Qdt\n"<<Q_dt;
      }
      double mg = (prior.P_mean() * Q_dt.gmean_i()).getvalue();
      double g_max = max(Q_dt.g());
      double g_min = min(Q_dt.g());
      double g_range = g_max - g_min;
      // double Neff=p_bi*q_bi/gSg*sqr(g_range);
      //    os<<" mg="<<mg;
      //     os<<" g_max="<<g_max<<" g_min="<<g_min;
      //    os<<"\t p_bi="<<p_bi<<" q_bi="<<q_bi<<"
      //    min(Nq,Np)="<<std::min(N*q_bi,N*p_bi)<<" Neff="<<Neff<<"
      //    NeffNp="<<std::min(N*q_bi,N*p_bi)*Neff;
      //   os<<" e="<<e<<" ms="<<ms<<" ms0="<<ms0<<" e_mu="<<e_mu<<"
      //   N*gSg="<<N*gSg; os<<" -N*zeta*sqr(sSg)="<<-N*zeta*sqr(sSg)<< "
      //   plogL="<<plogL<< " eplogL="<<eplogL<<"dif="<<plogL-eplogL<<"
      //   sSs="<<sSs;
      //  os<<"\nPcov \n"<<P__cov;
      //  os<<"\nP_mean \n"<<P_mean;
      //  os<<"\n Qdt\n"<<Q_dt;
      if (auto test = variance_value::test<true>(y_var, e, e * tolerance());
          !test.has_value()) {
        std::stringstream ss;
        ss << "\n step=" << p << " "
           << " y=" << y << " ymean=" << y_mean << " yvar=" << y_var
           << " e=" << e << " e_mu" << e_mu << " N*gSg" << N * gSg
           << " -N*zeta*sSg=" << -N * zeta * sqr(sSg);
        ss << "\n ms0=(ms-e/N)/2+std::sqrt(delta_emu)/2;\nms0=" << ms0
           << "\n ms=" << ms << "\n e/N/2=" << e / N;

        ss << "\n std::sqrt(delta_emu)/2=" << std::sqrt(delta_emu) / 2;
        ss << "\ndelta_emu=sqr(ms+e/N)-2.0/N*sSs \nsqr(ms+e/N)="
           << sqr(ms + e / N) << " 2.0/N*sSs=" << 2.0 / N * sSs
           << " sSs=" << sSs;
        ss << "\nP_mean \n" << P_mean;
        ss << " -N*zeta*sqr(sSg)=" << -N * zeta * sqr(sSg) << " plogL=" << plogL
           << " eplogL=" << eplogL << "dif=" << plogL - eplogL
           << " sSs=" << sSs;
        ss << "\nPcov \n" << P__cov;
        // ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;

        return Op(false, "error in variance!! " + test.error() + ss.str());
      }

      std::pair<Op_void, Op_void> test_Binomial;
      if constexpr (sizeof...(Aux) > 0) {
        if constexpr ((std::is_const_v<Aux> && ...)) {
          (test_Binomial = ... = aux);
          if (test_Binomial.first.has_value())
            os << " use 1 \t";
          else
            os << " use 0 \t" << test_Binomial.first.error() << "\t";
        } else {
          auto p_bi = (g_max - mg) / g_range;
          auto q_bi = (mg - g_min) / g_range;

          test_Binomial = mp_state_information::is_Binomial_Approximation_valid(
              N, p_bi, q_bi, 5, 0);
          (aux = ... = test_Binomial);
          if (test_Binomial.first.has_value())
            os << " set 1\t";
          else
            os << " set 0 \t" << test_Binomial.first.error() << "\t";
        }
      } else

      {
        auto p_bi = (g_max - mg) / g_range;
        auto q_bi = (mg - g_min) / g_range;
        test_Binomial = mp_state_information::is_Binomial_Approximation_valid(
            N, p_bi, q_bi, 5, 0);
        if (test_Binomial.first.has_value())
          os << " tes 1 \t";
        else
          os << " tes 0 \t" << test_Binomial.first.error() << "\t";
      }

      if (!test_Binomial.first.has_value()) {
        std::cerr << " \n\n----test----\n"
                  << test_Binomial.first.error() << "\n";
        P__cov = quadraticForm_BT_A_B(SmD, Q_dt.P());
        P_mean = prior.P_mean() * Q_dt.P();
        P__cov += diag(P_mean);
        //  std::cerr<<"\n SmD \n"<<SmD;

        //  std::cerr<<"\nPcov corr\n"<<P__cov<<"\nP_mean
        //  corr\n"<<P_mean<<"\nQ.P \n"<<Q_dt.P();
      }
      auto test = mp_state_information::test(P_mean, P__cov, y_mean, y_var,
                                             plogL, eplogL, e, tolerance());
      if (!test) {
        std::stringstream ss;
        ss << "\n step=" << p << " "
           << " y=" << y << " ymean=" << y_mean << " yvar=" << y_var
           << " e=" << e << " e_mu" << e_mu << " N*gSg" << N * gSg
           << " -N*zeta*sSg=" << -N * zeta * sqr(sSg);
        ss << "\n ms0=(ms-e/N)/2+std::sqrt(delta_emu)/2;\nms0=" << ms0
           << "\n ms=" << ms << "\n e/N/2=" << e / N;

        ss << "\n std::sqrt(delta_emu)/2=" << std::sqrt(delta_emu) / 2;
        ss << "\ndelta_emu=sqr(ms+e/N)-2.0/N*sSs \nsqr(ms+e/N)="
           << sqr(ms + e / N) << " 2.0/N*sSs=" << 2.0 / N * sSs
           << " sSs=" << sSs;
        ss << "\nP_mean \n" << P_mean;
        ss << " -N*zeta*sqr(sSg)=" << -N * zeta * sqr(sSg) << " plogL=" << plogL
           << " eplogL=" << eplogL << "dif=" << plogL - eplogL
           << " sSs=" << sSs;
        ss << "\nPcov \n" << P__cov;
        // ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;

        return Op(false,
                  "\nfails in trace!!!; error=" + test.error() + ss.str());
      } else
        return Op(mp_state_information::adjust(std::move(P_mean),
                                               std::move(P__cov), y_mean, y_var,
                                               plogL, eplogL, Q_dt.min_P(), e));
    }
  }

  template <class Model, class Step>
  mp_state_information run_old(const mp_state_information &prior, Model &m,
                               const Step &p) const {
    // double dt=Y.dt();
    //  double x=Y.x();
    // double y=Y.y();

    double eps = std::numeric_limits<double>::epsilon();
    Markov_Transition_step_double Q_dt = m.get_P(p, 0, eps);
    Markov_Transition_step_double Q_dt2 = m.get_P(p, 10, eps);
    Markov_Transition_step_double Q_dt3 = m.get_P(p, 20, eps);

    bool test, test2, test3;
    if (class_Invariants<Markov_Transition_step_double>::test<true>(Q_dt))
      test = true;
    else
      test = false;
    if (class_Invariants<Markov_Transition_step_double>::test<true>(Q_dt2))
      test2 = true;
    else
      test2 = false;

    if (class_Invariants<Markov_Transition_step_double>::test<true>(Q_dt3))
      test3 = true;
    else
      test3 = false;

    M_Matrix<double> Sm = prior.P_cov() - diag(prior.P_mean());

    M_Matrix<double> sSm = TranspMult(Q_dt.gvar_i(), Sm);
    M_Matrix<double> gSm = TranspMult(Q_dt.gmean_i(), Sm);

    double sms = 0;
    double smg = 0;
    double gmg = 0;
    std::size_t k_u = m.nstates();
    for (std::size_t i = 0; i < k_u; i++)
      for (std::size_t j = 0; j < k_u; j++) {
        sms += Q_dt.gvar_ij()(i, j) * prior.P_mean()[i] *
               Q_dt.gtotal_var_ij()(i, j);
        smg +=
            Q_dt.gvar_ij()(i, j) * prior.P_mean()[i] * Q_dt.gtotal_ij()(i, j);
        gmg +=
            Q_dt.gmean_ij()(i, j) * prior.P_mean()[i] * Q_dt.gtotal_ij()(i, j);
      };

    double gSg = (gSm * Q_dt.gmean_i())[0] + gmg;
    double sSg = (sSm * Q_dt.gmean_i())[0] + smg;
    double sSs = (sSm * Q_dt.gvar_i())[0] + sms;

    M_Matrix<double> gS = gSm * Q_dt.P() + prior.P_mean() * Q_dt.gtotal_ij();
    M_Matrix<double> sS =
        sSm * Q_dt.P() + prior.P_mean() * Q_dt.gtotal_var_ij();
    double N = m.AverageNumberOfChannels();
    double smean = N * (prior.P_mean() * Q_dt.gvar_i())[0];
    double e = m.noise_variance(p.nsamples());

    double sA2 = N * (prior.P_mean() * Q_dt.gvar_i())[0] +
                 m.noise_variance(p.nsamples());

    double sB4 = 2 * sA2 * sA2 - N * sSs;

    double y_var = sA2 + N * gSg + (N * N * sSg * sSg) / sB4;
    // sometimes sB4 might be negative, when P_cov_M is greter than one.

    if (y_var < 0) {
      std::cerr << " \n******** y_var_d negative!!!!!***********\n";
      std::cerr << "\n sA2+N*gSg+(N*N*sSg*sSg)/sB4 \n";
      std::cerr << "\t smean=" << smean << "\t e=" << e;
      std::cerr << "\n sA2=  " << sA2;
      std::cerr << "\n sB4=  " << sB4;

      std::cerr << "\n gSg =" << gSg;
      std::cerr << "\n sSs =" << sSs;

      std::cerr << "\n sB4=2*sA2*sA2-N*sSs  sB4 = " << sB4 << "\n";

      //        std::cerr<<*this;
      //        std::cerr<<this->model();
      //        //press_any_key_to_continue();
    }

    // auto y_std=std::sqrt(y_var);

    double y_mean = N * (prior.P_mean() * Q_dt.gmean_i())[0] - sA2 / sB4 * sSg;

    // product of S and g, used several places
    auto y = p.y();

    if (std::isnan(y)) {
      double plogL = std::numeric_limits<double>::quiet_NaN();
      double eplogL = std::numeric_limits<double>::quiet_NaN();
      double vplogL = 0;
      // auto chi2=std::numeric_limits<double>::quiet_NaN();
      /* auto
       * P_cov=TranspMult(Q_dt.P(),(prior.P_cov()-diag(prior.P_mean()))*Q_dt.P());*/
      auto P_cov = quadraticForm_BT_A_B((prior.P_cov() - diag(prior.P_mean())),
                                        Q_dt.P());

      auto P_mean = prior.P_mean() * Q_dt.P();
      P_cov += diag(prior.P_mean());
      return mp_state_information(std::move(P_mean), std::move(P_cov), y_mean,
                                  y_var, plogL, eplogL, vplogL);

    } else {
      M_Matrix<double> P_mean = prior.P_mean() * Q_dt.P();
      double dy = y - y_mean;
      P_mean += sS * (N * sSg / sB4 * (p.y() - y_mean) / y_var - sA2 / sB4) +
                gS * ((p.y() - y_mean) / y_var);
      if (!(P_mean >= 0.0)) {
        double summ = 0;
        for (std::size_t i = 0; i < size(P_mean); i++) {
          if (P_mean[i] < 0)
            P_mean[i] = 0;
          else if (P_mean[i] > 1)
            P_mean[i] = 1;
          summ += P_mean[i];
        }
        P_mean /= summ;
      }
      auto chi2 = dy * dy / y_var;

      double plogL;
      if (y_var > 0)
        plogL = -0.5 * log(2 * PI * y_var) - 0.5 * chi2;
      else
        plogL = std::numeric_limits<double>::infinity();

      double eplogL = -0.5 * (1.0 + log(2 * PI * y_var));
      double vplogL = 0.5;

      auto P_cov_ = TranspMult(Q_dt.P(), Sm) * Q_dt.P() +
                    diag(prior.P_mean() * Q_dt.P()) -
                    TranspMult(gS, gS) * (N / y_var) -
                    (TranspMult(sS, gS) + TranspMult(gS, sS)) *
                        (N * N * sSg / y_var / sB4) +
                    TranspMult(sS, sS) *
                        (N / sB4 - N / y_var * N * N * sSg * sSg / sB4 / sB4);

      auto P_cov =
          quadraticForm_BT_A_B(Sm, Q_dt.P()) + diag(prior.P_mean() * Q_dt.P()) -
          quadraticForm_XTX(gS) * (N / y_var) -
          TransposeSum(TranspMult(sS, gS)) * (N * N * sSg / y_var / sB4) +
          quadraticForm_XTX(sS) *
              (N / sB4 - N / y_var * N * N * sSg * sSg / sB4 / sB4);

      // assert((are_Equal<true,
      // M_Matrix<double>>(Q_dt.min_P(),Q_dt.min_P()).test_prod(P_cov_,P_cov,std::cerr)));

      if (!(diag(P_cov) >= 0.0) || !(diag(P_cov) <= 1.0)) {
        for (std::size_t i = 0; i < ncols(P_cov); i++) {
          if (P_cov(i, i) < 0.0)
            P_cov(i, i) = std::abs(P_cov(i, i));
          else if (P_cov(i, i) > 1.0)
            P_cov(i, i) = 1.0;
        }
      }
      return mp_state_information(std::move(P_mean), std::move(P_cov), y_mean,
                                  y_var, plogL, eplogL, vplogL);
    }
  }

  template <class Model, class Step>
  myOptional_t<mp_state_information> start(Model &m, const Step &p,
                                           double min_p) const {
    typedef myOptional_t<mp_state_information> Op;

    auto P_mean = m.Peq(p.begin()->x());
    //  std::cerr<<"gslreijgsorjgps INIT!!!"<<P_mean.value();
    if (!P_mean)
      return Op(false, "fails to get Peq :" + P_mean.error());
    else {
      auto P_cov = diag(P_mean.value()) - quadraticForm_XTX(P_mean.value());
      //     std::cerr<<"gslreijgsorjgps INIT!!!"<<P_cov;
      double nan = std::numeric_limits<double>::quiet_NaN();
      return Op(mp_state_information::adjust(std::move(P_mean).value(),
                                             std::move(P_cov), nan, nan, nan,
                                             nan, min_p, 0));
    }
  }
  double tolerance() const { return tolerance_; }

  MacroDR(double tolerance) : tolerance_{tolerance} {}
  MacroDR() = default;

private:
  double tolerance_ = 1e-5;
};


template<template<class...> class Tr>
struct logLikelihood_function_;
typedef logLikelihood_function_<C> logLikelihood_function;

template<template<class...> class Tr>
struct logLikelihood_function_ {
public:
    template<class T> using Tr_t=typename Tr<T>::type;


    template <class Experiment> Tr_t<double> operator()(const Experiment &) const
    {
        return Tr_t<double> (0.0);
    }
  template <class  mp_state_information>
  void operator()(const mp_state_information &mp, Tr_t<double> &logLsum,
                  std::size_t &i) const {
      if (std::isfinite(center(mp.plogL()))) {
      logLsum += mp.plogL();
      ++i;
    }
  }
};


struct partialLogLikelihood_function {
  template <class Experiment>
  std::vector<double> operator()(const Experiment &e) const {
    std::size_t n = e.num_measurements();

    return std::vector<double>(n);
  }

  template <class mp_state_information>
  void operator()(const mp_state_information &mp, std::vector<double> &v,
                  std::size_t &i) const {
    if (std::isfinite(mp.plogL())) {
      v[i] = mp.plogL();
      ++i;
    }
  }
};

struct partialDistribution_function {
  template <class Experiment> auto operator()(const Experiment &e) const {
    std::size_t n = e.num_measurements();

    return std::vector<Normal_Distribution<double>>(n);
  }

  // template<class mp_state_information>
  void operator()(const mp_state_information &mp,
                  std::vector<Normal_Distribution<double>> &v,
                  std::size_t &i) const {
    if (std::isfinite(mp.plogL())) {
      v[i] = Normal_Distribution<double>(mp.y_mean(), mp.y_var());
      ++i;
    } else // hack to avoid including intervals in the Jacobian calculation
    {
      v[i] = Normal_Distribution<double>(0, 1);
      ++i;
    }
  }
};

template <bool calc, class aux = markov::MACROR>
struct partialDistribution_function_aux {
  template <class Experiment>
  auto operator()(const Experiment &e, const std::vector<aux> &) const {
    std::size_t n = e.num_measurements();

    return std::vector<Normal_Distribution<double>>(n);
  }

  template <class Experiment>
  auto operator()(const Experiment &e, std::vector<aux> &a) const {
    std::size_t n = e.num_measurements();

    if constexpr (calc) {
      a.resize(n);
      return std::pair(std::vector<Normal_Distribution<double>>(n),
                       std::vector<aux>(n));
    } else

      return std::vector<Normal_Distribution<double>>(n);
  }

  // template<class mp_state_information>
  void operator()(
      const mp_state_information &mp,
      std::pair<std::vector<Normal_Distribution<double>>, std::vector<aux>> &v,
      std::size_t &i, aux &test) const {
    if (std::isfinite(mp.plogL())) {
      v.first[i] = Normal_Distribution<double>(mp.y_mean(), mp.y_var());
      v.second[i] = test;
      ++i;
    } else // hack to avoid including intervals in the Jacobian calculation
    {
      v.first[i] = Normal_Distribution<double>(0, 1);
      v.second[i] = test;
      ++i;
    }
  }
  void operator()(const mp_state_information &mp,
                  std::vector<Normal_Distribution<double>> &v, std::size_t &i,
                  const aux &) const {
    if (std::isfinite(mp.plogL())) {
      v[i] = Normal_Distribution<double>(mp.y_mean(), mp.y_var());
      ++i;
    } else // hack to avoid including intervals in the Jacobian calculation
    {
      v[i] = Normal_Distribution<double>(0, 1);
      ++i;
    }
  }
};

struct partialDistribution_function_mp {

  template <class Experiment>
  auto operator()(
      const Experiment &e,
      std::vector<std::tuple<markov::MACROR, markov::mp_state_information>> &v)
      const {
    std::size_t n = e.num_measurements();
    v.resize(n);
    return std::make_pair(
        std::vector<Normal_Distribution<double>>(n),
        std::vector<std::tuple<markov::MACROR, markov::mp_state_information>>(
            n));
  }
  template <class Experiment>
  auto operator()(
      const Experiment &e,
      const std::vector<
          std::tuple<markov::MACROR, markov::mp_state_information>> &) const {
    std::size_t n = e.num_measurements();
    return std::vector<Normal_Distribution<double>>(n);
  }

  // template<class mp_state_information>
  void operator()(
      const mp_state_information &mp,
      std::pair<
          std::vector<Normal_Distribution<double>>,
          std::vector<std::tuple<markov::MACROR, markov::mp_state_information>>>
          &v,
      std::size_t &i,
      std::tuple<markov::MACROR, markov::mp_state_information> &test) const {
    if (std::isfinite(mp.plogL())) {
      v.first[i] = Normal_Distribution<double>(mp.y_mean(), mp.y_var());
      std::get<0>(v.second[i]) = std::get<0>(test);
      std::get<1>(v.second[i]) = mp;
      ++i;
    } else // hack to avoid including intervals in the Jacobian calculation
    {
      v.first[i] = Normal_Distribution<double>(0, 1);
      std::get<0>(v.second[i]) = std::get<0>(test);
      std::get<1>(v.second[i]) = mp;
      ++i;
    }
  }
  void operator()(
      const mp_state_information &mp,
      std::vector<Normal_Distribution<double>> &v, std::size_t &i,
      const std::tuple<markov::MACROR, markov::mp_state_information> &) const {
    if (std::isfinite(mp.plogL())) {
      v[i] = Normal_Distribution<double>(mp.y_mean(), mp.y_var());
      ++i;
    } else // hack to avoid including intervals in the Jacobian calculation
    {
      v[i] = Normal_Distribution<double>(0, 1);
      ++i;
    }
  }
};



template <class F, class Calculator,class MacroDR, class Model,
          template <class, class> class Experiment, class Point, class Measure,
          class... Aux>
auto logLikelihood_experiment_calculation_new(const F &f, const MacroDR &a,
                                          Calculator& calc,
                                          const Model &m,
                                          const Experiment<Point, Measure> &e,
                                          std::ostream &os, Aux &... aux) {
    auto out = f(e, aux...);
    typedef myOptional_t<std::decay_t<decltype(out)>> Op;

    auto first_step = *e.begin_begin();
    auto prior = a.start(calc,m, first_step, m.min_P());
       if (!prior.has_value())
        return Op(false, "calculation interrupted at start :" + prior.error());
    std::size_t i = 0;
    auto end=e.end_end();
    --end;
    for (auto it = e.begin_begin(); it != end; ++it) {
        auto post = a.run(prior.value(),calc, m, *it, os, aux[i]...);
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





template <class F, class MacroDR, class Model,
          template <class, class> class Experiment, class Point, class Measure,
          class... Aux>
auto logLikelihood_experiment_calculation(const F &f, const MacroDR &a,
                                          Model &m,
                                          const Experiment<Point, Measure> &e,
                                          std::ostream &os, Aux &... aux) {
  auto out = f(e, aux...);
  typedef myOptional_t<std::decay_t<decltype(out)>> Op;

  auto first_step = *e.begin_begin();
  auto prior = a.start(m, first_step, m.min_P());
  if (!prior.has_value())
    return Op(false, "calculation interrupted at start :" + prior.error());
  std::size_t i = 0;
  auto end=e.end_end();
  --end;
  for (auto it = e.begin_begin(); it != end; ++it) {
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


template <class MacroDR,class Calculator, class Model, class Experiment, std::enable_if_t<!is_Derivative<Model>::value,int> =0>
auto logLikelihood_new(const MacroDR &a, Calculator& calc,const Model &m,
                                   const Experiment e, std::ostream &os) {

    typename MacroDR::template Tr_t<logLikelihood_function> lik;

    return logLikelihood_experiment_calculation_new(lik, a,calc, m, e,
                                                os);
}

template <class MacroDR,class Calculator, class Model, class Experiment, std::enable_if_t<is_Derivative<Model>::value,int> =0>
auto logLikelihood_new(const MacroDR &a, Calculator& calc,const Model &m,
                       const Experiment e, std::ostream &os) {

    typename MacroDR::template Tr_t<logLikelihood_function> lik(m.Qrun);

    return logLikelihood_experiment_calculation_new(lik, a,calc, m, e,
                                                    os);
}



template <class MacroDR, class Model, class Experiment>
myOptional_t<double> logLikelihood(const MacroDR &a, Model &m,
                                   const Experiment e, std::ostream &os) {
  return logLikelihood_experiment_calculation(logLikelihood_function(), a, m, e,
                                              os);
}

template <class MacroDR, class Model, class Experiment>
myOptional_t<std::vector<double>>
partialLogLikelihood(const MacroDR &a, Model &m, const Experiment e,
                     std::ostream &os) {
  return logLikelihood_experiment_calculation(partialLogLikelihood_function(),
                                              a, m, e, os);
}

template <class MacroDR, class Model, class Experiment>
myOptional_t<std::vector<Normal_Distribution<double>>>
partialDistribution(const MacroDR &a, Model &m, const Experiment e) {
  return logLikelihood_experiment_calculation(partialDistribution_function(), a,
                                              m, e, std::cerr);
}

template <class MacroDR, class Model, class Experiment>
myOptional_t<std::pair<std::vector<Normal_Distribution<double>>,
                       std::vector<markov::MACROR>>>
partialDistribution_aux(const MacroDR &a, Model &m, const Experiment e,
                        std::ostream &os) {
  std::vector<markov::MACROR> aux;
  return logLikelihood_experiment_calculation(
      partialDistribution_function_aux<true>(), a, m, e, os, aux);
}

template <class MacroDR, class Model, class Experiment>
myOptional_t<std::vector<Normal_Distribution<double>>>
partialDistribution_aux(const MacroDR &a, Model &m, const Experiment e,
                        std::ostream &os,
                        const std::vector<markov::MACROR> &aux) {
  return logLikelihood_experiment_calculation(
      partialDistribution_function_aux<false>(), a, m, e, os, aux);
}

template <class MacroDR, class Model, class Experiment>
myOptional_t<std::pair<
    std::vector<Normal_Distribution<double>>,
    std::vector<std::tuple<markov::MACROR, markov::mp_state_information>>>>
partialDistribution_aux_mp(const MacroDR &a, Model &m, const Experiment e,
                           std::ostream &os) {
  std::vector<std::tuple<markov::MACROR, markov::mp_state_information>> aux;
  return logLikelihood_experiment_calculation(partialDistribution_function_mp(),
                                              a, m, e, os, aux);
}

template <class MacroDR, class Model, class Experiment>
myOptional_t<std::vector<Normal_Distribution<double>>>
partialDistribution_aux_mp(
    const MacroDR &a, Model &m, const Experiment e, std::ostream &os,
    const std::vector<std::tuple<markov::MACROR, markov::mp_state_information>>
        &aux) {
  return logLikelihood_experiment_calculation(partialDistribution_function_mp(),
                                              a, m, e, os, aux);
}

template <class Y>
struct measure_likelihood : public experiment::measure_just_y<Y> {
  double y_mean;
  double y_var;
  double plogL;
  double eplogL;

  measure_likelihood() = default;

  template <class mp_state_information>
  measure_likelihood(const mp_state_information &mp)
      : y_mean{mp.y_mean()}, y_var{mp.y_var()}, plogL{mp.plogL()},
        eplogL{mp.eplogL()} {}

  measure_likelihood(const std::tuple<double, double, double, double> &data)
      : y_mean{std::get<0>(data)}, y_var{std::get<1>(data)},
        plogL{std::get<2>(data)}, eplogL{std::get<3>(data)} {}

  template <class DataFrame> static void insert_col(DataFrame &d) {
    d.insert_column("ymean", C<double>{});
    d.insert_column("y_var", C<double>{});
    d.insert_column("plogL", C<double>{});
    d.insert_column("eplogL", C<double>{});
  }
  std::tuple<double, double, double, double> data_row() const {
    return {y_mean, y_var, plogL, eplogL};
  }
};

template <class measure_likelihood>
struct partialLogLikelihood_monitor_function {
  template <class Experiment>
  std::vector<measure_likelihood> operator()(const Experiment &e) const {
    std::size_t n = e.end_end() - e.begin_begin();
    return std::vector<measure_likelihood>(n);
  }

  template <class mp_state_information>
  void operator()(const mp_state_information &mp,
                  std::vector<measure_likelihood> &v, std::size_t &i) const {
    v[i] = measure_likelihood(mp);
    ++i;
  }
};

template <class MacroDR, class Model, class Experiment>
auto monitorLikelihood(const MacroDR &a, Model &m, const Experiment e,
                       std::ostream &os) {
  typedef myOptional_t<experiment::basic_Experiment<
      experiment::point<double, double>, measure_likelihood<double>>>
      Op;
  auto v = logLikelihood_experiment_calculation(
      partialLogLikelihood_monitor_function<measure_likelihood<double>>(), a, m,
      e, os);
  if (!v.has_value())
    return Op(false,
              "logLikelihood_experiment_calculation fails :" + v.error());
  else
    return Op(experiment::basic_Experiment<experiment::point<double, double>,
                                           measure_likelihood<double>>(
        e, std::move(v).value()));
}

template <class MacroDR, class Calculator,class Model, class Experiment>
auto monitorLikelihood_new(const MacroDR &a, Calculator& calc, const Model &m, const Experiment e,
                       std::ostream &os) {
    typedef myOptional_t<experiment::basic_Experiment<
        experiment::point<double, double>, measure_likelihood<double>>>
        Op;
    auto v = logLikelihood_experiment_calculation_new(
        partialLogLikelihood_monitor_function<measure_likelihood<double>>(), a,calc, m,
        e, os);
    if (!v.has_value())
        return Op(false,
                  "logLikelihood_experiment_calculation fails :" + v.error());
    else
        return Op(experiment::basic_Experiment<experiment::point<double, double>,
                                               measure_likelihood<double>>(
            e, std::move(v).value()));
}


} // namespace markov

#endif // LIKELIHOOD_MARKOV_PROCESS_H
