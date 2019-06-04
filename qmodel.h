#ifndef QMODEL_H
#define QMODEL_H

#include "Matrix.h"
#include "maximum_entropy.h"
#include "myDistributions.h"
#include "mydynamicfunctions.h"
#include "myfields.h"
#include "mymath.h"
#include "myparameters.h"
#include "mytests.h"
#include "mytypetraits.h"
#include "matrixerror.h"

#include <cmath>
#include <map>
#include <set>
#include <vector>

struct i_micro_ensemble {
  constexpr static auto const title = my_static_string("i_micro_ensemble");
};

class Model_Parameter_label : public Parameter_label {
public:
  using Label::legal_chars;
  using Parameter_label::Parameter_label;

  constexpr static auto className = my_static_string("Model_Parameter_label");
};

class Number_of_Channels_Parameter_label : public Model_Parameter_label {
public:
  constexpr static auto className =
      my_static_string("Number_of_Channels_Parameter_label");
  Number_of_Channels_Parameter_label()
      : Model_Parameter_label("Number_of_Channels") {}
};

class Noise_Parameter_label : public Model_Parameter_label {
public:
  using Model_Parameter_label::Model_Parameter_label;

  constexpr static auto className = my_static_string("Noise_Parameter_label");
};

class gaussian_noise_Parameter_label : public Noise_Parameter_label {
public:
  constexpr static auto className =
      my_static_string("gaussian_noise_Parameter_label");
  gaussian_noise_Parameter_label() : Noise_Parameter_label("gaussian_noise") {}
};

class rate_Parameter_label : public Model_Parameter_label {
public:
  using Model_Parameter_label::Model_Parameter_label;
  constexpr static auto className = my_static_string("rate_Parameter_label");
};

class Equilibrium_Parameter_label : public Model_Parameter_label {
public:
  using Model_Parameter_label::Model_Parameter_label;
  constexpr static auto className =
      my_static_string("Equilibrium_Parameter_label");
};

class Time_Constant_Parameter_label : public Model_Parameter_label {
public:
  using Model_Parameter_label::Model_Parameter_label;
  constexpr static auto className =
      my_static_string("Time_Constant_Parameter_label");
};

class Coupling_factor_Parameter_label : public Model_Parameter_label {
public:
  using Model_Parameter_label::Model_Parameter_label;
  constexpr static auto className =
      my_static_string("Coupling_factor_Parameter_label");
};

class Coupling_coefficient_Parameter_label : public Model_Parameter_label {
public:
  using Model_Parameter_label::Model_Parameter_label;
  constexpr static auto className =
      my_static_string("Coupling_coefficient_Parameter_label");
};

class Conductance_Parameter_label : public Model_Parameter_label {
public:
  using Label::legal_chars;

  using Model_Parameter_label::Model_Parameter_label;
  constexpr static auto className =
      my_static_string("Conductance_Parameter_label");
};

class Conformational_change;
class Conformational_change_label : public Label<Conformational_change> {
public:
  using Label::Label;
  using Label::legal_chars;
  constexpr static auto className =
      my_static_string("Conformational_change_label");
};

class Conformational_change {
  int change_in_agonist_;
  int change_in_conductance_;
  Conformational_change_label label_;
  rate_Parameter_label par_on_;
  rate_Parameter_label par_off_;
  /*   Equilibrium_Parameter_label par_Eq_;
Time_Constant_Parameter_label par_tc_; */
public:
  typedef Conformational_change self_type;
  Conformational_change_label label() const { return label_; }
  const rate_Parameter_label &par_on() const { return par_on_; }
  const rate_Parameter_label &par_off() const { return par_off_; }
  /*  const Equilibrium_Parameter_label&  par_Eq()const {return par_Eq_;};
const Time_Constant_Parameter_label& par_tc()const {return par_tc_;}*/

  int change_in_agonist() const { return change_in_agonist_; }
  int change_in_conductance() const { return change_in_conductance_; }
  Conformational_change() = default;
  Conformational_change(int _change_in_agonist,
                          int _change_in_conductance,
                          Conformational_change_label _label,
                          rate_Parameter_label _par_on,
                          rate_Parameter_label _par_off /*,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                Equilibrium_Parameter_label&& _par_Eq,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                Time_Constant_Parameter_label&& _par_tc*/):
                              change_in_agonist_{_change_in_agonist},change_in_conductance_{_change_in_conductance},
                              label_{std::move(_label)},par_on_{std::move(_par_on)},par_off_{std::move(_par_off)}/*,par_Eq_{std::move(_par_Eq)},par_tc_{std::move(_par_tc)}*/{}

  constexpr static auto className = my_static_string("Conformational_change");

  static auto get_constructor_fields() {
    return std::make_tuple(grammar::field(C<self_type>{},"change_in_agonist", &self_type::change_in_agonist),
                               grammar::field(C<self_type>{},"change_in_conductance", &self_type::change_in_conductance),
                               grammar::field(C<self_type>{},"label",&self_type::label),
                               grammar::field(C<self_type>{},"par_on", &self_type::par_on),
                               grammar::field(C<self_type>{},"par_off", &self_type::par_off)
                               /*,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              grammar::field(C<self_type>{},"par_Eq", &self_type::par_Eq),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              grammar::field(C<self_type>{},"par_tc", &self_type::par_tc)*/);
  }
};
class Conformational_interaction;
class Conformational_interaction_label
    : public Label<Conformational_interaction> {
public:
  using Label::Label;
};

class Conformational_interaction {
  std::vector<Conformational_change_label> conf_changes_;
  Coupling_factor_Parameter_label factor_;
  std::vector<Coupling_coefficient_Parameter_label> par_inter_f_;

public:
  std::vector<Conformational_change_label> const &
  interacting_conformational_changes() const {
    return conf_changes_;
  }
  Coupling_factor_Parameter_label const &factor_label() const {
    return factor_;
  };
  std::vector<Coupling_coefficient_Parameter_label> const &
  coefficient_labels() const {
    return par_inter_f_;
  }

  typedef Conformational_interaction self_type;

  Conformational_interaction() = default;
  Conformational_interaction(
      std::vector<Conformational_change_label> _conf_changes,
      Coupling_factor_Parameter_label _factor,
      std::vector<Coupling_coefficient_Parameter_label> _par_inter_f)
      : conf_changes_{std::move(_conf_changes)}, factor_{std::move(_factor)},
        par_inter_f_{std::move(_par_inter_f)} {}
  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "interacting_conformational_changes",
                       &self_type::interacting_conformational_changes),
        grammar::field(C<self_type>{}, "factor_label",
                       &self_type::factor_label),
        grammar::field(C<self_type>{}, "coefficient_labels",
                       &self_type::coefficient_labels));
  }
  constexpr static auto className =
      my_static_string("Conformational_interaction");

  bool operator<(const Conformational_interaction &other) const {
    return factor_label() < other.factor_label();
  }
};

class State_Model {
  std::size_t numstates_;
  std::map<std::pair<std::size_t, std::size_t>,
           std::unique_ptr<Base_Function<double, State_Model>>>
      Q0s_;
  std::map<std::pair<std::size_t, std::size_t>,
           std::unique_ptr<Base_Function<double, State_Model>>>
      Qas_;
  std::map<std::size_t, std::unique_ptr<Base_Function<double, State_Model>>>
      gs_;

  std::vector<std::vector<std::size_t>> connections_to_0_;
  std::vector<std::vector<std::size_t>> connections_to_a_;
  std::vector<std::vector<std::size_t>> connections_to_x_;
  std::vector<std::vector<std::size_t>> connections_from_0_;
  std::vector<std::vector<std::size_t>> connections_from_x_;

  std::string errors_;
  bool isValid_;

  static auto getConnections(
      const std::map<std::pair<std::size_t, std::size_t>,
                     std::unique_ptr<Base_Function<double, State_Model>>> &Q0__,
      std::map<std::pair<std::size_t, std::size_t>,
               std::unique_ptr<Base_Function<double, State_Model>>> &Qa__,
      std::size_t k__) {
    std::vector<std::vector<std::size_t>> connections_to_0(k__);
    std::vector<std::vector<std::size_t>> connections_to_a(k__);
    std::vector<std::vector<std::size_t>> connections_to_x(k__);
    std::vector<std::vector<std::size_t>> connections_from_0(k__);
    std::vector<std::vector<std::size_t>> connections_from_x(k__);

    for (auto &e : Q0__) {
      connections_to_0[e.first.first].push_back(e.first.second);
      connections_to_x[e.first.first].push_back(e.first.second);
      connections_from_0[e.first.second].push_back(e.first.first);
      connections_from_x[e.first.second].push_back(e.first.first);
    }
    for (auto &e : Qa__) {
      connections_to_a[e.first.first].push_back(e.first.second);
      connections_to_x[e.first.first].push_back(e.first.second);
      connections_from_x[e.first.second].push_back(e.first.first);
    }
    return std::make_tuple(connections_to_0, connections_to_a, connections_to_x,
                           connections_from_0, connections_from_x);
  }

public:
  bool isValid() const { return isValid_; }
  std::string error() const { return errors_; }
  std::size_t k() const { return numstates_; }
  auto &transition_rates() const { return Q0s_; }
  std::map<std::pair<std::size_t, std::size_t>, std::string>
  transition_rates_text() const {
    return text_map<double, State_Model>(Q0s_);
  }

  auto &agonist_transitions_rates() const { return Qas_; }
  auto &connections_to_0() const { return connections_to_0_; }
  auto &connections_to_a() const { return connections_to_a_; }
  auto &connections_to_x(double x) const {
    if (x > 0)
      return connections_to_x_;
    else
      return connections_to_0_;
  }
  auto &connections_from_x(double x) const {
    if (x > 0)
      return connections_from_x_;
    else
      return connections_from_0_;
  }

  std::map<std::pair<std::size_t, std::size_t>, std::string>
  agonist_transitions_rates_text() const {
    return text_map<double, State_Model>(Qas_);
  }
  auto &conductances() const { return gs_; }
  std::map<std::size_t, std::string> conductances_text() const {
    return text_map(gs_);
  }

  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "number_of_states", &self_type::k),
        grammar::field(C<self_type>{}, "transition_rates",
                       &self_type::transition_rates_text),
        grammar::field(C<self_type>{}, "agonist_transitions_rates",
                       &self_type::agonist_transitions_rates_text),
        grammar::field(C<self_type>{}, "conductances",
                       &self_type::conductances_text));
  }
  typedef State_Model self_type;
  typedef Model_Parameter_label myParameter_label;
  constexpr static auto className = my_static_string("State_Model");

  static myOptional_t<State_Model>
  evaluate(std::size_t number_of_states,
           std::map<std::pair<std::size_t, std::size_t>, std::string>
               mytransition_rates,
           std::map<std::pair<std::size_t, std::size_t>, std::string>
               myagonist_transitions,
           std::map<std::size_t, std::string> myconductances) {
    typedef myOptional_t<State_Model> Op;
    auto transition_rates_r =
        compile_map<double, State_Model>(mytransition_rates);
    auto agonist_transitions_r =
        compile_map<double, State_Model>(myagonist_transitions);
    auto conductances_r = compile_map<double, State_Model>(myconductances);

    if (transition_rates_r.has_value() && agonist_transitions_r.has_value() &&
        conductances_r.has_value())
      return Op(State_Model(number_of_states,
                            std::move(transition_rates_r).value(),
                            std::move(agonist_transitions_r).value(),
                            std::move(conductances_r).value()));
    else
      return Op(false,
                "error building state model: " + transition_rates_r.error() +
                    agonist_transitions_r.error() + conductances_r.error());
  }

  State_Model() = default;

  State_Model(const State_Model &other)
      : numstates_{other.numstates_}, Q0s_{clone_map(other.Q0s_)},
        Qas_{clone_map(other.Qas_)}, gs_{clone_map(other.gs_)},
        connections_to_0_{other.connections_to_0_},
        connections_to_a_{other.connections_to_a_},
        connections_to_x_{other.connections_to_x_},
        connections_from_0_{other.connections_from_0_},
        connections_from_x_{other.connections_from_x_}, errors_{other.errors_},
        isValid_(other.isValid_) {}

  State_Model(State_Model &&other)
      : numstates_{other.numstates_}, Q0s_{std::move(other.Q0s_)},
        Qas_{std::move(other.Qas_)}, gs_{std::move(other.gs_)},
        connections_to_0_{std::move(other.connections_to_0_)},
        connections_to_a_{std::move(other.connections_to_a_)},
        connections_to_x_{std::move(other.connections_to_x_)},
        connections_from_0_{std::move(other.connections_from_0_)},
        connections_from_x_{std::move(other.connections_from_x_)},

        errors_{std::move(other.errors_)}, isValid_{other.isValid_} {}

  State_Model &operator=(const State_Model &other) {
    State_Model tmp(other);
    *this = std::move(tmp);
    return *this;
  }
  State_Model &operator=(State_Model &&other) {
    numstates_ = other.numstates_;
    Q0s_ = std::move(other.Q0s_);
    Qas_ = std::move(other.Qas_);
    gs_ = std::move(other.gs_);
    connections_to_0_ = std::move(other.connections_to_0_);
    connections_to_a_ = std::move(other.connections_to_a_);
    connections_to_x_ = std::move(other.connections_to_x_);
    connections_from_0_ = std::move(other.connections_from_0_);
    connections_from_x_ = std::move(other.connections_from_x_);

    errors_ = std::move(other.errors_);
    isValid_ = std::move(other.isValid_);
    return *this;
  }

  State_Model(std::size_t number_of_states,
              const std::map<std::pair<std::size_t, std::size_t>, std::string>
                  &mytransition_rates,
              const std::map<std::pair<std::size_t, std::size_t>, std::string>
                  &myagonist_transitions,
              const std::map<std::size_t, std::string> &myconductances) {
    auto transition_rates_r =
        compile_map<double, State_Model>(mytransition_rates);
    auto agonist_transitions_r =
        compile_map<double, State_Model>(myagonist_transitions);
    auto conductances_r = compile_map<double, State_Model>(myconductances);

    if (transition_rates_r.has_value() && agonist_transitions_r.has_value() &&
        conductances_r.has_value()) {
      numstates_ = number_of_states;
      Q0s_ = make_unique_map(std::move(transition_rates_r).value());
      Qas_ = make_unique_map(std::move(agonist_transitions_r).value());
      gs_ = make_unique_map(std::move(conductances_r).value());
      std::tie(connections_to_0_, connections_to_a_, connections_to_x_,
               connections_from_0_, connections_from_x_) =
          getConnections(Q0s_, Qas_, k());
      errors_.clear();
      isValid_ = true;

    } else {
      errors_ = transition_rates_r.error() + agonist_transitions_r.error() +
                conductances_r.error();
      isValid_ = false;
    }
  }

  State_Model(
      std::size_t number_of_states,
      std::map<std::pair<std::size_t, std::size_t>,
               Base_Function<double, State_Model> *> &&mytransition_rates,
      std::map<std::pair<std::size_t, std::size_t>,
               Base_Function<double, State_Model> *> &&myagonist_transitions,
      std::map<std::size_t, Base_Function<double, State_Model> *>
          &&myconductances)
      : numstates_{number_of_states}, Q0s_{make_unique_map(
                                          std::move(mytransition_rates))},
        Qas_{make_unique_map(std::move(myagonist_transitions))},
        gs_{make_unique_map(std::move(myconductances))}, errors_{},
        isValid_(true) {}

  template <class Parameters> auto Qs(const Parameters &p) const {

    M_Matrix<double> Q0(k(), k(), Matrix_TYPE::FULL, 0);
    M_Matrix<double> Qa(k(), k(), Matrix_TYPE::FULL, 0);

    for (auto &e : Q0s_) {
      Q0(e.first.first, e.first.second) = (*e.second)(p);
      Q0(e.first.first, e.first.first) -= Q0(e.first.first, e.first.second);
    }
    for (auto &e : Qas_) {
      Qa(e.first.first, e.first.second) = (*e.second)(p);
      Qa(e.first.first, e.first.first) -= Qa(e.first.first, e.first.second);
    }
    return std::pair(Q0, Qa);
  }

  template <class Parameters> auto g(const Parameters &p) const {
    M_Matrix<double> out(k(), 1, 0.0);
    for (auto &e : gs_)
      out[e.first] = (*e.second)(p);
    return out;
  }
};


template<template<class...> class Tr> class State_Model_new_;

typedef State_Model_new_<C> State_Model_new;


template<template<class...> class Tr>
class State_Model_new_ {
public:
    template<class T> using Tr_t=typename Tr<T>::type;
private:
    std::size_t numstates_;
    std::map<std::pair<std::size_t, std::size_t>,
             std::unique_ptr<Base_Function<double, State_Model_new>>>
        Q0s_;
    std::map<std::pair<std::size_t, std::size_t>,
             std::unique_ptr<Base_Function<double, State_Model_new>>>
        Qas_;
    std::map<std::size_t, std::unique_ptr<Base_Function<double, State_Model_new>>>
        gs_;

    std::vector<std::vector<std::size_t>> connections_to_0_;
    std::vector<std::vector<std::size_t>> connections_to_a_;
    std::vector<std::vector<std::size_t>> connections_to_x_;
    std::vector<std::vector<std::size_t>> connections_from_0_;
    std::vector<std::vector<std::size_t>> connections_from_x_;

    std::string errors_;
    bool isValid_;

    static auto getConnections(
        const std::map<std::pair<std::size_t, std::size_t>,
                       std::unique_ptr<Base_Function<double, State_Model_new>>> &Q0__,
        std::map<std::pair<std::size_t, std::size_t>,
                 std::unique_ptr<Base_Function<double, State_Model_new>>> &Qa__,
        std::size_t k__) {
        std::vector<std::vector<std::size_t>> connections_to_0(k__);
        std::vector<std::vector<std::size_t>> connections_to_a(k__);
        std::vector<std::vector<std::size_t>> connections_to_x(k__);
        std::vector<std::vector<std::size_t>> connections_from_0(k__);
        std::vector<std::vector<std::size_t>> connections_from_x(k__);

        for (auto &e : Q0__) {
            connections_to_0[e.first.first].push_back(e.first.second);
            connections_to_x[e.first.first].push_back(e.first.second);
            connections_from_0[e.first.second].push_back(e.first.first);
            connections_from_x[e.first.second].push_back(e.first.first);
        }
        for (auto &e : Qa__) {
            connections_to_a[e.first.first].push_back(e.first.second);
            connections_to_x[e.first.first].push_back(e.first.second);
            connections_from_x[e.first.second].push_back(e.first.first);
        }
        return std::make_tuple(connections_to_0, connections_to_a, connections_to_x,
                               connections_from_0, connections_from_x);
    }

public:
    bool isValid() const { return isValid_; }
    std::string error() const { return errors_; }
    std::size_t k() const { return numstates_; }
    auto &transition_rates() const { return Q0s_; }
    std::map<std::pair<std::size_t, std::size_t>, std::string>
    transition_rates_text() const {
        return text_map<double, State_Model_new>(Q0s_);
    }

    auto &agonist_transitions_rates() const { return Qas_; }
    auto &connections_to_0() const { return connections_to_0_; }
    auto &connections_to_a() const { return connections_to_a_; }
    auto &connections_to_x(double x) const {
        if (x > 0)
            return connections_to_x_;
        else
            return connections_to_0_;
    }
    auto &connections_from_x(double x) const {
        if (x > 0)
            return connections_from_x_;
        else
            return connections_from_0_;
    }

    std::map<std::pair<std::size_t, std::size_t>, std::string>
    agonist_transitions_rates_text() const {
        return text_map<double, State_Model_new>(Qas_);
    }
    auto &conductances() const { return gs_; }
    std::map<std::size_t, std::string> conductances_text() const {
        return text_map(gs_);
    }

    static auto get_constructor_fields() {
        return std::make_tuple(
            grammar::field(C<self_type>{}, "number_of_states", &self_type::k),
            grammar::field(C<self_type>{}, "transition_rates",
                           &self_type::transition_rates_text),
            grammar::field(C<self_type>{}, "agonist_transitions_rates",
                           &self_type::agonist_transitions_rates_text),
            grammar::field(C<self_type>{}, "conductances",
                           &self_type::conductances_text));
    }
    typedef State_Model_new self_type;
    typedef Model_Parameter_label myParameter_label;



    constexpr static auto className = my_static_string("State_Model_new");

    static myOptional_t<State_Model_new_<Tr> > evaluate(std::size_t number_of_states,
             std::map<std::pair<std::size_t, std::size_t>, std::string>
                 mytransition_rates,
             std::map<std::pair<std::size_t, std::size_t>, std::string>
                 myagonist_transitions,
             std::map<std::size_t, std::string> myconductances);

    State_Model_new_() = default;

    State_Model_new_(const State_Model_new_ &other)
        : numstates_{other.numstates_}, Q0s_{clone_map(other.Q0s_)},
          Qas_{clone_map(other.Qas_)}, gs_{clone_map(other.gs_)},
          connections_to_0_{other.connections_to_0_},
          connections_to_a_{other.connections_to_a_},
          connections_to_x_{other.connections_to_x_},
          connections_from_0_{other.connections_from_0_},
          connections_from_x_{other.connections_from_x_}, errors_{other.errors_},
          isValid_(other.isValid_) {}

    template<template<class...> class Tr_other>
    State_Model_new_(const State_Model_new_<Tr_other> &other)
        : numstates_{other.k()}, Q0s_{clone_map(other.transition_rates())},
          Qas_{clone_map(other.agonist_transitions_rates())}, gs_{clone_map(other.conductances())},
          connections_to_0_{other.connections_to_0()},
          connections_to_a_{other.connections_to_a()},
          connections_to_x_{other.connections_to_x(1.0)},
          connections_from_0_{other.connections_from_x(0.0)},
          connections_from_x_{other.connections_from_x(1.0)}, errors_{other.error()},
          isValid_(other.isValid()) {}


    State_Model_new_(State_Model_new_ &&other)
        : numstates_{other.numstates_}, Q0s_{std::move(other.Q0s_)},
          Qas_{std::move(other.Qas_)}, gs_{std::move(other.gs_)},
          connections_to_0_{std::move(other.connections_to_0_)},
          connections_to_a_{std::move(other.connections_to_a_)},
          connections_to_x_{std::move(other.connections_to_x_)},
          connections_from_0_{std::move(other.connections_from_0_)},
          connections_from_x_{std::move(other.connections_from_x_)},

          errors_{std::move(other.errors_)}, isValid_{other.isValid_} {}

    State_Model_new_ &operator=(const State_Model_new_ &other) {
        State_Model_new tmp(other);
        *this = std::move(tmp);
        return *this;
    }
    State_Model_new_ &operator=(State_Model_new_ &&other) {
        numstates_ = other.numstates_;
        Q0s_ = std::move(other.Q0s_);
        Qas_ = std::move(other.Qas_);
        gs_ = std::move(other.gs_);
        connections_to_0_ = std::move(other.connections_to_0_);
        connections_to_a_ = std::move(other.connections_to_a_);
        connections_to_x_ = std::move(other.connections_to_x_);
        connections_from_0_ = std::move(other.connections_from_0_);
        connections_from_x_ = std::move(other.connections_from_x_);

        errors_ = std::move(other.errors_);
        isValid_ = std::move(other.isValid_);
        return *this;
    }

    State_Model_new_(std::size_t number_of_states,
                const std::map<std::pair<std::size_t, std::size_t>, std::string>
                    &mytransition_rates,
                const std::map<std::pair<std::size_t, std::size_t>, std::string>
                    &myagonist_transitions,
                const std::map<std::size_t, std::string> &myconductances) {
        auto transition_rates_r =
            compile_map<double, State_Model_new>(mytransition_rates);
        auto agonist_transitions_r =
            compile_map<double, State_Model_new>(myagonist_transitions);
        auto conductances_r = compile_map<double, State_Model_new>(myconductances);

        if (transition_rates_r.has_value() && agonist_transitions_r.has_value() &&
            conductances_r.has_value()) {
            numstates_ = number_of_states;
            Q0s_ = make_unique_map(std::move(transition_rates_r).value());
            Qas_ = make_unique_map(std::move(agonist_transitions_r).value());
            gs_ = make_unique_map(std::move(conductances_r).value());
            std::tie(connections_to_0_, connections_to_a_, connections_to_x_,
                     connections_from_0_, connections_from_x_) =
                getConnections(Q0s_, Qas_, k());
            errors_.clear();
            isValid_ = true;

        } else {
            errors_ = transition_rates_r.error() + agonist_transitions_r.error() +
                conductances_r.error();
            isValid_ = false;
        }
    }

    State_Model_new_(
        std::size_t number_of_states,
        std::map<std::pair<std::size_t, std::size_t>,
                 Base_Function<double, State_Model_new> *> &&mytransition_rates,
        std::map<std::pair<std::size_t, std::size_t>,
                 Base_Function<double, State_Model_new> *> &&myagonist_transitions,
        std::map<std::size_t, Base_Function<double, State_Model_new> *>
            &&myconductances)
        : numstates_{number_of_states}, Q0s_{make_unique_map(
                                            std::move(mytransition_rates))},
          Qas_{make_unique_map(std::move(myagonist_transitions))},
          gs_{make_unique_map(std::move(myconductances))}, errors_{},
          isValid_(true) {}

    template <class Parameters> auto Qs(const Parameters &p) const {

        Tr_t<M_Matrix<double>> Q0(k(), k(), Matrix_TYPE::FULL, 0.0);
        Tr_t<M_Matrix<double>> Qa(k(), k(), Matrix_TYPE::FULL, 0.0);

        for (auto &e : Q0s_) {
            Q0.set(e.first.first, e.first.second,(*e.second)(p));
            Q0.add(e.first.first, e.first.first,- Q0(e.first.first, e.first.second));
        }
        for (auto &e : Qas_) {
            Qa.set(e.first.first, e.first.second,(*e.second)(p));
            Qa.add(e.first.first, e.first.first,-Qa(e.first.first, e.first.second));
        }
        return std::pair(Q0, Qa);
    }

    template <class Parameters> auto g(const Parameters &p) const {
        Tr_t<M_Matrix<double>> out(k(), 1, Matrix_TYPE::FULL,0.0);
        for (auto &e : gs_)
            out.set(e.first,(*e.second)(p));
        return out;
    }
};

template<template<class...> class Tr>
myOptional_t<State_Model_new_<Tr>> State_Model_new_<Tr>::evaluate(std::size_t number_of_states, std::map<std::pair<std::size_t, std::size_t>, std::string> mytransition_rates, std::map<std::pair<std::size_t, std::size_t>, std::string> myagonist_transitions, std::map<std::size_t, std::string> myconductances) {
    typedef myOptional_t<State_Model_new_<Tr>> Op;
    auto transition_rates_r =
        compile_map<double, State_Model_new>(mytransition_rates);
    auto agonist_transitions_r =
        compile_map<double, State_Model_new>(myagonist_transitions);
    auto conductances_r = compile_map<double, State_Model_new>(myconductances);

    if (transition_rates_r.has_value() && agonist_transitions_r.has_value() &&
        conductances_r.has_value())
        return Op(State_Model_new_(number_of_states,
                                  std::move(transition_rates_r).value(),
                                  std::move(agonist_transitions_r).value(),
                                  std::move(conductances_r).value()));
    else
        return Op(false,
                  "error building state model: " + transition_rates_r.error() +
                      agonist_transitions_r.error() + conductances_r.error());
}



/**
 * @brief The Allosteric_Model class build a kinetic rate model starting with a
 * set of conformational changes and their interactions.
 *
 *
 */
class Allosteric_Model {
public:
  typedef Allosteric_Model self_type;
  typedef Model_Parameter_label myParameter_label;
  struct transitions {
    bool on;
    bool agonist;
    std::string conformation;
    std::map<std::vector<std::pair<std::string, std::string>>, std::size_t>
        coupling;
  };
  struct model_definition {
    std::vector<std::string> conformational_changes;
    std::map<std::pair<std::string, bool>, std::string>
        conformational_changes_names;
    std::set<std::size_t> agonist_changes;
    std::set<std::size_t> conductance_changes;
    std::map<std::size_t, std::string> conductance_names;
    std::multimap<std::size_t, std::pair<std::set<std::size_t>,
                                         std::pair<std::string, std::string>>>
        conformational_inter_unit_cell;
  };

  struct new_model_definition {
    std::size_t number_of_units;
    std::map<Conformational_change_label, Conformational_change>
        conformational_changes;
    std::vector<Conformational_change_label> unit_of_conformational_changes;
    std::set<Conformational_interaction> conformational_interactions;
    std::map<std::size_t, Conductance_Parameter_label> conductance_names;
  };

  model_definition new_to_old(const new_model_definition &m) {
    model_definition out;
    for (auto &e : m.conductance_names)
      out.conductance_names[e.first] = e.second;
    auto k = m.unit_of_conformational_changes.size();
    std::size_t n = m.number_of_units * k;
    out.conformational_changes.resize(n);
    for (std::size_t i = 0; i < m.number_of_units; ++i)
      for (std::size_t j = 0; j < m.unit_of_conformational_changes.size();
           ++j) {
        auto ii = i * k + j;
        auto &cf = m.conformational_changes.at(
            m.unit_of_conformational_changes[j].name());
        out.conformational_changes[ii] = cf.label();
        if (cf.change_in_agonist() != 0)
          out.agonist_changes.insert(ii);
        if (cf.change_in_conductance() != 0)
          out.conductance_changes.insert(ii);
      }
    for (auto &e : m.conformational_changes) {
      out.conformational_changes_names[std::pair(e.second.label(), true)] =
          e.second.par_on();
      out.conformational_changes_names[std::pair(e.second.label(), false)] =
          e.second.par_off();
    }
    for (auto &e : m.conformational_interactions) {
      std::vector<std::size_t> cc;
      std::size_t current_i = 0;
      for (std::size_t j = 0; j < e.interacting_conformational_changes().size();
           ++j) {
        auto j_n = m.conformational_changes.at(
            e.interacting_conformational_changes()[j]);
        std::size_t i = current_i;
        while (j_n.label().name() != out.conformational_changes[i])
          ++i;
        cc.push_back(i);
        ++current_i;
      }

      for (std::size_t i = 0; i < cc.size(); ++i) {
        auto x = cc[i];
        auto ix = x % k;
        auto nx = x / k;
        auto shift = (n - nx) * k;
        std::set<std::size_t> s;
        for (std::size_t j = 0; j < cc.size(); ++j)
          if (j != i)
            s.insert(cc[j]);
        s = rotate(s, n * k, shift);
        out.conformational_inter_unit_cell.emplace(
            ix, std::pair(
                    s, std::pair(e.factor_label(), e.coefficient_labels()[i])));
      }
    }
    return out;
  }

  constexpr static auto className = my_static_string("Allosteric_Model");

protected:
  new_model_definition new_d_;
  std::vector<std::vector<std::string>> conformer_;
  std::vector<std::map<std::size_t, transitions>> transitions_;
  std::set<std::string> paramNames_;
  std::vector<std::string> conductances_;
  model_definition d_;

  std::vector<std::vector<std::size_t>> connections_to_0_;
  std::vector<std::vector<std::size_t>> connections_to_a_;
  std::vector<std::vector<std::size_t>> connections_to_x_;
  std::vector<std::vector<std::size_t>> connections_from_0_;
  std::vector<std::vector<std::size_t>> connections_from_x_;

  /**
   * @brief getParameterNamesFrom
   * @param conformational_changes_names
   * @param conductance_names
   * @param conformational_inter
   */
  static auto getParameterNamesFrom(
      const std::map<std::pair<std::string, bool>, std::string>
          conformational_changes_names,
      const std::map<std::size_t, std::string> &conductance_names,
      const std::multimap<
          std::size_t,
          std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>>
          &conformational_inter) {
    std::set<std::string> out;
    for (auto &e : conformational_changes_names)
      out.insert(e.second);
    for (auto &e : conductance_names)
      out.insert(e.second);
    for (auto &e : conformational_inter) {
      out.insert(e.second.second.first);
      out.insert(e.second.second.second);
    }
    return out;
  }

  static auto
  getConformers(const std::vector<std::string> &conformational_changes,
                std::size_t p) {
    std::vector<std::vector<std::string>> conformer;
    std::map<std::vector<std::string>, std::size_t> state_to_conformer;
    for (std::size_t i = 0; i < (1u << conformational_changes.size()); ++i) {
      auto c = index_to_conformational_change(conformational_changes, i);
      if (state_to_conformer.find(c) == state_to_conformer.end()) {
        state_to_conformer[c] = conformer.size();
        std::size_t n = 1;
        auto permute = rotate(c, p * n);
        while (c != permute) {
          state_to_conformer[permute] = conformer.size();
          ++n;
          permute = rotate(c, p * n);
        }
        conformer.push_back(c);
      }
    }
    return std::make_tuple(conformer, state_to_conformer);
  }

  static auto getTransitions(
      const std::vector<std::string> &conformational_changes,
      const std::vector<std::vector<std::string>> &conformer,
      const std::map<std::vector<std::string>, std::size_t> &state_to_conformer,
      const std::set<std::size_t> &agonist_changes,
      const std::map<std::pair<std::string, bool>, std::string>
          conformational_changes_names,
      std::multimap<std::size_t, std::pair<std::set<std::size_t>,
                                           std::pair<std::string, std::string>>>
          conformational_interactions) {
    std::vector<std::map<std::size_t, transitions>> transitions_out;
    for (std::size_t i = 0; i < conformer.size(); ++i) {
      std::map<std::size_t, transitions> myTransition;
      auto c = conformer[i];
      for (std::size_t k = 0; k < c.size(); ++k) {
        auto change = c;
        if (c[k] == "")
          change[k] = conformational_changes[k];
        else
          change[k] = "";
        auto j = state_to_conformer.at(change);
        if (myTransition.find(j) == myTransition.end()) {

          myTransition[j].agonist =
              agonist_changes.find(k) != agonist_changes.end();
          bool on = change[k] == conformational_changes[k];
          myTransition[j].on = on;
          auto name = conformational_changes_names.at(
              std::pair{conformational_changes[k], on});
          myTransition[j].conformation = name;
        }
        auto m = conformational_interactions.equal_range(k);
        std::vector<std::pair<std::string, std::string>> coupling;
        for (auto it = m.first; it != m.second; ++it) {
          std::set<std::size_t> s = it->second.first;
          bool all = true;
          for (auto e : s) {
            if (c[e].empty()) {
              all = false;
              break;
            }
          }
          if (all)
            coupling.push_back(it->second.second);
        }
        ++myTransition[j].coupling[coupling];
      }
      transitions_out.push_back(std::move(myTransition));
    }
    return transitions_out;
  }

  static auto getConections(
      const std::vector<std::map<std::size_t, transitions>> &mytransitions) {
    std::vector<std::vector<std::size_t>> connections_to_0(
        mytransitions.size());
    std::vector<std::vector<std::size_t>> connections_to_a(
        mytransitions.size());
    std::vector<std::vector<std::size_t>> connections_to_x(
        mytransitions.size());
    std::vector<std::vector<std::size_t>> connections_from_0(
        mytransitions.size());
    std::vector<std::vector<std::size_t>> connections_from_x(
        mytransitions.size());
    for (std::size_t i = 0; i < mytransitions.size(); ++i)
      for (auto it = mytransitions[i].begin(); it != mytransitions[i].end();
           ++it) {
        std::size_t j = it->first;
        connections_from_x[j].push_back(i);
        connections_to_x[i].push_back(j);
        if ((it->second.agonist) && (it->second.on)) {
          connections_to_a[i].push_back(j);
        } else {
          connections_to_0[i].push_back(j);
          connections_from_0[j].push_back(i);
        }
      }
    return std::make_tuple(connections_to_0, connections_to_a, connections_to_x,
                           connections_from_0, connections_from_x);
  }

  static

      std::string
      g_name_of_conformer(
          const std::vector<std::string> &c,
          const std::set<std::size_t> &conductance_changes,
          const std::map<std::size_t, std::string> &conductance_names) {
    std::size_t i = 0;
    for (auto e : conductance_changes)
      if (!c.at(e).empty())
        ++i;
    return conductance_names.at(i);
  }

  static std::vector<std::string>
  get_g_names(const std::vector<std::vector<std::string>> &conformers,
              const std::set<std::size_t> &conductance_changes,
              const std::map<std::size_t, std::string> &conductance_names) {
    std::vector<std::string> out(conformers.size());
    for (std::size_t i = 0; i < conformers.size(); ++i)
      out[i] = g_name_of_conformer(conformers[i], conductance_changes,
                                   conductance_names);
    return out;
  }

public:
  static auto get_constructor_fields_old() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "conformational_changes",
                       &self_type::get_conformational_changes),
        grammar::field(C<self_type>{}, "conformational_changes_names",
                       &self_type::get_conformational_changes_names),
        grammar::field(C<self_type>{}, "agonist_changes",
                       &self_type::get_agonist_changes),
        grammar::field(C<self_type>{}, "conductance_changes",
                       &self_type::get_conductance_changes),
        grammar::field(C<self_type>{}, "conductance_names",
                       &self_type::get_conductance_names),
        grammar::field(C<self_type>{}, "conformational_inter_unit_cell",
                       &self_type::get_conformational_inter_unit_cell));
  }

  const std::vector<std::string> &get_conformational_changes() const {
    return d_.conformational_changes;
  }
  const std::map<std::pair<std::string, bool>, std::string> &
  get_conformational_changes_names() const {
    return d_.conformational_changes_names;
  }
  const std::set<std::size_t> &get_agonist_changes() const {
    return d_.agonist_changes;
  }
  const std::set<std::size_t> &get_conductance_changes() const {
    return d_.conductance_changes;
  }
  const std::map<std::size_t, std::string> &get_conductance_names() const {
    return d_.conductance_names;
  }
  const std::multimap<
      std::size_t,
      std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>> &
  get_conformational_inter_unit_cell() const {
    return d_.conformational_inter_unit_cell;
  }

  std::size_t number_of_units() const { return new_d_.number_of_units; }
  std::map<Conformational_change_label, Conformational_change> const &
  conformational_changes() const {
    return new_d_.conformational_changes;
  }
  std::vector<Conformational_change_label> const &
  unit_of_conformational_changes() const {
    return new_d_.unit_of_conformational_changes;
  }
  std::set<Conformational_interaction> const &
  conformational_interactions() const {
    return new_d_.conformational_interactions;
  }
  std::map<std::size_t, Conductance_Parameter_label> const &
  conductance_names() const {
    return new_d_.conductance_names;
  }

  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "number_of_units",
                       &self_type::number_of_units),
        grammar::field(C<self_type>{}, "conformational_changes",
                       &self_type::conformational_changes),
        grammar::field(C<self_type>{}, "unit_of_conformational_changes",
                       &self_type::unit_of_conformational_changes),
        grammar::field(C<self_type>{}, "conformational_interactions",
                       &self_type::conformational_interactions),
        grammar::field(C<self_type>{}, "conductance_names",
                       &self_type::conductance_names));
  }

  Allosteric_Model() = default;

  Allosteric_Model(
      std::size_t number_of_units,
      std::map<Conformational_change_label, Conformational_change>
          conformational_changes,
      std::vector<Conformational_change_label> unit_of_conformational_changes,
      std::set<Conformational_interaction> conformational_interactions,
      std::map<std::size_t, Conductance_Parameter_label> conductance_names)
      : new_d_{number_of_units, conformational_changes,
               unit_of_conformational_changes, conformational_interactions,
               conductance_names},
        d_{new_to_old(new_d_)} {
    init();
  }

  Allosteric_Model(
      const std::vector<std::string> &conformational_changes,
      const std::map<std::pair<std::string, bool>, std::string>
          conformational_changes_names,
      const std::set<std::size_t> &agonist_changes,
      const std::set<std::size_t> &conductance_changes,
      const std::map<std::size_t, std::string> &conductance_names,
      const std::multimap<
          std::size_t,
          std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>>
          &conformational_inter_unit_cell)
      : d_{conformational_changes, conformational_changes_names,
           agonist_changes,        conductance_changes,
           conductance_names,      conformational_inter_unit_cell} {
    init();
  }

  template <class Parameters> auto Qs(const Parameters &p) const {

    M_Matrix<double> Q0(conformer_.size(), conformer_.size());
    M_Matrix<double> Qa(conformer_.size(), conformer_.size());
    for (std::size_t i = 0; i < Q0.nrows(); ++i) {
      Q0(i, i) = 0;
      Qa(i, i) = 0;
      for (auto it = transitions_[i].begin(); it != transitions_[i].end();
           ++it) {
        std::size_t j = it->first;
        if ((it->second.agonist) && (it->second.on)) {
          Qa(i, j) = rate(it->second, p);
          Qa(i, i) -= Qa(i, j);
        } else {
          Q0(i, j) = rate(it->second, p);
          Q0(i, i) -= Q0(i, j);
        }
      }
    }
    auto r = Q0 * ones<double>(Q0.ncols(), 1);
    auto q = Qa * ones<double>(Q0.ncols(), 1);

    return std::pair(Q0, Qa);
  }

  template <class Parameters> auto g(const Parameters &p) const {
    M_Matrix<double> out(conformer_.size(), 1);
    for (std::size_t i = 0; i < conformer_.size(); ++i)
      out(i, 0) = p.at(conductances_.at(i));
    return out;
  }

  std::size_t k() const { return conformer_.size(); }

  auto &connections_to_0() const { return connections_to_0_; }
  auto &connections_to_a() const { return connections_to_a_; }
  auto &connections_to_x(double x) const {
    if (x > 0)
      return connections_to_x_;
    else
      return connections_to_0_;
  }
  auto &connections_from_x(double x) const {
    if (x > 0)
      return connections_from_x_;
    else
      return connections_from_0_;
  }
  auto getParameterNames() const { return paramNames_; }

private:
  void init() {
    paramNames_ = getParameterNamesFrom(d_.conformational_changes_names,
                                        d_.conductance_names,
                                        d_.conformational_inter_unit_cell);
    std::cout << paramNames_;
    auto p = periodicity(d_.conformational_changes);

    auto conformational_interactions = fill_conformational_interactions(
        d_.conformational_inter_unit_cell, d_.conformational_changes.size(), p);

    std::map<std::vector<std::string>, std::size_t> state_to_conformer;

    std::tie(conformer_, state_to_conformer) =
        getConformers(d_.conformational_changes, p);

    transitions_ = getTransitions(d_.conformational_changes, conformer_,
                                  state_to_conformer, d_.agonist_changes,
                                  d_.conformational_changes_names,
                                  conformational_interactions);
    conductances_ =
        get_g_names(conformer_, d_.conductance_changes, d_.conductance_names);

    std::tie(connections_to_0_, connections_to_a_, connections_to_x_,
             connections_from_0_, connections_from_x_) =
        getConections(transitions_);
  }

  static std::vector<std::string>
  index_to_conformational_change(const std::vector<std::string> &cc,
                                 std::size_t index) {
    std::vector<std::string> out(cc.size(), "");
    for (std::size_t i = 0; i < cc.size(); ++i)
      if (((index >> i) & 1) == 1)
        out[i] = cc[i];
    return out;
  }

  static std::vector<std::string> rotate(const std::vector<std::string> &x,
                                         std::size_t n) {
    std::vector<std::string> out(x.size());
    for (std::size_t i = 0; i < out.size(); ++i)
      out[(i + n) % out.size()] = x[i];
    return out;
  }

  static std::size_t periodicity(const std::vector<std::string> &x) {
    std::size_t p = 1;
    auto t = rotate(x, p);
    while (t != x && p < x.size()) {
      ++p;
      t = rotate(x, p);
    }
    return p;
  }

  static std::set<std::size_t> rotate(const std::set<std::size_t> &c,
                                      std::size_t n, std::size_t i) {
    std::set<std::size_t> out;
    for (auto e : c)
      out.insert((e + i) % n);
    return out;
  }

  static std::multimap<
      std::size_t,
      std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>>
  fill_conformational_interactions(
      const std::multimap<
          std::size_t,
          std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>>
          &conformational_interactions,
      std::size_t n, std::size_t p) {
    std::multimap<std::size_t, std::pair<std::set<std::size_t>,
                                         std::pair<std::string, std::string>>>
        out;
    for (std::size_t i = 0; i < n / p; ++i) {
      for (auto &e : conformational_interactions) {
        out.insert({(e.first + i * p) % n,
                    {rotate(e.second.first, n, i * p), e.second.second}});
      }
    }
    return out;
  }

  template <class P> static double rate(const transitions &tr, const P &p) {
    double out = 0;
    if (tr.on) {
      for (auto &e : tr.coupling) {
        double b = e.second;
        for (auto &e2 : e.first)
          b *= std::pow(p.at(e2.first), p.at(e2.second));
        out += b;
      }
    } else {
      for (auto &e : tr.coupling) {
        double b = e.second;
        for (auto &e2 : e.first)
          b *= std::pow(p.at(e2.first), p.at(e2.second) - 1.0);
        out += b;
      }
    }
    return out * p.at(tr.conformation);
  }
};


template<template<class...> class Tr> class Allosteric_Model_new_;

typedef Allosteric_Model_new_<C> Allosteric_Model_new;


/**
 * @brief The Allosteric_Model class build a kinetic rate model starting with a
 * set of conformational changes and their interactions.
 *
 *
 */
struct new_model_definition {
    std::size_t number_of_units;
    std::map<Conformational_change_label, Conformational_change>
        conformational_changes;
    std::vector<Conformational_change_label> unit_of_conformational_changes;
    std::set<Conformational_interaction> conformational_interactions;
    std::map<std::size_t, Conductance_Parameter_label> conductance_names;
};
struct transitions {
    bool on;
    bool agonist;
    std::string conformation;
    std::map<std::vector<std::pair<std::string, std::string>>, std::size_t>
        coupling;
};
struct model_definition {
    std::vector<std::string> conformational_changes;
    std::map<std::pair<std::string, bool>, std::string>
        conformational_changes_names;
    std::set<std::size_t> agonist_changes;
    std::set<std::size_t> conductance_changes;
    std::map<std::size_t, std::string> conductance_names;
    std::multimap<std::size_t, std::pair<std::set<std::size_t>,
                                         std::pair<std::string, std::string>>>
        conformational_inter_unit_cell;
};

template<template<class...> class Tr>
class Allosteric_Model_new_ {
public:
    template<class T> using Tr_t=typename Tr<T>::type;


    typedef Allosteric_Model_new_ self_type;
    typedef Model_Parameter_label myParameter_label;


    model_definition new_to_old(const new_model_definition &m) {
        model_definition out;
        for (auto &e : m.conductance_names)
            out.conductance_names[e.first] = e.second;
        auto k = m.unit_of_conformational_changes.size();
        std::size_t n = m.number_of_units * k;
        out.conformational_changes.resize(n);
        for (std::size_t i = 0; i < m.number_of_units; ++i)
            for (std::size_t j = 0; j < m.unit_of_conformational_changes.size();
                 ++j) {
                auto ii = i * k + j;
                auto &cf = m.conformational_changes.at(
                    m.unit_of_conformational_changes[j].name());
                out.conformational_changes[ii] = cf.label();
                if (cf.change_in_agonist() != 0)
                    out.agonist_changes.insert(ii);
                if (cf.change_in_conductance() != 0)
                    out.conductance_changes.insert(ii);
            }
        for (auto &e : m.conformational_changes) {
            out.conformational_changes_names[std::pair(e.second.label(), true)] =
                e.second.par_on();
            out.conformational_changes_names[std::pair(e.second.label(), false)] =
                e.second.par_off();
        }
        for (auto &e : m.conformational_interactions) {
            std::vector<std::size_t> cc;
            std::size_t current_i = 0;
            for (std::size_t j = 0; j < e.interacting_conformational_changes().size();
                 ++j) {
                auto j_n = m.conformational_changes.at(
                    e.interacting_conformational_changes()[j]);
                std::size_t i = current_i;
                while (j_n.label().name() != out.conformational_changes[i])
                    ++i;
                cc.push_back(i);
                ++current_i;
            }

            for (std::size_t i = 0; i < cc.size(); ++i) {
                auto x = cc[i];
                auto ix = x % k;
                auto nx = x / k;
                auto shift = (n - nx) * k;
                std::set<std::size_t> s;
                for (std::size_t j = 0; j < cc.size(); ++j)
                    if (j != i)
                        s.insert(cc[j]);
                s = rotate(s, n * k, shift);
                out.conformational_inter_unit_cell.emplace(
                    ix, std::pair(
                            s, std::pair(e.factor_label(), e.coefficient_labels()[i])));
            }
        }
        return out;
    }

    constexpr static auto className = my_static_string("Allosteric_Model_new");
    auto& get_conductances()const{ return conductances_;}
    auto& get_paramNames()const {return paramNames_;}
    auto & get_transitions()const {return transitions_;}
    auto &get_conformers()const {return conformer_;}
    auto & get_new_model_definition()const {return new_d_;}

protected:
    new_model_definition new_d_;
    std::vector<std::vector<std::string>> conformer_;
    std::vector<std::map<std::size_t, transitions>> transitions_;
    std::set<std::string> paramNames_;
    std::vector<std::string> conductances_;
    model_definition d_;

    std::vector<std::vector<std::size_t>> connections_to_0_;
    std::vector<std::vector<std::size_t>> connections_to_a_;
    std::vector<std::vector<std::size_t>> connections_to_x_;
    std::vector<std::vector<std::size_t>> connections_from_0_;
    std::vector<std::vector<std::size_t>> connections_from_x_;

    /**
   * @brief getParameterNamesFrom
   * @param conformational_changes_names
   * @param conductance_names
   * @param conformational_inter
   */
    static auto getParameterNamesFrom(
        const std::map<std::pair<std::string, bool>, std::string>
            conformational_changes_names,
        const std::map<std::size_t, std::string> &conductance_names,
        const std::multimap<
            std::size_t,
            std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>>
            &conformational_inter) {
        std::set<std::string> out;
        for (auto &e : conformational_changes_names)
            out.insert(e.second);
        for (auto &e : conductance_names)
            out.insert(e.second);
        for (auto &e : conformational_inter) {
            out.insert(e.second.second.first);
            out.insert(e.second.second.second);
        }
        return out;
    }

    static auto
    getConformers(const std::vector<std::string> &conformational_changes,
                  std::size_t p) {
        std::vector<std::vector<std::string>> conformer;
        std::map<std::vector<std::string>, std::size_t> state_to_conformer;
        for (std::size_t i = 0; i < (1u << conformational_changes.size()); ++i) {
            auto c = index_to_conformational_change(conformational_changes, i);
            if (state_to_conformer.find(c) == state_to_conformer.end()) {
                state_to_conformer[c] = conformer.size();
                std::size_t n = 1;
                auto permute = rotate(c, p * n);
                while (c != permute) {
                    state_to_conformer[permute] = conformer.size();
                    ++n;
                    permute = rotate(c, p * n);
                }
                conformer.push_back(c);
            }
        }
        return std::make_tuple(conformer, state_to_conformer);
    }

    static auto getTransitions(
        const std::vector<std::string> &conformational_changes,
        const std::vector<std::vector<std::string>> &conformer,
        const std::map<std::vector<std::string>, std::size_t> &state_to_conformer,
        const std::set<std::size_t> &agonist_changes,
        const std::map<std::pair<std::string, bool>, std::string>
            conformational_changes_names,
        std::multimap<std::size_t, std::pair<std::set<std::size_t>,
                                             std::pair<std::string, std::string>>>
            conformational_interactions) {
        std::vector<std::map<std::size_t, transitions>> transitions_out;
        for (std::size_t i = 0; i < conformer.size(); ++i) {
            std::map<std::size_t, transitions> myTransition;
            auto c = conformer[i];
            for (std::size_t k = 0; k < c.size(); ++k) {
                auto change = c;
                if (c[k] == "")
                    change[k] = conformational_changes[k];
                else
                    change[k] = "";
                auto j = state_to_conformer.at(change);
                if (myTransition.find(j) == myTransition.end()) {

                    myTransition[j].agonist =
                        agonist_changes.find(k) != agonist_changes.end();
                    bool on = change[k] == conformational_changes[k];
                    myTransition[j].on = on;
                    auto name = conformational_changes_names.at(
                        std::pair{conformational_changes[k], on});
                    myTransition[j].conformation = name;
                }
                auto m = conformational_interactions.equal_range(k);
                std::vector<std::pair<std::string, std::string>> coupling;
                for (auto it = m.first; it != m.second; ++it) {
                    std::set<std::size_t> s = it->second.first;
                    bool all = true;
                    for (auto e : s) {
                        if (c[e].empty()) {
                            all = false;
                            break;
                        }
                    }
                    if (all)
                        coupling.push_back(it->second.second);
                }
                ++myTransition[j].coupling[coupling];
            }
            transitions_out.push_back(std::move(myTransition));
        }
        return transitions_out;
    }

    static auto getConections(
        const std::vector<std::map<std::size_t, transitions>> &mytransitions) {
        std::vector<std::vector<std::size_t>> connections_to_0(
            mytransitions.size());
        std::vector<std::vector<std::size_t>> connections_to_a(
            mytransitions.size());
        std::vector<std::vector<std::size_t>> connections_to_x(
            mytransitions.size());
        std::vector<std::vector<std::size_t>> connections_from_0(
            mytransitions.size());
        std::vector<std::vector<std::size_t>> connections_from_x(
            mytransitions.size());
        for (std::size_t i = 0; i < mytransitions.size(); ++i)
            for (auto it = mytransitions[i].begin(); it != mytransitions[i].end();
                 ++it) {
                std::size_t j = it->first;
                connections_from_x[j].push_back(i);
                connections_to_x[i].push_back(j);
                if ((it->second.agonist) && (it->second.on)) {
                    connections_to_a[i].push_back(j);
                } else {
                    connections_to_0[i].push_back(j);
                    connections_from_0[j].push_back(i);
                }
            }
        return std::make_tuple(connections_to_0, connections_to_a, connections_to_x,
                               connections_from_0, connections_from_x);
    }

    static

        std::string
        g_name_of_conformer(
            const std::vector<std::string> &c,
            const std::set<std::size_t> &conductance_changes,
            const std::map<std::size_t, std::string> &conductance_names) {
        std::size_t i = 0;
        for (auto e : conductance_changes)
            if (!c.at(e).empty())
                ++i;
        return conductance_names.at(i);
    }

    static std::vector<std::string>
    get_g_names(const std::vector<std::vector<std::string>> &conformers,
                const std::set<std::size_t> &conductance_changes,
                const std::map<std::size_t, std::string> &conductance_names) {
        std::vector<std::string> out(conformers.size());
        for (std::size_t i = 0; i < conformers.size(); ++i)
            out[i] = g_name_of_conformer(conformers[i], conductance_changes,
                                         conductance_names);
        return out;
    }

public:
    static auto get_constructor_fields_old() {
        return std::make_tuple(
            grammar::field(C<self_type>{}, "conformational_changes",
                           &self_type::get_conformational_changes),
            grammar::field(C<self_type>{}, "conformational_changes_names",
                           &self_type::get_conformational_changes_names),
            grammar::field(C<self_type>{}, "agonist_changes",
                           &self_type::get_agonist_changes),
            grammar::field(C<self_type>{}, "conductance_changes",
                           &self_type::get_conductance_changes),
            grammar::field(C<self_type>{}, "conductance_names",
                           &self_type::get_conductance_names),
            grammar::field(C<self_type>{}, "conformational_inter_unit_cell",
                           &self_type::get_conformational_inter_unit_cell));
    }


    auto& get_model_definition()const { return d_;}
    const std::vector<std::string> &get_conformational_changes() const {
        return d_.conformational_changes;
    }
    const std::map<std::pair<std::string, bool>, std::string> &
    get_conformational_changes_names() const {
        return d_.conformational_changes_names;
    }
    const std::set<std::size_t> &get_agonist_changes() const {
        return d_.agonist_changes;
    }
    const std::set<std::size_t> &get_conductance_changes() const {
        return d_.conductance_changes;
    }
    const std::map<std::size_t, std::string> &get_conductance_names() const {
        return d_.conductance_names;
    }
    const std::multimap<
        std::size_t,
        std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>> &
    get_conformational_inter_unit_cell() const {
        return d_.conformational_inter_unit_cell;
    }

    std::size_t number_of_units() const { return new_d_.number_of_units; }
    std::map<Conformational_change_label, Conformational_change> const &
    conformational_changes() const {
        return new_d_.conformational_changes;
    }
    std::vector<Conformational_change_label> const &
    unit_of_conformational_changes() const {
        return new_d_.unit_of_conformational_changes;
    }
    std::set<Conformational_interaction> const &
    conformational_interactions() const {
        return new_d_.conformational_interactions;
    }
    std::map<std::size_t, Conductance_Parameter_label> const &
    conductance_names() const {
        return new_d_.conductance_names;
    }

    static auto get_constructor_fields() {
        return std::make_tuple(
            grammar::field(C<self_type>{}, "number_of_units",
                           &self_type::number_of_units),
            grammar::field(C<self_type>{}, "conformational_changes",
                           &self_type::conformational_changes),
            grammar::field(C<self_type>{}, "unit_of_conformational_changes",
                           &self_type::unit_of_conformational_changes),
            grammar::field(C<self_type>{}, "conformational_interactions",
                           &self_type::conformational_interactions),
            grammar::field(C<self_type>{}, "conductance_names",
                           &self_type::conductance_names));
    }

    Allosteric_Model_new_() = default;

    Allosteric_Model_new_(const Allosteric_Model_new& other);


    Allosteric_Model_new_(
        std::size_t number_of_units,
        std::map<Conformational_change_label, Conformational_change>
            conformational_changes,
        std::vector<Conformational_change_label> unit_of_conformational_changes,
        std::set<Conformational_interaction> conformational_interactions,
        std::map<std::size_t, Conductance_Parameter_label> conductance_names)
        : new_d_{number_of_units, conformational_changes,
                 unit_of_conformational_changes, conformational_interactions,
                 conductance_names},
          d_{new_to_old(new_d_)} {
        init();
    }

    Allosteric_Model_new_(
        const std::vector<std::string> &conformational_changes,
        const std::map<std::pair<std::string, bool>, std::string>
            conformational_changes_names,
        const std::set<std::size_t> &agonist_changes,
        const std::set<std::size_t> &conductance_changes,
        const std::map<std::size_t, std::string> &conductance_names,
        const std::multimap<
            std::size_t,
            std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>>
            &conformational_inter_unit_cell)
        : d_{conformational_changes, conformational_changes_names,
             agonist_changes,        conductance_changes,
             conductance_names,      conformational_inter_unit_cell} {
        init();
    }

    template <class Parameters> auto Qs(const Parameters &p) const {

        Tr_t<M_Matrix<double>> Q0(conformer_.size(), conformer_.size());
        Tr_t<M_Matrix<double>> Qa(conformer_.size(), conformer_.size());
        for (std::size_t i = 0; i < Q0.nrows(); ++i) {
            Q0.set(i, i,0.0);
            Qa.set(i, i,0.0);
            for (auto it = transitions_[i].begin(); it != transitions_[i].end();
                 ++it) {
                std::size_t j = it->first;
                if ((it->second.agonist) && (it->second.on)) {
                    Qa.set(i, j,rate(it->second, p));
                    Qa.add(i, i,  - Qa(i, j));
                } else {
                    Q0.set(i, j, rate(it->second, p));
                    Q0.add(i, i, - Q0(i, j));
                }
            }
        }
        auto r = Q0 * ones<double>(Q0.ncols(), 1);
        auto q = Qa * ones<double>(Q0.ncols(), 1);

        return std::pair(Q0, Qa);
    }

    template <class Parameters> auto g(const Parameters &p) const {
        Tr_t<M_Matrix<double>> out(conformer_.size(), 1);
        for (std::size_t i = 0; i < conformer_.size(); ++i)
            out.set(i, 0, p.at(conductances_.at(i)));
        return out;
    }

    std::size_t k() const { return conformer_.size(); }

    auto &connections_to_0() const { return connections_to_0_; }
    auto &connections_to_a() const { return connections_to_a_; }
    auto &connections_to_x(double x) const {
        if (x > 0)
            return connections_to_x_;
        else
            return connections_to_0_;
    }
    auto &connections_from_x(double x) const {
        if (x > 0)
            return connections_from_x_;
        else
            return connections_from_0_;
    }
    auto getParameterNames() const { return paramNames_; }

private:
    void init() {
        paramNames_ = getParameterNamesFrom(d_.conformational_changes_names,
                                            d_.conductance_names,
                                            d_.conformational_inter_unit_cell);
        std::cout << paramNames_;
        auto p = periodicity(d_.conformational_changes);

        auto conformational_interactions = fill_conformational_interactions(
            d_.conformational_inter_unit_cell, d_.conformational_changes.size(), p);

        std::map<std::vector<std::string>, std::size_t> state_to_conformer;

        std::tie(conformer_, state_to_conformer) =
            getConformers(d_.conformational_changes, p);

        transitions_ = getTransitions(d_.conformational_changes, conformer_,
                                      state_to_conformer, d_.agonist_changes,
                                      d_.conformational_changes_names,
                                      conformational_interactions);
        conductances_ =
            get_g_names(conformer_, d_.conductance_changes, d_.conductance_names);

        std::tie(connections_to_0_, connections_to_a_, connections_to_x_,
                 connections_from_0_, connections_from_x_) =
            getConections(transitions_);
    }

    static std::vector<std::string>
    index_to_conformational_change(const std::vector<std::string> &cc,
                                   std::size_t index) {
        std::vector<std::string> out(cc.size(), "");
        for (std::size_t i = 0; i < cc.size(); ++i)
            if (((index >> i) & 1) == 1)
                out[i] = cc[i];
        return out;
    }

    static std::vector<std::string> rotate(const std::vector<std::string> &x,
                                           std::size_t n) {
        std::vector<std::string> out(x.size());
        for (std::size_t i = 0; i < out.size(); ++i)
            out[(i + n) % out.size()] = x[i];
        return out;
    }

    static std::size_t periodicity(const std::vector<std::string> &x) {
        std::size_t p = 1;
        auto t = rotate(x, p);
        while (t != x && p < x.size()) {
            ++p;
            t = rotate(x, p);
        }
        return p;
    }

    static std::set<std::size_t> rotate(const std::set<std::size_t> &c,
                                        std::size_t n, std::size_t i) {
        std::set<std::size_t> out;
        for (auto e : c)
            out.insert((e + i) % n);
        return out;
    }

    static std::multimap<
        std::size_t,
        std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>>
    fill_conformational_interactions(
        const std::multimap<
            std::size_t,
            std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>>
            &conformational_interactions,
        std::size_t n, std::size_t p) {
        std::multimap<std::size_t, std::pair<std::set<std::size_t>,
                                             std::pair<std::string, std::string>>>
            out;
        for (std::size_t i = 0; i < n / p; ++i) {
            for (auto &e : conformational_interactions) {
                out.insert({(e.first + i * p) % n,
                            {rotate(e.second.first, n, i * p), e.second.second}});
            }
        }
        return out;
    }

    template <class P> static Tr_t<double> rate(const transitions &tr, const P &p) {
        Tr_t<double> out(0);
        using std::pow;
        if (tr.on) {
            for (auto &e : tr.coupling) {
                Tr_t<double> b = e.second;
                for (auto &e2 : e.first)
                    b *= pow(p.at(e2.first), p.at(e2.second));
                out += b;
            }
        } else {
            for (auto &e : tr.coupling) {
                Tr_t<double> b = e.second;
                for (auto &e2 : e.first)
                    b *= pow(p.at(e2.first), p.at(e2.second) - Tr_t<double>(1.0));
                out += b;
            }
        }
        return out * p.at(tr.conformation);
    }
};



template<template<class...> class Tr>
Allosteric_Model_new_<Tr>::Allosteric_Model_new_(const Allosteric_Model_new &other)
    :new_d_{other.get_new_model_definition()},conformer_{other.get_conformers()},transitions_{other.get_transitions()},
      paramNames_{other.get_paramNames()},conductances_{other.get_conductances()},d_{other.get_model_definition()},
      connections_to_0_{other.connections_to_x(0.0)},connections_to_a_{other.connections_to_a()},
      connections_to_x_{other.connections_to_x(1.0)},connections_from_0_{connections_from_x(0.0)},
      connections_from_x_{other.connections_from_x(1.0)}
{}



inline std::vector<std::vector<std::size_t>>
get_connections_up(const M_Matrix<double> &Q) {
  std::vector<std::vector<std::size_t>> out(Q.nrows());
  for (std::size_t i = 0; i < Q.nrows(); ++i) {
    for (std::size_t j = 0; j < Q.ncols(); ++j)
      if (Q(i, j) > 0)
        out[i].push_back(j);
  }
  return out;
}

inline std::vector<std::vector<std::size_t>>
get_connections_down(const M_Matrix<double> &Q) {
  std::vector<std::vector<std::size_t>> out(Q.ncols());
  for (std::size_t i = 0; i < Q.nrows(); ++i) {
    for (std::size_t j = 0; j < Q.ncols(); ++j)
      if (Q(i, j) > 0)
        out[j].push_back(i);
  }
  return out;
}

class Microscopic_description {
public:
  typedef std::vector<std::size_t> mistate;

  M_Matrix<double> P_matrix() const {

    M_Matrix<double> out(size(), number_of_states());
    for (std::size_t i = 0; i < size(); ++i)
      for (std::size_t j = 0; j < number_of_states(); ++j)
        out(i, j) = ns_[i][j] * (1.0 / number_of_channels());
    return out;
  }

  M_Matrix<double> P_P2_matrix() const {
    std::size_t k = (number_of_states() * (number_of_states() + 3)) / 2;
    M_Matrix<double> out(size(), k);
    for (std::size_t i = 0; i < size(); ++i) {

      auto p = M_Matrix<double>(1, number_of_states());
      for (std::size_t j = 0; j < p.size(); ++j)
        p[j] = ns_[i][j] * (1.0 / number_of_channels());
      for (std::size_t j = 0; j < number_of_states(); ++j)
        out(i, j) = p[j];
      auto p2 = quadraticForm_XTX(p);
      //  auto p2=quadraticForm_XTX(p)*number_of_channels();
      for (auto j = 0ul; j < p2.size(); ++j) {
        out(i, number_of_states() + j) = p2[j];
      }
    }
    return out;
  }
  M_Matrix<double> P_P2_matrix(const M_Matrix<double> &Pmean,
                               const M_Matrix<double> &Pcov) const {
    //         auto P2=quadraticForm_XTX(Pmean)*number_of_channels()+Pcov;
    auto P2 = quadraticForm_XTX(Pmean) + Pcov / number_of_channels();
    std::size_t k = (number_of_states() * (number_of_states() + 3)) / 2;
    M_Matrix<double> out(1, k);
    for (std::size_t j = 0; j < number_of_states(); ++j)
      out[j] = Pmean[j];
    for (auto j = 0ul; j < P2.size(); ++j)
      out[number_of_states() + j] = P2[j];
    return out;
  }
  auto maximum_entropy(const M_Matrix<double> &Pmean, std::size_t maxiter,
                       double maxGradient, double max_landa,
                       double maxstep) const {
    typedef myOptional_t<M_Matrix<double>> Op;
    auto Pij = P_matrix();
    maxent::end_criteria end(maxiter, maxGradient, max_landa);
    auto res = maxent::optimize(Pmean, Pij, end, maxstep);
    if (!res)
      return Op(false, res.error());
    else
      return Op(res.value().pi());
  }
  auto maximum_entropy(const M_Matrix<double> &Pmean,
                       const M_Matrix<double> &Pcov, std::size_t maxiter,
                       double maxGradient, double maxlanda,
                       double maxstep) const {
    typedef myOptional_t<M_Matrix<double>> Op;
    auto Pij = P_P2_matrix();
    maxent::end_criteria end(maxiter, maxGradient, maxlanda);
    auto PP2 = P_P2_matrix(Pmean, Pcov);
    auto res = maxent::optimize(PP2, Pij, end, maxstep);
    if (!res)
      return Op(false, res.error());
    else
      return Op(res.value().pi());
  }

  typedef Microscopic_description self_type;
  struct Less {
    bool operator()(const mistate A, const mistate B) const {
      auto k = A.size();
      assert(A.size() == B.size());
      for (std::size_t i = k; i > 0; --i) {
        if ((A)[i - 1] < (B)[i - 1])
          return true;
        else if ((B)[i - 1] < (A)[i - 1])
          return false;
      }
      return false;
    }
  };
  auto &occupancy_vector_list() const { return ns_; }
  std::size_t
  number_of_common_states(const Microscopic_description &other) const {
    assert(number_of_channels() == other.number_of_channels());
    assert(number_of_states() == other.number_of_states());
    Microscopic_description const *bigger, *smaller;
    if (size() < other.size()) {
      bigger = &other;
      smaller = this;
    } else {
      bigger = this;
      smaller = &other;
    }
    std::size_t count = 0;
    for (auto &e : smaller->ns_)
      if (bigger->has_it(e))
        ++count;
    return count;
  }
  std::pair<std::size_t, std::size_t> number_of_common_states_closure(
      const Microscopic_description &other,
      const std::vector<std::vector<std::size_t>> &connections_to_x,
      const std::vector<std::vector<std::size_t>> &connections_from_x) const {
    assert(number_of_channels() == other.number_of_channels());
    assert(number_of_states() == other.number_of_states());
    Microscopic_description const *bigger, *smaller;
    if (size() < other.size()) {
      bigger = &other;
      smaller = this;
    } else {
      bigger = this;
      smaller = &other;
    }
    std::size_t count = 0;
    for (auto &e : smaller->ns_)
      if (bigger->has_it(e))
        ++count;
    auto closure = smaller->limit_set(connections_to_x, connections_from_x);
    for (auto &e : closure)
      if (bigger->has_it(e))
        ++count;
    auto n_ext = smaller->size() + closure.size();
    return std::make_pair(count, n_ext);
  }

  double connecting_ratio(const Microscopic_description &other) {
    auto n = number_of_common_states(other);
    return std::max((1.0 * n) / size(), (1.0 * n / other.size()));
  }

  double connecting_ratio_closure(
      const Microscopic_description &other,
      const std::vector<std::vector<std::size_t>> &connections_to_x,
      const std::vector<std::vector<std::size_t>> &connections_from_x) {
    auto [nc, nt] = number_of_common_states_closure(other, connections_to_x,
                                                    connections_from_x);
    return (1.0 * nc) / nt;
  }

  Microscopic_description() = default;

  template <class Model, class X>
  Microscopic_description(const Model &m, const X &x, std::size_t numChannels,
                          const M_Matrix<double> &Pmean,
                          const M_Matrix<double> &Pcov, double maxChi2)
      : N_{numChannels}, k_{Pmean.size()}, ns_{}, map_{} {
    extend(m, x, Pmean, Pcov, maxChi2);
  }

//  template <class Model, class X>
//  Microscopic_description
//  multinomial(const Model &/*m*/, [[maybe_unused]]const X &x, std::size_t numChannles,
//              const M_Matrix<double> &Pmean, [[maybe_unused]]double lost_p) const {
//    auto d = multinomial_distribution(numChannles, Pmean.toVector());
//    return d;
//  }

  template <class Model, class X>
  void extend(const Model &m, const X &x, const M_Matrix<double> &Pmean,
              const M_Matrix<double> &Pcov, double maxChi2) {
    auto &connections_to_x = m.connections_to_x(x);
    auto &connections_from_x = m.connections_from_x(x);
    auto origin = to_n(Pmean);
    push_back(origin);

    auto final_state = get_states_as_required_unconnected(Pmean, Pcov, maxChi2);

    add_states_to_fill_gaps(final_state, connections_to_x, connections_from_x);
  }

  template <class Model>
  M_Matrix<double> extend_P(const Model &m, const M_Matrix<double> &P,
                            const M_Matrix<double> &Pmean,
                            const M_Matrix<double> &Pcov, double maxChi2) {
    extend(m, Pmean, Pcov, maxChi2);
    M_Matrix<double> out(1, size(), 0.0);
    for (std::size_t i = 0; i < P.size(); ++i)
      out[i] = P[i];
    return out;
  }

  void limit_set_of_state(
      std::set<mistate> &out, const mistate &origin,
      const std::vector<std::vector<std::size_t>> &connections_x) const {

    for (std::size_t i_origin = 0; i_origin < origin.size(); ++i_origin) {
      if (origin[i_origin] > 0)
        for (auto &i_destination : connections_x[i_origin]) {
          auto destination = origin;
          --destination[i_origin];
          ++destination[i_destination];
          if (!has_it(destination)) {
            out.insert(destination);
          }
        }
    }
  }

  std::set<mistate> limit_set(
      const std::vector<std::vector<std::size_t>> &connections_to_x,
      const std::vector<std::vector<std::size_t>> &connections_from_x) const {
    std::set<mistate> out;
    for (auto &e : ns_) {
      limit_set_of_state(out, e, connections_to_x);
      limit_set_of_state(out, e, connections_from_x);
    }
    return out;
  }

  auto build_new_layer(
      const std::set<mistate> &boundary,
      const std::set<mistate> &already_included,
      std::map<mistate, std::set<mistate>> &layer_connections,
      const std::vector<std::vector<std::size_t>> &connections_from_x) const {
    std::set<mistate> new_layer;
    for (auto &destination_ : boundary) {

      for (std::size_t i_destination = 0; i_destination < destination_.size();
           ++i_destination) {
        if (destination_[i_destination] > 0)
          for (auto &i_origin : connections_from_x[i_destination]) {
            auto origin_ = destination_;
            ++origin_[i_origin];
            --origin_[i_destination];
            if (already_included.find(origin_) == already_included.end()) {
              layer_connections[origin_].insert(destination_);
              new_layer.insert(origin_);
            }
          }
      }
    }
    return new_layer;
  }

  auto build_initial_layer(
      const std::set<std::size_t> &component,
      const std::vector<std::vector<std::size_t>> &connections_x) const {
    std::set<mistate> new_layer;
    for (auto origin_index : component) {
      auto &origin = ns_[origin_index];
      for (std::size_t i_origin = 0; i_origin < origin.size(); ++i_origin) {
        if (origin[i_origin] > 0)
          for (auto &i_destination : connections_x[i_origin]) {
            auto destination = origin;
            --destination[i_origin];
            ++destination[i_destination];
            if (!has_it(destination)) {
              new_layer.insert(destination);
            }
          }
      }
    }
    return new_layer;
  }

  static std::set<std::size_t>
  find_new_connections(std::size_t origin,
                       const std::set<std::size_t> &connected,
                       std::vector<std::vector<std::size_t>> &graph) {
    std::set<std::size_t> new_connected;
    for (std::size_t j = 0; j < graph[origin].size(); ++j) {
      auto d = graph[origin][j];
      if (connected.find(d) == connected.end())
        new_connected.insert(d);
    }
    return new_connected;
  }
  static bool is_connected(std::vector<std::vector<std::size_t>> &graph) {
    if (graph[0].empty())
      return false;
    else {
      std::set<std::size_t> connected{0ul};
      std::set<std::size_t> scouts{0ul};
      while (!scouts.empty()) {
        std::set<std::size_t> new_scouts;
        for (auto &e : scouts) {
          auto s = find_new_connections(e, connected, graph);
          connected.insert(s.begin(), s.end());
          new_scouts.insert(s.begin(), s.end());
        }
        scouts = std::move(new_scouts);
      }
      return connected.size() == graph.size();
    }
  }

  std::pair<bool, std::set<mistate>> get_component_joining_states(
      std::vector<std::set<mistate>> &lastLayers_per_component) const {
    std::vector<std::vector<std::size_t>> components_connected(
        lastLayers_per_component.size());
    std::set<mistate> common_all;
    for (std::size_t i = 0; i < lastLayers_per_component.size(); ++i) {
      for (std::size_t j = i + 1; j < lastLayers_per_component.size(); ++j) {
        std::set<mistate> common;
        std::set_intersection(lastLayers_per_component[i].begin(),
                              lastLayers_per_component[i].end(),
                              lastLayers_per_component[j].begin(),
                              lastLayers_per_component[j].end(),
                              std::inserter(common, common.end()));
        if (!common.empty())
          components_connected[i].push_back(j);
        common_all.insert(common.begin(), common.end());
      }
    }
    bool res = is_connected(components_connected);
    return std::make_pair(res, common_all);
  }

  // to be removed
  std::pair<std::vector<std::set<mistate>>, std::set<mistate>> joint_components(
      const std::vector<std::set<mistate>> &lastLayers_per_component) const {
    std::vector<std::set<mistate>> joint_layers;

    std::set<mistate> common_all;
    std::set<std::size_t> unchecked_components;
    for (std::size_t i = 0; i < lastLayers_per_component.size(); ++i)
      unchecked_components.insert(i);
    while (!unchecked_components.empty()) {
      auto it = unchecked_components.begin();

      std::set<std::size_t> other_unchecked_components;
      auto tested_component = lastLayers_per_component[*it];
      ++it;
      other_unchecked_components.insert(it, unchecked_components.end());
      bool all_others_checked = false;
      while (!other_unchecked_components.empty() && (!all_others_checked)) {
        auto it2 = other_unchecked_components.begin();
        bool is_separate = true;
        while ((it2 != other_unchecked_components.end()) && is_separate) {
          std::set<mistate> common;
          auto &other = lastLayers_per_component[*it2];
          std::set_intersection(
              tested_component.begin(), tested_component.end(), other.begin(),
              other.end(), std::inserter(common, common.end()));
          if (!common.empty()) {
            tested_component.insert(other.begin(), other.end());
            other_unchecked_components.erase(it2);
            unchecked_components.erase(*it2);
            common_all.insert(common.begin(), common.end());
            is_separate = false;
          } else {
            ++it2;
          }
        }
        if (is_separate)
          all_others_checked = true;
      }
      joint_layers.push_back(tested_component);
      unchecked_components.erase(unchecked_components.begin());
    }
    return std::make_pair(joint_layers, common_all);
  }

  std::tuple<std::vector<std::set<mistate>>, std::vector<std::set<mistate>>,
             std::vector<std::set<mistate>>, std::set<mistate>>
  joint_components(
      const std::vector<std::set<mistate>> &layer_per_component,
      const std::vector<std::set<mistate>> &components_bulk,
      const std::vector<std::set<mistate>> &components_core) const {
    std::vector<std::set<mistate>> joint_layers;
    std::vector<std::set<mistate>> joint_bulks;
    std::vector<std::set<mistate>> joint_core;

    std::set<mistate> common_all;
    std::set<std::size_t> unchecked_components;
    for (std::size_t i = 0; i < components_core.size(); ++i)
      unchecked_components.insert(i);
    while (!unchecked_components.empty()) {
      auto it = unchecked_components.begin();

      std::set<std::size_t> other_unchecked_components;
      auto tested_core = components_core[*it];
      auto tested_bulk = components_bulk[*it];
      auto tested_layer = layer_per_component[*it];
      ++it;
      other_unchecked_components.insert(it, unchecked_components.end());
      bool all_others_checked = false;
      while (!other_unchecked_components.empty() && (!all_others_checked)) {
        auto it2 = other_unchecked_components.begin();
        bool is_separate = true;
        while ((it2 != other_unchecked_components.end()) && is_separate) {
          auto &other_core = components_core[*it2];
          auto other_bulk = components_bulk[*it];
          auto &other_layer = layer_per_component[*it2];
          std::set<mistate> new_layer;
          std::set<mistate> new_bulk;
          std::set<mistate> new_core;

          std::set<mistate> common;
          std::set_intersection(tested_layer.begin(), tested_layer.end(),
                                other_bulk.begin(), other_bulk.end(),
                                std::inserter(common, common.end()));
          if (!common.empty()) {
            new_layer = other_layer;
            new_bulk = tested_bulk;
            new_core = other_core;
            new_bulk.insert(other_bulk.begin(), other_bulk.end());
            common_all.insert(common.begin(), common.end());
            is_separate = false;
          }
          common.clear();
          std::set_intersection(other_layer.begin(), other_layer.end(),
                                tested_bulk.begin(), tested_bulk.end(),
                                std::inserter(common, common.end()));
          if (!common.empty()) {
            new_layer.insert(tested_layer.begin(), tested_layer.end());
            new_bulk.insert(tested_bulk.begin(), tested_bulk.end());
            new_bulk.insert(other_bulk.begin(), other_bulk.end());
            new_core.insert(tested_core.begin(), tested_core.end());
            common_all.insert(common.begin(), common.end());
            is_separate = false;
          }
          common.clear();
          std::set_intersection(other_layer.begin(), other_layer.end(),
                                tested_core.begin(), tested_core.end(),
                                std::inserter(common, common.end()));
          std::set_intersection(tested_layer.begin(), tested_layer.end(),
                                other_core.begin(), other_core.end(),
                                std::inserter(common, common.end()));
          std::set_intersection(other_layer.begin(), other_layer.end(),
                                tested_layer.begin(), tested_layer.end(),
                                std::inserter(common, common.end()));

          if (!common.empty()) {
            new_layer.insert(common.begin(), common.end());
            new_bulk.insert(tested_bulk.begin(), tested_bulk.end());
            new_bulk.insert(other_bulk.begin(), other_bulk.end());
            new_core.insert(common.begin(), common.end());
            common_all.insert(common.begin(), common.end());
            is_separate = false;
          }
          if (!is_separate) {
            other_unchecked_components.erase(it2);
            unchecked_components.erase(*it2);
            tested_core = std::move(new_core);
            tested_bulk = std::move(new_bulk);
            tested_layer = std::move(new_layer);
          } else
            ++it2;
        }
        if (is_separate)
          all_others_checked = true;
      }
      joint_layers.push_back(tested_layer);
      joint_core.push_back(tested_core);
      joint_bulks.push_back(tested_bulk);
      unchecked_components.erase(unchecked_components.begin());
    }
    return std::make_tuple(joint_layers, joint_bulks, joint_core, common_all);
  }

  std::set<mistate> get_component_connecting_states(
      const std::set<mistate> &joints,
      const std::map<mistate, std::set<mistate>> &connections_per_layer) const {
    if (connections_per_layer.size() == 0)
      return joints;
    else {
      std::set<mistate> out;
      auto end_points = joints;
      while (!end_points.empty()) {

        std::set<mistate> new_points;
        for (auto &e : end_points) {
          if (out.find(e) == out.end()) {
            auto it = connections_per_layer.find(e);
            if (it != connections_per_layer.end()) {
              auto new_found = it->second;
              new_points.insert(new_found.begin(), new_found.end());
            }
          }
        }
        out.insert(end_points.begin(), end_points.end());
        end_points = std::move(new_points);
      }
      return out;
    }
  }

  std::set<mistate> connect_components(
      const std::vector<std::set<mistate>> &components_core,
      const std::vector<std::set<mistate>> &components_bulk,
      const std::vector<std::vector<std::size_t>> &connections_to_x) const {
    std::vector<std::set<mistate>> all_accepted = components_bulk;
    auto components_core_run = components_core;
    auto components_bulk_run = components_bulk;
    std::vector<std::set<mistate>> old_layer_per_component = components_core;
    std::set<mistate> all_joints;
    std::map<mistate, std::set<mistate>> connections;
    for (std::size_t i = 0; i < old_layer_per_component.size(); ++i)
      all_accepted[i].insert(old_layer_per_component[i].begin(),
                             old_layer_per_component[i].end());
    while (old_layer_per_component.size() > 1) {
      std::vector<std::set<mistate>> new_layer_per_components(
          old_layer_per_component.size());
      for (std::size_t i = 0; i < old_layer_per_component.size(); ++i) {
        new_layer_per_components[i] =
            build_new_layer(old_layer_per_component[i], all_accepted[i],
                            connections, connections_to_x);
      }
      std::vector<std::set<mistate>> new_layers;
      std::vector<std::set<mistate>> new_bulks;
      std::vector<std::set<mistate>> new_cores;
      std::set<mistate> joints_n;

      std::tie(new_layers, new_bulks, new_cores, joints_n) = joint_components(
          new_layer_per_components, components_bulk_run, components_core_run);

      assert(([&components_core, &components_bulk, &new_layers, &new_bulks,
               &new_cores, &joints_n, &new_layer_per_components,
               &old_layer_per_component, &all_joints, &connections,
               &connections_to_x, &all_accepted, &components_bulk_run,
               &components_core_run]() {
        if ((new_layers.size() == 1) ||
            (std::count_if(new_layer_per_components.begin(),
                           new_layer_per_components.end(),
                           [](auto &x) { return x.empty(); }) <= 1))
          return true;
        else {
          std::cerr
              << "new_layer_per_components is empty for some components!! ";
          std::cerr << "components_core=\n" << components_core << "\n";
          std::cerr << "components_bulk=\n" << components_bulk << "\n";
          std::cerr << "new_layer_per_components=\n"
                    << new_layer_per_components << "\n";
          std::cerr << "old_layer_per_component\n"
                    << old_layer_per_component << "\n";
          std::cerr << "new_layers\n" << new_layers << "\n";
          std::cerr << "new_bulks\n" << new_bulks << "\n";
          std::cerr << "components_bulk_run\n" << components_bulk_run << "\n";
          std::cerr << "new_cores\n" << new_cores << "\n";
          std::cerr << "components_core_run\n" << components_core_run << "\n";
          std::cerr << "joints_n\n" << joints_n << "\n";

          std::cerr << "all_joints\n" << all_joints << "\n";
          std::cerr << "all_accepted\n" << all_accepted << "\n";
          std::cerr << "connections\n" << connections << "\n";
          std::cerr << "connections_to_x\n" << connections_to_x << "\n";
          return false;
        }
      }()));
      all_joints.insert(joints_n.begin(), joints_n.end());
      old_layer_per_component = std::move(new_layers);
      components_bulk_run = std::move(new_bulks);
      components_core_run = std::move(new_cores);
      for (std::size_t i = 0; i < components_core_run.size(); ++i)
        components_core_run[i].insert(old_layer_per_component[i].begin(),
                                      old_layer_per_component[i].end());
      all_accepted = components_bulk_run;
      for (std::size_t i = 0; i < all_accepted.size(); ++i)
        all_accepted[i].insert(components_core_run[i].begin(),
                               components_core_run[i].end());
    }
    return get_component_connecting_states(all_joints, connections);
  }

  static void connect(
      const mistate &from, const mistate &to,
      std::map<mistate, std::pair<std::set<mistate>, std::set<mistate>>> &net) {

    for (auto &parent : net[from].first) {
      net[parent].second.insert(to);
    }
    for (auto &child : net[to].second) {
      net[child].first.insert(from);
    }
  }

  std::set<mistate> get_connecting_states(
      const std::set<mistate> &start, const std::set<mistate> &end,
      const std::vector<std::vector<std::size_t>> &connections_to_x,
      std::size_t max_extra_steps) const {
    auto surface = start;
    std::map<mistate, std::pair<std::set<mistate>, std::set<mistate>>> net;
    std::set<mistate> orfan_start = start;
    std::set<mistate> orfan_end = end;
    for (auto &e : start) {
      net[e] = std::make_pair(std::set<mistate>{e}, std::set<mistate>{e});
    }
    auto n = number_of_states();
    std::size_t extra_steps = max_extra_steps;
    while (extra_steps > 0 || !orfan_end.empty() /*|| !orfan_start.empty()*/) {
      if (orfan_end.empty() && orfan_start.empty())
        --extra_steps;
      assert(([&net, &start, &end, &surface, &orfan_start, &orfan_end,
               &connections_to_x]() {
        if (surface.empty()) {
          std::cerr << "\n start \n" << start << "\n";
          std::cerr << "\nend \n" << end << "\n";
          std::cerr << "\n orfan_start \n" << orfan_start << "\n";
          std::cerr << "\n orfan_end \n" << orfan_end << "\n";
          std::cerr << "\nnet \n" << net << "\n";
          std::cerr << "\nsurface \n" << surface << "\n";
          std::cerr << "\nconnections_to_x \n" << connections_to_x << "\n";

          return false;
        } else
          return true;
      }()));
      std::set<mistate> new_surface;
      for (auto &origin : surface) {
        for (std::size_t i_origin = 0; i_origin < n; ++i_origin) {
          if (origin[i_origin] > 0) {
            for (auto i_destination : connections_to_x[i_origin]) {
              auto destination = origin;
              destination[i_origin]--;
              destination[i_destination]++;
              auto it = net.find(destination);
              if (it == net.end()) {
                net[destination] =
                    std::make_pair(std::set<mistate>{destination},
                                   std::set<mistate>{destination});
                for (auto &parent : net[origin].first) {
                  net[parent].second.insert(destination);
                }
                new_surface.insert(destination);
              } else {
                for (auto &parent : net[origin].first)
                  net[parent].second.insert(net[destination].second.begin(),
                                            net[destination].second.end());
                for (auto &child : net[destination].second)
                  net[child].first.insert(net[origin].first.begin(),
                                          net[origin].first.end());
              }
            }
          }
        }
      }
      std::set<mistate> new_orfan_start;
      for (auto &s : orfan_start)
        if (std::find_if(net[s].second.begin(), net[s].second.end(),
                         [&end](auto &x) {
                           return end.find(x) != end.end();
                         }) == net[s].second.end())
          new_orfan_start.insert(s);

      std::set<mistate> new_orfan_end;
      for (auto &e : orfan_end)
        if (net.find(e) == net.end())
          new_orfan_end.insert(e);
      orfan_start = std::move(new_orfan_start);
      orfan_end = std::move(new_orfan_end);

      surface = std::move(new_surface);
    }
    std::set<mistate> out;
    for (auto &end_state : end) {
      out.insert(net[end_state].first.begin(), net[end_state].first.end());
    }
    return out;
  }
  auto
  connected_set(const mistate &origin,
                const std::vector<std::vector<std::size_t>> &connections_from_x,
                std::set<mistate> &&already_connected) const {
    std::set<mistate> limit;
    already_connected.insert(origin);

    connected_search(origin, already_connected, limit, connections_from_x);

    while (!limit.empty()) {
      std::set<mistate> next_layer;
      for (auto &e : limit) {
        connected_search(e, already_connected, next_layer, connections_from_x);
      }
      limit = std::move(next_layer);
    }
    return std::move(already_connected);
  }

  auto connected_set(
      const mistate &origin,
      const std::vector<std::vector<std::size_t>> &connections_from_x) const {
    std::set<mistate> already_connected;
    return connected_set(origin, connections_from_x,
                         std::move(already_connected));
  }

  std::pair<std::vector<std::set<mistate>>, std::vector<std::set<mistate>>>
  connected_components(
      const std::vector<std::vector<std::size_t>> &connections_to_x,
      const std::vector<std::vector<std::size_t>> &connections_from_x) const {
    std::vector<std::set<mistate>> components_core;
    std::vector<std::set<mistate>> components_bulk;

    auto first_bulk = connected_set(ns_[0], connections_from_x);
    auto first_core = connected_set(ns_[0], connections_to_x);

    components_core.push_back(first_core);
    components_bulk.push_back(first_bulk);
    for (std::size_t origin_index = 1; origin_index < size(); ++origin_index) {
      auto origin = ns_[origin_index];
      bool already_in_bulk = false;
      for (auto &e : components_bulk)
        if (e.find(origin) != e.end()) {
          already_in_bulk = true;
          break;
        }
      if (!already_in_bulk) {
        bool already_in_core = false;
        std::vector<std::set<mistate>> components_core_new;
        std::vector<std::set<mistate>> components_bulk_new;
        for (std::size_t i = 0; i < components_core.size(); ++i) {
          if (components_core[i].find(origin) != components_core[i].end()) {
            if (!already_in_core) {
              auto new_bulk = connected_set(origin, connections_from_x,
                                            std::move(components_bulk[i]));
              auto new_core = connected_set(origin, connections_to_x);
              components_core_new.push_back(std::move(new_core));
              components_bulk_new.push_back(std::move(new_bulk));

              already_in_core = true;
            }
          } else {
            components_core_new.push_back(std::move(components_core[i]));
            components_bulk_new.push_back(std::move(components_bulk[i]));
          }
        }
        if (!already_in_core) {
          auto new_bulk = connected_set(origin, connections_from_x);
          auto new_core = connected_set(origin, connections_to_x);
          components_core_new.push_back(std::move(new_core));
          components_bulk_new.push_back(std::move(new_bulk));
        }
        components_bulk = std::move(components_bulk_new);
        components_core = std::move(components_core_new);
      }
    }
    return std::make_pair(components_core, components_bulk);
  }

  void add_states_to_fill_gaps(
      const std::set<mistate> final,
      const std::vector<std::vector<std::size_t>> &connections_to_x,
      [[maybe_unused]]const std::vector<std::vector<std::size_t>> &connections_from_x) {
    //    auto [components_core, components_bulk] =
    //        connected_components(connections_to_x, connections_from_x);
    //    if (components_core.size() > 1) {
    //            std::cerr << "components_core" << components_core << "\n";
    //            std::cerr << "components_bulk" << components_bulk << "\n";

    std::set<mistate> start(ns_.begin(), ns_.end());
    auto filling_states =
        get_connecting_states(start, final, connections_to_x, 0);
    for (auto &e : filling_states)
      push_back(e);
    //    assert(connected_components(connections_to_x, connections_from_x)
    //               .first.size() == 1);
  }

  void connected_search(
      const mistate &origin, std::set<mistate> &already_connected,
      std::set<mistate> &newly_connected,
      const std::vector<std::vector<std::size_t>> &connections_from_x) const {
    for (std::size_t i_origin = 0; i_origin < origin.size(); ++i_origin) {
      if (origin[i_origin] > 0)
        for (auto &i_destination : connections_from_x[i_origin]) {
          auto destination = origin;
          --destination[i_origin];
          ++destination[i_destination];
          if (has_it(destination)) {
            if (already_connected.find(destination) ==
                already_connected.end()) {
              newly_connected.insert(destination);
              already_connected.insert(destination);
            }
          }
        }
    }
  }

  M_Matrix<double>
  extend_P(const M_Matrix<double> &P, const M_Matrix<double> &Pmean,
           const M_Matrix<double> &Pcov,
           const std::vector<std::vector<std::size_t>> &connections_to_x,
           const std::vector<std::vector<std::size_t>> &connections_from_x,
           double maxChi2) {
    auto final = get_states_as_required_unconnected(Pmean, Pcov, maxChi2);
    add_states_to_fill_gaps(final, connections_to_x, connections_from_x);
    M_Matrix<double> out(1, size(), 0.0);
    for (std::size_t i = 0; i < P.size(); ++i)
      out[i] = P[i];
    return out;
  }

  Microscopic_description &merge(Microscopic_description &&other) {
    for (auto &&e : std::move(other.ns_))
      push_back(std::move(e));
    return *this;
  }
  Microscopic_description &merge(const Microscopic_description &other) {
    for (auto &e : other.ns_)
      push_back(e);
    return *this;
  }

  M_Matrix<double>
  merge_P(const M_Matrix<double> &P, const Microscopic_description &other,
          const std::vector<std::vector<std::size_t>> &connections_to_x,
          const std::vector<std::vector<std::size_t>> &connections_from_x) {
    std::set<mistate> final(other.ns_.begin(), other.ns_.end());
    add_states_to_fill_gaps(final, connections_to_x, connections_from_x);
    M_Matrix<double> out(1, size(), 0.0);
    for (std::size_t i = 0; i < P.size(); ++i)
      out[i] = P[i];

    return out;
  }
  auto search_neighbours_up(
      std::set<mistate> &rejected, std::set<mistate> &accepted,
      const mistate &origin,
      const std::vector<std::vector<std::size_t>> &connections_x,
      const std::vector<double> &Nmean, const M_Matrix<double> &pinvNCov,
      double maxChi2) const {
    std::vector<mistate> out;
    std::set<mistate> new_rejected;
    for (std::size_t i_origin = 0; i_origin < origin.size(); ++i_origin) {
      if (origin[i_origin] > 0)
        for (auto &i_destination : connections_x[i_origin]) {
          auto destination = origin;
          --destination[i_origin];
          ++destination[i_destination];
          if (!has_it(destination) &&
              (rejected.find(destination) == rejected.end()) &&
              (accepted.find(destination) == accepted.end())) {
            if (has_unallowed_states(destination, Nmean)) {
              new_rejected.insert(destination);
            } else {
              double logL = chi2sum(destination, Nmean, pinvNCov);
              if (logL < maxChi2) {
                out.push_back(destination);
                accepted.insert(destination);
              } else
                new_rejected.insert(destination);
            }
          }
        }
    }
    rejected.insert(new_rejected.begin(), new_rejected.end());
    return out;
  }

  std::set<mistate>
  search_on_normal(const mistate &origin,
                   const std::vector<std::vector<std::size_t>> &connections_x,
                   const std::vector<double> &Nmean,
                   const M_Matrix<double> &pinvNCov, double maxChi2) const {
    std::set<mistate> rejected;
    std::set<mistate> accepted;
    accepted.insert(origin);
    auto neighbours_list = search_neighbours_up(
        rejected, accepted, origin, connections_x, Nmean, pinvNCov, maxChi2);
    while (!neighbours_list.empty()) {

      std::vector<std::vector<std::size_t>> newNeigh;
      for (auto &e : neighbours_list) {
        auto new_neighbour_up = search_neighbours_up(
            rejected, accepted, e, connections_x, Nmean, pinvNCov, maxChi2);
        newNeigh.insert(newNeigh.end(), new_neighbour_up.begin(),
                        new_neighbour_up.end());
      }
      neighbours_list = newNeigh;
    }
    return accepted;
  }

  std::set<mistate> get_states_as_required(
      const M_Matrix<double> &Pmean, const M_Matrix<double> &Pcov,
      const std::vector<std::vector<std::size_t>> &connections_x,
      double maxChi2) const {
    std::vector<double> Nmean = (Pmean * number_of_channels()).toVector();
    auto origin = to_n(Pmean);
    std::set<mistate> accepted;
    accepted.insert(origin);
    if (maxAbs(Pcov) > std::numeric_limits<double>::epsilon() * 1000) {
      auto pinvNCov = pinv(Pcov * number_of_channels());
      if (pinvNCov.has_value() && (maxAbs(pinvNCov.value()) > 0)) {
        accepted = search_on_normal(origin, connections_x, Nmean,
                                    pinvNCov.value(), maxChi2);
        //  search_down(origin, connections_down, Nmean, pinvNCov.value(),
        //  maxChi2);
      }
    }
    assert(([&accepted, &Pmean, &Pcov, &maxChi2, this]() {
      auto mifull = add_states_as_required_full_normal(
          number_of_channels(), number_of_states(), Pmean, Pcov, maxChi2);
      auto set_dif = symmetric_difference(mifull, accepted);
      if (!set_dif.empty()) {
        std::cerr << "\nPmean =" << Pmean << "\nPcov= " << Pcov << "\nmaxChi2"
                  << maxChi2 << "\nset_dif" << set_dif << "\n smart=" << ns_
                  << "\n full=" << mifull.ns_ << "accepted \n"
                  << accepted << "\n";

        return false;
      }
      return true;
    }()));
    return accepted;
  }

  std::set<mistate>
  get_states_as_required_unconnected(const M_Matrix<double> &Pmean,
                                     const M_Matrix<double> &Pcov,
                                     double maxChi2) const {
    std::vector<std::vector<std::size_t>> connections_x(number_of_states());
    for (std::size_t i = 0; i < number_of_states(); ++i)
      for (std::size_t j = 0; j < number_of_states(); ++j)
        if (j != i)
          connections_x[i].push_back(j);
    auto out = get_states_as_required(Pmean, Pcov, connections_x, maxChi2);

    assert(([this, &Pmean, &Pcov, &maxChi2, &out]() {
      auto mifull = add_states_as_required_full_normal(
          number_of_channels(), number_of_states(), Pmean, Pcov, maxChi2);
      auto set_dif = symmetric_difference(mifull, out);
      if (!set_dif.empty()) {
        std::cerr << "\nPmean =" << Pmean << "\nPcov= " << Pcov << "\nmaxChi2"
                  << maxChi2 << "\nset_dif" << set_dif << "\n smart=" << ns_
                  << "\n full=" << mifull.ns_;
        return false;
      }
      return true;
    }()));
    return out;
  }

  static Microscopic_description add_states_as_required_full_normal(
      std::size_t number_of_channels__, std::size_t number_of_states__,
      const M_Matrix<double> &Pmean, const M_Matrix<double> &Pcov,
      double maxChi2) {
    Microscopic_description work =
        fullMap(number_of_channels__, number_of_states__);
    Microscopic_description out(number_of_channels__, number_of_states__);
    std::vector<double> Nmean = (Pmean * number_of_channels__).toVector();
    auto pinvNCov = pinv(Pcov * number_of_channels__);
    if (!pinvNCov) {
      auto ns = out.to_n(Pmean);
      out.push_back(ns);
      return out;
    } else {
      for (auto &e : work.ns_) {
        if ((chi2sum(e, Nmean, pinvNCov.value()) < maxChi2) &&
            !has_unallowed_states(e, Nmean))
          out.push_back(e);
      }
      return out;
    }
  }

  std::vector<std::vector<std::size_t>>
  symmetric_difference(const Microscopic_description &other) const {
    std::vector<std::vector<std::size_t>> out;
    for (auto &e : ns_) {
      if (!other.has_it(e))
        out.push_back(e);
    }

    for (auto &e : other.ns_) {
      if (!has_it(e))
        out.push_back(e);
    }
    return out;
  }
  std::vector<mistate>
  symmetric_difference(const Microscopic_description &other,
                       std::set<mistate> me) const {
    std::vector<mistate> out;
    for (auto &e : me) {
      if (!other.has_it(e))
        out.push_back(e);
    }

    for (auto &e : other.ns_) {
      if (me.find(e) == me.end())
        out.push_back(e);
    }
    return out;
  }

  auto Qs(const M_Matrix<double> &Q0m, const M_Matrix<double> &Qam) const {
    M_Matrix<double> Q0(size(), size(), Matrix_TYPE::FULL, 0);
    M_Matrix<double> Qa(size(), size(), Matrix_TYPE::FULL, 0);
    for (std::size_t i = 0; i < size(); ++i) {
      for (std::size_t j = i + 1; j < size(); ++j) {
        std::size_t i_origin = number_of_states();
        std::size_t i_destination = number_of_states();
        bool invalid = false;
        auto origin = ns_[i];
        auto destination = ns_[j];
        for (std::size_t k = 0; k < number_of_states(); ++k) {
          if (origin[k] + 1 == destination[k]) {
            if (i_destination == number_of_states())
              i_destination = k;
            else
              invalid = true;
          } else if (origin[k] == destination[k] + 1) {
            if (i_origin == number_of_states())
              i_origin = k;
            else
              invalid = true;
          } else if (origin[k] != destination[k])
            invalid = true;
          if (invalid)
            break;
        }
        if (!invalid) {
          Q0(i, j) = origin[i_origin] * Q0m(i_origin, i_destination);
          Q0(i, i) -= Q0(i, j);
          Q0(j, i) = destination[i_destination] * Q0m(i_destination, i_origin);
          Q0(j, j) -= Q0(j, i);
          Qa(i, j) = origin[i_origin] * Qam(i_origin, i_destination);
          Qa(i, i) -= Qa(i, j);
          Qa(j, i) = destination[i_destination] * Qam(i_destination, i_origin);
          Qa(j, j) -= Qa(j, i);
        }
      }
    }
    return std::make_tuple(Q0, Qa);
  }

  auto Qs(const M_Matrix<double> &Q0m, const M_Matrix<double> &Qam,
          const std::vector<std::vector<std::size_t>> &connections_to_0,
          const std::vector<std::vector<std::size_t>> &connections_to_a) const {
    M_Matrix<double> Q0(size(), size(), Matrix_TYPE::FULL, 0);
    M_Matrix<double> Qa(size(), size(), Matrix_TYPE::FULL, 0);
    for (std::size_t i = 0; i < size(); ++i) {
      auto &origin = ns_[i];
      auto n = number_of_states();
      for (std::size_t i_origin = 0; i_origin < n; ++i_origin) {
        if (origin[i_origin] > 0) {
          for (auto i_destination : connections_to_0[i_origin]) {
            auto destination = origin;
            destination[i_origin]--;
            destination[i_destination]++;
            auto j = index(destination);
            if (j < size()) {
              Q0(i, j) = origin[i_origin] * Q0m(i_origin, i_destination);
              Q0(i, i) -= Q0(i, j);
            }
          }

          for (auto i_destination : connections_to_a[i_origin]) {
            auto destination = origin;
            destination[i_origin]--;
            destination[i_destination]++;
            auto j = index(destination);
            if (j < size()) {
              Qa(i, j) = origin[i_origin] * Qam(i_origin, i_destination);
              Qa(i, i) -= Qa(i, j);
            }
          }
        }
      }
    }
    assert(([&Q0, &Qa, &Q0m, &Qam, this]() {
      auto [Q0t, Qat] = this->Qs(Q0m, Qam);
      return are_Equal_v(Q0, Q0t, std::numeric_limits<double>::epsilon(),
                         std::numeric_limits<double>::epsilon(), std::cerr) &&
             are_Equal_v(Qa, Qat, std::numeric_limits<double>::epsilon(),
                         std::numeric_limits<double>::epsilon(), std::cerr);
    }()));

    return std::make_pair(Q0, Qa);
  }

  auto g(const M_Matrix<double> &gm) const {
    assert(gm.size() == number_of_states());
    M_Matrix<double> out(size(), 1, 0.0);
    for (std::size_t i = 0; i < size(); ++i) {
      double sum = 0;
      for (std::size_t j = 0; j < gm.size(); ++j)
        sum += gm[j] * ns_[i][j];
      out[i] = sum;
    }
    return out;
  }

  // auto Qrun(const M_Matrix<double> &Qrun,
  //          const std::vector<std::vector<std::size_t>> &connections_x) {
  //  M_Matrix<double> Qr(size(), size(), Matrix_TYPE::FULL, 0);
  //  for (auto &e : map_) {
  //    std::size_t i = e.second;
  //    auto &origin = e.first;
  //    auto n = origin.size();
  //    for (std::size_t i_origin = 0; i_origin < n; ++i_origin) {
  //      if (origin[i_origin] > 0) {
  //        for (auto i_destination : connections_x[i_origin]) {
  //          auto destination = origin;
  //          destination[i_origin]--;
  //          destination[i_destination]++;
  //          auto j = index(destination);
  //          if (j < size()) {
  //            Qr(i, j) = origin[i_origin] * Qrun(i_origin, i_destination);
  //            Qr(i, i) -= Qr(i, j);
  //            Qr(j, i) =
  //                destination[i_destination] * Qrun(i_destination,
  //                i_origin);
  //            Qr(j, j) -= Qr(j, i);
  //          }
  //        }
  //      }
  //    }
  //  }
  //  return Qr;
  //}

  Microscopic_description(std::size_t number_of_channels,
                          std::vector<std::vector<std::size_t>> &&ns)
      : N_{number_of_channels}, k_{ns[0].size()}, ns_{std::move(ns)},
        map_{get_map(ns_)} {}
  Microscopic_description(std::size_t number_of_channels,
                          const std::vector<std::vector<std::size_t>> &ns)
      : N_{number_of_channels}, k_{ns[0].size()}, ns_{ns}, map_{get_map(ns_)} {}
  Microscopic_description(std::size_t number_of_channels,
                          std::size_t number_of_states_)
      : N_{number_of_channels}, k_{number_of_states_}, ns_{}, map_{} {}

  static auto get_constructor_fields() {
    return std::make_tuple(grammar::field(C<self_type>{},
                                          "occupancy_vector_list",
                                          &self_type::occupancy_vector_list));
  }
  void fillMap() {
    while (map_.size() < ns_.size()) {
      auto i = map_.size();
      map_.emplace(ns_[i], i);
    }
  }
  auto convert_to_Macroscopic(const M_Matrix<double> &P) const {
    assert(size() == P.size());
    M_Matrix<double> Nmean(1, number_of_states(), 0.0);
    M_Matrix<double> Nsqrmean(number_of_states(), number_of_states(),
                              Matrix_TYPE::SYMMETRIC, 0);

    for (std::size_t i = 0; i < size(); ++i) {
      M_Matrix<double> p =
          M_Matrix<std::size_t>(1ul, number_of_states(), ns_[i]);
      Nmean += p * P[i];
      Nsqrmean += quadraticForm_XTX(p) * P[i];
    }
    M_Matrix<double> Pmean = Nmean * (1.0 / number_of_channels());
    M_Matrix<double> Pcov =
        (Nsqrmean - quadraticForm_XTX(Nmean)) * (1.0 / number_of_channels());

    using namespace std::literals::string_literals;
    //  auto PP=convert_to_Microscopic(Pmean,Pcov);
    //    are_Equal_v(P,PP,
    //                       std::numeric_limits<double>::epsilon()*1e10,
    //                       std::numeric_limits<double>::epsilon()*1e10,std::cerr,
    //                       "P"s,P,"ns",M_Matrix<std::size_t>(ns_),"convert_to_Microscopic(Pmean,Pcov)"s,convert_to_Microscopic(Pmean,Pcov),
    //                       "Pmean"s,Pmean,"Pcov"s,Pcov);

    return std::make_tuple(Pmean, Pcov);
  }

  static void next(std::vector<std::size_t> &n) {
    std::size_t i_origin = 0;
    std::size_t i_destination = 1;
    if (n[i_origin] > 0) {
      --n[i_origin];
      ++n[i_destination];
    } else {
      ++i_origin;
      while (n[i_origin] == 0)
        ++i_origin;
      i_destination = i_origin + 1;
      --n[i_origin];
      ++n[i_destination];
      if (n[i_origin] > 0) {
        n[0] = n[i_origin];
        n[i_origin] = 0;
      }
    }
  }
  static Microscopic_description fullMap(std::size_t Number_of_Channels__,
                                         std::size_t number_of_states__) {

    auto k = number_of_states__;
    Microscopic_description out(Number_of_Channels__, number_of_states__);
    std::vector<std::vector<std::size_t>> ns;
    std::vector<std::size_t> start(k, 0);
    std::vector<std::size_t> end(k, 0);
    start[0] = Number_of_Channels__;
    end.back() = Number_of_Channels__;
    auto v = start;
    out.push_back(v);
    while (Less{}(v, end)) {
      next(v);
      out.push_back(v);
    }

    return out;
  }

  M_Matrix<double>
  convert_to_Microscopic_normal(const M_Matrix<double> &Pmean,
                                const M_Matrix<double> &Pcov) const {
    double sumlik = 0;
    if (maxAbs(Pcov) > 0) {
      auto pinvNCov = pinv(Pcov * number_of_channels()).value();
      if (maxAbs(pinvNCov)) {
        std::vector<double> Nmean = (Pmean * number_of_channels()).toVector();
        M_Matrix<double> lik(1, size());

        for (std::size_t i = 0; i < size(); ++i) {
          double l = std::exp(-0.5 * chi2sum(ns_[i], Nmean, pinvNCov));
          lik[i] = l;
          sumlik += l;
        }
        M_Matrix<double> P = lik / sumlik;
        return P;
      }
    }
    auto center = to_n(Pmean);
    M_Matrix<double> P(1, size(), 0.0);
    P[index(center)] = 1.0;
    return P;
  }
  auto convert_to_Microscopic_maxent(const M_Matrix<double> &Pmean,
                                     const M_Matrix<double> &Pcov) const {
    if (size() == 1)
      return myOptional_t<M_Matrix<double>>(M_Matrix<double>(1, 1, 1.0));
    else
      return maximum_entropy(Pmean, Pcov, 100, 1e-7, 16, 1.0);
  }

  Microscopic_description(
      std::size_t number_of_channels, const M_Matrix<double> &Pmean,
      const M_Matrix<double> &Pcov,
      [[maybe_unused]] const std::vector<std::vector<std::size_t>> &connections_to_x,
      [[maybe_unused]]const std::vector<std::vector<std::size_t>> &connections_from_x,
      double maxChi2)
      : N_{number_of_channels}, k_{Pmean.size()}, ns_{}, map_{} {
    auto origin = to_n(Pmean);
    push_back(origin);
    auto final = get_states_as_required_unconnected(Pmean, Pcov, maxChi2);
    for (auto &e : final)
      push_back(e);
    // add_states_to_fill_gaps(final,connections_to_x, connections_from_x);
    if constexpr (false) {
        assert((([this, &Pmean, &Pcov]() ->bool{
        auto P = convert_to_Microscopic_maxent(Pmean, Pcov);
        if (!P) {
          std::cerr << P.error();
          return false;
        } else {
          auto [Pmean_r, Pcov_r] = convert_to_Macroscopic(P.value());
          return are_Equal_v(Pmean, Pmean_r, 0.01, 0.01, std::cerr, "Pmean",
                             Pmean, "PCov", Pcov, "ns",
                             M_Matrix<std::size_t>(ns_)) &&
                 are_Equal_v(Pcov, Pcov_r, 0.05, 0.05, std::cerr, "Pmean",
                             Pmean, "Pmean_r", Pmean_r, "PCov", Pcov, "ns",
                             M_Matrix<std::size_t>(ns_))

              ;
        }
        }())));
    }
  }

  std::size_t size() const { return ns_.size(); }
  std::size_t number_of_channels() const { return N_; }
  bool has_it(const std::vector<std::size_t> &n) const {
    assert(n.size() == number_of_states());
    assert(sum(n) == number_of_channels());
    return map_.find(n) != map_.end();
  }
  std::size_t index(const std::vector<std::size_t> &n) const {
    auto it = map_.find(n);
    if (it != map_.end())
      return it->second;
    else
      return size();
  }

  void push_back(const mistate &n) {
    if (!has_it(n)) {
      ns_.push_back(n);
      map_[ns_.back()] = ns_.size() - 1;
    }
  }
  void push_back(mistate &&n) {
    if (!has_it(n)) {
      ns_.push_back(std::move(n));
      map_[ns_.back()] = ns_.size() - 1;
    }
  }

  std::size_t number_of_states() const { return k_; }
  static bool has_unallowed_states(
      const std::vector<std::size_t> &n, const std::vector<double> &Nmean,
      double eps = std::numeric_limits<double>::epsilon() * 100) {
    for (std::size_t i = 0; i < n.size(); ++i) {
      if ((Nmean[i] < eps) && (n[i] > 0))
        return true;
    }
    return false;
  }

  static double chi2sum(const std::vector<std::size_t> &n,
                        const std::vector<double> &Nmean,
                        const M_Matrix<double> &pinvNCov) {
    return xTSigmaX(Nmean - n, pinvNCov);
  }

  mistate to_n(const M_Matrix<double> &Pmean) const {
    std::vector<std::size_t> out(Pmean.size());
    std::multimap<double, std::size_t> remain;
    std::size_t sum = 0;
    for (std::size_t i = 0; i < out.size(); ++i) {
      double n = number_of_channels() * Pmean[i];
      out[i] = n;
      sum += out[i];
      remain.insert(std::make_pair(n - out[i], i));
    }
    auto it = remain.rbegin();
    while (sum < number_of_channels()) {
      auto i = it->second;
      ++out[i];
      ++sum;
      ++it;
    }
    return out;
  }

  void reduce(M_Matrix<double> &P, double remove_p) {
    assert(P.size() == size());
    remove_low_probabilities(ns_, P, remove_p);
    map_ = get_map(ns_);
    assert(P.size() == size());
  }

private:
  bool
  is_connected_up(const std::vector<std::size_t> &origin,
                  const std::vector<std::vector<std::size_t>> &connections_x) {
    for (std::size_t i_origin = 0; i_origin < origin.size(); ++i_origin) {
      if (origin[i_origin] > 0)
        for (auto &i_destination : connections_x[i_origin]) {
          auto destination = origin;
          --destination[i_origin];
          ++destination[i_destination];
          if (has_it(destination))
            return true;
        }
    }
    return false;
  }

  bool is_connected_down(
      const std::vector<std::size_t> &destination,
      const std::vector<std::vector<std::size_t>> &connections_down) {
    for (std::size_t i_destination = 0; i_destination + 1 < destination.size();
         ++i_destination) {
      if (destination[i_destination] < number_of_channels())
        for (auto &i_origin : connections_down[i_destination]) {
          if (destination[i_origin] > 0) {
            auto origin = destination;
            --origin[i_origin];
            ++origin[i_destination];
            if (has_it(origin)) {
              return true;
            }
          }
        }
    }
    return false;
  }

  std::size_t N_;
  std::size_t k_;
  std::vector<std::vector<std::size_t>> ns_;
  std::map<std::vector<std::size_t>, std::size_t, Less> map_;

  static std::map<std::vector<std::size_t>, std::size_t, Less>
  get_map(const std::vector<std::vector<std::size_t>> &ns) {
    std::map<std::vector<std::size_t>, std::size_t, Less> out;
    for (std::size_t i = 0; i < ns.size(); ++i) {
      out[ns[i]] = i;
    }
    return out;
  }
};

template <class Model> class Microscopic_Model_New {
public:
  template <class Parameters> auto Qs(const Parameters &p) const {
    auto [Q0m, Qam] = m_.Qs(p);
    auto &connections_to_0 = m_.connections_to_0();
    auto &connections_to_a = m_.connections_to_a();
    return mi_.Qs(Q0m, Qam, connections_to_0, connections_to_a);
  }
  auto &microscopic_description() const { return mi_; }

  template <class Parameters> auto g(const Parameters &p) const {
    auto gm = m_.g(p);
    return mi_.g(gm);
  }

  std::size_t k() const { return mi_.size(); }
  Microscopic_Model_New(const Model &m, const Microscopic_description &mi)
      : m_{m}, mi_{mi} {}
  Microscopic_Model_New(const Model &m, Microscopic_description &&mi)
      : m_{m}, mi_{std::move(mi)} {}

private:
  Model m_;
  Microscopic_description mi_;
};

template <class Model> class Microscopic_Model {
public:
  struct Less {
    bool operator()(const std::vector<std::size_t> &A,
                    const std::vector<std::size_t> &B) const {
      auto k = A.size();
      assert(A.size() == B.size());
      for (std::size_t i = k; i > 0; --i) {
        if (A[i - 1] < B[i - 1])
          return true;
        else if (B[i - 1] < A[i - 1])
          return false;
      }
      return false;
    }
  };

  typedef std::map<std::vector<std::size_t>, std::size_t, Less> myMap;

  template <class Parameters> auto Qs(const Parameters &p) const {

    auto [Q0m, Qam] = m_.Qs(p);

    auto gm = m_.g(p);

    auto &connections_to_0 = m_.connections_to_0();
    auto &connections_to_a = m_.connections_to_a();

    M_Matrix<double> Q0(k(), k(), Matrix_TYPE::FULL, 0);
    M_Matrix<double> Qa(k(), k(), Matrix_TYPE::FULL, 0);
    for (auto &e : map_) {
      std::size_t i = e.second;
      std::vector<std::size_t> origin = e.first;
      auto n = origin.size();
      for (std::size_t i_origin = 0; i_origin < n; ++i_origin) {
        if (origin[i_origin] > 0) {
          for (auto i_destination : connections_to_0[i_origin]) {
            auto destination = origin;
            destination[i_origin]--;
            destination[i_destination]++;
            auto it = map_.find(destination);
            if (it != map_.end()) {
              auto j = it->second;
              Q0(i, j) = origin[i_origin] * Q0m(i_origin, i_destination);
              Q0(i, i) -= Q0(i, j);
            }
          }

          for (auto i_destination : connections_to_a[i_origin]) {
            auto destination = origin;
            destination[i_origin]--;
            destination[i_destination]++;
            auto it = map_.find(destination);
            if (it != map_.end()) {
              auto j = it->second;
              Qa(i, j) = origin[i_origin] * Qam(i_origin, i_destination);
              Qa(i, i) -= Qa(i, j);
              Qa(j, i) =
                  destination[i_destination] * Qam(i_destination, i_origin);
              Qa(j, j) -= Qa(j, i);
            }
          }
        }
      }
    }
    return std::pair(Q0, Qa);
  }

  template <class Parameters> auto g(const Parameters &p) const {
    M_Matrix<double> out(k(), 1, 0.0);
    auto gm = m_.g(p);
    for (auto &e : map_) {
      double sum = 0;
      for (std::size_t i = 0; i < gm.size(); ++i)
        sum += gm[i] * e.first[i];
      out[e.second] = sum;
    }
    return out;
  }

  std::size_t k() const { return map_.size(); }

  Microscopic_Model(const Model &m, std::size_t Num_of_Channels)
      : m_{m}, ns_{}, map_{} {
    std::tie(ns_, map_) = fullMap(m, Num_of_Channels);
  }

  static void next(std::vector<std::size_t> &n) {
    std::size_t i_origin = 0;
    std::size_t i_destination = 1;
    if (n[i_origin] > 0) {
      --n[i_origin];
      ++n[i_destination];
    } else {
      ++i_origin;
      while (n[i_origin] == 0)
        ++i_origin;
      i_destination = i_origin + 1;
      --n[i_origin];
      ++n[i_destination];
      if (n[i_origin] > 0) {
        n[0] = n[i_origin];
        n[i_origin] = 0;
      }
    }
  }
  static auto fullMap(const Model &m, std::size_t Number_of_Channels) {
    auto k = m.k();
    myMap map;
    std::vector<std::vector<std::size_t>> ns;
    std::vector<std::size_t> start(k, 0);
    std::vector<std::size_t> end(k, 0);
    std::size_t i = 0;
    start[0] = Number_of_Channels;
    end.back() = Number_of_Channels;
    auto v = start;
    ns.push_back(v);
    map[v] = i;
    ++i;
    while (Less{}(v, end)) {
      next(v);
      map[v] = i;
      ns.push_back(v);
      ++i;
    }

    return std::make_pair(ns, map);
  }

  static double chi2sum(const std::vector<std::size_t> &n,
                        const std::vector<double> &Nmean,
                        const M_Matrix<double> &pinvNCov) {
    return xTSigmaX(Nmean - n, pinvNCov);
  }

  static std::vector<std::size_t> to_n(const M_Matrix<double> &Pmean,
                                       std::size_t N) {
    std::vector<std::size_t> out(Pmean.size());
    std::multimap<double, std::size_t> remain;
    std::size_t sum = 0;
    for (std::size_t i = 0; i < out.size(); ++i) {
      double n = N * Pmean[i];
      out[i] = n;
      sum += out[i];
      remain.insert(std::make_pair(n - out[i], i));
    }
    auto it = remain.rbegin();
    while (sum < N) {
      auto i = it->second;
      ++out[i];
      ++sum;
      ++it;
    }
    return out;
  }

  static auto
  new_neighbours_up(const std::vector<std::size_t> &n,
                    std::vector<std::vector<std::size_t>> &connections_x,
                    const std::vector<double> &Nmean,
                    const M_Matrix<double> &pinvNCov, double maxChi2,
                    std::map<std::vector<std::size_t>, double> &accepted,
                    std::set<std::vector<std::size_t>> &rejected) {
    std::vector<std::vector<std::size_t>> out;
    std::map<std::vector<std::size_t>, double> new_accepted;
    std::set<std::vector<std::size_t>> new_rejected;
    for (std::size_t i_origin = 0; i_origin < n.size(); ++i_origin) {
      if (n[i_origin] > 0)
        for (auto &i_destination : connections_x[i_origin]) {
          auto destination = n;
          --destination[i_origin];
          ++destination[i_destination];
          if ((accepted.find(destination) != accepted.end()) &&
              (rejected.find(destination) != rejected.end())) {
            double logL = chi2sum(destination, Nmean, pinvNCov);
            if (logL < maxChi2) {
              new_accepted[destination] = logL;
              out.push_back(destination);
            } else
              new_rejected.insert(destination);
          }
        }
    }
    accepted.insert(new_accepted.begin(), new_accepted.end());
    rejected.insert(new_rejected.begin(), new_rejected.end());
    return out;
  }

  static auto search_neighbours_up(
      myMap &map_accepted, std::set<std::vector<std::size_t>> &rejected,
      const std::vector<std::size_t> &origin,
      const std::vector<std::vector<std::size_t>> &connections_x,
      const std::vector<double> &Nmean, const M_Matrix<double> &pinvNCov,
      double maxChi2) {
    std::vector<std::vector<std::size_t>> out;
    myMap new_accepted;
    std::size_t i = map_accepted.size();
    std::set<std::vector<std::size_t>> new_rejected;
    for (std::size_t i_origin = 0; i_origin < origin.size(); ++i_origin) {
      if (origin[i_origin] > 0)
        for (auto &i_destination : connections_x[i_origin]) {
          auto destination = origin;
          --destination[i_origin];
          ++destination[i_destination];
          if ((map_accepted.find(destination) == map_accepted.end()) &&
              (rejected.find(destination) == rejected.end())) {
            double logL = chi2sum(destination, Nmean, pinvNCov);
            if (logL < maxChi2) {
              new_accepted[destination] = i;
              ++i;
              out.push_back(destination);
            } else
              new_rejected.insert(destination);
          }
        }
    }
    map_accepted.insert(new_accepted.begin(), new_accepted.end());
    rejected.insert(new_rejected.begin(), new_rejected.end());
    return out;
  }

  static auto justMap(const Model m, double x, const M_Matrix<double> &Pmean,
                      const M_Matrix<double> &Pcov, std::size_t numChannels,
                      double maxChi2) {
    double N = numChannels;
    auto connections_x = m.connections_x(x);
    auto pinvNCov = pinv(Pcov * N).value();
    std::vector<double> Nmean = (Pmean * numChannels).toVector();
    auto v = to_n(Pmean, numChannels);
    auto accepted_up = justMap_up(connections_x, Nmean, pinvNCov, maxChi2);
    accepted_up[v] = 1;
    return accepted_up;
  }

  static auto
  justMap_up(const std::vector<std::size_t> &v,
             const std::vector<std::vector<std::size_t>> &connections_x,
             const M_Matrix<double> &Nmean, const M_Matrix<double> &pinvNCov,
             double maxChi2) {
    std::map<std::vector<std::size_t>, std::size_t> accepted;
    std::set<std::vector<std::size_t>> rejected;
    auto neighbours_list = new_neighbours_up(v, connections_x, Nmean, pinvNCov,
                                             maxChi2, accepted, rejected);
    while (!neighbours_list.empty()) {

      std::vector<std::vector<std::size_t>> newNeigh;
      for (auto &e : neighbours_list) {
        auto new_neighbour = new_neighbours_up(
            e, connections_x, Nmean, pinvNCov, maxChi2, accepted, rejected);
        newNeigh.insert(new_neighbour.begin(), new_neighbour.end());
      }
      neighbours_list = newNeigh;
    }
  }

  static void
  search_up(myMap &map, const std::vector<std::size_t> &origin,
            const std::vector<std::vector<std::size_t>> &connections_x,
            const M_Matrix<double> &Nmean, const M_Matrix<double> &pinvNCov,
            double maxChi2) {
    std::set<std::vector<std::size_t>> rejected;
    auto neighbours_list = search_neighbours_up(
        map, rejected, origin, connections_x, Nmean, pinvNCov, maxChi2);
    while (!neighbours_list.empty()) {

      std::vector<std::vector<std::size_t>> newNeigh;
      for (auto &e : neighbours_list) {
        auto new_neighbour = new_neighbours_up(map, rejected, e, connections_x,
                                               Nmean, pinvNCov, maxChi2);
        newNeigh.insert(new_neighbour.begin(), new_neighbour.end());
      }
      neighbours_list = newNeigh;
    }
  }

  static void
  add_to_map(myMap &map_accepted, std::size_t numChannels,
             const M_Matrix<double> &Pmean, const M_Matrix<double> &Pcov,
             const std::vector<std::vector<std::size_t>> &connections_x,
             double maxChi2) {
    auto pinvNCov = pinv(Pcov * numChannels).value();
    std::vector<double> Nmean = (Pmean * numChannels).toVector();
    auto origin = to_n(Pmean, numChannels);
    auto i = map_accepted.size();
    map_accepted[origin] = i;
    search_up(map_accepted, origin, connections_x, Nmean, pinvNCov, maxChi2);
  }

  static myMap
  Macro_to_Micro(std::size_t numChannels, const M_Matrix<double> &Pmean0,
                 const M_Matrix<double> &Pcov0, const M_Matrix<double> &Pmean1,
                 const M_Matrix<double> &Pcov1,
                 const std::vector<std::vector<std::size_t>> &connections_x,
                 double maxChi2, myMap map_accepted) {

    add_to_map(map_accepted, numChannels, Pmean0, Pcov0, connections_x,
               maxChi2);
    add_to_map(map_accepted, numChannels, Pmean1, Pcov1, connections_x,
               maxChi2);
    return map_accepted;
  }

  static auto Macro_to_Micro(const std::vector<std::vector<std::size_t>> &ns,
                             const M_Matrix<double> &Pmean,
                             const M_Matrix<double> &Pcov,
                             std::size_t numChannels) {
    auto pinvNCov = pinv(Pcov * numChannels).value();
    std::vector<double> Nmean = (Pmean * numChannels).toVector();

    M_Matrix<double> lik(1, ns.size());
    double sumlik = 0;
    for (std::size_t i = 0; i < ns.size(); ++i) {
      double l = std::exp(-0.5 * chi2sum(ns[i], Nmean, pinvNCov));
      lik[i] = l;
      sumlik += l;
    }
    M_Matrix<double> P = lik / sumlik;

    return P;
  }

  static auto Micro_to_Macro(const std::vector<std::vector<std::size_t>> &ns,
                             const M_Matrix<double> &P,
                             std::size_t numChannels) {
    assert(ns.size() == P.size());
    std::size_t numStates = ns[0].size();

    M_Matrix<double> Nmean(1, numStates, 0.0);
    M_Matrix<double> Nsqrmean(numStates, numStates, Matrix_TYPE::SYMMETRIC, 0);

    for (std::size_t i = 0; i < ns.size(); ++i) {
      M_Matrix<double> p = M_Matrix<std::size_t>(1ul, numStates, ns[i]);
      p *= P[i];
      Nmean += p;
      Nsqrmean += quadraticForm_XTX(p);
    }
    M_Matrix<double> Pmean = Nmean * (1.0 / numChannels);
    M_Matrix<double> Pcov =
        (Nsqrmean - quadraticForm_XTX(Nmean)) * (1.0 / numChannels);
    return std::make_tuple(Pmean, Pcov);
  }
  template <class E = std::vector<std::size_t>>
  static void Reduce(std::vector<E> &ns, M_Matrix<double> &P,
                     double reduce_by_p) {
    std::multimap<double, std::size_t> pmap;
    for (std::size_t i = 0; i < ns.size(); ++i) {
      pmap.emplace(P[i], i);
    }
    std::set<std::size_t> to_be_removed;
    double psum = 0;
    for (auto it = pmap.begin(); it != pmap.end(); ++it) {
      psum += it->first;
      if (psum < reduce_by_p) {
        to_be_removed.insert(it->second);
      } else
        break;
    }
    if (!to_be_removed.empty()) {
      auto newN = ns.size() - to_be_removed.size();
      std::vector<E> new_ns(newN);
      M_Matrix<double> newP(1, newN);
      std::size_t start = 0;
      std::size_t i_destination = 0;
      to_be_removed.insert(ns.size()); // hack to force to include the
                                       // elements past the last removed.
      for (auto e : to_be_removed) {
        for (std::size_t i_origin = start; i_origin < e; ++i_origin) {
          new_ns[i_destination] = std::move(ns[i_origin]);
          newP[i_destination] = P[i_origin];
          ++i_destination;
        }
        start = e + 1;
      }

      ns = std::move(new_ns);
      P = std::move(newP);
    }
  }

private:
  Model m_;
  std::vector<std::vector<std::size_t>> ns_;
  myMap map_;
};



class Markov_Transition_rate {

  M_Matrix<double> Qrun_; // transition rate matrix at time zero
  M_Matrix<double> g_;

  M_Matrix<double> V_;     // eigenvector of Qrun
  M_Matrix<double> W_;     // eigenvector of Qrun
  M_Matrix<double> landa_; // eigenvalues
  M_Matrix<double> error_landa_;
  M_Matrix<double> error_V_;

  M_Matrix<double> Wg_;
  M_Matrix<double> WgV_;

  static void clean_landa(M_Matrix<double> &la) {
    double maxla = std::min(max(la), 0.0);
    for (std::size_t i = 1; i < la.size(); ++i)
      if (la[i] >= maxla) {
        la[i] = 0;
      }
  }

  void init(Matrix_Decompositions::eigensystem_type &&eig) {
    std::tie(V_, landa_, W_, error_landa_, error_V_) = std::move(eig);
    clean_landa(landa_);
    Wg_ = W_ * g_;
    WgV_ = W_ * Matrix_Unary_Transformations::diag(g_) * V_;
    // assert(Frobenius_test::test(Qrun(),V()*landa()*W(),1,
    // std::sqrt(std::numeric_limits<double>::epsilon())));
  }

public:
  struct landa : public invariant {
    static Op_void test(const M_Matrix<double> &landa, double tol) {
      std::size_t num_zeros = 0;
      // bool has_zero_value=false;
      bool has_positive_value = false;
      double lamax = maxAbs(landa);
      if (!std::isfinite(lamax))
        return Op_void(false, " landa has a nonfinite maximum value: " +
                                  ToString(lamax));
      double tolerance = tol;
      std::stringstream ss;
      are_zero<false, double> is_zero(tolerance);
      are_non_positive<true, double> is_neg(tolerance);
      for (std::size_t i = 0; i < landa.size(); ++i) {
        if (is_zero.test(landa[i], ss))
          num_zeros++;
        else if (!is_neg.test(landa[i], ss)) {
          has_positive_value = true;
          ss << " at i=" << i << "\n";
        }
      }

      if ((num_zeros == 1) && !has_positive_value) {
        //   std::cerr<<"\naccepted landa!!\n"<<landa;
        return Op_void(true, "");
      } else if ((num_zeros != 1) && !has_positive_value)
        return Op_void(false,
                       " landa has " + ToString(num_zeros) + " zero values");
      else if ((num_zeros != 1) && has_positive_value)
        return Op_void(
            false, " landa has " + ToString(num_zeros) +
                       " zero values and has positive value(s): " + ss.str());
      else
        return Op_void(false, " landa has positive value(s): " + ss.str());
    }
  };

  static Op_void test(const Matrix_Decompositions::eigensystem_type &e,
                      double tolerance) {
    return landa::test(std::get<1>(e), tolerance);
  }

  static myOptional_t<Markov_Transition_rate>
  evaluate(M_Matrix<double> &&_Qrun, const M_Matrix<double> &_g,
           double tolerance) {
    typedef myOptional_t<Markov_Transition_rate> Op;
    auto eig = Matrix_Decompositions::EigenSystem_full_real_eigenvalue(_Qrun);
    if (eig) {
      if (auto eigtest = test(eig.value(), tolerance); eigtest.has_value())
        return Op(Markov_Transition_rate(std::move(_Qrun), _g,
                                         std::move(eig).value()));
      else
        return Op(false, "invalid eigenvalues " + eigtest.error());
    } else
      return Op(false, "eigenvalue decomposition fails:" + eig.error());
  }

  Markov_Transition_rate() = default;
  /// virtual copy constructors

  Markov_Transition_rate(M_Matrix<double> &&_Qrun, const M_Matrix<double> &_g,
                         Matrix_Decompositions::eigensystem_type &&eig)
      : Qrun_{std::move(_Qrun)}, g_{_g} {
    init(std::move(eig));
  }

  const M_Matrix<double> &Qrun() const {
    return Qrun_;
  } // transition rate matrix at time zero
  const M_Matrix<double> &V() const { return V_; } // eigenvector of Qrun
  const M_Matrix<double> &W() const { return W_; } // eigenvector of Qrun
  const M_Matrix<double> &landa() const { return landa_; } // eigenvalues
  auto &error_landa() const { return error_landa_; }
  auto &error_vector() const { return error_V_; }
  const M_Matrix<double> &g() const { return g_; }
  const M_Matrix<double> &Wg() const { return Wg_; }
  const M_Matrix<double> &WgV() const { return WgV_; }

  M_Matrix<double> calc_Peq() const {
    M_Matrix<double> p0(1, Qrun().nrows(), 1.0 / Qrun().nrows());
    M_Matrix<double> laexp(landa_.size(), landa_.size(), Matrix_TYPE::DIAGONAL);
    for (std::size_t i = 0; i < landa_.size(); ++i) {
      if (landa_(i, i) == 0.0)
        laexp(i, i) = 1.0;
      else
        laexp(i, i) = 0.0;
    }
    //  std::cerr<<"\n
    //  ----calc_Peq-----\np0="<<p0<<"laexp="<<laexp<<"V()"<<V()<<"W"<<W();
    return p0 * V() * laexp * W();
  }
};

template<template<class...>class Tr>
class Markov_Transition_rate_new_;
typedef Markov_Transition_rate_new_<C> Markov_Transition_rate_new;

template<template<class...>class Tr>
class Markov_Transition_rate_new_ {

    template<class T> using Tr_t=typename Tr<T>::type;
    Tr_t<M_Matrix<double>> Qrun_; // transition rate matrix at time zero
    Tr_t<M_Matrix<double>> g_;

    Tr_t<M_Matrix<double>> V_;     // eigenvector of Qrun
    Tr_t<M_Matrix<double>> W_;     // eigenvector of Qrun
    Tr_t<M_Matrix<double>> landa_; // eigenvalues
    M_Matrix<double> error_landa_;
    M_Matrix<double> error_V_;

    Tr_t<M_Matrix<double>> Wg_;
    Tr_t<M_Matrix<double>> WgV_;

    static void clean_landa(Tr_t<M_Matrix<double>> &la) {
        double maxla = std::min(max(center(la)), 0.0);
        for (std::size_t i = 1; i < la.size(); ++i)
            if (center(la)[i] >= maxla) {
                la.set(i,0.0);
            }
    }

    void init(Tr_t<Matrix_Decompositions::eigensystem_type> &&eig) {
        std::tie(V_, landa_, W_, error_landa_, error_V_) = std::move(eig);
        clean_landa(landa_);
        Wg_ = W_ * g_;
        using Matrix_Unary_Transformations::diag;
        WgV_ = W_ * diag(g_) * V_;
        // assert(Frobenius_test::test(Qrun(),V()*landa()*W(),1,
        // std::sqrt(std::numeric_limits<double>::epsilon())));
    }

public:
    struct landa : public invariant {
        static Op_void test(const M_Matrix<double> &landa, double tol) {
            std::size_t num_zeros = 0;
            // bool has_zero_value=false;
            bool has_positive_value = false;
            double lamax = maxAbs(landa);
            if (!std::isfinite(lamax))
                return Op_void(false, " landa has a nonfinite maximum value: " +
                                   ToString(lamax));
            double tolerance = tol;
            std::stringstream ss;
            are_zero<false, double> is_zero(tolerance);
            are_non_positive<true, double> is_neg(tolerance);
            for (std::size_t i = 0; i < landa.size(); ++i) {
                if (is_zero.test(landa[i], ss))
                    num_zeros++;
                else if (!is_neg.test(landa[i], ss)) {
                    has_positive_value = true;
                    ss << " at i=" << i << "\n";
                }
            }

            if ((num_zeros == 1) && !has_positive_value) {
                //   std::cerr<<"\naccepted landa!!\n"<<landa;
                return Op_void(true, "");
            } else if ((num_zeros != 1) && !has_positive_value)
                return Op_void(false,
                               " landa has " + ToString(num_zeros) + " zero values");
            else if ((num_zeros != 1) && has_positive_value)
                return Op_void(
                    false, " landa has " + ToString(num_zeros) +
                        " zero values and has positive value(s): " + ss.str());
            else
                return Op_void(false, " landa has positive value(s): " + ss.str());
        }
    };

    static Op_void test(const Tr_t<Matrix_Decompositions::eigensystem_type> &e,
                        double tolerance) {
        return landa::test(std::get<1>(e), tolerance);
    }

    static myOptional<Markov_Transition_rate_new_<Tr>, reg_tag>
    evaluate(Tr_t<M_Matrix<double>> &&_Qrun, const Tr_t<M_Matrix<double>> &_g,
             double tolerance) {
        typedef myOptional_t<Markov_Transition_rate> Op;
        auto eig = Matrix_Decompositions::EigenSystem_full_real_eigenvalue(_Qrun);
        if (eig) {
            if (auto eigtest = test(eig.value(), tolerance); eigtest.has_value())
                return Op(Markov_Transition_rate(std::move(_Qrun), _g,
                                                 std::move(eig).value()));
            else
                return Op(false, "invalid eigenvalues " + eigtest.error());
        } else
            return Op(false, "eigenvalue decomposition fails:" + eig.error());
    }

    Markov_Transition_rate_new_() = default;
    /// virtual copy constructors

    Markov_Transition_rate_new_(Tr_t<M_Matrix<double>> &&_Qrun, const Tr_t<M_Matrix<double>> &_g,
                           Tr_t<Matrix_Decompositions::eigensystem_type> &&eig)
        : Qrun_{std::move(_Qrun)}, g_{_g} {
        init(std::move(eig));
    }

    const Tr_t<M_Matrix<double>> &Qrun() const {
        return Qrun_;
    } // transition rate matrix at time zero
    const Tr_t<M_Matrix<double>> &V() const { return V_; } // eigenvector of Qrun
    const Tr_t<M_Matrix<double>> &W() const { return W_; } // eigenvector of Qrun
    const Tr_t<M_Matrix<double>> &landa() const { return landa_; } // eigenvalues
    auto &error_landa() const { return error_landa_; }
    auto &error_vector() const { return error_V_; }
    const Tr_t<M_Matrix<double>> &g() const { return g_; }
    const Tr_t<M_Matrix<double>> &Wg() const { return Wg_; }
    const Tr_t<M_Matrix<double>> &WgV() const { return WgV_; }

    Tr_t<M_Matrix<double>> calc_Peq() const
    {
        M_Matrix<double> p0(1, Qrun().nrows(), 1.0 / Qrun().nrows());
        M_Matrix<double> laexp(landa_.size(), landa_.size(), Matrix_TYPE::DIAGONAL);
        for (std::size_t i = 0; i < landa_.size(); ++i) {
            if (center(landa_)(i, i) == 0.0)
                laexp(i, i) = 1.0;
            else
                laexp(i, i) = 0.0;
        }
        //  std::cerr<<"\n
        //  ----calc_Peq-----\np0="<<p0<<"laexp="<<laexp<<"V()"<<V()<<"W"<<W();
        return p0 * V() * laexp * W();
    }
};






class Markov_Transition_step_single;
struct Markov_Transition_step_single_minimum {
    M_Matrix<double> PPn;
    M_Matrix<double> PGn;
    M_Matrix<double> PG_n;
    M_Matrix<double> PGG_n;
    std::size_t n;
    Markov_Transition_step_single_minimum(
                const Markov_Transition_step_single &x_);
    Markov_Transition_step_single_minimum &
            operator*=(const Markov_Transition_step_single &x_);
};

template <template<class...> class Tr> class Markov_Transition_step_mean_;

typedef Markov_Transition_step_mean_<C> Markov_Transition_step_mean;

template <template<class...> class Tr>struct Markov_Transition_step_mean_minimum_ {
    template<class T> using Tr_t=typename Tr<T>::type;

    Tr_t<M_Matrix<double>> PPn;
    Tr_t<M_Matrix<double>> PGn;
    Tr_t<M_Matrix<double>> PG_n;
    Tr_t<M_Matrix<double>> PGGn;
    std::size_t n;
    Markov_Transition_step_mean_minimum_(
        const Tr_t<Markov_Transition_step_mean> &x_);
    Markov_Transition_step_mean_minimum_ &
    operator+=(const Tr_t<Markov_Transition_step_mean> &x_);
};



template <template<class...> class Tr> class Markov_Transition_step_new_;

typedef Markov_Transition_step_new_<C> Markov_Transition_step_new;

template <template<class...> class Tr> class Markov_Transition_step_new_ {
public:
     template<class T> using Tr_t=typename Tr<T>::type;
protected:
    Tr_t<M_Matrix<double>> P_; // transition matrix
    Tr_t<M_Matrix<double>> g_; // conductance matrix
  std::size_t nsamples_;
  double dt_;
  double min_p_;

public:

  double min_P() const { return min_p_; }
  Tr_t<M_Matrix<double>> const &myP() const { return P_; }  // transition matrix
  Tr_t<M_Matrix<double>> const &P() const & { return P_; }  // transition matrix
  Tr_t<M_Matrix<double>> &&P() && { return std::move(P_); } // transition matrix

  double dt() const { return dt_; }
  double fs() const { return 1.0 / dt(); }

  std::size_t nsamples() const { return nsamples_; }

  /// mean conductance for each starting state i
  Tr_t<M_Matrix<double>> const &g() const { return g_; } // conductance matrix


  Markov_Transition_step_new_(Tr_t<M_Matrix<double>> &&P, const Tr_t<M_Matrix<double>> &g,
                         std::size_t nsamples, double dt, double min_p)
      : P_{std::move(P)},  g_{g}, nsamples_{nsamples},
        dt_{dt}, min_p_{min_p} {}

  Markov_Transition_step_new getMinimum() const { return *this; }

  Markov_Transition_step_new &operator*=(const Markov_Transition_step_new &other) {
    P_ = Probability_transition::normalize((P_ * other.P()), min_P());
    g_ = other.g();
    nsamples_ += other.nsamples();
    return *this;
  }

  Markov_Transition_step_new_() = default;

  typedef Markov_Transition_step_new_<Tr> self_type;
  constexpr static auto className = my_static_string("Markov_Transition_step");
  static auto get_constructor_fields() {
    // M_Matrix<double> const & (self_type::*myP) ()const& =&self_type::P;

    return std::make_tuple(
        grammar::field(C<self_type>{}, "P", &self_type::myP),
        grammar::field(C<self_type>{}, "g", &self_type::g),
        grammar::field(C<self_type>{}, "nsamples", &self_type::nsamples),
        grammar::field(C<self_type>{}, "dt", &self_type::dt),
        grammar::field(C<self_type>{}, "min_p", &self_type::min_P));
  }

 };


class Markov_Transition_step;

class Markov_Transition_step {
protected:
    M_Matrix<double> P_; // transition matrix
    M_Matrix<double> g_; // conductance matrix
    std::size_t nsamples_;
    double dt_;
    double min_p_;

public:
    double min_P() const { return min_p_; }
    M_Matrix<double> const &myP() const { return P_; }  // transition matrix
    M_Matrix<double> const &P() const & { return P_; }  // transition matrix
    M_Matrix<double> &&P() && { return std::move(P_); } // transition matrix

    double dt() const { return dt_; }
    double fs() const { return 1.0 / dt(); }

    std::size_t nsamples() const { return nsamples_; }

    /// mean conductance for each starting state i
    M_Matrix<double> const &g() const { return g_; } // conductance matrix

    Markov_Transition_step(const Markov_Transition_rate &Qx, std::size_t nsamples,
                           double fs, double min_p)
        : P_{}, /*ladt_{},exp_ladt_{},*/ g_{Qx.g()}, nsamples_{nsamples},
          dt_{1.0 / fs * nsamples}, min_p_{min_p} {
        init(Qx);
    }

    Markov_Transition_step(M_Matrix<double> &&P, const M_Matrix<double> &g,
                           std::size_t nsamples, double fs, double min_p)
        : P_{std::move(P)}, /*ladt_{},exp_ladt_{},*/ g_{g}, nsamples_{nsamples},
          dt_{1.0 / fs * nsamples}, min_p_{min_p} {}

    Markov_Transition_step(std::size_t n_sub_samples,
                           const Markov_Transition_rate &Qx, std::size_t nsamples,
                           double fs, double min_p)
        : P_{}, /*ladt_{},exp_ladt_{},*/ g_{Qx.g()}, nsamples_{nsamples},
          dt_{1.0 / (fs * n_sub_samples) * nsamples}, min_p_{min_p} {
        init(Qx);
    }

    Markov_Transition_step(M_Matrix<double> &&P, M_Matrix<double> &&g,
                           std::size_t nsamples, double fs, double min_p)
        : P_{std::move(P)},
          //          P_{Probability_transition::normalize(std::move(P),min_p)},
          g_{std::move(g)}, nsamples_{nsamples}, dt_{(1.0 * nsamples) / fs},
          min_p_{min_p} {}

    Markov_Transition_step getMinimum() const { return *this; }

    Markov_Transition_step &operator*=(const Markov_Transition_step &other) {
        P_ = Probability_transition::normalize((P_ * other.P()), min_P());
        g_ = other.g();
        nsamples_ += other.nsamples();
        return *this;
    }

    Markov_Transition_step() = default;

    typedef Markov_Transition_step self_type;
    constexpr static auto className = my_static_string("Markov_Transition_step");
    static auto get_constructor_fields() {
        // M_Matrix<double> const & (self_type::*myP) ()const& =&self_type::P;

        return std::make_tuple(
            grammar::field(C<self_type>{}, "P", &self_type::myP),
            grammar::field(C<self_type>{}, "g", &self_type::g),
            grammar::field(C<self_type>{}, "y_mean", &self_type::nsamples),
            grammar::field(C<self_type>{}, "min_p", &self_type::min_P));
    }

private:
    void init(const Markov_Transition_rate &Qx) {
        auto ladt = Qx.landa() * dt();
        auto exp_ladt = ladt.apply([](double x) { return std::exp(x); });
        P_ = Probability_transition::normalize(Qx.V() * exp_ladt * Qx.W(), min_P());
        g_ = Qx.g();
    }
};



inline M_Matrix<double> calc_P_taylor(const M_Matrix<double> &Qrun,
                                      std::size_t nsamples, double fs) {
  auto x = Qrun * (1.0 / fs * nsamples);
  return Matrix_Unary_Transformations::expm_taylor(x);
}

template <template<class...> class Tr>
class Markov_Transition_step_mean_ : public Markov_Transition_step_new_<Tr> {
public:
    template<class T> using Tr_t=typename Tr<T>::type;
private:
    Tr_t<M_Matrix<double>> gmean_i_; // conductance matrix


    Tr_t<M_Matrix<double>> gsqr_i_;

    Tr_t<M_Matrix<double>> gvar_i_;

    Tr_t<M_Matrix<double>> gmean_ij_;

    Tr_t<M_Matrix<double>> gtotal_ij_; // conductance matrix
    // M_Matrix<double> ladt_; // exp(la dt)

public:
    typedef Tr_t<Markov_Transition_step_new> base_type;

    /// total conductance for each starting state i and ending state j
    Tr_t<M_Matrix<double>> const &gtotal_ij() const {
        return gtotal_ij_;
    } // conductance matrix

    /// squared mean conductance for each starting state i
    Tr_t<M_Matrix<double>> const &gsqr_i() const {
        return gsqr_i_;
    } // conductance matrix

    std::size_t nsamples() const { return base_type::nsamples_; }

    /// mean conductance for each starting state i
    Tr_t<M_Matrix<double>> const &gmean_i() const {
        return gmean_i_;
    } // conductance matrix

    Tr_t<M_Matrix<double>> const &gmean_ij() const {
        return gmean_ij_;
    } // conductance matrix

    /// variance of the mean conductance for each starting state i
    Tr_t<M_Matrix<double>> const &gvar_i() const {
        return gvar_i_;
    }


        Markov_Transition_step_mean_(Tr_t<M_Matrix<double>> &&P, Tr_t<M_Matrix<double>> &&g,
                                      std::size_t nsamples, double fs, double min_p,
                                      Tr_t<M_Matrix<double>> &&gmean_i,
                                      Tr_t<M_Matrix<double>> &&gsqr_i,
                                      Tr_t<M_Matrix<double>> &&gvar_i,
                                      Tr_t<M_Matrix<double>> &&gtotal_ij,
                                      Tr_t<M_Matrix<double>> &&gmean_ij
                                      )
        : base_type(std::move(P), std::move(g), nsamples, nsamples/fs, min_p),
           gmean_i_{std::move(gmean_i)},  gsqr_i_{std::move(gsqr_i)}, gvar_i_{std::move(gvar_i)},gmean_ij_{std::move(gmean_ij)},gtotal_ij_{std::move(gtotal_ij)}{}

        Markov_Transition_step_mean_(base_type&& base,
                                    Tr_t<M_Matrix<double>> &&gmean_i,
                                    Tr_t<M_Matrix<double>> &&gsqr_i,
                                    Tr_t<M_Matrix<double>> &&gvar_i,
                                    Tr_t<M_Matrix<double>> &&gtotal_ij,
                                    Tr_t<M_Matrix<double>> &&gmean_ij
                                    )
            : base_type(std::move(base)),
              gmean_i_{gmean_i},  gsqr_i_{gsqr_i}, gvar_i_{gvar_i}, gmean_ij_{gmean_ij},gtotal_ij_{gtotal_ij}{}


    Markov_Transition_step_mean_() = default;

};



class Markov_Transition_step_single {
protected:
  M_Matrix<double> P_; // transition matrix

  M_Matrix<double> gmean_i_; // conductance matrix

  M_Matrix<double> gtotal_ij_; // conductance matrix

  M_Matrix<double> gsqr_i_;
  // M_Matrix<double> ladt_; // exp(la dt)

  //  M_Matrix<double> exp_ladt_; // exp(la dt)

  M_Matrix<double> g_; // conductance matrix

  std::size_t nsamples_;
  double dt_;

public:
  Markov_Transition_step_single_minimum getMinimum() const {
    return Markov_Transition_step_single_minimum(*this);
  }
  M_Matrix<double> const &P() const { return P_; } // transition matrix

  /// total conductance for each starting state i and ending state j
  M_Matrix<double> const &gtotal_ij() const {
    return gtotal_ij_;
  } // conductance matrix


  /// squared mean conductance for each starting state i
  M_Matrix<double> const &gsqr_i() const {
    return gsqr_i_;
  } // conductance matrix

  double dt() const { return dt_; }

  // M_Matrix<double> const & ladt(){return ladt_;}; // exp(la dt)

  //  M_Matrix<double> const& exp_ladt(){return exp_ladt_;} // exp(la dt)

  std::size_t nsamples() const { return nsamples_; }

  /// mean conductance for each starting state i
  M_Matrix<double> const &gmean_i() const {
    return gmean_i_;
  } // conductance matrix

  Markov_Transition_step_single(const Markov_Transition_rate &Qx,
                                std::size_t nsamples, double fs)
      : P_{}, /*ladt_{},exp_ladt_{},*/ g_{Qx.g()}, nsamples_{nsamples},
        dt_{1.0 / fs * nsamples} {
    init(Qx);
  }

  Markov_Transition_step_single(M_Matrix<double> &&P, M_Matrix<double> &&g,
                                std::size_t nsamples, double fs)
      : P_{std::move(P)}, g_{std::move(g)}, nsamples_{nsamples},
        dt_{(1.0 * nsamples) / fs} {}

  Markov_Transition_step_single() = default;

private:
  void init(const Markov_Transition_rate &Qx) {

    auto ladt = Qx.landa() * dt();
    auto exp_ladt = exp(ladt);
    P_ = Qx.V() * exp_ladt * Qx.W();
    g_ = Qx.g();
  }
};

inline Markov_Transition_step_single_minimum::
    Markov_Transition_step_single_minimum(
        const Markov_Transition_step_single &x) {
  n = x.nsamples();
  PPn = x.P();
  PG_n = x.gmean_i() * n;
  PG_n = x.gtotal_ij() * n;
  PGG_n = x.gsqr_i() * (n * n);
}

template <template<class...> class Tr>
inline Markov_Transition_step_mean_minimum_<Tr>::
    Markov_Transition_step_mean_minimum_(const Tr_t<Markov_Transition_step_mean> &x) {
    n = x.nsamples();
    PPn = x.P();
    PG_n = x.gtotal_ij() * n;
}



template <typename E>
auto calc_Gtot_taylor(const M_Matrix<E> &Qx, const M_Matrix<E> &P,
                      const M_Matrix<E> &G, std::size_t order) {
  assert(G.isDiagonal());
  assert(G.ncols() == Qx.ncols());
  assert(G.ncols() == P.ncols());
  assert(Qx.nrows() == P.nrows());

  M_Matrix<E> dGn = G * P;
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

template <typename E>
M_Matrix<E> calc_Gvar_taylor(const M_Matrix<E> &Qx, const M_Matrix<E> &P,
                             const M_Matrix<E> &G, std::size_t order) {
  assert(G.isDiagonal());
  assert(G.ncols() == Qx.ncols());
  assert(G.ncols() == P.ncols());
  assert(Qx.nrows() == P.nrows());

  M_Matrix<E> Tau_n = G;
  std::vector<double> coeff(1, 1.0);
  auto GP = G * P;
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

inline Markov_Transition_step_single_minimum &
Markov_Transition_step_single_minimum::
operator*=(const Markov_Transition_step_single &x) {
  auto n1 = x.nsamples();
  PGn += ((PPn * x.gmean_i()) * n1);
  PGG_n += ((PPn * x.gsqr_i()) * (n1 * n1) + (PG_n * x.gmean_i()) * n1);
  PG_n = (PG_n * x.P()) + (PPn * x.gtotal_ij()) * n1;
  PPn = PPn * x.P();
  n += n1;
  return *this;
}

template <template<class...> class Tr>
inline Markov_Transition_step_mean_minimum_<Tr> &
Markov_Transition_step_mean_minimum_<Tr>::
operator+=(const Tr_t<Markov_Transition_step_mean> &x) {
    auto n1 = x.nsamples();
    PGn += ((PPn * x.gmean_i()) * n1);
    PG_n = (PG_n * x.P()) + (PPn * x.gtotal_ij()) * n1;
    PPn = PPn * x.P();
    n += n1;
    return *this;
}


template <template<class...> class Tr> class Markov_Transition_step_variance_;

typedef Markov_Transition_step_variance_<C> Markov_Transition_step_variance;

template <template<class...> class Tr>
struct Markov_Transition_step_variance_minimum_;
typedef Markov_Transition_step_variance_minimum_<C> Markov_Transition_step_variance_minimum;

template <template<class...> class Tr>
struct Markov_Transition_step_variance_minimum_ {
    template<class T> using Tr_t=typename Tr<T>::type;

    std::size_t n;
    Tr_t<M_Matrix<double>> PPn;
    Tr_t<M_Matrix<double>> PG_n;
    Tr_t<M_Matrix<double>> PGG_n;
    double min_p;
    double tolerance;
    Markov_Transition_step_variance_minimum_(
        const Tr_t<Markov_Transition_step_variance> &x_);

    Markov_Transition_step_variance_minimum_(const Tr_t<M_Matrix<double>> &P,
                                          const Tr_t<M_Matrix<double>> &PG,
                                          const Tr_t<M_Matrix<double>> &PGG,
                                          std::size_t n, double min_p__,
                                          double tolerance__)
        : PPn{P}, PG_n{PG * n}, PGG_n{PGG * n}, min_p{min_p__},
          tolerance{tolerance__} {}

    Markov_Transition_step_variance_minimum_ &
    operator+=(const Tr_t<Markov_Transition_step_variance> &x_);
};



class Markov_Transition_step_double;

struct Markov_Transition_step_double_minimum {
  std::size_t n;
  M_Matrix<double> PPn;
  M_Matrix<double> PG_n;
  M_Matrix<double> PGG_n;
  double min_p;
  double tolerance;
  Markov_Transition_step_double_minimum(
      const Markov_Transition_step_double &x_);

  Markov_Transition_step_double_minimum(const M_Matrix<double> &P,
                                        const M_Matrix<double> &PG,
                                        const M_Matrix<double> &PGG,
                                        std::size_t n, double min_p__,
                                        double tolerance__)
      : PPn{P}, PG_n{PG * n}, PGG_n{PGG * n}, min_p{min_p__},
        tolerance{tolerance__} {}

  Markov_Transition_step_double_minimum &
  operator*=(const Markov_Transition_step_double &x_);
};





template <bool output> struct are_Equal<output, Markov_Transition_step_double> {
public:
  template <class ostream>
  bool test(const Markov_Transition_step_double &one,
            const Markov_Transition_step_double &other, double eps,
            ostream &os = std::cerr) const;
};

template <template<class...> class Tr> class Markov_Transition_step_variance_;
typedef Markov_Transition_step_variance_<C> Markov_Transition_step_variance;

template <template<class...> class Tr> class Markov_Transition_step_variance_ : public Markov_Transition_step_mean_<Tr> {
public:
    template<class T> using Tr_t=typename Tr<T>::type;
    typedef Markov_Transition_step_mean_<Tr> base_type;
private:
    Tr_t<M_Matrix<double>> gtotal_sqr_ij_; // conductance matrix
    Tr_t<M_Matrix<double>> gtotal_var_ij_; // variance of the conductance matrix
    Tr_t<M_Matrix<double>> gvar_ij_; // variance of the conductance matrix

public:

    Markov_Transition_step_variance_(Tr_t<M_Matrix<double>> &&P, Tr_t<M_Matrix<double>> &&g,
                                  std::size_t nsamples, double fs, double min_p,
                                    Tr_t<M_Matrix<double>> &&gmean_i,
                                  Tr_t<M_Matrix<double>> &&gsqr_i,
                                  Tr_t<M_Matrix<double>> &&gvar_i,
                                  Tr_t<M_Matrix<double>> &&gtotal_ij,
                                  Tr_t<M_Matrix<double>> &&gmean_ij,
                                  Tr_t<M_Matrix<double>> &&gtotal_sqr_ij,
                                  Tr_t<M_Matrix<double>> &&gtotal_var_ij,
                                  Tr_t<M_Matrix<double>> &&gvar_ij)
        : Markov_Transition_step_mean_<Tr>(std::move(P), std::move(g), nsamples, fs,
                                      min_p,
                                      std::move(gmean_i),std::move(gsqr_i),std::move(gvar_i),std::move(gtotal_ij),std::move(gmean_ij)),
          gtotal_sqr_ij_{std::move(gtotal_sqr_ij)},
          gtotal_var_ij_{std::move(gtotal_var_ij)}, gvar_ij_{std::move(gvar_ij)} {}


    Markov_Transition_step_variance_(Tr_t<Markov_Transition_step_new>&& P,
                                    Tr_t<M_Matrix<double>> &&gmean_i,
                                    Tr_t<M_Matrix<double>> &&gsqr_i,
                                    Tr_t<M_Matrix<double>> &&gvar_i,
                                    Tr_t<M_Matrix<double>> &&gtotal_ij,
                                    Tr_t<M_Matrix<double>> &&gmean_ij,
                                    Tr_t<M_Matrix<double>> &&gtotal_sqr_ij,
                                    Tr_t<M_Matrix<double>> &&gtotal_var_ij,
                                    Tr_t<M_Matrix<double>> &&gvar_ij)
        : Tr_t<Markov_Transition_step_mean>(std::move(P),std::move(gmean_i),std::move(gsqr_i),std::move(gvar_i),std::move(gtotal_ij),std::move(gmean_ij)),
          gtotal_sqr_ij_{std::move(gtotal_sqr_ij)},
          gtotal_var_ij_{std::move(gtotal_var_ij)}, gvar_ij_{std::move(gvar_ij)} {}


    Markov_Transition_step_variance_(Tr_t<Markov_Transition_step_mean>&& m,
                                        Tr_t<M_Matrix<double>> &&gtotal_sqr_ij,
                                        Tr_t<M_Matrix<double>> &&gtotal_var_ij,
                                        Tr_t<M_Matrix<double>> &&gvar_ij)
        : Tr_t<Markov_Transition_step_mean>(std::move(m)),
          gtotal_sqr_ij_{std::move(gtotal_sqr_ij)},
          gtotal_var_ij_{std::move(gtotal_var_ij)}, gvar_ij_{std::move(gvar_ij)} {}


    Tr_t<Markov_Transition_step_variance_minimum> getMinimum() const {
        return Tr_t<Markov_Transition_step_variance_minimum>(*this);
    }


    /// squared mean conductance for each starting state i and ending state j
    auto &gtotal_sqr_ij() const {
        return gtotal_sqr_ij_;
    } // conductance matrix

    /// variance of the mean conductance for each starting state i contributed
    /// by the ones ending at state j
    auto const &gtotal_var_ij() const {
        return gtotal_var_ij_;
    } // variance of the conductance matrix

    /// variance of the mean conductance for each starting state i and ending
    /// state j
    Tr_t<M_Matrix<double>> const &gvar_ij() const {
        return gvar_ij_;
    } // variance of the conductance matrix



    Markov_Transition_step_variance_() = default;

    typedef Markov_Transition_step_variance_ self_type;
    constexpr static auto className =my_template_trait<Tr>::className+
        my_static_string("Markov_Transition_step_variance");
    static auto get_constructor_fields() {
        // M_Matrix<double> const & (self_type::*myP) ()const& =&self_type::P;
        /*
 *   Markov_Transition_step_double(M_Matrix<double> &&P,
 * M_Matrix<double> &&g,
 * std::size_t nsamples,
 * double fs,
 *  double min_p,
 * M_Matrix<double>
     &&gmean_i, M_Matrix<double> &&gtotal_ij, M_Matrix<double> &&gmean_ij,

                                    M_Matrix<double> &&gtotal_sqr_ij,

                                    M_Matrix<double> &&gsqr_i,

                                    M_Matrix<double> &&gvar_i,
                                    M_Matrix<double> &&gtotal_var_ij,
                                    M_Matrix<double> &&gvar_ij)

 * */
        return std::tuple_cat(
            base_type::get_constructor_fields(),
            std::make_tuple(
            grammar::field(C<self_type>{}, "gtotal_sqr_ij",
                           &self_type::gtotal_sqr_ij),
            grammar::field(C<self_type>{}, "gtotal_var_ij",
                           &self_type::gtotal_var_ij),
            grammar::field(C<self_type>{}, "gvar_ij", &self_type::gvar_ij))
                );
    }


};
template <template<class...> class Tr>
class Markov_Transition_step_calculations_;

typedef Markov_Transition_step_calculations_<C> Markov_Transition_step_calculations;
template <template<class...> class Tr>
class Markov_Transition_step_calculations_ {

public:
    template<typename T> using Tr_t=typename Tr<T>::type;

    Tr_t<Markov_Transition_step_new> calc_step(const Tr_t<Markov_Transition_rate_new> &Qx,
                                                       std::size_t nsamples, double dt, double min_p) const  {

            auto ladt = Qx.landa() * dt;
            auto exp_ladt = ladt.apply([](double x) { return std::exp(x); });
            auto P = Probability_transition::normalize(Qx.V() * exp_ladt * Qx.W(), min_p);
            auto g = Qx.g();
            return Tr_t<Markov_Transition_step_new>(std::move(P),std::move(g),nsamples,dt,min_p);

        }


        myOptional_t<Tr_t<Markov_Transition_step_variance>> calc_step_variance(const Tr_t<Markov_Transition_rate_new> &Qx,
                                                           std::size_t nsamples, double fs, double min_p,
                                double max_landa_error) const  {
    if (Qx.error_landa() > max_landa_error)
      return  calc_step_variance(Qx, nsamples, fs,min_p);
   else
    return calc_step_variance(Qx, nsamples, fs, min_p, true);
  }


  myOptional_t<Tr_t<Markov_Transition_step_variance>> calc_step_variance(const Tr_t<Markov_Transition_rate_new> &Qx,
                                                     std::size_t nsamples,double fs, double min_P)const
{
      using std::exp;
      typedef myOptional_t<Tr_t<Markov_Transition_step_variance>> Op;
       std::size_t N=Qx.Qrun().ncols();
       auto ladt = Qx.landa() * (1.0*nsamples/fs);
       auto exp_ladt = ladt.apply([](Tr_t<double> x) { return exp(x); });
       auto P = Tr_t<Probability_transition_new>::normalize(Qx.V() * exp_ladt * Qx.W(), min_P);
       auto g = Qx.g();

           double C_V = max(Qx.error_vector());
           double C_l = max(Qx.error_landa());
           double err = C_V * C_V * C_l * std::numeric_limits<double>::epsilon();

           assert(are_Equal_v(Qx.Qrun(), Qx.V() * Qx.landa() * Qx.W(), err, err,
                              std::cerr, "Qx.Qrun()=", Qx.Qrun(), "Qx.V()", Qx.V(),
                              "Qx.landa()=", Qx.landa(), "Qx.W()=", Qx.W()));
           [[maybe_unused]] double RV = Matrix_Unary_Functions::norm_1(center(Qx.V())) * Matrix_Unary_Functions::norm_1(center(Qx.W()));


           Tr_t<M_Matrix<double>> E2m(N, N, Matrix_TYPE::SYMMETRIC);
           // M_Matrix<double> E2mb(N, N, Matrix_TYPE::SYMMETRIC);
           for (std::size_t i = 0; i < N; ++i)
               for (std::size_t j = 0; j < i + 1; ++j)
                   E2m.set(i, j,
                           Ee(ladt[i], ladt[j], exp_ladt[i], exp_ladt[j] /*, sqr(min_P())*/));

           // build E2
           Tr_t<M_Matrix<double>> WgV_E2(N, N,Matrix_TYPE::FULL);

           for (std::size_t i = 0; i < N; ++i)
               for (std::size_t j = 0; j < N; ++j)
                   WgV_E2.set(i, j,Qx.WgV()(i, j) * E2m(i, j));
           auto G = diag(g);
           auto gtotal_ij_p = (G * P + P * G) * 0.5;
           using std::abs;
           auto gtotal_var_ij_p =
               (G * P - P * G).apply([](auto x) { return abs(x); }) * 0.5;
           auto gmean_ij_p = elemDivSafe(gtotal_ij_p, P);
           auto gtotal_sqr_ij_p = gtotal_var_ij_p + elemMult(gtotal_ij_p, gmean_ij_p);

           auto gtotal_ij = (Qx.V() * WgV_E2 * Qx.W());
           auto gmean_ij = elemDivSafe(gtotal_ij, P);
           Tr_t<double> zero(0);

           Tr_t<M_Matrix<double>> WgV_E3(N, N,Matrix_TYPE::FULL, 0.0);
           for (std::size_t n1 = 0; n1 < N; n1++)
               for (std::size_t n3 = 0; n3 < N; n3++)
                   for (std::size_t n2 = 0; n2 < N; n2++) {
                       Tr_t<double> e3 = E3(ladt[n1], ladt[n2], ladt[n3], exp_ladt[n1],
                                      exp_ladt[n2], exp_ladt[n3] /*, sqr(min_P())*/);
                       WgV_E3.add(n1, n3, Qx.WgV()(n1, n2) * Qx.WgV()(n2, n3) * e3);
                       // optimizable
                   }


           double minP = N * N * 20 * 1e-6;
           // std::pow(R_V, 3) * C_l *std::numeric_limits<double>::epsilon();
           auto gtotal_sqr_ij = Qx.V() * WgV_E3 * Qx.W() * 2.0;
           for (std::size_t i = 0; i < gtotal_sqr_ij.size(); ++i) {
               if (std::abs(center(gtotal_ij)[i]) < minP)
                   gtotal_ij.set(i, zero);
               if (std::abs(center(gtotal_sqr_ij)[i]) < minP)
                   gtotal_sqr_ij.set(i, 0);
           }
           auto gtotal_var_ij = gtotal_sqr_ij - elemMult(gtotal_ij, gmean_ij);
           for (std::size_t i = 0; i < gtotal_sqr_ij.size(); ++i) {
               if (std::abs(center(gtotal_sqr_ij)[i]) < minP)
                   gtotal_var_ij.set(i, 0);
           }

           auto gvar_ij = elemDivSafe(gtotal_var_ij, P);
           for (std::size_t i = 0; i < gvar_ij.size(); ++i)
               if (std::abs(center(gvar_ij)[i]) < minP)
                   gvar_ij.set(i, 0);

           M_Matrix<double> u(N, 1, 1.0);
           auto gmean_i = gtotal_ij * u;
           auto gsqr_i = gtotal_sqr_ij * u;
           auto gvar_i = gtotal_var_ij * u;
           auto r = elemDivSafe(gvar_ij, gmean_ij);
           if (false) {
               //      auto out2 =
               //          Markov_Transition_step_double(4, Qx, nsamples(), fs(),
               //          1e-7);
               assert(([&gvar_ij]() {
                   for (std::size_t i = 0; i < gvar_ij.size(); ++i)
                       if (center(gvar_ij)[i] < 0)
                           return false;
                   return true;
               }()));
               //      assert(([this, &Qx, &check, &out2]() {
               //        return are_Equal<true, Markov_Transition_step_double>().test(
               //            *this, out2, 1e-7, std::cerr);
               //      }()));
           }
           return Op( Tr_t<Markov_Transition_step_variance>(std::move(P),std::move(g),nsamples,fs,min_P,std::move(gmean_i),std::move(gsqr_i),std::move(gvar_i),std::move(gtotal_ij),std::move(gmean_ij),std::move(gtotal_sqr_ij),std::move(gtotal_var_ij),std::move(gvar_ij)));
       }


       inline Tr_t<M_Matrix<double>> calc_P_taylor(const Tr_t<M_Matrix<double>> &Qrun,
                                             std::size_t nsamples, double fs) {
           auto x = Qrun * (1.0 / fs * nsamples);
           return Matrix_Unary_Transformations::expm_taylor(x);
       }


       Tr_t<Markov_Transition_step_variance> calc_step_variance(Tr_t<Markov_Transition_step_new> &&s)
       {

           M_Matrix<double> G = diag(s.g());
           std::size_t N = s.P().ncols();
           auto gtotal_ij = (G * s.P() + s.P() * G) * 0.5;
           auto gtotal_var_ij =
               (G * s.P() - s.P() * G).apply([](double x) { return std::abs(x); }) * 0.5;
           auto gmean_ij = elemDivSafe(gtotal_ij, s.P());
           auto gtotal_sqr_ij = gtotal_var_ij + elemMult(gtotal_ij, gmean_ij);
           auto gvar_ij = elemDivSafe(gtotal_var_ij, s.P());
           M_Matrix<double> u(N, 1, 1.0);
           auto gmean_i = gtotal_ij * u;
           auto gsqr_i = gtotal_sqr_ij * u;
           auto gvar_i = gtotal_var_ij * u;
           return Markov_Transition_step_variance(std::move(s),
                                                  std::move(gmean_i),std::move(gsqr_i),std::move(gvar_i),std::move(gtotal_ij),
                                                  std::move(gmean_ij),std::move(gtotal_sqr_ij),std::move(gtotal_var_ij),std::move(gvar_ij));
       }


       Tr_t<Markov_Transition_step_variance> calc_step_variance(
           const Tr_t<M_Matrix<double>> &Qrun,
           const Tr_t<M_Matrix<double>> &grun,
                                std::size_t nsamples, double fs, double min_p,std::size_t recursion_number_extra) {
      int recursion_number =
          2 + std::log2(maxAbs(Qrun) * 1.0 / fs * nsamples) + recursion_number_extra;
      recursion_number = std::max(4, recursion_number);
      std::size_t times = std::pow(2, recursion_number);
      Tr_t<Markov_Transition_step_variance> s=calc_step_variance(Markov_Transition_step_new(
          calc_P_taylor(Qrun, nsamples, fs * times), grun, nsamples, fs* times, min_p));
      Tr_t<Markov_Transition_step_variance> step_ini(s);
      for (int i = 0; i < recursion_number; ++i) {
          Tr_t<Markov_Transition_step_variance_minimum> step(step_ini);
          step += step_ini;
          step_ini = calc_step_variance(std::move(step), step_ini.g(),
                                                   fs * times, min_p);
      }
      return step_ini;

     }




     Tr_t<Markov_Transition_step_variance> calc_step_variance(Tr_t<Markov_Transition_step_variance_minimum> &&x,
                                                            Tr_t<M_Matrix<double>> g, double fs, double min_P)
   //   : Markov_Transition_step(std::move(x.PPn), std::move(g), x.n, fs, min_p) {
  {
      M_Matrix<double> u = Matrix_Generators::ones<double>(x.PPn.ncols(), 1ul);

      auto gtotal_sqr_ij = x.PGG_n / x.n / x.n * 2;
      auto gtotal_ij = x.PG_n / x.n;

      auto gmean_ij= elemDivSafe(gtotal_ij, x.PPn, min_P);

      auto gtotal_var_ij = gtotal_sqr_ij - elemMult(gtotal_ij, gmean_ij);


      auto gvar_ij = elemDivSafe(gtotal_var_ij, x.PPn, min_P);
      auto gmean_i = gtotal_ij * u;
      auto gsqr_i = gtotal_sqr_ij * u;
      auto gvar_i = gtotal_var_ij * u;

      return Tr_t<Markov_Transition_step_variance>(std::move(x.PPn),std::move(g),x.n,fs,min_P,
                                             std::move(gmean_i),std::move(gsqr_i),std::move(gvar_i),
                                             std::move(gtotal_ij),std::move(gmean_ij), std::move(gtotal_sqr_ij),
                                             std::move(gtotal_var_ij), std::move(gvar_ij));
  }

  Markov_Transition_step_calculations_() = default;

  typedef Markov_Transition_step_calculations_<Tr> self_type;
  constexpr static auto className =
      my_static_string("Markov_Transition_step_calculations");
  static auto get_constructor_fields() {
    // M_Matrix<double> const & (self_type::*myP) ()const& =&self_type::P;
    /*
 *   Markov_Transition_step_double(M_Matrix<double> &&P,
 * M_Matrix<double> &&g,
 * std::size_t nsamples,
 * double fs,
 *  double min_p,
 * M_Matrix<double>
     &&gmean_i, M_Matrix<double> &&gtotal_ij, M_Matrix<double> &&gmean_ij,

                                    M_Matrix<double> &&gtotal_sqr_ij,

                                    M_Matrix<double> &&gsqr_i,

                                    M_Matrix<double> &&gvar_i,
                                    M_Matrix<double> &&gtotal_var_ij,
                                    M_Matrix<double> &&gvar_ij)

 * */
  }

  static double E1(double x) {
    if (std::abs(x) < std::numeric_limits<double>::epsilon() * 100)
      return 1.0;
    else if (std::abs(x) < 1e-2)
      return std::expm1(x) / x;
    else
      return (std::exp(x) - 1.0) / x;
  }

  static double E2(double x, double y) {
    const double eps = std::numeric_limits<double>::epsilon();
    if (x * x < eps) {
      if (y * y < eps)
        return 0.5;
      else
        return (E1(y) - 1.0) / y;
    } else if (y * y < eps)
      return (E1(x) - 1.0) / x;
    else if ((y - x) * (y - x) < eps)
      return (std::exp(x) - E1(x)) / x;
    else
      return (E1(y) - E1(x)) / (y - x);
  }

  static Tr_t<double> Ee(Tr_t<double> x, Tr_t<double> y, Tr_t<double> exp_x, Tr_t<double> exp_y,
                   double eps = std::numeric_limits<double>::epsilon()) {
      if (sqr(center(x) - center(y)) < eps)
      return exp_x;
    else
      return (exp_x - exp_y) / (x - y);
  }

  static Tr_t<double> EX_111(Tr_t<double> x, Tr_t<double> y, Tr_t<double> z, Tr_t<double> exp_x) {
    return exp_x / ((x - y) * (x - z));
  }
  static Tr_t<double> E12_exmp1(Tr_t<double> x, Tr_t<double> y, Tr_t<double> z, Tr_t<double> exp_x,
                          Tr_t<double> exp_y) {
      using std::expm1;
    Tr_t<double> dz = z - x;
    return (dz * (exp_x - exp_y) + (y - x) * exp_x * expm1(dz)) /
           ((x - y) * (y - z) * dz);
  }
  static Tr_t<double> E111(Tr_t<double> x, Tr_t<double> y, Tr_t<double> z, Tr_t<double> exp_x, Tr_t<double> exp_y,
                     Tr_t<double> exp_z, double min2) {
      if (sqr(center(x) - center(y)) < min2) // x==y
    {
          if (sqr(center(y) - center(z)) < min2) // y==z
        return E3(x, y, z);  // x==y==z
      else
        return E12_exmp1(z, x, y, exp_z, exp_x); // x==y!=z
      } else if (sqr(center(y) - center(z)) < min2)                // x!=y==z
    {
      return E12_exmp1(x, y, z, exp_x, exp_y);
      } else if (sqr(center(x) - center(z)) < min2) // y!=z==x!=y
    {
      return E12_exmp1(y, x, z, exp_y, exp_x);
    } else
      return EX_111(x, y, z, exp_x) + EX_111(y, x, z, exp_y) +
             EX_111(z, y, x, exp_z);
  }
  static Tr_t<double> expm2(Tr_t<double> x) {
      assert(std::abs(center(x)) < 0.001);
    std::int8_t nmax = 7;
    Tr_t<double> xr = 1.0 / 2.0;
    Tr_t<double> sum = xr;
    for (std::int8_t n = 3; n < nmax; ++n) {
      xr *= x / n;
      sum += xr;
    }
    return sum;
  }

  static Tr_t<double> E12(Tr_t<double> x, Tr_t<double> y, Tr_t<double> exp_x, Tr_t<double> exp_y,
                    double min) {
      if (std::abs(center(y) - center(x)) < min) {
      return exp_x * expm2(y - x);
    }
    return EX_111(x, y, y, exp_x) + exp_y / (y - x) * (Tr_t<double>(1.0) - Tr_t<double>(1.0) / (y - x));
  }

  static Tr_t<double> E3(const Tr_t<double>& x, const Tr_t<double>& y, const Tr_t<double>& z, const Tr_t<double>& exp_x, const Tr_t<double>& exp_y,
                   const Tr_t<double>& exp_z, double min = 1e-3,
                   double eps = std::numeric_limits<double>::epsilon()) {
      if (sqr(center(x) - center(y)) < eps) // x==y
    {
          if (sqr(center(y) -center( z)) < eps) // y==z
        return exp_x / 2.0; // x==y==z
      else
        return E12(z, x, exp_z, exp_x, min); // x==y!=z
      } else if (sqr(center(y) - center(z)) < eps)             // x!=y==z
    {
      return E12(x, y, exp_x, exp_y, min);
      } else if (sqr(center(x) - center(z)) < eps) // y!=z==x!=y
    {
      return E12(y, x, exp_y, exp_x, min);
    } else
      return E111(x, y, z, exp_x, exp_y, exp_z, min); // x!=y!=z!=x
  }
  static Tr_t<double> expm2(Tr_t<double> x, Tr_t<double> y) {
      using std::expm1;
      if (std::abs(center(x) -center(y)) < std::numeric_limits<double>::epsilon() * 100)
          return Tr_t<double>(0.5) + (x + y) / 3.0;
      else if ((std::abs(center(x)) > 0.001) || (std::abs(center(y)) > 0.001)) {
          if (std::abs(center(x)) < std::numeric_limits<double>::epsilon() * 100) {
              if (std::abs(center(y)) < std::numeric_limits<double>::epsilon() * 100)
          return 0.5;
        else
          return (y - expm1(y)) / (y * (x - y));
          } else if (std::abs(center(y)) < std::numeric_limits<double>::epsilon() * 100) {
        return (x - expm1(x)) / (x * (y - x));
      } else
        return (y * expm1(x) - x * expm1(y)) / (x * y * (x - y));
    } else {
      const int8_t nmax = 10;
      Tr_t<double> xr = x / 2;
      Tr_t<double> yr = y / 2;
      Tr_t<double> sum (0);
      for (int8_t n = 3; n < nmax; ++n) {
        xr *= x / n;
        yr *= y / n;
        sum += xr - yr;
      }
      return Tr_t<double>(0.5) + sum / (x - y);
    }
  }

  static Tr_t<double> E3_expm2_order(const Tr_t<double>& x, const Tr_t<double>& y, const Tr_t<double>& z) {
      assert(center(x) <= center(y));
      assert(center(y) <= center(z));
    using std::exp;
    Tr_t<double> dx = x - y;
    Tr_t<double> dz = z - y;
    Tr_t<double> out = exp(y) * expm2(dx, dz);
    auto out2 = E3_expm1_order(x, y, z);
    // auto out3=E3(x,y,z,std::exp(x),std::exp(y),std::exp(z));
    assert(([&out, &out2]() { return are_Equal_v(out, out2, std::cerr); }()));
    return out;
  }

  static Tr_t<double>
  E3_expm1_order(Tr_t<double> x, Tr_t<double> y,Tr_t<double> z,
                 double eps = std::numeric_limits<double>::epsilon()) {
      assert(center(x) <= center(y));
      assert(center(y) <= center(z));
      using std::expm1;
    Tr_t<double> dx = x - y;
    Tr_t<double> dz = z - y;
    if (center(dx) * center(dx) < eps) {
        if (center(dz) * center(dz) < eps)
        return exp(y) / 2.0;
      else
          return exp(y) * (Tr_t<double>(1.0) - expm1(dz) / dz) / (-dz);
    } else if (center(dz) * center(dz) < eps)
      return exp(y) * (Tr_t<double>(1.0) - expm1(dx) / dx) / (-dx);
    else
      return exp(y) * (dx * expm1(dz) - dz * expm1(dx)) /
             (dx * dz * (dz - dx));
  }

  static Tr_t<double> E3(Tr_t<double> x, Tr_t<double> y, Tr_t<double> z) {
      if (center(y) < center(z)) {
          if (center(x) <center( y))
        return E3_expm2_order(x, y, z);
          else if (center(x) < center(z))
        return E3_expm2_order(y, x, z);
      else
        return E3_expm2_order(y, z, x);
      } else if (center(x) < center(z))
      return E3_expm2_order(x, z, y);
      else if (center(x) < center(y))
      return E3_expm2_order(z, x, y);
    else
      return E3_expm2_order(z, y, x);
  }

};




class Markov_Transition_step_double : public Markov_Transition_step {
  M_Matrix<double> gmean_i_;   // conductance matrix
  M_Matrix<double> gtotal_ij_; // conductance matrix
  M_Matrix<double> gmean_ij_;  // conductance matrix

  M_Matrix<double> gtotal_sqr_ij_; // conductance matrix

  M_Matrix<double> gsqr_i_; // conductance matrix

  M_Matrix<double> gvar_i_;        // variance of the conductance matrix
  M_Matrix<double> gtotal_var_ij_; // variance of the conductance matrix

  M_Matrix<double> gvar_ij_; // variance of the conductance matrix

public:
  Markov_Transition_step_double(M_Matrix<double> &&P, M_Matrix<double> &&g,
                                std::size_t nsamples, double fs, double min_p,
                                M_Matrix<double> &&gmean_i,
                                M_Matrix<double> &&gtotal_ij,
                                M_Matrix<double> &&gmean_ij,

                                M_Matrix<double> &&gtotal_sqr_ij,

                                M_Matrix<double> &&gsqr_i,

                                M_Matrix<double> &&gvar_i,
                                M_Matrix<double> &&gtotal_var_ij,
                                M_Matrix<double> &&gvar_ij)
      : Markov_Transition_step(std::move(P), std::move(g), nsamples, fs, min_p),
        gmean_i_{gmean_i}, gtotal_ij_{gtotal_ij}, gmean_ij_{gmean_ij},
        gtotal_sqr_ij_{gtotal_sqr_ij}, gsqr_i_{gsqr_i}, gvar_i_{gvar_i},
        gtotal_var_ij_{gtotal_var_ij}, gvar_ij_{gvar_ij} {}

  Markov_Transition_step_double_minimum getMinimum() const {
    return Markov_Transition_step_double_minimum(*this);
  }

  /// mean conductance for each starting state i
  M_Matrix<double> const &gmean_i() const {
    return gmean_i_;
  } // conductance matrix
  /// total conductance for each starting state i and ending state j
  M_Matrix<double> const &gtotal_ij() const {
    return gtotal_ij_;
  } // conductance matrix
  /// mean conductance for each starting state i and ending state j
  M_Matrix<double> const &gmean_ij() const {
    return gmean_ij_;
  } // conductance matrix

  /// squared mean conductance for each starting state i and ending state j
  M_Matrix<double> const &gtotal_sqr_ij() const {
    return gtotal_sqr_ij_;
  } // conductance matrix

  /// squared mean conductance for each starting state i
  M_Matrix<double> const &gsqr_i() const {
    return gsqr_i_;
  } // conductance matrix

  /// variance of the mean conductance for each starting state i
  M_Matrix<double> const &gvar_i() const {
    return gvar_i_;
  } // variance of the conductance matrix
  /// variance of the mean conductance for each starting state i contributed
  /// by the ones ending at state j
  M_Matrix<double> const &gtotal_var_ij() const {
    return gtotal_var_ij_;
  } // variance of the conductance matrix

  /// variance of the mean conductance for each starting state i and ending
  /// state j
  M_Matrix<double> const &gvar_ij() const {
    return gvar_ij_;
  } // variance of the conductance matrix

  Markov_Transition_step_double(const Markov_Transition_rate &Qx,
                                std::size_t nsamples, double fs, double min_p,
                                double max_landa_error) {
    if (Qx.error_landa() > max_landa_error)
      *this = Markov_Transition_step_double(1, Qx.Qrun(), Qx.g(), nsamples, fs,
                                            min_p);
    else
      *this = Markov_Transition_step_double(Qx, nsamples, fs, min_p, true);
  }

  Markov_Transition_step_double(const Markov_Transition_rate &Qx,
                                std::size_t nsamples, double fs, double min_p,
                                bool check = true)
      : Markov_Transition_step(Qx, nsamples, fs, min_p) {
    init(Qx, check);
  }

  Markov_Transition_step_double(std::size_t recursion_number,
                                const M_Matrix<double> &Qx,
                                const M_Matrix<double> &gx,
                                std::size_t nsamples, double fs, double min_p) {
    init_by_step(Qx, gx, recursion_number, nsamples, fs, min_p);
  }

  Markov_Transition_step_double(  Microscopic_description mi [[maybe_unused]],
                                [[maybe_unused]]const M_Matrix<double> &Qx,
                                [[maybe_unused]]const M_Matrix<double> &gx,
                                [[maybe_unused]]std::size_t nsamples, [[maybe_unused]]double fs, [[maybe_unused]]double min_p) {
    //   init_by_step(mi, Qx, gx, nsamples, fs, min_p);
  }

  Markov_Transition_step_double(std::mt19937_64 &mt,
                                const Markov_Transition_rate &Qx,
                                std::size_t nsamples, double fs, double min_p,
                                std::size_t numsteps, std::size_t numreplicates)
      : Markov_Transition_step(Qx, nsamples, fs, min_p) {
    measure(Qx, mt, numsteps, numreplicates);
  }

  Markov_Transition_step_double(Markov_Transition_step &&step)
      : Markov_Transition_step(step) {
    init_by_step_init();
  }

  Markov_Transition_step_double(Markov_Transition_step_double_minimum &&x,
                                M_Matrix<double> g, double fs, double min_p)
      : Markov_Transition_step(std::move(x.PPn), std::move(g), x.n, fs, min_p) {
    init(std::move(x), fs);
  }

  Markov_Transition_step_double() = default;

  typedef Markov_Transition_step_double self_type;
  typedef Markov_Transition_step base_type;
  constexpr static auto className =
      my_static_string("Markov_Transition_step_double");
  static auto get_constructor_fields() {
    // M_Matrix<double> const & (self_type::*myP) ()const& =&self_type::P;
    /*
 *   Markov_Transition_step_double(M_Matrix<double> &&P,
 * M_Matrix<double> &&g,
 * std::size_t nsamples,
 * double fs,
 *  double min_p,
 * M_Matrix<double>
     &&gmean_i, M_Matrix<double> &&gtotal_ij, M_Matrix<double> &&gmean_ij,

                                    M_Matrix<double> &&gtotal_sqr_ij,

                                    M_Matrix<double> &&gsqr_i,

                                    M_Matrix<double> &&gvar_i,
                                    M_Matrix<double> &&gtotal_var_ij,
                                    M_Matrix<double> &&gvar_ij)

 * */
    return std::make_tuple(
        grammar::field(C<self_type>{}, "P", &self_type::myP),
        grammar::field(C<self_type>{}, "g", &self_type::g),

        grammar::field(C<self_type>{}, "nsamples", &self_type::nsamples),
        grammar::field(C<self_type>{}, "fs", &self_type::fs),
        grammar::field(C<self_type>{}, "min_p", &self_type::min_P),
        grammar::field(C<self_type>{}, "gmean_i", &self_type::gmean_i),
        grammar::field(C<self_type>{}, "gtotal_ij", &self_type::gtotal_ij),
        grammar::field(C<self_type>{}, "gmean_ij", &self_type::gmean_ij),
        grammar::field(C<self_type>{}, "gtotal_sqr_ij",
                       &self_type::gtotal_sqr_ij),
        grammar::field(C<self_type>{}, "gsqr_i", &self_type::gsqr_i),
        grammar::field(C<self_type>{}, "gvar_i", &self_type::gvar_i),
        grammar::field(C<self_type>{}, "gtotal_var_ij",
                       &self_type::gtotal_var_ij),
        grammar::field(C<self_type>{}, "gvar_ij", &self_type::gvar_ij));
  }

  static double E1(double x) {
    if (std::abs(x) < std::numeric_limits<double>::epsilon() * 100)
      return 1.0;
    else if (std::abs(x) < 1e-2)
      return std::expm1(x) / x;
    else
      return (std::exp(x) - 1.0) / x;
  }

  static double E2(double x, double y) {
    const double eps = std::numeric_limits<double>::epsilon();
    if (x * x < eps) {
      if (y * y < eps)
        return 0.5;
      else
        return (E1(y) - 1.0) / y;
    } else if (y * y < eps)
      return (E1(x) - 1.0) / x;
    else if ((y - x) * (y - x) < eps)
      return (std::exp(x) - E1(x)) / x;
    else
      return (E1(y) - E1(x)) / (y - x);
  }

  static double Ee(double x, double y, double exp_x, double exp_y,
                   double eps = std::numeric_limits<double>::epsilon()) {
    if (sqr(x - y) < eps)
      return exp_x;
    else
      return (exp_x - exp_y) / (x - y);
  }

  static double EX_111(double x, double y, double z, double exp_x) {
    return exp_x / ((x - y) * (x - z));
  }
  static double E12_exmp1(double x, double y, double z, double exp_x,
                          double exp_y) {
    double dz = z - x;
    return (dz * (exp_x - exp_y) + (y - x) * exp_x * std::expm1(dz)) /
           ((x - y) * (y - z) * dz);
  }
  static double E111(double x, double y, double z, double exp_x, double exp_y,
                     double exp_z, double min2) {
    if (sqr(x - y) < min2) // x==y
    {
      if (sqr(y - z) < min2) // y==z
        return E3(x, y, z);  // x==y==z
      else
        return E12_exmp1(z, x, y, exp_z, exp_x); // x==y!=z
    } else if (sqr(y - z) < min2)                // x!=y==z
    {
      return E12_exmp1(x, y, z, exp_x, exp_y);
    } else if (sqr(x - z) < min2) // y!=z==x!=y
    {
      return E12_exmp1(y, x, z, exp_y, exp_x);
    } else
      return EX_111(x, y, z, exp_x) + EX_111(y, x, z, exp_y) +
             EX_111(z, y, x, exp_z);
  }
  static double expm2(double x) {
    assert(std::abs(x) < 0.001);
    std::int8_t nmax = 7;
    double xr = 1.0 / 2.0;
    double sum = xr;
    for (std::int8_t n = 3; n < nmax; ++n) {
      xr *= x / n;
      sum += xr;
    }
    return sum;
  }

  static double E12(double x, double y, double exp_x, double exp_y,
                    double min) {
    if (std::abs(y - x) < min) {
      return exp_x * expm2(y - x);
    }
    return EX_111(x, y, y, exp_x) + exp_y / (y - x) * (1.0 - 1.0 / (y - x));
  }

  static double E3(double x, double y, double z, double exp_x, double exp_y,
                   double exp_z, double min = 1e-3,
                   double eps = std::numeric_limits<double>::epsilon()) {
    if (sqr(x - y) < eps) // x==y
    {
      if (sqr(y - z) < eps) // y==z
        return exp_x / 2.0; // x==y==z
      else
        return E12(z, x, exp_z, exp_x, min); // x==y!=z
    } else if (sqr(y - z) < eps)             // x!=y==z
    {
      return E12(x, y, exp_x, exp_y, min);
    } else if (sqr(x - z) < eps) // y!=z==x!=y
    {
      return E12(y, x, exp_y, exp_x, min);
    } else
      return E111(x, y, z, exp_x, exp_y, exp_z, min); // x!=y!=z!=x
  }
  static double expm2(double x, double y) {
    if (std::abs(x - y) < std::numeric_limits<double>::epsilon() * 100)
      return 0.5 + (x + y) / 3.0;
    else if ((std::abs(x) > 0.001) || (std::abs(y) > 0.001)) {
      if (std::abs(x) < std::numeric_limits<double>::epsilon() * 100) {
        if (std::abs(y) < std::numeric_limits<double>::epsilon() * 100)
          return 0.5;
        else
          return (y - std::expm1(y)) / (y * (x - y));
      } else if (std::abs(y) < std::numeric_limits<double>::epsilon() * 100) {
        return (x - std::expm1(x)) / (x * (y - x));
      } else
        return (y * std::expm1(x) - x * std::expm1(y)) / (x * y * (x - y));
    } else {
      const int8_t nmax = 10;
      double xr = x / 2;
      double yr = y / 2;
      double sum = 0;
      for (int8_t n = 3; n < nmax; ++n) {
        xr *= x / n;
        yr *= y / n;
        sum += xr - yr;
      }
      return 0.5 + sum / (x - y);
    }
  }

  static double E3_expm2_order(double x, double y, double z) {
    assert(x <= y);
    assert(y <= z);
    double dx = x - y;
    double dz = z - y;
    double out = std::exp(y) * expm2(dx, dz);
    auto out2 = E3_expm1_order(x, y, z);
    // auto out3=E3(x,y,z,std::exp(x),std::exp(y),std::exp(z));
    assert(([&out, &out2]() { return are_Equal_v(out, out2, std::cerr); }()));
    return out;
  }

  static double
  E3_expm1_order(double x, double y, double z,
                 double eps = std::numeric_limits<double>::epsilon()) {
    assert(x <= y);
    assert(y <= z);
    double dx = x - y;
    double dz = z - y;
    if (dx * dx < eps) {
      if (dz * dz < eps)
        return std::exp(y) / 2;
      else
        return std::exp(y) * (1.0 - std::expm1(dz) / dz) / (-dz);
    } else if (dz * dz < eps)
      return std::exp(y) * (1.0 - std::expm1(dx) / dx) / (-dx);
    else
      return std::exp(y) * (dx * std::expm1(dz) - dz * std::expm1(dx)) /
             (dx * dz * (dz - dx));
  }

  static double E3(double x, double y, double z) {
    if (y < z) {
      if (x < y)
        return E3_expm2_order(x, y, z);
      else if (x < z)
        return E3_expm2_order(y, x, z);
      else
        return E3_expm2_order(y, z, x);
    } else if (x < z)
      return E3_expm2_order(x, z, y);
    else if (x < y)
      return E3_expm2_order(z, x, y);
    else
      return E3_expm2_order(z, y, x);
  }

private:
  void init(Markov_Transition_step_double_minimum &&x, double) {
    M_Matrix<double> u = ones<double>(P().ncols(), 1);

    gtotal_sqr_ij_ = x.PGG_n / x.n / x.n * 2;
    gtotal_ij_ = x.PG_n / x.n;
    gmean_ij_ = elemDivSafe(gtotal_ij_, P_, min_P());

    gtotal_var_ij_ = gtotal_sqr_ij_ - elemMult(gtotal_ij_, gmean_ij_);

    gmean_i_ = gtotal_ij_ * u;

    gsqr_i_ = gtotal_sqr_ij_ * u;
    gvar_ij_ = elemDivSafe(gtotal_var_ij_, P_, min_P());
    gmean_i_ = gtotal_ij_ * u;
    gsqr_i_ = gtotal_sqr_ij_ * u;
    gvar_i_ = gtotal_var_ij_ * u;
  }
  void normalize() {}

  void init(const Markov_Transition_rate &Qx, [[maybe_unused]]bool check) {
    // const double eps=std::numeric_limits<double>::epsilon();
    double dt = Markov_Transition_step::dt();

    std::size_t N = P().ncols();

    auto ladt = Qx.landa() * dt;
    double C_V = max(Qx.error_vector());
    double C_l = max(Qx.error_landa());
    double err = C_V * C_V * C_l * std::numeric_limits<double>::epsilon();

    assert(are_Equal_v(Qx.Qrun(), Qx.V() * Qx.landa() * Qx.W(), err, err,
                       std::cerr, "Qx.Qrun()=", Qx.Qrun(), "Qx.V()", Qx.V(),
                       "Qx.landa()=", Qx.landa(), "Qx.W()=", Qx.W()));
    [[maybe_unused]]double RV = Matrix_Unary_Functions::norm_1(Qx.V()) * Matrix_Unary_Functions::norm_1(Qx.W());

    auto exp_ladt = ladt.apply([](double x) { return std::exp(x); });

    M_Matrix<double> E2m(N, N, Matrix_TYPE::SYMMETRIC);
    // M_Matrix<double> E2mb(N, N, Matrix_TYPE::SYMMETRIC);
    for (std::size_t i = 0; i < N; ++i)
      for (std::size_t j = 0; j < i + 1; ++j)
        E2m(i, j) =
            Ee(ladt[i], ladt[j], exp_ladt[i], exp_ladt[j] /*, sqr(min_P())*/);

    // build E2
    M_Matrix<double> WgV_E2(N, N);

    for (std::size_t i = 0; i < N; ++i)
      for (std::size_t j = 0; j < N; ++j)
        WgV_E2(i, j) = Qx.WgV()(i, j) * E2m(i, j);
    M_Matrix<double> G = diag(g());
    auto gtotal_ij_p = (G * P() + P() * G) * 0.5;
    auto gtotal_var_ij_p =
        (G * P() - P() * G).apply([](double x) { return std::abs(x); }) * 0.5;
    auto gmean_ij_p = elemDivSafe(gtotal_ij_p, P_);
    auto gtotal_sqr_ij_p = gtotal_var_ij_p + elemMult(gtotal_ij_p, gmean_ij_p);

    gtotal_ij_ = (Qx.V() * WgV_E2 * Qx.W());
    gmean_ij_ = elemDivSafe(gtotal_ij_, P());

    M_Matrix<double> WgV_E3(N, N, 0.0);
    for (std::size_t n1 = 0; n1 < N; n1++)
      for (std::size_t n3 = 0; n3 < N; n3++)
        for (std::size_t n2 = 0; n2 < N; n2++) {
          double e3 = E3(ladt[n1], ladt[n2], ladt[n3], exp_ladt[n1],
                         exp_ladt[n2], exp_ladt[n3] /*, sqr(min_P())*/);
          WgV_E3(n1, n3) += Qx.WgV()(n1, n2) * Qx.WgV()(n2, n3) * e3;
          // optimizable
        }

    double minP = N * N * 20 * 1e-6;
    // std::pow(R_V, 3) * C_l *std::numeric_limits<double>::epsilon();
    gtotal_sqr_ij_ = Qx.V() * WgV_E3 * Qx.W() * 2.0;
    for (std::size_t i = 0; i < gtotal_sqr_ij_.size(); ++i) {
      if (std::abs(gtotal_ij_[i]) < minP)
        gtotal_ij_[i] = 0;
      if (std::abs(gtotal_sqr_ij_[i]) < minP)
        gtotal_sqr_ij_[i] = 0;
    }
    gtotal_var_ij_ = gtotal_sqr_ij_ - elemMult(gtotal_ij_, gmean_ij_);
    for (std::size_t i = 0; i < gtotal_sqr_ij_.size(); ++i) {
      if (std::abs(gtotal_sqr_ij_[i]) < minP)
        gtotal_var_ij_[i] = 0;
    }

    gvar_ij_ = elemDivSafe(gtotal_var_ij_, P());
    for (std::size_t i = 0; i < gvar_ij_.size(); ++i)
      if (std::abs(gvar_ij_[i]) < minP)
        gvar_ij_[i] = 0;

    M_Matrix<double> u(N, 1, 1.0);
    gmean_i_ = gtotal_ij_ * u;
    gsqr_i_ = gtotal_sqr_ij_ * u;
    gvar_i_ = gtotal_var_ij_ * u;
    auto r = elemDiv(gvar_ij_, gmean_ij_);
    if (false) {
      //      auto out2 =
      //          Markov_Transition_step_double(4, Qx, nsamples(), fs(),
      //          1e-7);
      assert(([this]() {
        for (std::size_t i = 0; i < gvar_ij_.size(); ++i)
          if (gvar_ij_[i] < 0)
            return false;
        return true;
      }()));
      //      assert(([this, &Qx, &check, &out2]() {
      //        return are_Equal<true, Markov_Transition_step_double>().test(
      //            *this, out2, 1e-7, std::cerr);
      //      }()));
    }
    // normalize();
  }

  void measure(const Markov_Transition_rate &Qx, std::mt19937_64 &mt,
               std::size_t numsteps, std::size_t nrep)

  {
    // const double eps=std::numeric_limits<double>::epsilon();
    double dt = Markov_Transition_step::dt() / numsteps;

    std::size_t N = P().ncols();

    auto ladt = Qx.landa() * dt;

    auto exp_ladt = ladt.apply([](double x) { return std::exp(x); });

    auto P =
        Probability_transition::normalize(Qx.V() * exp_ladt * Qx.W(), min_P());

    M_Matrix<double> gmean_ij_sum(N, N, 0.0);
    M_Matrix<double> gsqr_ij_sum(N, N, 0.0);
    M_Matrix<std::size_t> Ns(N, N, 0);
    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t r = 0; r < nrep; ++r) {
        double y = 0;
        markov_process<std::size_t> mp(i, P);
        for (std::size_t t = 0; t < numsteps; ++t) {
          auto y0 = g()[mp.N()];
          auto Nr = mp(mt);
          auto y1 = g()[Nr];
          y += (y0 + y1) / 2;
          mp.set_N(Nr);
        }
        y /= numsteps;
        gmean_ij_sum(i, mp.N()) += y;
        gsqr_ij_sum(i, mp.N()) += sqr(y);
        Ns(i, mp.N()) += 1;
      }
    }
    P_ = M_Matrix<double>(Ns) * (1.0 / nrep);
    gtotal_sqr_ij_ = gsqr_ij_sum * (1.0 / nrep);
    gtotal_ij_ = gmean_ij_sum * (1.0 / nrep);
    gmean_ij_ = elemDivSafe(gmean_ij_sum, Ns);
    gtotal_var_ij_ = gtotal_ij_ - elemMult(gmean_ij_, gtotal_ij_);
    gvar_ij_ = elemDivSafe(gsqr_ij_sum, Ns) - elemMult(gmean_ij_, gmean_ij_);
    M_Matrix<double> U = Matrix_Generators::ones<double>(1, g().size());
    M_Matrix<double> u(N, 1, 1.0);
    gmean_i_ = gtotal_ij_ * u;
    gsqr_i_ = gtotal_sqr_ij_ * u;
    gvar_i_ = gtotal_var_ij_ * u;
  }

  void init_by_step_init() {

    M_Matrix<double> G = diag(g());
    std::size_t N = P().nrows();
    gtotal_ij_ = (G * P() + P() * G) * 0.5;
    gtotal_var_ij_ =
        (G * P() - P() * G).apply([](double x) { return std::abs(x); }) * 0.5;
    gmean_ij_ = elemDivSafe(gtotal_ij_, P_);
    gtotal_sqr_ij_ = gtotal_var_ij_ + elemMult(gtotal_ij_, gmean_ij_);
    gvar_ij_ = elemDivSafe(gtotal_var_ij_, P_);
    M_Matrix<double> u(N, 1, 1.0);
    gmean_i_ = gtotal_ij_ * u;
    gsqr_i_ = gtotal_sqr_ij_ * u;
    gvar_i_ = gtotal_var_ij_ * u;
  }

  void init_by_step(const M_Matrix<double> &Qrun, const M_Matrix<double>& grun,
                    std::size_t recursion_number_extra, std::size_t n,
                    double fs, double min_p) {

    int recursion_number =
        2 + std::log2(maxAbs(Qrun) * 1.0 / fs * n) + recursion_number_extra;
    recursion_number = std::max(4, recursion_number);
    std::size_t times = std::pow(2, recursion_number);
    Markov_Transition_step_double s(Markov_Transition_step(
        calc_P_taylor(Qrun, n, fs * times), grun, n, fs * times, min_p));
    Markov_Transition_step_double step_ini(s);
    for (int i = 0; i < recursion_number; ++i) {
      Markov_Transition_step_double_minimum step(step_ini);
      step *= step_ini;
      step_ini = Markov_Transition_step_double(std::move(step), step_ini.g(),
                                               fs * times, min_p);
    }
    *this = step_ini;
  }
  void init_by_step([[maybe_unused]]Microscopic_description mi, const M_Matrix<double> &Qrun,
                    const M_Matrix<double>& grun,
                    std::size_t recursion_number_extra, std::size_t n,
                    double fs, double min_p) {

    int recursion_number =
        2 + std::log2(maxAbs(Qrun) * 1.0 / fs * n) + recursion_number_extra;
    recursion_number = std::max(4, recursion_number);
    std::size_t times = std::pow(2, recursion_number);
    Markov_Transition_step_double s(Markov_Transition_step(
        calc_P_taylor(Qrun, n, fs * times), grun, n, fs * times, min_p));
    Markov_Transition_step_double step_ini(s);
    for (int i = 0; i < recursion_number; ++i) {
      Markov_Transition_step_double_minimum step(step_ini);
      step *= step_ini;
      step_ini = Markov_Transition_step_double(std::move(step), step_ini.g(),
                                               fs * times, min_p);
    }
    *this = step_ini;
  }
};


template <bool output>
template <class ostream>
bool are_Equal<output, Markov_Transition_step_double>::test(
    const Markov_Transition_step_double &one,
    const Markov_Transition_step_double &other, double eps, ostream &os) const {
  are_Equal<output, M_Matrix<double>> mytest(eps, eps);
  assert(([&one]() {
    for (std::size_t i = 0; i < one.gvar_ij().size(); ++i)
      if (one.gtotal_var_ij()[i] < -0.1)
        return false;
    return true;
  }()));
  bool P_equal = mytest.test_prod(one.P(), other.P(), os);
  if constexpr (output)
    if (!P_equal)
      os << "\n differ in  P_ij !!  \n";
  bool gtotal_equal = mytest.test_prod(one.gtotal_ij(), other.gtotal_ij(), os);
  if constexpr (output)
    if (!gtotal_equal)
      os << "\n differ in  gtotal_ij !!\n";
  bool gsqr_equal =
      mytest.test_prod(one.gtotal_sqr_ij(), other.gtotal_sqr_ij(), os);
  if constexpr (output)
    if (!gsqr_equal)
      os << "\n differ in gtotal_sqr_ij!!\n";
  return gtotal_equal && gsqr_equal;
}

template <>
struct class_Invariants<Markov_Transition_step_double> : public invariant {
  template <bool output>
  static bool test_variances(
      const Markov_Transition_step_double &x,
      double eps = std::sqrt(std::numeric_limits<double>::epsilon())) {
    M_Matrix<double> u(x.gmean_i().nrows(), 1, 1.0);
    are_non_negative<output, M_Matrix<double>> non_negative(eps);
    are_Equal<output, M_Matrix<double>> areEqual(eps, eps);
    if (!non_negative.test_sum(x.gtotal_sqr_ij())) {
      if constexpr (output)
        std::cerr << "\n gtotal_sqr_ij is negative!!!! \n";
      return false;
    } else if (!non_negative.test_sum(x.gtotal_var_ij())) {
      if constexpr (output)
        std::cerr << "\n gtotal_var_ij is negative!!!! \n";
      return false;
    } else if (!areEqual.test_sum(x.gvar_i(), x.gtotal_var_ij() * u)) {
      if constexpr (output)
        std::cerr << "\n gvar_i is not equalt to  gtotal_var_ij*u !!!! \n";
      return false;

    } else if (!areEqual.test_sum(elemMult(x.gvar_ij(), x.P()),
                                  x.gtotal_var_ij())) {
      if constexpr (output)
        std::cerr << "\n gvar_ij is not equalt to  gtotal_var_ij*u !!!! \n";
      return false;
    } else if (!areEqual.test_sum(x.gtotal_var_ij() +
                                      elemMult(x.gmean_ij(), x.gtotal_ij()),
                                  x.gtotal_sqr_ij())) {
      if constexpr (output)
        std::cerr << "\n gtotal_var_ij+   elemMult(x.gmean_ij(),x.gtotal_ij()) "
                     "is not equalt to x.gtotal_sqr_ij() !!!! \n";
      return false;
    } else
      return true;
  }

  template <bool output>
  static bool
  test_means(const Markov_Transition_step_double &x,
             double eps = std::sqrt(std::numeric_limits<double>::epsilon())) {
    double gmin = min(x.g());
    double gmax = max(x.g());
    are_in_range<output, M_Matrix<double>> in_range(true, gmin, gmax, eps);
    bool out = true;
    are_Equal<output, M_Matrix<double>> equal(eps, eps);
    if (!are_not_more<output, M_Matrix<double>>(true, x.P() * gmax, eps)
             .test(x.gtotal_ij())) {
      if constexpr (output)
        std::cerr << "\n gtotal_ij is more than gmax=" << gmax << "!!! \n";
      out = false;
    }
    if (!are_not_less<output, M_Matrix<double>>(true, x.P() * gmin, eps)
             .test(x.gtotal_ij())) {
      if constexpr (output)
        std::cerr << "\n gtotal_ij is less than gmin=" << gmin << "!!! \n";
      out = false;
    }
    if (!in_range.test(x.gmean_i())) {
      if constexpr (output)
        std::cerr << "\n gmean_i is not in range!! \n";
      out = false;
    }
    if (!equal.test_prod(elemMult(x.gmean_ij(), x.P()), x.gtotal_ij())) {
      if constexpr (output)

        std::cerr << "\n gtotal_ij is not equal to  "
                     "elemMult(x.gmean_ij(),x.P())!! \n";
      out = false;
    }
    return out;
  }

  template <bool output>
  static bool
  test(const Markov_Transition_step_double &x,
       double eps = std::sqrt(std::numeric_limits<double>::epsilon())) {
    bool out = true;
    if (!test_variances<output>(x, eps)) {
      out = false;
    }
    if (!test_means<output>(x, eps)) {
      out = false;
    }
    if (out)
      return true;
    else {
      if constexpr (output)
        std::cerr << "<\n Markov_Transition_step_double=\n" << x;
      return false;
    }
  }
};

template<template<class...> class Tr>
inline Markov_Transition_step_variance_minimum_<Tr>::
    Markov_Transition_step_variance_minimum_(
        const Tr_t<Markov_Transition_step_variance> &x)
    : n{x.nsamples()}, PPn{x.P()}, PG_n{x.gtotal_ij() * x.nsamples()},
      PGG_n{x.gtotal_sqr_ij() * (x.nsamples() * x.nsamples() * 0.5)},
      min_p(x.min_P()) {}

template<template<class...> class Tr>
inline Markov_Transition_step_variance_minimum_<Tr> &
Markov_Transition_step_variance_minimum_<Tr>::
operator+=(const Tr_t<Markov_Transition_step_variance> &x) {
    auto n1 = x.nsamples();
    PGG_n = (PGG_n * x.P()) + (PG_n * x.gtotal_ij()) * n1 +
        (PPn * x.gtotal_sqr_ij()) * (0.5 * n1 * n1);
    PG_n = (PG_n * x.P()) + (PPn * x.gtotal_ij()) * n1;
    PPn = Probability_transition::normalize(PPn * x.P(), min_p);
    n = n + n1;
    return *this;
}






inline Markov_Transition_step_double_minimum::
    Markov_Transition_step_double_minimum(
        const Markov_Transition_step_double &x)
    : n{x.nsamples()}, PPn{x.P()}, PG_n{x.gtotal_ij() * x.nsamples()},
      PGG_n{x.gtotal_sqr_ij() * (x.nsamples() * x.nsamples() * 0.5)},
      min_p(x.min_P()) {}

inline Markov_Transition_step_double_minimum &
Markov_Transition_step_double_minimum::
operator*=(const Markov_Transition_step_double &x) {
  auto n1 = x.nsamples();
  PGG_n = (PGG_n * x.P()) + (PG_n * x.gtotal_ij()) * n1 +
          (PPn * x.gtotal_sqr_ij()) * (0.5 * n1 * n1);
  PG_n = (PG_n * x.P()) + (PPn * x.gtotal_ij()) * n1;
  PPn = Probability_transition::normalize(PPn * x.P(), min_p);
  n = n + n1;
  return *this;
}



class Microscopic_Transition_rate {
  M_Matrix<double> Qrun_; // transition rate matrix at time zero
  M_Matrix<double> g_;
  Microscopic_description mi_;
  std::vector<bool> surface_state_;

  M_Matrix<double> p0_;

public:
  auto &Qrun() const { return Qrun_; }
  auto &g() const { return g_; }
  auto &mi() const { return mi_; }
  auto &p0() const { return p0_; }
  Microscopic_Transition_rate(M_Matrix<double> &&Qrun__, M_Matrix<double> &&g__)
      : Qrun_{std::move(Qrun__)}, g_{std::move(g__)} {}
};


template<template<class...> class>
class Microscopic_description_new_;

typedef Microscopic_description_new_<C> Microscopic_description_new;


template<template<class...> class Tr>
class Microscopic_description_new_
{
public:
    template<typename T> using Tr_t=typename Tr<T>::type;
  typedef std::vector<std::size_t> mistate;
  typedef Microscopic_description_new_<Tr> self_type;
  std::size_t number_of_states() const { return k_; }
  std::size_t size() const { return ns_.size(); }
  std::size_t number_of_channels() const { return N_; }
  std::size_t index(const mistate &n) const {
    auto it = map_.find(n);
    if (it != map_.end())
      return it->second;
    else
      return size();
  }
  bool has_it(const mistate &n) const {
    assert(n.size() == number_of_states());
    assert(sum(n) == number_of_channels());
    return map_.find(n) != map_.end();
  }

  void push_back(const mistate &n) {
    if (!has_it(n)) {
      ns_.push_back(n);
      map_[ns_.back()] = ns_.size() - 1;
    }
  }
  void push_back(mistate &&n) {
    if (!has_it(n)) {
      ns_.push_back(std::move(n));
      map_[ns_.back()] = ns_.size() - 1;
    }
  }
  auto convert_to_Macroscopic(const Tr_t<M_Matrix<double>> &P) const {
    assert(size() == P.size());
    M_Matrix<double> Nmean(1, number_of_states(), 0.0);
    M_Matrix<double> Nsqrmean(number_of_states(), number_of_states(),
                              Matrix_TYPE::SYMMETRIC, 0);

    for (std::size_t i = 0; i < size(); ++i) {
      M_Matrix<double> p =
          M_Matrix<std::size_t>(1ul, number_of_states(), ns_[i]);
      Nmean += p * P[i];
      Nsqrmean += quadraticForm_XTX(p) * P[i];
    }
    M_Matrix<double> Pmean = Nmean * (1.0 / number_of_channels());
    M_Matrix<double> Pcov =
        (Nsqrmean - quadraticForm_XTX(Nmean)) * (1.0 / number_of_channels());

    using namespace std::literals::string_literals;
    //  auto PP=convert_to_Microscopic(Pmean,Pcov);
    //    are_Equal_v(P,PP,
    //                       std::numeric_limits<double>::epsilon()*1e10,
    //                       std::numeric_limits<double>::epsilon()*1e10,std::cerr,
    //                       "P"s,P,"ns",M_Matrix<std::size_t>(ns_),"convert_to_Microscopic(Pmean,Pcov)"s,convert_to_Microscopic(Pmean,Pcov),
    //                       "Pmean"s,Pmean,"Pcov"s,Pcov);

    return std::make_tuple(Pmean, Pcov);
  }
  auto maximum_entropy(const Tr_t<M_Matrix<double>> &Pmean, std::size_t maxiter,
                       double maxGradient, double max_landa,
                       double maxstep) const {
    typedef myOptional_t<M_Matrix<double>> Op;
    auto Pij = P_matrix();
    maxent::end_criteria end(maxiter, maxGradient, max_landa);
    auto res = maxent::optimize(Pmean, Pij, end, maxstep);
    if (!res)
      return Op(false, res.error());
    else
      return Op(res.value().pi());
  }
  auto maximum_entropy(const M_Matrix<double> &Pmean,
                       const M_Matrix<double> &Pcov, std::size_t maxiter,
                       double maxGradient, double maxlanda,
                       double maxstep) const {
    typedef myOptional_t<M_Matrix<double>> Op;
    auto Pij = P_P2_matrix();
    maxent::end_criteria end(maxiter, maxGradient, maxlanda);
    auto PP2 = P_P2_matrix(Pmean, Pcov);
    auto res = maxent::optimize(PP2, Pij, end, maxstep);
    if (!res)
      return Op(false, res.error());
    else
      return Op(res.value().pi());
  }
  M_Matrix<double> P_matrix() const {

    M_Matrix<double> out(size(), number_of_states());
    for (std::size_t i = 0; i < size(); ++i)
      for (std::size_t j = 0; j < number_of_states(); ++j)
        out(i, j) = ns_[i][j] * (1.0 / number_of_channels());
    return out;
  }

  M_Matrix<double> P_P2_matrix() const {
    std::size_t k = (number_of_states() * (number_of_states() + 3)) / 2;
    M_Matrix<double> out(size(), k);
    for (std::size_t i = 0; i < size(); ++i) {

      auto p = M_Matrix<double>(1, number_of_states());
      for (std::size_t j = 0; j < p.size(); ++j)
        p[j] = ns_[i][j] * (1.0 / number_of_channels());
      for (std::size_t j = 0; j < number_of_states(); ++j)
        out(i, j) = p[j];
      auto p2 = quadraticForm_XTX(p);
      //  auto p2=quadraticForm_XTX(p)*number_of_channels();
      for (auto j = 0ul; j < p2.size(); ++j) {
        out(i, number_of_states() + j) = p2[j];
      }
    }
    return out;
  }
  M_Matrix<double> P_P2_matrix(const M_Matrix<double> &Pmean,
                               const M_Matrix<double> &Pcov) const {
    //         auto P2=quadraticForm_XTX(Pmean)*number_of_channels()+Pcov;
    auto P2 = quadraticForm_XTX(Pmean) + Pcov / number_of_channels();
    std::size_t k = (number_of_states() * (number_of_states() + 3)) / 2;
    M_Matrix<double> out(1, k);
    for (std::size_t j = 0; j < number_of_states(); ++j)
      out[j] = Pmean[j];
    for (auto j = 0ul; j < P2.size(); ++j)
      out[number_of_states() + j] = P2[j];
    return out;
  }

  auto convert_to_Microscopic_maxent(const M_Matrix<double> &Pmean,
                                     const M_Matrix<double> &Pcov) const {
    if (size() == 1)
      return myOptional_t<M_Matrix<double>>(M_Matrix<double>(1, 1, 1.0));
    else
      return maximum_entropy(Pmean, Pcov, 100, 1e-7, 16, 1.0);
  }

  Microscopic_description_new_() = default;

  Microscopic_description_new_(std::size_t numChannels,
                               const Tr_t<M_Matrix<double>> &Pmean,
                               const Tr_t<M_Matrix<double>> &Pcov, double maxChi2)
      : N_{numChannels}, k_{Pmean.size()}, ns_{}, map_{} {
    auto s = get_states_as_required(Pmean, Pcov, maxChi2);
    for (auto &e : s)
      push_back(e);
  }

  Microscopic_description_new_(std::size_t N__, std::size_t k__,
                              std::vector<std::vector<std::size_t>> &&ns__)
      : N_{N__}, k_{k__}, ns_{std::move(ns__)}, map_{get_map(ns_)} {}

  void reduce(Tr_t<M_Matrix<double>> &P, double remove_p) {
    assert(P.size() == size());
    remove_low_probabilities(ns_, P, remove_p);
    map_ = get_map(ns_);
    assert(P.size() == size());
  }

  void reduce(const std::set<std::size_t> &to_be_removed_) {
    if (!to_be_removed_.empty()) {
      auto to_be_removed = to_be_removed_;
      auto newN = ns_.size() - to_be_removed.size();
      std::vector<mistate> new_ns(newN);
      std::size_t start = 0;
      std::size_t i_destination = 0;
      to_be_removed.insert(ns_.size()); // hack to force to include the elements
                                        // past the last removed.
      for (auto e : to_be_removed) {
        for (std::size_t i_origin = start; i_origin < e; ++i_origin) {
          new_ns[i_destination] = ns_[i_origin];
          ++i_destination;
        }
        start = e + 1;
      }

      ns_ = std::move(new_ns);
      map_ = get_map(ns_);
    }
  }

  static bool has_unallowed_states(
      const std::vector<std::size_t> &n, const std::vector<double> &Nmean,
      double eps = std::numeric_limits<double>::epsilon() * 100) {
    for (std::size_t i = 0; i < n.size(); ++i) {
      if ((Nmean[i] < eps) && (n[i] > 0))
        return true;
    }
    return false;
  }

  static double chi2sum(const std::vector<std::size_t> &n,
                        const std::vector<double> &Nmean,
                        const M_Matrix<double> &pinvNCov) {
    return xTSigmaX(Nmean - n, pinvNCov);
  }

  auto search_neighbours_up(
      std::set<mistate> &rejected, std::set<mistate> &accepted,
      const mistate &origin,
      const std::vector<std::vector<std::size_t>> &connections_x,
      const std::vector<double> &Nmean, const M_Matrix<double> &pinvNCov,
      double maxChi2) const {
    std::vector<mistate> out;
    std::set<mistate> new_rejected;
    for (std::size_t i_origin = 0; i_origin < origin.size(); ++i_origin) {
      if (origin[i_origin] > 0)
        for (auto &i_destination : connections_x[i_origin]) {
          auto destination = origin;
          --destination[i_origin];
          ++destination[i_destination];
          if (!has_it(destination) &&
              (rejected.find(destination) == rejected.end()) &&
              (accepted.find(destination) == accepted.end())) {
            if (has_unallowed_states(destination, Nmean)) {
              new_rejected.insert(destination);
            } else {
              double logL = chi2sum(destination, Nmean, pinvNCov);
              if (logL < maxChi2) {
                out.push_back(destination);
                accepted.insert(destination);
              } else
                new_rejected.insert(destination);
            }
          }
        }
    }
    rejected.insert(new_rejected.begin(), new_rejected.end());
    return out;
  }

  std::set<mistate>
  search_on_normal(const mistate &origin,
                   const std::vector<std::vector<std::size_t>> &connections_x,
                   const std::vector<double> &Nmean,
                   const M_Matrix<double> &pinvNCov, double maxChi2) const {
    std::set<mistate> rejected;
    std::set<mistate> accepted;
    accepted.insert(origin);
    auto neighbours_list = search_neighbours_up(
        rejected, accepted, origin, connections_x, Nmean, pinvNCov, maxChi2);
    while (!neighbours_list.empty()) {

      std::vector<std::vector<std::size_t>> newNeigh;
      for (auto &e : neighbours_list) {
        auto new_neighbour_up = search_neighbours_up(
            rejected, accepted, e, connections_x, Nmean, pinvNCov, maxChi2);
        newNeigh.insert(newNeigh.end(), new_neighbour_up.begin(),
                        new_neighbour_up.end());
      }
      neighbours_list = newNeigh;
    }
    return accepted;
  }

  mistate to_n(const M_Matrix<double> &Pmean) const {
    std::vector<std::size_t> out(Pmean.size());
    std::multimap<double, std::size_t> remain;
    std::size_t sum = 0;
    for (std::size_t i = 0; i < out.size(); ++i) {
      double n = number_of_channels() * Pmean[i];
      out[i] = n;
      sum += out[i];
      remain.insert(std::make_pair(n - out[i], i));
    }
    auto it = remain.rbegin();
    while (sum < number_of_channels()) {
      auto i = it->second;
      ++out[i];
      ++sum;
      ++it;
    }
    return out;
  }

  std::set<mistate> get_states_as_required(const M_Matrix<double> &Pmean,
                                           const M_Matrix<double> &Pcov,
                                           double maxChi2) const {
    std::vector<std::vector<std::size_t>> connections_x(number_of_states());
    for (std::size_t i = 0; i < number_of_states(); ++i)
      for (std::size_t j = 0; j < number_of_states(); ++j)
        if (j != i)
          connections_x[i].push_back(j);

    std::vector<double> Nmean = (Pmean * number_of_channels()).toVector();
    auto origin = to_n(Pmean);
    std::set<mistate> accepted;
    accepted.insert(origin);
    if (maxAbs(Pcov) > std::numeric_limits<double>::epsilon() * 1000) {
      auto pinvNCov = pinv(Pcov * number_of_channels());
      if (pinvNCov.has_value() && (maxAbs(pinvNCov.value()) > 0)) {
        accepted = search_on_normal(origin, connections_x, Nmean,
                                    pinvNCov.value(), maxChi2);
        //  search_down(origin, connections_down, Nmean, pinvNCov.value(),
        //  maxChi2);
      }
    }
    //      assert(([&accepted, &Pmean, &Pcov, &maxChi2, this]() {
    //          auto mifull = add_states_as_required_full_normal(
    //              number_of_channels(), number_of_states(), Pmean, Pcov,
    //              maxChi2);
    //          auto set_dif = symmetric_difference(mifull, accepted);
    //          if (!set_dif.empty()) {
    //              std::cerr << "\nPmean =" << Pmean << "\nPcov= " << Pcov <<
    //              "\nmaxChi2"
    //                        << maxChi2 << "\nset_dif" << set_dif << "\n
    //                        smart=" << ns_
    //                        << "\n full=" << mifull.ns_ << "accepted \n"
    //                        << accepted << "\n";

    //              return false;
    //          }
    //          return true;ay
    //      }()));
    return accepted;
  }

  auto &ns(std::size_t i) const { return ns_[i]; }
  auto &ns() const { return ns_; }

  struct Less {
    bool operator()(const mistate &A, const mistate &B) const {
      assert(A.size() == B.size());
      auto k = A.size();
      for (std::size_t i = k; i > 0; --i) {
        if ((A)[i - 1] < (B)[i - 1])
          return true;
        else if ((B)[i - 1] < (A)[i - 1])
          return false;
      }
      return false;
    }
  };
  static std::map<std::vector<std::size_t>, std::size_t, Less>
  get_map(const std::vector<std::vector<std::size_t>> &ns) {
    std::map<std::vector<std::size_t>, std::size_t, Less> out;
    for (std::size_t i = 0; i < ns.size(); ++i) {
      out[ns[i]] = i;
    }
    return out;
  }

private:
  std::size_t N_;
  std::size_t k_;
  std::vector<std::vector<std::size_t>> ns_;
  std::map<std::vector<std::size_t>, std::size_t, Less> map_;
};

class Microscopic_Transition_step;

struct Microscopic_Transition_step_double_minimum {
  std::size_t n;
  Microscopic_description_new mi0;
  Microscopic_description_new mi1;
  M_Matrix<double> PPn;
  M_Matrix<double> PG_n;
  M_Matrix<double> PGG_n;
  double min_p;
  double tolerance;
  Microscopic_Transition_step_double_minimum(
      const Microscopic_Transition_step &x_);

  Microscopic_Transition_step_double_minimum &
  operator*=(const Microscopic_Transition_step &x_);
};

class Microscopic_Transition_step_double {
  Microscopic_description_new mi0_;
  Microscopic_description_new mi1_;
  double min_P_;
  M_Matrix<double> P_;
  M_Matrix<double> gtotal_ij_;
  M_Matrix<double> gtotal_sqr_ij_;

  M_Matrix<double> gmean_ij_;
  M_Matrix<double> gtotal_var_ij_;

  M_Matrix<double> gmean_i_;

  M_Matrix<double> gvar_ij_;
  M_Matrix<double> gvar_i_;

  void calc() {
    gmean_ij_ = elemDivSafe(gtotal_ij_, P_, min_P());

    gtotal_var_ij_ = gtotal_sqr_ij_ - elemMult(gtotal_ij_, gmean_ij_);
    std::size_t N = P_.ncols();
    M_Matrix<double> u(N, 1, 1.0);

    gmean_i_ = gtotal_ij_ * u;

    gvar_ij_ = elemDivSafe(gtotal_var_ij_, P_, min_P());
    gvar_i_ = gtotal_var_ij_ * u;
  }

public:
  Microscopic_Transition_step_double(Microscopic_description_new mi0,
                                     Microscopic_description_new &&mi1,
                                     M_Matrix<double> &&P,
                                     M_Matrix<double> &&PG,
                                     M_Matrix<double> &&PGG)
      : mi0_(std::move(mi0)), mi1_(std::move(mi1)), P_{std::move(P)},
        gtotal_ij_{std::move(PG)}, gtotal_sqr_ij_{std::move(PGG)} {
    calc();
  }

  double min_P() const { return min_P_; }
  auto &P() const { return P_; }
  auto &gmean_i() const { return gmean_i_; }
  auto &gmean_ij() const { return gmean_ij_; }
  auto &gvar_ij() const { return gvar_ij_; }
  auto &gvar_i() const { return gvar_i_; }
  Microscopic_description_new const &mi0() const { return mi0_; }
  Microscopic_description_new const &mi1() const { return mi1_; }
};

template<template<class...> class Tr>
class Microscopic_Transition_step_mean_;

typedef Microscopic_Transition_step_mean_<C> Microscopic_Transition_step_mean;


template <template<class...> class > class Microscopic_Transition_step_mean_minimum_;

template <template<class...> class Tr>
class Microscopic_Transition_step_mean_minimum_ {
public:
    template <typename T> using Tr_t=typename Tr<T>::type;
    std::size_t n;
    Tr_t<Microscopic_description_new> mi0;
    Tr_t<Microscopic_description_new> mi1;
    Tr_t<M_Matrix<double>> PPn;
    Tr_t<M_Matrix<double>> PG_n;
    Tr_t<M_Matrix<double>> PGG_n;
    double min_p;
    double tolerance;
    Microscopic_Transition_step_mean_minimum_(
        const Tr_t<Microscopic_Transition_step_mean> &x_);

    Microscopic_Transition_step_mean_minimum_ &
    operator+=(const Tr_t<Microscopic_Transition_step_mean> &x_);
};
template<template<class...> class Tr>
    class Microscopic_Transition_step_mean_ :public Markov_Transition_step_mean_<Tr>
{
public:
    template <typename T> using Tr_t=typename Tr<T>::type;
private:
    Tr_t<Microscopic_description_new> mi0_;
    Tr_t<Microscopic_description_new> mi1_;

public:
    Tr_t<Microscopic_description_new> const &mi0() const { return mi0_; }
    Tr_t<Microscopic_description_new> const &mi1() const { return mi1_; }
};

template<template<class...> class Tr>
class Microscopic_Transition_step_variance_;

typedef Microscopic_Transition_step_variance_<C> Microscopic_Transition_step_variance;

template<template<class...> class Tr>
class Microscopic_Transition_step_variance_minimum_;

typedef Microscopic_Transition_step_variance_minimum_<C> Microscopic_Transition_step_variance_minimum;

template<template<class...> class Tr>
class Microscopic_Transition_step_variance_minimum_ {
public:
    template<typename T> using Tr_t=typename Tr<T>::type;
    std::size_t n;
    Tr_t<Microscopic_description_new> mi0;
    Tr_t<Microscopic_description_new> mi1;
    Tr_t<M_Matrix<double>> PPn;
    Tr_t<M_Matrix<double>> PG_n;
    Tr_t<M_Matrix<double>> PGG_n;
    double min_p;
    double tolerance;
    Microscopic_Transition_step_variance_minimum_(
        const Tr_t<Microscopic_Transition_step_variance> &x_);

    Microscopic_Transition_step_variance_minimum_ &
    operator+=(const Tr_t<Microscopic_Transition_step_variance> &x_);
};

template<template<class...> class Tr>
class Microscopic_Transition_step_variance_: public Markov_Transition_step_variance_<Tr>
{
    template <typename T> using Tr_t=typename Tr<T>::type;
    Tr_t<Microscopic_description_new> mi0_;
    Tr_t<Microscopic_description_new> mi1_;


public:

    Tr_t<Microscopic_description_new> const &mi0() const { return mi0_; }
    Tr_t<Microscopic_description_new> const &mi1() const { return mi1_; }
};

template<template<class...> class Tr>
struct mi_Q_;

typedef mi_Q_<C> mi_Q;

template<template<class...> class Tr>
struct mi_Q_ {
    template<typename T> using Tr_t=typename Tr<T>::type;
  std::vector<std::size_t> mi_start;
  std::vector<std::size_t> mi_end;
  std::vector<std::size_t> num_rows;
  std::size_t mi_current;
  Tr_t<Microscopic_description_new> mi;
  Tr_t<M_Matrix<double>> prun;

  Tr_t<M_Matrix<M_Matrix<double>>> Q_s;
  Tr_t<M_Matrix<M_Matrix<double>>> P_s;
  Tr_t<M_Matrix<M_Matrix<double>>> Od_s;

  mi_Q_() = default;
  mi_Q_(Tr_t<Microscopic_description_new> const &start)
      : mi_start(), mi_end(), num_rows(),
        mi_current(0), mi(start), Q_s(),
        P_s(),
        Od_s() {}
};

inline auto
calc_Q(Microscopic_description_new const &start, const M_Matrix<double> &Qm,
       const std::vector<std::vector<std::size_t>> &connections_to_x,
       std::size_t iter_num) {
  Microscopic_description_new mi_out(start);

  std::vector<std::vector<double>> Qo;

  std::vector<M_Matrix<double>> Q_s(iter_num);
  std::size_t initial_row = 0;
  std::vector<std::size_t> nrows(iter_num);

  for (std::size_t ncycle = 0; ncycle < iter_num; ++ncycle) {
    nrows[ncycle] = mi_out.size();
    auto final_row = nrows[ncycle];
    std::vector<std::vector<double>> Qapp(
        final_row - initial_row, std::vector<double>(mi_out.size(), 0.0));
    Qo.insert(Qo.end(), Qapp.begin(), Qapp.end());
    for (std::size_t i = initial_row; i < final_row; ++i) {
      auto &origin = mi_out.ns(i);
      auto n = mi_out.number_of_states();
      for (std::size_t i_origin = 0; i_origin < n; ++i_origin) {
        if (origin[i_origin] > 0) {
          for (auto i_destination : connections_to_x[i_origin]) {
            auto destination = origin;
            destination[i_origin]--;
            destination[i_destination]++;
            auto j = mi_out.index(destination);
            if (j == mi_out.size()) {
              mi_out.push_back(destination);
              for (auto &e : Qo)
                e.push_back(0.0);
            }
            double q = origin[i_origin] * Qm(i_origin, i_destination);
            Qo[i][j] = q;
            Qo[i][i] -= q;
          }
        }
      }
    }
    initial_row = final_row;
    Q_s[iter_num - 1 - ncycle] = M_Matrix<double>(Qo);
    for (std::size_t i = 0; i < ncycle; ++i)
      Q_s[iter_num - 1 - i] =
          Q_s[iter_num - 1 - i] * Q_s[iter_num - 1 - ncycle];
  }

  //      assert(([&Q0, &Qa, &Q0m, &Qam, this]() {
  //          auto [Q0t, Qat] = this->Qs(Q0m, Qam);
  //          return are_Equal_v(Q0, Q0t,
  //          std::numeric_limits<double>::epsilon(),
  //                             std::numeric_limits<double>::epsilon(),
  //                             std::cerr) &&
  //              are_Equal_v(Qa, Qat, std::numeric_limits<double>::epsilon(),
  //                          std::numeric_limits<double>::epsilon(),
  //                          std::cerr);
  //      }()));

  return std::make_tuple(M_Matrix<double>(Qo), mi_out, nrows);
}

inline void reduce_P_mi(Microscopic_description_new &mi1, M_Matrix<double> &P1,
                        double reduce_by_p) {
  auto to_be_removed = get_low_probabilities_indexes(mi1.ns(), P1, reduce_by_p);
  if (!to_be_removed.empty()) {
    auto newN = mi1.size() - to_be_removed.size();
    M_Matrix<double> newP(1, newN);
    std::size_t start = 0;
    std::size_t i_destination = 0;
    to_be_removed.insert(mi1.size());
    for (auto e : to_be_removed) {
      for (std::size_t i_origin = start; i_origin < e; ++i_origin) {
        newP[i_destination] = P1[i_origin];
        ++i_destination;
      }
      start = e + 1;
    }
    mi1.reduce(to_be_removed);
    P1 = std::move(newP);
  }
}

inline void reduce_P_by_p(Microscopic_description_new &mi1,
                             M_Matrix<double> &Prun, M_Matrix<double> &p1,
                             double reduce_by_p, std::size_t i_begin = 0,
                             std::size_t i_end = std::string::npos) {
  if (i_end == std::string::npos)
    i_end = mi1.size();
  if (i_end > i_begin) {

    auto to_be_removed = get_low_probabilities_indexes(
        mi1.ns(), p1, reduce_by_p, i_begin, i_end);
    if (!to_be_removed.empty()) {
      auto newN = p1.size() - to_be_removed.size();
      mi1.reduce(to_be_removed);
      if (newN>0)
      {
      M_Matrix<double> newP(1, newN);
      M_Matrix<double> newPrun(Prun.nrows(), newN);
      std::size_t start = 0;
      std::size_t i_destination = 0;
      to_be_removed.insert(mi1.size());
      for (auto e : to_be_removed) {
        for (std::size_t i_origin = start; i_origin < e-i_begin; ++i_origin) {
          newP[i_destination] = p1[i_origin];
          ++i_destination;
        }
        start = e + 1;
      }
      for (std::size_t i = 0; i < Prun.nrows(); ++i) {
        std::size_t start = 0;
        std::size_t i_destination = 0;
        for (auto e : to_be_removed) {
          for (std::size_t i_origin = start; i_origin < e-i_begin; ++i_origin) {
            newPrun(i, i_destination) = Prun(i, i_origin);
            ++i_destination;
          }
          start = e + 1;
        }
      }
      p1 = std::move(newP);
      Prun = std::move(newPrun);
      }
      else {
          p1=M_Matrix<double>();
          Prun=M_Matrix<double>();
      }
    }
  }
}


inline void reduce_P_G_by_p(Microscopic_description_new &mi1,M_Matrix<double> &p1,
                          M_Matrix<double> &Prun, M_Matrix<double> &Gtot,M_Matrix<double> &Gvar,
                          double reduce_by_p, std::size_t i_begin = 0,
                          std::size_t i_end = std::string::npos) {
    if (i_end == std::string::npos)
        i_end = mi1.size();
    if (i_end > i_begin) {

        auto to_be_removed = get_low_probabilities_indexes(
            mi1.ns(), p1, reduce_by_p, i_begin, i_end);
        if (!to_be_removed.empty()) {
            auto newN = p1.size() - to_be_removed.size();
            mi1.reduce(to_be_removed);
            if (newN>0)
            {
                M_Matrix<double> newP(1, newN);
                M_Matrix<double> newPrun(Prun.nrows(), newN);
                M_Matrix<double> newGtot(Prun.nrows(), newN);
                M_Matrix<double> newGvar(Prun.nrows(), newN);
                std::size_t start = 0;
                std::size_t i_destination = 0;
                to_be_removed.insert(mi1.size());
                for (auto e : to_be_removed) {
                    for (std::size_t i_origin = start; i_origin < e-i_begin; ++i_origin) {
                        newP[i_destination] = p1[i_origin];
                        ++i_destination;
                    }
                    start = e + 1;
                }
                for (std::size_t i = 0; i < Prun.nrows(); ++i) {
                    std::size_t start = 0;
                    std::size_t i_destination = 0;
                    for (auto e : to_be_removed) {
                        for (std::size_t i_origin = start; i_origin < e-i_begin; ++i_origin) {
                            newPrun(i, i_destination) = Prun(i, i_origin);
                            newGtot(i, i_destination) = Gtot(i, i_origin);
                            newGvar(i, i_destination) = Gvar(i, i_origin);
                            ++i_destination;
                        }
                        start = e + 1;
                    }
                }
                p1 = std::move(newP);
                Prun = std::move(newPrun);
                Gtot=std::move(newGtot);
                Gvar=std::move(newGvar);
            }
            else {
                p1=M_Matrix<double>();
                Prun=M_Matrix<double>();
                Gtot=M_Matrix<double>();
                Gvar=M_Matrix<double>();
            }
        }
    }
}



inline void pre_Qmi(std::size_t n, mi_Q &m,
                    std::vector<std::vector<M_Matrix<double>>> &Q_sm) {
    m.mi_start.push_back(m.mi_current);
    m.mi_end.push_back(m.mi.size());
    m.num_rows.push_back(m.mi_end[n] - m.mi_start[n]);
  Q_sm.push_back(std::vector<M_Matrix<double>>(n + 1));
  for (std::size_t nj = 0; nj + 1 < n; ++nj)
    Q_sm[n][nj] =
        M_Matrix<double>(m.num_rows[n], m.num_rows[nj], Matrix_TYPE::ZERO);

  if (n > 0)
    Q_sm[n][n - 1] = M_Matrix<double>(m.num_rows[n], m.num_rows[n - 1]);
  Q_sm[n][n] = M_Matrix<double>(m.num_rows[n], m.num_rows[n], 0.0);
}

inline void post_Qmi(std::size_t n, double reduce_p_by, mi_Q &m,
                     std::vector<std::vector<M_Matrix<double>>> &Q_sm,
                     const std::vector<std::vector<double>>& Q_01) {

  if (!Q_01[0].empty())
  {
      M_Matrix<double> Q01=M_Matrix<double>(Q_01);
  if (n < 2)
    m.prun = m.prun * Q01;
  else
    m.prun = (m.prun * Q01) * (1.0 / n);

  reduce_P_by_p(m.mi, Q01, m.prun, reduce_p_by, m.mi_end[n], m.mi.size());

  if (Q01.size()>0)
      Q_sm[n].push_back(std::move(Q01));
  }
  m.mi_current = m.mi_end[n];

}

inline void
calc_Qmi(std::size_t n, std::size_t nsteps,
         const std::vector<std::vector<std::size_t>> &connections_to_x,
         const M_Matrix<double> &Qm, mi_Q &m,
         std::vector<std::vector<M_Matrix<double>>> &Q_sm,
         std::vector<std::vector<double>> &Q_01) {

  for (std::size_t i = 0; i < m.num_rows[n]; ++i) {
    auto origin = m.mi.ns(i + m.mi_start[n]);
    auto nu = m.mi.number_of_states();

    for (std::size_t i_origin = 0; i_origin < nu; ++i_origin) {
      if (origin[i_origin] > 0) {
        for (auto i_destination : connections_to_x[i_origin]) {
          auto destination = origin;
          destination[i_origin]--;
          destination[i_destination]++;
          auto j = m.mi.index(destination);
          if (j == m.mi.size()) {
            m.mi.push_back(destination);
            for (auto &e : Q_01)
              e.push_back(0.0);
          }
          double q = origin[i_origin] * Qm(i_origin, i_destination) / nsteps;
          if (j >= m.mi_end[n]) {
            Q_01[i][j - m.mi_end[n]] = q;
            Q_sm[n][n](i, i) -= q;

          } else
            for (std::size_t nj = 0; nj <= n; ++nj)
              if (j < m.mi_end[nj]) {
                if ((nj + 1 < n) && Q_sm[n][nj].type() == Matrix_TYPE::ZERO)
                  Q_sm[n][nj] = M_Matrix<double>(m.num_rows[n], m.num_rows[nj]);
                Q_sm[n][nj](i, j - m.mi_start[nj]) = q;
                Q_sm[n][n](i, i) -= q;
                break;
              }
        }
      }
    }
  }
}

inline std::size_t suggested_steps_number(const M_Matrix<double>& Q,
                                          double max_tol = 0.1) {
  double maxq = maxAbs(Matrix_Unary_Transformations::diag(Q));
  if (maxq<max_tol) return 1;
  else
  return std::max(1ul, 2ul << std::size_t(std::ceil(
                           std::log2(maxq / max_tol ))));
}

inline auto
calc_Q_for_taylor(std::size_t &nsteps, Microscopic_description_new const &start,
                  const M_Matrix<double> &Qm, const M_Matrix<double>& P0,
                  const std::vector<std::vector<std::size_t>> &connections_to_x,
                  std::size_t order_num, double reduce_p_by) {
  typedef myOptional_t<mi_Q> Op;
  mi_Q m(start);
  std::vector<std::vector<double>> Q_01;
  std::vector<std::vector<M_Matrix<double>>> Q_sm;

  std::size_t suggested_steps;
  m.prun = P0;
  std::size_t n = 0;
  while ((n < order_num) && (m.mi_current < m.mi.size())) {
    pre_Qmi(n, m, Q_sm);
    Q_01 = std::vector<std::vector<double>>(m.num_rows[n]);
    calc_Qmi(n, nsteps, connections_to_x, Qm, m, Q_sm, Q_01);
    double tol=0.1;
    suggested_steps = suggested_steps_number(Q_sm[n][n],tol);
    if (suggested_steps > nsteps) {
      nsteps = suggested_steps;
      return Op(false, "step_too_big");
    }
    //        if (Q01[0].empty())
    //            break;
    post_Qmi(n, reduce_p_by, m, Q_sm, Q_01);
    ++n;
  }
  m.Q_s = to_Matrix(Q_sm);
  m.P_s = M_Matrix<M_Matrix<double>>(1, m.Q_s.nrows());
  std::vector<M_Matrix<double>> Ods(m.Q_s.nrows());

  for (std::size_t n = 0; n < m.Q_s.nrows(); ++n) {
    if (n == 0) {
      Ods[n] = Matrix_Generators::eye<double>(m.num_rows[n]);
      m.P_s[n] = P0;
    } else {
      Ods[n] =
          M_Matrix<double>(m.num_rows[n], m.num_rows[n], Matrix_TYPE::ZERO);
      m.P_s[n] = M_Matrix<double>(1, m.num_rows[n], Matrix_TYPE::ZERO);
    }
  }
  m.Od_s = M_Matrix<M_Matrix<double>>(m.Q_s.nrows(),m.Q_s.nrows(),Matrix_TYPE::DIAGONAL,Ods);
  return Op(std::move(m));
}

inline auto calc_exmp_taylor_step(
    Microscopic_description_new const &start, const M_Matrix<double> &Qm,
    const M_Matrix<double> &gm, const M_Matrix<double> &P0,
    const std::vector<std::vector<std::size_t>> &connections_to_x,
    std::size_t order_num, double reduce_p_by, std::size_t nsteps) {

  auto Qo = calc_Q_for_taylor(nsteps, start, Qm, P0, connections_to_x,
                              order_num, reduce_p_by);
  while (!Qo) {
    Qo = calc_Q_for_taylor(nsteps, start, Qm, P0, connections_to_x, order_num,
                           reduce_p_by);
  }
  auto m = std::move(Qo).value();
  M_Matrix<M_Matrix<double>> xr = m.Od_s * m.Q_s;
  M_Matrix<M_Matrix<double>> expmQ = m.Od_s + xr;
  double a = 1.0;
  for (std::size_t n = 2; n + 1 < order_num; ++n) {
    a /= n;
    xr = xr * m.Q_s;
    expmQ = expmQ + xr * a;
  }


  std::vector<M_Matrix<double>> G_sv(m.Q_s.nrows());
  for (std::size_t n = 0; n < m.Q_s.nrows(); ++n) {
    G_sv[n] = M_Matrix<double>(m.num_rows[n], m.num_rows[n],
                              Matrix_TYPE::DIAGONAL, 0.0);
    for (std::size_t i = 0; i < m.num_rows[n]; ++i)
      for (std::size_t k = 0; k < m.mi.number_of_states(); ++k)
        G_sv[n](i, i) += m.mi.ns(i + m.mi_start[n])[k] * gm[k];
  }
  M_Matrix<M_Matrix<double>> G_s(m.Q_s.nrows(), m.Q_s.ncols(),
                                 Matrix_TYPE::DIAGONAL,G_sv);

  auto Gtot = calc_Gtot_taylor(m.Q_s, expmQ, G_s, 4) * (1.0 / nsteps);
  auto Gvar = calc_Gvar_taylor(m.Q_s, expmQ, G_s, 4) * (1.0 / nsteps / nsteps);
  M_Matrix<double> expQ = full(expmQ, {0ul, 1ul});
  M_Matrix<double> P1 = P0 * expQ;
  M_Matrix<double> Gt=full(Gtot,{0ul,1ul});
  M_Matrix<double> Gv=full(Gvar,{0ul,1ul});

  reduce_P_G_by_p(m.mi,P1,expQ,Gt,Gv,reduce_p_by);


  return std::make_tuple(nsteps, std::move(m.mi), std::move(P1),
                         std::move(expQ), std::move(Gt),
                         std::move(Gv));
}

inline auto
calc_exmp_taylor(Microscopic_description_new const &start,
                 const M_Matrix<double> &Qm, const M_Matrix<double> &gm,
                 const M_Matrix<double> &p0,
                 const std::vector<std::vector<std::size_t>> &connections_to_x,
                 std::size_t order_num, double reduce_p_by) {
  std::size_t min_steps = 1;
  auto [nsteps, mi, p, P, PG, PGG] = calc_exmp_taylor_step(
      start, Qm, gm, p0, connections_to_x, order_num, reduce_p_by, min_steps);
  double f = 1 / nsteps;
  while (f < 1) {
    min_steps = std::ceil(1 / (1 - f));
    auto [nsteps1, mi1, p1, P1, PG1, PGG1] = calc_exmp_taylor_step(
        mi, Qm, gm, p, connections_to_x, order_num, reduce_p_by, min_steps);
    PGG = (PGG * P1) + 2*(PG * PG1) + (P * PGG1);
    PG = (PG * P1) + (P * PG1);
    P = P * P1;
    p = p1;
    f = f + 1.0 / nsteps1;
    mi = mi1;
  }

  return std::make_tuple(f, p, P, mi, PG, PGG);
}

class SingleLigandModel {
public:
  SingleLigandModel() {}
  auto &Q0() const { return Q0_; }
  auto &Qa() const { return Qa_; }

  auto Q(double x) const { return Q0_ + Qa_ * x; }

  auto &g(double) const { return g_; }
  auto &myg() const { return g_; }

  auto Qx(double x, double tolerance) const {
    return Markov_Transition_rate::evaluate(Q(x), g(x), tolerance);
  }

  auto P(const Markov_Transition_rate &aQx, std::size_t n, double fs) const {
    return Markov_Transition_step(aQx, n, fs, min_P_);
  }

  auto Pdt(const Markov_Transition_rate &x, std::size_t n, double fs) const {
    return Markov_Transition_step_double(x, n, fs, min_P_);
  }

  std::size_t nstates() const { return Q0_.nrows(); }

  double noise_variance(std::size_t nsamples, double fs) const {
    return noise_variance_ * fs / nsamples;
  }

  double noise_variance_1Hz() const { return noise_variance_; }

  double AverageNumberOfChannels() const { return N_channels_; }

  std::size_t N_channels() const { return N_channels_; }

  SingleLigandModel(std::pair<M_Matrix<double>, M_Matrix<double>> &&Qs,
                    M_Matrix<double> &&g, double Vm, double N, double noise,
                    double min_P)
      : Q0_{std::move(Qs.first)}, Qa_{std::move(Qs.second)}, g_{g * Vm},
        N_channels_{N}, noise_variance_{noise}, min_P_{min_P} {}
  SingleLigandModel(const std::pair<M_Matrix<double>, M_Matrix<double>> &Qs,
                    const M_Matrix<double> &g, double Vm, double N,
                    double noise, double min_P)
      : Q0_{std::move(Qs.first)}, Qa_{std::move(Qs.second)}, g_{g * Vm},
        N_channels_{N}, noise_variance_{noise}, min_P_{min_P} {}

  SingleLigandModel(const M_Matrix<double> &Q0, const M_Matrix<double> &Qa,
                    const M_Matrix<double> &g, double N, double noise,
                    double min_P)
      : Q0_{Q0}, Qa_{Qa}, g_{g}, N_channels_{N},
        noise_variance_{noise}, min_P_{min_P} {}

  double min_P() const { return min_P_; }

  typedef SingleLigandModel self_type;

  constexpr static auto className = my_static_string("SingleLigandModel");
  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "Q0", &self_type::Q0),
        grammar::field(C<self_type>{}, "Qa", &self_type::Qa),
        grammar::field(C<self_type>{}, "g", &self_type::myg),
        grammar::field(C<self_type>{}, "AverageNumberOfChannels",
                       &self_type::AverageNumberOfChannels),
        grammar::field(C<self_type>{}, "plogL", &self_type::noise_variance_1Hz),
        grammar::field(C<self_type>{}, "eplogL", &self_type::min_P));
  }

private:
  M_Matrix<double> Q0_;
  M_Matrix<double> Qa_;
  M_Matrix<double> g_;
  double N_channels_;
  double noise_variance_;
  double min_P_;
};


template<template<class...> class Tr> class SingleLigandModel_new_;

typedef SingleLigandModel_new_<C> SingleLigandModel_new;


template<template<class...> class Tr>
class SingleLigandModel_new_ {
public:
    template<class T> using Tr_t=typename Tr<T>::type;

    SingleLigandModel_new_() {}
    auto &Q0() const { return Q0_; }
    auto &Qa() const { return Qa_; }

    auto Q(double x) const { return Q0_ + Qa_ * x; }

    auto &g(double) const { return g_; }
    auto &myg() const { return g_; }

    auto Qx(double x, double tolerance) const {
        return Markov_Transition_rate::evaluate(Q(x), g(x), tolerance);
    }

    auto P(const Markov_Transition_rate &aQx, std::size_t n, double fs) const {
        return Markov_Transition_step(aQx, n, fs, min_P_);
    }

    auto Pdt(const Markov_Transition_rate &x, std::size_t n, double fs) const {
        return Markov_Transition_step_double(x, n, fs, min_P_);
    }

    std::size_t nstates() const { return Q0_.nrows(); }

    Tr_t<double> noise_variance(std::size_t nsamples, double fs) const {
        return noise_variance_ * fs / nsamples;
    }

    Tr_t<double> noise_variance_1Hz() const { return noise_variance_; }

    Tr_t<double> AverageNumberOfChannels() const { return N_channels_; }

    std::size_t N_channels() const { return N_channels_; }

    SingleLigandModel_new_(std::pair<Tr_t<M_Matrix<double>>, Tr_t<M_Matrix<double>>> &&Qs,
                           Tr_t<M_Matrix<double>> &&g, double Vm, Tr_t<double> N, Tr_t<double> noise,
                      double min_P)
        : Q0_{std::move(Qs.first)}, Qa_{std::move(Qs.second)}, g_{g * Vm},
          N_channels_{N}, noise_variance_{noise}, min_P_{min_P} {}
    SingleLigandModel_new_(const std::pair<Tr_t<M_Matrix<double>>, Tr_t<M_Matrix<double>>> &Qs,
                           const Tr_t<M_Matrix<double>> &g, double Vm, Tr_t<double> N,
                           Tr_t<double> noise, double min_P)
        : Q0_{std::move(Qs.first)}, Qa_{std::move(Qs.second)}, g_{g * Vm},
          N_channels_{N}, noise_variance_{noise}, min_P_{min_P} {}

    SingleLigandModel_new_(const Tr_t<M_Matrix<double>> &Q0, const Tr_t<M_Matrix<double>> &Qa,
                           const Tr_t<M_Matrix<double>> &g, Tr_t<double> N, Tr_t<double> noise,
                      double min_P)
        : Q0_{Q0}, Qa_{Qa}, g_{g}, N_channels_{N},
          noise_variance_{noise}, min_P_{min_P} {}

    double min_P() const { return min_P_; }

    typedef SingleLigandModel_new_ self_type;

    constexpr static auto className = my_static_string("SingleLigandModel_new");
    static auto get_constructor_fields() {
        return std::make_tuple(
            grammar::field(C<self_type>{}, "Q0", &self_type::Q0),
            grammar::field(C<self_type>{}, "Qa", &self_type::Qa),
            grammar::field(C<self_type>{}, "g", &self_type::myg),
            grammar::field(C<self_type>{}, "AverageNumberOfChannels",
                           &self_type::AverageNumberOfChannels),
            grammar::field(C<self_type>{}, "plogL", &self_type::noise_variance_1Hz),
            grammar::field(C<self_type>{}, "eplogL", &self_type::min_P));
    }

private:
    Tr_t<M_Matrix<double>> Q0_;
    Tr_t<M_Matrix<double>> Qa_;
    Tr_t<M_Matrix<double>> g_;
    Tr_t<double> N_channels_;
    Tr_t<double> noise_variance_;
    double min_P_;
};





template <class Markov_Transition_step, class Markov_Transition_rate, class M,
          class Experiment, class X>
class Markov_Model_calculations {

public:
  myOptional_t<Markov_Transition_rate> calc_Qx(X x) const {
    return m_.Qx(x, tolerance());
  }
  myOptional_t<Markov_Transition_step> calc_sub_P(X x, std::size_t nsamples,
                                                  std::size_t n_sub_samples) {
    typedef myOptional_t<Markov_Transition_step> Op;
    auto Qx = calc_Qx(x);
    if (!Qx)
      return Op(false, "fails to calculate Qx: " + Qx.error());
    else
      return Op(Markov_Transition_step(n_sub_samples, Qx.value(), nsamples,
                                       fs(), Model().min_P()));
  }

  myOptional_t<Markov_Transition_step> calc_P(X x, std::size_t nsamples,
                                              std::size_t num_recursion = 0) {
    //typedef myOptional_t<Markov_Transition_step> Op;
    if (num_recursion == 0) {
      return Markov_Transition_step(num_recursion, Model().Q(x), Model().g(x),
                                    nsamples, fs(), Model().min_P() * 0.0001);

    } else {
      auto Qx = calc_Qx(x);
      if (!Qx) {
        // eigenvalues are too close or something like that, try another
        // approach.

        return Markov_Transition_step(2, Model().Q(x), Model().g(x), nsamples,
                                      fs(), Model().min_P() * 0.0001);
      } else {

        if constexpr (std::is_same_v<Markov_Transition_step,
                                     Markov_Transition_step_double>) {
          Markov_Transition_step out(Qx.value(), nsamples, fs(),
                                     Model().min_P());
          if constexpr (false) {
            std::random_device rd;
            auto initseed = rd();
            std::mt19937_64 mt(initseed);
            Markov_Transition_step out2(mt, Qx.value(), nsamples, fs(),
                                        Model().min_P(), 2000, 40000);
            are_Equal<true, Markov_Transition_step_double>().test(out, out2,
                                                                  std::cerr);
          }
          assert(([&out, &Qx, nsamples, this]() {
            auto out2 = Markov_Transition_step(2, Qx.value().Qrun(),
                                               Qx.value().g(), nsamples, fs(),
                                               Model().min_P() * 0.0001);
            auto out01 = Markov_Transition_step(1, Qx.value().Qrun(),
                                                Qx.value().g(), nsamples, fs(),
                                                Model().min_P() * 0.0001);
            auto out4 = Markov_Transition_step(4, Qx.value().Qrun(),
                                               Qx.value().g(), nsamples, fs(),
                                               Model().min_P() * 0.0001);
            auto out8 = Markov_Transition_step(8, Qx.value().Qrun(),
                                               Qx.value().g(), nsamples, fs(),
                                               Model().min_P() * 0.0001);
            assert((are_Equal<true, Markov_Transition_step_double>().test(
                out, out2, 1e-7, std::cerr)));
            return true;
          }()));

          return out;

        } else
          return Markov_Transition_step(Qx.value(), nsamples, fs(),
                                        Model().min_P());
        //     return Markov_Transition_step(get_Qx(x),nsamples,fs());
      }
    }
  }
  auto Peq(const X &x) {

    typedef myOptional_t<std::invoke_result_t<
        decltype(&Markov_Transition_rate::calc_Peq), Markov_Transition_rate>>
        Op;
    //  auto Qx=get_Qx(x);
    auto Qx = calc_Qx(x);
    // std::cerr<<"calcQx"<<Qx.value().Qrun();
    if (!Qx)
      return Op(false, "cannot calculate Qx: " + Qx.error());
    else
      return Op(Qx.value().calc_Peq());
  };

  auto AverageNumberOfChannels() const { return m_.AverageNumberOfChannels(); }

  auto noise_variance(std::size_t nsamples) const {
    return m_.noise_variance(nsamples, fs());
  }

  auto g(const X &x) const { return m_.g(x); }

  std::size_t nstates() const { return m_.nstates(); }

  const auto &get_Qx(X x) {
    auto it = rate_map_.find(x);
    if (it != rate_map_.end())
      return it->second;
    else {
      rate_map_[x] = calc_Qx(x);
      return get_Qx(x);
    }
  }

  Markov_Transition_step const &get_P(X x, std::size_t nsamples) {
    if ((x == current_x_) && (nsamples == current_nsamples_))
      return *current_P_;

    else {
      auto it = step_map_.find(x);
      if (it == step_map_.end()) {
        current_x_ = x;
        current_nsamples_ = nsamples;

        rate_map_[x] = calc_Qx(x);
        step_map_[x][nsamples] = calc_P(x, nsamples);
        current_P_ = &step_map_[x][nsamples];
        return *current_P_;
      } else {
        auto it2 = it->second.find(nsamples);
        if (it2 != it->second.end()) {
          current_x_ = x;
          current_nsamples_ = nsamples;
          current_P_ = &(it2->second);
          return *current_P_;
        } else {
          current_x_ = x;
          current_nsamples_ = nsamples;
          step_map_[x][nsamples] = calc_P(x, nsamples);
          current_P_ = &step_map_[x][nsamples];
          return *current_P_;
        }
      }
    }
  }

  template <typename step>
  //    Markov_Transition_step const& get_P(const step& s)
  auto get_P(const step &s, std::size_t num_recursion = 0) {
    typedef myOptional_t<Markov_Transition_step> Op;
    auto b = s.begin();
    auto e = s.end();
    if ((e - b) > 1) {
      auto PP = calc_P(b->x(), b->nsamples());
      if (!PP.has_value())
        return Op(false, "fails at start :" + PP.error());
      // auto& P=get_P(b->x(),b->nsamples());
      //   auto PP=P.getMinimum();
      ++b;
      for (auto it = b; it != e; ++it) {
        auto pc = calc_P(it->x(), it->nsamples());
        if (!pc)
          return Op(false,
                    "fails at interval :" + ToString(*it) + " : " + pc.error());

        else
          PP.value() *= std::move(pc).value();

        //    PP*=get_P(it->x(), it->nsamples());
      }
      // current_step_.reset(new Markov_Transition_step(std::move(PP),fs()));
      // return *current_step_;
      return PP;
    } else
      return calc_P(b->x(), b->nsamples(), num_recursion);
    //        return get_P(b->x(), b->nsamples());
  }

  template <typename step>
  auto get_Pmi(const Microscopic_description &mi, const step &s,
               std::size_t num_recursion = 0) {

      typedef myOptional_t<Microscopic_Transition_step_double> Op;
    auto Qs = mi.Qs(Model().Q0(), Model().Qa(), connections_to_0(),
                    connections_to_a());
    auto g = mi.g(Model().g(0.0));

    SingleLigandModel SM(std::move(Qs), std::move(g), 1.0,
                         AverageNumberOfChannels(), noise_variance(1), min_P());

    Markov_Model_calculations<Markov_Transition_step, Markov_Transition_rate,
                              SingleLigandModel, Experiment, X>
        MiCalc(SM, get_Experiment(), n_sub_samples(), tolerance());
    auto res=MiCalc.get_P(s, num_recursion);
  //  if (! res)
        return Op(false, res.error());
  //  else
  //          return Op(Microscopic_Transition_step_double(res.value()));
  }

  template <typename step>
  auto get_Pmi_new(const Microscopic_description_new &mi0,
                   const M_Matrix<double>& p0, const step &s,
                   std::size_t order_num, double reduce_by_p) {

    typedef myOptional_t<Microscopic_Transition_step_double> Op;
    double x = s.x().x();
    double dt=s.x().nsamples()/fs();
    auto Qm = Model().Q(x)*dt;
    auto g = Model().g(x);
    auto connections_x = connections_to_x(x);
    auto [f, p1, P, mi1, PG, PGG] =
        calc_exmp_taylor(mi0, Qm, g, p0, connections_x, order_num, reduce_by_p);
    assert(f == 1);
    return Op(Microscopic_Transition_step_double(
        mi0, std::move(mi1), std::move(P), std::move(PG), std::move(PGG)));
  }

  const M &Model() const { return m_; }

  double min_P() const { return Model().min_P(); }
  double tolerance() const { return tolerance_; }
  double n_sub_samples() const { return n_sub_samples_; }
  std::size_t N(const X &) const { return N_; }
  double fs() const { return fs_; }
  auto &get_Experiment() const { return e_; }

  auto &connections_to_0() const { return connections_to_0_; }
  auto &connections_to_a() const { return connections_to_a_; }
  auto &connections_to_x(double x) const {
    if (x > 0)
      return connections_to_x_;
    else
      return connections_to_0_;
  }
  auto &connections_from_x(double x) const {
    if (x > 0)
      return connections_from_x_;
    else
      return connections_from_0_;
  }
  Markov_Model_calculations(
      const M &m, std::vector<std::vector<std::size_t>> connections_to_0,
      std::vector<std::vector<std::size_t>> connections_to_a,
      std::vector<std::vector<std::size_t>> connections_to_x,
      std::vector<std::vector<std::size_t>> connections_from_0,
      std::vector<std::vector<std::size_t>> connections_from_x,
      const Experiment &e, std::size_t n_sub_samples, double tolerance)
      : m_{m}, connections_to_0_{connections_to_0},
        connections_to_a_{connections_to_a},
        connections_to_x_{connections_to_x},
        connections_from_0_{connections_from_0},
        connections_from_x_{connections_from_x},

        N_{std::size_t(m.N_channels())}, fs_{e.frequency_of_sampling()},
        n_sub_samples_{n_sub_samples}, tolerance_{tolerance}, e_{e},
        current_x_{}, current_step_{}, rate_map_{}, map_{}, step_map_{} {}
  Markov_Model_calculations(const M &m, const Experiment &e,
                            std::size_t n_sub_samples, double tolerance)
      : m_{m}, connections_to_0_{}, connections_to_a_{}, connections_to_x_{},
        connections_from_0_{}, connections_from_x_{},

        N_{std::size_t(m.N_channels())}, fs_{e.frequency_of_sampling()},
        n_sub_samples_{n_sub_samples}, tolerance_{tolerance}, e_{e},
        current_x_{}, current_step_{}, rate_map_{}, map_{}, step_map_{} {}

private:
  M m_;
  std::vector<std::vector<std::size_t>> connections_to_0_;
  std::vector<std::vector<std::size_t>> connections_to_a_;
  std::vector<std::vector<std::size_t>> connections_to_x_;
  std::vector<std::vector<std::size_t>> connections_from_0_;
  std::vector<std::vector<std::size_t>> connections_from_x_;
  size_t N_;
  double fs_;
  std::size_t n_sub_samples_;
  double tolerance_;
  Experiment e_;
  mutable X current_x_;
  mutable std::size_t current_nsamples_;
  Markov_Transition_step *current_P_;

  std::unique_ptr<Markov_Transition_step> current_step_;

  mutable std::map<X, Markov_Transition_rate> rate_map_;

  std::map<X, std::map<std::size_t, std::size_t>> map_;

  std::map<X, std::map<std::size_t, Markov_Transition_step>> step_map_;
};

template <bool recursive, int averaging, bool variance> class Macro_R_new;




#endif // QMODEL_H



