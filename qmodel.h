#ifndef QMODEL_H
#define QMODEL_H

#include "Matrix.h"
#include "myfields.h"
#include "mymath.h"
#include "mytypetraits.h"
#include "myparameters.h"
#include "myDistributions.h"
#include "mytests.h"
#include "mydynamicfunctions.h"
#include <vector>
#include <map>
#include <set>
#include <cmath>

class Model_Parameter_label : public Parameter_label {
public:
  using Parameter_label::Parameter_label;
  using Label::legal_chars;

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
  std::map<std::pair<std::size_t, std::size_t>, std::unique_ptr<Base_Function<double,State_Model>>> Q0s_;
  std::map<std::pair<std::size_t, std::size_t>, std::unique_ptr<Base_Function<double,State_Model>>> Qas_;
  std::map<std::size_t, std::unique_ptr<Base_Function<double,State_Model>>> gs_;
  std::string errors_;
  bool isValid_;


public:
  bool isValid()const {return isValid_;}
  std::string error()const { return errors_;}
  std::size_t k() const { return numstates_; }
  auto &transition_rates() const { return Q0s_; }
  std::map<std::pair<std::size_t,std::size_t>,std::string> transition_rates_text() const { return text_map<double,State_Model>(Q0s_); }

  auto &agonist_transitions_rates() const { return Qas_; }

  std::map<std::pair<std::size_t,std::size_t>,std::string> agonist_transitions_rates_text() const { return text_map<double,State_Model>(Qas_); }
  auto &conductances() const { return gs_; }
  std::map<std::size_t,std::string> conductances_text() const { return text_map(gs_); }

  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "number_of_states", &self_type::k),
        grammar::field(C<self_type>{}, "transition_rates", &self_type::transition_rates_text),
        grammar::field(C<self_type>{}, "agonist_transitions_rates",&self_type::agonist_transitions_rates_text),
        grammar::field(C<self_type>{}, "conductances", &self_type::conductances_text));
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
           std::map<std::size_t, std::string> myconductances)
  {
   typedef myOptional_t<State_Model> Op;
    auto transition_rates_r=compile_map<double,State_Model>(mytransition_rates);
    auto agonist_transitions_r=compile_map<double,State_Model>(myagonist_transitions);
    auto conductances_r=compile_map<double,State_Model>(myconductances);

    if(transition_rates_r.has_value()&&agonist_transitions_r.has_value()&&conductances_r.has_value())
      return Op(State_Model(number_of_states,std::move(transition_rates_r).value(),std::move(agonist_transitions_r).value(),std::move(conductances_r).value()));
    else return Op(false,"error building state model: "+transition_rates_r.error()+agonist_transitions_r.error()+conductances_r.error());
  }

      State_Model() = default;

      State_Model(const State_Model& other):numstates_{other.numstates_}, Q0s_{clone_map(other.Q0s_)},
                                              Qas_{clone_map(other.Qas_)}, gs_{clone_map(other.gs_)},errors_{other.errors_}, isValid_(other.isValid_) {}

      State_Model( State_Model &&other)
          : numstates_{other.numstates_}, Q0s_{std::move(other.Q0s_)},
            Qas_{std::move(other.Qas_)}, gs_{std::move(other.gs_)},
            errors_{std::move(other.errors_)}, isValid_{other.isValid_} {}

      State_Model& operator=(const State_Model& other)
      {
        State_Model tmp(other);
        *this=std::move(tmp);
        return *this;
      }
      State_Model &operator=(State_Model &&other) {
        numstates_=other.numstates_;
        Q0s_=std::move(other.Q0s_);
        Qas_=std::move(other.Qas_);
        gs_=std::move(other.gs_);
        errors_=std::move(other.errors_);
        isValid_=std::move(other.isValid_);
        return *this;
      }

      State_Model(std::size_t number_of_states,
               const std::map<std::pair<std::size_t, std::size_t>, std::string>&
                   mytransition_rates,
               const std::map<std::pair<std::size_t, std::size_t>, std::string>&
                   myagonist_transitions,
                 const  std::map<std::size_t, std::string>& myconductances)
      {
        auto transition_rates_r =
            compile_map<double, State_Model>(mytransition_rates);
        auto agonist_transitions_r =
            compile_map<double, State_Model>(myagonist_transitions);
        auto conductances_r = compile_map<double, State_Model>(myconductances);

        if (transition_rates_r.has_value() &&
            agonist_transitions_r.has_value() && conductances_r.has_value())
        {
          numstates_=number_of_states;
          Q0s_=make_unique_map(std::move(transition_rates_r).value());
          Qas_=make_unique_map(std::move(agonist_transitions_r).value());
          gs_=make_unique_map(std::move(conductances_r).value());
          errors_.clear();
          isValid_=true;

        } else {
          errors_=transition_rates_r.error()+agonist_transitions_r.error()+conductances_r.error();
          isValid_=false;
        }
      }

  State_Model(
      std::size_t number_of_states,
      std::map<std::pair<std::size_t, std::size_t>, Base_Function<double,State_Model>*>&&
          mytransition_rates,
      std::map<std::pair<std::size_t, std::size_t>, Base_Function<double,State_Model>*>&&
          myagonist_transitions,
      std::map<std::size_t, Base_Function<double,State_Model>*>&& myconductances)
          : numstates_{number_of_states}, Q0s_{make_unique_map(std::move(mytransition_rates))},
            Qas_{make_unique_map(std::move(myagonist_transitions))}, gs_{make_unique_map(std::move(myconductances))},errors_{}, isValid_(true) {}




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

class Markov_Transition_rate {

  M_Matrix<double> Qrun_; // transition rate matrix at time zero
  M_Matrix<double> g_;

  M_Matrix<double> V_;     // eigenvector of Qrun
  M_Matrix<double> W_;     // eigenvector of Qrun
  M_Matrix<double> landa_; // eigenvalues

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
    std::tie(V_, landa_, W_) = std::move(eig);
    clean_landa(landa_);
    Wg_ = W_ * g_;
    WgV_ = W_ * Matrix_Unary_Transformations::diag(g_) * V_;
    // assert(Frobenius_test::test(Qrun(),V()*landa()*W(),1,
    // std::sqrt(std::numeric_limits<double>::epsilon())));
  }

public:
  struct landa : public invariant {
    static Op_void test(const M_Matrix<double> landa, double tol) {
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

Markov_Transition_step_single_minimum::Markov_Transition_step_single_minimum(
    const Markov_Transition_step_single &x) {
  n = x.nsamples();
  PPn = x.P();
  PG_n = x.gmean_i() * n;
  PG_n = x.gtotal_ij() * n;
  PGG_n = x.gsqr_i() * (n * n);
}

Markov_Transition_step_single_minimum &Markov_Transition_step_single_minimum::
operator*=(const Markov_Transition_step_single &x) {
  auto n1 = x.nsamples();
  PGn += ((PPn * x.gmean_i()) * n1);
  PGG_n += ((PPn * x.gsqr_i()) * (n1 * n1) + (PG_n * x.gmean_i()) * n1);
  PG_n = (PG_n * x.P()) + (PPn * x.gtotal_ij()) * n1;
  PPn = PPn * x.P();
  n += n1;
  return *this;
}

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
  Markov_Transition_step_double_minimum
  operator*=(const Markov_Transition_step_double &x_);
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
  /// variance of the mean conductance for each starting state i contributed by
  /// the ones ending at state j
  M_Matrix<double> const &gtotal_var_ij() const {
    return gtotal_var_ij_;
  } // variance of the conductance matrix

  /// variance of the mean conductance for each starting state i and ending
  /// state j
  M_Matrix<double> const &gvar_ij() const {
    return gvar_ij_;
  } // variance of the conductance matrix

  Markov_Transition_step_double(const Markov_Transition_rate &Qx,
                                std::size_t nsamples, double fs, double min_p)
      : Markov_Transition_step(Qx, nsamples, fs, min_p) {
    init(Qx);
  }

  Markov_Transition_step_double(std::size_t recursion_number,
                                const Markov_Transition_rate &Qx,
                                std::size_t nsamples, double fs, double min_p)
      : Markov_Transition_step(Qx, nsamples, fs, min_p) {
    init_by_step(Qx, recursion_number, fs);
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
  };

  static double EX_111(double x, double y, double z, double exp_x) {
    return exp_x / ((x - y) * (x - z));
  }

  static double E111(double x, double y, double z, double exp_x, double exp_y,
                     double exp_z) {
    return EX_111(x, y, z, exp_x) + EX_111(y, x, z, exp_y) +
           EX_111(z, y, x, exp_z);
  }
  static double E12(double x, double y, double exp_x, double exp_y) {
    return EX_111(x, y, y, exp_x) + exp_y / (y - x) * (1.0 - 1.0 / (y - x));
  }

  static double E3(double x, double y, double z, double exp_x, double exp_y,
                   double exp_z,
                   double eps = std::numeric_limits<double>::epsilon()) {
    if (sqr(x - y) < eps) // x==y
    {
      if (sqr(y - z) < eps) // y==z
        return exp_x / 2.0; // x==y==z
      else
        return E12(z, x, exp_z, exp_x); // x==y!=z
    } else if (sqr(y - z) < eps)        // x!=y==z
    {
      return E12(x, y, exp_x, exp_y);
    } else if (sqr(x - z) < eps) // y!=z==x!=y
    {
      return E12(y, x, exp_y, exp_x);
    } else
      return E111(x, y, z, exp_x, exp_y, exp_z); // x!=y!=z!=x
  }

private:
  void init(Markov_Transition_step_double_minimum &&x, double) {
    M_Matrix<double> u = ones<double>(P().ncols(), 1);

    gtotal_sqr_ij_ = x.PGG_n / (x.n * x.n * 0.5);
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

  void init(const Markov_Transition_rate &Qx) {
    // const double eps=std::numeric_limits<double>::epsilon();
    double dt = Markov_Transition_step::dt();

    std::size_t N = P().ncols();

    auto ladt = Qx.landa() * dt;

    auto exp_ladt = ladt.apply([](double x) { return std::exp(x); });

    M_Matrix<double> E2m(N, N, Matrix_TYPE::SYMMETRIC);
    M_Matrix<double> E2mb(N, N, Matrix_TYPE::SYMMETRIC);
    for (std::size_t i = 0; i < N; ++i)
      for (std::size_t j = 0; j < i + 1; ++j)
        E2m(i, j) = Ee(ladt[i], ladt[j], exp_ladt[i], exp_ladt[j], min_P());

    // build E2
    M_Matrix<double> WgV_E2(N, N);

    for (std::size_t i = 0; i < N; ++i)
      for (std::size_t j = 0; j < N; ++j)
        WgV_E2(i, j) = Qx.WgV()(i, j) * E2m(i, j);

    gtotal_ij_ = Qx.V() * WgV_E2 * Qx.W();

    M_Matrix<double> WgV_E3(N, N, 0.0);
    for (std::size_t n1 = 0; n1 < N; n1++)
      for (std::size_t n3 = 0; n3 < N; n3++)
        for (std::size_t n2 = 0; n2 < N; n2++) {
          WgV_E3(n1, n3) +=
              Qx.WgV()(n1, n2) * Qx.WgV()(n2, n3) *
              E3(ladt[n1], ladt[n2], ladt[n3], exp_ladt[n1], exp_ladt[n2],
                 exp_ladt[n3], min_P()); // optimizable
        }

    gtotal_sqr_ij_ = Qx.V() * WgV_E3 * Qx.W() * 2.0;
    for (std::size_t i = 0; i < N; ++i)
      for (std::size_t j = 0; j < N; ++j)
        if (P()(i, j) == 0) {
          gtotal_ij_(i, j) = 0;
          gtotal_sqr_ij_(i, j) = 0;
        }

    M_Matrix<double> U = Matrix_Generators::ones<double>(1, g().size());
    M_Matrix<double> UU = M_Matrix<double>(g().size(), g().size(), 1);
    auto gmean_ij_p = TransposeSum(g() * U) * (0.5);
    auto gvar_ij_p = (g() * U - Transpose(g() * U)).apply([](double x) {
      return std::abs(x);
    }) * (0.5);

    // std::cerr<<"\ngmean_ij_p=\n"<<gmean_ij_p<<"\ngvar_ij_p=\n"<<gvar_ij_p<<"\n";
    // std::cerr<<"\n UU=ņ"<<UU<<"\n";
    auto gmean_ij_tot = gtotal_ij_ + gmean_ij_p * min_P();
    auto P_p = P() + UU * min_P();
    gmean_ij_ = elemDiv(gmean_ij_tot, P_p);
    gtotal_var_ij_ = gtotal_sqr_ij_ - elemMult(gtotal_ij_, gmean_ij_);
    auto gvar_ij_tot = gtotal_var_ij_ + gvar_ij_p * min_P();
    gvar_ij_ = elemDiv(gvar_ij_tot, P_p);
    M_Matrix<double> u(N, 1, 1.0);
    gmean_i_ = gtotal_ij_ * u;
    gsqr_i_ = gtotal_sqr_ij_ * u;
    gvar_i_ = gtotal_var_ij_ * u;
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

  void init_by_step(const Markov_Transition_rate &Qx,
                    std::size_t recursion_number, double fs) {
    std::size_t times = std::pow(2, recursion_number);
    std::size_t n = nsamples();
    Markov_Transition_step_double s(Qx, n, fs * times, min_P());
    Markov_Transition_step_double step_ini(s);
    for (std::size_t i = 0; i < recursion_number; ++i) {
      Markov_Transition_step_double_minimum step(step_ini);
      step *= step_ini;
      step_ini = Markov_Transition_step_double(std::move(step), step_ini.g(),
                                               fs, min_P());
    }
    *this = step_ini;
  }
};

template <bool output> class are_Equal<output, Markov_Transition_step_double> {
public:
  template <class ostream>
  bool test(const Markov_Transition_step_double &one,
            const Markov_Transition_step_double &other,
            ostream &os = std::cerr) const {
    are_Equal<output, M_Matrix<double>> mytest(one.min_P(), one.min_P());
    bool P_equal = mytest.test_prod(one.P(), other.P(), os);
    if constexpr (output)
      if (!P_equal)
        os << "\n differ in  P_ij !!  \n";
    bool gtotal_equal =
        mytest.test_prod(one.gtotal_ij(), other.gtotal_ij(), os);
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
};

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

Markov_Transition_step_double_minimum::Markov_Transition_step_double_minimum(
    const Markov_Transition_step_double &x)
    : n{x.nsamples()}, PPn{x.P()}, PG_n{x.gtotal_ij() * x.nsamples()},
      PGG_n{x.gtotal_sqr_ij() * (x.nsamples() * x.nsamples() * 0.5)},
      min_p(x.min_P()) {}

Markov_Transition_step_double_minimum Markov_Transition_step_double_minimum::
operator*=(const Markov_Transition_step_double &x) {
  auto n1 = x.nsamples();
  PGG_n = (PGG_n * x.P()) + (PG_n * x.gtotal_ij()) * n1 +
          (PPn * x.gtotal_sqr_ij()) * (0.5 * n1 * n1);
  PG_n = (PG_n * x.P()) + (PPn * x.gtotal_ij()) * n1;
  PPn = Probability_transition::normalize(PPn * x.P(), min_p);
  n = n + n1;
  return *this;
}

class SingleLigandModel {
public:
  SingleLigandModel(){}
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

  double noise_variance_1Hz() const {
    return noise_variance_ ;
  }


  double AverageNumberOfChannels() const { return N_channels_; }

  std::size_t N_channels() const { return N_channels_; }



  SingleLigandModel(std::pair<M_Matrix<double>, M_Matrix<double>> &&Qs,
                    M_Matrix<double> &&g, double Vm, double N, double noise,
                    double min_P)
      : Q0_{std::move(Qs.first)}, Qa_{std::move(Qs.second)}, g_{g * Vm},
        N_channels_{N}, noise_variance_{noise}, min_P_{min_P} {}
  SingleLigandModel(const std::pair<M_Matrix<double>, M_Matrix<double>> &Qs,
                    const M_Matrix<double> &g, double Vm, double N, double noise,
                    double min_P)
      : Q0_{std::move(Qs.first)}, Qa_{std::move(Qs.second)}, g_{g * Vm},
        N_channels_{N}, noise_variance_{noise}, min_P_{min_P} {}


  SingleLigandModel(const M_Matrix<double>& Q0, const M_Matrix<double>& Qa,
                    const M_Matrix<double> &g,  double N, double noise,
                    double min_P)
      : Q0_{Q0}, Qa_{Qa}, g_{g },
        N_channels_{N}, noise_variance_{noise}, min_P_{min_P} {}

  double min_P() const { return min_P_; }

  typedef SingleLigandModel self_type;

  constexpr static auto className = my_static_string("SingleLigandModel");
  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "Q0", &self_type::Q0),
        grammar::field(C<self_type>{}, "Qa", &self_type::Qa),
        grammar::field(C<self_type>{}, "g", &self_type::myg),
        grammar::field(C<self_type>{}, "AverageNumberOfChannels", &self_type::AverageNumberOfChannels),
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
    typedef myOptional_t<Markov_Transition_step> Op;
    auto Qx = calc_Qx(x);
    if (!Qx)
      return Op(false, "fails to calculate Qx: " + Qx.error());
    else {
      if constexpr (std::is_same_v<Markov_Transition_step,
                                   Markov_Transition_step_double>) {
        if (num_recursion == 0) {

          Markov_Transition_step out(Qx.value(), nsamples, fs(),
                                     Model().min_P());
          if (0) {
            std::random_device rd;
            auto initseed = rd();
            std::mt19937_64 mt(initseed);
            Markov_Transition_step out2(mt, Qx.value(), nsamples, fs(),
                                        Model().min_P(), 2000, 40000);
            are_Equal<true, Markov_Transition_step_double>().test(out, out2,
                                                                  std::cerr);
          }
          return out;

        } else
          return Markov_Transition_step(num_recursion, Qx.value(), nsamples,
                                        fs(), Model().min_P());
      } else
        return Markov_Transition_step(Qx.value(), nsamples, fs(),
                                      Model().min_P());
      //     return Markov_Transition_step(get_Qx(x),nsamples,fs());
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

  const M &Model() const { return m_; }

  double min_P() const { return Model().min_P(); }
  double tolerance() const { return tolerance_; }
  double n_sub_samples() const { return n_sub_samples_; }
  std::size_t N(const X &) const { return N_; }
  double fs() const { return fs_; }
  auto& get_Experiment()const { return e_;}

  Markov_Model_calculations(const M &m, const Experiment &e,
                            std::size_t n_sub_samples, double tolerance)
      : m_{m}, N_{std::size_t(m.N_channels())}, fs_{e.frequency_of_sampling()},
        n_sub_samples_{n_sub_samples},
        tolerance_{tolerance}, e_{e}, current_x_{},current_step_{},rate_map_{}, map_{}, step_map_{} {}

private:
  M m_;
  size_t N_;
  double fs_;
  std::size_t n_sub_samples_;
  double tolerance_;
  const Experiment &e_;
  mutable X current_x_;
  mutable std::size_t current_nsamples_;
  Markov_Transition_step *current_P_;

  std::unique_ptr<Markov_Transition_step> current_step_;

  mutable std::map<X, Markov_Transition_rate> rate_map_;

  std::map<X, std::map<std::size_t, std::size_t>> map_;

  std::map<X, std::map<std::size_t, Markov_Transition_step>> step_map_;
};

#endif // QMODEL_H
