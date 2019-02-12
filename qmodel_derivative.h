#ifndef QMODEL_DERIVATIVE_H
#define QMODEL_DERIVATIVE_H
#include "qmodel.h"
#include "matrixderivative.h"
#include "myparameters_derivative.h"

template <> class Derivative<State_Model> : public State_Model {
public:
  typedef Derivative self_type;
  typedef Model_Parameter_label myParameter_label;
  typedef State_Model base_type;
  using base_type::k;
  using base_type::transition_rates_text;
  using base_type::agonist_transitions_rates_text;
  using base_type::conductances_text;

  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "number_of_states", &self_type::k),
        grammar::field(C<self_type>{}, "transition_rates", &self_type::transition_rates_text),
        grammar::field(C<self_type>{}, "agonist_transitions_rates",&self_type::agonist_transitions_rates_text),
        grammar::field(C<self_type>{}, "conductances", &self_type::conductances_text));
  }
  Derivative(std::size_t number_of_states,
               const std::map<std::pair<std::size_t, std::size_t>, std::string>&
                   mytransition_rates,
                        const std::map<std::pair<std::size_t, std::size_t>, std::string>&
                   myagonist_transitions,
                 const  std::map<std::size_t, std::string>& myconductances)

      : base_type(number_of_states,mytransition_rates,myagonist_transitions,myconductances){}


public:
  Derivative() = default;
  Derivative(const base_type &model) : base_type{model} {}

  Derivative(
      std::size_t number_of_states,
      std::map<std::pair<std::size_t, std::size_t>,
              Base_Function<double, State_Model>*>
          mytransition_rates,
      std::map<std::pair<std::size_t, std::size_t>,
               Base_Function<double, State_Model>*>
          myagonist_transitions,
      std::map<std::size_t, Base_Function<double, State_Model>*>
          myconductances)
      : base_type(number_of_states, std::move(mytransition_rates), std::move(myagonist_transitions),
                  std::move(myconductances)) {}

  template <class Parameters> auto Qs(const Derivative<Parameters> &p) const {

    M_Matrix<Derivative<double>> Q0(k(), k(), Matrix_TYPE::FULL,
                                    Derivative<double>(p.x()));
    M_Matrix<Derivative<double>> Qa(k(), k(), Matrix_TYPE::FULL,
                                    Derivative<double>(p.x()));

    for (auto &e : transition_rates()) {
      Q0(e.first.first, e.first.second) = (*e.second)(p);
      Q0(e.first.first, e.first.first) -= Q0(e.first.first, e.first.second);
    }
    for (auto &e : agonist_transitions_rates()) {
      Qa(e.first.first, e.first.second) = (*e.second)(p);
      Qa(e.first.first, e.first.first) -= Qa(e.first.first, e.first.second);
    }
    //    auto r=Q0*ones<Constant<double>>(Q0.ncols(),1);
    //    auto q=Qa*ones<Constant<double>>(Q0.ncols(),1);

    return std::pair(Derivative<M_Matrix<double>>(std::move(Q0)),
                     Derivative<M_Matrix<double>>(std::move(Qa)));
  }

  template <class Parameters> auto g(const Derivative<Parameters> &p) const {
    M_Matrix<Derivative<double>> out(k(), 1, Derivative<double>(p.x()));
    for (auto &e : conductances())
      out[e.first] = (*e.second)(p);
    return Derivative<M_Matrix<double>>(out);
  }
};

template <class F, class Model, class Parameters>
auto Incremental_ratio(double eps, const F &fun, const Model &M,
                       const Derivative<Parameters_values<State_Model>> &p) {
  auto f = std::invoke(fun, M, p.f());
  M_Matrix<M_Matrix<double>> dfdx(p.x().nrows(), p.x().ncols(), p.x().type());

  for (std::size_t i = 0; i < p.x().size(); ++i) {
    double e = eps;
    if (p.x()[i] != 0)
      e *= std::abs(p.x()[i]);
    auto ppos = p.f_dfdx(i, e);
    auto pneg = p.f_dfdx(i, -e);

    auto fpos = fun(ppos);
    auto fneg = fun(pneg);
    dfdx[i] = (fpos - fneg) / (2.0 * e);
  }
  return Derivative<M_Matrix<double>>(f, p.x(), dfdx);
}

template <class F, class Model, class Parameters>
auto Incremental_ratio_pair(
    double eps, const F &fun, const Model &M,
    const Derivative<Parameters_values<State_Model>> &p) {
  auto [f1, f2] = std::invoke(fun, M, p.f());
  M_Matrix<M_Matrix<double>> df1dx(p.x().nrows(), p.x().ncols(), p.x().type());
  M_Matrix<M_Matrix<double>> df2dx(p.x().nrows(), p.x().ncols(), p.x().type());

  for (std::size_t i = 0; i < p.x().size(); ++i) {
    double e = eps;
    if (p.x()[i] != 0)
      e *= std::abs(p.x()[i]);
    auto ppos = p.f_dfdx(i, e);
    auto pneg = p.f_dfdx(i, -e);

    auto [fpos1, fpos2] = fun(ppos);
    auto [fneg1, fneg2] = fun(pneg);
    df1dx[i] = (fpos1 - fneg1) / (2.0 * e);
    df2dx[i] = (fpos2 - fneg2) / (2.0 * e);
  }
  return std::pair(Derivative<M_Matrix<double>>(f1, p.x(), df1dx),
                   Derivative<M_Matrix<double>>(f2, p.x(), df2dx));
}

template <> class Derivative<Allosteric_Model>: public Allosteric_Model
{
public:
  typedef Derivative self_type;
  typedef Model_Parameter_label myParameter_label;
  typedef Allosteric_Model base_type;
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

public:
  Derivative() = default;
  Derivative(const base_type &model) : base_type{model} {}

  Derivative(
      std::size_t number_of_units,
      std::map<Conformational_change_label, Conformational_change>
          conformational_changes,
      std::vector<Conformational_change_label> unit_of_conformational_changes,
      std::set<Conformational_interaction> conformational_interactions,
      std::map<std::size_t, Conductance_Parameter_label> conductance_names)
      : base_type(number_of_units, conformational_changes,
                  unit_of_conformational_changes, conformational_interactions,
                  conductance_names) {}

  Derivative(const std::vector<std::string> &conformational_changes,
             const std::map<std::pair<std::string, bool>, std::string>
                 conformational_changes_names,
             const std::set<std::size_t> &agonist_changes,
             const std::set<std::size_t> &conductance_changes,
             const std::map<std::size_t, std::string> &conductance_names,
             const std::multimap<std::size_t,
                                 std::pair<std::set<std::size_t>,
                                           std::pair<std::string, std::string>>>
                 &conformational_inter_unit_cell)
      : base_type{conformational_changes, conformational_changes_names,
                  agonist_changes,        conductance_changes,
                  conductance_names,      conformational_inter_unit_cell} {}

  template <class Parameter>
  static Derivative<double> rate(const transitions &tr,
                                 const Derivative<Parameter> &p) {
    Derivative<double> out(p.x());
    if (tr.on) {
      for (auto &e : tr.coupling) {
        Derivative<double> b(e.second, p.x());
        for (auto &e2 : e.first)
          b = b * pow(p.at(e2.first), p.at(e2.second));
        out += b;
      }
    } else {
      for (auto &e : tr.coupling) {
        Derivative<double> b(e.second, p.x());
        for (auto &e2 : e.first)
          b = b * pow(p.at(e2.first), p.at(e2.second) + Constant(-1.0));
        out += b;
      }
    }
    return out * p.at(tr.conformation);
  }

  template <class Parameters> auto Qs(const Derivative<Parameters> &p) const {

    M_Matrix<Derivative<double>> Q0(conformer_.size(), conformer_.size(),
                                    Matrix_TYPE::FULL,
                                    Derivative<double>(p.x()));
    M_Matrix<Derivative<double>> Qa(conformer_.size(), conformer_.size(),
                                    Matrix_TYPE::FULL,
                                    Derivative<double>(p.x()));

    for (std::size_t i = 0; i < Q0.nrows(); ++i) {
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
    //    auto r=Q0*ones<Constant<double>>(Q0.ncols(),1);
    //    auto q=Qa*ones<Constant<double>>(Q0.ncols(),1);

    return std::pair(Derivative<M_Matrix<double>>(std::move(Q0)),
                     Derivative<M_Matrix<double>>(std::move(Qa)));
  }

  template <class Parameters> auto g(const Derivative<Parameters> &p) const {
    M_Matrix<Derivative<double>> out(conformer_.size(), 1, Matrix_TYPE::FULL,
                                     Derivative<double>(p.x()));
    for (std::size_t i = 0; i < conformer_.size(); ++i)
      out(i, 0) = p.at(conductances_.at(i));
    return Derivative<M_Matrix<double>>(out);
  }
};

template <> class Derivative<Markov_Transition_rate> {

  Derivative<M_Matrix<double>> Qrun_; // transition rate matrix at time zero
  Derivative<M_Matrix<double>> g_;

  Derivative<M_Matrix<double>> V_;     // eigenvector of Qrun
  Derivative<M_Matrix<double>> W_;     // eigenvector of Qrun
  Derivative<M_Matrix<double>> landa_; // eigenvalues

  Derivative<M_Matrix<double>> Wg_;
  Derivative<M_Matrix<double>> WgV_;

  static void clean_landa(Derivative<M_Matrix<double>> &la) {
    double maxla = std::min(max(la.f()), 0.0);
    for (std::size_t i = 1; i < la.f().size(); ++i)
      if (la.f()[i] >= maxla) {
        la.f()[i] = 0;
        la.dfdx().transform([i](auto &x) { x[i] = 0; });
      }
  }

  void init(Derivative_t<Matrix_Decompositions::eigensystem_type> &&eig) {
    std::tie(V_, landa_, W_) = std::move(eig);
    clean_landa(landa_);
    Wg_ = W_ * g_;
    assert((are_Equal<true, Derivative<M_Matrix<double>>>().test_prod(
        Wg_,
        Incremental_ratio(
            1e-6, [](auto const &W, auto const &g) { return W * g; }, W_, g_),
        std::cerr)));

    WgV_ = W_ * diag(g_) * V_;

    // auto WV=W_*V_;
    // auto dWV=Incremental_ratio([](1e-6,auto const &W,
    //                                auto const &V) { return W * V; },
    //                            W_, V_);

    // are_Equal_v(WV,dWV,std::cerr);

    // auto VWg = V_*W_ * g_;
    // auto dVWg = Incremental_ratio(1e-6,
    //   [](auto const &V, auto const &W, auto const & g) { return  V*W*g; },
    //   V_,W_,g_);

    // are_Equal_v(VWg, dVWg, std::cerr);

    assert((are_Equal<true, Derivative<M_Matrix<double>>>().test_prod(
        Incremental_ratio(1e-5,
                          [](auto const &W, auto const &g, auto const &V) {
                            return W * diag(g) * V;
                          },
                          W_, g_, V_),
        Incremental_ratio(1e-6,
                          [](auto const &W, auto const &g, auto const &V) {
                            return W * diag(g) * V;
                          },
                          W_, g_, V_),
        std::cerr)));

    assert((are_Equal<true, Derivative<M_Matrix<double>>>().test_prod(
        WgV_,
        Incremental_ratio(1e-6,
                          [](auto const &W, auto const &g, auto const &V) {
                            return W * diag(g) * V;
                          },
                          W_, g_, V_),
        std::cerr)));

    // assert(Frobenius_test::test(Qrun(),V()*landa()*W(),1,
    // std::sqrt(std::numeric_limits<double>::epsilon())));
  }

public:
  static Op_void
  test(const Derivative_t<Matrix_Decompositions::eigensystem_type> &e,
       double tolerance) {
    return Markov_Transition_rate::landa::test(std::get<1>(e).f(), tolerance);
  }

  static myOptional_t<Derivative<Markov_Transition_rate>>
  evaluate(Derivative<M_Matrix<double>> &&_Qrun,
           const Derivative<M_Matrix<double>> &_g, double tolerance) {
    typedef myOptional_t<Derivative<Markov_Transition_rate>> Op;
    auto eig = Matrix_Decompositions::EigenSystem_full_real_eigenvalue(_Qrun);
    if (eig) {
      // std::cerr<<"\n --------1e-5  vs 1e-6-------\n";

      if constexpr (false){
  are_Equal<true, Derivative_t<
                          typename Matrix_Decompositions::eigensystem_type>>()
          .test(Incremental_ratio_tuple_3(
                    1e-5,
                    [](auto &x) {
                      return Matrix_Decompositions::
                          EigenSystem_full_real_eigenvalue(x)
                              .value();
                    },
                    _Qrun),
                Incremental_ratio_tuple_3(
                    1e-6,
                    [](auto &x) {
                      return Matrix_Decompositions::
                          EigenSystem_full_real_eigenvalue(x)
                              .value();
                    },
                    _Qrun),
                std::cerr);
      // std::cerr<<"\n --------1e-6-------\n";

      //bool test_e6 =
          are_Equal<
              true,
              Derivative_t<typename Matrix_Decompositions::eigensystem_type>>()
              .test(eig.value(),
                    Incremental_ratio_tuple_3(
                        1e-6,
                        [](auto &x) {
                          return Matrix_Decompositions::
                              EigenSystem_full_real_eigenvalue(x)
                                  .value();
                        },
                        _Qrun),
                    std::cerr);
      }
      if (auto eigtest = test(eig.value(), tolerance); eigtest.has_value())
        return Op(Derivative(std::move(_Qrun), _g, std::move(eig).value()));
      else
        return Op(false, "invalid eigenvalues " + eigtest.error());
    } else
      return Op(false, "eigenvalue decomposition fails:" + eig.error());
  }

  Derivative() = default;
  /// virtual copy constructors

  Derivative(Derivative<M_Matrix<double>> &&_Qrun,
             const Derivative<M_Matrix<double>> &_g,
             Derivative_t<Matrix_Decompositions::eigensystem_type> &&eig)
      : Qrun_{std::move(_Qrun)}, g_{_g} {
    init(std::move(eig));
  }

  auto &Qrun() const { return Qrun_; }   // transition rate matrix at time zero
  auto &V() const { return V_; }         // eigenvector of Qrun
  auto &W() const { return W_; }         // eigenvector of Qrun
  auto &landa() const { return landa_; } // eigenvalues

  auto &g() const { return g_; }
  auto &Wg() const { return Wg_; }
  auto &WgV() const { return WgV_; }
  auto &x() const { return g().x(); }
  auto calc_Peq() const {
    M_Matrix<double> p0(1, Qrun().f().nrows(), 1.0 / Qrun().f().nrows());
    M_Matrix<double> laexp(landa().f().size(), landa().f().size(),
                           Matrix_TYPE::DIAGONAL);
    for (std::size_t i = 0; i < landa().f().size(); ++i) {
      if (landa().f()(i, i) == 0.0)
        laexp(i, i) = 1.0;
      else
        laexp(i, i) = 0.0;
    }
    //  std::cerr<<"\n
    //  ----calc_Peq-----\np0="<<p0<<"laexp="<<laexp<<"V()"<<V()<<"W"<<W();

    return Constant(p0) * V() * Constant(laexp) * W();
  }
  Derivative(
      const Derivative<M_Matrix<double>>& Qrun, // transition rate matrix at time zero
      Derivative<M_Matrix<double>> const& g,

      Derivative<M_Matrix<double>> const&V,     // eigenvector of Qrun
      Derivative<M_Matrix<double>> const&W,     // eigenvector of Qrun
      Derivative<M_Matrix<double>> const&landa, // eigenvalues

      Derivative<M_Matrix<double>> const&Wg, Derivative<M_Matrix<double>> const&WgV)
      : Qrun_{Qrun}, g_{g}, V_{V}, W_{W}, landa_{landa}, Wg_{Wg}, WgV_{WgV} {}





  typedef Markov_Transition_rate self_type;

constexpr static auto const className =
    my_static_string("Markov_Transition_rate");
// std::string myClass()const  { return className.str();}

static auto get_constructor_fields() {
  return std::make_tuple(
      grammar::field(C<self_type>{}, "Qrun", &self_type::Qrun),
      grammar::field(C<self_type>{}, "g", &self_type::g),
      grammar::field(C<self_type>{}, "V", &self_type::V),
      grammar::field(C<self_type>{}, "W", &self_type::W),
      grammar::field(C<self_type>{}, "landa", &self_type::landa),
      grammar::field(C<self_type>{}, "Wg", &self_type::Wg),
      grammar::field(C<self_type>{}, "WgV", &self_type::WgV));
}
}
;

template <> class Derivative<Markov_Transition_step> {
protected:
  Derivative<M_Matrix<double>> P_; // transition matrix

  Derivative<M_Matrix<double>> g_; // conductance matrix
  std::size_t nsamples_;
  double dt_;
  double min_p_;

public:
  double min_P() const { return min_p_; }
  auto const &myP() const { return P_; }  // transition matrix
  auto const &P() const & { return P_; }  // transition matrix
  auto &&P() && { return std::move(P_); } // transition matrix

  double dt() const { return dt_; }
  double fs() const { return 1.0 / dt(); }

  std::size_t nsamples() const { return nsamples_; }

  /// mean conductance for each starting state i
  auto const &g() const { return g_; } // conductance matrix

  Derivative(const Derivative<Markov_Transition_rate> &Qx, std::size_t nsamples,
             double fs, double min_p)
      : P_{}, /*ladt_{},exp_ladt_{},*/ g_{Qx.g()}, nsamples_{nsamples},
        dt_{1.0 / fs * nsamples}, min_p_{min_p} {
    init(Qx);
  }

  Derivative(std::size_t n_sub_samples,
             const Derivative<Markov_Transition_rate> &Qx, std::size_t nsamples,
             double fs, double min_p)
      : P_{}, /*ladt_{},exp_ladt_{},*/ g_{Qx.g()}, nsamples_{nsamples},
        dt_{1.0 / (fs * n_sub_samples) * nsamples}, min_p_{min_p} {
    init(Qx);
  }

  Derivative(Derivative<M_Matrix<double>> &&P, Derivative<M_Matrix<double>> &&g,
             std::size_t nsamples, double fs, double min_p)
      : P_{std::move(P)},
        //          P_{Probability_transition::normalize(std::move(P),min_p)},
        g_{std::move(g)}, nsamples_{nsamples}, dt_{(1.0 * nsamples) / fs},
        min_p_{min_p} {}

  Derivative<Markov_Transition_step> getMinimum() const { return *this; }

  Derivative<Markov_Transition_step> &
  operator*=(const Derivative<Markov_Transition_step> &other) {
    P_ = P_ * other.P();
    P_.f() = Probability_transition::normalize(std::move(P_.f()), min_P());
    g_ = other.g();
    nsamples_ += other.nsamples();
    return *this;
  }

  Derivative() = default;

  typedef Derivative self_type;
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
  void init(const Derivative<Markov_Transition_rate> &Qx) {
    auto ladt = Qx.landa() * dt();
    auto exp_ladt = exp(ladt);
    P_ = Qx.V() * exp_ladt * Qx.W();
    P_.f() = Probability_transition::normalize(std::move(P_.f()), min_P());
    g_ = Qx.g();
  }
};

template <> class Derivative<Markov_Transition_step_double_minimum> {
public:
    std::size_t n;
  Derivative<M_Matrix<double>> PPn;
  Derivative<M_Matrix<double>> PG_n;
  Derivative<M_Matrix<double>> PGG_n;
  double min_p;
  double tolerance;
  Derivative(const Derivative<Markov_Transition_step_double> &x_);
  Derivative operator*=(const Derivative<Markov_Transition_step_double> &x_);
};

template <>
class Derivative<Markov_Transition_step_double>
    : public Derivative<Markov_Transition_step> {
  Derivative<M_Matrix<double>> gmean_i_;   // conductance matrix
  Derivative<M_Matrix<double>> gtotal_ij_; // conductance matrix
  Derivative<M_Matrix<double>> gmean_ij_;  // conductance matrix

  Derivative<M_Matrix<double>> gtotal_sqr_ij_; // conductance matrix

  Derivative<M_Matrix<double>> gsqr_i_; // conductance matrix

  Derivative<M_Matrix<double>> gvar_i_; // variance of the conductance matrix
  Derivative<M_Matrix<double>>
      gtotal_var_ij_; // variance of the conductance matrix

  Derivative<M_Matrix<double>> gvar_ij_; // variance of the conductance matrix

public:
  auto &x() const { return g().x(); }
  typedef Derivative<Markov_Transition_step> base_type;
  Derivative(Derivative<M_Matrix<double>> &&P, Derivative<M_Matrix<double>> &&g,
             std::size_t nsamples, double fs, double min_p,
             Derivative<M_Matrix<double>> &&gmean_i,
             Derivative<M_Matrix<double>> &&gtotal_ij,
             Derivative<M_Matrix<double>> &&gmean_ij,

             Derivative<M_Matrix<double>> &&gtotal_sqr_ij,

             Derivative<M_Matrix<double>> &&gsqr_i,

             Derivative<M_Matrix<double>> &&gvar_i,
             Derivative<M_Matrix<double>> &&gtotal_var_ij,
             Derivative<M_Matrix<double>> &&gvar_ij)
      : base_type(std::move(P), std::move(g), nsamples, fs, min_p),
        gmean_i_{gmean_i}, gtotal_ij_{gtotal_ij}, gmean_ij_{gmean_ij},
        gtotal_sqr_ij_{gtotal_sqr_ij}, gsqr_i_{gsqr_i}, gvar_i_{gvar_i},
        gtotal_var_ij_{gtotal_var_ij}, gvar_ij_{gvar_ij} {}

  Derivative<Markov_Transition_step_double_minimum> getMinimum() const {
    return Derivative<Markov_Transition_step_double_minimum>(*this);
  }

  /// mean conductance for each starting state i
  auto const &gmean_i() const { return gmean_i_; } // conductance matrix
  /// total conductance for each starting state i and ending state j
  auto const &gtotal_ij() const { return gtotal_ij_; } // conductance matrix
  /// mean conductance for each starting state i and ending state j
  auto const &gmean_ij() const { return gmean_ij_; } // conductance matrix

  /// squared mean conductance for each starting state i and ending state j
  auto const &gtotal_sqr_ij() const {
    return gtotal_sqr_ij_;
  } // conductance matrix

  /// squared mean conductance for each starting state i
  auto const &gsqr_i() const { return gsqr_i_; } // conductance matrix

  /// variance of the mean conductance for each starting state i
  auto const &gvar_i() const {
    return gvar_i_;
  } // variance of the conductance matrix
  /// variance of the mean conductance for each starting state i contributed by
  /// the ones ending at state j
  auto const &gtotal_var_ij() const {
    return gtotal_var_ij_;
  } // variance of the conductance matrix

  /// variance of the mean conductance for each starting state i and ending
  /// state j
  auto const &gvar_ij() const {
    return gvar_ij_;
  } // variance of the conductance matrix

  Derivative(const Derivative<Markov_Transition_rate> &Qx, std::size_t nsamples,
             double fs, double min_p)
      : base_type(Qx, nsamples, fs, min_p) {
    init(Qx);
  }

  Derivative(std::size_t recursion_number,
             const Derivative<Markov_Transition_rate> &Qx, std::size_t nsamples,
             double fs, double min_p)
      : base_type(Qx, nsamples, fs, min_p) {
    init_by_step(Qx, recursion_number, fs);
  }

  Derivative(Derivative<Markov_Transition_step> &&step) : base_type(step) {
    init_by_step_init();
  }

  Derivative(Derivative<Markov_Transition_step_double_minimum> &&x,
             Derivative<M_Matrix<double>> g, double fs, double min_p, double)
      : base_type(std::move(x.PPn), std::move(g), x.n, fs, min_p) {
    init(std::move(x), fs);
  }

  Derivative() = default;

  typedef Derivative<Markov_Transition_step_double> self_type;
  constexpr static auto className =
      my_static_string("Markov_Transition_step_double_derivative");
  static auto get_constructor_fields() {
    // M_Matrix<double> const & (self_type::*myP) ()const& =&self_type::P;

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

private:
  static double E1(double x) { return Markov_Transition_step_double::E1(x); }
  static double dE1(double x) {
    if (std::abs(x) < std::numeric_limits<double>::epsilon() * 100)
      return 0.5;
    else
      return (std::exp(x) * x - std::exp(x) + 1) / sqr(x);
  }

  static double dE2dx(double x, double y) {

    const double eps = std::numeric_limits<double>::epsilon();
    if (x * x < eps) {
      return 0;
    } else if (y * y < eps)
      return (dE1(x) * x - E1(x) + 1.0) / sqr(x);
    else if ((y - x) * (y - x) < eps)
      return ((std::exp(x) * x - dE1(x) * x) - (std::exp(x) - E1(x))) / sqr(x);
    else
      return (dE1(x) * (x - y) - (E1(x) - E1(y))) / sqr(x - y);
  }
  static double Ee(double x, double y, double exp_x, double exp_y,
                   double eps = std::numeric_limits<double>::epsilon()) {
    return Markov_Transition_step_double::Ee(x, y, exp_x, exp_y, eps);
  }
  static double dEedx(double x, double y, double exp_x, double exp_y,
                      double eps = std::numeric_limits<double>::epsilon()) {
    if (sqr(x - y) < eps)
      return exp_x / 2.0;
    else
      return (exp_x * (x - y) - (exp_x - exp_y)) / sqr(x - y);
  };
  static Derivative<double>
  Ee(Derivative<double> x, Derivative<double> y, double exp_x, double exp_y,
     double eps = std::numeric_limits<double>::epsilon()) {
    return Derivative<double>(
        Ee(x.f(), y.f(), exp_x, exp_y, eps), x.x(),
        dEedx(x.f(), y.f(), exp_x, exp_y, eps) * x.dfdx() +
            dEedx(y.f(), x.f(), exp_y, exp_x, eps) * y.dfdx());
  }

  static double EX_111(double x, double y, double z, double exp_x) {
    return Markov_Transition_step_double::EX_111(x, y, z, exp_x);
  }
  //  static double EX_111(double x, double y, double z, double exp_x) {
  //    return exp_x / ((x - y) * (x - z));
  //  }

  static double dEX_111_dx(double x, double y, double z, double exp_x) {
    // (u'v-uv')/v2
    // return   (exp_x*((x - y) * (x - z))-(exp_x*( (x - z)+(x - y) )))/sqr((x -
    // y) * (x - z));

    return exp_x * ((x - y) * (x - z) - (2 * x - z - y)) /
           sqr((x - y) * (x - z));
  }
  static double dEX_111_dy(double x, double y, double z, double exp_x) {
    // -v'/v2

    return exp_x / (x - z) / sqr(x - y);
  }

  static double E111(double x, double y, double z, double exp_x, double exp_y,
                     double exp_z) {
    return EX_111(x, y, z, exp_x) + EX_111(y, x, z, exp_y) +
           EX_111(z, y, x, exp_z);
  }
  static double dE111_dx(double x, double y, double z, double exp_x,
                         double exp_y, double exp_z) {
    return dEX_111_dx(x, y, z, exp_x) + dEX_111_dy(y, x, z, exp_y) +
           dEX_111_dy(z, x, y, exp_z);
  }
  static double dEX_111_dyy(double x, double y, double exp_x) {

    //    return exp_x / sqr(x - y);
    return 2 * exp_x / std::pow(x - y, 3);
  }

  static double E12(double x, double y, double exp_x, double exp_y) {
    return EX_111(x, y, y, exp_x) + exp_y / (y - x) * (1.0 - 1.0 / (y - x));
  }

  static double dE12_dx(double x, double y, double exp_x, double exp_y) {

    // exp_y / (y - x) * (1.0 - 1.0 / (y - x));

    //   +exp_y / sqr(y - x) * (1.0 - 1.0 / (y - x))
    //  -exp_y / (y - x) * ( 1.0 / sqr(y - x))
    return dEX_111_dx(x, y, y, exp_x) +
           exp_y / sqr(y - x) * (1.0 - 1.0 / (y - x)) -
           exp_y / (y - x) * (1.0 / sqr(y - x));
  }
  static double dE12_dy(double x, double y, double exp_x, double exp_y) {
    return dEX_111_dyy(x, y, exp_x) + exp_y / (y - x) * (1.0 - 1.0 / (y - x)) -
           exp_y / sqr(y - x) * (1.0 - 1.0 / (y - x)) +
           exp_y / (y - x) * (1.0 / sqr(y - x));
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

  static double dE3dx(double x, double y, double z, double exp_x, double exp_y,
                      double exp_z,
                      double eps = std::numeric_limits<double>::epsilon()) {
    if (sqr(x - y) < eps) // x==y
    {
      if (sqr(y - z) < eps)       // y==z
        return exp_x / 2.0 / 3.0; // x==y==z
      else
        return dE12_dy(z, x, exp_z, exp_x) / 2.0; // x==y!=z
    } else if (sqr(y - z) < eps)                  // x!=y==z
    {
      return dE12_dx(x, y, exp_x, exp_y);
    } else if (sqr(x - z) < eps) // y!=z==x!=y
    {
      return dE12_dy(y, x, exp_y, exp_x) / 2.0;
    } else
      return dE111_dx(x, y, z, exp_x, exp_y, exp_z); // x!=y!=z!=x
  }

  static Derivative<double>
  E3(Derivative<double> x, Derivative<double> y, Derivative<double> z,
     double exp_x, double exp_y, double exp_z,
     double eps = std::numeric_limits<double>::epsilon()) {
    auto df = dE3dx(x.f(), y.f(), z.f(), exp_x, exp_y, exp_z, eps) * x.dfdx() +
              dE3dx(y.f(), z.f(), x.f(), exp_y, exp_z, exp_x, eps) * y.dfdx() +
              dE3dx(z.f(), x.f(), y.f(), exp_z, exp_x, exp_y, eps) * z.dfdx();

    return Derivative<double>(E3(x.f(), y.f(), z.f(), exp_x, exp_y, exp_z, eps),
                              x.x(), std::move(df));
  }

  void init(Derivative<Markov_Transition_step_double_minimum> &&x, double) {
    M_Matrix<double> u = ones<double>(P().f().ncols(), 1);

    gtotal_sqr_ij_ = x.PGG_n * (1.0 / (x.n * x.n * 0.5));
    gtotal_ij_ = x.PG_n * (1.0 / x.n);
    gmean_ij_ = elemDivSafe(gtotal_ij_, P_, min_P());

    gtotal_var_ij_ = gtotal_sqr_ij_ - elemMult(gtotal_ij_, gmean_ij_);

    gmean_i_ = gtotal_ij_ * u;

    gsqr_i_ = gtotal_sqr_ij_ * u;
    gvar_ij_ = elemDivSafe(gtotal_var_ij_, P(), min_P());
    gmean_i_ = gtotal_ij_ * u;
    gsqr_i_ = gtotal_sqr_ij_ * u;
    gvar_i_ = gtotal_var_ij_ * u;
  }

  void init(const Derivative<Markov_Transition_rate> &Qx) {
    // const double eps=std::numeric_limits<double>::epsilon();
    double dt = base_type::dt();

    std::size_t N = P().f().ncols();

    auto ladt = Qx.landa().to_Matrix() * dt;

    auto exp_ladt = exp(Qx.landa().f() * dt);

    //auto minP = min_P();
    M_Matrix<Derivative<double>> E2m(N, N, Matrix_TYPE::SYMMETRIC);
    M_Matrix<Derivative<double>> E2mb(N, N, Matrix_TYPE::SYMMETRIC);
    for (std::size_t i = 0; i < N; ++i)
      for (std::size_t j = 0; j < i + 1; ++j) {
        E2m(i, j) = Ee(ladt[i], ladt[j], exp_ladt[i], exp_ladt[j], min_P());
        assert(are_Equal_v(E2m(i, j),
                           Incremental_ratio(1e-4,
                                             [&minP](auto ladt1, auto ladt2) {
                                               return Ee(ladt1, ladt2,
                                                         std::exp(ladt1),
                                                         std::exp(ladt2), minP);
                                             },
                                             ladt[i], ladt[j]),
                           std::cerr));

        assert(Derivative_correctness_mean_value_test(1e-2,1e4,
                                             [&minP](auto ladt1, auto ladt2) {
                                               return Ee(ladt1, ladt2,
                  std::exp(ladt1.f()),
                  std::exp(ladt2.f()), minP);
                                             },
                                                      std::forward_as_tuple(ladt[i], ladt[j]), E2m(i, j),
                                                      std::cerr, "i= ", i, " j=", j, "ladt1= ", ladt[i], "ladt2= ", ladt[j]));

      }

    // build E2
    M_Matrix<Derivative<double>> WgV_E2(N, N);

    auto WgV = Qx.WgV().to_Matrix();
    for (std::size_t i = 0; i < N; ++i)
      for (std::size_t j = 0; j < N; ++j)
        WgV_E2(i, j) = WgV(i, j) * E2m(i, j);

    gtotal_ij_ = Qx.V() * Derivative<M_Matrix<double>>(WgV_E2) * Qx.W();

    M_Matrix<Derivative<double>> WgV_E3(N, N,
                                        Derivative<double>(Qx.Qrun().x()));
    //bool succed = true;
    for (std::size_t n1 = 0; n1 < N; n1++)
      for (std::size_t n3 = 0; n3 < N; n3++)
        for (std::size_t n2 = 0; n2 < N; n2++) {
          auto E3v=E3(ladt[n1], ladt[n2], ladt[n3], exp_ladt[n1], exp_ladt[n2],
                        exp_ladt[n3], min_P());
          WgV_E3(n1, n3) +=
              WgV(n1, n2) * WgV(n2, n3) * E3v;
//              E3(ladt[n1], ladt[n2], ladt[n3], exp_ladt[n1], exp_ladt[n2],
//                 exp_ladt[n3], min_P()); // optimizable
          /*if (!are_Equal_v(E3(ladt[n1], ladt[n2], ladt[n3], exp_ladt[n1],
                              exp_ladt[n2], exp_ladt[n3], min_P()),
                           Incremental_ratio(
                               1e-2,
                               [&minP](auto ladt1, auto ladt2, auto ladt3) {
                                 return E3(ladt1, ladt2, ladt3, std::exp(ladt1),
                                           std::exp(ladt2), std::exp(ladt3),
                                           minP);
                               },
                               ladt[n1], ladt[n2], ladt[n3]),
                           std::cerr, " n1=", n1, " n2=", n2, " n3=", n3))
            succed = false; */
          assert(Derivative_correctness_mean_value_test(
              1e-2,1e6,
                               [&minP](auto ladt1, auto ladt2, auto ladt3) {
          return E3(ladt1, ladt2, ladt3, std::exp(ladt1.f()), std::exp(ladt2.f()),
                    std::exp(ladt3.f()), minP);
                               },
                               std::forward_as_tuple(ladt[n1], ladt[n2], ladt[n3]),E3v, std::cerr,
              "n1= ", n1, " n2=", n2, " n3=", n3, "\nladt1= ", ladt[n1], "ladt2= ", ladt[n2], "ladt3= ", ladt[n3]));
        //    succed = false;
        }
    assert(succed);
    gtotal_sqr_ij_ =
        Qx.V() * Derivative<M_Matrix<double>>(WgV_E3) * Qx.W() * 2.0;
    for (std::size_t i = 0; i < N; ++i)
      for (std::size_t j = 0; j < N; ++j)
        if (P().f()(i, j) == 0) {
          gtotal_ij_.f()(i, j) = 0;
          gtotal_sqr_ij_.f()(i, j) = 0;
        }
    if constexpr (true)
    {

    M_Matrix<double> U = Matrix_Generators::ones<double>(1, g().f().size());
    M_Matrix<double> UU = M_Matrix<double>(g().f().size(), g().f().size(), 1);
    auto gmean_ij_p = TransposeSum(g() * Constant(U)) * (0.5);
    auto gvar_ij_p = (g() * Constant(U) - Transpose(g() * Constant(U))) * 0.5;
    gvar_ij_p =
        gvar_ij_p.apply([](Derivative<double> const &x) { return abs(x); });

    // std::cerr<<"\ngmean_ij_p=\n"<<gmean_ij_p<<"\ngvar_ij_p=\n"<<gvar_ij_p<<"\n";
    // std::cerr<<"\n UU=ņ"<<UU<<"\n";
    auto gmean_ij_tot = gtotal_ij_ + gmean_ij_p * min_P();
    auto P_p = P() + Constant(UU * min_P());
    gmean_ij_ = elemDiv(gmean_ij_tot, P_p);
    gtotal_var_ij_ = gtotal_sqr_ij_ - elemMult(gtotal_ij_, gmean_ij_);
    auto gvar_ij_tot = gtotal_var_ij_ + gvar_ij_p * min_P();
    gvar_ij_ = elemDiv(gvar_ij_tot, P_p);
    M_Matrix<double> u(N, 1, 1.0);
    gmean_i_ = gtotal_ij_ * Constant(u);
    gsqr_i_ = gtotal_sqr_ij_ * Constant(u);
    gvar_i_ = gtotal_var_ij_ * Constant(u);
  } else
    {

      gmean_ij_ = elemDiv(gtotal_ij_, P_);
      gtotal_var_ij_ = gtotal_sqr_ij_ - elemMult(gtotal_ij_, gmean_ij_);
      gvar_ij_ = elemDiv(gtotal_var_ij_, P_);
      M_Matrix<double> u(N, 1, 1.0);
      gmean_i_ = gtotal_ij_ * Constant(u);
      gsqr_i_ = gtotal_sqr_ij_ * Constant(u);
      gvar_i_ = gtotal_var_ij_ * Constant(u);
    }
  }

  void init_by_step_init() {

    auto G = diag(g());
    std::size_t N = P().f().nrows();
    gtotal_ij_ = (G * P() + P() * G) * 0.5;
    gtotal_var_ij_ = (G * P() - P() * G) * 0.5;
    gtotal_var_ij_.f() =
        gtotal_var_ij_.f().apply([](double x) { return std::abs(x); });
    gmean_ij_ = elemDivSafe(gtotal_ij_, P_);
    gtotal_sqr_ij_ = gtotal_var_ij_ + elemMult(gtotal_ij_, gmean_ij_);
    gvar_ij_ = elemDivSafe(gtotal_var_ij_, P_);
    M_Matrix<double> u(N, 1, 1.0);
    gmean_i_ = gtotal_ij_ * u;
    gsqr_i_ = gtotal_sqr_ij_ * u;
    gvar_i_ = gtotal_var_ij_ * u;
  }

  void init_by_step(const Derivative<Markov_Transition_rate> &Qx,
                    std::size_t recursion_number, double fs) {
    std::size_t times = std::pow(2, recursion_number);
    std::size_t n = nsamples();
    Derivative<Markov_Transition_step_double> s(Qx, n, fs * times, min_P());
    Derivative<Markov_Transition_step_double> step_ini(s);
    for (std::size_t i = 0; i < recursion_number; ++i) {
      Derivative<Markov_Transition_step_double_minimum> step(step_ini);
      step *= step_ini;
      step_ini = Derivative<Markov_Transition_step_double>(
          std::move(step), step_ini.g(), fs, min_P(), 0);
    }
    *this = step_ini;
  }
};


std::ostream& operator<<(std::ostream& os, const Derivative<Markov_Transition_step_double>& d ){return io::output_operator_on_Object(os,d);}

template <> class Derivative<SingleLigandModel> : public SingleLigandModel {
public:
  typedef SingleLigandModel base_type;
  auto &x() const { return Q0_.x(); }

  Derivative(){};
  SingleLigandModel const &f() const

  {
    return *this;
  }
  auto Q(double x) const { return Q0_ + Qa_ * x; }

  auto &g(double) const { return g_; }
  auto Qx(double x, double tolerance) const {
    auto out =
        Derivative_t<Markov_Transition_rate>::evaluate(Q(x), g(x), tolerance);
    if (out)
    assert((Derivative_correctness_mean_value_test(
        1e-2,1e4,
        [&tolerance](auto Qx, auto gx) {
          return std::move(Derivative_t<Markov_Transition_rate>::evaluate(
                               std::move(Qx), gx, tolerance))
              .value();
        },
        std::forward_as_tuple(Q(x), g(x)), out.value(), std::cerr, "x= ", x,
        "  tolerance=", tolerance)));
    return out;
  }

  auto P(const Derivative_t<Markov_Transition_rate> &aQx, std::size_t n,
         double fs) const {
    return Derivative<Markov_Transition_step>(aQx, n, fs, min_P_);
  }

  auto Pdt(const Markov_Transition_rate &x, std::size_t n, double fs) const {
    return Markov_Transition_step_double(x, n, fs, min_P_);
  }

  std::size_t nstates() const { return Q0_.f().nrows(); }

  auto noise_variance(std::size_t nsamples, double fs) const {
    return noise_variance_ * (fs / nsamples);
  }

  auto &AverageNumberOfChannels() const { return N_channels_; }

  std::size_t N_channels() const { return N_channels_.f(); }

  Derivative(std::pair<Derivative<M_Matrix<double>>,
                       Derivative<M_Matrix<double>>> &&Qs,
             Derivative<M_Matrix<double>> &&g, double Vm, Derivative<double> N,
             Derivative<double> noise, double min_P)
      : SingleLigandModel({Qs.first.f(), Qs.second.f()}, g.f(), Vm, N.f(),
                          noise.f(), min_P),
        Q0_{std::move(Qs.first)}, Qa_{std::move(Qs.second)}, g_{g * Vm},
        N_channels_{N}, noise_variance_{noise}, min_P_{min_P} {}


  Derivative(const Derivative<M_Matrix<double>>& Q0,
                       const Derivative<M_Matrix<double>>&Qa,
             const Derivative<M_Matrix<double>> &g, const Derivative<double>& N,
             const Derivative<double>& noise, double min_P)
      : SingleLigandModel(Q0.f(), Qa.f(), g.f(), N.f(),
                          noise.f(), min_P),
        Q0_{Q0}, Qa_{Qa}, g_{g },
        N_channels_{N}, noise_variance_{noise}, min_P_{min_P} {}



  double min_P() const { return min_P_; }

  auto &Q0() const { return Q0_; }
  auto &Qa() const { return Qa_; }
  auto &g0() const { return g_; }
  auto &noise_std() const { return noise_variance_; }

  typedef Derivative self_type;

  constexpr static auto className = my_static_string("SingleLigandModel_derivative");
  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "Q0", &self_type::Q0),
        grammar::field(C<self_type>{}, "Qa", &self_type::Qa),
        grammar::field(C<self_type>{}, "g", &self_type::myg),
        grammar::field(C<self_type>{}, "AverageNumberOfChannels", &self_type::N_channels),
        grammar::field(C<self_type>{}, "plogL", &self_type::noise_variance_1Hz),
        grammar::field(C<self_type>{}, "eplogL", &self_type::min_P));
  }

 private:
  Derivative<M_Matrix<double>> Q0_;
  Derivative<M_Matrix<double>> Qa_;
  Derivative<M_Matrix<double>> g_;
  Derivative<double> N_channels_;
  Derivative<double> noise_variance_;
  double min_P_;
};



static_assert (is_field_Object<Derivative<SingleLigandModel>>::value,"" );

inline SingleLigandModel Taylor_first(const Derivative<SingleLigandModel> &dy,
                                      std::size_t i, double eps) {
  auto Q0 = Taylor_first(dy.Q0(), i, eps);
  auto Qa = Taylor_first(dy.Qa(), i, eps);
  auto g = Taylor_first(dy.g0(), i, eps);
  auto N = Taylor_first(dy.AverageNumberOfChannels(), i, eps);
  auto noise = Taylor_first(dy.noise_std(), i, eps);
  return SingleLigandModel(std::pair(std::move(Q0), std::move(Qa)),
                           std::move(g), 1.0, N, noise, dy.min_P());
}

inline Derivative<SingleLigandModel>
Directional_Derivative(const Derivative<SingleLigandModel> &dy,
                       const M_Matrix<double> &new_x, std::size_t i,
                       double eps) {
  auto Q0 = Directional_Derivative(dy.Q0(), new_x, i, eps);
  auto Qa = Directional_Derivative(dy.Qa(), new_x, i, eps);
  auto g = Directional_Derivative(dy.g0(), new_x, i, eps);
  auto N = Directional_Derivative(dy.AverageNumberOfChannels(), new_x, i, eps);
  auto noise = Directional_Derivative(dy.noise_std(), new_x, i, eps);
  return Derivative<SingleLigandModel>(std::pair(std::move(Q0), std::move(Qa)),
                                       std::move(g), 1.0, N, noise, dy.min_P());
}

template <class Markov_Transition_step, class Markov_Transition_rate, class M,
          class Experiment, class X>
Markov_Model_calculations<Markov_Transition_step, Markov_Transition_rate, M,
                          Experiment, X>
Primitive(const Markov_Model_calculations<Derivative<Markov_Transition_step>,
                                          Derivative<Markov_Transition_rate>,
                                          Derivative<M>, Experiment, X> &dm) {
  return Markov_Model_calculations<Markov_Transition_step,
                                   Markov_Transition_rate, M, Experiment, X>(
      dm.Model().f(), dm.get_Experiment(), dm.n_sub_samples(), dm.tolerance());
}

template <class Markov_Transition_step, class Markov_Transition_rate, class M,
          class Experiment, class X>
Markov_Model_calculations<Markov_Transition_step, Markov_Transition_rate, M,
                          Experiment, X>
Taylor_first(const Markov_Model_calculations<Derivative<Markov_Transition_step>,
                                             Derivative<Markov_Transition_rate>,
                                             Derivative<M>, Experiment, X> &dm,
             std::size_t i, double eps) {
  return Markov_Model_calculations<Markov_Transition_step,
                                   Markov_Transition_rate, M, Experiment, X>(
      Taylor_first(dm.Model(), i, eps), dm.get_Experiment(), dm.n_sub_samples(),
      dm.tolerance());
}

template <class Markov_Transition_step, class Markov_Transition_rate, class M,
          class Experiment, class X>
Markov_Model_calculations<Derivative<Markov_Transition_step>,
                          Derivative<Markov_Transition_rate>, Derivative<M>,
                          Experiment, X>
Directional_Derivative(
    const Markov_Model_calculations<Derivative<Markov_Transition_step>,
                                    Derivative<Markov_Transition_rate>,
                                    Derivative<M>, Experiment, X> &dm,
    const M_Matrix<double> &new_x, std::size_t i, double eps) {
  return Markov_Model_calculations<Derivative<Markov_Transition_step>,
                                   Derivative<Markov_Transition_rate>,
                                   Derivative<M>, Experiment, X>(
      Directional_Derivative(dm.Model(), new_x, i, eps), dm.get_Experiment(),
      dm.n_sub_samples(), dm.tolerance());
}

template <class F, class ModelCalculationsDerivative, class... Ds>
auto Incremental_ratio_model(double eps, const F &fun,
                             const ModelCalculationsDerivative &dm,
                             const Ds &... ys) {

  auto m = Primitive(dm);
  auto const &x = dm.Model().x();
  auto f = std::invoke(fun, m, Primitive(ys)...);
  std::vector<std::decay_t<decltype(f)>> fpos(x.size());
  std::vector<std::decay_t<decltype(f)>> fneg(x.size());
  std::vector<double> e(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    e[i] = eps;
    if (x[i] != 0)
      e[i] *= std::abs(x[i]);
    auto mpos = Taylor_first(dm, i, e[i]);
    auto mneg = Taylor_first(dm, i, -e[i]);

    fpos[i] = std::invoke(fun, mpos, Taylor_first(ys, i, e[i])...);
    fneg[i] = std::invoke(fun, mneg, Taylor_first(ys, i, -e[i])...);
  }
  return Incremental_ratio_calc(e, f, x, fpos, fneg);
}

#endif // QMODEL_DERIVATIVE_H
