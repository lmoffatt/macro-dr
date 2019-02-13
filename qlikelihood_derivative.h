#ifndef QLIKELIHOOD_DERIVATIVE_H
#define QLIKELIHOOD_DERIVATIVE_H
#include "qmodel.h"
#include "likelihood_markov_process.h"
#include "measure_markov_process.h"
#include "myevidence.h"
#include "mytests.h"
#include "qlikelihood.h"
#include "likelihood_markov_process_derivative.h"
#include "myparameters_derivative.h"



template<class Model, class Parameters_distribution>
class Derivative<Markov_Model_Likelihood<Model,Parameters_distribution>>

{
public:
    typedef  Derivative<Markov_Model_Likelihood<Model,Parameters_distribution>> self_type;
    typedef  Cs<Model> template_types;
    typedef Parameters_distribution Parameters_Distribution;
    typedef M_Matrix<double> Parameters;
    constexpr static auto const className=my_static_string("Markov_Model_Likelihood")+my_trait<template_types>::className;

    typedef markov::MACROR aux_type;

    template< class Experiment>
    auto compute_Distribution(const Experiment& e, const M_Matrix<double>& parameters) const
    {
        typedef myOptional_t<std::vector<Derivative<Normal_Distribution<double>>>> Op;
        auto p_opt=p_.tr_to_Parameter(parameters);
        if (!p_opt)
          return Op(false, "compute_Distribution_aux error in Parameter conversion: "+p_opt.error());
        else{
          auto p=std::move(p_opt).value();
        Derivative<SingleLigandModel> SM(m.Qs(p),m.g(p),e.Vm(), p.at(Number_of_Channels_Parameter_label()), p.at(gaussian_noise_Parameter_label()), min_P_);
        Markov_Model_calculations<Derivative<Markov_Transition_step_double>,Derivative<Markov_Transition_rate>,
                Derivative<SingleLigandModel>,Experiment,double> MC(SM,e,1,tolerance_);
          if(algorithm_==my_trait<markov::MacroDVR>::className.str())
        {
            auto out= markov::partialDistribution_derivative(Derivative<markov::MacroDVR>(tolerance_),MC,e);
            return out;
        }
        else  if(algorithm_==my_trait<markov::MacroDMR>::className.str())
        {
            auto out= markov::partialDistribution_derivative(Derivative<markov::MacroDMR>(tolerance_),MC,e);
            return out;
        }
        else  if(algorithm_==my_trait<markov::MacroDVNR>::className.str())
        {
            auto out= markov::partialDistribution_derivative(Derivative<markov::MacroDVNR>(tolerance_),MC,e);
            return out;
        }

        else  if(algorithm_==my_trait<markov::MacroDMNR>::className.str())
        {
            auto out= markov::partialDistribution_derivative(Derivative<markov::MacroDMNR>(tolerance_),MC,e);
            return out;
        }

        else
            return Op(false,"algoritm "+algorithm_+ " not found");
        }
    };



        template <class Experiment>
    auto compute_DLikelihood(const M_Matrix<double> &parameters,std::ostream& os,const Experiment &e
                                                                                                                           ) const {
      typedef myOptional_t<evidence::PartialDLogLikelihood<markov::MACROR>> Op;
      auto p_opt=p_.tr_to_Parameter_derivative(parameters);
      if (!p_opt)
        return Op(false, "compute_Distribution_aux error in Parameter conversion: "+p_opt.error());
      else{
        auto p=p_opt.value();


        Derivative<SingleLigandModel> SM(
            m.Qs(p), m.g(p), e.Vm(), p.at(Number_of_Channels_Parameter_label()),
            p.at(gaussian_noise_Parameter_label()), min_P_);

        Markov_Model_calculations<Derivative<Markov_Transition_step_double>,
                                Derivative<Markov_Transition_rate>,
                                  Derivative<SingleLigandModel>, Experiment,
                                  double>
            MC(SM, e, 1, tolerance_);
    assert((Derivative_correctness_mean_value_test(
        1e-2,1e4,
        [this,&e](auto p) {
            return Derivative<SingleLigandModel>( m.Qs(p), m.g(p), e.Vm(), p.at(Number_of_Channels_Parameter_label()),
                                                   p.at(gaussian_noise_Parameter_label()), min_P_);
            },
            p,SM, std::cerr,
            "  tolerance=", tolerance_)));



    if (algorithm_ == my_trait<markov::MacroDVR>::className.str())
    {
        return markov::partialDlikelihood_derivative<markov::MACROR,evidence::PartialDLogLikelihood<markov::MACROR>>(
          Derivative<markov::MacroDVR>(tolerance_,BiNumber_,VaNumber_), MC, e,os);
    }
    else  if (algorithm_ == my_trait<markov::MacroDMR>::className.str()) {
        return markov::partialDlikelihood_derivative<
          markov::MACROR, evidence::PartialDLogLikelihood<markov::MACROR>>(
          Derivative<markov::MacroDMR>(tolerance_,BiNumber_,VaNumber_), MC, e,os);
    } else if (algorithm_ == my_trait<markov::MacroDMNR>::className.str()) {
        return markov::partialDlikelihood_derivative<
          markov::MACROR, evidence::PartialDLogLikelihood<markov::MACROR>>(
          Derivative<markov::MacroDMNR>(tolerance_,BiNumber_,VaNumber_), MC, e,os);
    }

    else if (algorithm_ == my_trait<markov::MacroDVNR>::className.str()) {
        return markov::partialDlikelihood_derivative<
          markov::MACROR, evidence::PartialDLogLikelihood<markov::MACROR>>(
          Derivative<markov::MacroDVNR>(tolerance_,BiNumber_,VaNumber_), MC, e,os);
    } else
      return Op(false, "algoritm " + algorithm_ + " not found");
      }
    };




    template <class Experiment>
    auto compute_PartialDLikelihood(const M_Matrix<double> &parameters,std::ostream& os,const Experiment &e
                             ) const {
        typedef myOptional_t<evidence::PartialDLogLikelihood<markov::MACROR,
                                                                        Derivative_t<markov::mp_state_information>>> Op;
      auto p_opt=p_.tr_to_Parameter_derivative(parameters);
      if (!p_opt)
        return Op(false, "compute_Distribution_aux error in Parameter conversion: "+p_opt.error());
      else{
        auto p=p_opt.value();


      Derivative<SingleLigandModel> SM(
          m.Qs(p), m.g(p), e.Vm(), p.at(Number_of_Channels_Parameter_label()),
          p.at(gaussian_noise_Parameter_label()), min_P_);

      Markov_Model_calculations<Derivative<Markov_Transition_step_double>,
                                Derivative<Markov_Transition_rate>,
                                Derivative<SingleLigandModel>, Experiment,
                                double>
          MC(SM, e, 1, tolerance_);
    assert((Derivative_correctness_mean_value_test(
        1e-2,1e4,
        [this,&e](auto p) {
            return Derivative<SingleLigandModel>( m.Qs(p), m.g(p), e.Vm(), p.at(Number_of_Channels_Parameter_label()),
                                                 p.at(gaussian_noise_Parameter_label()), min_P_);
          },
          p,SM, std::cerr,
          "  tolerance=", tolerance_)));



      if (algorithm_ == my_trait<markov::MacroDVR>::className.str())
      {
        return markov::partialDlikelihood_derivative_mp
              <markov::MACROR,evidence::PartialDLogLikelihood<markov::MACROR, Derivative<markov::mp_state_information>>>(
            Derivative<markov::MacroDVR>(tolerance_,BiNumber_,VaNumber_), MC, e,os);
      }
      else  if (algorithm_ == my_trait<markov::MacroDMR>::className.str()) {
          return markov::partialDlikelihood_derivative_mp
              <markov::MACROR,evidence::PartialDLogLikelihood<markov::MACROR, Derivative<markov::mp_state_information>>>(
            Derivative<markov::MacroDMR>(tolerance_,BiNumber_,VaNumber_), MC, e,os);
      } else if (algorithm_ == my_trait<markov::MacroDMNR>::className.str()) {
          return markov::partialDlikelihood_derivative_mp
              <markov::MACROR,evidence::PartialDLogLikelihood<markov::MACROR, Derivative<markov::mp_state_information>>>(
            Derivative<markov::MacroDMNR>(tolerance_,BiNumber_,VaNumber_), MC, e,os);
      }

      else if (algorithm_ == my_trait<markov::MacroDVNR>::className.str()) {
          return markov::partialDlikelihood_derivative_mp
              <markov::MACROR,evidence::PartialDLogLikelihood<markov::MACROR, Derivative<markov::mp_state_information>>>(
            Derivative<markov::MacroDVNR>(tolerance_,BiNumber_,VaNumber_), MC, e,os);
      } else
        return Op(false, "algoritm " + algorithm_ + " not found");
      }
    }

    Derivative(const Derivative<Model>& m, const Derivative<Parameters_distribution>& p, const std::string& algorithm, double min_P, double tolerance, double BiNumber,double VaNumber): m{m},p_{p}, algorithm_{algorithm}, min_P_{min_P}, tolerance_{tolerance}, BiNumber_(BiNumber),VaNumber_{VaNumber}{}
    Derivative()=default;

    auto& get_Model()const {return m;}
    auto& get_Parameters_Distribution()const {return p_;}
    auto get_Algorithm()const { return algorithm_;}
    auto min_P()const {return min_P_;}
    double tolerance()const {return tolerance_;}
    double BiNumber()const {return BiNumber_;}
    double VaNumber()const {return VaNumber_;}


    static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "Model", &self_type::get_Model),
          grammar::field(C<self_type>{}, "Parameters_Distribution", &self_type::get_Parameters_Distribution),
          grammar::field(C<self_type>{}, "Algorithm", &self_type::get_Algorithm),
          grammar::field(C<self_type>{}, "min_P", &self_type::min_P),
          grammar::field(C<self_type>{}, "tolerance", &self_type::tolerance),
          grammar::field(C<self_type>{}, "BiNumber", &self_type::BiNumber),
          grammar::field(C<self_type>{}, "VaNumber", &self_type::VaNumber)
                                                  );
    }

  private:
    Derivative<Model> m;
    Derivative<Parameters_distribution> p_;
    std::string algorithm_;
    double min_P_;
    double tolerance_;
    double BiNumber_;
    double VaNumber_;
};


template<class Model, class Parameters_distribution, class Experiment>
class Derivative<Markov_Model_Likelihood<Model,Parameters_distribution,Experiment>>
{
public:
  typedef  Derivative self_type;
  typedef  Cs<Model> template_types;
  typedef Parameters_distribution Parameters_Distribution;
  typedef M_Matrix<double> Parameters;
  constexpr static auto const className=my_trait< Derivative<Markov_Model_Likelihood<Model,Parameters_distribution>>>::className+
                                          my_trait<Experiment>::className;


  Derivative(const Experiment& e,Derivative<Markov_Model_Likelihood<Model,Parameters_distribution>>&& d
                        ):d_{std::move(d)},e_{e}{}
  Derivative(const Experiment& e,const Derivative<Markov_Model_Likelihood<Model,Parameters_distribution>>& d
                                                              ):d_{d},e_{e}{}




  Derivative()=default;
  template <class mySample>
  auto compute_DLikelihood(const Parameters &p, const mySample &,
                               std::ostream &os) const
  {
    return d_.compute_DLikelihood( p, os,e_);

  }


  auto& get_ModelLikelihood()const {return d_;}
  auto& get_Experiment()const {return e_;}


  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "Experiment", &self_type::get_Experiment),
          grammar::field(C<self_type>{}, "Model_DLikelihood", &self_type::get_ModelLikelihood));
  }


      auto compute_DLikelihood_init(const Parameters &p,
                          std::ostream &os) const
  {
    return d_.compute_DLikelihood( p, os,e_);

  }
  template<class State>
  auto compute_PartialDLikelihood(const M_Matrix<double> &parameters,const State&,std::ostream& os) const
  {
      return d_.compute_PartialDLikelihood(parameters,os,e_);
  }
private:
  Derivative<Markov_Model_Likelihood<Model,Parameters_distribution>> d_;
  Experiment e_;
};


template<class Model, class Parameters_distribution, class Experiment>
template<class Sample>
auto Markov_Model_Likelihood<Model,Parameters_distribution,Experiment>::compute_PartialDLikelihood(const M_Matrix<double>& p,const Sample& current,std::ostream & os) const
{
    auto dm=Derivative<Model>(this->get_Model().get_Model());
    auto dprior=Derivative<Parameters_distribution>(this->get_Model().get_ParametersDistribution());

    Derivative<Markov_Model_Likelihood<Model,Parameters_distribution,Experiment>>
        lik(this->get_Experiment(),
            Derivative<Markov_Model_Likelihood<Model,Parameters_distribution>>(dm,dprior,this->get_Model().get_Algorithm(),this->get_Model().min_P(), this->get_Model().tolerance(),this->get_Model().BiNumber(),this->get_Model().VaNumber()));
    return lik.compute_PartialDLikelihood(p,current,os);
}




#endif // QLIKELIHOOD_DERIVATIVE_H
