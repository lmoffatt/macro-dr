#ifndef QLIKELIHOOD_DERIVATIVE_H
#define QLIKELIHOOD_DERIVATIVE_H
#include "qmodel.h"
#include "likelihood_markov_process.h"
#include "measure_markov_process.h"
#include "myevidence.h"
#include "mytests.h"
#include "qlikelihood.h"
#include "likelihood_markov_process_derivative.h"



template<class Model, class Parameters_distribution>
class Derivative<Markov_Model_Likelihood<Model,Parameters_distribution>>

{
public:
    typedef  Derivative<Markov_Model_Likelihood<Model,Parameters_distribution>> self_type;
    typedef  Cs<Model> template_types;
    typedef Parameters_distribution Parameters_Distribution;
    constexpr static auto const className=my_static_string("Markov_Model_Likelihood")+my_trait<template_types>::className;

    typedef markov::MACROR aux_type;

    template< class Experiment>
    auto compute_Distribution(const Experiment& e, const M_Matrix<double>& parameters) const
    {
        typedef myOptional_t<std::vector<Normal_Distribution<double>>> Op;
        auto p=p_.tr_to_Parameter_derivative(parameters);
        SingleLigandModel SM(m.Qs(p),m.g(p),e.Vm(), p.at(Number_of_Channels_Parameter_label()), p.at(gaussian_noise_Parameter_label()), min_P_);
        Markov_Model_calculations<Markov_Transition_step_double,SingleLigandModel,Experiment,double> MC(SM,e,1,tolerance_);
          if(algorithm_==my_trait<markov::MacroDVR>::className.str())
        {
            auto out= markov::partialDistribution(markov::MacroDVR(tolerance_),MC,e);
            return out;
        }
        else  if(algorithm_==my_trait<markov::MacroDMR>::className.str())
        {
            auto out= markov::partialDistribution(markov::MacroDMR(tolerance_),MC,e);
            return out;
        }
        else  if(algorithm_==my_trait<markov::MacroDVNR>::className.str())
        {
            auto out= markov::partialDistribution(markov::MacroDVNR(tolerance_),MC,e);
            return out;
        }

        else  if(algorithm_==my_trait<markov::MacroDMNR>::className.str())
        {
            auto out= markov::partialDistribution(markov::MacroDMNR(tolerance_),MC,e);
            return out;
        }

        else
            return Op(false,"algoritm "+algorithm_+ " not found");
    };

    template< class Experiment>
    auto compute_Distribution_aux(const Experiment& e, const M_Matrix<double>& parameters, std::ostream& os) const
    {
        typedef myOptional_t<std::pair<std::vector<Normal_Distribution<double>>,std::vector<markov::MACROR>>> Op;
        auto p=p_.tr_to_Parameter(parameters);
        SingleLigandModel SM(m.Qs(p),m.g(p),e.Vm(), p.at(Number_of_Channels_Parameter_label()), p.at(gaussian_noise_Parameter_label()), min_P_);
        Markov_Model_calculations<Markov_Transition_step_double,SingleLigandModel,Experiment,double> MC(SM,e,1,tolerance_);
          if(algorithm_==my_trait<markov::MacroDVR>::className.str())
        {
            auto out= markov::partialDistribution_aux(markov::MacroDVR(tolerance_,BiNumber_,VaNumber_),MC,e,os);
            return out;
        }
          else  if(algorithm_==my_trait<markov::MacroDMR>::className.str())
          {
              auto out= markov::partialDistribution_aux(markov::MacroDMR(tolerance_,BiNumber_),MC,e,os);
              return out;
          }
          else  if(algorithm_==my_trait<markov::MacroDVNR>::className.str())
          {
              auto out= markov::partialDistribution_aux(markov::MacroDVNR(tolerance_,VaNumber_),MC,e,os);
              return out;
          }
          else  if(algorithm_==my_trait<markov::MacroDMNR>::className.str())
          {
              auto out= markov::partialDistribution_aux(markov::MacroDMNR(tolerance_),MC,e,os);
              return out;
          }

        else return Op(false,"algoritm "+algorithm_+ " not found");
    };


    template< class Experiment>
    auto compute_Distribution_aux(const Experiment& e, const M_Matrix<double>& parameters, std::ostream& os,const std::vector<markov::MACROR> & aux) const
    {
        typedef myOptional_t<std::vector<Normal_Distribution<double>>> Op;
        auto p=p_.tr_to_Parameter(parameters);
        SingleLigandModel SM(m.Qs(p),m.g(p),e.Vm(), p.at(Number_of_Channels_Parameter_label()), p.at(gaussian_noise_Parameter_label()), min_P_);
        Markov_Model_calculations<Markov_Transition_step_double,SingleLigandModel,Experiment,double> MC(SM,e,1,tolerance_);
       if(algorithm_==my_trait<markov::MacroDVR>::className.str())
       {
           auto out= markov::partialDistribution_aux(markov::MacroDVR(tolerance_,BiNumber_,VaNumber_),MC,e,os,aux);
           return out;
       }
       else  if(algorithm_==my_trait<markov::MacroDMR>::className.str())
       {
           auto out= markov::partialDistribution_aux(markov::MacroDMR(tolerance_,BiNumber_),MC,e,os,aux);
           return out;
       }
       else  if(algorithm_==my_trait<markov::MacroDVNR>::className.str())
       {
           auto out= markov::partialDistribution_aux(markov::MacroDVNR(tolerance_,VaNumber_),MC,e,os,aux);
           return out;
       }
       else  if(algorithm_==my_trait<markov::MacroDMNR>::className.str())
       {
           auto out= markov::partialDistribution_aux(markov::MacroDMNR(tolerance_),MC,e,os,aux);
           return out;
       }


        else return Op(false,"algoritm "+algorithm_+ " not found");
    };



    Derivative(const Derivative<Model>& m, const Derivative<Parameters_distribution>& p, const std::string& algorithm, double min_P, double tolerance, double BiNumber,double VaNumber): m{m},p_{p}, algorithm_{algorithm}, min_P_{min_P}, tolerance_{tolerance}, BiNumber_(BiNumber),VaNumber_{VaNumber}{}
private:
    Derivative<Model> m;
    Derivative<Parameters_distribution> p_;
    std::string algorithm_;
    double min_P_;
    double tolerance_;
    double BiNumber_;
    double VaNumber_;
};



#endif // QLIKELIHOOD_DERIVATIVE_H
