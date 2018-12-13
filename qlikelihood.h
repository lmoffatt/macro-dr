#ifndef QLIKELIHOOD_H
#define QLIKELIHOOD_H
#include "qmodel.h"
#include "likelihood_markov_process.h"
#include "measure_markov_process.h"
#include "myevidence.h"
#include "mytests.h"
template<class ...>
class Markov_Model_Likelihood;

template<class Model, class Parameters_distribution>
class Markov_Model_Likelihood<Model,Parameters_distribution>
{
public:
  typedef  Markov_Model_Likelihood self_type;
  typedef  Cs<Model> template_types;
  typedef Parameters_distribution Parameters_Distribution;
  constexpr static auto const className=my_static_string("Markov_Model_Likelihood")+my_trait<template_types>::className;
  typedef M_Matrix<double> Parameters;

  typedef markov::MACROR aux_type;

  template< class Experiment>
  auto compute_Distribution(const Experiment& e, const M_Matrix<double>& parameters) const
  {
    typedef myOptional_t<std::vector<Normal_Distribution<double>>> Op;
    auto p_opt=p_.tr_to_Parameter(parameters);
    if (!p_opt)
      return Op(false, "compute_Distribution_aux error in Parameter conversion: "+p_opt.error());
    else{
      auto p=std::move(p_opt).value();
      SingleLigandModel SM(m.Qs(p),m.g(p),e.Vm(), p.at(Number_of_Channels_Parameter_label()), p.at(gaussian_noise_Parameter_label()), min_P_);
      Markov_Model_calculations<Markov_Transition_step_double,Markov_Transition_rate,SingleLigandModel,Experiment,double> MC(SM,e,1,tolerance_);
      if(algorithm_==my_trait<markov::MacroDVR>::className.str())
      {
        auto out= markov::partialDistribution(markov::MacroDVR(tolerance_,BiNumber_,VaNumber_),MC,e);
        return out;
      }
      else  if(algorithm_==my_trait<markov::MacroDMR>::className.str())
      {
        auto out= markov::partialDistribution(markov::MacroDMR(tolerance_,BiNumber_,VaNumber_),MC,e);
        return out;
      }
      else  if(algorithm_==my_trait<markov::MacroDVNR>::className.str())
      {
        auto out= markov::partialDistribution(markov::MacroDVNR(tolerance_,BiNumber_,VaNumber_),MC,e);
        return out;
      }

      else  if(algorithm_==my_trait<markov::MacroDMNR>::className.str())
      {
        auto out= markov::partialDistribution(markov::MacroDMNR(tolerance_,BiNumber_,VaNumber_),MC,e);
        return out;
      }

      else
        return Op(false,"algoritm "+algorithm_+ " not found");
    }
  };

  template< class Experiment>
  auto compute_Distribution_aux(const Experiment& e, const M_Matrix<double>& parameters, std::ostream& os) const
  {
    typedef myOptional_t<std::pair<std::vector<Normal_Distribution<double>>,std::vector<markov::MACROR>>> Op;
    auto p_opt=p_.tr_to_Parameter(parameters);
    if (!p_opt)
      return Op(false, "compute_Distribution_aux error in Parameter conversion: "+p_opt.error());
    else{
      auto p=std::move(p_opt).value();
      SingleLigandModel SM(m.Qs(p),m.g(p),e.Vm(), p.at(Number_of_Channels_Parameter_label()), p.at(gaussian_noise_Parameter_label()), min_P_);
      Markov_Model_calculations<Markov_Transition_step_double,Markov_Transition_rate,SingleLigandModel,Experiment,double> MC(SM,e,1,tolerance_);
      if(algorithm_==my_trait<markov::MacroDVR>::className.str())
      {
        auto out= markov::partialDistribution_aux(markov::MacroDVR(tolerance_,BiNumber_,VaNumber_),MC,e,os);
        return out;
      }
      else  if(algorithm_==my_trait<markov::MacroDMR>::className.str())
      {
        auto out= markov::partialDistribution_aux(markov::MacroDMR(tolerance_,BiNumber_,VaNumber_),MC,e,os);
        return out;
      }
      else  if(algorithm_==my_trait<markov::MacroDVNR>::className.str())
      {
        auto out= markov::partialDistribution_aux(markov::MacroDVNR(tolerance_,BiNumber_,VaNumber_),MC,e,os);
        return out;
      }
      else  if(algorithm_==my_trait<markov::MacroDMNR>::className.str())
      {
        auto out= markov::partialDistribution_aux(markov::MacroDMNR(tolerance_,BiNumber_,VaNumber_),MC,e,os);
        return out;
      }

      else return Op(false,"algoritm "+algorithm_+ " not found");
    }
  };


  template< class Experiment>
  auto compute_Distribution_aux(const Experiment& e, const M_Matrix<double>& parameters, std::ostream& os,const std::vector<markov::MACROR> & aux) const
  {
    typedef myOptional_t<std::vector<Normal_Distribution<double>>> Op;
    auto p_opt=p_.tr_to_Parameter(parameters);
    if (!p_opt)
      return Op(false, "compute_Distribution_aux error in Parameter conversion: "+p_opt.error());
    else{
      auto p=std::move(p_opt).value();
      SingleLigandModel SM(m.Qs(p),m.g(p),e.Vm(), p.at(Number_of_Channels_Parameter_label()), p.at(gaussian_noise_Parameter_label()), min_P_);
      Markov_Model_calculations<Markov_Transition_step_double,Markov_Transition_rate,SingleLigandModel,Experiment,double> MC(SM,e,1,tolerance_);
      if(algorithm_==my_trait<markov::MacroDVR>::className.str())
      {
        auto out= markov::partialDistribution_aux(markov::MacroDVR(tolerance_,BiNumber_,VaNumber_),MC,e,os,aux);
        return out;
      }
      else  if(algorithm_==my_trait<markov::MacroDMR>::className.str())
      {
        auto out= markov::partialDistribution_aux(markov::MacroDMR(tolerance_,BiNumber_,VaNumber_),MC,e,os,aux);
        return out;
      }
      else  if(algorithm_==my_trait<markov::MacroDVNR>::className.str())
      {
        auto out= markov::partialDistribution_aux(markov::MacroDVNR(tolerance_,BiNumber_,VaNumber_),MC,e,os,aux);
        return out;
      }
      else  if(algorithm_==my_trait<markov::MacroDMNR>::className.str())
      {
        auto out= markov::partialDistribution_aux(markov::MacroDMNR(tolerance_,BiNumber_,VaNumber_),MC,e,os,aux);
        return out;
      }


      else return Op(false,"algoritm "+algorithm_+ " not found");
    }
  };


  template<class Experiment>
  auto compute_Likelihood(const Experiment& e,const M_Matrix<double>& p)const
  {
    typedef myOptional_t<evidence::logLikelihood> Op;
    auto D0 = compute_Distribution(e, p);
    if (D0.has_value()) {
      auto logL = calculate_Likelihood(D0, e);
      std::stringstream ss;
      if (are_finite<true, double>().test(logL, ss))
        return Op(logLikelihood(logL));
      else
        return Op(false, ss.str());
    } else
      return Op(false, "fails to compute model data" + D0.error());

  }



  Markov_Model_Likelihood(const Model& m, const Parameters_distribution& p, const std::string& algorithm, double min_P, double tolerance, double BiNumber,double VaNumber): m{m},p_{p}, algorithm_{algorithm}, min_P_{min_P}, tolerance_{tolerance}, BiNumber_(BiNumber),VaNumber_{VaNumber}{}
private:
  Model m;
  Parameters_distribution p_;
  std::string algorithm_;
  double min_P_;
  double tolerance_;
  double BiNumber_;
  double VaNumber_;
};

template<class Model, class Parameters_distribution, class Experiment>
class Markov_Model_Likelihood<Model,Parameters_distribution, Experiment>
{
public:
  typedef  Markov_Model_Likelihood self_type;
  typedef  Cs<Model, Experiment> template_types;
  typedef Parameters_distribution Parameters_Distribution;
  constexpr static auto const className=my_static_string("Markov_Model_Likelihood")+my_trait<template_types>::className;
  typedef M_Matrix<double> Parameters;

  typedef markov::MACROR aux_type;

  auto compute_Distribution(const Experiment& e, const M_Matrix<double>& parameters) const
  {
    return l_.compute_Distribution(e,parameters);
  };

  auto compute_Distribution_aux(const Experiment& e, const M_Matrix<double>& parameters, std::ostream& os) const
  {
    return l_.compute_Distribution_aux(e,parameters,os);
  }

  auto compute_Distribution_aux(const Experiment& e, const M_Matrix<double>& parameters, std::ostream& os,const std::vector<markov::MACROR> & aux) const
  {
    return compute_Distribution_aux(e,parameters,os,aux);
  }
  auto compute_Likelihood(const M_Matrix<double>& p)const
  {
    typedef myOptional_t<evidence::logLikelihood> Op;
    auto D0 = compute_Distribution(e_, p);
    if (D0.has_value()) {
      auto logL = evidence::calculate_Likelihood(D0.value(), e_);
      std::stringstream ss;
      if (are_finite<true, double>().test(std::get<0>(logL), ss))
        return Op(evidence::logLikelihood(logL));
      else
        return Op(false, ss.str());
    } else
      return Op(false, "fails to compute model data" + D0.error());

  }


  Markov_Model_Likelihood(const Model& m, const Parameters_distribution& p,const Experiment& e, const std::string& algorithm, double min_P, double tolerance, double BiNumber,double VaNumber): e_{e},l_{m,p,algorithm,min_P,tolerance,BiNumber,VaNumber}{}
private:
  Experiment  const &e_;
  Markov_Model_Likelihood<Model,Parameters_distribution> l_;

};










template<class Model, class Experiment, class Parameters_Distribution>
class Markov_Model_DLikelihood: public evidence::FIM_Model<Markov_Model_Likelihood<Model, Parameters_Distribution>,Experiment>
{
public:
  typedef  Markov_Model_DLikelihood self_type;
  typedef  Cs<Model,Experiment> template_types;
  constexpr static auto const className=my_static_string("Markov_Model_DLikelihood")+my_trait<template_types>::className;


    Markov_Model_DLikelihood(const Model& m, const Parameters_Distribution& p, const Experiment& e, const std::string& algorithm, double eps_G,double min_P, double tolerance, double BiNumber, double VaNumber, double epsf)
        :evidence::FIM_Model<Markov_Model_Likelihood<Model,Parameters_Distribution>,Experiment>(Markov_Model_Likelihood<Model,Parameters_Distribution>(m,p,algorithm,min_P,tolerance,BiNumber,VaNumber),e, eps_G, epsf){}

    using evidence::FIM_Model<Markov_Model_Likelihood<Model,Parameters_Distribution>,Experiment>::compute_DLikelihood;
    using evidence::FIM_Model<Markov_Model_Likelihood<Model,Parameters_Distribution>,Experiment>::compute_Likelihood;
    using evidence::FIM_Model<Markov_Model_Likelihood<Model,Parameters_Distribution>,Experiment>::set_Data;
    //using evidence::FIM_Model<Markov_Model_Likelihood<Model>,Experiment>::getikelihood;

};








#endif // QLIKELIHOOD_H
