#ifndef QLIKELIHOOD_H
#define QLIKELIHOOD_H
#include "qmodel.h"
#include "likelihood_markov_process.h"
#include "measure_markov_process.h"
#include "myevidence.h"
#include "mytests.h"


template<class Model, class Parameters_distribution>
class Markov_Model_Likelihood
{
public:
    typedef  Markov_Model_Likelihood self_type;
    typedef  Cs<Model> template_types;
    typedef Parameters_distribution Parameters_Distribution;
    constexpr static auto const className=my_static_string("Markov_Model_Likelihood")+my_trait<template_types>::className;


    template< class Experiment>
    auto compute_Distribution(const Experiment& e, const M_Matrix<double>& parameters) const
    {
        typedef myOptional_t<std::vector<Normal_Distribution<double>>> Op;
        auto p=p_.tr_to_Parameter(parameters);
        SingleLigandModel SM(m.Qs(p),m.g(p),e.Vm(), p.at(Number_of_Channels_Parameter_label()), p.at(gaussian_noise_Parameter_label()), min_P_);
        Markov_Model_calculations<Markov_Transition_step_double,SingleLigandModel,Experiment,double> MC(SM,e,1,tolerance_);
        if(algorithm_=="MacroDR")
        {
        auto out= markov::partialDistribution(markov::MacroDR(tolerance_),MC,e);
        return out;
        }
        else return Op(false,"algoritm "+algorithm_+ " not found");
    };

    template< class Experiment>
    auto compute_Distribution_aux(const Experiment& e, const M_Matrix<double>& parameters, std::ostream& os) const
    {
        typedef myOptional_t<std::pair<std::vector<Normal_Distribution<double>>,std::vector<Op_void>>> Op;
        auto p=p_.tr_to_Parameter(parameters);
        SingleLigandModel SM(m.Qs(p),m.g(p),e.Vm(), p.at(Number_of_Channels_Parameter_label()), p.at(gaussian_noise_Parameter_label()), min_P_);
        Markov_Model_calculations<Markov_Transition_step_double,SingleLigandModel,Experiment,double> MC(SM,e,1,tolerance_);
        if(algorithm_=="MacroDR")
        {
        auto out= markov::partialDistribution_aux(markov::MacroDR(tolerance_),MC,e,os);
        return out;
        }
        else return Op(false,"algoritm "+algorithm_+ " not found");
    };


    template< class Experiment>
    auto compute_Distribution_aux(const Experiment& e, const M_Matrix<double>& parameters, std::ostream& os,const std::vector<Op_void> & aux) const
    {
        typedef myOptional_t<std::vector<Normal_Distribution<double>>> Op;
        auto p=p_.tr_to_Parameter(parameters);
        SingleLigandModel SM(m.Qs(p),m.g(p),e.Vm(), p.at(Number_of_Channels_Parameter_label()), p.at(gaussian_noise_Parameter_label()), min_P_);
        Markov_Model_calculations<Markov_Transition_step_double,SingleLigandModel,Experiment,double> MC(SM,e,1,tolerance_);
        if(algorithm_=="MacroDR")
        {
        auto out= markov::partialDistribution_aux(markov::MacroDR(tolerance_),MC,e,os,aux);
        return out;
        }
        else return Op(false,"algoritm "+algorithm_+ " not found");
    };



    Markov_Model_Likelihood(const Model& m, const Parameters_distribution& p, const std::string& algorithm, double min_P, double tolerance): m{m},p_{p}, algorithm_{algorithm}, min_P_{min_P}, tolerance_{tolerance}{}
private:
   Model m;
   Parameters_distribution p_;
   std::string algorithm_;
   double min_P_;
   double tolerance_;
};


template<class Model, class Experiment, class Parameters_Distribution>
class Markov_Model_DLikelihood: public evidence::FIM_Model<Markov_Model_Likelihood<Model, Parameters_Distribution>,Experiment>
{
public:
    typedef  Markov_Model_DLikelihood self_type;
    typedef  Cs<Model,Experiment> template_types;
    constexpr static auto const className=my_static_string("Markov_Model_DLikelihood")+my_trait<template_types>::className;


    Markov_Model_DLikelihood(const Model& m, const Parameters_Distribution& p, const Experiment& e, const std::string& algorithm, double eps_G,double min_P, double tolerance)
        :evidence::FIM_Model<Markov_Model_Likelihood<Model,Parameters_Distribution>,Experiment>(Markov_Model_Likelihood<Model,Parameters_Distribution>(m,p,algorithm,min_P,tolerance),e, eps_G){}

    using evidence::FIM_Model<Markov_Model_Likelihood<Model,Parameters_Distribution>,Experiment>::compute_DLikelihood;
    using evidence::FIM_Model<Markov_Model_Likelihood<Model,Parameters_Distribution>,Experiment>::compute_Likelihood;
    using evidence::FIM_Model<Markov_Model_Likelihood<Model,Parameters_Distribution>,Experiment>::set_Data;
    //using evidence::FIM_Model<Markov_Model_Likelihood<Model>,Experiment>::getikelihood;

};








#endif // QLIKELIHOOD_H
