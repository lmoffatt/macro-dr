#ifndef QLIKELIHOOD_H
#define QLIKELIHOOD_H
#include "qmodel.h"
#include "likelihood_markov_process.h"
#include "myevidence.h"



template<class Model>
class Markov_Model_Likelihood
{
public:
    template< class Experiment>
    std::vector<Normal_Distribution<double>> getDistribution(const Experiment& e, const M_Matrix<double>& parameters) const
    {
        auto p=p_.tr_to_Parameter(parameters);
        SingleLigandModel SM(m.Qs(p),m.g(p), p.at(Number_of_Channels_Parameter_label()), p.at(gaussian_noise_Parameter_label()));
        Markov_Model_calculations<Markov_Transition_step_double,SingleLigandModel,Experiment,double> MC(SM,e);
        if(algorithm_=="macroDR")
        {
        auto out= markov::partialDistribution(markov::MacroDR<true>(),MC,e);
        return out;
        }
        else return {};
    };
    Markov_Model_Likelihood(const Model& m, Parameters_distribution<Model> p, const std::string algorithm): m{m},p_{p}, algorithm_{algorithm}{}
private:
   Model m;
   Parameters_distribution<Model> p_;
   std::string algorithm_;
};


template<class Model, class Experiment>
class Markov_Model_DLikelihood: public evidence::FIM_Model<Markov_Model_Likelihood<Model>,Experiment>
{
public:
    Markov_Model_DLikelihood(const Model& m, const Parameters_distribution<Model>& p, const Experiment& e, const std::string& algorithm, double eps)
        :evidence::FIM_Model<Markov_Model_Likelihood<Model>,Experiment>(Markov_Model_Likelihood<Model>(m,p,algorithm),e,eps){}

    using evidence::FIM_Model<Markov_Model_Likelihood<Model>,Experiment>::getDLikelihood;
    using evidence::FIM_Model<Markov_Model_Likelihood<Model>,Experiment>::getLikelihood;
    //using evidence::FIM_Model<Markov_Model_Likelihood<Model>,Experiment>::getikelihood;

};



#endif // QLIKELIHOOD_H
