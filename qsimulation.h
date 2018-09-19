#ifndef QSIMULATION_H
#define QSIMULATION_H

#include "qmodel.h"
#include "measure_markov_process.h"

struct get_current
{
    template <class Model, class X>
    double operator()(markov_process<M_Matrix<std::size_t>> mp, const Model& m, const X& x)const
    {
        return (mp.N()*m.g(x)).getvalue();
    }


};


template <class> class Simulator;

template<>
class Simulator<SingleLigandModel>
{
public:
    template<class Model,class Experiment, class Parameters>
auto compute_simulation(const Model& m ,const Experiment& e,   const Parameters& p,std::mt19937_64& mt)
   {
       SingleLigandModel SM(m.Qs(p),m.g(p),e.Vm(),p.at(Number_of_Channels_Parameter_label()),p.at(gaussian_noise_Parameter_label()), min_P_);

       Markov_Model_calculations<Markov_Transition_step,SingleLigandModel,Experiment,double> MC(SM,e,n_sub_intervals_,tolerance_);

       return markov::measure_experiment(get_current{},mt,MC,e,n_sub_intervals_);
   }

Simulator(std::size_t n_sub_intervals, double min_P, double tolerance)
    : n_sub_intervals_{n_sub_intervals},min_P_{min_P},tolerance_{tolerance}
{}

Simulator()=default;

private:
   std::size_t n_sub_intervals_;
   double min_P_;
   double tolerance_;

};
#endif // QSIMULATION_H
