#ifndef SIMULATION_H
#define SIMULATION_H

#include "Markov.h"
#include "Experiment.h"

#include <type_traits>
namespace Simulation {

using namespace experiment;

         template <class Model>
        struct get_current
        {
            static double operator()(markov_process<std::size_t> mp, const Model& m)
            {
                return mp.N()*m.g();
            }
        };


        template<class Model>
        auto Simulate(std::mt19937_64& mt,const Experiment<double,double>& e, const Model m , std::size_t n)
        {
            auto points= measure_experiment(get_current{},mt,m,e,n);
            return e.set_Points(std::move(points),e);

        }

} //namespace Simulation


#endif // SIMULATION_H
