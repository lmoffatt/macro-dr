#ifndef SIMULATION_H
#define SIMULATION_H

#include "Markov.h"
#include "Experiment.h"
#include "myDistributions.h"

#include "measure_markov_process.h"
#include <type_traits>
#include "myfields.h"

 namespace markov{

        struct get_current
        {
            template <class Model, class X>
             double operator()(markov_process<std::size_t> mp, const Model& m, const X& x)const
            {
                 return (mp.N()*m.g(x));
            }
        };


        template<bool real_calc,class Experiment,class Model>
        auto simulate(std::mt19937_64::result_type initseed,const Experiment& e,  Model& m , std::size_t n)
        {
            if (initseed==0)
            {
                std::random_device rd;
                initseed=rd();
            }
            std::mt19937_64 mt(initseed);
            return measure_experiment(get_current{},mt,m,e,n);
        }

        template<class Experiment,class Model>
         auto get_arguments_Simulate()
         {

             return std::make_tuple(
                    grammar::argument(C<typename std::mt19937_64::result_type>{},"initseed",typename std::mt19937_64::result_type(0)),
                         grammar::argument(C<const Experiment&>{},my_trait<Experiment>::name),
                     grammar::argument(C<Model&>{},my_trait<Model>::name),
                         grammar::argument(C<std::size_t>{},"number_of_sub_intervals",10ul));
         }



 }



 //namespace experiment


#endif // SIMULATION_H
