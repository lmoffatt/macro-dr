#ifndef COMMANDS_H
#define COMMANDS_H


#include "qmodel.h"

#include "Experiment.h"
#include "myDistributions.h"

#include "measure_markov_process.h"
#include <type_traits>
#include "myfields.h"
#include "qmodel.h"
#include "myparameters.h"

using namespace experiment;
struct get_current
{
    template <class Model, class X>
    double operator()(markov_process<std::size_t> mp, const Model& m, const X& x)const
    {
        return (mp.N()*m.g(x));
    }
};


struct model_tag;





struct to_experiment
{
    static constexpr auto className=my_static_string("to_experiment");

   static auto run(const io::myDataFrame<double>& da
                                            ,const std::string& colname_time,
                                            const std::string& colname_nsample,
                                            const std::string& colname_x,
                                            const std::string& colname_y,
                                            double frequency_of_sampling)
    {
       return experiment::DataFrame_to_Experiment(da,colname_time,colname_nsample,colname_x,colname_y,frequency_of_sampling);
    }


    static auto get_arguments()
    {
        return std::make_tuple(grammar::argument(C<const io::myDataFrame<double>&>{},"data_frame"),
                               grammar::argument(C<std::string>{},"colname_time"),
                               grammar::argument(C<std::string>{},"colname_nsample"),
                               grammar::argument(C<std::string>{},"colname_x"),
                               grammar::argument(C<std::string>{},"colname_y"),
                    grammar::argument(C<double>{},"frequency_of_sampling"));
    }


};







template<class Experiment,class Model>
struct simulate{

    typedef Cs<experiment_tag,model_tag> template_tags;


    static constexpr auto className=my_static_string("simulate");
    static auto run(std::mt19937_64::result_type initseed,const Experiment& e, const  Model& m , const Parameters_values<Model>& p,std::size_t N,std::size_t n)
    {
        SingleLigandModel SM(m.Qs(p),m.g(p));

        Markov_Model_calculations<SingleLigandModel,Experiment,double> MC(SM,N,e);

        if (initseed==0)
        {
            std::random_device rd;
            initseed=rd();
        }
        std::mt19937_64 mt(initseed);
        return markov::measure_experiment(get_current{},mt,MC,e,n);
    }

    static auto get_arguments()
    {
        return std::make_tuple(
                    grammar::argument(C<typename std::mt19937_64::result_type>{},"initseed",typename std::mt19937_64::result_type(0)),
                    grammar::argument(C<const Experiment&>{},my_trait<Experiment>::className.c_str()),
                    grammar::argument(C<const Model&>{},my_trait<Model>::className.c_str()),
                    grammar::argument(C<const Parameters_values<Model>&>{},"model_parameters"),
                    grammar::argument(C<std::size_t>{},"number_channels"),
                    grammar::argument(C<std::size_t>{},"number_of_sub_intervals",10ul));
    }
};




typedef typename experiment::basic_Experiment<point<double,double>> singleLigandExperiment;

template<>
struct my_trait < singleLigandExperiment >
{
    static constexpr  auto className=my_static_string("singleLigandExperiment");

};

struct Objects
{

    typedef Cs<Allosteric_Model,singleLigandExperiment> types;
    typedef Cs<simulate<singleLigandExperiment,Allosteric_Model>, to_experiment> commands;
};








#endif // COMMANDS_H
