#ifndef COMMANDS_H
#define COMMANDS_H


#include "qmodel.h"

#include "Experiment.h"
#include "myDistributions.h"

#include "measure_markov_process.h"
#include "likelihood_markov_process.h"
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

struct to_DataFrame
{
    static constexpr auto className=my_static_string("to_dataframe");

    static auto run(const basic_Experiment<point<double,double>>& e,
                            const std::string& colname_time="time",
                            const std::string& colname_nsample="nsample",
                            const std::string& colname_x="x",
                            const std::string& colname_y="y")
    {
       return experiment::Experiment_to_DataFrame(e,colname_time,colname_nsample,colname_x,colname_y);
    }


    static auto get_arguments()
    {
        return std::make_tuple(grammar::argument(C< const basic_Experiment<point<double,double>>&>{},"experiment"),
                               grammar::argument(C<std::string>{},"colname_time"),
                               grammar::argument(C<std::string>{},"colname_nsample"),
                               grammar::argument(C<std::string>{},"colname_x"),
                               grammar::argument(C<std::string>{},"colname_y"));
    }


};


using namespace  io;


template<class T>
struct save{
    static constexpr auto className=my_static_string("save");

    template<class type_candidate>
    struct uses
    {
        //typedef decltype (operator<< (std::declval<std::ostream&>(),std::declval<type_candidate const &>())) test;
       // static_assert(!has_global_extractor<type_candidate>::value,"test");
        typedef has_global_extractor<type_candidate> type;
    };


    static std::string run(const T& x, const std::string& filename)
    {
         std::ofstream of;
         of.open(filename.c_str());
         if (!of.good())
             return filename+" could not be created";
         else
         {
             if (!(of<<x))
                 return "variable could not be written";
             else
                 return "written successfully";
         }
    }

    static auto get_arguments()
    {
        return std::make_tuple(
                    grammar::argument(C<const T&>{},"variable"),
                    grammar::argument(C<std::string>{},"filename"));
    }
};

template std::ostream& io::write(std::ostream&, const experiment::basic_Experiment<point<double,double>>&);

template<class T>
struct write_variable{
    static constexpr auto className=my_static_string("write");

    template<class type_candidate>
    struct uses
    {
        //typedef decltype (operator<< (std::declval<std::ostream&>(),std::declval<type_candidate const &>())) test;
       // static_assert(!has_global_extractor<type_candidate>::value,"test");
        typedef is_write_Object<type_candidate> type;

       //static_assert (type::value,"true" );
       // static_assert (!type::value,"false" );

    };


    static std::string run(const T& x, const std::string& filename)
    {
         std::ofstream of;
         of.open(filename.c_str());
         if (!of.good())
             return filename+" could not be created";
         else
         {
             if (!x.write(of))
                 return "variable could not be written";
             else
                 return "written successfully";
         }
    }

    static auto get_arguments()
    {
        return std::make_tuple(
                    grammar::argument(C<const T&>{},"variable"),
                    grammar::argument(C<std::string>{},"filename"));
    }
};




template<class Experiment,class Model>
struct simulate{

    typedef Cs<experiment_tag,model_tag> template_tags;


    static constexpr auto className=my_static_string("simulate");
    static auto run(std::mt19937_64::result_type initseed,const Experiment& e,   Model& m , const Parameters_values<Model>& p,std::size_t n)
    {
        SingleLigandModel SM(m.Qs(p),m.g(p),p.at(Number_of_Channels_Parameter_label()),p.at(gaussian_noise_Parameter_label()));

        Markov_Model_calculations<Markov_Transition_step,SingleLigandModel,Experiment,double> MC(SM,e);

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
                    grammar::argument(C< Model&>{},my_trait<Model>::className.c_str()),
                    grammar::argument(C<const Parameters_values<Model>&>{},"model_parameters"),
                    grammar::argument(C<std::size_t>{},"number_of_sub_intervals",10ul));
    }
};




template<class Experiment,class Model>
struct likelihood{


    static constexpr auto className=my_static_string("likelihood");

    static auto run(const Experiment& e,   Model& m , const Parameters_values<Model>& p,const std::string algorithm)
    {
        SingleLigandModel SM(m.Qs(p),m.g(p), p.at(Number_of_Channels_Parameter_label()), p.at(gaussian_noise_Parameter_label()));


        Markov_Model_calculations<Markov_Transition_step_double,SingleLigandModel,Experiment,double> MC(SM,e);

        if (algorithm==markov::MacroDR<true>::className.str())
        {
           auto out= markov::logLikelihood(markov::MacroDR<true>(),MC,e);
           std::cerr<<"logLikelihodd = "<<out<<std::endl;
           return out;

        }
        else return mynan;
    }

    static auto get_arguments()
    {
        return std::make_tuple(
                    grammar::argument(C<const Experiment&>{},my_trait<Experiment>::className.c_str()),
                    grammar::argument(C< Model&>{},my_trait<Model>::className.c_str()),
                    grammar::argument(C<const Parameters_values<Model>&>{},"model_parameters"),
                    grammar::argument(C<std::string>{},"algorithm"));
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
    typedef Cs<simulate<singleLigandExperiment,Allosteric_Model>, likelihood<singleLigandExperiment,Allosteric_Model>,to_experiment, to_DataFrame> commands;
    typedef CCs<save,write_variable> templateCommands;
};








#endif // COMMANDS_H
