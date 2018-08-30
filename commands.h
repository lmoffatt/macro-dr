#ifndef COMMANDS_H
#define COMMANDS_H


#include "qmodel.h"

#include "Experiment.h"
#include "myDistributions.h"

#include "measure_markov_process.h"
#include "likelihood_markov_process.h"
#include "myfields.h"
#include "qmodel.h"
#include "myparameters.h"
//#include "myoptimization.h"
#include "qlikelihood.h"
#include <type_traits>

//using opt::Template_Tempering_mcmc;
using namespace experiment;
struct get_current
{
    template <class Model, class X>
    double operator()(markov_process<std::size_t> mp, const Model& m, const X& x)const
    {
        return (mp.N()*m.g(x)).getvalue();
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
                    double Vm,
                    double frequency_of_sampling)
    {
        return experiment::DataFrame_to_Experiment(da,colname_time,colname_nsample,colname_x,colname_y,Vm,frequency_of_sampling);
    }


    static auto get_arguments()
    {
        return std::make_tuple(grammar::argument(C<const io::myDataFrame<double>&>{},"data_frame"),
                               grammar::argument(C<std::string>{},"colname_time"),
                               grammar::argument(C<std::string>{},"colname_nsample"),
                               grammar::argument(C<std::string>{},"colname_x"),
                               grammar::argument(C<std::string>{},"colname_y"),
                               grammar::argument(C<double>{},"holding_potential"),
                               grammar::argument(C<double>{},"frequency_of_sampling"));
    }


};

struct function_log10
{
    static constexpr auto className=my_static_string("log10");
    static double run(double x){ return std::log10(x);}
    static   auto get_arguments()
    {
        return std::make_tuple(grammar::argument(C< double>{},"x"));
    }
};



template<class measure>
struct to_DataFrame
{
    static constexpr auto className=my_static_string("to_dataframe");

    static auto run(const basic_Experiment<point<double,double>,measure>& e,
                    const std::string& colname_trace="trace",
                    const std::string& colname_time="time",
                    const std::string& colname_nsample="nsample",
                    const std::string& colname_x="x",
                    const std::string& colname_y="y")
    {
        if constexpr (std::is_same_v<measure,measure_just_y<double >>)
                return experiment::Experiment_to_DataFrame(e,colname_trace,colname_time,colname_nsample,colname_x,colname_y);
        else
        return experiment::Experiment_Steps_to_DataFrame(e,colname_trace,colname_time,colname_nsample,colname_x,colname_y);


    }


    static auto get_arguments()
    {
        return std::make_tuple(grammar::argument(C< const basic_Experiment<point<double,double>,measure>&>{},"experiment"),
                               grammar::argument(C<std::string>{},"colname_trace", std::string("trace")),
                               grammar::argument(C<std::string>{},"colname_time", std::string("time")),
                               grammar::argument(C<std::string>{},"colname_nsample",std::string("nsample")),
                               grammar::argument(C<std::string>{},"colname_x",std::string("x")),
                               grammar::argument(C<std::string>{},"colname_y",std::string("y")));
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

template std::ostream& io::write(std::ostream&, const experiment::basic_Experiment<point<double,double>,measure_just_y<double>>&);

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
    static auto run(std::mt19937_64::result_type initseed,const Experiment& e,  const Model& m , const Parameters_values<Model>& p,std::size_t n_sub_intervals
                    ,double min_P,double tolerance)
    {
        SingleLigandModel SM(m.Qs(p),m.g(p),e.Vm(),p.at(Number_of_Channels_Parameter_label()),p.at(gaussian_noise_Parameter_label()), min_P, tolerance);

        Markov_Model_calculations<Markov_Transition_step,SingleLigandModel,Experiment,double> MC(SM,e,n_sub_intervals);

        if (initseed==0)
        {
            std::random_device rd;
            initseed=rd();
        }
        std::mt19937_64 mt(initseed);
        return markov::measure_experiment(get_current{},mt,MC,e,n_sub_intervals);
    }

    static auto get_arguments()
    {
        return std::make_tuple(
                    grammar::argument(C<typename std::mt19937_64::result_type>{},"initseed",typename std::mt19937_64::result_type(0)),
                    grammar::argument(C<const Experiment&>{},my_trait<Experiment>::className.c_str()),
                    grammar::argument(C<const Model&>{},my_trait<Model>::className.c_str()),
                    grammar::argument(C<const Parameters_values<Model>&>{},"model_parameters"),
                    grammar::argument(C<std::size_t>{},"number_of_sub_intervals",10ul),
                    grammar::argument(C<double>{},"min_probability",1e-9),
                    grammar::argument(C<double>{},"tolerance_error",1e-7));
    }
};




template<class Experiment,class Model>
struct likelihood{


    static constexpr auto className=my_static_string("likelihood");

    static auto run(const Experiment& e,   Model& m , const Parameters_values<Model>& p,const std::string algorithm, double min_P, double tolerance)
    {
        SingleLigandModel SM(m.Qs(p),m.g(p),e.Vm(), p.at(Number_of_Channels_Parameter_label()), p.at(gaussian_noise_Parameter_label()),min_P,tolerance);


        Markov_Model_calculations<Markov_Transition_step_double,SingleLigandModel,Experiment,double> MC(SM,e,1);

        if (algorithm==my_trait<markov::MacroDR<true>>::className.str())
        {
            auto out= markov::logLikelihood(markov::MacroDR<true>(),MC,e);
            std::cerr<<"logLikelihodd = "<<out<<std::endl;
            return out;

        }
        else if (algorithm==my_trait<markov::MacroDR<false>>::className.str())
        {
            auto out= markov::logLikelihood(markov::MacroDR<false>(),MC,e);
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
                    grammar::argument(C<std::string>{},"algorithm"),
                    grammar::argument(C<double>{},"min_probability",1e-9),
                    grammar::argument(C<double>{},"tolerance_error",1e-7));
    }
};



template<class Experiment,class Model>
struct likelihood_detail{


    static constexpr auto className=my_static_string("likelihood_detail");

    static auto run(const Experiment& e,   Model& m , const Parameters_values<Model>& p,const std::string algorithm, double min_P, double tolerance)
    {
        SingleLigandModel SM(m.Qs(p),m.g(p), e.Vm(),p.at(Number_of_Channels_Parameter_label()), p.at(gaussian_noise_Parameter_label()),min_P,tolerance);


        Markov_Model_calculations<Markov_Transition_step_double,SingleLigandModel,Experiment,double> MC(SM,e,1);

        if (algorithm==my_trait<markov::MacroDR<true>>::className.str())
        {
            auto out= markov::monitorLikelihood(markov::MacroDR<true>(),MC,e);
            return out;

        }
        else         if (algorithm==my_trait<markov::MacroDR<false>>::className.str())
        {
            auto out= markov::monitorLikelihood(markov::MacroDR<false>(),MC,e);
            return out;

        }
        else return experiment::basic_Experiment<experiment::point<double,double>,markov::measure_likelihood<double>>{};
    }

    static auto get_arguments()
    {
        return std::make_tuple(
                    grammar::argument(C<const Experiment&>{},my_trait<Experiment>::className.c_str()),
                    grammar::argument(C< Model&>{},my_trait<Model>::className.c_str()),
                    grammar::argument(C<const Parameters_values<Model>&>{},"model_parameters"),
                    grammar::argument(C<std::string>{},"algorithm"),
                    grammar::argument(C<double>{},"min_probability",1e-9),
                    grammar::argument(C<double>{},"tolerance_error",1e-7));
    }
};



typedef typename experiment::basic_Experiment<point<double,double>, measure_just_y<double>> singleLigandExperiment;
typedef typename experiment::basic_Experiment<point<double,double>, markov::measure_likelihood<double>> singleLigandLikelihood;

template<>
struct my_trait < singleLigandExperiment >
{
    static constexpr  auto className=my_static_string("singleLigandExperiment");
};

template<>
struct my_trait < singleLigandLikelihood >
{
    static constexpr  auto className=my_static_string("singleLigandLikelihood");
};



template<class Experiment,class Model>
struct Evidence{


    static constexpr auto className=my_static_string("evidence");

    static auto run(const Experiment& e,
                    const Model& m ,
                    const Parameters_distribution<Model>& p,
                    std::string algorithm,
                    double min_P,
                    double tolerance,
                    std::mt19937_64::result_type initseed,
                    std::vector<double> betas,
                    std::vector<double>landa,
                    std::vector<std::vector<double>>landa_50_hill,
                    double gain_moment,
                    std::size_t nSamples,
                    bool parameters,
                    bool gradient)
    {

        Markov_Model_DLikelihood<Model,Experiment> lik(m,p,e,algorithm,min_P, tolerance);
        if (initseed==0)
        {
            std::random_device rd;
            initseed=rd();
        }
        std::mt19937_64 mt(initseed);
        evidence::OutputGenerator out(std::cerr,parameters,gradient);

        return evidence::run_Thermo_Levenberg_ProbVel(evidence::Prior_Model<Model>(p),lik,mt,betas,landa,landa_50_hill,gain_moment,nSamples,out);
    }

    static auto get_arguments()
    {
        return std::make_tuple(
                    grammar::argument(C<const Experiment&>{},my_trait<Experiment>::className.c_str()),
                    grammar::argument(C< Model&>{},my_trait<Model>::className.c_str()),
                    grammar::argument(C<const Parameters_distribution<Model>&>{},"model_parameters_distribution"),
                    grammar::argument(C<std::string>{},"algorithm"),
                    grammar::argument(C<double>{},"min_probability"),
                    grammar::argument(C<double>{},"tolerance_error"),
                    grammar::argument(C<std::mt19937_64::result_type>{},"initseed"),
                    grammar::argument(C<std::vector<double>>{},"betas"),
                    grammar::argument(C<std::vector<double>>{},"landas"),
                    grammar::argument(C<std::vector<std::vector<double>>>{},"landa_50_hill"),
                    grammar::argument(C<double>{},"gain_moment"),
                    grammar::argument(C<std::size_t>{},"nSamples"),
                    grammar::argument(C<bool>{},"parameters_output"),
                    grammar::argument(C<double>{},"gradient_output")

                    );
    }
};




struct Objects
{

    typedef Cs<Allosteric_Model,singleLigandExperiment> types;
    typedef Cs<
    function_log10,
    simulate<singleLigandExperiment,Allosteric_Model>, likelihood<singleLigandExperiment,Allosteric_Model>,
    likelihood_detail<singleLigandExperiment,Allosteric_Model>,
    Evidence<singleLigandExperiment,Allosteric_Model>,

    to_experiment, to_DataFrame<measure_just_y<double>>, to_DataFrame<markov::measure_likelihood<double>>> commands;
    typedef CCs<save,write_variable> templateCommands;
};








#endif // COMMANDS_H
