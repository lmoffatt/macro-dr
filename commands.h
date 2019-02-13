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
#include "qsimulation.h"
#include <type_traits>
#include "likelihood_markov_process_derivative.h"
#include "qlikelihood_derivative.h"
#include "myparameters_derivative.h"
#include "myevidence.h"
//using opt::Template_Tempering_mcmc;
using namespace experiment;


struct model_tag;





struct to_experiment
{
    static constexpr auto className=my_static_string("to_experiment");

    static
        basic_Experiment<point<double,double>,measure_just_y<double>>
        run(const io::myDataFrame<double>& da
                    ,const std::string& colname_time,
                    const std::string& colname_nsample,
                    const std::string& colname_x,
                    const std::string& colname_y,
                    double Vm,
                    double frequency_of_sampling);


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

    static io::myDataFrame<double,std::size_t, std::string> run(const basic_Experiment<point<double,double>,measure>& e);


    static auto get_arguments()
    {
        return std::make_tuple(grammar::argument(C< const basic_Experiment<point<double,double>,measure>&>{},"experiment"));
    }


};
typedef typename experiment::basic_Experiment<point<double,double>, measure_just_y<double>> singleLigandExperiment;
typedef typename experiment::basic_Experiment<point<double,double>, markov::measure_likelihood<double>> singleLigandLikelihood;


//template<class Parameters_Distribution>
//struct to_DataFrame_index
//{
//    static constexpr auto className=my_static_string("to_dataframe");

//    static auto run(const evidence::Likelihood_Test::sample<singleLigandExperiment,markov::MACROR>& l, const Parameters_Distribution& prior)
//    {
//        std::vector<io::myDataFrame<double,std::size_t,std::string>> d(l.getLikelihoods().size());
//        for (auto i=0lu; i<l.getLikelihoods().size(); ++i )
//        {
//            d[i]=experiment::Experiment_steps_to_DataFrame(l.getSimulations()[i],l.getLikelihoods()[i].partial_DlogL(),prior);
//        }
//        return myDataFrame<double,std::size_t,std::string>::consolidate(d,"i_simul");

//    }


//    static auto get_arguments()
//    {
//        return std::make_tuple(
//                    grammar::argument(C< const evidence::Likelihood_Test::sample<singleLigandExperiment>& >{},"samples"),
//                    grammar::argument(C< const Parameters_Distribution& >{},"parameters")

//                    );
//    }


//};

template<class Parameters_Distribution,class Parameters_Values,class ExperimentData, class logLikelihood>
struct to_DataFrame<evidence::Likelihood_Analisis<Parameters_Distribution,Parameters_Values,ExperimentData,logLikelihood>>
{
    static constexpr auto className=my_static_string("to_dataframe");

    static io::myDataFrame<double, std::size_t, std::string> run(const evidence::Likelihood_Analisis<Parameters_Distribution,Parameters_Values,ExperimentData,logLikelihood>& l);


    static auto get_arguments()
    {
        return std::make_tuple(grammar::argument(C< const evidence::Likelihood_Analisis<Parameters_Distribution,Parameters_Values,ExperimentData,logLikelihood>&>{},"analysis"));
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


    static std::string run(const T& x, const std::string& filename) {
        std::ofstream of;
        of.open(filename.c_str());
        if (!of.good())
            return filename + " could not be created";
        else {
            if (!(of << x))
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


    static std::string run(const T& x, const std::string& filename) {
        std::ofstream of;
        of.open(filename.c_str());
        if (!of.good())
            return filename + " could not be created";
        else {
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
    static myOptional_t<Experiment>
 run(std::mt19937_64::result_type initseed,const Experiment& e,  const Model& m , const Parameters_values<Model>& p,std::size_t n_sub_intervals,
                    double max_dt,double min_P,double tolerance);

    static auto get_arguments()
    {
        return std::make_tuple(
            grammar::argument(C<typename std::mt19937_64::result_type>{},"initseed",typename std::mt19937_64::result_type(0)),
            grammar::argument(C<const Experiment&>{},my_trait<Experiment>::className.c_str()),
            grammar::argument(C<const Model&>{},my_trait<Model>::className.c_str()),
            grammar::argument(C<const Parameters_values<Model>&>{},"model_parameters"),
            grammar::argument(C<std::size_t>{},"number_of_sub_intervals",10ul),
            grammar::argument(C<double>{},"max_dt",1e-4),
            grammar::argument(C<double>{},"min_probability",1e-9),
            grammar::argument(C<double>{},"tolerance_error",1e-7));
    }
};




template<class Experiment,class Model>
struct likelihood{


    static constexpr auto className=my_static_string("likelihood");

    static myOptional_t<double> run(const Experiment& e,   Model& m , const Parameters_values<Model>& p,const std::string algorithm, double min_P, double tolerance, double biNumber,double Vanumber);

    static auto get_arguments()
    {
        return std::make_tuple(
            grammar::argument(C<const Experiment&>{},my_trait<Experiment>::className.c_str()),
            grammar::argument(C< Model&>{},my_trait<Model>::className.c_str()),
            grammar::argument(C<const Parameters_values<Model>&>{},"model_parameters"),
            grammar::argument(C<std::string>{},"algorithm"),
            grammar::argument(C<double>{},"min_probability",1e-9),
            grammar::argument(C<double>{},"tolerance_error",1e-7),
            grammar::argument(C<double>{},"Binomial_threshold",5.0),
            grammar::argument(C<double>{},"Variance_threshold",1.0)
                                    );
    }
};


template<class Experiment,class Model, class ParametersDistribution>
struct likelihoodtest{


    static constexpr auto className=my_static_string("likelihoodtest");
    static auto run(std::size_t initseed,
                    const Experiment& e,
                    const Model& m ,
                    const Parameters_values<Model>& p,
                    const ParametersDistribution& prior,
                    const std::string algorithm,
                    double BiNumber,
                    double VaNumber,
                    double min_P,
                    double tolerance,
                    double eps_Gradient,
                    bool eps_adjust,
                    bool center_Gradient,
                    std::size_t n_sub_intervals,
                    double max_dt,
                    std::size_t nsamples,
                    double pvalue,
                    double epsf)
    {
        std::cerr<<"\nparameters\n"<<p;
        std::cerr<<"\n initseed="<<initseed<<"\n";
        std::cerr<<"\n n_sub_intervals="<<n_sub_intervals<<"\n";
        std::cerr<<"\n nsamples="<<nsamples<<"\n";
        std::cerr<<"\n center_Gradient="<<center_Gradient<<"\n";


        std::mt19937_64 mt=init_mt(initseed, std::cerr);
        Markov_Model_DLikelihood<Model,Experiment,ParametersDistribution> lik(m,prior,e,algorithm,eps_Gradient,min_P, tolerance,BiNumber,VaNumber, epsf);
        Simulator<Model> sim(m,n_sub_intervals, max_dt,min_P,tolerance);
        return evidence::Likelihood_Test::compute_test(std::cerr,sim,lik,e,prior,p,mt,nsamples,pvalue,!eps_adjust,center_Gradient);
    }

    static auto get_arguments()
    {   return std::make_tuple(
            grammar::argument(C<std::size_t>{},"initseed"),
            grammar::argument(C<const Experiment&>{},my_trait<Experiment>::className.c_str()),
            grammar::argument(C< Model&>{},my_trait<Model>::className.c_str()),
            grammar::argument(C<const Parameters_values<Model>&>{},"model_parameters"),
            grammar::argument(C<const ParametersDistribution&>{},"model_parameters_distribution"),
            grammar::argument(C<std::string>{},"algorithm"),
            grammar::argument(C<double>{},"Binomial_threshold",5.0),
            grammar::argument(C<double>{},"Variance_threshold",1.0),
            grammar::argument(C<double>{},"min_probability",1e-9),
            grammar::argument(C<double>{},"tolerance_error",1e-7),
            grammar::argument(C<double>{},"eps_G"),
            grammar::argument(C<bool>{},"eps_G_adjust"),
            grammar::argument(C<bool>{},"Center_gradient"),
            grammar::argument(C<std::size_t>{},"number_of_sub_intervals"),
            grammar::argument(C<double>{},"max_dt"),
            grammar::argument(C<std::size_t>{},"nsamples"),
            grammar::argument(C<double>{},"p_value"),
            grammar::argument(C<double>{},"eps_factor")
                                    );
    }

};


template<class Experiment,class Model, class ParametersDistribution>
struct likelihoodtest_derivative{


    static constexpr auto className=my_static_string("likelihoodtest_derivative");
    static auto run(std::size_t initseed,
                    const Experiment& e,
                    const Model& m ,
                    const Parameters_values<Model>& p,
                    const ParametersDistribution& prior,
                    const std::string algorithm,
                    double BiNumber,
                    double VaNumber,
                    double min_P,
                    double tolerance,
                    std::size_t n_sub_intervals,
                    double max_dt,
                    std::size_t nsamples,
                    double pvalue)
    {
        std::cerr<<"\nparameters\n"<<p;
        std::cerr<<"\n initseed="<<initseed<<"\n";
        std::cerr<<"\n n_sub_intervals="<<n_sub_intervals<<"\n";
        std::cerr<<"\n nsamples="<<nsamples<<"\n";


        std::mt19937_64 mt=init_mt(initseed, std::cerr);
        auto dm=Derivative<Model>(m);
        auto dprior=Derivative<ParametersDistribution>(prior);
        Derivative<Markov_Model_Likelihood<Model,ParametersDistribution>> lik(dm,dprior,algorithm,min_P, tolerance,BiNumber,VaNumber);
        Simulator<Model> sim(m,n_sub_intervals, max_dt,min_P,tolerance);
        return evidence::Likelihood_Test::compute_test_derivative(std::cerr,sim,lik,e,dprior,p,mt,nsamples,pvalue);
    }

    static auto get_arguments()
    {   return std::make_tuple(
            grammar::argument(C<std::size_t>{},"initseed"),
            grammar::argument(C<const Experiment&>{},my_trait<Experiment>::className.c_str()),
            grammar::argument(C< Model&>{},my_trait<Model>::className.c_str()),
            grammar::argument(C<const Parameters_values<Model>&>{},"model_parameters"),
            grammar::argument(C<const ParametersDistribution&>{},"model_parameters_distribution"),
            grammar::argument(C<std::string>{},"algorithm"),
            grammar::argument(C<double>{},"Binomial_threshold",5.0),
            grammar::argument(C<double>{},"Variance_threshold",1.0),
            grammar::argument(C<double>{},"min_probability",1e-9),
            grammar::argument(C<double>{},"tolerance_error",1e-7),
            grammar::argument(C<std::size_t>{},"number_of_sub_intervals"),
            grammar::argument(C<double>{},"max_dt"),
            grammar::argument(C<std::size_t>{},"nsamples"),
            grammar::argument(C<double>{},"p_value")
                                    );
    }

};



template<class Experiment,class Model>
struct likelihood_detail{


    static constexpr auto className=my_static_string("likelihood_detail");

    static auto run(const Experiment& e,   Model& m , const Parameters_values<Model>& p,const std::string algorithm, double min_P, double tolerance, double biNumber,double Vanumber)
    {
        typedef myOptional_t<experiment::basic_Experiment<experiment::point<double,double>,markov::measure_likelihood<double>>> Op;

        SingleLigandModel SM(m.Qs(p),m.g(p), e.Vm(),p.at(Number_of_Channels_Parameter_label()), p.at(gaussian_noise_Parameter_label()),min_P);


        Markov_Model_calculations<Markov_Transition_step_double,Markov_Transition_rate,SingleLigandModel,Experiment,double> MC(SM,e,1,tolerance);

        if (algorithm==my_trait<markov::MacroDVR>::className.str())
        {
            return markov::monitorLikelihood(markov::MacroDVR(tolerance,biNumber,Vanumber),MC,e, std::cerr);
        }
        else if (algorithm==my_trait<markov::MacroDMR>::className.str())
        {
            return markov::monitorLikelihood(markov::MacroDMR(tolerance,biNumber,Vanumber),MC,e, std::cerr);
        }
        else if (algorithm==my_trait<markov::MacroDVNR>::className.str())
        {
            return markov::monitorLikelihood(markov::MacroDVNR(tolerance,biNumber,Vanumber),MC,e, std::cerr);
        }
        else if (algorithm==my_trait<markov::MacroDMNR>::className.str())
        {
            return markov::monitorLikelihood(markov::MacroDMNR(tolerance,biNumber,Vanumber),MC,e, std::cerr);
        }
        else return Op(false,"algorithm "+algorithm+" is not recognized");
    }

    static auto get_arguments()
    {
        return std::make_tuple(
            grammar::argument(C<const Experiment&>{},my_trait<Experiment>::className.c_str()),
            grammar::argument(C< Model&>{},my_trait<Model>::className.c_str()),
            grammar::argument(C<const Parameters_values<Model>&>{},"model_parameters"),
            grammar::argument(C<std::string>{},"algorithm"),
            grammar::argument(C<double>{},"min_probability",1e-9),
            grammar::argument(C<double>{},"tolerance_error",1e-7),
            grammar::argument(C<double>{},"Binomial_threshold",5.0),
            grammar::argument(C<double>{},"Variance_threshold",1.0));
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



template<class Experiment,class Model, class ParametersDistribution>
struct Evidence{


    static constexpr auto className=my_static_string("evidence");

    static std::string run(const Experiment& e,
                    const Model& m ,
                    const ParametersDistribution& p,
                    std::string algorithm,
                    double pjump,
                    double eps_Gradient,
                    double min_P,
                    double tolerance,
                    double BiNumber,
                    double VaNumber,
                    std::mt19937_64::result_type initseed,
                    std::vector<double> betas,
                    std::vector<double>landa,
                    std::vector<std::vector<double>>landa_50_hill,
                    double gain_moment,
                    std::size_t nSamples,
                    std::size_t ntrials,
                    bool parameters,
                    bool gradient,
                    double epsf,
                    std::string idfile
                                                );

    static auto get_arguments()
    {
        return std::make_tuple(
            grammar::argument(C<const Experiment&>{},my_trait<Experiment>::className.c_str()),
            grammar::argument(C< Model&>{},my_trait<Model>::className.c_str()),
            grammar::argument(C<const ParametersDistribution&>{},"model_parameters_distribution"),
            grammar::argument(C<std::string>{},"algorithm"),
            grammar::argument(C<double>{},"p_jump"),
            grammar::argument(C<double>{},"eps_G"),
            grammar::argument(C<double>{},"min_probability"),
            grammar::argument(C<double>{},"tolerance_error"),
            grammar::argument(C<double>{},"Binomial_threshold",5.0),
            grammar::argument(C<double>{},"Variance_threshold",1.0),
            grammar::argument(C<std::mt19937_64::result_type>{},"initseed"),
            grammar::argument(C<std::vector<double>>{},"betas"),
            grammar::argument(C<std::vector<double>>{},"landas"),
            grammar::argument(C<std::vector<std::vector<double>>>{},"landa_50_hill"),
            grammar::argument(C<double>{},"gain_moment"),
            grammar::argument(C<std::size_t>{},"nSamples"),
            grammar::argument(C<std::size_t>{},"ntrials"),
            grammar::argument(C<bool>{},"parameters_output"),
            grammar::argument(C<double>{},"gradient_output"),
            grammar::argument(C<double>{},"eps_factor"),
            grammar::argument(C<std::string>{},"id_file")

                                    );
    }
};



template<class Experiment,class Model, class ParametersDistribution>
struct Evidence_prob{


    static constexpr auto className=my_static_string("evidence_prob");

    static std::string run(const Experiment& e,
                    const Model& m ,
                    const ParametersDistribution& p,
                    std::string algorithm,
                    double pjump,
                    double eps_Gradient,
                    double min_P,
                    double tolerance,
                    double BiNumber,
                    double VaNumber,
                    std::mt19937_64::result_type initseed,
                    std::vector<double> betas,
                    std::vector<double>landa,
                    double target_probability,
                    std::size_t nSamples,
                    std::size_t ntrials,
                    bool parameters,
                    bool gradient,
                    double epsf,
                    std::string idfile
                                                                        );

    static auto get_arguments()
    {
        return std::make_tuple(
            grammar::argument(C<const Experiment&>{},my_trait<Experiment>::className.c_str()),
            grammar::argument(C< Model&>{},my_trait<Model>::className.c_str()),
            grammar::argument(C<const ParametersDistribution&>{},"model_parameters_distribution"),
            grammar::argument(C<std::string>{},"algorithm"),
            grammar::argument(C<double>{},"p_jump"),
            grammar::argument(C<double>{},"eps_G"),
            grammar::argument(C<double>{},"min_probability"),
            grammar::argument(C<double>{},"tolerance_error"),
            grammar::argument(C<double>{},"Binomial_threshold",5.0),
            grammar::argument(C<double>{},"Variance_threshold",1.0),
            grammar::argument(C<std::mt19937_64::result_type>{},"initseed"),
            grammar::argument(C<std::vector<double>>{},"betas"),
            grammar::argument(C<std::vector<double>>{},"landas"),
            grammar::argument(C<double>{},"target_probability"),
            grammar::argument(C<std::size_t>{},"nSamples"),
            grammar::argument(C<std::size_t>{},"ntrials"),
            grammar::argument(C<bool>{},"parameters_output"),
            grammar::argument(C<double>{},"gradient_output"),
            grammar::argument(C<double>{},"eps_factor"),
            grammar::argument(C<std::string>{},"id_file")

                                                    );
    }
};





template<class Experiment,class Model, class ParametersDistribution>
struct Evidence_Derivative{


    static constexpr auto className=my_static_string("evidence_derivative");

    static std::string run(const Experiment& e,
                    const Model& m ,
                    const ParametersDistribution& prior,
                    std::string algorithm,
                    double pjump,
                    double min_P,
                    double tolerance,
                    double BiNumber,
                    double VaNumber,
                    std::mt19937_64::result_type initseed,
                    std::vector<double> betas,
                    std::vector<double>landa,
                    std::vector<std::vector<double>>landa_50_hill,
                    double gain_moment,
                    std::size_t nSamples,
                    std::size_t n_trials,
                    bool parameters,
                    bool gradient,
                    std::string id_file);

    static auto get_arguments()
    {
        return std::make_tuple(
            grammar::argument(C<const Experiment&>{},my_trait<Experiment>::className.c_str()),
            grammar::argument(C< Model&>{},my_trait<Model>::className.c_str()),
            grammar::argument(C<const ParametersDistribution&>{},"model_parameters_distribution"),
            grammar::argument(C<std::string>{},"algorithm"),
            grammar::argument(C<double>{},"p_jump"),
            grammar::argument(C<double>{},"min_probability"),
            grammar::argument(C<double>{},"tolerance_error"),
            grammar::argument(C<double>{},"Binomial_threshold",5.0),
            grammar::argument(C<double>{},"Variance_threshold",1.0),
            grammar::argument(C<std::mt19937_64::result_type>{},"initseed"),
            grammar::argument(C<std::vector<double>>{},"betas"),
            grammar::argument(C<std::vector<double>>{},"landas"),
            grammar::argument(C<std::vector<std::vector<double>>>{},"landa_50_hill"),
            grammar::argument(C<double>{},"gain_moment"),
            grammar::argument(C<std::size_t>{},"nSamples"),
            grammar::argument(C<std::size_t>{},"n_trials"),
            grammar::argument(C<bool>{},"parameters_output"),
            grammar::argument(C<double>{},"gradient_output"),
            grammar::argument(C<std::string>{},"id_file")

                                                  );
    }
};


template<class Experiment,class Model, class ParametersDistribution>
struct Evidence_Derivative_prob{


    static constexpr auto className=my_static_string("evidence_derivative");

    static std::string run(const Experiment& e,
                    const Model& m ,
                    const ParametersDistribution& prior,
                    std::string algorithm,
                    double pjump,
                    double min_P,
                    double tolerance,
                    double BiNumber,
                    double VaNumber,
                    std::mt19937_64::result_type initseed,
                    std::vector<double> betas,
                    std::vector<double>landa,
                    double target_probability,
                    std::size_t nSamples,
                    std::size_t n_trials,
                    bool parameters,
                    bool gradient,
                    std::string id_file);

    static auto get_arguments()
    {
        return std::make_tuple(
            grammar::argument(C<const Experiment&>{},my_trait<Experiment>::className.c_str()),
            grammar::argument(C< Model&>{},my_trait<Model>::className.c_str()),
            grammar::argument(C<const ParametersDistribution&>{},"model_parameters_distribution"),
            grammar::argument(C<std::string>{},"algorithm"),
            grammar::argument(C<double>{},"p_jump"),
            grammar::argument(C<double>{},"min_probability"),
            grammar::argument(C<double>{},"tolerance_error"),
            grammar::argument(C<double>{},"Binomial_threshold",5.0),
            grammar::argument(C<double>{},"Variance_threshold",1.0),
            grammar::argument(C<std::mt19937_64::result_type>{},"initseed"),
            grammar::argument(C<std::vector<double>>{},"betas"),
            grammar::argument(C<std::vector<double>>{},"landas"),
            grammar::argument(C<double>{},"target_probability"),
            grammar::argument(C<std::size_t>{},"nSamples"),
            grammar::argument(C<std::size_t>{},"n_trials"),
            grammar::argument(C<bool>{},"parameters_output"),
            grammar::argument(C<double>{},"gradient_output"),
            grammar::argument(C<std::string>{},"id_file")

                                                                  );
    }
};



template<class Experiment,class Model, class ParametersDistribution>
struct Evidence_emcee{


    static constexpr auto className=my_static_string("evidence_emcee");

    static std::string run(const Experiment& e,
                    const Model& m ,
                    const ParametersDistribution& p,
                    std::string algorithm,
                    double pjump,
                    double min_P,
                    double tolerance,
                    double BiNumber,
                    double VaNumber,
                    std::mt19937_64::result_type initseed,
                    std::vector<double> betas,
                    std::vector<double>alfas,
                    std::size_t nSamples,
                    bool parameters,
                    bool gradient,
                    std::size_t numWalkers,
                    double target_prob,
                    std::size_t ntrials,
                    std::string n_file);

    static auto get_arguments()
    {
        return std::make_tuple(
            grammar::argument(C<const Experiment&>{},my_trait<Experiment>::className.c_str()),
            grammar::argument(C< Model&>{},my_trait<Model>::className.c_str()),
            grammar::argument(C<const ParametersDistribution&>{},"model_parameters_distribution"),
            grammar::argument(C<std::string>{},"algorithm"),
            grammar::argument(C<double>{},"p_jump"),
            grammar::argument(C<double>{},"min_probability"),
            grammar::argument(C<double>{},"tolerance_error"),
            grammar::argument(C<double>{},"Binomial_threshold",5.0),
            grammar::argument(C<double>{},"Variance_threshold",1.0),
            grammar::argument(C<std::mt19937_64::result_type>{},"initseed"),
            grammar::argument(C<std::vector<double>>{},"betas"),
            grammar::argument(C<std::vector<double>>{},"alfas"),
            grammar::argument(C<std::size_t>{},"nSamples"),
            grammar::argument(C<bool>{},"parameters_output"),
            grammar::argument(C<bool>{},"gradient_output"),
            grammar::argument(C<std::size_t>{},"numWalkers"),
            grammar::argument(C<double>{},"target_prob"),
            grammar::argument(C<std::size_t>{},"n_trials_at_init"),
            grammar::argument(C<std::string>{},"id_file")
                                                  );
        //evidence_works = evidence_emcee ( singleLigandExperiment = mySimulation  State_Model = Model_1   model_parameters_distribution = paramPrior_1   algorithm = "MacroDMR"  Binomial_threshold =5.0 Variance_threshold =1.0  p_jump = 0.5   min_probability = 1e-14 tolerance_error=1e-2 initseed = 3034446629   betas = { 1.0 0.5 0.3 0.1 1e-2 1e-3 1e-4 0}  alfas = {2 1.5 1.2 1.1 1.05 1.02 1.01 1.005  1.002  1.001 1.0005  1.0002 1.0001 }  nSamples = 10000  parameters_output = 0  gradient_output = 0  numWalkers = 8  target_prob = 0.2  n_trials_at_init = 100 )


    }
};





struct Objects
{

    typedef Cs<State_Model,Allosteric_Model,singleLigandExperiment> types;
    typedef Cs<
        function_log10,
        simulate<singleLigandExperiment,Allosteric_Model>,
        likelihood<singleLigandExperiment,Allosteric_Model>,
      //  likelihood_detail<singleLigandExperiment,Allosteric_Model>,
        Evidence<singleLigandExperiment,Allosteric_Model,Parameters_distribution<Allosteric_Model>>,
        Evidence_prob<singleLigandExperiment,Allosteric_Model,Parameters_distribution<Allosteric_Model>>,
        Evidence_Derivative<singleLigandExperiment,Allosteric_Model,Parameters_distribution<Allosteric_Model>>,
        Evidence_Derivative_prob<singleLigandExperiment,Allosteric_Model,Parameters_distribution<Allosteric_Model>>,
        Evidence_emcee<singleLigandExperiment,Allosteric_Model,Parameters_distribution<Allosteric_Model>>,
      //  likelihoodtest<singleLigandExperiment,Allosteric_Model,Parameters_distribution<Allosteric_Model>>,
      //  likelihoodtest<singleLigandExperiment,Allosteric_Model,Parameters_partial_distribution<Allosteric_Model>>,
      //  likelihoodtest_derivative<singleLigandExperiment,Allosteric_Model,Parameters_distribution<Allosteric_Model>>,
        // not yet  likelihoodtest_derivative<singleLigandExperiment,Allosteric_Model,Parameters_partial_distribution<Allosteric_Model>>,
        simulate<singleLigandExperiment, State_Model>,
        likelihood<singleLigandExperiment, State_Model>,
        likelihood_detail<singleLigandExperiment, State_Model>,
        Evidence<singleLigandExperiment, State_Model,
                 Parameters_distribution<State_Model>>,
        Evidence_prob<singleLigandExperiment, State_Model,
                 Parameters_distribution<State_Model>>,
        Evidence_emcee<singleLigandExperiment, State_Model,Parameters_distribution<State_Model>>,
        Evidence_Derivative<singleLigandExperiment, State_Model,Parameters_distribution<State_Model>>,
        Evidence_Derivative_prob<singleLigandExperiment, State_Model,Parameters_distribution<State_Model>>,
        Evidence<singleLigandExperiment, State_Model,
                 Parameters_partial_distribution<State_Model>>,
        Evidence_prob<singleLigandExperiment, State_Model,
                 Parameters_partial_distribution<State_Model>>,
//        likelihoodtest<singleLigandExperiment, State_Model,
//                       Parameters_distribution<State_Model>>,
//        likelihoodtest<singleLigandExperiment, State_Model,
//                       Parameters_partial_distribution<State_Model>>,
//        likelihoodtest_derivative<singleLigandExperiment,State_Model,Parameters_distribution<State_Model>>,

        to_experiment,
        to_DataFrame<measure_just_y<double>>,
        to_DataFrame<markov::measure_likelihood<double>>,
        to_DataFrame<evidence::Likelihood_Analisis<Parameters_distribution<Allosteric_Model>,M_Matrix<double>,singleLigandExperiment,evidence::PartialDLogLikelihood<markov::MACROR>>>,
        to_DataFrame<evidence::Likelihood_Analisis<Parameters_distribution<State_Model>,M_Matrix<double>,singleLigandExperiment,evidence::PartialDLogLikelihood<markov::MACROR>>>
        > commands;
    typedef CCs<save,write_variable> templateCommands;
};








#endif // COMMANDS_H
