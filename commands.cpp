#include "commands.h"
#include "likelihood_markov_process_derivative.h"

basic_Experiment<point<double, double>, measure_just_y<double>>
to_experiment::run(const io::myDataFrame<double> &da,
                   const std::string &colname_time,
                   const std::string &colname_nsample,
                   const std::string &colname_x, const std::string &colname_y,
                   double Vm, double frequency_of_sampling) {
    return experiment::DataFrame_to_Experiment(
        da, colname_time, colname_nsample, colname_x, colname_y, Vm,
        frequency_of_sampling);
}

template <class measure>
io::myDataFrame<double, std::size_t, std::string> to_DataFrame<measure>::run(
    const basic_Experiment<point<double, double>, measure> &e) {
    return experiment::Experiment_steps_to_DataFrame(e);
}

template <class Parameters_Distribution, class Parameters_Values,
          class ExperimentData, class logLikelihood>
io::myDataFrame<double, std::size_t, std::string> to_DataFrame<
    evidence::Likelihood_Analisis<Parameters_Distribution, Parameters_Values,
                                  ExperimentData, logLikelihood>>::
    run(const evidence::Likelihood_Analisis<Parameters_Distribution,
                                            Parameters_Values, ExperimentData,
                                            logLikelihood> &l) {
    return l.make_Data_Frame();
}

template <class Experiment, class Model>
myOptional_t<Experiment>
simulate<Experiment, Model>::run(typename std::mt19937_64::result_type initseed,
                                 const Experiment &e, const Model &m,
                                 const Parameters_values<Model> &p,
                                 std::size_t n_sub_intervals, double max_dt,
                                 double min_P, double tolerance) {

    SingleLigandModel SM(m.Qs(p), m.g(p), e.Vm(),
                         p.at(Number_of_Channels_Parameter_label()),
                         p.at(gaussian_noise_Parameter_label()), min_P);

    Markov_Model_calculations<Markov_Transition_step, Markov_Transition_rate,
                              SingleLigandModel, Experiment, double>
        MC(SM, e, n_sub_intervals, tolerance);

    if (initseed == 0) {
        std::random_device rd;
        initseed = rd();
    }
    std::mt19937_64 mt(initseed);
    return markov::measure_experiment(get_current{}, mt, MC, e, n_sub_intervals,
                                      max_dt);
}

template <class Experiment, class Model>
myOptional_t<double> likelihood<Experiment, Model>::run(
    const Experiment &e, Model &m, const Parameters_values<Model> &p,
    const std::string algorithm, double min_P, double tolerance,
    double biNumber, double Vanumber) {
    std::cerr << "\nparameters\n" << p;
    SingleLigandModel SM(m.Qs(p), m.g(p), e.Vm(),
                         p.at(Number_of_Channels_Parameter_label()),
                         p.at(gaussian_noise_Parameter_label()), min_P);

    Markov_Model_calculations<Markov_Transition_step_double,
                              Markov_Transition_rate, SingleLigandModel,
                              Experiment, double>
        MC(SM, e, 1, tolerance);

    if (algorithm == my_trait<markov::MacroDVR>::className.str()) {
        auto out = markov::logLikelihood(
            markov::MacroDVR(tolerance, biNumber, Vanumber), MC, e, std::cerr);
        if (out.has_value())
            std::cerr << "logLikelihodd = " << out.value() << std::endl;
        else
            std::cerr << out.error();
        return out;

    }

    else if (algorithm == my_trait<markov::MacroDMR>::className.str()) {
        auto out = markov::logLikelihood(
            markov::MacroDMR(tolerance, biNumber, Vanumber), MC, e, std::cerr);
        if (out.has_value())
            std::cerr << "logLikelihodd = " << out.value() << std::endl;
        else
            std::cerr << out.error();
        return out;

    } else if (algorithm == my_trait<markov::MacroDVNR>::className.str()) {
        auto out = markov::logLikelihood(
            markov::MacroDVNR(tolerance, biNumber, Vanumber), MC, e, std::cerr);
        if (out.has_value())
            std::cerr << "logLikelihodd = " << out.value() << std::endl;
        else
            std::cerr << out.error();
        return out;

    } else if (algorithm == my_trait<markov::MacroDMNR>::className.str()) {
        auto out = markov::logLikelihood(
            markov::MacroDMNR(tolerance, biNumber, Vanumber), MC, e, std::cerr);
        if (out.has_value())
            std::cerr << "logLikelihodd = " << out.value() << std::endl;
        else
            std::cerr << out.error();
        return out;

    }

    else
        return myOptional_t<double>(false,
                                    " algorithm " + algorithm + "not found");
}

template struct simulate<singleLigandExperiment, Allosteric_Model>;
template struct likelihood<singleLigandExperiment, Allosteric_Model>;
template struct simulate<singleLigandExperiment, State_Model>;
template struct likelihood<singleLigandExperiment, State_Model>;

template struct to_DataFrame<measure_just_y<double>>;
template struct to_DataFrame<markov::measure_likelihood<double>>;
template struct to_DataFrame<evidence::Likelihood_Analisis<
    Parameters_distribution<Allosteric_Model>, M_Matrix<double>,
    singleLigandExperiment, evidence::PartialDLogLikelihood<markov::MACROR>>>;
template struct to_DataFrame<evidence::Likelihood_Analisis<
    Parameters_distribution<State_Model>, M_Matrix<double>,
    singleLigandExperiment, evidence::PartialDLogLikelihood<markov::MACROR>>>;

template struct to_DataFrame<typename likelihoodtest<
    singleLigandExperiment, State_Model,
    Parameters_distribution<State_Model>>::return_type>;
