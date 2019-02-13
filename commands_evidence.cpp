#include "commands.h"



template <class Experiment, class Model, class ParametersDistribution>
std::string
Evidence<Experiment, Model, ParametersDistribution>::
    run(
        const Experiment &e, const Model &m, const ParametersDistribution &p,
        std::string algorithm, double pjump, double eps_Gradient, double min_P,
        double tolerance, double BiNumber, double VaNumber,
        std::mt19937_64::result_type initseed,
        std::vector<double> betas, std::vector<double> landa,
        std::vector<std::vector<double>> landa_50_hill, double gain_moment,
        std::size_t nSamples, std::size_t ntrials, bool parameters, bool gradient,
        double epsf, std::string idfile) {

    typedef Markov_Model_DLikelihood<Model, Experiment, ParametersDistribution>
        Likelihood_Model;
    typedef evidence::Prior_Model<Model, ParametersDistribution> PriorModel;
    typedef evidence::Thermodynamic_Model<PriorModel, Likelihood_Model> ThModel;

    typedef std::vector<std::mt19937> RG;
    typedef std::vector<evidence::Adaptive_Parameterized_Distribution_Generator<
        evidence::LevenbergMarquardt<ThModel>>>
        Adaptive;
    typedef evidence::Thermodynamic_Model_Series<PriorModel, Likelihood_Model>
        Th_Models;
    typedef evidence::Parallel_Tempering<true> MCMC;

    Likelihood_Model lik(m, p, e, algorithm, eps_Gradient, min_P, tolerance,
                         BiNumber, VaNumber, epsf);
    std::mt19937_64 mt = init_mt(initseed, std::cerr);
    std::string info = "";

    std::size_t samples_interval = 10;
    auto out = evidence::OutputGenerator(Cs<RG>{}, Cs<MCMC>{}, Cs<Th_Models>{},
                                         Cs<Adaptive>{},
                                         [samples_interval](std::size_t isample) {
                                             return isample % samples_interval == 0;
                                         },
                                         std::cerr, parameters, gradient);
    std::string id_f = idfile + time_now() + "_" + std::to_string(initseed);

    return evidence::run_Thermo_Levenberg_ProbVel(
               id_f, info,
               evidence::Prior_Model<Model, ParametersDistribution>(p), lik, mt,
               betas, landa, landa_50_hill, pjump, gain_moment, nSamples, ntrials,
               out)
        .error();
}

template <class Experiment, class Model, class ParametersDistribution>
std::string

Evidence_prob<Experiment, Model, ParametersDistribution>::run(
    const Experiment &e, const Model &m, const ParametersDistribution &p,
    std::string algorithm, double pjump, double eps_Gradient, double min_P,
    double tolerance, double BiNumber, double VaNumber,
    std::mt19937_64::result_type initseed,
    std::vector<double> betas, std::vector<double> landa,
    double target_probability, std::size_t nSamples, std::size_t ntrials,
    bool parameters, bool gradient, double epsf, std::string idfile) {

    typedef Markov_Model_DLikelihood<Model, Experiment, ParametersDistribution>
        Likelihood_Model;
    typedef evidence::Prior_Model<Model, ParametersDistribution> PriorModel;
    typedef evidence::Thermodynamic_Model<PriorModel, Likelihood_Model> ThModel;

    typedef std::vector<std::mt19937> RG;
    typedef std::vector<evidence::Adaptive_Probability_Distribution_Generator<
        evidence::LevenbergMarquardt<ThModel>>>
        Adaptive;
    typedef evidence::Thermodynamic_Model_Series<PriorModel, Likelihood_Model>
        Th_Models;
    typedef evidence::Parallel_Tempering<true> MCMC;

    Likelihood_Model lik(m, p, e, algorithm, eps_Gradient, min_P, tolerance,
                         BiNumber, VaNumber, epsf);
    std::mt19937_64 mt = init_mt(initseed, std::cerr);
    std::string info = "";
    std::size_t sample_interval = 10;
    auto out = evidence::OutputGenerator(Cs<RG>{}, Cs<MCMC>{}, Cs<Th_Models>{},
                                         Cs<Adaptive>{},
                                         [sample_interval](std::size_t i_sample) {
                                             return i_sample % sample_interval == 0;
                                         },
                                         std::cerr, parameters, gradient);
    std::string id_f = idfile + time_now() + "_" + std::to_string(initseed);

    return evidence::run_Thermo_Levenberg_Prob(
               id_f, info,
               evidence::Prior_Model<Model, ParametersDistribution>(p), lik, mt,
               betas, landa, target_probability, pjump, nSamples, ntrials, out)
        .error();
}


template struct Evidence<singleLigandExperiment, Allosteric_Model,
                         Parameters_distribution<Allosteric_Model>>;
template struct Evidence_prob<singleLigandExperiment, Allosteric_Model,
                              Parameters_distribution<Allosteric_Model>>;

template struct Evidence<singleLigandExperiment, State_Model,
                         Parameters_distribution<State_Model>>;
template struct Evidence_prob<singleLigandExperiment, State_Model,
                              Parameters_distribution<State_Model>>;
template struct Evidence<singleLigandExperiment, State_Model,
                         Parameters_partial_distribution<State_Model>>;
template struct Evidence_prob<singleLigandExperiment, State_Model,
                              Parameters_partial_distribution<State_Model>>;



