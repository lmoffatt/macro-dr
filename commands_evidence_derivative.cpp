#include "commands.h"

template <class Experiment, class Model, class ParametersDistribution>
std::string Evidence_Derivative<Experiment, Model, ParametersDistribution>::run(
    const Experiment &e, const Model &m, const ParametersDistribution &prior,
    std::string algorithm, double pjump, double min_P, double tolerance,
    double BiNumber, double VaNumber,
    std::mt19937_64::result_type initseed,
    std::vector<double> betas, std::vector<double> landa,
    std::vector<std::vector<double>> landa_50_hill, double gain_moment,
    std::size_t nSamples, std::size_t n_trials, bool parameters, bool gradient,
    std::string id_file) {
    typedef Derivative<
        Markov_Model_Likelihood<Model, ParametersDistribution, Experiment>>
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

    std::mt19937_64 mt = init_mt(initseed, std::cerr);
    auto dm = Derivative<Model>(m);
    auto dprior = Derivative<ParametersDistribution>(prior);
    Derivative<Markov_Model_Likelihood<Model, ParametersDistribution, Experiment>>
        lik(e, Derivative<Markov_Model_Likelihood<Model, ParametersDistribution>>(
                   dm, dprior, algorithm, min_P, tolerance, BiNumber, VaNumber));

    std::size_t sample_interval = 10;
    auto out = evidence::OutputGenerator(Cs<RG>{}, Cs<MCMC>{}, Cs<Th_Models>{},
                                         Cs<Adaptive>{},
                                         [sample_interval](std::size_t i_sample) {
                                             return i_sample % sample_interval == 0;
                                         },
                                         std::cerr, parameters, gradient);

    std::string info;
    std::string id_f = id_file + time_now() + "_" + std::to_string(initseed);

    return evidence::run_Thermo_Levenberg_ProbVel(
               id_f, info,
               evidence::Prior_Model<Model, ParametersDistribution>(prior), lik,
               mt, betas, landa, landa_50_hill, pjump, gain_moment, nSamples,
               n_trials, out)
        .error();
}

template <class Experiment, class Model, class ParametersDistribution>
std::string
Evidence_Derivative_prob<Experiment, Model, ParametersDistribution>::run(
    const Experiment &e, const Model &m, const ParametersDistribution &prior,
    std::string algorithm, double pjump, double min_P, double tolerance,
    double BiNumber, double VaNumber,
    std::mt19937_64::result_type initseed,
    std::vector<double> betas, std::vector<double> landa,
    double target_probability, std::size_t nSamples, std::size_t n_trials,
    bool parameters, bool gradient, std::string id_file) {
    typedef Derivative<
        Markov_Model_Likelihood<Model, ParametersDistribution, Experiment>>
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

    std::mt19937_64 mt = init_mt(initseed, std::cerr);
    auto dm = Derivative<Model>(m);
    auto dprior = Derivative<ParametersDistribution>(prior);
    Derivative<Markov_Model_Likelihood<Model, ParametersDistribution, Experiment>>
        lik(e, Derivative<Markov_Model_Likelihood<Model, ParametersDistribution>>(
                   dm, dprior, algorithm, min_P, tolerance, BiNumber, VaNumber));

    std::size_t sample_interval = 10;
    auto out = evidence::OutputGenerator(Cs<RG>{}, Cs<MCMC>{}, Cs<Th_Models>{},
                                         Cs<Adaptive>{},
                                         [sample_interval](std::size_t i_sample) {
                                             return i_sample % sample_interval == 0;
                                         },
                                         std::cerr, parameters, gradient);

    std::string info;

    std::string id_f = id_file + time_now() + "_" + std::to_string(initseed);
    return evidence::run_Thermo_Levenberg_Prob(
               id_f, info,
               evidence::Prior_Model<Model, ParametersDistribution>(prior), lik,
               mt, betas, landa, target_probability, pjump, nSamples, n_trials,
               out)
        .error();
}

template struct Evidence_Derivative<singleLigandExperiment, Allosteric_Model,
                                    Parameters_distribution<Allosteric_Model>>;
template struct Evidence_Derivative_prob<
    singleLigandExperiment, Allosteric_Model,
    Parameters_distribution<Allosteric_Model>>;

template struct Evidence_Derivative<singleLigandExperiment, State_Model,
                                    Parameters_distribution<State_Model>>;
template struct Evidence_Derivative_prob<singleLigandExperiment, State_Model,
                                         Parameters_distribution<State_Model>>;
