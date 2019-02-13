#include "commands.h"

template <class Experiment, class Model, class ParametersDistribution>
std::string Evidence_emcee<Experiment, Model, ParametersDistribution>::run(
    const Experiment &e, const Model &m, const ParametersDistribution &p,
    std::string algorithm, double pjump, double min_P, double tolerance,
    double BiNumber, double VaNumber,
    std::mt19937_64::result_type initseed,
    std::vector<double> betas, std::vector<double> alfas, std::size_t nSamples,
    bool parameters, bool gradient, std::size_t numWalkers, double target_prob,
    std::size_t ntrials, std::string n_file) {
    typedef Markov_Model_Likelihood<Model, ParametersDistribution, Experiment>
        Likelihood_Model;
    typedef evidence::Ensemble_Parallel_Tempering MCMC;

    typedef evidence::Prior_Model<Model, ParametersDistribution> PriorModel;
    typedef evidence::Thermodynamic_Model<PriorModel, Likelihood_Model> ThModel;

    typedef std::vector<std::mt19937> RG;
    typedef evidence::Thermodynamic_Model_Series<PriorModel, Likelihood_Model>
        Th_Models;
    typedef std::vector<typename evidence::Ensemble_Metropolis_Hastings::
                            adaptive_stretch_mover<ThModel>>
        Adaptive;

    Likelihood_Model lik(m, p, e, algorithm, min_P, tolerance, BiNumber,
                         VaNumber);
    std::mt19937_64 mt = init_mt(initseed, std::cerr);
    std::cerr << "\np.tr_to_Parameter(p.sample(mt))\n"
              << p.tr_to_Parameter(p.sample(mt)).value();
    std::size_t sample_interval = 10;
    auto out = evidence::OutputGenerator(Cs<RG>{}, Cs<MCMC>{}, Cs<Th_Models>{},
                                         Cs<Adaptive>{},
                                         [sample_interval](std::size_t i_sample) {
                                             return i_sample % sample_interval == 0;
                                         },
                                         std::cerr, parameters, gradient);

    std::string info = "";
    return evidence::run_Thermo_emcee(n_file, info, PriorModel(p), lik, mt, betas,
                                      numWalkers, alfas, pjump, target_prob,
                                      nSamples, ntrials, out)
        .error();
}

template struct Evidence_emcee<singleLigandExperiment, State_Model,
                               Parameters_distribution<State_Model>>;
template struct Evidence_emcee<singleLigandExperiment, Allosteric_Model,
                               Parameters_distribution<Allosteric_Model>>;