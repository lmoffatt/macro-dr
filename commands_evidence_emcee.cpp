#include "commands.h"

template <class Experiment, class Model, class ParametersDistribution>
std::string Evidence_emcee<Experiment, Model, ParametersDistribution>::run(
    const Experiment &e,  Model &m, const ParametersDistribution &p,
    std::string algorithm, double pjump, double min_P, double tolerance,
    double BiNumber, double VaNumber,
    std::mt19937_64::result_type initseed,
    std::vector<double> betas, std::vector<double> alfas, std::size_t nSamples,
     std::size_t numWalkers, double target_prob,
    std::size_t ntrials, std::string n_file ,  const  std::map<std::size_t, std::size_t>& state_sampling_cycles,
    const std::map<std::size_t, std::size_t>& gen_sampling_cycles,
    const std::map<std::size_t, std::size_t>& ana_sampling_cycles)
 {
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
    LinearIndexSampling state_sampling(state_sampling_cycles);
    LinearIndexSampling gen_sampling(gen_sampling_cycles);
    LinearIndexSampling ana_sampling(ana_sampling_cycles);
    auto out = evidence::OutputGenerator(Cs<RG>{}, Cs<MCMC>{}, Cs<Th_Models>{},
                                         Cs<Adaptive>{},
                                         state_sampling,gen_sampling,ana_sampling,
                                         std::cerr);

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
