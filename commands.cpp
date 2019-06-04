#include "commands.h"

basic_Experiment<point<double, double>, measure_just_y<double>>
to_experiment::run(const io::myDataFrame<double> &da,
                   const std::string &colname_time,
                   const std::string &colname_nsample,
                   const std::string &colname_x, const std::string &colname_y,
                   double Vm, double frequency_of_sampling) {
  return experiment::DataFrame_to_Experiment(da, colname_time, colname_nsample,
                                             colname_x, colname_y, Vm,
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
    double biNumber, double Vanumber, double maxChi2, double reduce_by_p, std::size_t order_number,[[maybe_unused]]double min_connection_ratio) {
  std::cerr << "\nparameters\n" << p;

  if (algorithm == my_trait<markov::MicroDVR>::className.str()) {

        Microscopic_Model mi(m,p.at(Number_of_Channels_Parameter_label{}));

    SingleLigandModel SM(mi.Qs(p), mi.g(p), e.Vm(),
                         1.0,
                         p.at(gaussian_noise_Parameter_label()), min_P);

    Markov_Model_calculations<Markov_Transition_step_double,
                              Markov_Transition_rate, SingleLigandModel,
                              Experiment, double>
        MC(SM, e, 1, tolerance);

      auto out = markov::logLikelihood(
          markov::MicroDVR(), MC, e, std::cerr);
      if (out.has_value())
        std::cerr << "logLikelihodd = " << out.value() << std::endl;
      else
        std::cerr << out.error();
      return out;
  }

//  else if (algorithm == my_trait<markov::RMicroDVR>::className.str()) {

//      SingleLigandModel SM(m.Qs(p), m.g(p), e.Vm(),
//                           p.at(Number_of_Channels_Parameter_label()),
//                           p.at(gaussian_noise_Parameter_label()), min_P);

//      Markov_Model_calculations<Markov_Transition_step_double,
//                                Markov_Transition_rate, SingleLigandModel,
//                                Experiment, double>
//          MC(SM,m.connections_to_0(),m.connections_to_a(),m.connections_to_x(1.0),m.connections_from_x(0.0),m.connections_from_x(1.0), e, 1, tolerance);

//      auto out = markov::logLikelihood(
//          markov::RMicroDVR(maxChi2,reduce_by_p,min_connection_ratio), MC, e, std::cerr);
//      if (out.has_value())
//          std::cerr << "logLikelihodd = " << out.value() << std::endl;
//      else
//          std::cerr << out.error();
//      return out;
//  }
  else if (algorithm == my_trait<markov::RMicroDVR_new>::className.str()) {

      SingleLigandModel SM(m.Qs(p), m.g(p), e.Vm(),
                           p.at(Number_of_Channels_Parameter_label()),
                           p.at(gaussian_noise_Parameter_label()), min_P);

      Markov_Model_calculations<Markov_Transition_step_double,
                                Markov_Transition_rate, SingleLigandModel,
                                Experiment, double>
          MC(SM,m.connections_to_0(),m.connections_to_a(),m.connections_to_x(1.0),m.connections_from_x(0.0),m.connections_from_x(1.0), e, 1, tolerance);

      auto out = markov::logLikelihood(
          markov::RMicroDVR_new(min_P,reduce_by_p, maxChi2, order_number), MC, e, std::cerr);
      if (out.has_value())
          std::cerr << "logLikelihodd = " << out.value() << std::endl;
      else
          std::cerr << out.error();
      return out;
  }

  else {
    SingleLigandModel SM(m.Qs(p), m.g(p), e.Vm(),
                         p.at(Number_of_Channels_Parameter_label()),
                         p.at(gaussian_noise_Parameter_label()), min_P);

    Markov_Model_calculations<Markov_Transition_step_double,
                              Markov_Transition_rate, SingleLigandModel,
                              Experiment, double>
        MC(SM, e, 1, tolerance);

     markov::Markov_Model_calculations_new<true> MC_new(e.frequency_of_sampling(),min_P);
     E_t<markov::Markov_Model_calculations_new<true>> MC_new_Et(e.frequency_of_sampling(),min_P);


    if (algorithm == my_trait<markov::MacroDVR>::className.str()) {
      auto out = markov::logLikelihood(
          markov::MacroDVR(tolerance, biNumber, Vanumber), MC, e, std::cerr);
      if (out.has_value())
        std::cerr << "logLikelihodd = " << out.value() << std::endl;
      else
        std::cerr << out.error();
      return out;
    } else if (algorithm == my_trait<markov::MacroDMR>::className.str()) {
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
    else if (algorithm == my_trait<markov::MacroSDMR>::className.str()) {
        auto out = markov::logLikelihood_new(
      //      markov::MacroSDMR(tolerance,biNumber,Vanumber), MC_new,SM, e, std::cerr);
            E_t<markov::MacroSDMR>(tolerance,biNumber,Vanumber), MC_new_Et,SM, e, std::cerr);
        if (out.has_value())
            std::cerr << "logLikelihodd = " << out.value() << std::endl;
        else
            std::cerr << out.error();
        return myOptional_t<double>(std::move(out).value().center());
    }
    else if (algorithm == my_trait<markov::MacroSDVR>::className.str()) {
        auto out = markov::logLikelihood_new(
            markov::MacroSDVR(tolerance,biNumber,Vanumber), MC_new,SM, e, std::cerr);
        if (out.has_value())
            std::cerr << "logLikelihodd = " << out.value() << std::endl;
        else
            std::cerr << out.error();
        return out;
    } else
      return myOptional_t<double>(false,
                                  " algorithm " + algorithm + "not found");
  }
}


template <class Experiment, class Model, class ParametersDistribution>
myOptional_t<Der_t<double>> likelihood_der<Experiment, Model, ParametersDistribution>::run(
    const Experiment &e, Model &m, const Parameters_values_new<Model> &p,
    const ParametersDistribution &prior,
    const std::string algorithm, double min_P, double tolerance,
    double biNumber, double Vanumber, double maxChi2 [[maybe_unused]], double reduce_by_p [[maybe_unused]], std::size_t order_number [[maybe_unused]],double min_connection_ratio [[maybe_unused]]) {
    std::cerr << "\nparameters\n" << p;

    typedef myOptional_t<Der_t<double>> Op;
    auto dm = Der_t<Model>(m);
    auto dprior = Der_t<ParametersDistribution>(prior);


//    if (algorithm == my_trait<markov::MicroDVR>::className.str()) {

//        Microscopic_Model mi(m,p.at(Number_of_Channels_Parameter_label{}));

//        SingleLigandModel SM(mi.Qs(p), mi.g(p), e.Vm(),
//                             1.0,
//                             p.at(gaussian_noise_Parameter_label()), min_P);

//        Markov_Model_calculations<Markov_Transition_step_double,
//                                  Markov_Transition_rate, SingleLigandModel,
//                                  Experiment, double>
//            MC(SM, e, 1, tolerance);

//        auto out = markov::logLikelihood(
//            markov::MicroDVR(), MC, e, std::cerr);
//        if (out.has_value())
//            std::cerr << "logLikelihodd = " << out.value() << std::endl;
//        else
//            std::cerr << out.error();
//        return out;
//    }

    //  else if (algorithm == my_trait<markov::RMicroDVR>::className.str()) {

    //      SingleLigandModel SM(m.Qs(p), m.g(p), e.Vm(),
    //                           p.at(Number_of_Channels_Parameter_label()),
    //                           p.at(gaussian_noise_Parameter_label()), min_P);

    //      Markov_Model_calculations<Markov_Transition_step_double,
    //                                Markov_Transition_rate, SingleLigandModel,
    //                                Experiment, double>
    //          MC(SM,m.connections_to_0(),m.connections_to_a(),m.connections_to_x(1.0),m.connections_from_x(0.0),m.connections_from_x(1.0), e, 1, tolerance);

    //      auto out = markov::logLikelihood(
    //          markov::RMicroDVR(maxChi2,reduce_by_p,min_connection_ratio), MC, e, std::cerr);
    //      if (out.has_value())
    //          std::cerr << "logLikelihodd = " << out.value() << std::endl;
    //      else
    //          std::cerr << out.error();
    //      return out;
    //  }
    /*else if (algorithm == my_trait<markov::RMicroDVR_new>::className.str()) {

        SingleLigandModel SM(m.Qs(p), m.g(p), e.Vm(),
                             p.at(Number_of_Channels_Parameter_label()),
                             p.at(gaussian_noise_Parameter_label()), min_P);

        Markov_Model_calculations<Markov_Transition_step_double,
                                  Markov_Transition_rate, SingleLigandModel,
                                  Experiment, double>
            MC(SM,m.connections_to_0(),m.connections_to_a(),m.connections_to_x(1.0),m.connections_from_x(0.0),m.connections_from_x(1.0), e, 1, tolerance);

        auto out = markov::logLikelihood(
            markov::RMicroDVR_new(min_P,reduce_by_p, maxChi2, order_number), MC, e, std::cerr);
        if (out.has_value())
            std::cerr << "logLikelihodd = " << out.value() << std::endl;
        else
            std::cerr << out.error();
        return out;
    }

    else*/ {
        auto x=dprior.Parameter_to_tr(p);
        auto odp=dprior.tr_to_Parameter(x);
        if (!odp.has_value())
            return Op(false,odp.error());
        else{
            auto dp=std::move(odp).value();
        Der_t<SingleLigandModel_new> SM(dm.Qs(dp), dm.g(dp), e.Vm(),
                             dp.at(Number_of_Channels_Parameter_label()),
                             dp.at(gaussian_noise_Parameter_label()), min_P);


        Der_t<markov::Markov_Model_calculations_new<true>> MC_new_Der(e.frequency_of_sampling(),min_P);


         if (algorithm == my_trait<markov::MacroSDMR>::className.str()) {
            auto out = markov::logLikelihood_new(
                //      markov::MacroSDMR(tolerance,biNumber,Vanumber), MC_new,SM, e, std::cerr);
                Der_t<markov::MacroSDMR>(tolerance,biNumber,Vanumber), MC_new_Der,SM, e, std::cerr);
            if (out.has_value())
                std::cerr << "logLikelihodd = " << out.value() << std::endl;
            else
                std::cerr << out.error();
            return out;
        }
        else if (algorithm == my_trait<markov::MacroSDVR>::className.str()) {
            auto out = markov::logLikelihood_new(
                Der_t<markov::MacroSDVR>(tolerance,biNumber,Vanumber), MC_new_Der,SM, e, std::cerr);
            if (out.has_value())
                std::cerr << "logLikelihodd = " << out.value() << std::endl;
            else
                std::cerr << out.error();
            return out;
        } else
            return Op(false,
                                        " algorithm " + algorithm + "not found");
        }}
}



template struct simulate<singleLigandExperiment, Allosteric_Model>;
template struct likelihood<singleLigandExperiment, Allosteric_Model>;
template struct likelihood_der<singleLigandExperiment, Allosteric_Model_new,Parameters_distribution_new<Allosteric_Model_new>>;
template struct simulate<singleLigandExperiment, State_Model>;
template struct likelihood<singleLigandExperiment, State_Model>;
template struct likelihood_der<singleLigandExperiment, State_Model_new,Parameters_distribution_new<State_Model_new>>;

template struct to_DataFrame<measure_just_y<double>>;
template struct to_DataFrame<markov::measure_likelihood<double>>;
template struct to_DataFrame<evidence::Likelihood_Analisis<
    Parameters_distribution<Allosteric_Model>, M_Matrix<double>,
    singleLigandExperiment, evidence::PartialDLogLikelihood<markov::MACROR>>>;
template struct to_DataFrame<evidence::Likelihood_Analisis<
    Parameters_distribution<State_Model>, M_Matrix<double>,
    singleLigandExperiment, evidence::PartialDLogLikelihood<markov::MACROR>>>;
