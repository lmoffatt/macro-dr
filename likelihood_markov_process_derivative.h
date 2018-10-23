#ifndef LIKELIHOOD_MARKOV_PROCESS_DERIVATIVE_H
#define LIKELIHOOD_MARKOV_PROCESS_DERIVATIVE_H

#include "likelihood_markov_process.h"
#include "matrixderivative.h"
#include "qmodel_derivative.h"
#include "mydistributions_derivative.h"

template<>
class Derivative<markov::mp_state_information>
{
    Derivative<M_Matrix<double>> P_mean_;
    Derivative<M_Matrix<double>> P_cov_;

    Derivative<double> y_mean_;
    Derivative<double> y_var_;
    Derivative<double> plogL_;
    Derivative<double> eplogL_;


public:


    Derivative(Derivative<M_Matrix<double>>&& P_mean__,
                         Derivative<M_Matrix<double>>&& P_cov__,
                         Derivative<double>&& y_mean__,
                         Derivative<double>&& y_var__,
                         Derivative<double>&& plogL__,
                         Derivative<double>&& eplogL__
                         ):
        P_mean_{std::move(P_mean__)},
        P_cov_{std::move( P_cov__)},
        y_mean_{std::move(y_mean__)},
        y_var_{std::move(y_var__)},
        plogL_{std::move(plogL__)},eplogL_{std::move(eplogL__)}{
    }

    Derivative()=default;
    auto& P_mean()const{return P_mean_;}
    auto& P_cov()const {return P_cov_;}

    auto& y_mean()const {return y_mean_;}
    auto& y_var()const {return y_var_;}
    auto& plogL()const {return plogL_;}
    auto& eplogL()const{return eplogL_;}

    typedef  Derivative self_type;
    constexpr static auto  className=my_static_string("mp_state_information_Derivative");
    static auto get_constructor_fields()
    {
        return std::make_tuple(
                    grammar::field(C<self_type>{},"P_mean",&self_type::P_mean),
                    grammar::field(C<self_type>{},"P_cov",&self_type::P_cov),
                    grammar::field(C<self_type>{},"y_mean",&self_type::y_mean),
                    grammar::field(C<self_type>{},"y_var",&self_type::y_var),
                    grammar::field(C<self_type>{},"plogL",&self_type::plogL),
                    grammar::field(C<self_type>{},"eplogL",&self_type::eplogL)
                    );
    }

    static Derivative<markov::mp_state_information> adjust(Derivative<M_Matrix<double>>&& P_mean__,
                                       Derivative<M_Matrix<double>>&& P_cov__,
                                       Derivative<double> y_mean__,
                                       Derivative<double> y_var__,
                                       Derivative<double> plogL__,
                                       Derivative<double> eplogL__,
                                       double min_p,
                                       double min_var)
    {
        return Derivative<markov::mp_state_information>( Derivative<Probability_distribution>::normalize(std::move(P_mean__),min_p),
                                     Derivative<Probability_distribution_covariance>::normalize(std::move( P_cov__),min_p),
                                     std::move(y_mean__),
                                     Derivative<variance_value>::adjust(std::move(y_var__),min_var),
                                     std::move(plogL__),std::move(eplogL__));
    }


};

using namespace markov;

template<>
class Derivative<markov::MacroDMR>:public markov::MacroDMR
{
public:
    inline constexpr static auto const className=my_static_string("MacroDMR");

    typedef markov::MacroDMR base_type;
    template< class Model, class Step>
    myOptional_t<Derivative<mp_state_information>>
    run(const Derivative<mp_state_information>& prior,  Model& m,const Step& p,std::ostream& os)const
    {
        MACROR alg;
        return run(prior,m,p,os,alg);
    }

    template< class Model, class Step, class MACROR>
    myOptional_t<Derivative<mp_state_information>> run(const Derivative<mp_state_information>& prior,  Model& m,const Step& p,std::ostream& os , MACROR& alg)const
    {
        typedef   myOptional_t<Derivative<mp_state_information>> Op;
        auto Q_dto=m.get_P(p,0);
        if (!Q_dto)
            return Op(false,"fails in auto Q_dt=m.get_P(p,0) :"+Q_dto.error());
        Derivative<Markov_Transition_step_double> Q_dt=std::move(Q_dto).value();
        return run(prior,Q_dt,m,p,os, alg);
    }

    template< class Model, class Step>
    myOptional_t<mp_state_information> run(const Derivative<mp_state_information>& prior,  const Derivative<Markov_Transition_step_double> Q_dt,Model& m,const Step& p,std::ostream& os ,const MACROR& alg)const
    {
        typedef   myOptional_t<mp_state_information> Op;
        switch (alg) {
        case MACRO_DMNR: return MacroDMNR(tolerance()).run(prior,Q_dt,m,p,os);
        case MACRO_DVNR: return  Op(false,"Not a contemplated algorithm :" +MACROR_string[MACRO_DVNR]+" is not valid for "+className.str());
        case MACRO_DMR: return MacroDMR(tolerance(),Binomial_magic_number()).run(prior,Q_dt,m,p,os);
        case MACRO_DVR: default:return  Op(false,"Not a contemplated algorithm :" +MACROR_string[MACRO_DVR]+" is not valid for "+className.str());
        }

    }


    template< class Model, class Step>
    myOptional_t<Derivative<mp_state_information>>
    run(const Derivative<mp_state_information>& prior,
        const Derivative<Markov_Transition_step_double>& Q_dt,
        Model& m,const Step& p,std::ostream& os , MACROR& alg)const
    {
        auto y=p.y();
        if (std::isnan(y))
            alg= MACRO_DMNR;
        else
        {
            double mg=(prior.P_mean().f()*Q_dt.gmean_i().f()).getvalue();
            double g_max=max(Q_dt.g().f());
            double g_min=min(Q_dt.g().f());
            double g_range=g_max-g_min;
            double N=m.AverageNumberOfChannels().f();
            auto p_bi=(g_max-mg)/g_range;
            auto q_bi=(mg-g_min)/g_range;
            auto test_Binomial=mp_state_information::is_Binomial_Approximation_valid(N,p_bi,q_bi,Binomial_magic_number());
            if (!test_Binomial.has_value())
                alg=MACRO_DMNR;
            else  alg=MACRO_DMR;
        }
        const MACROR alg2=alg;
        return run(prior,Q_dt,m,p,os,alg2);
    }






    template< class Model, class Step>
    myOptional_t<Derivative<mp_state_information>>
    run(const Derivative<mp_state_information>& prior,
        const Derivative<Markov_Transition_step_double>& Q_dt,
        Model& m,const Step& p,std::ostream& /*os*/)const
    {
        typedef   myOptional_t<Derivative<mp_state_information>> Op;
        auto y=p.y();

       // double mg=(prior.P_mean()*Q_dt.gmean_i()).getvalue();
        //double g_max=max(Q_dt.g());
       // double g_min=min(Q_dt.g());
   //     double g_range=g_max-g_min;
        Derivative<double> N=m.AverageNumberOfChannels();
        Derivative<double> e=m.noise_variance(p.nsamples());
        M_Matrix<double> u(prior.P_mean().f().size(),1,1.0);
        auto SmD=prior.P_cov()-diag(prior.P_mean());

        Derivative<double> gSg=(TranspMult(Q_dt.gmean_i(),SmD)*Q_dt.gmean_i()).getvalue()+(prior.P_mean()*(elemMult(Q_dt.gtotal_ij(),Q_dt.gmean_ij())*u)).getvalue();

        //   auto mu_n=;
        auto gS=TranspMult(Q_dt.gmean_i(),SmD)*Q_dt.P()+prior.P_mean()*Q_dt.gtotal_ij();

        Derivative<double> ms=(prior.P_mean()*Q_dt.gvar_i()).getvalue();

        auto e_mu=e+N*ms;
        auto y_mean=N*(prior.P_mean()*Q_dt.gmean_i()).getvalue();
        auto y_var=e_mu+N*gSg;
        auto dy=y-y_mean;
        auto chi=dy/y_var;
        auto P_mean=prior.P_mean()*Q_dt.P()+chi*gS;

        auto P__cov=quadraticForm_BT_A_B(SmD,Q_dt.P())+diag(prior.P_mean()*Q_dt.P())
                -(N/y_var)*quadraticForm_XTX(gS);

        auto chi2=dy*chi;

        double plogL;
        if (y_var.f()>0)
            plogL=-0.5*log(2.0*PI*y_var)-0.5*chi2;
        else
            plogL=std::numeric_limits<double>::infinity();

        Derivative<double> eplogL=-0.5*log(2.0*PI*y_var)+Constant(-0.5);  //e_mu+N*gSg"-N*zeta*sqr(sSg)"
        //double chilogL=(eplogL-plogL)/std::sqrt(0.5);

        auto test=mp_state_information::test(P_mean,P__cov,y_mean,y_var,plogL,eplogL,e,tolerance());
        if (!test)
        {
            std::stringstream ss;

            ss<<"\nP_mean \n"<<P_mean;
            ss<<"\nPcov \n"<<P__cov;
            //ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;

            return Op(false,"\nfails in trace!!!; error="+test.error()+ss.str());
        }
        else
            return Op(Derivative<mp_state_information>::adjust(std::move(P_mean),std::move(P__cov),std::move(y_mean),std::move(y_var),std::move(plogL),std::move(eplogL), Q_dt.min_P(),e));


    }




    template< class Model, class Step>
    myOptional_t<mp_state_information> start(Model& m,const Step& p, double min_p )const
    {
        typedef myOptional_t<mp_state_information> Op;

        auto P_mean=m.Peq(p.begin()->x());
        //  std::cerr<<"gslreijgsorjgps INIT!!!"<<P_mean.value();
        if (! P_mean)
            return Op(false,"fails to get Peq :"+P_mean.error());
        else
        {
            auto P_cov=diag(P_mean.value())-quadraticForm_XTX(P_mean.value());
            //     std::cerr<<"gslreijgsorjgps INIT!!!"<<P_cov;
            double nan=std::numeric_limits<double>::quiet_NaN();
            return Op(mp_state_information::adjust(std::move(P_mean).value(),std::move(P_cov),nan,nan,nan,nan, min_p,0));
        }
    }
    Derivative(double tolerance, double binomial_magical):base_type(tolerance,binomial_magical){}
    Derivative()=default;

    template< class Model, class Step>
    myOptional_t<mp_state_information> start(Derivative<Model>& m,const Step& p, double min_p )const
    {
        typedef myOptional_t<mp_state_information> Op;

        auto P_mean=m.Peq(p.begin()->x());
        //  std::cerr<<"gslreijgsorjgps INIT!!!"<<P_mean.value();
        if (! P_mean)
            return Op(false,"fails to get Peq :"+P_mean.error());
        else
        {
            auto P_cov=diag(P_mean.value())-quadraticForm_XTX(P_mean.value());
            //     std::cerr<<"gslreijgsorjgps INIT!!!"<<P_cov;
            double nan=std::numeric_limits<double>::quiet_NaN();
            return Op(mp_state_information::adjust(std::move(P_mean).value(),std::move(P_cov),nan,nan,nan,nan, min_p,0));
        }
    }




};





#endif // LIKELIHOOD_MARKOV_PROCESS_DERIVATIVE_H
