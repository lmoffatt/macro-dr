#ifndef LIKELIHOOD_MARKOV_PROCESS_H
#define LIKELIHOOD_MARKOV_PROCESS_H
#include "Matrix.h"
#include "myDistributions.h"
#include "myoptimization.h"
#include "Experiment.h"
#include "qmodel.h"
#include "mymath.h"


#include <type_traits>

namespace markov {


class mp_state_information
{
    M_Matrix<double> P_mean_;
    M_Matrix<double> P_cov_;
    
    double y_mean_;
    double y_var_;
    double plogL_;
    double eplogL_;

    
public:
    static Op_void close_to_zero_test(const M_Matrix<double>& P_mean__,
                        double tolerance)
    {
        std::stringstream ss;
        auto ck_mean=are_in_range<true,M_Matrix<double>>(false,0,1,tolerance).test(P_mean__,ss);
        if (ck_mean)
            return Op_void(true,"");
        else
        {
            return Op_void(false,"a probability value crossed zero: "+ss.str());
        }
    }
    static Op_void test(const M_Matrix<double>& P_mean__,
                        const M_Matrix<double>& P_cov__,
                       double tolerance)
    {
        auto ck_mean=Probability_distribution::test<true>(P_mean__,tolerance);
        auto ck_cov=Probability_distribution_covariance::test<true>(P_cov__,tolerance);
        if (ck_mean.has_value()&&ck_cov.has_value())
            return Op_void(true,"");
        else
        {
            std::stringstream ss;
            ss<<" Pmean test: "<<ck_mean.error()<<" Pcov test: "<<ck_cov.error();
            return Op_void(false,ss.str());
        }
    }

    static Op_void test(const M_Matrix<double>& P_mean__,
                        const M_Matrix<double>& P_cov__,
                        double y_mean__,
                        double y_var__,
                        double plogL__,
                        double eplogL__,
                        double minvariance,
                        double tolerance)
    {
        auto ck_mean=Probability_distribution::test<true>(P_mean__,tolerance);
        auto ck_cov=Probability_distribution_covariance::test<true>(P_cov__,tolerance);
        auto ck_ymean=variable_value::test<true>(y_mean__);
        auto ck_yvar=variance_value::test<true>(y_var__,minvariance,tolerance);
        auto ck_plogL=logLikelihood_value::test<true>(plogL__);
        auto ck_eplogL=logLikelihood_value::test<true>(eplogL__);
        if (ck_mean.has_value()&&ck_cov.has_value()&&ck_ymean.has_value()&&ck_yvar.has_value()&&ck_plogL.has_value()&&ck_eplogL.has_value())
            return Op_void(true,"");
        else
        {
            std::stringstream ss;
            ss<<" Pmean test: "<<ck_mean.error()<<" Pcov test: "<<ck_cov.error()<< " yvar test:"<<ck_yvar.error();
            ss<<" pLogL test: "<<ck_plogL.error()<<" eplogL test"<<ck_eplogL.error();
            return Op_void(false,ss.str());
        }
    }


    static mp_state_information adjust(M_Matrix<double>&& P_mean__,
                                       M_Matrix<double>&& P_cov__,
                                       double y_mean__,
                                       double y_var__,
                                       double plogL__,
                                       double eplogL__,
                                       double min_p,
                                       double min_var)
    {
        return mp_state_information( Probability_distribution::normalize(std::move(P_mean__),min_p),
                                     Probability_distribution_covariance::normalize(std::move( P_cov__),min_p),
                                     y_mean__,
                                     variance_value::adjust(y_var__,min_var),
                                     plogL__,eplogL__);
    }
    mp_state_information(M_Matrix<double>&& P_mean__,
                         M_Matrix<double>&& P_cov__,
                         double y_mean__,
                         double y_var__,
                         double plogL__,
                         double eplogL__
                         ):
        P_mean_{std::move(P_mean__)},
        P_cov_{std::move( P_cov__)},
        y_mean_{y_mean__},
        y_var_{y_var__},
        plogL_{plogL__},eplogL_{eplogL__}{
    }

    mp_state_information()=default;
    const M_Matrix<double>& P_mean()const{return P_mean_;}
    const M_Matrix<double>& P_cov()const {return P_cov_;}
    
    double y_mean()const {return y_mean_;}
    double y_var()const {return y_var_;}
    double plogL()const {return plogL_;}
    double eplogL()const{return eplogL_;}
    
    typedef  mp_state_information self_type;
    constexpr static auto  className=my_static_string("mp_state_information");
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
};




class MacroDR
{
public:
    inline constexpr static auto const className=my_static_string("MacroDR");

    template< class Model, class Step>
    myOptional_t<mp_state_information> run(const mp_state_information& prior,  Model& m,const Step& p )const
    {
        typedef   myOptional_t<mp_state_information> Op;
        //double eps=std::numeric_limits<double>::epsilon();
        auto Q_dto=m.get_P(p,0);
        if (!Q_dto)
            return Op(false,"fails in auto Q_dt=m.get_P(p,0) :"+Q_dto.error());

        Markov_Transition_step_double Q_dt=std::move(Q_dto).value();
        //  Markov_Transition_step_double Q_dt2=m.get_P(p,20,eps);
        //  assert(class_Invariants<Markov_Transition_step_double>::test<true>(Q_dt, prior.tolerance()));
        //  assert(class_Invariants<Markov_Transition_step_double>::test<true>(Q_dt2, prior.tolerance()));
        //  areEqual(Q_dt,Q_dt2,1e-6);

        double e=m.noise_variance(p.nsamples());
        //      double dt=p.nsamples()/m.fs();
        double N=m.AverageNumberOfChannels();
        M_Matrix<double> u(prior.P_mean().size(),1,1.0);

        auto SmD=prior.P_cov()-diag(prior.P_mean());

        double gSg=(TranspMult(Q_dt.gmean_i(),SmD)*Q_dt.gmean_i()).getvalue()+(prior.P_mean()*(elemMult(Q_dt.gtotal_ij(),Q_dt.gmean_ij())*u)).getvalue();

        double sSg=(TranspMult(Q_dt.gvar_i(),SmD)*Q_dt.gmean_i()).getvalue()+(prior.P_mean()*(elemMult(Q_dt.gtotal_var_ij(),Q_dt.gmean_ij())*u)).getvalue();
        double sSs=(TranspMult(Q_dt.gvar_i(),SmD)*Q_dt.gvar_i()).getvalue()+(prior.P_mean()*(elemMult(Q_dt.gtotal_var_ij(),Q_dt.gvar_ij())*u)).getvalue();
        //   auto mu_n=;
        auto sS=TranspMult(Q_dt.gvar_i(),SmD)*Q_dt.P()+prior.P_mean()*Q_dt.gtotal_var_ij();
        auto gS=TranspMult(Q_dt.gmean_i(),SmD)*Q_dt.P()+prior.P_mean()*Q_dt.gtotal_ij();


        double ms=(prior.P_mean()*Q_dt.gvar_i()).getvalue();
        double ms0;
        double delta_emu=sqr(ms+e/N)-2.0/N*sSs;
        if (delta_emu>0)
        {
            ms0=(ms-e/N)/2+std::sqrt(delta_emu)/2;
           // std::cerr<<"(ms0==ms-sSs/(2*(e+N*ms0)))="<<(ms0==ms-sSs/(2*(e+N*ms0)))<< "ms0="<<ms0;
            assert((are_Equal<true,double>(std::numeric_limits<double>::epsilon()*1000).test(ms0,ms-sSs/(2*(e+N*ms0)))));

        }
        else
        {
            ms0=ms;
            //  assert(false);
        }
        if (false && (ms0<0))
        {
            std::stringstream ss;
            ss<<"\n ms0 is negative!! \n ms0=(ms-e/N)/2+std::sqrt(delta_emu)/2;\nms0="<<ms0<<"\n ms="<<ms<<"\n e/N/2="<<e/N;
            ss<<"\n N="<<N;
            ss<<"\n std::sqrt(delta_emu)/2="<<std::sqrt(delta_emu)/2;
            ss<<"\ndelta_emu=sqr(ms+e/N)-2.0/N*sSs \nsqr(ms+e/N)="<<sqr(ms+e/N)<<" 2.0/N*sSs="<<2.0/N*sSs<<" sSs="<<sSs;
            ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;
            return Op(false,ss.str());
        }


        auto e_mu=e+N*ms0;
        auto y_mean=N*(prior.P_mean()*Q_dt.gmean_i()).getvalue()-N*0.5/e_mu*sSg;
        auto zeta=N/(2*sqr(e_mu)+N*sSs);
        if (zeta<0)
        {
            std::cerr<<"\n****************************************   negative e_4!!!!****\n";
            std::cerr<<"\n t= "<<p.x().t()<<" x="<<p.x().x()<<" y="<<p.y();
            std::cerr<<"\t N*gSg= "<<N*gSg;
            std::cerr<<"\t e_4= 2*sqr(e_mu)- N*sSs="<<zeta<<"\t";
            std::cerr<<"\t N*sSs=\t"<<N*sSs;
            std::cerr<<"\t 2*sqr(e_mu)="<<2*sqr(e_mu);
            zeta=-zeta;
            std::cerr<<"\n****************************************   negative e_4!!!!****\n";

        }
        auto y_var=e_mu+N*gSg-N*zeta*sqr(sSg);
        if (!variance_value::test<true>(y_var,e,Q_dt.min_P()))
        {
            y_var=e;
        }
            auto y=p.y();


        if constexpr(true)
        {
          //  std::cerr<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;
        //    std::cerr<<"\n step="<<p<<" "<<" y="<<y<<" ymean="<<y_mean<<" yvar="<<y_var<<" e="<<e<<" e_mu"<<e_mu<<" N*gSg"<<N*gSg<<" -N*zeta*sSg="<<-N*zeta*sqr(sSg);

           // std::cerr<<"delta_emu="<<delta_emu<<"  sqr(ms+e/N)="<<sqr(ms+e/N)<<" -2/N*sSs="<<-2/N*sSs<<" ms="<<ms<<" sSs="<<sSs<<" sSg="<<sSg<<" gSg="<<gSg<<"\n";
        }

        if (std::isnan(y))
        {
            double plogL=std::numeric_limits<double>::quiet_NaN();
            double eplogL=std::numeric_limits<double>::quiet_NaN();
            auto P__cov=quadraticForm_BT_A_B(SmD,Q_dt.P());
            auto P_mean = prior.P_mean()*Q_dt.P();
            P__cov+=diag(P_mean);
            // std::cerr<<"\nPcov nana corr\n"<<P__cov<<"\nP_mean nana corr\n"<<P_mean<<"\nQ.P \n"<<Q_dt.P();
            auto test=mp_state_information::test(P_mean,P__cov,tolerance_);
            if (test.has_value())
               return Op(mp_state_information::adjust(std::move(P_mean),std::move(P__cov),y_mean,y_var,plogL,eplogL, Q_dt.min_P(),e));
            else
                return Op(false,"fails at intertrace prediction!!: "+test.error());
        }
        else
        {
            auto dy=y-y_mean;
            auto chi=dy/y_var;
            auto P_mean=prior.P_mean()*Q_dt.P()+chi*gS
                    -(chi*zeta*sSg+0.5/sqr(e_mu))*sS;

            auto P__cov=quadraticForm_BT_A_B(SmD,Q_dt.P())+diag(prior.P_mean()*Q_dt.P())
                    -(zeta+N/y_var*sqr(zeta*sSg))*quadraticForm_XTX(sS)
                    +(2.0*N/y_var*zeta*sSg)*TransposeSum(TranspMult(sS,gS))
                    -(N/y_var)*quadraticForm_XTX(gS);

            auto chi2=dy*chi;

            double plogL;
            if (y_var>0)
                plogL=-0.5*log(2*PI*y_var)-0.5*chi2;
            else
                plogL=std::numeric_limits<double>::infinity();

            double eplogL=-0.5*log(2*PI*y_var)-0.5;  //e_mu+N*gSg"-N*zeta*sqr(sSg)"
            //  std::cerr<<"\n t= "<<p.x().t()<<" nsamples="<<p.nsamples()<<" dt="<<dt<<" x="<<p.x().x()<<" y="<<p.y()<<" y_mean="<<y_mean<<" chi2="<<chi2<<" y_var"<<y_var;
            //   std::cerr<<" e="<<e<<" ms="<<ms<<" ms0="<<ms0<<" e_mu="<<e_mu<<" N*gSg="<<N*gSg;
            //   std::cerr<<" -N*zeta*sqr(sSg)="<<-N*zeta*sqr(sSg)<< " plogL="<<plogL<< " eplogL="<<eplogL<<"dif="<<plogL-eplogL<<" sSs="<<sSs;
            //  std::cerr<<"\nPcov \n"<<P__cov;
            //  std::cerr<<"\nP_mean \n"<<P_mean;
            //  std::cerr<<"\n Qdt\n"<<Q_dt;
            if (auto test=variance_value::test<true>(y_var,e,e*tolerance()); !test.has_value())
            {
                std::stringstream ss;
                ss<<"\n step="<<p<<" "<<" y="<<y<<" ymean="<<y_mean<<" yvar="<<y_var<<" e="<<e<<" e_mu"<<e_mu<<" N*gSg"<<N*gSg<<" -N*zeta*sSg="<<-N*zeta*sqr(sSg);
                ss<<"\n ms0=(ms-e/N)/2+std::sqrt(delta_emu)/2;\nms0="<<ms0<<"\n ms="<<ms<<"\n e/N/2="<<e/N;

                ss<<"\n std::sqrt(delta_emu)/2="<<std::sqrt(delta_emu)/2;
                ss<<"\ndelta_emu=sqr(ms+e/N)-2.0/N*sSs \nsqr(ms+e/N)="<<sqr(ms+e/N)<<" 2.0/N*sSs="<<2.0/N*sSs<<" sSs="<<sSs;
                ss<<"\nP_mean \n"<<P_mean;
                ss<<" -N*zeta*sqr(sSg)="<<-N*zeta*sqr(sSg)<< " plogL="<<plogL<< " eplogL="<<eplogL<<"dif="<<plogL-eplogL<<" sSs="<<sSs;
                ss<<"\nPcov \n"<<P__cov;
                ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;

                return Op(false,"error in variance!! "+test.error() +ss.str());
            }

            auto testzero=mp_state_information::close_to_zero_test(P_mean,0);
            if (!testzero.has_value())
            {
          //      std::cerr<<" \n\n----test----\n"<<test.error()<<"\n";
                P__cov=quadraticForm_BT_A_B(SmD,Q_dt.P());
                P_mean = prior.P_mean()*Q_dt.P();
                P__cov+=diag(P_mean);
                //  std::cerr<<"\n SmD \n"<<SmD;

                //  std::cerr<<"\nPcov corr\n"<<P__cov<<"\nP_mean corr\n"<<P_mean<<"\nQ.P \n"<<Q_dt.P();



            }
            auto test=mp_state_information::test(P_mean,P__cov,y_mean,y_var,plogL,eplogL,e,tolerance());
            if (!test)
            {
                std::stringstream ss;
                ss<<"\n step="<<p<<" "<<" y="<<y<<" ymean="<<y_mean<<" yvar="<<y_var<<" e="<<e<<" e_mu"<<e_mu<<" N*gSg"<<N*gSg<<" -N*zeta*sSg="<<-N*zeta*sqr(sSg);
                ss<<"\n ms0=(ms-e/N)/2+std::sqrt(delta_emu)/2;\nms0="<<ms0<<"\n ms="<<ms<<"\n e/N/2="<<e/N;

                ss<<"\n std::sqrt(delta_emu)/2="<<std::sqrt(delta_emu)/2;
                ss<<"\ndelta_emu=sqr(ms+e/N)-2.0/N*sSs \nsqr(ms+e/N)="<<sqr(ms+e/N)<<" 2.0/N*sSs="<<2.0/N*sSs<<" sSs="<<sSs;
                ss<<"\nP_mean \n"<<P_mean;
                ss<<" -N*zeta*sqr(sSg)="<<-N*zeta*sqr(sSg)<< " plogL="<<plogL<< " eplogL="<<eplogL<<"dif="<<plogL-eplogL<<" sSs="<<sSs;
                ss<<"\nPcov \n"<<P__cov;
                ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;

                return Op(false,"\nfails in trace!!!; error="+test.error()+ss.str());
            }
            else
            return Op(mp_state_information::adjust(std::move(P_mean),std::move(P__cov),y_mean,y_var,plogL,eplogL, Q_dt.min_P(),e));


        }

    }


    template< class Model, class Step>
    mp_state_information run_old(const mp_state_information& prior,  Model& m,const Step& p )const
    {
        //double dt=Y.dt();
        //  double x=Y.x();
        // double y=Y.y();
        
        double eps=std::numeric_limits<double>::epsilon();
        Markov_Transition_step_double Q_dt=m.get_P(p,0,eps);
        Markov_Transition_step_double Q_dt2=m.get_P(p,10,eps);
        Markov_Transition_step_double Q_dt3=m.get_P(p,20,eps);

        bool test, test2, test3;
        if (class_Invariants<Markov_Transition_step_double>::test<true>(Q_dt))
            test=true;
        else
            test=false;
        if (class_Invariants<Markov_Transition_step_double>::test<true>(Q_dt2))
            test2=true;
        else
            test2=false;

        if (class_Invariants<Markov_Transition_step_double>::test<true>(Q_dt3))
            test3=true;
        else
            test3=false;

        M_Matrix<double> Sm=prior.P_cov()-diag(prior.P_mean());
        
        M_Matrix<double> sSm=TranspMult(Q_dt.gvar_i(),Sm);
        M_Matrix<double> gSm=TranspMult(Q_dt.gmean_i(),Sm);
        



        double sms=0;
        double smg=0;
        double gmg=0;
        std::size_t k_u=m.nstates();
        for (std::size_t i=0; i<k_u; i++)
            for (std::size_t j=0; j<k_u; j++)
            {
                sms+=Q_dt.gvar_ij()(i,j)*prior.P_mean()[i]*Q_dt.gtotal_var_ij()(i,j);
                smg+=Q_dt.gvar_ij()(i,j)*prior.P_mean()[i]*Q_dt.gtotal_ij()(i,j);
                gmg+=Q_dt.gmean_ij()(i,j)*prior.P_mean()[i]*Q_dt.gtotal_ij()(i,j);
            };
        
        double gSg=(gSm*Q_dt.gmean_i())[0]+gmg;
        double sSg=(sSm*Q_dt.gmean_i())[0]+smg;
        double sSs=(sSm*Q_dt.gvar_i())[0]+sms;
        
        
        M_Matrix<double> gS=gSm*Q_dt.P()+prior.P_mean()*Q_dt.gtotal_ij();
        M_Matrix<double> sS=sSm*Q_dt.P()+prior.P_mean()*Q_dt.gtotal_var_ij();
        double N=m.AverageNumberOfChannels();
        double smean=N*(prior.P_mean()*Q_dt.gvar_i())[0];
        double e=m.noise_variance(p.nsamples());

        double sA2=N*(prior.P_mean()*Q_dt.gvar_i())[0]+m.noise_variance(p.nsamples());

        double sB4=2*sA2*sA2-N*sSs;



        double y_var=sA2+N*gSg+(N*N*sSg*sSg)/sB4;
        //sometimes sB4 might be negative, when P_cov_M is greter than one.
        
        if (y_var<0)
        {
            std::cerr<<" \n******** y_var_d negative!!!!!***********\n";
            std::cerr<<"\n sA2+N*gSg+(N*N*sSg*sSg)/sB4 \n";
            std::cerr<<"\t smean="<<smean<<"\t e="<<e;
            std::cerr<<"\n sA2=  "<<sA2;
            std::cerr<<"\n sB4=  "<<sB4;

            std::cerr<<"\n gSg ="<<gSg;
            std::cerr<<"\n sSs ="<<sSs;

            std::cerr<<"\n sB4=2*sA2*sA2-N*sSs  sB4 = "<<sB4<<"\n";
            
            //        std::cerr<<*this;
            //        std::cerr<<this->model();
            //        //press_any_key_to_continue();
        }
        
        //auto y_std=std::sqrt(y_var);
        
        double y_mean=N*(prior.P_mean()*Q_dt.gmean_i())[0]-sA2/sB4*sSg;
        
        // product of S and g, used several places
        auto y=p.y();
        

        
        if (std::isnan(y))
        {
            double plogL=std::numeric_limits<double>::quiet_NaN();
            double eplogL=std::numeric_limits<double>::quiet_NaN();
            //auto chi2=std::numeric_limits<double>::quiet_NaN();
            /* auto P_cov=TranspMult(Q_dt.P(),(prior.P_cov()-diag(prior.P_mean()))*Q_dt.P());*/
            auto P_cov=quadraticForm_BT_A_B((prior.P_cov()-diag(prior.P_mean())),Q_dt.P());

            auto P_mean = prior.P_mean()*Q_dt.P();
            P_cov+=diag(prior.P_mean());
            return mp_state_information(std::move(P_mean),std::move(P_cov),y_mean,y_var,plogL,eplogL);
            
        }
        else
        {
            M_Matrix<double> P_mean=prior.P_mean()*Q_dt.P();
            double dy=y-y_mean;
            P_mean+=sS*(N*sSg/sB4*(p.y()-y_mean)/y_var-sA2/sB4)
                    +gS*((p.y()-y_mean)/y_var);
            if (!(P_mean>=0.0))
            {
                double summ=0;
                for (std::size_t i=0;i<size(P_mean);i++)
                {
                    if (P_mean[i]<0)
                        P_mean[i]=0;
                    else if (P_mean[i]>1)
                        P_mean[i]=1;
                    summ+=P_mean[i];
                }
                P_mean/=summ;
            }
            auto chi2=dy*dy/y_var;
            
            double plogL;
            if (y_var>0)
                plogL=-0.5*log(2*PI*y_var)-0.5*chi2;
            else
                plogL=std::numeric_limits<double>::infinity();
            
            double eplogL=-0.5*(1.0+log(2*PI*y_var));
            
            auto P_cov_=TranspMult(Q_dt.P(),Sm)*Q_dt.P()+diag(prior.P_mean()*Q_dt.P())-
                    TranspMult(gS,gS)*(N/y_var)-
                    (TranspMult(sS,gS)+TranspMult(gS,sS))*(N*N*sSg/y_var/sB4)+
                    TranspMult(sS,sS)*(N/sB4-N/y_var*N*N*sSg*sSg/sB4/sB4);
            
            auto P_cov=quadraticForm_BT_A_B(Sm,Q_dt.P())+
                    diag(prior.P_mean()*Q_dt.P())-
                    quadraticForm_XTX(gS)*(N/y_var)-
                    TransposeSum(TranspMult(sS,gS))*(N*N*sSg/y_var/sB4)+
                    quadraticForm_XTX(sS)*(N/sB4-N/y_var*N*N*sSg*sSg/sB4/sB4);

            //assert((are_Equal<true, M_Matrix<double>>(Q_dt.min_P(),Q_dt.min_P()).test_prod(P_cov_,P_cov,std::cerr)));

           if (!(diag(P_cov)>=0.0)||!(diag(P_cov)<=1.0))
            {
                for (std::size_t i=0;i<ncols(P_cov);i++)
                {
                    if (P_cov(i,i)<0.0)
                        P_cov(i,i)=std::abs(P_cov(i,i));
                    else if (P_cov(i,i)>1.0)
                        P_cov(i,i)=1.0;
                }
            }
            return mp_state_information(std::move(P_mean),std::move(P_cov),y_mean,y_var,plogL,eplogL);
            
            
        }
        
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
    double tolerance()const {return tolerance_;}

    MacroDR(double tolerance):tolerance_{tolerance}{}
    MacroDR()=default;
private:
   double tolerance_=1e-5;
    
};


struct logLikelihood_function
{
    template<class Experiment>
    double operator()(const Experiment& )const
    {
        return 0;
    }
    template<class mp_state_information>
    void operator()(const mp_state_information& mp, double& logLsum, std::size_t& i)const
    {
        if (std::isfinite(mp.plogL()))
        {
            logLsum+=mp.plogL();
            ++i;
        }
    }

};

struct partialLogLikelihood_function
{
    template<class Experiment>
    std::vector<double> operator()(const Experiment& e)const
    {
        std::size_t n=e.num_measurements();

        return std::vector<double>(n);
    }

    template<class mp_state_information>
    void operator()(const mp_state_information& mp,std::vector<double>& v, std::size_t& i)const
    {
        if (std::isfinite(mp.plogL()))
        {
            v[i]=mp.plogL();
            ++i;
        }
    }

};

struct partialDistribution_function
{
    template<class Experiment>
    auto operator()(const Experiment& e)const
    {
        std::size_t n=e.num_measurements();

        return std::vector<Normal_Distribution<double>>(n);
    }

    // template<class mp_state_information>
    void operator()(const mp_state_information& mp,std::vector<Normal_Distribution<double>>& v, std::size_t& i)const
    {
        if (std::isfinite(mp.plogL()))
        {
            v[i]=Normal_Distribution<double>(mp.y_mean(),mp.y_var());
            ++i;
        }
        else // hack to avoid including intervals in the Jacobian calculation
        {
            v[i]=Normal_Distribution<double>(0,1);
            ++i;
        }
    }

};




template<class F, class MacroDR,class Model,template<class, class > class Experiment, class Point, class Measure>
auto logLikelihood_experiment_calculation(const F& f,const MacroDR& a,  Model& m,const Experiment<Point,Measure>& e)
{
    auto out=f(e);
    typedef myOptional_t<std::decay_t<decltype (out)>> Op;

    auto first_step=*e.begin_begin();
    auto prior=a.start(m,first_step, m.min_P());
    if (!prior.has_value())
        return Op(false,"calculation interrupted at start :"+prior.error());
    std::size_t i=0;
    for (auto it=e.begin_begin(); it!=e.end_end(); ++it)
    {
        auto post=a.run(prior.value(),m,*it);
        if (!post.has_value())
            return Op(false,"calculation interrupted at "+ToString(*it)+"  :"+post.error());
        else
        {
         f(post.value(),out,i);
          prior.value()=std::move(post).value();
        }
    }
    return Op(out);
    
}
template<class MacroDR,class Model,class Experiment>
myOptional_t<double> logLikelihood(const MacroDR& a, Model& m, const Experiment e )
{
    return logLikelihood_experiment_calculation(logLikelihood_function(),a,m,e);
}

template<class MacroDR,class Model,class Experiment>
myOptional_t<std::vector<double>> partialLogLikelihood(const MacroDR& a, Model& m, const Experiment e )
{
    return logLikelihood_experiment_calculation(partialLogLikelihood_function(),a,m,e);
}

template<class MacroDR,class Model,class Experiment>
myOptional_t<std::vector<Normal_Distribution<double>>> partialDistribution(const MacroDR& a, Model& m, const Experiment e )
{
    return logLikelihood_experiment_calculation(partialDistribution_function(),a,m,e);
}


template<class Y>
struct measure_likelihood:public experiment::measure_just_y<Y>
{
    double y_mean;
    double y_var;
    double plogL;
    double eplogL;

    measure_likelihood()=default;
    measure_likelihood(const mp_state_information& mp): y_mean{mp.y_mean()},y_var{mp.y_var()},plogL{mp.plogL()},eplogL{mp.eplogL()}{}

    measure_likelihood(const std::tuple<double,double,double,double>& data):
        y_mean{std::get<0>(data)},y_var{std::get<1>(data)},plogL{std::get<2>(data)},eplogL{std::get<3>(data)}{}

    static void insert_col(io::myDataFrame<double,std::size_t>& d)
    {
        d.insert_column("ymean",C<double>{});
        d.insert_column("y_var",C<double>{});
        d.insert_column("plogL",C<double>{});
        d.insert_column("eplogL",C<double>{});
    }
    std::tuple<double,double,double,double> data()const
    {
        return {y_mean,y_var,plogL,eplogL};
    }
};


template< class measure_likelihood>
struct partialLogLikelihood_monitor_function
{
    template<class Experiment>
    std::vector<measure_likelihood> operator()(const Experiment& e)const
    {
        std::size_t n=e.end_end()-e.begin_begin();
        return std::vector<measure_likelihood>(n);
    }

    template<class mp_state_information>
    void operator()(const mp_state_information& mp,std::vector<measure_likelihood>& v, std::size_t& i)const
    {
        v[i]=measure_likelihood(mp);
        ++i;
    }

};


template<class MacroDR,class Model,class Experiment>
auto monitorLikelihood(const MacroDR& a, Model& m, const Experiment e )
{
    typedef  myOptional_t< experiment::basic_Experiment<experiment::point<double,double>,measure_likelihood<double>>> Op;
    auto v= logLikelihood_experiment_calculation(partialLogLikelihood_monitor_function<measure_likelihood<double>>(),a,m,e);
    if (!v.has_value())
        return Op(false,"logLikelihood_experiment_calculation fails :"+v.error());
    else
     return Op(experiment::basic_Experiment<experiment::point<double,double>,measure_likelihood<double>>(e,std::move(v).value()));
}










} // namespace markov

#endif // LIKELIHOOD_MARKOV_PROCESS_H
