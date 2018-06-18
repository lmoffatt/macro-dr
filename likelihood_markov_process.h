#ifndef LIKELIHOOD_MARKOV_PROCESS_H
#define LIKELIHOOD_MARKOV_PROCESS_H
#include "Matrix.h"
#include "myDistributions.h"
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
    mp_state_information(M_Matrix<double>&& P_mean__,
                         M_Matrix<double>&& P_cov__,
                         double y_mean__,
                         double y_var__,
                         double plogL__,
                         double eplogL__):
        P_mean_{P_mean__},P_cov_{P_cov__}, y_mean_{y_mean__}, y_var_{y_var__},plogL_{plogL__},eplogL_{eplogL__}{}
    
    const M_Matrix<double>& P_mean()const{return P_mean_;}
    const M_Matrix<double>& P_cov()const {return P_cov_;}
    
    double y_mean()const {return y_mean_;}
    double y_var()const {return y_var_;}
    double plogL()const {return plogL_;}
    double eplogL()const{return eplogL_;}
    
    
};




template<bool zero_guard>
struct MacroDR
{
    inline constexpr static auto const className=my_static_string("MacroDR");
    template< class Model, class Step>
    mp_state_information run(const mp_state_information& prior,  Model& m,const Step& p )const
    {
        //double dt=Y.dt();
        //  double x=Y.x();
        // double y=Y.y();
        
        
        Markov_Transition_step_double Q_dt=m.get_P(p);
        M_Matrix<double> Sm=prior.P_cov()-diag(prior.P_mean());
        
        M_Matrix<double> sSm=TranspMult(Q_dt.gtotal_var_i(),Sm);
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
        double sSs=(sSm*Q_dt.gtotal_var_i())[0]+sms;
        
        
        M_Matrix<double> gS=gSm*Q_dt.P()+prior.P_mean()*Q_dt.gtotal_ij();
        M_Matrix<double> sS=sSm*Q_dt.P()+prior.P_mean()*Q_dt.gtotal_var_ij();
        double N=m.AverageNumberOfChannels();
        double sA2=N*(prior.P_mean()*Q_dt.gtotal_var_i())[0]+m.noise_variance(p.nsamples());
        double sB4=2*sA2*sA2-N*sSs;
        
        double y_var=sA2+N*gSg+(N*N*sSg*sSg)/std::abs(sB4);
        //sometimes sB4 might be negative, when P_cov_M is greter than one.
        
        if (y_var<0)
        {
            //        std::cerr<<" \n******** y_var_d negative!!!!!***********\n";
            //        std::cerr<<"\n sA2+N*gSg+(N*N*sSg*sSg)/sB4 \n";
            //        std::cerr<<"\n sA2=  "<<sA2;
            //        std::cerr<<"\n gSg ="<<gSg;
            //        std::cerr<<"\n sB4=2*sA2*sA2-N*sSs  sB4 = "<<sB4<<"\n";
            
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
            auto P_cov=TranspMult(Q_dt.P(),(prior.P_cov()-diag(prior.P_mean()))*Q_dt.P());
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
            if constexpr(zero_guard)
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
            
            auto P_cov=TranspMult(Q_dt.P(),Sm)*Q_dt.P()+diag(prior.P_mean()*Q_dt.P())-
                    TranspMult(gS,gS)*(N/y_var)-
                    (TranspMult(sS,gS)+TranspMult(gS,sS))*(N*N*sSg/y_var/sB4)+
                    TranspMult(sS,sS)*(N/sB4-N/y_var*N*N*sSg*sSg/sB4/sB4);
            
            if constexpr (zero_guard)
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
    mp_state_information start(Model& m,const Step& p )const {
        auto P_mean=m.Peq(p.begin()->x());
        auto P_cov=diag(P_mean)-TranspMult(P_mean,P_mean);
        double nan=std::numeric_limits<double>::quiet_NaN();
       return mp_state_information(std::move(P_mean),std::move(P_cov),nan,nan,nan,nan);
    }
    
    
    
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





template<class F, class MacroDR,class Model,template<class, class > class Experiment, class Point, class Measure>
auto logLikelihood_experiment_calculation(const F& f,const MacroDR& a,  Model& m,const Experiment<Point,Measure>& e)
{
    
    auto first_step=*e.begin_begin();
    auto prior=a.start(m,first_step);

    auto out=f(e);
    std::size_t i=0;
    for (auto it=e.begin_begin(); it!=e.end_end(); ++it)
    {
        auto post=a.run(prior,m,*it);
         f(post,out,i);
         prior=std::move(post);
    }
    return out;
    
}
template<class MacroDR,class Model,class Experiment>
double logLikelihood(const MacroDR& a, Model& m, const Experiment e )
{
    return logLikelihood_experiment_calculation(logLikelihood_function(),a,m,e);
}

template<class MacroDR,class Model,class Experiment>
std::vector<double> partialLogLikelihood(const MacroDR& a, Model& m, const Experiment e )
{
    return logLikelihood_experiment_calculation(partialLogLikelihood_function(),a,m,e);
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
    auto v= logLikelihood_experiment_calculation(partialLogLikelihood_monitor_function<measure_likelihood<double>>(),a,m,e);
    return experiment::basic_Experiment<experiment::point<double,double>,measure_likelihood<double>>(e,std::move(v));
}







} // namespace markov

#endif // LIKELIHOOD_MARKOV_PROCESS_H
