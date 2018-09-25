#ifndef MYEVIDENCE_H
#define MYEVIDENCE_H
#include "Matrix.h"
#include "myDistributions.h"
#include "myparameters.h"
#include <iomanip>
#include "mylikelihood.h"

namespace evidence {
template<class Model, class Parameters_distribution>
class Prior_Model
{
public:
    typedef  Prior_Model self_type;
    typedef  Cs<Model> template_types;
    constexpr static auto const className=my_static_string("Prior_Model")+my_trait<template_types>::className;

    typedef M_Matrix<double> Parameters;

    Parameters sample(std::mt19937_64& mt)const { return prior_.sample(mt);}
    myOptional_t<logLikelihood> compute_Likelihood(const M_Matrix<double>& x)const
    {
        typedef myOptional_t<logLikelihood> Op;
        double logL=prior_.logP(x);
        if (std::isfinite(logL))
        {
            double elogL=prior_.expected_logP();
            double vlogL=prior_.variance_logP();
            return Op(logLikelihood(logL,elogL,vlogL));
        }
        else return Op(false,"not a finite value="+ToString(logL));
    }
    myOptional_t<DlogLikelihood> compute_DLikelihood(const M_Matrix<double>& x)const
    {
        typedef myOptional_t<DlogLikelihood> Op;
        double logL=prior_.logP(x);
        double elogL=prior_.expected_logP();
        double vlogL=prior_.variance_logP();

        auto G=prior_.dlogL_dx(x);
        auto H=prior_.dlogL_dx2(x);
        std::stringstream ss;
        if (are_finite<true,double>().test(logL,ss)&&
                are_finite<true,M_Matrix<double>>().test(G,ss)&&
                are_finite<true,M_Matrix<double>>().test(H,ss)
                )
            return Op(DlogLikelihood(logL,elogL,vlogL,std::move(G),std::move(H)));
        else
            return Op(false, ss.str() );
    }



    Prior_Model(const Parameters_distribution& prior): prior_{prior}{}
private:
    Parameters_distribution prior_;

};


class ThlogLikelihood: public logLikelihood
{
public:
    typedef   logLikelihood base_type;
    typedef  ThlogLikelihood self_type;
    constexpr static auto const className=my_static_string("ThlogLikelihood")+base_type::className;

    logLikelihood prior()const { return prior_;}
    logLikelihood likelihood()const { return lik_;}
    ThlogLikelihood(logLikelihood prior, logLikelihood lik, double beta):
        logLikelihood(prior.logL()+beta*lik.logL(),prior.elogL()+beta*lik.elogL(),prior.vlogL()+beta*lik.vlogL())
      ,prior_{prior},lik_{lik}{}
    ThlogLikelihood()=default;
    void set_beta(double beta)
    { set_logL(prior().logL()+beta*likelihood().logL(),
               prior().elogL()+beta*likelihood().elogL(),
               prior().vlogL()+beta*likelihood().vlogL());

    }
private:
    logLikelihood prior_;
    logLikelihood lik_;

};

class ThDlogLikelihood: public DlogLikelihood
{
public:
    typedef  DlogLikelihood base_type;

    constexpr static auto const className=my_static_string("ThDlogLikelihood_")+base_type::className;
    //std::string myClass()const  { return className.str();}

    typedef   ThDlogLikelihood self_type ;
    static auto get_constructor_fields()
    {
        return std::make_tuple(
                    grammar::field(C<self_type>{},"prior",&self_type::prior),
                    grammar::field(C<self_type>{},"likelihood",&self_type::likelihood),
                    grammar::field(C<self_type>{},"beta",&self_type::beta)
                    );
    }

    DlogLikelihood prior()const { return prior_;}
    DlogLikelihood likelihood()const { return lik_;}
    double beta()const {return beta_;}
    ThDlogLikelihood(const DlogLikelihood& prior,const  DlogLikelihood& lik, double beta):
        DlogLikelihood(prior.logL()+beta*lik.logL(),
                       prior.elogL()+beta*lik.elogL(),
                       prior.vlogL()+beta*lik.vlogL(),
                       prior.G()+lik.G()*beta,prior.H()+lik.H()*beta),prior_{prior},lik_{lik}, beta_{beta}{}
    ThDlogLikelihood(DlogLikelihood&& prior,  DlogLikelihood&& lik, double beta):
        DlogLikelihood(prior.logL()+beta*lik.logL(),
                       prior.elogL()+beta*lik.elogL(),
                       prior.vlogL()+beta*lik.vlogL(),
                       prior.G()+lik.G()*beta,
                       prior.H()+lik.H()*beta),prior_{std::move(prior)},lik_{std::move(lik)}, beta_{beta}{}
    ThDlogLikelihood()=default;
    void set_beta(double beta){
        beta_=beta;
        set_logL(prior().logL()+beta*likelihood().logL(),
                 prior().elogL()+beta*likelihood().elogL(),
                 prior().vlogL()+beta*likelihood().vlogL());
        set_G(prior().G()+likelihood().G()*beta);
        set_H(prior().H()+likelihood().H()*beta);
    }


private:
    DlogLikelihood prior_;
    DlogLikelihood lik_;
    double beta_;

};


template<class PriorModel, class LikelihoodModel>
class Thermodynamic_Model
{
    static_assert (std::is_same_v<typename PriorModel::Parameters,typename LikelihoodModel::Parameters > );
public:

    typedef  Thermodynamic_Model self_type;
    typedef  Cs<PriorModel, LikelihoodModel> template_types;
    constexpr static auto const className=my_static_string("Thermodynamic_Model")+my_trait<template_types>::className;



    typedef typename PriorModel::Parameters Parameters;

    Parameters sample(std::mt19937_64& mt)const {return prior().sample(mt);}

    typedef ThDlogLikelihood DLikelihoodResult;
    typedef ThlogLikelihood LikelihoodResult;
    Thermodynamic_Model(const PriorModel& prior, const LikelihoodModel& lik, double beta):prior_{prior}, lik_{lik}, beta_{beta}{}
    Thermodynamic_Model(const Thermodynamic_Model& other):prior_{other.prior_}, lik_{other.lik_}, beta_{other.beta_}{};
    Thermodynamic_Model(Thermodynamic_Model&& other):prior_{other.prior_}, lik_{other.lik_}, beta_{std::move(other.beta_)}{};
    Thermodynamic_Model& operator=(Thermodynamic_Model&&)=default;
    Thermodynamic_Model& operator=(const Thermodynamic_Model&other){
        if (&other!=this)
        {
            Thermodynamic_Model tmp(other);
            std::swap(*this,tmp);
            return *this;
        };
    }


    double beta()const {return beta_;}

    void set_beta(double _beta){beta_=_beta;}
    auto   compute_DLikelihood(const Parameters& x, std::ostream& os) const
    {
        typedef myOptional_t<ThDlogLikelihood> Op;
        auto p=prior_.compute_DLikelihood(x);
        if (!p.has_value()) return Op(false, "invalid prior DlogLik :"+p.error());
        auto l=lik_.compute_DLikelihood(x,os);
        if (!l.has_value())
            return Op(false,"fails in likelihood :"+l.error());
        else
            return Op(ThDlogLikelihood (std::move(p).value(),std::move(l).value(), beta()));
    }

    auto compute_Likelihood(const Parameters& x) const
    {
        typedef myOptional_t<ThlogLikelihood> Op;
        auto p=prior_.compute_Likelihood(x);
        if (!p.has_value()) return Op(false, "invalid prior logLik :"+p.error());
        auto l=lik_.compute_Likelihood(x);
        if (!l.has_value)
            return Op(false,"fails in likelihood :"+l.error());
        else
            return Op(ThlogLikelihood (std::move(p).value(),std::move(l).value(), beta()));
    }

    PriorModel const & prior()const { return prior_;};
    LikelihoodModel const & likelihood() const { return lik_;}

private:
    PriorModel  prior_;
    LikelihoodModel lik_;
    double beta_;
};


template<class PriorModel, class LikelihoodModel>
class Thermodynamic_Model_Series
{
    static_assert (std::is_same_v<typename PriorModel::Parameters,typename LikelihoodModel::Parameters > );

    typedef  Thermodynamic_Model_Series self_type;
    typedef  Cs<PriorModel, LikelihoodModel> template_types;
    constexpr static auto const className=my_static_string("Thermodynamic_Model_Series")+my_trait<template_types>::className;



public:
    typedef typename PriorModel::Parameters Parameters;

    typedef ThDlogLikelihood DLikelihoodResult;
    typedef ThlogLikelihood LikelihoodResult;
    Thermodynamic_Model_Series(const PriorModel& prior, const LikelihoodModel& lik, std::vector<double> beta):
        prior_{prior}, lik_{lik}, models_{getModels(prior_,lik_,beta)}{}

    double beta(std::size_t i)const {return models_[i].beta();}

    std::size_t size()const {return models_.size();}

    void set_beta(double abeta, std::size_t i){models_[i].set_beta(abeta);}

    Thermodynamic_Model<PriorModel,LikelihoodModel>& model(std::size_t i){ return models_[i];}
    Thermodynamic_Model<PriorModel,LikelihoodModel>const & model(std::size_t i)const { return models_[i];}

    Thermodynamic_Model_Series(const Thermodynamic_Model_Series&)=delete ;




private:
    PriorModel prior_;
    LikelihoodModel lik_;
    std::vector<Thermodynamic_Model<PriorModel,LikelihoodModel>> models_;

    static std::vector<Thermodynamic_Model<PriorModel,LikelihoodModel>>
    getModels(const PriorModel& prior, const LikelihoodModel& lik, std::vector<double> beta)
    {
        std::vector<Thermodynamic_Model<PriorModel,LikelihoodModel>> out;
        for (std::size_t i=0; i<beta.size(); ++i)
            out.push_back(Thermodynamic_Model<PriorModel,LikelihoodModel>(prior,lik,beta[i]));
        return out;
    }
};







template<class Model>
class  LevenbergMarquardt
{
public:
    typedef typename Model::Parameters Parameters;

    class myDistribution: public Normal_Distribution<Parameters>
    {
    public:
        typedef Normal_Distribution<Parameters> base_type;
        constexpr static auto const className=my_static_string("Normal_Distribution_landa")+my_trait<Parameters>::className;
        std::string myClass()const override { return className.str();}
        virtual myDistribution* clone()const override{ return new myDistribution(*this);};

        typedef   myDistribution self_type ;
        static auto get_constructor_fields()
        {
            return std::make_tuple(
                        grammar::field(C<self_type>{},"mean",&base_type::mean),
                        grammar::field(C<self_type>{},"Cov",&base_type::Cov),
                        grammar::field(C<self_type>{},"landa",&self_type::landa)

                        );
        }
        myDistribution(Normal_Distribution<Parameters>&& d,
                       double landa):
            base_type{std::move(d)}, landa_{landa}{}

        double landa()const {return landa_;}
        myDistribution()=default;
    private:
        double landa_;
    };



    typedef myDistribution Distribution;
    typedef typename Model::DLikelihoodResult LikelihoodResult;
    typedef  LevenbergMarquardt self_type;
    typedef  Cs<Model> template_types;
    constexpr static auto const className=my_static_string("LevenbergMarquardt")+my_trait<template_types>::className;

    //std::string myClass()const  { return className.str();}



    struct myAcceptProb
    {
        constexpr static auto const className=my_trait<LevenbergMarquardt>::className+my_static_string("_myAcceptProb");
        static std::tuple<> get_constructor_fields() {return std::tuple<>();}

        typedef  std::pair<double, double> Parameters;
        double operator()(const LevenbergMarquardt<Model>& LM,const Parameters& param)const
        {
            double landa50=param.first;
            double h=param.second;
            return 1.0/(1.0+std::pow(landa50/(LM.landa()+1.0),h));
        }
        double operator()(const Parameters& param, const LevenbergMarquardt<Model>& LM)const
        {
            return operator()(LM,param);
        }
        static  std::map<Parameters, double> uniform_parameter_prior(const  std::vector<std::vector<double>>& v, double p=-1)
        {
            std::map<Parameters, double> out;
            if (p==-1)
                p=1.0/(v[0].size()*v[1].size());
            for (std::size_t i=0; i<v[0].size(); ++i)
                for (std::size_t j=0; j<v[1].size(); ++j)
                    out[{v[0][i],v[1][j]}]+=p;
        return out;
    }
};
    typedef myAcceptProb AcceptanceProbability;

    struct myExpectVelocity
    {
        constexpr static auto const className=my_trait<LevenbergMarquardt>::className+my_static_string("_myExpectVelocity");
        static std::tuple<> get_constructor_fields() {return std::tuple<>();}

        double operator()(const LevenbergMarquardt<Model>& LM)const
        {
            return (1.0/(1.0+LM.landa()));
        }
    };

    typedef myExpectVelocity ExpectedVelocity;


    double landa()const { return landa_;}



    myOptional<Distribution,reg_tag> calculate_Distribution(const Model& ,const Parameters x, const LikelihoodResult& logL ) const
    {
        typedef myOptional<Distribution,reg_tag> Op;
        auto& H=logL.H();
        auto& G=logL.G();
        auto inv_cov=-(H+diag(H)*landa_);
        auto cov=inv(inv_cov);
        if (cov)
        {
            Parameters d=(G*cov.value());
            Parameters newpoint=x+d;
            auto chol_cov=chol(cov.value(),"lower");
            if (chol_cov.has_value())
            {
                return Op(Distribution(typename Distribution::base_type(newpoint,cov.value(),inv_cov,chol_cov.value()),landa_));
            }
            else return Op(false," invalid cholesky decomposition: "+chol_cov.error());

        }
        else return Op(false,"cannot invert Hessian landa matrix "+ cov.error());
    }

    auto compute_Likelihood(const Model& m, const Parameters& x, std::ostream& os) const
    {
        return m.compute_DLikelihood(x,os);
    }

    LevenbergMarquardt(double _landa):landa_{_landa}{}
    LevenbergMarquardt()=default;

    bool operator<(const LevenbergMarquardt& other )const { return landa()<other.landa();}

    static auto get_constructor_fields()
    {
        return std::make_tuple(
                    grammar::field(C<self_type>{},"landa",&self_type::landa)
                    );
    }



private:
    double landa_;

};


template<class Model>
class  AdaptiveCovariance
{
public:
    typedef typename Model::Parameters Parameters;

    typedef Normal_Distribution<Parameters> Distribution;
    typedef typename Model::LikelihoodResult LikelihoodResult;

    typedef  AdaptiveCovariance self_type;
    typedef  Cs<Model> template_types;
    constexpr static auto const className=my_static_string("AdaptiveCovariance")+my_trait<template_types>::className;


    Distribution getDistribution(const Model& ,const Parameters& x, const LikelihoodResult&  ) const
    {
        return Normal_Distribution<Parameters>(x,cov_);
    }

    LikelihoodResult compute_Likelihood(const Model& m, const Parameters& x)
    {
        return m.compute_Likelihood(x);
    }

    AdaptiveCovariance(double r,const M_Matrix<double>& c):cov_{c*r}, ratio_(r){}

    double ratio()const {return ratio_;}
private:
    M_Matrix<double> cov_;
    double ratio_;

};



template<class Model,class aLikelihoodResult,class  aParameters,class aDistribution>
class  mcmc_sample: public aLikelihoodResult
{
public:
    typedef   aParameters  Parameters;
    typedef aDistribution Distribution;
    typedef aLikelihoodResult LikelihoodResult;

    typedef  aLikelihoodResult base_type;

    typedef  mcmc_sample self_type;
    typedef  Cs<Model,aLikelihoodResult,aParameters,aDistribution> template_types;
    constexpr static auto const className=my_static_string("AdaptiveCovariance")+my_trait<base_type>::className+my_trait<template_types>::className;



    std::string myClass()const  { return className.str();}

    static auto get_constructor_fields()
    {
        return std::make_tuple(
                    grammar::field(C<self_type>{},"Parameters",&self_type::x),
                    grammar::field(C<self_type>{},"likelihood_result",&self_type::likelihood_result),
                    grammar::field(C<self_type>{},"Distribution",&self_type::g)
                    );
    }

    LikelihoodResult const & likelihood_result()const { return *this;}
    Parameters const & x() const {return x_;}
    Distribution const & g()const {return g_;}
    bool accept()const { return accepted_;}
    void push_rejected(){accepted_=false;}
    void push_accepted(){accepted_=true;}

    void set_g(Distribution&& new_g){g_=std::move(new_g);}
    mcmc_sample(const Parameters& par,LikelihoodResult&& lik,Distribution&& dist)
        :LikelihoodResult(std::move(lik)),x_{par},g_{std::move(dist)},accepted_{true}{}


    mcmc_sample()=default;



private:
    Parameters x_;
    Distribution g_;
    bool accepted_;

};


template<class Model , class Distribution_Generator>
using mcmc_sample_t=mcmc_sample<Model,typename Distribution_Generator::LikelihoodResult,typename Model::Parameters, typename Distribution_Generator::Distribution>;


template<class Model,class Distribution_Generator, class Parameters=typename Model::Parameters>
static
myOptional_t<mcmc_sample_t<Model,Distribution_Generator>>
calculate(const Model& M, const Distribution_Generator& G, Parameters&& x, std::ostream& os)
{
    typedef myOptional_t<mcmc_sample_t<Model,Distribution_Generator>> Op;

    auto logL=G.compute_Likelihood(M,x,os);
    os<<"\nlogL.likelihood().H()\n"<<logL.value().likelihood().H();
    os<<"\nlogL.H()\n"<<logL.value().H();
    //   std::cerr<<" gets a logL!!----------="<<logL<<"\n";
    if (!logL.has_value())
        return Op(false,"cannot compute logLikelihood because "+logL.error());
    auto g=G.calculate_Distribution(M,x,logL.value());
    if (!g)
        return Op(false,"cannot compute distribution because "+g.error());
    return Op(mcmc_sample_t<Model,Distribution_Generator>(std::move(x),std::move(logL).value(),std::move(g).value()));
}

template<bool Hastings>
constexpr static auto str_Hastings(){
    if constexpr (Hastings) return my_static_string("_Hastings");
    else return my_static_string("");
}

template<bool Hastings=true>
class Metropolis_H
{
public:


    typedef  Metropolis_H self_type;
    constexpr static auto const className=my_static_string("Metropolis")+str_Hastings<Hastings>();

    template<class mcmcsample, class stream>
    static double Acceptance(const mcmcsample& current, const mcmcsample& candidate, stream& s)
    {
        double logL_current=double(current.logL());
        double logL_candidate=double(candidate.logL());
        if constexpr (Hastings)
        {
            double logP_forward=current.g().logP(candidate.x());
            double logP_backward=candidate.g().logP(current.x());
            double logA=(logL_candidate+logP_backward)-(logL_current+logP_forward);
            double A=std::min(1.0,exp(logA));
            s<<"\n logA=(logL_candidate+logP_backward)-(logL_current+logP_forward)\n logA="<<logA<<"\n";
            s<<"\nlogL_candidate="<<logL_candidate<<"\n";
            s<<"\nlogL_current="<<logL_current<<"\n";

            s<<"\nlogP_foward="<<logP_forward<<"\n";
            s<<"\nlogP_backward="<<logP_backward<<"\n";
            s<<"\ncandidate.beta()="<<candidate.beta()<<"\n";
            s<<"\ncurrent.beta()="<<current.beta()<<"\n";

            s<<"\ncandidate.g().landa()="<<candidate.g().landa()<<"\n";
            s<<"\ncurrent.g().landa()="<<current.g().landa()<<"\n";


            s<<"\n logDetCov_forward ="<<current.g().logDetCov()<<"\n";
            s<<"\n logDetCov_backward ="<<candidate.g().logDetCov()<<"\n";

            auto dcandidate=candidate.x()-candidate.g().mean();
            auto dforward=candidate.x()-current.g().mean();
            auto dback=current.x()-candidate.g().mean();
            auto dcurrent=current.x()-current.g().mean();

            s<<"\ndcandidate\n"<<dcandidate;
            s<<"\ndforward\n"<<dforward;
            s<<"\n dback\n"<<dback;
            s<<"\n dcurrent\n"<<dcurrent;

            s<<"\nxTSigmaX(dforward, current.g().CovInv())=\n"<<xTSigmaX(current.g().CovInv(),dforward)<<"\n";

            s<<"\nxTSigmaX(dforward, current.g().CovInv())=\n"<<xTSigmaX(dforward, current.g().CovInv())<<"\n";
            s<<"\nxTSigmaX(dcurrent, current.g().CovInv())=\n"<<xTSigmaX(dcurrent, current.g().CovInv())<<"\n";
            s<<"\n INV xTSigmaX(dcandidate, current.g().CovInv())=\n"<<xTSigmaX(dcandidate, current.g().CovInv())<<"\n";
            s<<"\n INV xTSigmaX(dback, current.g().CovInv())=\n"<<xTSigmaX(dback, current.g().CovInv())<<"\n";

            s<<"\nxTSigmaX(dforward, candidate.g().CovInv())=\n"<<xTSigmaX(dforward, candidate.g().CovInv())<<"\n";
            s<<"\nxTSigmaX(dcurrent, candidate.g().CovInv())=\n"<<xTSigmaX(dcurrent, candidate.g().CovInv())<<"\n";
            s<<"\n INV xTSigmaX(dcandidate, candidate.g().CovInv())=\n"<<xTSigmaX(dcandidate, candidate.g().CovInv())<<"\n";
            s<<"\n INV xTSigmaX(dback, candidate.g().CovInv())=\n"<<xTSigmaX(dback, candidate.g().CovInv())<<"\n";


            s<<"\nquadraticForm_B_A_BT( current.g().CovInv(),dforward)=\n"<<quadraticForm_B_A_BT( current.g().CovInv(),dforward)<<"\n";

            s<<"\nquadraticForm_B_A_BT(candidate.g().CovInv(),dback)=\n"<<quadraticForm_B_A_BT(candidate.g().CovInv(),dback)<<"\n";


            s<<"\n current.x() \n"<<current.x();
            s<<"\n current.g().mean() \n"<<current.g().mean();


            s<<"\n candidate.x() \n"<<candidate.x();
            s<<"\n candidate.g().mean() \n"<<candidate.g().mean();

            s<<"\n current.x()-candidate.x() \n"<<current.x()-candidate.x();
            s<<"\n ds=-current.g().mean()+candidate.x() \n"<<-current.g().mean()+candidate.x();
            s<<"\n current.g().mean()-candidate.g().mean() \n"<<current.g().mean()-candidate.g().mean();


            s<<"\ncurrent.g().Chol()\n"<<current.g().Chol();


            s<<"\ncandidate.g().Chol()\n"<<candidate.g().Chol();

            s<<"\ncurrent.g().Chol()-candidate.g().Chol()\n"<<current.g().Chol()-candidate.g().Chol();

            s<<"\ncurrent.g().CovInv()\n"<<current.g().CovInv();
            s<<"\ncandidate.g().CovInv()\n"<<candidate.g().CovInv();


            s<<"\n current \n"<<current;
            s<<"\n candidate \n"<<candidate;
            s<<"-------------------------------------------------------------------------\n";




            return A;
        }
        else
        {
            double logA=(logL_candidate)-(logL_current);
            double A=std::min(1.0,exp(logA));
            return A;
        }
    }

    template<class Model, class Distribution_Generator, class Parameters=typename Model::Parameters>
    myOptional_t<mcmc_sample_t<Model,Distribution_Generator>>
    start(const Model& M,std::mt19937_64& mt, const Distribution_Generator& G, std::ostream& os) const
    {
        std::size_t n=0;
        typedef  myOptional_t<mcmc_sample_t<Model,Distribution_Generator>> Op;

        //  std::cerr<<"\nmax number of trials="<<max_number_of_trials();
        while(n<max_number_of_trials())
        {

            Parameters x=M.sample(mt);
            //    std::cerr<<"\n x="<<x;
            auto out=calculate(M,G,x,os);
            //  std::cerr<<" finishela!!";
            if (out.has_value())
            {
                //      std::cerr<<" and succeds"<<out.value();
                return out;
            }
            else
            {
                std::cerr<<out.error();
                ++n;
            }
        }
        return Op(false,"fails to start after n="+ToString(n)+" trials");
    }

    template<class Model, class Distribution_Generator, class stream,class mySample=mcmc_sample_t<Model,Distribution_Generator>>
    Op_void next(const Model& M,std::mt19937_64& mt,const Distribution_Generator& G,mySample& current, stream& s)const
    {
        auto candidate_point=current.g().sample(mt);
        // current.g().autoTest(mt,900);
        myOptional_t<mySample> candidate=calculate(M,G,candidate_point,s);
        if (!candidate)
        {
            current.push_rejected();
        }
        else
        {
            // auto gg=G.calculate_Distribution(M,candidate.value().x(), candidate.value());
            //            s<<"gg.mean()"<<gg.value().mean();
            //            s<<"\ncandidate.g().mean()\n"<<candidate.value().g().mean();
            //            s<<"gg.Cov()"<<gg.value();
            //            s<<"\ncandidate.g()\n"<<candidate.value().g();
            //            s<<"\ncandidate.g().Cov()- gg.value().Cov()\n"<<candidate.value().g().Cov()- gg.value().Cov();
            double A=Acceptance(current,candidate.value(),s);
            std::uniform_real_distribution<double> u(0,1);
            double r=u(mt);
            bool accept=r<A;
            if (accept)
            {
                current=std::move(candidate.value());
                current.push_accepted();
            }
            else current.push_rejected();
        }
        return Op_void(true,"");

    }
    std::size_t max_number_of_trials()const {return ntrials_;}
    Metropolis_H<Hastings>(std::size_t max_trials):ntrials_{max_trials}{}
private:
    std::size_t ntrials_;

};


typedef   Metropolis_H<true> Metropolis_Hastings;

typedef Metropolis_H<false> Metropolis;

template< bool Hastings=true>
class Adaptive_Metropolis_H
{
public:

    typedef  Adaptive_Metropolis_H self_type;
    constexpr static auto const className=my_static_string("Adaptive_Metropolis")+str_Hastings<Hastings>();


    template<class Model, class Adapative_Distribution_Generator,
             class mySample=myOptional_t<mcmc_sample_t<Model,Adapative_Distribution_Generator>>>
    mySample start(const Model& M,std::mt19937_64& mt, Adapative_Distribution_Generator& adaptive, std::ostream& os)const
    {
        auto G=adaptive.sample(mt);
        std::cerr<<" start "<<G;
        auto   out=get_Metropolis().start(M,mt,G,os);
        return out;
    }




    template<class Model, class Adapative_Distribution_Generator,class stream,
             class mySample=mcmc_sample_t<Model,Adapative_Distribution_Generator>>
    myOptional_t<void> next(const Model& M,std::mt19937_64& mt, Adapative_Distribution_Generator& adaptive, mySample& current, stream& s)const
    {
        auto G=adaptive.sample(mt);
        auto gg=G.calculate_Distribution(M,current.x(), current.likelihood_result());
        if (!gg)
        {
            return myOptional_t<void>(false,"cannot calculate distribution "+gg.error());
        }
        else
        {
            //    std::cerr<<"\nold current.g()\n"<<current.g();

            //  auto gold=gg.value();
            current.set_g(std::move(gg).value());
            //   s<<"\n new current.g().mean()-ggold.mean()\n"<<current.g().mean()-gold.mean();
            //   s<<"\n new current.g().Cov()-ggold.COv()\n"<<current.g().Cov()-gold.Cov();

            get_Metropolis().next(M,mt,G,current,s);
            adaptive.push_outcome(G,current, current.accept());
            return myOptional_t<void>(true,"");
        }
    }

    Metropolis_H<Hastings> const & get_Metropolis()const {return mh_;}

    Adaptive_Metropolis_H(Metropolis_H<Hastings> mh): mh_{std::move(mh)}{}
private:
    Metropolis_H<Hastings> mh_;
};


typedef Adaptive_Metropolis_H<true> Adaptive_Metropolis_Hastring;
typedef Adaptive_Metropolis_H<false> Adaptive_Metropolis;


template<class Model, class Adaptive>
class emcee_sample
{
public:


    typedef  emcee_sample self_type;
    typedef  Cs<Model,Adaptive> template_types;
    constexpr static auto const className=my_static_string("emcee_sample")+my_trait<template_types>::className;


    typedef typename Adaptive::Mover Mover;
    typedef mcmc_sample_t<Model,Mover> mcmc_s;
    mcmc_s& Walker(std::size_t i){ return walkers_[i];}
    mcmc_s const & Walker(std::size_t i)const { return walkers_[i];}
    std::size_t numWalkers()const { return walkers_.size();}
    std::mt19937_64& mt(std::size_t i){ return mt_[i];}


    emcee_sample(const Model& m,std::mt19937_64& mt,Adaptive& mov,std::size_t nwalkers):
        mt_{mts(mt,nwalkers)}, walkers_{getWalkers(m,mt,mov,nwalkers)}{}

private:
    std::vector<std::mt19937_64> mt_;
    std::vector<mcmc_s> walkers_;
    static std::vector<std::mt19937_64> mts(std::mt19937_64& mt,std::size_t n)
    {
        std::uniform_int_distribution<typename std::mt19937_64::result_type> useed;
        std::vector<std::mt19937_64> out(n);
        for (std::size_t i=0; i<n; ++i)
            out[i].seed(useed(mt));
        return out;
    }

    static std::vector<mcmc_s> getWalkers(const Model& model,std::mt19937_64& mt, const Adaptive& adaptive,std::size_t n)
    {
        std::vector<mcmc_s> out(n);
        for (std::size_t i=0; i<n; ++i)
        {
            auto x=model.sample(mt);
            auto logL=Mover::compute_Likelihood(model,x);
            auto d=adaptive.sample(mt);
            out[i]=mcmc_s(logL,x,d);
        }

    }


};


class Ensemble_Metropolis_Hastings
{
private:
    std::size_t numWalkers_;
public:
    typedef  Ensemble_Metropolis_Hastings self_type;
    constexpr static auto const className=my_static_string("Ensemble_Metropolis_Hastings");


    std::size_t numWalkers()const {return numWalkers_;}
    Ensemble_Metropolis_Hastings(std::size_t numberWalkers):numWalkers_{numberWalkers}{}
    template<class Model>
    class stretch_move
    {
    public:
        typedef logLikelihood LikelihoodResult;
        typedef stretch_move_Distribution Distribution;
        typedef  typename Model::Parameters Parameters;


        static
        myOptional_t<LikelihoodResult> compute_Likelihood(const Model& m, Parameters& x){ return m.compute_Likelihood(x);}

        static double Acceptance(double Z, const logLikelihood& candidateLogLik, const logLikelihood& currentLogLik, std::size_t numParam)
        {
            double logA=candidateLogLik.logL()-currentLogLik.logL()+log(Z)*(numParam-1);
            return std::min(0.0,logA);
        }


        template<class Adaptive_Strecth_Move>
        static mcmc_sample_t<Model,stretch_move> move(const Model& model,
                                                      std::mt19937_64& mt,
                                                      const stretch_move_Distribution& d,
                                                      const emcee_sample<Model,Adaptive_Strecth_Move>& s,
                                                      const std::vector<std::size_t>& index,
                                                      mcmc_sample_t<Model,stretch_move>& current
                                                      )
        {
            std::uniform_int_distribution<std::size_t> u(0,index.size());
            std::size_t i=index[u(mt)];
            double z=d.sample(mt);
            Parameters candidate=(current.x()-(s.Walker(i).x()))*z+s.Walker(i).x();
            LikelihoodResult lik=compute_Likelihood(model,candidate);
            auto logA=test(z,lik,current,candidate.size());
            double A=std::exp(logA);
            std::uniform_real_distribution<double> real(0,1);
            double r=real(mt);
            bool accept=r<A;
            if (accept)
            {
                current=mcmc_sample_t<Model, stretch_move>(std::move(lik), std::move(candidate), std::move(current));
                current.push_accepted();
            }
            else
                current.push_rejected();



        }


    };

    template< class Model>
    class adaptive_stretch_mover
    {
    public:

        typedef stretch_move<Model> Mover;
        stretch_move_Distribution sample(std::mt19937_64& mt)
        {
            return stretch_move_Distribution(alfa_.sample(mt));
        }

        void push_outcome(const stretch_move_Distribution& d, bool accept)
        {
            if (accept) alfa_map_[d.alpha()].push_accept();
            else alfa_map_[d.alpha()].push_reject();

        }

        void actualize()
        {
            alfa_=alfa_map_.Distribute_on_p(target_);
        }

    private:
        double target_;

        Beta_map<double> alfa_map_;
        Probability_map<double> alfa_;



    };

    static
    std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
    random_split(std::mt19937_64& mt, std::size_t n)
    {
        std::vector<std::size_t> index(n);
        for (std::size_t i=0; i<n; ++i) index[i]=i;

        std::shuffle(index.begin(),index.end(),mt);
        std::size_t half=n/2;

        std::vector<std::size_t> one(index.begin(), index.begin()+half);
        std::vector<std::size_t> two(index.begin()+half, index.end());
        return {one, two};
    }

    template<class Model, class Adaptive_Mover, class Mover=typename Adaptive_Mover::Mover>
    emcee_sample<Model,Mover> start(const Model& model,std::mt19937_64& mt, Adaptive_Mover& adaptive)
    {
        return emcee_sample<Model,Mover>(model,mt,adaptive,numWalkers());
    }


    template<class Model, template<class >class Adaptive_Mover, class Mover=typename Adaptive_Mover<Model>::Mover>
    bool next(const Model& model,std::mt19937_64& mt, adaptive_stretch_mover<Model>& adaptive,emcee_sample<Model,stretch_move<Model>>& current)
    {

        //   typedef typename Adaptive_Mover<Model>::Mover stretch_move;
        auto [one,two]=random_split(current.numWalkers());
        #pragma omp parallel for
                for (std::size_t ii=0; ii<one.size(); ++ii)
        {
            auto i=one[ii];
            current.Walker(i).set_g(adaptive.sample(mt));
            Mover::move(model,current.mt(i),current.Walker(i).g(),current,one,current.Walker(i));
        }
        for (std::size_t ii=0; ii<one.size(); ++ii)
        {
            auto i=one[ii];
            adaptive.push_outcome(current.Walker(i).g(),current.accept(i));
        }

#pragma omp parallel for
        for (std::size_t ii=0; ii<two.size(); ++ii)
        {
            auto i=two[ii];
            current.Walker(i).g()=adaptive.sample(mt);
            stretch_move<Model>::move(model,current.mt(i),current.Walker(i).g(),current,two,current.Walker(i));
        }
        for (std::size_t ii=0; ii<two.size(); ++ii)
        {
            auto i=two[ii];
            adaptive.push_outcome(current.Walker(i).g(),current.accept(i));
        }
        adaptive.actualize();
        return current;
    }


};



template <class Parameterized_Distribution/*=Landa*/,
          class Gain=typename Parameterized_Distribution::ExpectedVelocity,
          class Likelihood=typename Parameterized_Distribution::myAcceptProb>
class Adaptive_Parameterized_Distribution_Generator
{
public:
    typedef typename Parameterized_Distribution::LikelihoodResult LikelihoodResult;
    typedef typename Parameterized_Distribution::Distribution Distribution;

    typedef  Adaptive_Parameterized_Distribution_Generator self_type;
    typedef  Cs<Parameterized_Distribution,Gain,Likelihood> template_types;
    constexpr static auto const className=my_static_string("Adaptive_Parameterized_Distribution_Generator")+my_trait<template_types>::className;

    static auto get_constructor_fields()
    {
        return std::make_tuple(
                    grammar::field(C<self_type>{},"likelihoodFunction",&self_type::getLikelihood),
                    grammar::field(C<self_type>{},"gainFunction",&self_type::getGain),
                    grammar::field(C<self_type>{},"gainMoment",&self_type::getGainMoment),
                    grammar::field(C<self_type>{},"probabilityMap",&self_type::get_ProbabilityMap),
                    grammar::field(C<self_type>{},"likelihoodProbabilityMap",&self_type::get_Parameters_Map),
                    grammar::field(C<self_type>{},"logEvidence",&self_type::getLogEvidence)

                    );
    }


    Parameterized_Distribution sample(std::mt19937_64& mt)const{return p_.sample(mt);}


    template<class mcmc>
    void push_outcome(const Parameterized_Distribution& landa, const mcmc&, bool accept){
        if (accept) push_acceptance(landa); else push_rejection(landa);
        actualize();
    }

    void push_acceptance(const Parameterized_Distribution& landa){
        parDist_=logBayes_rule(Log_of(lik_),landa,parDist_);}


    void push_rejection(const Parameterized_Distribution& landa){
        parDist_=logBayes_rule(Log_of(Complement_prob(lik_)),landa,parDist_);}

    void actualize(double nmax){
        actualize();
        parDist_.reduce(nmax);}

    static std::pair<Probability_map<Parameterized_Distribution>, double>
    Distribute_on_gain(const Gain& g, const Likelihood& lik,const  logLikelihood_map<typename Likelihood::Parameters>& par,
                       const Probability_map<Parameterized_Distribution> & landas, double moment)
    {
        auto out=landas.p();
        auto p_par=par.p();
        for (auto& e:out)
        {
            double meangain=Expectance(g,lik,p_par,e.first);
            e.second=std::pow(meangain,moment);
        }
        auto o=Probability_map<Parameterized_Distribution>::normalize(out,landas.nsamples());
        double sum=0;
        for (auto& e:o.first.p())
        {
            sum+=e.second*std::pow(out[e.first],1.0/moment);
        }
        return {o.first,sum};
    }

    void actualize(){
        auto pold=p_;
        auto o=Distribute_on_gain(g_,lik_,parDist_,p_,gainMoment_);
        p_=o.first;
    }


    /*Adaptive_Parameterized_Distribution_Generator(const std::map<Parameterized_Distribution,double>& prior_landa,
                                                  const std::map<typename Likelihood::Parameters, double>& prior_par,
                                                  double nsamples,std::size_t gainMoment):
        gainMoment_(gainMoment),
        p_{prior_landa,nsamples},
        parDist_(prior_par,nsamples){}
*/

    template<template<typename...>class V>
    Adaptive_Parameterized_Distribution_Generator(const V<Parameterized_Distribution>& landa,
                                                  const  std::vector<std::vector<double>>& par,
                                                  double gainMoment):
        lik_{},g_{},
        gainMoment_(gainMoment),
        p_{landa},
        parDist_{Likelihood::uniform_parameter_prior(par,0.0),0},logEvidence{0}{}


    Adaptive_Parameterized_Distribution_Generator(){}

    Adaptive_Parameterized_Distribution_Generator(Likelihood lik, Gain g, std::size_t gainMoment,const Probability_map<Parameterized_Distribution>& pmap,
                                                  const logLikelihood_map<typename Likelihood::Parameters>& pardist, double logEvide)
        : lik_{lik}, g_{g}, gainMoment_{gainMoment}, p_{pmap},parDist_{pardist}, logEvidence{logEvide}{}


    Likelihood const & getLikelihood()const {return lik_;}

    Gain const & getGain()const { return g_;}

    std::size_t getGainMoment()const {return gainMoment_;}
    Probability_map<Parameterized_Distribution> get_ProbabilityMap()const  {return p_;}
    logLikelihood_map<typename Likelihood::Parameters> get_Parameters_Map() const {return parDist_;}
    double getLogEvidence()const {return logEvidence;}



private:
    Likelihood lik_;
    Gain g_;
    std::size_t gainMoment_;
    Probability_map<Parameterized_Distribution> p_;
    logLikelihood_map<typename Likelihood::Parameters> parDist_;
    double logEvidence;
};


template <class Model>
class AdaptiveCovarianceGenerator

{
    typedef  AdaptiveCovarianceGenerator self_type;
    typedef  Cs<Model> template_types;
    constexpr static auto const className=my_static_string("AdaptiveCovarianceGenerator")+my_trait<template_types>::className;


    typedef typename Model::Parameters Parameters;
    AdaptiveCovariance<Model> sample(std::mt19937_64& mt)
    {
        double r=r_.sample(mt);
        return AdaptiveCovariance<Model>(r,cov_);
    }


    void actualize()
    {
        auto m=Sx_*(1.0/nsamples_);
        cov_=Sxx_*(1.0/nsamples_)-TranspMult(m,m);
        r_=rmap_.Distribute_on_p(target_);
    }

    void reset(std::size_t nmax)
    {
        double f=1.0*nmax/nsamples_;
        if (f<1)
        {
            nsamples_=nmax;
            Sx_*=f;
            Sxx_*=f;
            rmap_.reduce(nmax);
            r_=rmap_.Distribute_on_p(target_);
        }
    }
    void push_outcome(const AdaptiveCovariance<Model>& cov,const Parameters& p, bool accept)
    {
        Sxx_+=TranspMult(p,p);
        Sx_+=p;
        ++nsamples_;
        if (accept)
            rmap_[cov.ratio()].push_accept();
        else
            rmap_[cov.ratio()].push_reject();
        if (nsamples_ % nkip_ ==0)
        {
            actualize();
        }
    }

private:
    double target_;
    std::size_t nkip_;
    Beta_map<double> rmap_;
    std::size_t nsamples_;
    M_Matrix<double> Sxx_;
    M_Matrix<double> Sx_;
    M_Matrix<double> cov_;
    Probability_map<double> r_;
};


template < class EV,class Tp,class DistGen, class DistAdapt>
class Adaptive_Probability_Distribution_Generator
{
public:
    typedef  Adaptive_Probability_Distribution_Generator self_type;
    typedef  Cs< EV, Tp, DistGen,  DistAdapt> template_types;
    constexpr static auto const className=my_static_string("Adaptive_Probability_Distribution_Generator")+my_trait<template_types>::className;


    DistGen sample(std::mt19937_64& mt)const
    { return sample_rev_map(rev_,mt);}

    void push_outcome(const DistGen& landa, const M_Matrix<double>& ,bool accept)
    {
        if (accept) push_acceptance(landa);
        else push_rejection(landa);
        actualize();

    }

    void push_acceptance(DistGen landa){
        landaDist_[landa].push_accept();}


    void push_rejection(DistGen landa){
        landaDist_[landa].push_reject();}

    void actualize(double nmax)
    {
        actualize();
        landaDist_.reduce(nmax);
    }

    void actualize()
    {
        auto pnew=this->p_;

        double sum=0;
        for (auto it=pnew.begin(); it!=pnew.end(); ++it)
        {
            auto ns=landaDist_[it->first].Parameters();
            double l=f_(it->first)*tp_(ns);
            sum+=l;
            it->second=l;
        }
        double expectedGain=0;
        for (auto it=pnew.begin(); it!=pnew.end(); ++it)
        {
            auto ns=landaDist_[it->first].Parameters();
            it->second*=1.0/sum;
            expectedGain+=it->second*f_(it->first)*tp_(ns);
        }

        p_=pnew;
        this->rev_=cumulative_reverse_map(this->p_);

    }



    Adaptive_Probability_Distribution_Generator (const std::map<DistGen,double>& prior_landa):
        f_(),tp_(),
        p_{prior_landa},
        rev_{cumulative_reverse_map(p_)},
        landaDist_{Beta_map<DistGen>::UnInformativePrior(prior_landa)}{}

    Adaptive_Probability_Distribution_Generator(const Tp& tp,const std::map<DistGen,double>& prior_landa):
        f_(),tp_(tp),
        p_{prior_landa},
        rev_{cumulative_reverse_map(p_)},
        landaDist_{Beta_map<DistGen>::UnInformativePrior(prior_landa)}{}

    template<template<typename>class V>
    Adaptive_Probability_Distribution_Generator(const V<DistGen>& landa):
        Adaptive_Probability_Distribution_Generator(uniform_prior(landa)){}

    template<template<typename>class V>
    Adaptive_Probability_Distribution_Generator(const Tp tp,const V<DistGen>& landa):
        Adaptive_Probability_Distribution_Generator(tp,uniform_prior(landa)){}

    Adaptive_Probability_Distribution_Generator(){}




private:
    EV f_;
    Tp tp_;
    DistAdapt A_;
    std::map<DistGen,double> p_;

    std::map<double,DistGen> rev_;

    Beta_map<DistGen> landaDist_;

    template<template<typename>class V>
    static std::map<DistGen,double> uniform_prior(const V<DistGen>& v)
    {
        std::map<DistGen,double> out;
        double p=1.0/v.size();
        for (auto it=v.begin(); it!=v.end(); ++it)
            out[*it]+=p;
        return out;
    }



};


template <bool Hastings>
class Parallel_Tempering
{
public:

    typedef  Parallel_Tempering self_type;
    constexpr static auto const className=my_static_string("Parallel_Tempering")+str_Hastings<Hastings>();


    template<class priorModel, class LikelihoodModel,class Adapative_Distribution_Generator>
    class samples
    {
    public:
        samples()=default;
        typedef Thermodynamic_Model<priorModel,LikelihoodModel> Model;
        typedef Thermodynamic_Model_Series<priorModel,LikelihoodModel> ModelSeries;
        auto& Scout(std::size_t i){ return scouts_[i_to_scout[i]];}
        auto& Scout(std::size_t i)const {return scouts_[i_to_scout.at(i)];}
        void swap(std::size_t i, const ModelSeries& model)
        {
            std::swap(i_to_scout[i], i_to_scout[i+1]);
            Scout(i).set_beta(model.beta(i));
            Scout(i+1).set_beta(model.beta(i+1));

        }

        double Accept(std::size_t i, const ModelSeries& model)const
        {
            double logA=(model.beta(i)-model.beta(i+1))*
                    (Scout(i+1).likelihood().logL()-Scout(i).likelihood().logL());
            std::cerr<<"\n thermologA"<<logA<<"  "<<model.beta(i)<<" "<<model.beta(i+1)<<"  "<<Scout(i).likelihood().logL()<<" "<<Scout(i+1).likelihood().logL()<<"\n";
            return  std::min(1.0,exp(logA));

        }
        std::size_t size()const { return scouts_.size();}

        std::mt19937_64& mt(std::size_t i){return mt_[i];}


        static myOptional<samples,reg_tag> evaluate(const Adaptive_Metropolis_H<Hastings>& adm,
                                                    const Thermodynamic_Model_Series<priorModel,LikelihoodModel>& m,std::mt19937_64& _mt,
                                                    std::vector<Adapative_Distribution_Generator>& adaptive)
        {
            typedef myOptional<samples,reg_tag> Op;
            auto mymt=mts(_mt,m.size());
            auto itoscout=init_map(m.size());
            std::vector<mcmc_sample_t<Model,Adapative_Distribution_Generator>> myscouts(m.size());
            std::vector<Op_void> res(m.size());
            std::cerr<<"\n por empezar\n";
            std::vector<std::stringstream> os(m.size());

#pragma omp parallel for
            for (std::size_t i=0; i<m.size(); ++i)
            {
                std::cerr<<"\n"<<i <<" aqui\n";
                auto s=adm.start(m.model(i),mymt[i],adaptive[i],os[i]);
                //     std::cerr<<"\n"<<i <<" next\n";
                //     std::cerr<<"\n"<<i <<s.value()<<" value\n";

                res[i]=s;
                //     std::cerr<<"\n"<<i <<" res\n";
                myscouts[i]=std::move(s).value();
                //     std::cerr<<"\n"<<i <<" scouts\n";
            }
            for (auto& ss: os)
                std::cerr<<ss.str();
            auto ops=consolidate(std::move(res));
            if (ops.has_value())
                return Op(samples(std::move(mymt),std::move(itoscout), std::move(myscouts)));
            else
                return Op(false,"fails to build a scout :"+ops.error());


        }

        samples(std::vector<std::mt19937_64>&& mt,
                std::map<std::size_t, std::size_t>&& i_to_scouts,
                std::vector<mcmc_sample_t<Model,Adapative_Distribution_Generator>>&& scouts):
            mt_{std::move(mt)},i_to_scout{std::move(i_to_scouts)},scouts_{std::move(scouts)}
        {}

        typedef  samples self_type;
        typedef Parallel_Tempering enclosing_type;
        constexpr static auto  className=enclosing_type::className+my_static_string("_samples");
        static auto get_constructor_fields()
        {
            return std::make_tuple(
                        grammar::field(C<self_type>{},"mts",&self_type::get_mt),
                        grammar::field(C<self_type>{},"i_to_scout",&self_type::get_i_to_scout),
                        grammar::field(C<self_type>{},"scouts",&self_type::get_scouts)
                        );
        }
        std::vector<std::mt19937_64>const & get_mt()const {return mt_;}
        std::map<std::size_t, std::size_t>const & get_i_to_scout()const {return i_to_scout;}
        std::vector<mcmc_sample_t<Model,Adapative_Distribution_Generator>>const & get_scouts()const {return scouts_;}


    private:
        std::vector<std::mt19937_64> mt_;
        std::map<std::size_t, std::size_t> i_to_scout;
        std::vector<mcmc_sample_t<Model,Adapative_Distribution_Generator>> scouts_;

        static std::vector<std::mt19937_64> mts(std::mt19937_64& mt,std::size_t n)
        {
            std::uniform_int_distribution<typename std::mt19937_64::result_type> useed;
            std::vector<std::mt19937_64> out(n);
            for (std::size_t i=0; i<n; ++i)
                out[i].seed(useed(mt));
            return out;
        }

        static std::map<std::size_t, std::size_t> init_map(std::size_t n)
        {
            std::map<std::size_t, std::size_t> out;
            for (std::size_t i=0; i<n; ++i) out[i]=i;
            return out;
        }


    };

    template<class Model,class samples>
    static double Evidence(const Model& m,const samples& s)
    {
        auto n=m.size();
        double beta0=m.beta(0);
        double sum=0;
        double sumdb=0;
        double logLik0=s.Scout(0).likelihood().logL();
        for (std::size_t i=1; i<m.size(); ++i)
        {
            double beta=m.beta(i);
            double db=beta0-beta;
            double logLik=s.Scout(i).likelihood().logL();
            sum+=db*(logLik0+logLik)/2;
            sumdb+=db;
            if ((i==n-1)&&(beta>0))
            {
                double db0=beta;
                sum+=(db0*(1+db0/db/2))*logLik-sqr(db0)/db/2*logLik0;
                sumdb+=db0;
            }
            beta0=beta;
            logLik0=logLik;
        }
        return sum;
    }


    template<class priorModel, class LikelihoodModel,class Adapative_Distribution_Generator>
    myOptional<samples<priorModel,LikelihoodModel,Adapative_Distribution_Generator>,reg_tag>
    start(const Thermodynamic_Model_Series<priorModel,LikelihoodModel>& model,std::mt19937_64& mt,std::vector<Adapative_Distribution_Generator>& g)const
    {
        return  samples<priorModel,LikelihoodModel,Adapative_Distribution_Generator>::evaluate(get_Metropolis(),model,mt,g);


    }




    template<class priorModel, class LikelihoodModel,class Adapative_Distribution_Generator>
    Op_void next(const Thermodynamic_Model_Series<priorModel,LikelihoodModel>& model,std::mt19937_64& mt, std::vector<Adapative_Distribution_Generator>& adaptives,samples<priorModel,LikelihoodModel,Adapative_Distribution_Generator>& current) const
    {
        std::vector<Op_void> res(model.size());
        std::vector<std::stringstream> ss(model.size());
#pragma omp parallel for
        for (std::size_t i=0; i<model.size(); ++i)
        {
            res[i]=get_Metropolis().next(model.model(i),current.mt(i),adaptives[i],current.Scout(i),ss[i]);
        }

        for (auto& s: ss) std::cerr<<s.str();
        auto ops=consolidate(std::move(res));
        if (!ops.has_value())
            return Op_void(false,"get_Metropolis fails: "+ops.error());
        std::uniform_real_distribution<> U;
        double r=U(mt);

        std::cerr<<" test thermo r="<<r<< "P_jump_="<<P_jump_;
        if (r<P_jump_)
        {
            std::size_t j;
            if (U(mt)<0.5)
                j=0;
            else
                j=1;

            std::cerr<<" j="<<j<<"\n";
            //#pragma omp parallel for
            for (std::size_t i=j; i<model.size()-1; i+=2)
            {
                double A=current.Accept(i,model);
                double r=std::uniform_real_distribution<>(0,1)(current.mt(i));
                if (r<A)
                    current.swap(i,model);
            }
        }
        return Op_void(true,"");
    }


    Parallel_Tempering(Adaptive_Metropolis_H<Hastings>amh, double P_jump):amh_{std::move(amh)},P_jump_{P_jump}{}

    Adaptive_Metropolis_H<Hastings> const & get_Metropolis()const {return amh_;}

private:
    Adaptive_Metropolis_H<Hastings> amh_;
    double P_jump_;

};


class Ensemble_Parallel_Tempering
{

public:


    template<class Model,class Adaptive_Mover>
    class samples
    {
    public:
        emcee_sample<Model,Adaptive_Mover>& Scout(std::size_t i){ return scouts_[i];}
        emcee_sample<Model,Adaptive_Mover> const & Scout(std::size_t i)const {return scouts_[i];}

        void swap(std::size_t i, std::size_t j, std::size_t k,const Model& model)
        {
            std::swap(Scout(i).Walker(j),Scout(i+1).Walker(k));
            std::swap(ij_history_[i][j], ij_history_[i+1][k]);
            Scout(i).Walker(j).set_beta(model.beta(i));
            Scout(i+1).Walker(k).set_beta(model.beta(i+1));

        }

        double Accept(std::size_t i, std::size_t j, std::size_t k,const Model& model)
        {
            double logA=(model.beta(i)-model.beta(i+1))*(Scout(i).Walker(j).likelihood().logL()-Scout(i+1).Walker(k).likelihood().logL());
            return  std::min(1.0,exp(logA));

        }

        std::mt19937_64& mt(std::size_t i){return mt_[i];}

        samples(const Ensemble_Metropolis_Hastings& emh,const Model& m,std::mt19937_64& _mt,std::vector<Adaptive_Mover>& adaptives):
            mt_{mts(_mt,m.size())}, ij_history_{init_history(m.size(),emh.numWalkers())}, scouts_(m.size()){
            for (std::size_t i=0; i<m.size(); ++i)
            {
                scouts_[i]=emh.start(m.model(i),mt(i),adaptives[i]);
            }

        }


    private:
        std::vector<std::mt19937_64>& mt_;
        std::vector<std::vector<std::pair<std::size_t, std::size_t>>> ij_history_;
        std::vector<emcee_sample<Model,Adaptive_Mover>> scouts_;

        std::vector<std::mt19937_64> mts(std::mt19937_64& mt,std::size_t n)
        {
            std::uniform_int_distribution<typename std::mt19937_64::result_type> useed;
            std::vector<std::mt19937_64> out(n);
            for (std::size_t i=0; i<n; ++i)
                out[i].seed(useed(mt));
            return out;
        }

        static std::vector<std::vector<std::pair<std::size_t, std::size_t>>> init_history(std::size_t nScouts,std::size_t nWalkers)
        {
            std::vector<std::vector<std::pair<std::size_t, std::size_t>>> out(nScouts,std::vector<std::pair<std::size_t, std::size_t>>(nWalkers));
            for (std::size_t i=0; i<nScouts; ++i)
                for(std::size_t j=0; j<nWalkers; ++j)
                    out[i][j]=std::pair(i,j);
            return out;

        }





    };



    template<class Model, class AdaptiveMover>
    myOptional_t<samples<Model,AdaptiveMover>> start(const Model& model,std::mt19937_64& mt, std::vector<AdaptiveMover> adaptives)
    {
        return  samples<Model,AdaptiveMover>(get_Ensamble_Metropolis(),model,mt,adaptives);
    }




    template<class Model, class AdaptiveMover>
    void next(const Model& model,std::mt19937_64& mt, std::vector<AdaptiveMover>& adaptives,samples<Model,AdaptiveMover>& current)
    {
#pragma omp parallel for
        for (std::size_t i=0; i<model.size(); ++i)
        {
            get_Ensamble_Metropolis().next(model.model(i),current.mt(i),adaptives(i),current.Scout(i));
        }

        std::uniform_real_distribution<> U;
        double r=U(mt);

        if (r<P_jump_)
        {
            std::size_t j;
            if (U(mt)<0.5) j=0; else j=1;

#pragma omp parallel for
            for (std::size_t i=j; i<model.size()-1; i+=2)
            {
                std::size_t j=std::uniform_int_distribution<std::size_t>(0,current.Scout(i).numWalkers()-1)(current.mt(i));
                std::size_t k=std::uniform_int_distribution<std::size_t>(0,current.Scout(i+1).numWalkers()-1)(current.mt(i));
                double A=current.Accept(i,j,k,model);
                double r=std::uniform_real_distribution<>(0,1)(current.mt(i));
                if (r<A)
                    current.swap(i,j,k,model);
            }
        }
    }

    Ensemble_Metropolis_Hastings const & get_Ensamble_Metropolis(){return emh_;}
    Ensemble_Parallel_Tempering(Ensemble_Metropolis_Hastings ensamble_mh, double probability_jump):emh_{ensamble_mh},P_jump_{probability_jump}{}

private:

    Ensemble_Metropolis_Hastings emh_;
    double P_jump_;



};


struct OutputGenerator
{


    template<class priorModel, class LikelihoodModel>
    using ALM=
    std::vector<Adaptive_Parameterized_Distribution_Generator<LevenbergMarquardt<Thermodynamic_Model<priorModel,LikelihoodModel>>>>;

    template<template <class, class, class>class samples, class priorModel, class LikelihoodModel,class A>
    void print_title(const samples<priorModel,LikelihoodModel,A>&  )
    {
        os<<"nsample"<<sep<<"beta"<<sep<<"Field"<<sep<<"par"<<sep<<"value"<<end_of_line{};
    };


    template<class Model, class Samples>
    void print_sample(std::size_t isample,const Model& m,const Samples& state )
    {
        os<<std::setprecision(std::numeric_limits<double>::digits10 + 1);

        os<<isample<<sep<<""<<sep<<"Evidence"<<sep<<""<<sep<<Parallel_Tempering<true>::Evidence(m,state)<<end_of_line{};
        for (std::size_t i=0; i<state.size(); ++i)
        {
            os<<isample<<sep<<m.beta(i)<<sep<<"Prob"<<sep<<""<<sep<<state.Scout(i).logL()<<end_of_line{};
            os<<isample<<sep<<m.beta(i)<<sep<<"Prior"<<sep<<""<<sep<<state.Scout(i).prior().logL()<<end_of_line{};
            os<<isample<<sep<<m.beta(i)<<sep<<"LogLik"<<sep<<""<<sep<<state.Scout(i).likelihood().logL()<<end_of_line{};
            if (parameter_)
                for (std::size_t k=0; k<state.Scout(i).x().size(); ++k)
                    os<<isample<<sep<<m.beta(i)<<sep<<"X"<<sep<<k<<sep<<state.Scout(i).x()[k]<<end_of_line{};
            if (gradient_)
            {
                for (std::size_t k=0; k<state.Scout(i).x().size(); ++k)
                    os<<isample<<sep<<m.beta(i)<<sep<<"bG"<<sep<<k<<sep<<state.Scout(i).G()[k]<<end_of_line{};
                for (std::size_t k=0; k<state.Scout(i).x().size(); ++k)
                    os<<isample<<sep<<m.beta(i)<<sep<<"PG"<<sep<<k<<sep<<state.Scout(i).prior().G()[k]<<end_of_line{};
                for (std::size_t k=0; k<state.Scout(i).x().size(); ++k)
                    os<<isample<<sep<<m.beta(i)<<sep<<"LG"<<sep<<k<<sep<<state.Scout(i).likelihood().G()[k]<<end_of_line{};
            }
        }
        os.flush();

    };


    OutputGenerator(std::ostream& ost, bool parameter,bool gradient, std::string separator="\t"):os(ost),parameter_{parameter}, gradient_{gradient},sep{separator}{}
private:
    std::ostream& os;
    bool parameter_;
    bool gradient_;
    std::string sep;
};




template<class Metropolis_algorithm,class Model, class Distribution_Generator,class OutputGenerator>
Op_void run_Montecarlo_Markov_Chain(const Metropolis_algorithm& mcmc,std::mt19937_64& mt,const Model& M,  Distribution_Generator& G,std::size_t nsamples, OutputGenerator& O)
{
    auto state=mcmc.start(M,mt, G);
    if (!state)
    {
        std::cerr<<state.error();
        return Op_void(false," does not start at all "+state.error());
    }

    std::cerr<<"\nstate=mcmc.start(M,mt, G)=\n"<<state.value();
    O.print_title(state.value());
    //  O.print_title(G);
    std::size_t isample=0;

    O.print_sample(isample,M,state.value());

    for (std::size_t i=0; i<nsamples; ++i)
    {
        auto res=mcmc.next(M,mt,G,state.value());
        if (!res) return Op_void(false," interrupted at sample i="+ToString(i)+res.error());
        std::cerr<<"\nmcmc.next(M,mt,G,state)  state=\n"<<state.value();

        std::cerr<<"\n ------Adaptive GENERATOR----------\n"<<G<<"\n ------Adaptive GENERATOR-FIN---------\n";

        ++isample;
        O.print_sample(isample,M,state.value());
        //O.print_row(G);
    }
    return Op_void(true,"run for "+ToString(nsamples)+" samples");
}

template<class PriorModel, class Likelihood_Model, class OutputGenerator>
Op_void run_Thermo_Levenberg_ProbVel(const PriorModel& prior,const Likelihood_Model& lik, std::mt19937_64& mt, const std::vector<double>& betas,
                                     const std::vector<double>& landa, const std::vector<std::vector<double>>& landa_50_hill,double pjump,double gain_moment,std::size_t nSamples, OutputGenerator& output)
{
    typedef  Thermodynamic_Model<PriorModel,Likelihood_Model> Model;

    Thermodynamic_Model_Series<PriorModel,Likelihood_Model> model(prior,lik,betas);
    Parallel_Tempering<true> PT(Adaptive_Metropolis_H<true>(100),pjump);
    std::vector<LevenbergMarquardt<Model>> las(landa.size());
    for (std::size_t i=0; i<las.size(); ++i)
    {
        las[i]=LevenbergMarquardt<Model>(landa[i]);
    }
    std::vector<Adaptive_Parameterized_Distribution_Generator<LevenbergMarquardt<Model>>>
            ala
            (betas.size(),
             Adaptive_Parameterized_Distribution_Generator<LevenbergMarquardt<Model>>(las,landa_50_hill,gain_moment));

    std::cerr<<"\n  std::vector<Adaptive_Parameterized_Distribution_Generator<LevenbergMarquardt<Model>>>"
               "ala (betas.size(),    Adaptive_Parameterized_Distribution_Generator<LevenbergMarquardt<Model>>(las,landa_50_hill,gain_moment)); \n";

    //  std::cerr<<ala;

    return run_Montecarlo_Markov_Chain(PT,mt,model,ala,nSamples,output);

}














} // namespace evidence


#endif // MYEVIDENCE_H
