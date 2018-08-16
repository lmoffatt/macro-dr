#ifndef MYEVIDENCE_H
#define MYEVIDENCE_H
#include "Matrix.h"
#include "myDistributions.h"
#include "myparameters.h"
namespace evidence {





/*
 *
 *  supongamos en el medio que tenemos un Modelo que tira
 *
 *
 *   Modelo + vector de Parametros  +  vector de Datos -> vector de distribucion de los parametros  y de los datos.
 *
 *
 *
 *
 * Como obtengo lo de arriba a partir de esto?
 * Posibilidad 1, a partir de derivada por diferencia
 *
 * Posibilidad 2 a partir de derivada explicita. (esta posibilidad es una pesadilla para macrodr no tiene sentido en esta etapa)
 *
 *
 *
 * */

template <class E, class D>
double getLikelihood(const Base_Distribution<E>& p, const D& data)
{
    return p.logP(data);
}




template <template<class...>class V,template <class>class Distribution,class T, class D>
double getLikelihood(const V<Distribution<T>>& P, const D& data)
{
    assert(P.size()==data.num_measurements());
    double logL=0;
    for (std::size_t i=0; i<data.num_measurements(); ++i)
    {
        logL+=getLikelihood(P[i],data[i]);
    }
    return logL;
}


template <class E, class D>
double getGradient(const Base_Distribution<E>& d,const Base_Distribution<E>& d0, const D& data)
{
    return d.logP(data)-d0.logP(data);
}


template <template <class...>class V,template <class...>class Distributions, typename T,class D>
double getGradient(const V<Distributions<T>>& d,const V<Distributions<T>>& d0, const D& data)
{
    assert(d.size()==d0.size());
    double out=0;
    for (std::size_t i=0; i<d.size(); ++i)
    {
        out+=getGradient(d[i],d0[i],data[i]);
    }
    return out;
}

template <template <class>class Distributions, typename T,class D>
double getGradient(const Distributions<T>& d,const Distributions<T>& d0, const D& data,double eps)
{
    return getGradient(d,d0,data)/eps;
}

template <template <class...>class Distributions, typename T,class D>
M_Matrix<double> getGradient(const Distributions<T>& d0,const std::vector<Distributions<T>>& d, const D& data,double eps)
{
    M_Matrix<double> out(1,d.size());
    for (std::size_t i=0; i<d.size();)
        out[i]= getGradient(d[i],d0,data)/eps;
    return out;
}

template<class E>
auto getParameter(const Base_Distribution<E>& d)
{
    return d.param();
}

template<class E>
auto getParameter(const std::unique_ptr<Base_Distribution<E>>& d)
{
    return d->param();
}

template<class E>
auto getParameter(const Base_Distribution<E>* d)
{
    return d->param();
}

template<class E>
auto  getFIM(const Base_Distribution<E>& d)
{
    return d.Fisher_Information();
}

template<class E>
auto getFIM(const std::unique_ptr<Base_Distribution<E>>& d)
{
    return d->Fisher_Information();
}

template<class E>
auto getFIM(const Base_Distribution<E>* d)
{
    return d->Fisher_Information();
}


template <class Distributions, class D>
M_Matrix<double> getHessian(const std::vector<Distributions>& d0, const std::vector<std::vector<Distributions>>& d,const D& data,double eps)
{
    auto k=d.size();
    auto n=d0.size();
    assert(n==data.size());

    M_Matrix<double> out(k,k,0.0);
    for (std::size_t i=0; i<n; ++i)
    {
        auto& d_i=d0[i];
        auto param=getParameter(d_i);
        auto FIM=getFIM(d_i);
        auto npar=param.size();
        M_Matrix<double> J(npar,k);
        for (std::size_t j=0; j<k; ++j)
        {
            J(j,":")=(getParameter(d[j][i])-param)*(1.0/eps);
        }
        auto H=TranspMult(J,FIM)*J;
        out+=H;
    }
    return out;
}



class logLikelihood{
public:
    double logL()const { return logL_;}
    logLikelihood(double value):logL_{value}{}
    logLikelihood():logL_{std::numeric_limits<double>::quiet_NaN()}{}
    operator bool()const { return std::isfinite(logL_);}
    void logL(double logLik){ logL_=logLik;}
private:
    double logL_;
};

class DlogLikelihood: public logLikelihood
{
public:
    const M_Matrix<double>& G()const {return G_;}
    const M_Matrix<double>& H()const {return H_;}
    operator bool()const { return (G_.size()>0)&&(H_.size()>0)&&logLikelihood::operator bool();}
    void  G(M_Matrix<double>&& gradient) { G_=std::move(gradient);}
    void  H(M_Matrix<double>&& hessian) { H_=std::move(hessian);}



    DlogLikelihood(double logL, const M_Matrix<double>& Gradient, const M_Matrix<double>& Hessian)
        :logLikelihood(logL), G_{Gradient},H_{Hessian}{}
    DlogLikelihood(double logL,  M_Matrix<double>&& Gradient,  M_Matrix<double>&& Hessian)
        :logLikelihood(logL), G_{std::move(Gradient)},H_{std::move(Hessian)}{}
    DlogLikelihood()=default;
private:
    M_Matrix<double> G_;
    M_Matrix<double> H_;
};


class ThlogLikelihood: public logLikelihood
{
public:
    logLikelihood prior()const { return prior_;}
    logLikelihood likelihood()const { return lik_;}
    ThlogLikelihood(logLikelihood prior, logLikelihood lik, double beta):
        logLikelihood(prior.logL()+beta*lik.logL()),prior_{prior},lik_{lik}{}
    ThlogLikelihood()=default;
    void set_beta(double beta){ logL(prior().logL()+beta*likelihood().logL());

                              }
private:
    logLikelihood prior_;
    logLikelihood lik_;

};

class ThDlogLikelihood: public DlogLikelihood
{
public:
    DlogLikelihood prior()const { return prior_;}
    DlogLikelihood likelihood()const { return lik_;}
    ThDlogLikelihood(const DlogLikelihood& prior,const  DlogLikelihood& lik, double beta):
        DlogLikelihood(prior.logL()+beta*lik.logL(),prior.G()+lik.G()*beta,prior.H()+lik.H()*beta),prior_{prior},lik_{lik}{}
    ThDlogLikelihood(DlogLikelihood&& prior,  DlogLikelihood&& lik, double beta):
        DlogLikelihood(prior.logL()+beta*lik.logL(),prior.G()+lik.G()*beta,prior.H()+lik.H()*beta),prior_{std::move(prior)},lik_{std::move(lik)}{}
    ThDlogLikelihood()=default;
    void set_beta(double beta){
        logL(prior().logL()+beta*likelihood().logL());
        G(prior().G()+likelihood().G()*beta);
        H(prior().H()+likelihood().H()*beta);
    }


private:
    DlogLikelihood prior_;
    DlogLikelihood lik_;

};





/// Aproximacion por Fisher Information Matrix al Hessiano
///
///

template<class Model>
class Prior_Model
{
public:
    typedef M_Matrix<double> Parameters;

    Parameters sample(std::mt19937_64& mt)const { return prior_.sample(mt);}
    logLikelihood getLikelihood(const M_Matrix<double>& x)
    {
        return logLikelihood(prior_.logP(x));
    }
    DlogLikelihood getDLikelihood(const M_Matrix<double>& x)const
    {
        return DlogLikelihood(prior_.logP(x), prior_.dlogL_dx(x), prior_.dlogL_dx2(x));
    }
    Prior_Model(const Parameters_distribution<Model>& prior): prior_{prior}{}
private:
    Parameters_distribution<Model> prior_;

};





template <class Distribution_Model, class Data>
class Likelihood_Model
{
public:
    typedef M_Matrix<double> Parameters;


    logLikelihood getLikelihood(const Parameters& p) const
    {
        auto D0=l_.getDistribution(d_,p);
        return getLikelihood(D0,d_);
    }
    Likelihood_Model(const Distribution_Model& l, const Data& d):l_{l},d_{d}{}

    const Distribution_Model& model() const {return l_;}
    const Data&  data() const {return d_;}


protected:
    const Distribution_Model& l_;
    const Data& d_;
};




template <class Distribution_Model, class Data>
class FIM_Model: public Likelihood_Model<Distribution_Model,Data>
{
public:
    typedef M_Matrix<double> Parameters;

    typedef Likelihood_Model<Distribution_Model,Data> L;

    DlogLikelihood getDLikelihood(const Parameters& p)const
    {
        auto D0=L::model().getDistribution(L::data(),p);
        std::vector<decltype (D0)> D(p.size());
        for (std::size_t i=0; i<p.size(); ++i)
        {
            Parameters x(p);
            x[i]+=eps_;
            D[i]=L::model().getDistribution(L::data(),x);
        }

        auto logL=getLikelihood(D0,L::d_);

        auto G=getGradient(D0,D,L::data(), eps_);

        auto H=getHessian(D0,D,L::data(), eps_);
        return DlogLikelihood(logL,std::move(G), std::move(H));
    }
    FIM_Model(const Distribution_Model& l, const Data& d, double eps):
        Likelihood_Model<Distribution_Model, Data> (l,d), eps_{eps}{}

private:
    double eps_;
};







template<class PriorModel, class LikelihoodModel>
class Thermodynamic_Model
{
    static_assert (std::is_same_v<typename PriorModel::Parameters,typename LikelihoodModel::Parameters > );


public:
    typedef typename PriorModel::Parameters Parameters;

    Parameters sample(std::mt19937_64& mt)const {return prior().sample(mt);}

    typedef ThDlogLikelihood DLikelihoodResult;
    typedef ThlogLikelihood LikelihoodResult;
    Thermodynamic_Model(const PriorModel& prior, const LikelihoodModel& lik, double beta):prior_{prior}, lik_{lik}, beta_{beta}{}
    Thermodynamic_Model(const Thermodynamic_Model& other):prior_{other.prior_}, lik_{other.lik_}, beta_{other.beta_}{};
    Thermodynamic_Model(Thermodynamic_Model&& other):prior_{other.prior_}, lik_{other.lik_}, beta_{std::move(other.beta_)}{};
    Thermodynamic_Model& operator=(const Thermodynamic_Model&)=default;
    Thermodynamic_Model& operator=(Thermodynamic_Model&&other){
        if (&other!=this)
        {
            Thermodynamic_Model tmp(other);
            std::swap(*this,tmp);
            return *this;
        };
    }

    double beta()const {return beta_;}

    void set_beta(double _beta){beta_=_beta;}
    ThDlogLikelihood getDLikelihood(const Parameters& x) const
    {
        auto p=prior_.getDLikelihood(x);
        auto l=lik_.getDLikelihood(x);
        ThDlogLikelihood out(std::move(p),std::move(l), beta());
        return out;
    }

    ThlogLikelihood getLikelihood(const Parameters& x) const
    {
        auto p=prior_.getLikelihood(x);
        auto l=lik_.getLikelihood(x);
        ThlogLikelihood out(std::move(p),std::move(l), beta());
        return out;
    }

    PriorModel const & prior()const { return prior_;};
    LikelihoodModel const & likelihood() const { return lik_;}

private:
    PriorModel const & prior_;
    LikelihoodModel const & lik_;
    double beta_;
};


template<class PriorModel, class LikelihoodModel>
class Thermodynamic_Model_Series
{
    static_assert (std::is_same_v<typename PriorModel::Parameters,typename LikelihoodModel::Parameters > );


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

    typedef Normal_Distribution<Parameters> Distribution;
    typedef typename Model::DLikelihoodResult LikelihoodResult;


    struct myAcceptProb
    {
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
        double operator()(const LevenbergMarquardt<Model>& LM)const
        {
            return (1.0/(1.0+LM.landa()));
        }
    };

    typedef myExpectVelocity ExpectedVelocity;


    double landa()const { return landa_;}



    Distribution getDistribution(const Model& ,const Parameters x, const LikelihoodResult& logL ) const
    {
        auto& H=logL.H();
        auto& G=logL.G();
        auto Hl=H+diag(H)*landa_;
        auto cov=inv(Hl);
        if (cov.second.empty())
        {
            Parameters d=-(G*cov.first);
            Parameters newpoint=x+d;
            return Normal_Distribution<Parameters>(newpoint,cov.first);
        }
        else return {};
    }

    LikelihoodResult getLikelihood(const Model& m, const Parameters& x) const
    {
        return m.getDLikelihood(x);
    }

    LevenbergMarquardt(double _landa):landa_{_landa}{}
    LevenbergMarquardt()=default;

    bool operator<(const LevenbergMarquardt& other )const { return landa()<other.landa();}

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

    Distribution getDistribution(const Model& ,const Parameters& x, const LikelihoodResult&  ) const
    {
        return Normal_Distribution<Parameters>(x,cov_);
    }

    LikelihoodResult getLikelihood(const Model& m, const Parameters& x)
    {
        return m.getLikelihood(x);
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
static mcmc_sample_t<Model,Distribution_Generator> calculate(const Model& M, const Distribution_Generator& G, Parameters&& x)
{
    auto logL=G.getLikelihood(M,x);
    auto g=G.getDistribution(M,x,logL);
    return mcmc_sample_t<Model,Distribution_Generator>(std::move(x),std::move(logL),std::move(g));
}



template<bool Hastings=true>
class Metropolis_H
{
public:
    template<class mcmcsample>
    static double Acceptance(const mcmcsample& current, const mcmcsample& candidate)
    {
        double logL_current=double(current.logL());
        double logL_candidate=double(candidate.logL());
        if constexpr (Hastings)
        {
            double logP_forward=current.g().logP(candidate.x());
            double logP_backward=candidate.g().logP(current.x());
            double logA=(logL_candidate+logP_forward)-(logL_current+logP_backward);
            double A=std::min(1.0,exp(logA));
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
    mcmc_sample_t<Model,Distribution_Generator>
    start(const Model& M,std::mt19937_64& mt, const Distribution_Generator& G) const
    {
        mcmc_sample_t<Model,Distribution_Generator> out;
        std::size_t n=0;
        while(!out &&n<max_number_of_trials())
        {
            Parameters x=M.sample(mt);
            out=calculate(M,G,x);
            ++n;
        }
        return out;
    }

    template<class Model, class Distribution_Generator, class mySample=mcmc_sample_t<Model,Distribution_Generator>>
    bool next(const Model& M,std::mt19937_64& mt,const Distribution_Generator& G,mySample& current)const
    {
        auto candidate_point=current.g().sample(mt);
        mySample candidate=calculate(M,G,candidate_point);
        double A=Acceptance(current,candidate);
        std::uniform_real_distribution<double> u(0,1);
        double r=u(mt);
        bool accept=r<A;
        if (accept)
        {
            current=std::move(candidate);
            current.push_accepted();
        }
        else current.push_rejected();
        return bool(current);

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
    template<class Model, class Adapative_Distribution_Generator,
             class mySample=mcmc_sample_t<Model,Adapative_Distribution_Generator>>
    mySample start(const Model& M,std::mt19937_64& mt, Adapative_Distribution_Generator& adaptive)const
    {
        auto G=adaptive.sample(mt);
        auto   out=get_Metropolis().start(M,mt,G);
        return out;
    }




    template<class Model, class Adapative_Distribution_Generator,
             class mySample=mcmc_sample_t<Model,Adapative_Distribution_Generator>>
    bool next(const Model& M,std::mt19937_64& mt, Adapative_Distribution_Generator& adaptive, mySample& current)const
    {
        auto G=adaptive.sample(mt);
        current.set_g(G.getDistribution(M,current.x(), current));
        get_Metropolis().next(M,mt,G,current);
        adaptive.push_outcome(G,current, current.accept());
        return bool(current);
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
            auto logL=Mover::getLikelihood(model,x);
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
        LikelihoodResult getLikelihood(const Model& m, Parameters& x){ return m.getLikelihood(x);}

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
            LikelihoodResult lik=getLikelihood(model,candidate);
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


    Parameterized_Distribution sample(std::mt19937_64& mt)const{return p_.sample(mt);}


    template<class mcmc>
    void push_outcome(const Parameterized_Distribution& landa, const mcmc&, bool accept){ if (accept) push_acceptance(landa); else push_rejection(landa);}

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
    DistGen sample(std::mt19937_64& mt)const
    { return sample_rev_map(rev_,mt);}

    void push_outcome(const DistGen& landa, const M_Matrix<double>& ,bool accept)
    { if (accept) push_acceptance(landa);
        else push_rejection(landa);


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

    template<class priorModel, class LikelihoodModel,class Adapative_Distribution_Generator>
    class samples
    {
    public:
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
            double logA=(model.beta(i)-model.beta(i+1))*(Scout(i).likelihood().logL()-Scout(i+1).likelihood().logL());
            return  std::min(1.0,exp(logA));

        }

        std::mt19937_64& mt(std::size_t i){return mt_[i];}

        samples(const Adaptive_Metropolis_H<Hastings>& adm,Thermodynamic_Model_Series<priorModel,LikelihoodModel> m,std::mt19937_64& _mt, std::vector<Adapative_Distribution_Generator>& adaptive):
            mt_{mts(_mt,m.size())}, i_to_scout{init_map(m.size())}, scouts_(m.size())
        {
            for (std::size_t i=0; i<m.size(); ++i)
            {
                scouts_[i]=adm.start(m.model(i),mt(i),adaptive[i]);
            }

        }

        operator bool()const { for (auto&e: scouts_) if (!e) return  false; return true;}
    private:
        std::vector<std::mt19937_64> mt_;
        std::map<std::size_t, std::size_t> i_to_scout;
        std::vector<mcmc_sample_t<Model,Adapative_Distribution_Generator>> scouts_;

        std::vector<std::mt19937_64> mts(std::mt19937_64& mt,std::size_t n)
        {
            std::uniform_int_distribution<typename std::mt19937_64::result_type> useed;
            std::vector<std::mt19937_64> out(n);
            for (std::size_t i=0; i<n; ++i)
                out[i].seed(useed(mt));
            return out;
        }

        std::map<std::size_t, std::size_t> init_map(std::size_t n)
        {
            std::map<std::size_t, std::size_t> out;
            for (std::size_t i=0; i<n; ++i) out[i]=i;
            return out;
        }


    };



    template<class priorModel, class LikelihoodModel,class Adapative_Distribution_Generator>
    samples<priorModel,LikelihoodModel,Adapative_Distribution_Generator>
    start(const Thermodynamic_Model_Series<priorModel,LikelihoodModel>& model,std::mt19937_64& mt,std::vector<Adapative_Distribution_Generator>& g)const
    {
        return  samples<priorModel,LikelihoodModel,Adapative_Distribution_Generator>(get_Metropolis(),model,mt,g);


    }




    template<class priorModel, class LikelihoodModel,class Adapative_Distribution_Generator>
    void next(const Thermodynamic_Model_Series<priorModel,LikelihoodModel>& model,std::mt19937_64& mt, std::vector<Adapative_Distribution_Generator>& adaptives,samples<priorModel,LikelihoodModel,Adapative_Distribution_Generator>& current) const
    {
#pragma omp parallel for
        for (std::size_t i=0; i<model.size(); ++i)
        {
            get_Metropolis().next(model.model(i),current.mt(i),adaptives[i],current.Scout(i));
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
                double A=current.Accept(i,model);
                double r=std::uniform_real_distribution<>(0,1)(current.mt(i));
                if (r<A)
                    current.swap(i,model);
            }
        }
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
    samples<Model,AdaptiveMover> start(const Model& model,std::mt19937_64& mt, std::vector<AdaptiveMover> adaptives)
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



template<class Metropolis_algorithm,class Model, class Distribution_Generator>
bool run_Montecarlo_Markov_Chain(const Metropolis_algorithm& mcmc,std::mt19937_64& mt,const Model& M,  Distribution_Generator& G,std::size_t nsamples)
{
    auto state=mcmc.start(M,mt, G);
    std::size_t isample=0;
    for (std::size_t i=0; i<nsamples; ++i)
    {
        mcmc.next(M,mt,G,state);
        if (!bool(state)) return false;
        else
            ++isample;
    }
    return true;
}

template<class PriorModel, class Likelihood_Model>
bool run_Thermo_Levenberg_ProbVel(const PriorModel& prior,const Likelihood_Model& lik, std::mt19937_64& mt, std::vector<double> betas,
                                  const std::vector<double>& landa, const std::vector<std::vector<double>>& landa_50_hill,double gain_moment,std::size_t nSamples)
{
    typedef  Thermodynamic_Model<PriorModel,Likelihood_Model> Model;

    Thermodynamic_Model_Series<PriorModel,Likelihood_Model> model(prior,lik,betas);
    Parallel_Tempering<true> PT(Adaptive_Metropolis_H<true>(100),0.5);
    std::vector<LevenbergMarquardt<Model>> las(landa.size());
    for (std::size_t i=0; i<las.size(); ++i)
    {
        las[i]=LevenbergMarquardt<Model>(landa[i]);
    }
    std::vector<Adaptive_Parameterized_Distribution_Generator<LevenbergMarquardt<Model>>>
            ala
            (betas.size(),
             Adaptive_Parameterized_Distribution_Generator<LevenbergMarquardt<Model>>(las,landa_50_hill,gain_moment));
    return run_Montecarlo_Markov_Chain(PT,mt,model,ala,nSamples);

}













} // namespace evidence


#endif // MYEVIDENCE_H
