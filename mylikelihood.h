#ifndef MYLIKELIHOOD_H
#define MYLIKELIHOOD_H

#include "Matrix.h"
#include "myDistributions.h"
#include "myparameters.h"
#include <iomanip>

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
double calculate_Likelihood(const Base_Distribution<E>& p, const D& data)
{
    if (std::isfinite(data))
        return p.logP(data);
    else return 0;
}



template <template<class...>class V,template <class>class Distribution,class T, class D>
double calculate_Likelihood(const V<Distribution<T>>& P, const D& data)
{

    assert(P.size()==data.num_measurements());
    double logL=0;
    for (std::size_t i=0; i<data.num_measurements(); ++i)
    {
        double lik=calculate_Likelihood(P[i],data[i]);
        logL+=lik;
    }
    return logL;
}


template <class E, class D>
double calculate_Gradient(const Base_Distribution<E>& d,const Base_Distribution<E>& d0, const D& data)
{
    if (std::isfinite(data))
        return d.logP(data)-d0.logP(data);
    else return 0;
}


template <template <class...>class V,template <class...>class Distributions, typename T,class D>
double calculate_Gradient(const V<Distributions<T>>& d,const V<Distributions<T>>& d0, const D& data)
{
    assert(d.size()==d0.size());
    double out=0;
    for (std::size_t i=0; i<d.size(); ++i)
    {
        out+=calculate_Gradient(d[i],d0[i],data[i]);
    }
    return out;
}

template <template <class>class Distributions, typename T,class D>
double calculate_Gradient(const Distributions<T>& d,const Distributions<T>& d0, const D& data,double eps)
{
    return calculate_Gradient(d,d0,data)/eps;
}

template <template <class...>class Distributions, typename T,class D>
M_Matrix<double> calculate_Gradient(const Distributions<T>& d0,const std::vector<Distributions<T>>& d, const D& data,double eps)
{
    //typedef myOptional_t<M_Matrix<double>> Op;
    M_Matrix<double> out(1,d.size());
    for (std::size_t i=0; i<d.size();++i)
    {
        out[i]= calculate_Gradient(d[i],d0,data)/eps;
    }
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
auto calculate_Hessian(const std::vector<Distributions>& d0, const std::vector<std::vector<Distributions>>& d,const D& ,double eps, std::ostream& os)
{
    auto k=d.size();
    auto n=d0.size();
  //  assert(n==data.size());

    M_Matrix<double> out(k,k,M_Matrix<double>::SYMMETRIC,0.0);
    for (std::size_t i=0; i<n; ++i)
    {
        auto& d_i=d0[i];
        auto param=getParameter(d_i);
     //   os<<"\n -------------i="<<i<<"------------\n";
     //   os<<"\nparam "<<param;
        auto FIM=getFIM(d_i);
    //    os<<"\nFIM\n"<<FIM;
    //    os<<"eps="<<eps;
        auto npar=param.size();
        M_Matrix<double> J(npar,k);
        for (std::size_t j=0; j<k; ++j)
        {
            auto dpar=(getParameter(d.at(j).at(i))-param)*1.0/eps;
            for(std::size_t jj=0; jj<npar; ++jj)
                J(jj,j)=dpar[jj];
        }
      //  os<<"\nJ\n"<<J;

        /*auto H_=TranspMult(J,FIM)*J;*/
        auto H=quadraticForm_BT_A_B(FIM,J);
//        are_Equal<true,M_Matrix<double>>().test_sum(H,Transpose(J)*FIM*J, std::cerr);
     //   os<<"\nTranspose(J)*FIM*J-H\n"<<Transpose(J)*FIM*J-H;



        out-=H;
      //  os<<"\nout\n"<<out;

    }
   // os<<"\nout\n"<<out;
    return out;

}



class logLikelihood{
public:

    constexpr static auto const className=my_static_string("logLikelihood");
    //std::string myClass()const  { return className.str();}

    typedef   logLikelihood self_type ;
    static auto get_constructor_fields()
    {
        return std::make_tuple(
                    grammar::field(C<self_type>{},"logL",&self_type::logL)
                    );
    }

    double logL()const { return logL_;}
    logLikelihood(double value):logL_{value}{}
    logLikelihood():logL_{std::numeric_limits<double>::quiet_NaN()}{}
    operator bool()const { return std::isfinite(logL_);}
    void set_logL(double logLik){ logL_=logLik;}
private:
    double logL_;
};

class DlogLikelihood: public logLikelihood
{
public:
    typedef  logLikelihood base_type;

    constexpr static auto const className=my_static_string("DlogLikelihood_")+base_type::className;
    //std::string myClass()const  { return className.str();}

    typedef   DlogLikelihood self_type ;
    static auto get_constructor_fields()
    {
        double (logLikelihood::*myLogL) () const=&base_type::logL;
        M_Matrix<double> const & (self_type::*myG) () const=&self_type::G;
        const M_Matrix<double>& (self_type::*myH) () const=&self_type::H;

        return std::make_tuple(
                    grammar::field(C<self_type>{},"logL",myLogL),
                    grammar::field(C<self_type>{},"Gradient",myG),
                    grammar::field(C<self_type>{},"Hessian",myH)
                    );
    }

    const M_Matrix<double>& G()const {return G_;}
    const M_Matrix<double>& H()const {return H_;}
    operator bool()const { return (G_.size()>0)&&(H_.size()>0)&&logLikelihood::operator bool();}
    void  set_G(M_Matrix<double>&& gradient) { G_=std::move(gradient);}
    void  set_H(M_Matrix<double>&& hessian) {
        assert(hessian.isSymmetric());
        H_=std::move(hessian);}




    DlogLikelihood(double logL, const M_Matrix<double>& Gradient, const M_Matrix<double>& Hessian)
        :logLikelihood(logL), G_{Gradient},H_{Hessian}{assert(H_.isSymmetric());}
    DlogLikelihood(double logL,  M_Matrix<double>&& Gradient,  M_Matrix<double>&& Hessian)
        :logLikelihood(logL), G_{std::move(Gradient)},H_{std::move(Hessian)}{assert(H_.isSymmetric());}
    DlogLikelihood()=default;
private:
    M_Matrix<double> G_;
    M_Matrix<double> H_;
};






/// Aproximacion por Fisher Information Matrix al Hessiano
///
///

template<class Parameters_distribution>
class Prior_Model
{
public:
    typedef  Cs<Parameters_distribution> template_types;
    constexpr static auto const className=my_static_string("Prior_Model_")+my_trait<template_types>::className;

    typedef  typename Parameters_distribution::Parameters Parameters;
    auto sample(std::mt19937_64& mt)const { return prior_.sample(mt);}

    myOptional_t<logLikelihood> compute_Likelihood(const Parameters& x)const
    {
        typedef myOptional_t<logLikelihood> Op;
        double logL=prior_.logP(x);
        if (std::isfinite(logL))
            return Op(logLikelihood(logL));
        else return Op(false,"not a finite value="+ToString(logL));
    }
    myOptional_t<DlogLikelihood> compute_DLikelihood(const Parameters& x)const
    {
        typedef myOptional_t<DlogLikelihood> Op;
        double logL=prior_.logP(x);
        auto G=prior_.dlogL_dx(x);
        auto H=prior_.dlogL_dx2(x);
        std::stringstream ss;
        if (are_finite<true,double>().test(logL,ss)&&
                are_finite<true,M_Matrix<double>>().test(G,ss)&&
                are_finite<true,M_Matrix<double>>().test(H,ss)
                )
            return Op(DlogLikelihood(std::move(logL),std::move(G),std::move(H)));
        else
            return Op(false, ss.str() );
    }
    Prior_Model(const Parameters_distribution& prior):prior_{prior}{}

private:
    Parameters_distribution prior_;
};





template <class Distribution_Model>
class Likelihood_Model
{
public:
    typedef M_Matrix<double> Parameters;

    typedef  Likelihood_Model self_type;
    constexpr static auto const className=my_static_string("Likelihood_Model")+my_trait<Distribution_Model>::className;

    template <class Data>
    myOptional_t<logLikelihood> compute_Likelihood(std::ostream& os,const Parameters& p, const Data& d_) const
    {
        typedef myOptional_t<logLikelihood> Op;
        auto D0=l_.compute_Distribution(d_,p);
        if (D0.has_value())
        {
            double logL=calculate_Likelihood(D0,d_);
            std::stringstream ss;
            if (are_finite<true,double>().test(logL,ss))
                return Op(logLikelihood(logL));
            else
                return Op(false,ss.str());
        }
        else return Op(false,"fails to compute model data"+D0.error());

    }
    Likelihood_Model(const Distribution_Model& l):l_{l}{}

    const Distribution_Model& model() const {return l_;}


protected:
    Distribution_Model l_;
 };




template <class Distribution_Model>
class FIM_Model: public Likelihood_Model<Distribution_Model>
{
public:
    typedef  FIM_Model self_type;
    typedef  Likelihood_Model<Distribution_Model> base_type;
    constexpr static auto const className=my_static_string("FIM_Model")+my_trait<base_type>::className;



    typedef M_Matrix<double> Parameters;

    typedef Likelihood_Model<Distribution_Model> L;

    template <class Data>
    auto compute_DLikelihood(std::ostream& os,const Parameters& p, const Data& data)const
    {
        typedef myOptional_t<DlogLikelihood> Op;
        auto D0res=L::model().compute_Distribution_aux(data,p,os);
        if (!D0res.has_value())
            return Op(false,"fails to compute the model :"+D0res.error());
        else
        {

           auto [D0,aux]=std::move(D0res).value();

            std::vector<std::decay_t<decltype (D0)>> D(p.size());
            for (std::size_t i=0; i<p.size(); ++i)
            {
                Parameters x(p);
                x[i]=x[i]+eps_;
                auto OpD=L::model().compute_Distribution_aux(data,x,os,aux);
                if (OpD.has_value())
                    D[i]=std::move(OpD).value();
                else
                    return Op(false," getDLikelihood error at i="+ToString(i)+"th parameter  :"+OpD.error());
            }

            auto logL=calculate_Likelihood(D0,data);
            std::stringstream ss;
            if (!are_finite<true,double>().test(logL,ss)) return Op(false," getDLikelihood error for getLikelihood "+ ss.str());
            auto G=calculate_Gradient(D0,D,data, eps_);
            if (!are_finite<true,M_Matrix<double>>().test(G,ss)) return Op(false," getDLikelihood error for gradient "+ ss.str());

            auto H=calculate_Hessian(D0,D,data, eps_,os);
            if (!are_finite<true,M_Matrix<double>>().test(H,ss)) return Op(false," getDLikelihood error for Hessian "+ ss.str());
            return Op(DlogLikelihood(logL,std::move(G), std::move(H)));
        }
    }
    FIM_Model(const Distribution_Model& l,  double eps):
        Likelihood_Model<Distribution_Model> (l), eps_{eps}{}

private:
    double eps_;
};



class Likelihood_Test
{
public:
    template<class Simulation_Model, class Experiment, class Distribution_Model, class Parameters>
    auto compute_DLikelihood(std::ostream& os,const Simulation_Model& m, const Experiment& e,const Distribution_Model& lik,  const Parameters& p, double eps,std::mt19937_64& mt  )const
    {
        typedef myOptional_t<DlogLikelihood> Op;
        auto Data=m(e,p,mt);
        if (!Data.has_value())
            return Op(false,"Simulation failed "+Data.error());
        else
        {
            FIM_Model<Distribution_Model> fim(lik,eps);

            return fim.compute_DLikelihood(os,p,Data.value());

        }
    }


};


}










#endif // MYLIKELIHOOD_H
