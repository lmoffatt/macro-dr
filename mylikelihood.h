#ifndef MYLIKELIHOOD_H
#define MYLIKELIHOOD_H
#include "Matrix.h"
#include "myDistributions.h"
#include "myparameters.h"
#include "mydataframe.h"
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
std::tuple<double,double,double> calculate_Likelihood(const Base_Distribution<E>& p, const D& data)
{
    if (std::isfinite(data))
        return {p.logP(data),p.expected_logP(),p.variance_logP()};
    else return {0.0,0.0,0.0};
}



template <template<class...>class V,template <class>class Distribution,class T, class D>
std::tuple<double,double,double> calculate_Likelihood(const V<Distribution<T>>& P, const D& data)
{

    assert(P.size()==data.num_measurements());
    double logL=0;
    double elogL=0;
    double vlogL=0;

    for (std::size_t i=0; i<data.num_measurements(); ++i)
    {
        auto [lik,elik,vlik]=calculate_Likelihood(P[i],data[i]);
                logL+=lik;
                elogL+=elik;
                vlogL+=vlik;
    }
                return {logL,elogL,vlogL};
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


    template<class E>
    auto  get_elogL(const Base_Distribution<E>& d)
    {
        return d.expected_logP();
    }

    template<class E>
    auto get_elogL(const std::unique_ptr<Base_Distribution<E>>& d)
    {
        return d->expected_logP();
    }

    template<class E>
    auto get_elogL(const Base_Distribution<E>* d)
    {
        return d->expected_logP();
    }

    template <class Distributions>
    auto calculate_Hessian(std::size_t i,const std::vector<Distributions>& d0, const std::vector<std::vector<Distributions>>& d ,double eps)
    {
        auto k=d.size();
        //auto n=d0.size();
        //  assert(n==data.size());

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
        auto H=quadraticForm_BT_A_B(FIM,J);
        return -H;

    }



    template <class Distributions, class D>
    auto calculate_Hessian(const std::vector<Distributions>& d0, const std::vector<std::vector<Distributions>>& d,const D& ,double eps, std::ostream& )
    {
        auto k=d.size();
        auto n=d0.size();

        M_Matrix<double> out(k,k,M_Matrix<double>::SYMMETRIC,0.0);
        for (std::size_t i=0; i<n; ++i)
        {
            out+=calculate_Hessian(i,d0,d,eps);

        }
        return out;

    }


    class logLikelihood
    {
    public:
        //std::string myClass()const  { return className.str();}

        constexpr static auto const className=my_static_string("logLikelihood");

        typedef   logLikelihood self_type ;
        static auto get_constructor_fields()
        {
            double (logLikelihood::*myLogL) () const=&self_type::logL;
            return std::make_tuple(
                        grammar::field(C<self_type>{},"logL",myLogL),
                        grammar::field(C<self_type>{},"elogL",&self_type::elogL),
                        grammar::field(C<self_type>{},"vlogL",&self_type::vlogL)
                        );
        }

        double logL()const { return logL_;}
        double elogL()const { return elogL_;}
        double vlogL()const { return vlogL_;}

        double chilogL()const {return (logL()-elogL())/std::sqrt(vlogL());}

        double chi2logL()const {return sqr(logL()-elogL())/vlogL();}


        logLikelihood(double value, double elogL , double vlogL):logL_{value}, elogL_{elogL},vlogL_{vlogL}{}
        logLikelihood(std::tuple<double,double,double> logL):logL_{std::get<0>(logL)}, elogL_{std::get<1>(logL)},vlogL_{std::get<2>(logL)}{}

        logLikelihood():logL_{std::numeric_limits<double>::quiet_NaN()}{}
        operator bool()const { return std::isfinite(logL_);}
        void set_logL(double logLik, double elogL, double vlogL){ logL_=logLik;elogL_=elogL; vlogL_=vlogL;}


        template<class DataFrame>
        static void insert_col(DataFrame& d, const std::string& pre)
        {
            d.insert_column(pre+"logL",C<double>{});
            d.insert_column(pre+"elogL",C<double>{});
            d.insert_column(pre+"vlogL",C<double>{});
        }
        auto data_row(std::size_t)const {return std::tuple(logL(),elogL(),vlogL());}
        std::size_t data_size()const { return 1;}



    private:
        double logL_;
        double elogL_;
        double vlogL_;
    };
}
template<>
class moments<evidence::logLikelihood>
{
public:
    //std::string myClass()const  { return className.str();}

    constexpr static auto const className=my_static_string("logLikelihood_moments");

    typedef   moments<evidence::logLikelihood> self_type ;
    static auto get_constructor_fields()
    {
        return std::make_tuple(
                    grammar::field(C<self_type>{},"logL",&self_type::logL),
                    grammar::field(C<self_type>{},"elogL",&self_type::elogL),
                    grammar::field(C<self_type>{},"vlogL",&self_type::vlogL)
                    );
    }

    moments<double> logL()const { return logL_;}
    moments<double> elogL()const { return elogL_;}
    moments<double> vlogL()const { return vlogL_;}
    moments<double> chilogL()const { return chilogL_;}

    moments(const std::vector<evidence::logLikelihood>& l):
        logL_{l,&evidence::logLikelihood::logL},
        elogL_{l,&evidence::logLikelihood::elogL},
        vlogL_{l,&evidence::logLikelihood::vlogL},
        chilogL_{l,&evidence::logLikelihood::chilogL}
    {}

    template<class logL>
    moments(const std::vector<logL>& l):
        logL_{l,&logL::logL},
        elogL_{l,&logL::elogL},
        vlogL_{l,&logL::vlogL},
        chilogL_{l,&logL::chilogL}
    {}

    template<class logL, class Op, typename...Ts>
    moments(const std::vector<logL>& l, const Op& u, Ts...t):
        logL_{l,[&u,&t...](const logL& d){return u(d,t...).logL();}},
        elogL_{l,[&u,&t...](const logL& d){return u(d,t...).elogL();}},
        vlogL_{l,[&u,&t...](const logL& d){return u(d,t...).vlogL();}},
        chilogL_{l,[&u,&t...](const logL& d){return u(d,t...).chilogL();}}
    {}



    moments(moments<double> value, moments<double> elogL , moments<double> vlogL, moments<double> chilogL):logL_{value}, elogL_{elogL},vlogL_{vlogL},chilogL_{chilogL}{}



    moments()=default;
    operator bool()const { return std::isfinite(logL_.mean());}
    void set_logL(moments<double> logLik, moments<double> elogL, moments<double> vlogL, moments<double> chilogL){ logL_=logLik;elogL_=elogL; vlogL_=vlogL;chilogL_=chilogL;}


    template<class DataFrame>
    static void insert_col(DataFrame& d, const std::string& pre)
    {
        moments<double>::insert_col(d,"logL");
        moments<double>::insert_col(d,"elogL");
        moments<double>::insert_col(d,"vlogL");
        moments<double>::insert_col(d,"chilogL");

    }
    auto data_row(std::size_t i)const {
        return std::tuple_cat(logL().data_row(i),elogL().data_row(i),vlogL().data_row(i),chilogL().data_row(i));}
    std::size_t data_size()const { return 1;}

private:
    moments<double> logL_;
    moments<double> elogL_;
    moments<double> vlogL_;
    moments<double> chilogL_;
};

namespace evidence {
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
                    grammar::field(C<self_type>{},"elogL",&self_type::elogL),
                    grammar::field(C<self_type>{},"vlogL",&self_type::vlogL),
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




    DlogLikelihood(double logL, double elogL, double vlogL, const M_Matrix<double>& Gradient, const M_Matrix<double>& Hessian)
        :logLikelihood(logL,elogL,vlogL), G_{Gradient},H_{Hessian}{assert(H_.isSymmetric());}
    DlogLikelihood(std::tuple<double,double,double> logL, M_Matrix<double>&& Gradient,  M_Matrix<double>&& Hessian)
        :DlogLikelihood(std::get<0>(logL),std::get<1>(logL),std::get<2>(logL),std::move(Gradient),std::move(Hessian)){}
    DlogLikelihood(double logL, double elogL, double vlogL, M_Matrix<double>&& Gradient,  M_Matrix<double>&& Hessian)
        :logLikelihood(logL,elogL,vlogL), G_{std::move(Gradient)},H_{std::move(Hessian)}{assert(H_.isSymmetric());}
    DlogLikelihood()=default;
    template<class DataFrame>
    static void insert_col(DataFrame& d, const std::string pre)
    {
        base_type::insert_col(d,pre);
        d.insert_column(pre+"Gradient",C<double>{});
        d.insert_column(pre+"Hessian",C<double>{});

    }
    auto data_row(std::size_t i, std::size_t j) const
    {
        return std::tuple_cat(base_type::data_row(i),std::tuple(G_[i],H_(i,j)));
    }
    template<class Parameters_Distribution>
    auto data_row(std::size_t k, const Parameters_Distribution& prior) const
    {
        std::size_t n=G_.size();
        auto [i,j]=M_Matrix<double>::pos_to_ij_Symmetric(k,n);
                return std::tuple_cat(base_type::data_row(k),std::tuple(prior.name(j),prior.name(i),G_[j],H_(j-i,i+j)));
    }
                std::size_t data_size() const
        { return G_.size();}
        std::size_t data_big_size()const
        { return H_.size();}

        private:
        M_Matrix<double> G_;
        M_Matrix<double> H_;
    };


}

template<>
class moments<evidence::DlogLikelihood>: public moments<evidence::logLikelihood>
{
public:
    typedef  moments<evidence::logLikelihood> base_type;

    constexpr static auto const className=my_static_string("moments_DlogLikelihood_")+base_type::className;
    //std::string myClass()const  { return className.str();}

    typedef   moments self_type ;
    static auto get_constructor_fields()
    {

        return  std::tuple_cat(
                    base_type::get_constructor_fields(),
                    std::tuple(
                        grammar::field(C<self_type>{},"Gradient",&self_type::G),
                        grammar::field(C<self_type>{},"Hessian",&self_type::H)
                        ));
    }

    const moments<M_Matrix<double>>& G()const {return G_;}
    const M_Matrix<double>& H()const {return H_;}
    operator bool()const { return (G().mean().size()>0)&&(H_.size()>0)&&base_type::operator bool();}
    void  set_G(moments<M_Matrix<double>>&& gradient) { G_=std::move(gradient);}
    void  set_H(M_Matrix<double>&& hessian) {
        assert(hessian.isSymmetric());
        H_=std::move(hessian);}




    moments(moments<double> logL, moments<double> elogL, moments<double> vlogL,
            moments<double> chilogL,
            const moments<M_Matrix<double>>& Gradient, const M_Matrix<double>& Hessian)
        :base_type(logL,elogL,vlogL,chilogL), G_{Gradient},H_{Hessian}{assert(H_.isSymmetric());}


    moments(const std::vector<evidence::DlogLikelihood>& l):
        base_type(l),
        G_{l,&evidence::DlogLikelihood::G},
        H_{op::mean(l,&evidence::DlogLikelihood::H)}{}




    template<class DlogL>
    moments(const std::vector<DlogL>& l):
        base_type(l),
        G_{l,&DlogL::G},
        H_{op::mean(l,&DlogL::H)}{}

    template<class logL, class Op, typename...Ts>
    moments(const std::vector<logL>& l, const Op& u, Ts...t):
        base_type(l,u,t...),
        G_{l,[&u,&t...](const logL& d){return u(d,t...).G();}},
        H_{op::mean(l,[&u,&t...](const logL& d){return u(d,t...).H();})}{}




    moments()=default;
    template<class DataFrame>
    static void insert_col(DataFrame& d, const std::string& pre)
    {
        base_type::insert_col(d,pre);
        d.insert_column(pre+"parameters", C<std::string>());
        d.insert_column(pre+"parameters_T", C<std::string>());

        moments<M_Matrix<double>>::insert_col(d,pre+"Gradient");
        d.insert_column(pre+"mean_Hessian",C<double>{});

    }
    template<class Parameters_Distribution>
    auto data_row(std::size_t k, const Parameters_Distribution& prior) const
    {
        std::size_t n=G().mean().size();
        auto [i,j]=M_Matrix<double>::pos_to_ij_Symmetric(k,n);
                return std::tuple_cat(
                base_type::data_row(k),
                std::tuple(
                prior.name(i),
                prior.name(j)),
                G().data_row(i,j),
                std::tuple(H_(j-i,i+j)))
                ;
    }
                std::size_t data_size() const
        { return G().mean().size();}
        std::size_t data_big_size()const
        { return H_.size();}

        private:
        moments<M_Matrix<double>> G_;
        M_Matrix<double> H_;
    };


    namespace evidence {

    class PartialDLogLikelihood: public DlogLikelihood
    {
    public:
        typedef  DlogLikelihood base_type;

        constexpr static auto const className=my_static_string("DlogLikelihood_")+base_type::className;
        //std::string myClass()const  { return className.str();}

        typedef   PartialDLogLikelihood self_type ;
        static auto get_constructor_fields()
        {
            double (logLikelihood::*myLogL) () const=&base_type::logL;
            M_Matrix<double> const & (base_type::*myG) () const=&self_type::G;
            const M_Matrix<double>& (base_type::*myH) () const=&self_type::H;

            return std::make_tuple(
                        grammar::field(C<self_type>{},"logL",myLogL),
                        grammar::field(C<self_type>{},"elogL",&self_type::elogL),
                        grammar::field(C<self_type>{},"vlogL",&self_type::vlogL),
                        grammar::field(C<self_type>{},"Gradient",myG),
                        grammar::field(C<self_type>{},"Hessian",myH),
                        grammar::field(C<self_type>{},"Partial_DlogL",&self_type::partial_DlogL)
                        );
        }

        template<class DataFrame>
        static void insert_col(DataFrame& d)
        {
            base_type::insert_col(d,"");
            DlogLikelihood::insert_col(d,"partial_");
        }
        auto data_row(std::size_t isample,std::size_t i, std::size_t j) const
        {
            return std::tuple_cat(base_type::data_row(i,j), partial_DlogL()[isample].data_row(i,j));
        }




        std::vector<DlogLikelihood> const& partial_DlogL()const {return  partial_;}

        PartialDLogLikelihood(DlogLikelihood&& dlogL, std::vector<DlogLikelihood>&& dist )
            :DlogLikelihood(std::move(dlogL)),partial_{std::move(dist)}{}

        PartialDLogLikelihood(double logL, double elogL, double vlogL,const M_Matrix<double>& Gradient, const M_Matrix<double>& Hessian,const std::vector<DlogLikelihood>& dist )
            :DlogLikelihood(logL,elogL,vlogL,Gradient,Hessian),partial_{dist}{}
        PartialDLogLikelihood(double logL, double elogL, double vlogL, M_Matrix<double>&& Gradient,  M_Matrix<double>&& Hessian, std::vector<DlogLikelihood>&& dist )
            :DlogLikelihood(logL,elogL,vlogL,std::move(Gradient),std::move(Hessian)),partial_{std::move(dist)}{}
        PartialDLogLikelihood()=default;
    private:
        std::vector<DlogLikelihood> partial_;

    };



    }

    template <>
    class moments<evidence::PartialDLogLikelihood>: public moments<evidence::DlogLikelihood>
    {
    public:
        typedef  moments<evidence::DlogLikelihood> base_type;

        constexpr static auto const className=my_static_string("moments_DlogLikelihood_")+base_type::className;
        //std::string myClass()const  { return className.str();}

        typedef   moments self_type ;
        static auto get_constructor_fields()
        {

            return std::tuple_cat(
                        base_type::get_constructor_fields(),
                        std::tuple(
                            grammar::field(C<self_type>{},"Partial_DlogL",&self_type::partial_DlogL)
                            ));
        }

        std::vector<moments<evidence::DlogLikelihood>> const& partial_DlogL()const {return  partial_;}


        moments(const std::vector<evidence::PartialDLogLikelihood>& data):
            base_type(data),
            partial_(calc_moments(data))
        {}

        moments()=default;
    private:
        std::vector<moments<evidence::DlogLikelihood>> partial_;

        std::vector<moments<evidence::DlogLikelihood>>
        calc_moments(const std::vector<evidence::PartialDLogLikelihood>& data)
        {
            auto n=data[0].partial_DlogL().size();
            std::vector<moments<evidence::DlogLikelihood>> out(n);
            for (std::size_t i=0; i<n;++i)
            {
                out[i]=moments<evidence::DlogLikelihood>
                        (data,[](const evidence::PartialDLogLikelihood& p,std::size_t ii)
                {return p.partial_DlogL()[ii];},i);
            }
            return out;

        }

    };



    namespace evidence {


    /// Aproximacion por Fisher Information Matrix al Hessiano
    ///
    ///





    template <class Distribution_Model, class Data>
    class Likelihood_Model
    {
    public:
        typedef M_Matrix<double> Parameters;

        typedef  Likelihood_Model self_type;
        typedef  Cs<Distribution_Model, Data> template_types;
        constexpr static auto const className=my_static_string("Likelihood_Model")+my_trait<template_types>::className;

        myOptional_t<logLikelihood> compute_Likelihood(const Parameters& p) const
        {
            typedef myOptional_t<logLikelihood> Op;
            auto D0=l_.compute_Distribution(d_,p);
            if (D0.has_value())
            {
                auto logL=calculate_Likelihood(D0,d_);
                std::stringstream ss;
                if (are_finite<true,double>().test(logL,ss))
                    return Op(logLikelihood(logL));
                else
                    return Op(false,ss.str());
            }
            else return Op(false,"fails to compute model data"+D0.error());

        }
        Likelihood_Model(const Distribution_Model& l, const Data& d):l_{l},d_{d}{}

        const Distribution_Model& model() const {return l_;}
        const Data&  data() const {return d_;}
        void set_Data(Data&& d){ d_=std::move(d);}

    protected:
        Distribution_Model l_;
        Data d_;
    };




    template <class Distribution_Model, class Data>
    class FIM_Model: public Likelihood_Model<Distribution_Model,Data>
    {
    public:
        typedef  FIM_Model self_type;
        typedef  Likelihood_Model<Distribution_Model,Data> base_type;
        constexpr static auto const className=my_static_string("FIM_Model")+my_trait<base_type>::className;



        typedef M_Matrix<double> Parameters;

        typedef Likelihood_Model<Distribution_Model,Data> L;

        auto compute_DLikelihood(const Parameters& p, std::ostream& os)const
        {
            return compute_DLikelihood(p,os,L::data());
        }


        auto compute_Distributions(const Parameters& p, std::ostream& os, const Data& data)const
        {
            std::string error;
            auto D0res=L::model().compute_Distribution_aux(data,p,os);
            typedef decltype (D0res.value().first) myDistr;
            typedef myOptional_t<std::pair<myDistr,std::vector<myDistr>>> Op;

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
                return Op(std::pair(D0,D));

            }
        }

        auto compute_DLikelihood(const Parameters& p, std::ostream& os, const Data& data)const

        {
            typedef myOptional_t<DlogLikelihood> Op;
            auto res=compute_Distributions(p,os,data);
            if (!res) return Op(false,res.error());
            else
            {

                auto [D0,D]=std::move(res).value();
                        return  compute_DLikelihood(D0,D,os,data);

            }
            }

                        template<class Distribution>
                        auto compute_DLikelihood(const Distribution& D0, const std::vector<Distribution>& D,std::ostream& os, const Data& data)const
                {
                    typedef myOptional_t<DlogLikelihood> Op;

                    auto logL=calculate_Likelihood(D0,data);
                    std::stringstream ss;
                    if (!are_finite<true,double>().test(std::get<0>(logL),ss)) return Op(false," getDLikelihood error for getLikelihood "+ ss.str());
                    auto G=calculate_Gradient(D0,D,data, eps_);
                    if (!are_finite<true,M_Matrix<double>>().test(G,ss)) return Op(false," getDLikelihood error for gradient "+ ss.str());

                    auto H=calculate_Hessian(D0,D,data, eps_,os);
                    if (!are_finite<true,M_Matrix<double>>().test(H,ss)) return Op(false," getDLikelihood error for Hessian "+ ss.str());
                    return Op(DlogLikelihood(logL,std::move(G), std::move(H)));
                }




                auto compute_PartialDLikelihood(const Parameters& p, std::ostream& os, const Data& data)const
                {
                    typedef myOptional_t<PartialDLogLikelihood> Op;
                    auto res=compute_Distributions(p,os,data);
                    if (!res) return Op(false,res.error());
                    else
                    {
                        auto [D0,D]=std::move(res).value();

                                auto Dlik=  compute_DLikelihood(D0,D,os,data);
                                if (!Dlik) return Op(false,Dlik.error());

                                std::vector<DlogLikelihood> partials(D0.size());

                                for (std::size_t n=0; n<D0.size(); ++n)
                        {

                            auto logL=calculate_Likelihood(D0[n],data[n]);
                            std::stringstream ss;
                            if (!are_finite<true,double>().test(std::get<0>(logL),ss)) return Op(false," getDLikelihood error for getLikelihood "+ ss.str());
                            M_Matrix<double> G(1,p.size());
                            for (std::size_t i=0; i<p.size(); ++i)
                            {
                                G[i]= calculate_Gradient(D[i][n],D0[n],data[n],eps_);
                            }
                            if (!are_finite<true,M_Matrix<double>>().test(G,ss)) return Op(false," getDLikelihood error for gradient "+ ss.str());

                            auto H=calculate_Hessian(n,D0,D,eps_);
                            if (!are_finite<true,M_Matrix<double>>().test(H,ss)) return Op(false," getDLikelihood error for Hessian "+ ss.str());
                            partials[n]=DlogLikelihood(logL,std::move(G), std::move(H));
                        }

                        return Op(PartialDLogLikelihood(std::move(Dlik).value(),std::move(partials)));

                    }
                }


                FIM_Model(const Distribution_Model& l, const Data& d, double eps):
                    Likelihood_Model<Distribution_Model, Data> (l,d), eps_{eps}{}

                void set_Data(Data&& d){ base_type::set_Data(std::move(d));}

                private:
                double eps_;
            };
            template<class Parameters_Distribution,class Parameters_Values,class ExperimentData, class logLikelihood>
            class Likelihood_Analisis
            {

            public:
                constexpr static auto const className=my_static_string("Likelihood_Analisis");

                typedef   Likelihood_Analisis self_type ;
                static auto get_constructor_fields()
                {
                    return std::make_tuple(
                                grammar::field(C<self_type>{},"prior",&self_type::get_Parameters_Distribution),
                                grammar::field(C<self_type>{},"simulations",&self_type::get_Parameters),
                                grammar::field(C<self_type>{},"simulations",&self_type::get_Simulations),
                                grammar::field(C<self_type>{},"likelihoods",&self_type::get_Likelihoods)
                                );
                }

                auto size()const {return s_.size();}

                double mean_logL()const
                {
                    double out=s_[0].logL();
                    for (std::size_t i=1; i<s_.size(); ++i)
                        out+=s_[i].logL();
                    return out/s_.size();
                }



                M_Matrix<double> mean_Gradient()const
                {
                    M_Matrix<double> out=s_[0].G();
                    for (std::size_t i=1; i<s_.size(); ++i)
                        out+=s_[i].G();
                    return out/s_.size();
                }
                std::vector<M_Matrix<double>> partial_mean_Gradient()const
                {
                    auto ns=s_[0].partial_DlogL().size();
                    std::vector<M_Matrix<double>> out(ns);
                    for (std::size_t n=0; n<ns; ++n)
                    {
                        M_Matrix<double> m=s_[0].partial_DlogL()[n].G();
                        for (std::size_t i=1; i<s_.size(); ++i)
                            m+=s_[i].partial_DlogL()[n].G();

                        out[n]=m/s_.size();
                    }
                    return out;
                }
                M_Matrix<double> mean_sqr_Gradient()const
                {
                    M_Matrix<double> out=quadraticForm_XTX(s_[0].G());
                    for (std::size_t i=1; i<s_.size(); ++i)
                        out+=quadraticForm_XTX(s_[i].G());
                    return out/s_.size();
                }
                std::vector<M_Matrix<double>> partial_sqr_Gradient()const
                {
                    auto ns=s_[0].partial_DlogL().size();
                    std::vector<M_Matrix<double>> out(ns);
                    for (std::size_t n=0; n<ns; ++n)
                    {
                        M_Matrix<double> m=quadraticForm_XTX(s_[0].partial_DlogL()[n].G());
                        for (std::size_t i=1; i<s_.size(); ++i)
                            m+=quadraticForm_XTX(s_[i].partial_DlogL()[n].G());
                        out[n]=m/s_.size();
                    }
                    return out;
                }

                M_Matrix<double> mean_Hessian()const
                {
                    M_Matrix<double> out=s_[0].H();
                    for (std::size_t i=1; i<s_.size(); ++i)
                        out+=s_[i].H();
                    return out/s_.size();
                }


                std::vector<M_Matrix<double>> partial_mean_Hessian()const
                {
                    auto ns=s_[0].partial_DlogL().size();
                    std::vector<M_Matrix<double>> out(ns);
                    for (std::size_t n=0; n<ns; ++n)
                    {
                        M_Matrix<double> m=s_[0].partial_DlogL()[n].H();
                        for (std::size_t i=1; i<s_.size(); ++i)
                            m+=s_[i].partial_DlogL()[n].H();

                        out[n]=m/s_.size();
                    }
                    return out;
                }


                auto& get_Parameters_Distribution()const{ return prior_;}
                auto& get_Parameters()const{ return p_;}


                std::vector<ExperimentData> const & get_Simulations()const {return  e_;}

                std::vector<logLikelihood>const & get_Likelihoods()const { return s_;}


                moments<logLikelihood> get_Likelihood_Moments()const
                { return moments<logLikelihood>(get_Likelihoods());}


                Likelihood_Analisis(Parameters_Distribution&& prior,
                                    std::vector<Parameters_Values>&& par,
                                    std::vector<ExperimentData>&& e,
                                    std::vector<logLikelihood>&& s):
                    prior_{std::move(prior)},
                    p_{std::move(par)},
                    e_{std::move(e)},
                    s_{std::move(s)}{}

                Likelihood_Analisis(const Parameters_Distribution& prior,
                                    const std::vector<Parameters_Values>& par,
                                    const std::vector<ExperimentData>& e,
                                    const std::vector<logLikelihood>& s):
                    prior_{prior},
                    p_{par},
                    e_{e},
                    s_{s}{}


                Likelihood_Analisis()=default;




                io::myDataFrame<double, std::size_t,std::string>
                        make_Data_Frame()const
                {
                    io::myDataFrame<double, std::size_t,std::string> d;
                    d.insert_column("statistic", C<std::string>());
                    d.insert_column("n_sample", C<std::size_t>());


                    ExperimentData::insert_col(d);

                    d.insert_column("ParName", C<std::string>());
                    d.insert_column("Par_trValue", C<double>());
                    d.insert_column("Par_Value", C<double>());
                    d.insert_column("ParName_T", C<std::string>());
                    d.insert_column("Par_trValue_T", C<double>());
                    d.insert_column("Par_Value_T", C<double>());

                    PartialDLogLikelihood::insert_col(d);
                    auto nsamples=s_.size();
                    for (std::size_t i=0; i<s_.size(); ++i)
                    {
                        auto data_sample=std::tuple(std::string("value"),i);
                        auto& e=e_[e_.size()==nsamples?i:0];
                        auto& p=p_[p_.size()==nsamples?i:0];
                        auto& l=s_[i];
                        assert(e.steps().size()==l.partial_DlogL().size());
                        auto nsteps=e.steps().size();
                        for (std::size_t i_step=0; i_step<nsteps; ++i_step)
                        {
                            auto data_exp=e.data_row(i_step);
                            auto npar=p.size();
                            assert(npar==l.G().size());
                            for(std::size_t i_par=0; i_par<npar; ++i_par)
                            {
                                for (std::size_t i_parT=0; i_parT<=i_par; ++i_parT)
                                {
                                    auto data_lik=l.data_row(i_step,i_par,i_parT);
                                    auto data_par=
                                            std::tuple(prior_.name(i_par),
                                                       p[i_par],
                                                       prior_.tr_to_Parameter(p[i_par],i_par),
                                                       prior_.name(i_parT),
                                                       p[i_parT],
                                                       prior_.tr_to_Parameter(p[i_parT],i_parT));
                                    auto data_t=std::tuple_cat(data_sample,data_exp,data_par,data_lik);
                                    d.push_back_t(data_t);
                                }
                            }
                        }
                    }
                    return d;

                }


            private:
                Parameters_Distribution prior_;
                std::vector<Parameters_Values> p_;
                std::vector<ExperimentData> e_;
                std::vector<PartialDLogLikelihood> s_;
            };


            class Likelihood_Test
            {
            public:


                static    myOptional_t<std::pair<double,std::size_t>>
                        chitest(const M_Matrix<double>&  mean,const M_Matrix<double>& Cov, std::size_t nsamples, std::ostream& os )
                {
                    typedef  myOptional_t<std::pair<double,std::size_t>> Op;
                    auto Cov_inv=inv(Cov);
                    if (!Cov_inv)
                    {
                        return Op(false, "Covariance matrix is not inversible");
                    }
                    else
                    {
                        double chi2=xTSigmaX(mean,Cov_inv.value())*nsamples;
                        os<<"\nchi2 ="<<chi2<<" df="<<mean.size();

                        return  Op({chi2,nsamples*mean.size()});

                    }
                }


                static   std::pair<std::vector<double>,std::size_t>
                        t_test(const M_Matrix<double>&  mean,const M_Matrix<double>& Cov, std::size_t nsamples, std::ostream& os )
                {

                    std::vector<double> out(mean.size());
                    for (std::size_t i=0; i<mean.size(); ++i)
                        if (mean[i]==0) out[i]=0;
                        else
                            out[i]=mean[i]/std::sqrt(Cov(i,i)/nsamples);
                    os<<"\nt test ="<<out<<" df="<<nsamples;
                    return {out,nsamples};
                }



                static    myOptional_t<std::pair<double,std::size_t>>
                        Cov_test_1(const M_Matrix<double>&  Cov,const M_Matrix<double>& CovExp, std::size_t n, std::ostream& os )
                {
                    typedef  myOptional_t<std::pair<double,std::size_t>> Op;
                    auto p=CovExp.ncols();
                    auto Cov_inv=inv(CovExp);
                    if (!Cov_inv.has_value()) return Op(false, "singular target covariance");
                    auto S=Cov*Cov_inv.value();
                    auto chol_S=chol(S,"lower");
                    double logdet=logDiagProduct(chol_S.value());


                    double TLR=n*(p*std::log(1)-logdet+Trace(S)-p);
                    os<<"\n Trace(S)="<<Trace(S)<<" logdet="<<logdet;
                    os<<" TLR="<<TLR;
                    TLR=(1.0-1.0/(6*n-1)*(2*p+1-2.0/(p+1)))*TLR;
                    os<<" TLRc="<<TLR<<" df="<< (p*(p+1))/2;
                    return {{TLR,(p*(p+1))/2}};

                }

                static  myOptional_t<std::pair<double,std::size_t>>
                        Cov_test_2(const M_Matrix<double>&  Cov,const M_Matrix<double>& CovExp, std::size_t  , std::ostream& os)
                {
                    typedef  myOptional_t<std::pair<double,std::size_t>> Op;
                    auto p=CovExp.ncols();
                    auto Cov_inv=inv(CovExp);
                    if (!Cov_inv.has_value()) return Op(false, "singular target covariance");
                    auto S=Cov*Cov_inv.value();
                    auto TJn= 1.0/p *Trace(sqr(S/(1.0/p*Trace(S))-Matrix_Generators::eye<double>(p)));
                    os<<"\nTJn= "<<TJn-1<<" df="<< (p*(p+1))/2;

                    return {{TJn,(p*(p+1))/2-1}};

                }



                static auto sample_Cov(std::mt19937_64& mt,const M_Matrix<double>& Cov,std::size_t nsamples)
                {
                    typedef myOptional_t<std::pair<M_Matrix<double>,M_Matrix<double>>> Op;
                    auto k=Cov.ncols();
                    M_Matrix<double> mean(1,k,0.0);
                    auto normal=Normal_Distribution<M_Matrix<double>>::make(mean,Cov);
                    if (!normal.has_value())
                    {
                        return Op(false,"\nsingular Covariance!!\n");
                    }
                    else
                    {
                        M_Matrix<double> S(k,k,M_Matrix<double>::SYMMETRIC,0.0);
                        M_Matrix<double> m(1,k,0.0);
                        for (std::size_t i=0; i<nsamples; ++i)
                        {
                            auto x=normal.value().sample(mt);
                            S+=quadraticForm_XTX(x);
                            m+=x;
                        }
                        return Op(std::pair(m/nsamples,S/nsamples));
                    }
                }


                template<class ExperimentSimulation>
                class sample
                {
                public:
                    constexpr static auto const className=my_static_string("Likelihood_Test_sample");

                    typedef   sample self_type ;
                    static auto get_constructor_fields()
                    {
                        return std::make_tuple(
                                    grammar::field(C<self_type>{},"simulations",&self_type::getSimulations),
                                    grammar::field(C<self_type>{},"likelihoods",&self_type::getLikelihoods)
                                    );
                    }

                    auto size()const {return s_.size();}
                    M_Matrix<double> mean_Gradient()const
                    {
                        M_Matrix<double> out=s_[0].G();
                        for (std::size_t i=1; i<s_.size(); ++i)
                            out+=s_[i].G();
                        return out/s_.size();
                    }
                    std::vector<M_Matrix<double>> partial_mean_Gradient()const
                    {
                        auto ns=s_[0].partial_DlogL().size();
                        std::vector<M_Matrix<double>> out(ns);
                        for (std::size_t n=0; n<ns; ++n)
                        {
                            M_Matrix<double> m=s_[0].partial_DlogL()[n].G();
                            for (std::size_t i=1; i<s_.size(); ++i)
                                m+=s_[i].partial_DlogL()[n].G();

                            out[n]=m/s_.size();
                        }
                        return out;
                    }
                    M_Matrix<double> mean_sqr_Gradient()const
                    {
                        M_Matrix<double> out=quadraticForm_XTX(s_[0].G());
                        for (std::size_t i=1; i<s_.size(); ++i)
                            out+=quadraticForm_XTX(s_[i].G());
                        return out/s_.size();
                    }
                    std::vector<M_Matrix<double>> partial_sqr_Gradient()const
                    {
                        auto ns=s_[0].partial_DlogL().size();
                        std::vector<M_Matrix<double>> out(ns);
                        for (std::size_t n=0; n<ns; ++n)
                        {
                            M_Matrix<double> m=quadraticForm_XTX(s_[0].partial_DlogL()[n].G());
                            for (std::size_t i=1; i<s_.size(); ++i)
                                m+=quadraticForm_XTX(s_[i].partial_DlogL()[n].G());
                            out[n]=m/s_.size();
                        }
                        return out;
                    }

                    M_Matrix<double> mean_Hessian()const
                    {
                        M_Matrix<double> out=s_[0].H();
                        for (std::size_t i=1; i<s_.size(); ++i)
                            out+=s_[i].H();
                        return out/s_.size();
                    }


                    std::vector<M_Matrix<double>> partial_mean_Hessian()const
                    {
                        auto ns=s_[0].partial_DlogL().size();
                        std::vector<M_Matrix<double>> out(ns);
                        for (std::size_t n=0; n<ns; ++n)
                        {
                            M_Matrix<double> m=s_[0].partial_DlogL()[n].H();
                            for (std::size_t i=1; i<s_.size(); ++i)
                                m+=s_[i].partial_DlogL()[n].H();

                            out[n]=m/s_.size();
                        }
                        return out;
                    }



                    std::vector<ExperimentSimulation> const & getSimulations()const {return  e_;}

                    std::vector<PartialDLogLikelihood>const & getLikelihoods()const { return s_;}


                    sample(std::vector<ExperimentSimulation>&& e,std::vector<PartialDLogLikelihood>&& s):
                        e_{std::move(e)},
                        s_{std::move(s)}{}

                    sample(const std::vector<ExperimentSimulation>& e,const std::vector<PartialDLogLikelihood>& s):
                        e_{e},
                        s_{s}{}


                    sample()=default;

                private:
                    std::vector<ExperimentSimulation> e_;
                    std::vector<PartialDLogLikelihood> s_;
                };


                template<class Simulation_Model,class FIM_Model,  class Data, class ParametersDistribution,class Parameters>
                static    auto compute_PartialDLikelihood
                        (std::ostream& os,const Simulation_Model& sim, const FIM_Model& fim, const Data& e, const ParametersDistribution& prior, const Parameters& p,std::mt19937_64& mt  )
                {
                    auto data=sim.compute_simulation(e,p,mt);
                    typedef std::decay_t<decltype (data.value())> sim_type;
                    typedef myOptional_t<std::pair<sim_type,PartialDLogLikelihood>> Op;
                    if (!data.has_value())
                        return Op(false,"Simulation failed "+data.error());
                    else
                    {
                        auto Par=prior.Parameter_to_tr(p);
                        auto Dlik=fim.compute_PartialDLikelihood(Par,os,data.value());
                        if (!Dlik)
                            return Op(false,"Likelihood failed!!: "+Dlik.error());

                        return Op(std::pair(std::move(data).value(),std::move(Dlik).value()));
                    }
                }


                template<class Simulation_Model,class FIM_Model,  class Data, class ParametersDistribution,class Parameters>
                static auto
                        get_sample
                        (std::ostream& os,const Simulation_Model& sim, const FIM_Model& fim, const Data& e,  const ParametersDistribution& prior,const Parameters& p,std::mt19937_64& mt  , std::size_t nsamples)
                {
                    typedef std::decay_t<decltype (sim.compute_simulation(e,p,mt).value())>  sim_type;
                    typedef myOptional_t<sample<sim_type>> Op;
                    std::stringstream ss;
                    std::vector<sim_type> simuls(nsamples);
                    std::vector<PartialDLogLikelihood> s(nsamples);

                    bool succed=true;
#pragma omp parallel for
                    for (std::size_t i=0; i<nsamples; ++i)
                    {
                        auto logL=compute_PartialDLikelihood(os,sim,fim,e,prior,p,mt);
                        if (!logL.has_value())
                            succed=false;
                        else
                        {
                            auto pair=std::move(logL).value();
                            s[i]=std::move(pair.second);
                            simuls[i]=std::move(pair.first);
                            //             os<<"\n in sample\n"<<s[i];
                        }
                    }
                    if (!succed)
                        return  Op(false,"failed in one of the samples");
                    os<<"\n successfully obtained "+ToString(nsamples)+" samples for Gradient_Expectancy_Test \n";
                    for (std::size_t i=0; i<nsamples; ++i)
                    {
                        os<<"\n"<<s[i].logL()<<"\t"<<sqr(s[i].logL()-s[i].elogL())/s[i].vlogL()<<"\t"<<s[i].G();
                    }

                    auto samples=sample(std::move(simuls),std::move(s));
                    return Op(samples);

                }


                static Op_void Gradient_Expectancy_Test
                        (const M_Matrix<double>& Gradient,
                         const M_Matrix<double>& GradientVariance,
                         std::size_t nsamples, double pvalue, std::ostream& os)

                {


                    auto t_stud=t_test(Gradient,GradientVariance,nsamples,os);
                    auto chi2_FIM=chitest(Gradient,GradientVariance,nsamples, os);
                    if (!chi2_FIM) return Op_void(false, "Gradient_Expectancy_Test could not be performed: "+chi2_FIM.error());
                    else
                    {
                        auto chivalue=1.0-chi2_cdf(chi2_FIM.value().first, chi2_FIM.value().second);
                        if (chivalue<pvalue)
                            return Op_void(false,"Gradient_Expectancy_Test indicates a low probability of being true: chi2="+
                                           std::to_string(chi2_FIM.value().first)+" : p("+
                                           std::to_string(chivalue)+")<"+ToString(pvalue)+" : "+chi2_FIM.error());
                        else
                            return Op_void(true,"Gradient_Expectancy_Test fall within the expected values: chi2="+
                                           std::to_string(chi2_FIM.value().first)+" p("+
                                           std::to_string(chivalue)+")>"+ToString(pvalue)+" : "+chi2_FIM.error());

                    }
                }

                static Op_void Gradient_Variance_Expectancy_Test
                        (const M_Matrix<double>& GradientVariance, const M_Matrix<double>& ExpectedGradientVariance, std::size_t nsamples, double pvalue, std::ostream& os)
                {


                    auto chi2=Cov_test_1(GradientVariance,ExpectedGradientVariance,nsamples,os);
                    auto chi22=Cov_test_2(GradientVariance,ExpectedGradientVariance,nsamples,os);
                    if (!chi2) return Op_void(false, "Gradient_Variance_Expectancy_Test could not be performed: "+chi2.error());
                    else
                    {
                        auto chivalue=chi2_cdf(chi2.value().first, chi2.value().second);
                        if (chivalue<pvalue)
                            return Op_void(false,"Gradient_Variance_Expectancy_Test indicates a low probability of being true: p("+
                                           std::to_string(chivalue)+")<"+ToString(pvalue)+" : "+chi2.error());
                        else
                            return Op_void(true,"Gradient_Variance_Expectancy_Test fall within the expected values: p("+
                                           std::to_string(chivalue)+")>"+ToString(pvalue)+" : "+chi2.error());

                    }

                }


                template<class Simulation_Model,class FIM_Model,  class Data, class ParametersDistribution,class Parameters>
                static auto
                        compute_test
                        (std::ostream& os,const Simulation_Model& sim, const FIM_Model& fim, const Data& e, const ParametersDistribution& prior,const Parameters& p,std::mt19937_64& mt  , std::size_t nsamples, double pvalue)
                {
                    auto mysamples=get_sample(os,sim,fim,e,prior,p,mt,nsamples);
                    auto Gmean=mysamples.value().mean_Gradient();
                    os<<"\n Gmean \n"<<Gmean;

                    auto partialGmean=mysamples.value().partial_mean_Gradient();
                    os<<"\n partialGmean \n";
                    for (std::size_t n=0; n<partialGmean.size(); ++n)
                        os<<"\n"<<partialGmean[n];


                    auto partialGsqr=mysamples.value().partial_sqr_Gradient();
                    auto partialH=mysamples.value().partial_mean_Hessian();

                    auto Gsqr=mysamples.value().mean_sqr_Gradient();
                    os<<"\n Gsqr \n"<<Gsqr;

                    auto H=-mysamples.value().mean_Hessian();
                    os<<"\n H \n"<<H;
                    os<<"\n elemDiv(Gsqr,H)\n"<<elemDiv(Gsqr,H);

                    // lets simulate results...
                    auto resres=sample_Cov(mt,H,nsamples);
                    M_Matrix<double> Smean=resres.value().first;
                    M_Matrix<double> Ssqr=resres.value().second;

                    os<<"\n TEST gradient_information_test \n";
                    auto gradient_information_test_test=Gradient_Expectancy_Test(Smean,H,nsamples,pvalue,os);

                    os<<gradient_information_test_test.error()<<"\n\n";

                    os<<"\n gradient_information_test \n";
                    auto gradient_information_test=Gradient_Expectancy_Test(Gmean,H,nsamples,pvalue,os);
                    os<<gradient_information_test.error()<<"\n\n";
                    os<<"\n partial gradient_information_test \n";
                    for (std::size_t n=0; n<partialGmean.size(); ++n)
                    {
                        auto partial_gradient_information_test=Gradient_Expectancy_Test(partialGmean[n],partialH[n],nsamples,pvalue,os);
                        os<<partial_gradient_information_test.error()<<"\n\n";

                    }

                    os<<"\n TEST gradient_variance_test \n";
                    auto gradient_variance_test_test=Gradient_Expectancy_Test(Smean,Ssqr,nsamples,pvalue,os);
                    os<<gradient_variance_test_test.error()<<"\n\n";

                    os<<"\n gradient_variance_test \n";
                    auto gradient_variance_test=Gradient_Expectancy_Test(Gmean,Gsqr,nsamples,pvalue,os);
                    os<<gradient_variance_test.error()<<"\n\n";
                    os<<"\n partial gradient_variance_test \n";
                    for (std::size_t n=0; n<partialGmean.size(); ++n)
                    {
                        auto partial_gradient_variance_test=Gradient_Expectancy_Test(partialGmean[n],partialGsqr[n],nsamples,pvalue,os);
                        os<<partial_gradient_variance_test.error()<<"\n\n";

                    }


                    os<<"\n TEST gradient_variance_information_test \n";
                    auto gradient_variance_information_test_test=Gradient_Variance_Expectancy_Test(Ssqr,H,nsamples,pvalue,os);
                    os<<gradient_variance_information_test_test.error()<<"\n\n";

                    os<<"\n gradient_variance_information_test \n";
                    auto gradient_variance_information_test=Gradient_Variance_Expectancy_Test(Gsqr,H,nsamples,pvalue,os);
                    os<<gradient_variance_information_test.error()<<"\n\n";;

                    Op_void all_tests=gradient_information_test<<gradient_variance_information_test<<gradient_variance_test;
                    auto par=std::vector({prior.Parameter_to_tr(p)});
                    return Likelihood_Analisis(prior,par,mysamples.value().getSimulations(),mysamples.value().getLikelihoods());
                }



            };


        }

#endif // MYLIKELIHOOD_H
