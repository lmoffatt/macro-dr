#ifndef MYOPTIMIZATION_H
#define MYOPTIMIZATION_H

#include "myoptional.h"
#include <cmath>
#include <iostream>
#include <map>
struct golden_section_search
{

    template<class F>
    static double opt(const F& f,double a, double c,std::size_t itermax)
    {
        double goldenRatio = (1.0 + std::sqrt(5.0)) / 2.0;
        std::size_t iter=0;
        double fa=f(a); ++iter;
        double fc=f(c); ++iter;
        double b= c - (2.0 - goldenRatio) * (c - a);
        double  fb=f(b); ++iter;
        return step(f,a,fa,b,fb,c,fc,iter,itermax);
    }
    template<class F>
    static double step(const F& f,
                       double a, double fa,
                       double b, double fb,
                       double c, double fc,
                       std::size_t iter,
                       std::size_t itermax)
    {

        double goldenRatio = (1.0 + std::sqrt(5.0)) / 2.0;

        double x;
        if (b < c)
            x = b + (2.0 - goldenRatio) * (c - b);
        else
            x = b - (2.0 - goldenRatio) * (b - a);
        if (iter>itermax)
            return (a+c)/2;
        double  fx=f(x);
        if (fx < fb)
        {if (b < c)
                return step(f, b, fb,x,fx, c,fc, ++iter,itermax);
            else
                return step(f, a, fa,x, fx,b, fb,++iter,itermax);
        }
        else
        {
            if (b < c)

                return step(f, a, fa,b,fb, x,fx, ++iter,itermax);
            else
                return step(f, x, fx, b, fb,c, fc,++iter,itermax);

        }


    };

};



struct zero_binary_search
{
    template<class F>
    static myOptional_t<double> zero(const F& f,double a, double c,double min_x, double min_f,std::size_t itermax)
    {
        std::size_t iter=0;
        double fa=f(a); ++iter;
        double fc=f(c); ++iter;
        //std::cerr<<"\n zero   fa="<<fa<<" fc"<<fc;
        if (std::isfinite(fa)&& std::isfinite(fc))
            return  myOptional_t<double>(false,"invalid starting point");
        if (fa>0)
        {
            if (fc<0)
                return step(f,c,fc,a,fa,min_x,min_f,iter,itermax);
            else if (fc>0)
                return myOptional_t<double>(false,"both starting points positive");
            else return c;
        }
        else if (fa<0)
        {
            if (fc>0)
                return step(f,a,fa,c,fc,min_x,min_f,iter,itermax);
            else if (fc<0)

                return myOptional_t<double>(false,"both starting points negative");
            else return c;
        }
        else return a;
    }
    template<class F>
    static double step(const F& f,
                       double a, double fa,
                       double c, double fc,
                       double min_x, double min_f,
                       std::size_t iter,
                       std::size_t itermax)
    {
        if ((std::abs(a-c)<min_x)||(std::abs(fa-fc)<min_f)||iter>itermax)
            return 0.5*(a+c);
        else
        {
            double b=(a+c)/2;

            double fb=f(b); ++iter;
            //std::cerr<<" b="<<b<<" fb"<<fb<<"\t next::";
            if (!std::isfinite(fb))
                return myOptional_t<double>(false,"invalid value");
            if (fb>0)
                return step(f,a,fa,b,b,min_x, min_f,iter,itermax);
            else if (fb<0)
                return step(f,b,fb,c,fc,min_x, min_f,iter,itermax);
            else return b;
        }

    }

};


struct zero_secant_method
{

    template<class F>
    static double zero(const F& f,double x0, double x1,double min_x, double min_f,std::size_t itermax)
    {
        std::size_t iter=0;
        double f0=f(x0); ++iter;
        double f1=f(x1); ++iter;
        while ((std::abs(x0-x1)>min_x)&&(std::abs(f1)>min_f)&&(itermax>iter))
        {
            double x2=x1-f1*(x1-x0)/(f1-f0);
            f0=f1; x0=x1;
            x1=x2; f1=f(x1);++iter;
            //std::cerr<<" \n zero iter="<<iter<<" f0="<<f0<<" f1="<<f1<<" x0="<<x0<<" x1="<<x1;
        }
        return 0.5*(x1+x0);
    }

};



class nelder_mead_method {

private :
    double alpha=1;
    double beta=0.5;
    double gamma=2;
    double delta=0.5;



    template<class Vector>
    Vector centroid(std::map<double,Vector>& P)
    {
        assert(P.size()>1);
         auto n=P.size()-1;
         auto it=P.begin();
         Vector out=it->second;
         ++it;
         auto end=--P.end();
         for (; it!=end; ++it)  out+=it->second;
         out/=n;
         return out;
    }

    template<class Vector>
    Vector reflect(Vector xh, Vector c)
    {
        return c+(c-xh)*alpha;
    }
    template<class Vector>
    Vector expand(Vector xr, Vector c)
    {
        return c+(xr-c)*gamma;
    }


    template<class Vector>
    Vector contract_outside(Vector xr, Vector c)
    {
        return c+(xr-c)*beta;
    }

    template<class Vector>
    Vector contract_inside(Vector xh, Vector c)
    {
        return c+(xh-c)*beta;
    }

    template<class Vector>
    void replace (std::map<double, Vector>& P, double fc, Vector&& Xc)
    {
        P.erase(--P.end());
        P.emplace(fc,std::move(Xc));
    }

    template<class Function,class Vector>
    std::map<double, Vector> shrink (const Function& f,const std::map<double, Vector>& P)
    {
        std::map<double, Vector> out;
        auto it=P.begin();
        auto xl=it->second;

        out.insert(*it);
        for (++it; it!=P.end(); ++it)
        {
            auto xj=xl+(it->second - xl)*delta;
            double fj=f(xj);
            out.insert(fj,std::move(xj));
        }
        return out;
    }

public :
    template<class Function, class Vector>
    void next(const Function& f, std::map<double, Vector>& P)
    {
        auto c=centroid(P);
        auto it_h=--P.end(); auto fh=it_h->first; auto xh=it_h->second;
        --it_h;
        auto fs=it_h->first;
        it_h=P.begin();
        auto fl=it_h->first;

        auto xr=reflect(xh,c);
        auto fr=f(xr);

        if (fr>=fs)
        {
            Vector xc;
            if (fr>=fh)
                xc=contract_inside(xh,c);
            else
                xc=contract_outside(xr,c);
            auto fc=f(xc);
            if (fc<fh)
                replace(P,fc,xc);
            else
                P=shrink(f,P);

        }
        else if (fr>=fl)
            replace(P,fr,xr);
        else
        {
            auto xe=expand(xr,c);
            auto fe=f(xe);
            if (fe<fr)
                replace(P,fe,std::move(xe));
            else
                replace(P,fr,std::move(xr));

        }
     }






};








#endif // MYOPTIMIZATION_H
