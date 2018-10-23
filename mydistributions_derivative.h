#ifndef MYDISTRIBUTIONS_DERIVATIVE_H
#define MYDISTRIBUTIONS_DERIVATIVE_H

#include "myDistributions.h"
#include "matrixderivative.h"

template <>
struct Derivative<Probability_distribution>: public Probability_distribution
{

    template<class dP>
    dP static normalize_derivative(dP&& dp)
    {
        double dpos=0;
        double dneg=0;
        for (std::size_t i=0; i<dp.size(); ++i)
        {
            if(dp[i]>0) dpos+=dp[i];
            else dneg-=dp[i];
        }
        double fpos=dneg/(dpos+dneg);
        double fneg=dpos/(dpos+dneg);
        for (std::size_t i=0; i<dp.size(); ++i)
        {
            if(dp[i]>0) dp[i]*=fpos;
            else dp[i]*=fneg;
        }
        return dp;
    }


    template<class P>
    Derivative<P> static normalize(Derivative<P>&& p, double min_p)
    {
        double sum=0;
        for (std::size_t i=0; i<p.f().size(); ++i)
        {
            if (p.f()[i]<min_p)
            {
                p.f()[i]=0;
                p.dfdx().transform([&i](auto &m){m[i]=0;});
            }
            else
            {
                if (p.f()[i]+min_p>1)
                {
                    p.f()[i]=1;
                    p.dfdx().transform([&i](auto &m){m[i]=0;});
                }
                sum+=p.f()[i];
            }
        }
        for (std::size_t i=0; i<p.f().size(); ++i)
            p.f()[i]=p.f()[i]/sum;
        p.dfdx().transform([](auto &m ){m=normalize_derivative(std::move(m));});
        return p;
    }
};

template<>
struct Derivative<Probability_distribution_covariance>: public Probability_distribution_covariance
{
    static M_Matrix<double> normalize_derivative(M_Matrix<double>&& dp, const std::set<std::size_t>& non_zero_i)
    {
        for (auto i:non_zero_i)
        {

            double dpos=0;
            double dneg=0;
            for (auto j: non_zero_i)
            {
                if(dp(i,j)>0) dpos+=dp(i,j);
                else dneg-=dp(i,j);
            }
            double fpos=dneg/(dpos+dneg);
            double fneg=dpos/(dpos+dneg);
            for (auto j: non_zero_i)
            {
                if(dp(i,j)>0) dp(i,j)*=fpos;
                else dp(i,j)*=fneg;
            }

        }
        return dp;

    }


    static Derivative<M_Matrix<double>> normalize(Derivative<M_Matrix<double>>&& p,double min_p)
    {
        std::set<std::size_t> non_zero_i;

        for (std::size_t i=0; i<p.f().nrows(); ++i) {
            if(p.f()(i,i)<min_p)
            {
                for (std::size_t j=0; j<p.f().ncols(); ++j)
                {
                    p.f()(i,j)=0;
                    p.dfdx().transform([&i,&j](auto& m){m(i,j)=0;});
                }
            }
            else if(p.f()(i,i)+min_p>1.0)
            {
                for (std::size_t n=0; n<p.f().size(); ++n)
                {
                    p.f()[n]=0;
                    p.dfdx().transform([&n](auto &m){m[n]=0;});
                }
                p.f()(i,i)=1;

                return p;
            }
            else
                non_zero_i.insert(i);
        }
        for (auto i:non_zero_i)
        {

            double sum=0;
            for (auto j: non_zero_i)
                if (j!=i)
                    sum+=p.f()(i,j);

            if (-sum!=p.f()(i,i))
            {
                auto sum_new=(sum-p.f()(i,i))*0.5;
                double f=sum_new/sum;
                p.f()(i,i)=-sum_new;
                for (auto j: non_zero_i)
                    if (i!=j)
                        p.f()(i,j)*=f;

            }
        }
        p.dfdx().transform([&non_zero_i](auto &m){m=normalize_derivative(std::move(m),non_zero_i);});

        return p;
    }



};


template<>
struct Derivative<variance_value>: public variable_value
{


    static Derivative<double> adjust(Derivative<double>&& value, double min_variance)
    {
        value.f()= std::max(value.f(),min_variance);
        return value;
    }

};



#endif // MYDISTRIBUTIONS_DERIVATIVE_H
