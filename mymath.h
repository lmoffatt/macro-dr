#ifndef MYMATH_H
#define MYMATH_H

#include <cmath>
#include <cassert>


inline std::size_t base2_floor(std::size_t x)
{
    unsigned n=0;
    while (x>>=1) ++n;
    return 1<<n;
}

inline constexpr double PI = 3.14159265358979323846;

inline constexpr double mynan = std::numeric_limits<double>::quiet_NaN();


inline double log_beta_f(double a,double b) {return std::lgamma(a)+lgamma(b)-lgamma(a+b);}

inline double digamma(double x){
    return std::log(x)
            -0.5*x
            -(1./12.)*std::pow(x,-2)
            +(1./120.)*std::pow(x,-4)
            -(1./252.)*std::pow(x,-6)
            +(1./240.)*std::pow(x,-8)
            -(5./660.)*std::pow(x,-10)
            +(691./32760.)*std::pow(x,-12)
            -(1./12.)*std::pow(x,-14);
    ;

}

inline double normal_pdf(double x) {return 1.0/std::sqrt(2*PI)*std::exp(-0.5*x*x);}




inline double normal_cdf(double x)
{
     return    0.5*(1+std::erf(x/std::sqrt(2)));
}

inline double normal_cdf(double x, double mean, double stddev)
{
    return normal_cdf((x-mean)/stddev);
}

inline double normal_pdf(double x, double mean, double stddev)
{
    return normal_pdf((x-mean)/stddev);
}


inline double chi2_cdf(double x, std::size_t k)
{
    //Chernoff bound
    if constexpr(false)
    {
    double z=x/k;
    if (z<1)
    return std::pow(z*std::exp(1-z),0.5*k);
    else return 1.0-std::pow(z*std::exp(1-z),0.5*k);

    }
    else
        return normal_cdf(std::pow(x/k,1.0/3),1.0-2.0/9.0/k,std::sqrt(2.0/9.0/k));


}
inline double chi2_pdf(double x, std::size_t k)
{
        return normal_pdf(std::pow(x/k,1.0/3),1.0-2.0/9.0/k,std::sqrt(2.0/9.0/k));
}


#endif // MYMATH_H
