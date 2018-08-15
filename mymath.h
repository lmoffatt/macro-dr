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



#endif // MYMATH_H
