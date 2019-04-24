#ifndef ERROR_ANALYSIS_H
#define ERROR_ANALYSIS_H

#include<cmath>
#include<limits>

template<typename T> class E;


class E<double>
{
private:
    double value_;
    double variance_;
  public:
      E(double val): value_{val}, variance_{std::numeric_limits<double>::epsilon()*val}
    {
          variance_*=variance_;
    }
    E(double val, double variance):value_{val}, variance_{variance}
    {
    }
    double value()const {return value_;}
    double variance()const { return variance_;}
    double relative_variance()const { return variance()/value()/value();}
    double stddev()const {return std::sqrt(variance_);}

};


inline E<double>
operator+(const E<double>& one, const E<double>& two)
{
    return E(one.value()+two.value(),one.variance()+two.variance());
}

inline E<double>
operator-(const E<double>& one, const E<double>& two)
{
    return E(one.value()-two.value(),one.variance()+two.variance());
}

inline E<double>
operator*(const E<double>& one, const E<double>& two)
{
    double f=one.value()*two.value();
    double varf=f*(one.relative_variance()+two.relative_variance());
    return E(f,varf);
}

inline E<double>
operator/(const E<double>& one, const E<double>& two)
{
    double f=one.value()/two.value();
    double varf=f*(one.relative_variance()+two.relative_variance());
    return E(f,varf);
}







#endif // ERROR_ANALYSIS_H
