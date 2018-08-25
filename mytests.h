#ifndef MYTESTS_H
#define MYTESTS_H
#include <cmath>
#include <limits>

struct invariant{};


template<class... > struct class_Invariants;

struct are_Equal

{

    static bool test(double x, double y, double eps=std::sqrt(std::numeric_limits<double>::epsilon()*100))
    {
       if  (std::abs(x-y)/std::max(std::abs(x), 1.0)<eps)
            return true;
        else
            return false;
    }

    template<class Matrix>
    static bool test(const Matrix& one, const Matrix& other, double eps=std::sqrt(std::numeric_limits<double>::epsilon())*100)
    {
        auto dif=one-other;
        if ((one.ncols()!=other.ncols())||(one.nrows()!=other.nrows()))
            return false;
        for (std::size_t i=0; i<one.nrows(); ++i)
        {
            for (std::size_t j=0; j<one.ncols(); ++j)
            if (!test(one(i,j),other(i,j),eps))
                return false;
        }
        return true;
    }

};



struct are_zero: public invariant

{

    static bool test(double x, double conditionNumber=1,double eps=std::numeric_limits<double>::epsilon()*10)
    {
       if  (std::abs(x)*conditionNumber<eps)
            return true;
        else
            return false;
    }

    template<class C>
    static bool test(const C& one, double conditionNumber=1, double eps=std::numeric_limits<double>::epsilon()*10)
    {
        for (std::size_t i=0; i<one.size(); ++i)
        {
            if (!test(one[i],conditionNumber,eps))
                return false;
        }
        return true;
    }

};




struct are_non_negative: public invariant
{
    static bool test(double x, double eps=std::numeric_limits<double>::epsilon())
    {
        if (x+eps>0)
            return true;
        else
            return false;

   }

    template<class C>
    static auto test(const C& x, double eps=std::numeric_limits<double>::epsilon())->decltype (x.size(),std::declval<C>()[0],true)
    {
        for (std::size_t i=0; i<x.size(); ++i)
            if (!test(x[i],eps))
                return false;
        return true;
    }

};

struct are_in_range: public invariant
{
    static bool test(double x, double min, double max, double eps=std::numeric_limits<double>::epsilon())
    {
        if ((x+eps>min)&&(x-eps<max))
            return true;
        else
            return false;

   }

    template<class C>
    static auto test(const C& x, double min, double max,double eps=std::numeric_limits<double>::epsilon())->decltype (x.size(),std::declval<C>()[0],true)
    {
        for (std::size_t i=0; i<x.size(); ++i)
            if (!test(x[i],min,max,eps))
                return false;
        return true;
    }

};




#endif // MYTESTS_H
