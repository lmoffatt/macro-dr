#ifndef MYTESTS_H
#define MYTESTS_H
#include <cmath>
#include <limits>
#include <iostream>
struct invariant{};


template<class... > struct class_Invariants;

template<bool, class > struct are_Equal;
template<bool, class > struct are_zero;
template<bool, class > struct are_non_negative;
template<bool, class > struct are_in_range;
template<bool, class > struct are_not_less;
template<bool, class > struct are_not_more;


template <bool output>
class are_Equal<output,double>
{
public:
    bool test(double x, double y)const
    {
        if  ((std::abs(x-y)<absolute_error())||(std::abs(x-y)/std::abs(x+y))<relative_error())
            return true;
        else
        {
            {
                if constexpr(output){
                    std::cerr<<"\n not equal!! absolute="<<absolute_error()<<"x="<<x<<" y="<<y<<" abs(x-y)="<<std::abs(x-y);
                    std::cerr<<" relative="<<relative_error()<<"  std::abs(x-y)/std::abs(x+y)"<<std::abs(x-y)/std::abs(x+y);
                }
                return false;

            }
        }
    }

    double relative_error()const {return relative_;}
    double absolute_error()const { return absolute_;}
    are_Equal(double absoluteError=std::numeric_limits<double>::epsilon(), double relativeError=std::numeric_limits<double>::epsilon()):absolute_{absoluteError},relative_{relativeError}{}
private:
    double absolute_=std::numeric_limits<double>::epsilon();
    double relative_=std::numeric_limits<double>::epsilon();
};





template <bool output>
class are_zero<output,double>

{
public:
    bool test(double x)const
    {
        if  (std::abs(x)<absolute_error())
            return true;
        else {
            if constexpr (output)
                    std::cerr<<" not zero!!! absolute="<<absolute_error()<<" x="<<x;
            return false;
        }
    }

    double absolute_error()const { return absolute_;}
    are_zero(double absoluteError):absolute_{absoluteError}{}
private:
    double absolute_=std::numeric_limits<double>::epsilon();
};


template <bool output>
struct are_not_less<output,double>
{

     bool test(double x, double y)const
    {
        if (x+absolute_error()>y)
        {
            return true;

        }
        else
        {
            if constexpr(output)
                    std::cerr<<"\nfails are_not_less test  absolute="<<absolute_error()<<" x="<<x<<" y="<<y<<" x-y="<<x-y<<"\n";
            return false;
        }
    }

    double absolute_error()const { return absolute_;}
    are_not_less(double absoluteError):absolute_{absoluteError}{}
private:
    double absolute_=std::numeric_limits<double>::epsilon();

};


template <bool output>
struct are_not_more<output,double>
{

     bool test(double x, double y)const
    {
        if (x-absolute_error()<y)
        {
            return true;

        }
        else
        {
            if constexpr(output)
                    std::cerr<<"\nfails are_not_more test  absolute="<<absolute_error()<<" x="<<x<<" y="<<y<<" x-y="<<x-y<<"\n";
            return false;
        }
    }

    double absolute_error()const { return absolute_;}
    are_not_more(double absoluteError):absolute_{absoluteError}{}
private:
    double absolute_=std::numeric_limits<double>::epsilon();

};




template <bool output>
struct are_non_negative<output,double>
{

     bool test(double x)const
    {
        if (x+absolute_error()>0)
        {
            return true;

        }
        else
        {
            if constexpr(output)
                    std::cerr<<"\nfails are_non_negative test  absolute="<<absolute_error()<<" x="<<x<<"\n";
            return false;
        }
    }

    double absolute_error()const { return absolute_;}
    are_non_negative(double absoluteError):absolute_{absoluteError}{}
private:
    double absolute_=std::numeric_limits<double>::epsilon();

};

template <bool output>
struct are_in_range<output,double>
{
     bool test(double x)const
    {
         if ((missing_) && (x==0)) return true;
         if ((x+absolute_error()>min())&&(x-absolute_error()<max()))
            return true;
        else
        {
            if constexpr(output)
                    std::cerr<<"\nfails are_in_range test  absolute_error="<<absolute_error()<<" min="<<min()<<" max="<<max()<<" x="<<x<<" x-min="<<x-min()<<" x-max="<<x-max()<<"\n";
            return false;
        }

    }
    double min()const {return min_;}
    double max()const {return max_;}

    double absolute_error()const { return absolute_;}
    are_in_range(bool missing,double min, double max,double absoluteError):missing_{missing},min_{min},max_{max},absolute_{absoluteError}{}
private:
    bool missing_;
    double min_;
    double max_;
    double absolute_=std::numeric_limits<double>::epsilon();


};






#endif // MYTESTS_H
