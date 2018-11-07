#ifndef MYTESTS_H
#define MYTESTS_H
#include <cmath>
#include <limits>
#include <iostream>

#include "myoptional.h"
#include "myfields.h"
struct invariant{};


template<class... > struct class_Invariants;

template<bool, class, typename = void> struct are_Equal;
template<bool, class > struct are_zero;
template<bool, class > struct are_finite;
template<bool, class > struct are_non_negative;
template<bool, class > struct are_non_positive;

template<bool, class > struct are_in_range;
template<bool, class > struct are_not_less;
template<bool, class > struct are_not_more;


template <bool output>
class are_Equal<output,double>
{
public:
    template<class ostream=std::ostream>
    bool test(double x, double y,ostream& s=std::cerr)const
    {
      if  ((std::isnan(x)==std::isnan(y))
          &&((std::isnan(x)==std::isnan(y)
               ||(x==y)
               ||(std::abs(x-y)<=absolute_error())
               ||(std::abs(x-y)/std::abs(x+y))<relative_error())))
            return true;
        else
        {
            {
                if constexpr(output){
                    s<<"\n not equal!!";
                    if (std::isnan(x)!=std::isnan(y))
                      s<<" nan value: "<<"  x="<<x<<"  y="<<y;
                    if  ((std::abs(x-y)>=absolute_error()))
                        s<<"\tabsolute="<<absolute_error()<<"  x="<<x<<"  y="<<y<<"  abs(x-y)="<<std::abs(x-y);
                    if  ((std::abs(x-y)/std::abs(x+y))>=relative_error())

                        s<<"\t relative="<<relative_error()<<"  std::abs(x-y)/std::abs(x+y)"<<std::abs(x-y)/std::abs(x+y);
                }
                return false;

            }
        }
    }



    double relative_error()const {return relative_;}
    double absolute_error()const { return absolute_;}
    are_Equal(double absoluteError=std::numeric_limits<double>::epsilon()*1000, double relativeError=std::numeric_limits<double>::epsilon()*1000):absolute_{absoluteError},relative_{relativeError}{}
private:
    double absolute_=std::numeric_limits<double>::epsilon()*100;
    double relative_=std::numeric_limits<double>::epsilon()*100;
};










template <bool output>
class are_finite<output,double>

{
public:
    template<class ostream>
    bool test(double x,ostream& s=std::cerr)const
    {
        if  (std::isfinite(x))
            return true;
        else {
            if constexpr (output)
                    s<<" not finite!!!  x="<<x;
            return false;
        }
    }
};


template <bool output>
class are_zero<output,double>

{
public:
    template<class ostream>
    bool test(double x,ostream& s=std::cerr)const
    {
        if  (std::isfinite(x)&&std::abs(x)<absolute_error())
            return true;
        else {
            if constexpr (output)
                    s<<" not zero!!! absolute="<<absolute_error()<<" x="<<x;
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
    template<class ostream>
    bool test(double x, double y,ostream& s=std::cerr)const
    {
        if (x+absolute_error()>=y)
        {
            return true;

        }
        else
        {
            if constexpr(output)
                    s<<"\nfails are_not_less test  absolute="<<absolute_error()<<" x="<<x<<" y="<<y<<" x-y="<<x-y<<"\n";
            return false;
        }
    }
    bool test(double x, double y)const
    {
        return test(std::cerr,x,y);
    }


    double absolute_error()const { return absolute_;}
    are_not_less(double absoluteError):absolute_{absoluteError}{}
private:
    double absolute_=std::numeric_limits<double>::epsilon();

};


template <bool output>
struct are_not_more<output,double>
{
    template<class ostream>
    bool test(double x, double y,ostream& os=std::cerr)const
    {
        if (x-absolute_error()<=y)
        {
            return true;

        }
        else
        {
            if constexpr(output)
                    os<<"\nfails are_not_more test  absolute="<<absolute_error()<<" x="<<x<<" y="<<y<<" x-y="<<x-y<<"\n";
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
    template<class ostream>
    bool test(double x,ostream& s=std::cerr)const
    {
        if (x+absolute_error()>=0)
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
struct are_non_positive<output,double>
{
    template<class ostream>
    bool test(double x,ostream& os)const
    {
        if (std::isfinite(x)&&(x-absolute_error()<=0))
        {
            return true;

        }
        else
        {
            if constexpr(output)
                    os<<"\nfails are_non_negative test  absolute="<<absolute_error()<<" x="<<x<<"\n";
            return false;
        }
    }


    double absolute_error()const { return absolute_;}
    are_non_positive(double absoluteError):absolute_{absoluteError}{}
private:
    double absolute_=std::numeric_limits<double>::epsilon();

};



template <bool output>
struct are_in_range<output,double>
{
    template<class ostream>
    bool test(double x, ostream& s=std::cerr)const
    {
        if ((missing_) && (x==0)) return true;
        if ((x+absolute_error()>=min())&&(!(x-absolute_error()>max())))
            return true;
        else
        {
            if constexpr(output)
                    s<<"\nfails are_in_range test  absolute_error="<<absolute_error()<<" min="<<min()<<" max="<<max()<<" x="<<x<<" x-min="<<x-min()<<" x-max="<<x-max()<<"\n";
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


template <class T>
bool are_Equal_v(const T &one, const T &other, std::ostream &os) {
  return are_Equal<true, T>().test(one, other, os);
}

template <bool output, class Object, class method>
bool field_test(const grammar::field<Object, method> &m, const Object &one,
                const Object &other, std::ostream &os) {
  std::stringstream ss;
  bool res = are_Equal_v(std::invoke(m.access_method, one),
                         std::invoke(m.access_method, other), ss);
  if constexpr (output)
    if (!res) {
      os << "\nEquality test failed for field " << m.idField
         << " details: " << ss.str();
    }
  return res;
}

template <bool output, class Object>
class are_Equal<
    output, Object,
    std::void_t<decltype(std::declval<Object &>().get_constructor_fields())>> {

public:
  bool test(const Object &one, const Object &other, std::ostream &os) {
    auto fields = one.get_constructor_fields();
    std::stringstream ss;
    bool res = std::apply(
        [&one, &other, &ss](auto &... f) {
          return (field_test<output>(f, one, other, ss) * ...);
        },
        fields);
    if constexpr (output)
      if (!res) {
        os << "\n Equality test failed for "
           << my_trait<Object>::className.str() << " details: " << ss.str();
      }
    return res;
  }
};

#endif // MYTESTS_H
