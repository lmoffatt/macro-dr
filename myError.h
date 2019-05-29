#ifndef MYERROR_H
#define MYERROR_H

#include <cmath>
#include <limits>

#include "mytypetraits.h"
template <typename T, class Norm, bool diff> class Error;


struct norm_1
{
    static double pow(double x) {return std::abs(x);}
    static double pow_inv(double x) {return x;}
};

struct norm_2
{
    static double pow(double x) {return x*x;}
    static double pow_inv(double x) {return std::sqrt(x);}
};




template <typename T> struct E { typedef Error<T, ::norm_1, true> type; };

template<template<template<class...>class,auto...> class myClass,auto...Ts>
struct E<myClass<C,Ts...>>
{
    typedef myClass<E,Ts...> type;
};

template<>
struct my_template_trait<E>
{
    constexpr static auto className = my_static_string("");

};



template <typename...Ts>
using E_t=typename E<Ts...>::type;


template <> struct E<double> { typedef Error<double, ::norm_1, true> type; };




template <class Norm, bool diff> class Error<double, Norm, diff> {
private:
  double center_;
  double norm_;

public:
  Error operator-() && {
    center_ = -center_;
    return *this;
  }
  Error operator-() const & {
    Error out(*this);
    out.center_ = -center_;
    return out;
  }

  double center() const { return center_; }
  double norm() const { return norm_; }
  double error() const {
      return Norm::pow_inv(norm_);
  }
  Error()=default;
  Error(double val)
      : center_{val}, norm_{Norm::pow(std::numeric_limits<double>::epsilon() * val)} {}
  Error(double val, double err) : center_{val}, norm_{err} {
      assert(std::isnan(norm())||error() >= 0);
  }

  Error &operator+=(const Error &other) {
    center_ += other.center();
    norm_ = norm()+other.norm();
    return *this;
  }
  Error &operator-=(const Error &other) {
    center_ -= other.center();
    norm_ = norm()+other.norm();
    return *this;
  }
  Error &operator*=(const Error &other) {
    center_ *= other.center();
    norm_ = norm() * Norm::pow(other.center())+ Norm::pow(center()) * other.norm();
    return *this;
  }
  Error &operator/=(const Error &other) {
    center_ /= other.center();
    norm_ = norm() / Norm::pow(other.center())+other.norm()*Norm::pow(center() / (other.center()* other.center()));
    return *this;
  }

  template <typename Number> Error &operator+=(double other) {
    center_ += other;
    return *this;
  }
  template <typename Number> Error &operator-=(double other) {
    center_ -= other;
    return *this;
  }
  template <typename Number> Error &operator*=(double other) {
    center_ *= other;
    norm_ *= other;
    return *this;
  }

  template <typename Number> Error &operator/=(double other) {
    center_ /= other;
    norm_ /= other;
    return *this;
  }
};


template <class Norm, bool diff>
bool operator==(const Error<double, Norm,diff> &one, const Error<double, Norm, diff> &two) {
    return Norm::pow(one.center()-two.center())/(one.norm()+two.norm())<=1.0;
}
template <class Norm, bool diff>
bool operator==(const Error<double, Norm,diff> &one, double two) {
    if (one.center()==two)
        return true;
    else if (Norm::pow(one.center()-two)/(one.norm())<=1.0)
        return true;
    else
        return false;
}

template <class Norm, bool diff>
bool operator==(double one,const Error<double, Norm,diff> &two) {
    return two==one;
}

template <class Norm, bool diff ,typename...Ts>
bool are_Equal_v(const Error<double,Norm,diff> &one, const Error<double,Norm,diff> &other, std::ostream &os, Ts...context) {
    if (!(one==other))
    {
        (os<<...<<context);
        return false;
    }
    else return true;
}



template <class Norm, bool diff>
 Error<double, Norm,diff> operator+(const Error<double, Norm,diff> &one,
                               const Error<double, Norm,diff> &two) {
  Error<double, Norm,diff> out(one);
  out += two;
  return out;
}

template <class Norm, bool diff>
Error<double, Norm,diff> operator+(const Error<double, Norm,diff> &one,
                                     double two) {
    return Error<double, Norm,diff> (one.center()+two,one.norm());
}

template <class Norm, bool diff>
Error<double, Norm,diff> operator+(double two,const Error<double, Norm,diff> &one
                                    ) {
    return Error<double, Norm,diff> (two+one.center(),one.norm());
}

template <class Norm, bool diff >
Error<double, Norm,diff> operator*(const Error<double, Norm,diff> &one,
                                  double two) {
    return Error<double, Norm,diff> (one.center()*two,one.norm()*Norm::pow(two));
}

template <class Norm, bool diff >
Error<double, Norm,diff> operator*(double two,const Error<double, Norm,diff> &one
                                    ) {
    return Error<double, Norm,diff> (two*one.center(),Norm::pow(two)*one.norm());
}

template <class Norm, bool diff >
Error<double, Norm,diff> operator/(const Error<double, Norm,diff> &one,
                                    double two) {
    return Error<double, Norm,diff> (one.center()/two,one.norm()/Norm::pow(two));
}

template <class Norm, bool diff >
Error<double, Norm,diff> operator/(double two,const Error<double, Norm,diff> &one
                                    ) {
    return Error<double, Norm,diff> (two/one.center(),one.norm()*two/Norm::pow(one.center())/Norm::pow(one.center()));
}





template <class Norm, bool diff>
Error<double, Norm,diff> operator-(const Error<double, Norm,diff> &one,
                                    const Error<double, Norm,diff> &two) {
    Error<double, Norm,diff> out(one);
    out -= two;
    return out;
}
template <class Norm, bool diff>
Error<double, Norm,diff> operator*(const Error<double, Norm,diff> &one,
                                    const Error<double, Norm,diff> &two) {
    Error<double, Norm,diff> out(one);
    out *= two;
    return out;
}
template <class Norm, bool diff>
Error<double, Norm,diff> operator/(const Error<double, Norm,diff> &one,
                                    const Error<double, Norm,diff> &two) {
    Error<double, Norm,diff> out(one);
    out /= two;
    return out;
}




template <class T,class Norm, bool diff>
auto const& center(const Error<T,Norm,diff>& x)
{
    return x.center();
}

template <class T>
auto const& center(const T& x)
{
    return x;
}

template<typename T, class Norm, bool diff>
std::ostream& operator<<(std::ostream& os, const Error<T,Norm,diff>& x)
{
    os<<x.center()<<"+/-"<<x.norm()<<"\n";
    return os;
}


template <class Norm, bool diff>
Error<double, Norm,diff> exp(const Error<double, Norm,diff> &one) {
    auto c=std::exp(one.center());
    double n=Norm::pow(c)*one.norm();
    return Error<double, Norm,diff> (c,n);

}

template <class Norm, bool diff>
Error<double, Norm,diff> sqr(const Error<double, Norm,diff> &one) {
    auto c=sqr(one.center());
    double n=Norm::pow(2*one.center())*one.norm();
    return Error<double, Norm,diff> (c,n);

}

template <class Norm, bool diff>
Error<double, Norm,diff> sqrt(const Error<double, Norm,diff> &one) {
    auto c=std::sqrt(one.center());
    double n=one.norm()/Norm::pow(2.0*std::sqrt(one.center()));
    return Error<double, Norm,diff> (c,n);

}


template <class Norm, bool diff>
Error<double, Norm,diff> log(const Error<double, Norm,diff> &one) {
    auto log=std::log(one.center());
    double n=Norm::pow(1.0/one.center())*one.norm();
    return Error<double, Norm,diff> (log,n);

}



template <class Norm, bool diff>
Error<double, Norm,diff> abs(const Error<double, Norm,diff> &one) {
    auto c=std::abs(one.center());
    return Error<double, Norm,diff> (c,one.norm());

}

template <class Norm, bool diff>
Error<double, Norm,diff> expm1(const Error<double, Norm,diff> &one) {
    auto c=std::expm1(one.center());
    return Error<double, Norm,diff> (c,one.norm());

}

template <class Norm, bool diff>
Error<double, Norm,diff> max(const Error<double, Norm,diff> &one,const Error<double, Norm,diff> &two) {
    if (one.center()>two.center())
    {
        if (one==two)
            return Error<double,Norm,diff>(one.center(),one.norm()+two.norm());
        else
            return one;
    }
    else
    {
        if (one==two)
            return Error<double,Norm,diff>(two.center(),one.norm()+two.norm());
        else
            return two;
    }


}

template <class Norm, bool diff >
Error<double, Norm,diff> max(const Error<double, Norm,diff> &one,double two) {
    if (one.center()>two)
    {
            return one;
    }
    else
    {
        if (one==two)
            return Error<double,Norm,diff>(two,one.norm());
        else
            return Error<double,Norm,diff>(two);
    }


}


#endif // MYERROR_H
