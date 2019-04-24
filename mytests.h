#ifndef MYTESTS_H
#define MYTESTS_H
#include <cmath>
#include <limits>
#include <iostream>

#include "myoptional.h"
#include "myfields.h"
#include "mySerializer.h"
struct invariant {};

template <class...> struct class_Invariants;

template <bool, class, class =void> struct are_Equal;
template <bool, class> struct are_zero;
template <bool, class> struct are_finite;
template <bool, class> struct are_non_negative;
template <bool, class> struct are_non_positive;

template <bool, class> struct are_in_range;
template <bool, class> struct are_not_less;
template <bool, class> struct are_not_more;

template <bool output> struct are_Equal<output, double> {
public:
  template <class ostream = std::ostream>
  bool test(double x, double y, ostream &s = std::cerr) const {
    if ((std::isnan(x) == std::isnan(y)) &&
        ((std::isnan(x) && std::isnan(y)) || (x == y) ||
         (std::abs(x - y) <= absolute_error()) ||
         (std::abs(x - y) / std::abs(x + y)) < relative_error()))
      return true;
    else {
      {
        if constexpr (output) {
          s << "\n not equal!!";
          if (std::isnan(x) != std::isnan(y))
            s << " nan value: "
              << "  x=" << x << "  y=" << y;
          if ((std::abs(x - y) >= absolute_error()))
            s << "\tabsolute=" << absolute_error() << "  x=" << x << "  y=" << y
              << "  abs(x-y)=" << std::abs(x - y);
          if ((std::abs(x - y) / std::abs(x + y)) >= relative_error())

            s << "\t relative=" << relative_error()
              << "  std::abs(x-y)/std::abs(x+y)"
              << std::abs(x - y) / std::abs(x + y);
        }
        return false;
      }
    }
  }

  double relative_error() const { return relative_; }
  double absolute_error() const { return absolute_; }
  are_Equal(double absoluteError = std::numeric_limits<double>::epsilon() *
                                   1000,
            double relativeError = std::numeric_limits<double>::epsilon() *
                                   1000)
      : absolute_{absoluteError}, relative_{relativeError} {}

private:
  double absolute_ = std::numeric_limits<double>::epsilon() * 100;
  double relative_ = std::numeric_limits<double>::epsilon() * 100;
};

template <bool output>
struct are_finite<output, double>

{
public:
  template <class ostream> bool test(double x, ostream &s = std::cerr) const {
    if (std::isfinite(x))
      return true;
    else {
      if constexpr (output)
        s << " not finite!!!  x=" << x;
      return false;
    }
  }
};

template <bool output>
struct are_zero<output, double>

{
public:
  template <class ostream> bool test(double x, ostream &s = std::cerr) const {
    if (std::isfinite(x) && std::abs(x) < absolute_error())
      return true;
    else {
      if constexpr (output)
        s << " not zero!!! absolute=" << absolute_error() << " x=" << x;
      return false;
    }
  }

  double absolute_error() const { return absolute_; }
  are_zero(double absoluteError) : absolute_{absoluteError} {}

private:
  double absolute_ = std::numeric_limits<double>::epsilon();
};

template <bool output> struct are_not_less<output, double> {
  template <class ostream>
  bool test(double x, double y, ostream &s = std::cerr) const {
    if (x + absolute_error() >= y) {
      return true;

    } else {
      if constexpr (output)
        s << "\nfails are_not_less test  absolute=" << absolute_error()
          << " x=" << x << " y=" << y << " x-y=" << x - y << "\n";
      return false;
    }
  }
  bool test(double x, double y) const { return test(std::cerr, x, y); }

  double absolute_error() const { return absolute_; }
  are_not_less(double absoluteError) : absolute_{absoluteError} {}

private:
  double absolute_ = std::numeric_limits<double>::epsilon();
};

template <bool output> struct are_not_more<output, double> {
  template <class ostream>
  bool test(double x, double y, ostream &os = std::cerr) const {
    if (x - absolute_error() <= y) {
      return true;

    } else {
      if constexpr (output)
        os << "\nfails are_not_more test  absolute=" << absolute_error()
           << " x=" << x << " y=" << y << " x-y=" << x - y << "\n";
      return false;
    }
  }

  double absolute_error() const { return absolute_; }
  are_not_more(double absoluteError) : absolute_{absoluteError} {}

private:
  double absolute_ = std::numeric_limits<double>::epsilon();
};

template <bool output> struct are_non_negative<output, double> {
  template <class ostream> bool test(double x, ostream &s = std::cerr) const {
    if (x + absolute_error() >= 0) {
      return true;

    } else {
      if constexpr (output)
        s << "\nfails are_non_negative test  absolute="
                  << absolute_error() << " x=" << x << "\n";
      return false;
    }
  }

  double absolute_error() const { return absolute_; }
  are_non_negative(double absoluteError) : absolute_{absoluteError} {}

private:
  double absolute_ = std::numeric_limits<double>::epsilon();
};
template <bool output> struct are_non_positive<output, double> {
  template <class ostream> bool test(double x, ostream &os) const {
    if (std::isfinite(x) && (x - absolute_error() <= 0)) {
      return true;

    } else {
      if constexpr (output)
        os << "\nfails are_non_negative test  absolute=" << absolute_error()
           << " x=" << x << "\n";
      return false;
    }
  }

  double absolute_error() const { return absolute_; }
  are_non_positive(double absoluteError) : absolute_{absoluteError} {}

private:
  double absolute_ = std::numeric_limits<double>::epsilon();
};

template <bool output> struct are_in_range<output, double> {
  template <class ostream> bool test(double x, ostream &s = std::cerr) const {
    if ((missing_) && (x == 0))
      return true;
    if ((x + absolute_error() >= min()) && (!(x - absolute_error() > max())))
      return true;
    else {
      if constexpr (output)
        s << "\nfails are_in_range test  absolute_error=" << absolute_error()
          << " min=" << min() << " max=" << max() << " x=" << x
          << " x-min=" << x - min() << " x-max=" << x - max() << "\n";
      return false;
    }
  }
  double min() const { return min_; }
  double max() const { return max_; }

  double absolute_error() const { return absolute_; }
  are_in_range(bool missing, double min, double max, double absoluteError)
      : missing_{missing}, min_{min}, max_{max}, absolute_{absoluteError} {}

private:
  bool missing_;
  double min_;
  double max_;
  double absolute_ = std::numeric_limits<double>::epsilon();
};

template <class T, typename...Ts>
bool are_Equal_v(const T &one, const T &other, std::ostream &os, Ts...context) {
  if (!are_Equal<true, T>(std::sqrt(std::numeric_limits<double>::epsilon()*1e10),std::sqrt(std::numeric_limits<double>::epsilon()*1e10)).test(one, other, os))
  {
    (os<<...<<context);
    return false;
  }
  else return true;
}
template <class T, typename...Ts>
auto are_Equal_v(const T &one, const T &other, double eps, double epsf,std::ostream &os, Ts...context)->std::enable_if_t<!std::is_pointer_v<T>,bool>
{
  if (!are_Equal<true, T>(eps,epsf).test(one, other, os))
  {
      using namespace io;
    (os<<...<<context);
    return false;
  }
  else return true;
}

template <class T, typename...Ts>
bool are_Equal_v(const T *one, const T * other, double eps, double epsf,std::ostream &os, Ts...context) {
  if (!are_Equal<true, T>(eps,epsf).test(one, other, os))
  {
    (os<<...<<context);
    return false;
  }
  else return true;
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


template <bool output, class Object, class method>
bool field_test(const grammar::field<Object, method> &m, const Object &one,
                const Object &other, double eps, double epsf,std::ostream &os) {
  std::stringstream ss;
  bool res = are_Equal_v(std::invoke(m.access_method, one),
                         std::invoke(m.access_method, other), eps,epsf,ss);
  if constexpr (output)
    if (!res) {
      os << "\nEquality test failed for field " << m.idField
         << " details: " << ss.str();
    }
  return res;
}


template <bool output, class Object>
struct are_Equal<
    output, Object,
    std::void_t<decltype(std::decay_t<Object>::get_constructor_fields())>> {
private :
  double absolute_; double relative_;
public:
  are_Equal(double absolute, double relative):absolute_{absolute},relative_{relative}{}
  bool test(const Object &one, const Object &other, std::ostream &os) {
    return test(one,other,absolute_,relative_,os);
  }
  bool test(const Object &one, const Object &other, double abs,double rel,std::ostream &os) {
    auto fields = one.get_constructor_fields();
    std::stringstream ss;
    bool res = std::apply(
        [&one, &other, &ss,&abs,&rel](auto &... f) {
          return (field_test<output>(f, one, other,abs,rel, ss) * ...);
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


template <bool output, class Int>
struct are_Equal<
    output, Int,std::enable_if_t<std::is_integral_v<Int>,void>> {

public:
  are_Equal(double, double ) {}
  bool test(const Int &one, const Int &other, std::ostream &os) {
    bool res = one==other;
    if constexpr (output)
      if (!res) {
        os << "\n Equality test failed for "
           << my_trait<Int>::className.str() << " one = " << one<<" other= "<<other;
      }
    return res;
  }
};

template <bool output, class C>
struct are_Equal<
    output, C,std::enable_if_t<is_std_container<C>::value,void>> {

  double absolute_;
  double relative_;

public:
  are_Equal(double abs, double rel) : absolute_{abs},relative_{rel}{}
  bool test(const C &one, const C &two, double abs, double rel,std::ostream &os) {
    bool res=true;
    auto it1=one.begin();
    auto it2=two.begin();
    std::size_t pos=0;
    for (; (it1!=one.end())&&(it2!=two.end()); ++it1,++it2)
    {

      if (!are_Equal_v(*it1,*it2,abs,rel,os))
      {
        os<<"differ at pos: "<<pos<<"\t";
        res=false;
      }
      ++pos;
    }
    if (it1!=one.end()|| it2!=two.end())
    {
      os<<"Differ in size!!: "<<one.size()<<" vs "<<two.size()<<"\n";
      res=false;
    }
    return res;
  }
  bool test(const C &one, const C &two, std::ostream &os) {
    return test(one,two,absolute_,relative_,os);
  }
};



template <bool output, class...Ts>
struct are_Equal<output,std::tuple<Ts...>>
{
    double absolute_;
double relative_;

public:
  are_Equal(double abs, double rel) : absolute_{abs},relative_{rel}{}
  template<std::size_t...Is>
  bool test(const std::tuple<Ts...> &one, const std::tuple<Ts...> &two, std::ostream &os, std::index_sequence<Is...>)
  {
    return (are_Equal_v(std::get<Is>(one),std::get<Is>(two),absolute_,relative_,os)&&...);
  }

  bool test(const std::tuple<Ts...> &one, const std::tuple<Ts...> &two, std::ostream &os) {

    return test(one,two,os,std::index_sequence_for<Ts...>());
  }
};

template <bool output, class T, class K>
struct are_Equal<output,std::pair<T,K>,void>
{
  double absolute_;
  double relative_;

public:
  are_Equal(double abs, double rel) : absolute_{abs},relative_{rel}{}

  bool test(const std::pair<T,K> &one, const  std::pair<T,K> &two, std::ostream &os) {

    return are_Equal_v(one.first,two.first,absolute_,relative_,os)&&are_Equal_v(one.second,two.second,absolute_,relative_,os);
  }
};


#endif // MYTESTS_H
