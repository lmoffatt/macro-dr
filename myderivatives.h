#ifndef MYDERIVATIVES_H
#define MYDERIVATIVES_H
#include <utility>
#include <map>
#include "mytypetraits.h"
template< class...>
class Derivative;

template<class X>
struct myDerivative
{
    typedef Derivative<X> type;
};


template <class X>
using Derivative_t=typename myDerivative<X>::type;


template <typename K, typename T, typename ...X>
struct myDerivative<std::map<K,T,X...>>
{
    typedef std::map<K, Derivative_t<T>> type;
};



class D{};


template <class T>
struct Constant
{
  constexpr static auto className=my_static_string("Constant")+my_trait<T>::className;

   T value;
   Constant(T&& x):value{std::move(x)}{}
   Constant(const T& x):value{x}{}
   Constant()=default;
   
};

namespace std{
template <class T> struct numeric_limits<Constant<T>>: public numeric_limits<T> {
  typedef numeric_limits<T> base_type;
  using base_type::epsilon;
};
template <class T>
struct numeric_limits<Derivative<T>> : public numeric_limits<T> {
  typedef numeric_limits<T> base_type;
//  using base_type::epsilon;
};
}
#endif // MYDERIVATIVES_H
