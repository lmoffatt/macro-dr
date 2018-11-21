#ifndef MYDERIVATIVES_H
#define MYDERIVATIVES_H
#include <utility>
#include <map>
#include <functional>
#include "mytypetraits.h"
template< class...>
class Derivative;

template<class X,class=void>
struct myDerivative
{
    typedef Derivative<X> type;
};


template <class X>
      struct myDerivative<X, std::void_t<typename X::Derivative>> {
  typedef typename X::Derivative type;
};

template <class X>
using Derivative_t=typename myDerivative<std::decay_t<X>>::type;


template <typename K, typename T, typename ...X>
struct myDerivative<std::map<K,T,X...>>
{
    typedef std::map<K, Derivative_t<T>> type;
};



template <class> struct is_Derivative : public std::false_type {};

template <typename T>
struct is_Derivative<Derivative<T>> : public std::true_type {};

template <class D>
constexpr static bool is_Derivative_v = is_Derivative<D>::value;


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



template <class C>
auto Primitive(const Derivative<C> &dy)->decltype (std::declval<const Derivative<C>&>().f())
{
  return dy.f();
}

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
