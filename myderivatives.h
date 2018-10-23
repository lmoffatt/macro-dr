#ifndef MYDERIVATIVES_H
#define MYDERIVATIVES_H
#include <utility>
#include <map>
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
   T value;
   Constant(T&& x):value{std::move(x)}{}
   Constant(const T& x):value{x}{}
   
};



#endif // MYDERIVATIVES_H
