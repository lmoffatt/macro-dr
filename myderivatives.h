#ifndef MYDERIVATIVES_H
#define MYDERIVATIVES_H
#include <utility>
template< class...>
class Derivative;

template<class X>
struct myDerivative
{
    typedef Derivative<X> type;
};


template <class X>
using Derivative_t=typename myDerivative<X>::type;

class D{};


template <class T>
struct Constant
{
   T value;
   Constant(T&& x):value{std::move(x)}{}
 };



#endif // MYDERIVATIVES_H
