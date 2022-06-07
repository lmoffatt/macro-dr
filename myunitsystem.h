#ifndef MYUNITSYSTEM_H
#define MYUNITSYSTEM_H

#include "mytypetraits.h"
#include <cstddef>

struct m{};
struct s{};
struct kg{};
struct A{};

template <class m, int N>
struct u{};

template<class ,class>
struct add_if_same_unit;

template <class m, int N, int N2>
struct add_if_same_unit<u<m,N>,u<m,N2>>{
    typedef u<m,N+N2> type ;
};

template <class m,class m2, int N, int N2>
struct add_if_same_unit<u<m,N>,u<m2,N2>>{
    typedef u<m,N> type ;
};

template<class ...>
struct substr_if_same_unit;

template <class m, int N, int N2>
struct substr_if_same_unit<u<m,N>,u<m,N2>>{
    typedef u<m,N-N2> type ;
};

template <class m,class m2, int N, int N2>
struct substr_if_same_unit<u<m,N>,u<m2,N2>>{
    typedef u<m,N> type ;
};



template<class ...>
struct wrap_if_not_unit;

template <class m>
struct wrap_if_not_unit<u<m,0>>{
    typedef std::tuple<> type ;
};

template <class m, int N>
struct wrap_if_not_unit<u<m,N>>{
    typedef std::tuple<u<m,N>> type ;
};



template<class ...> struct p;





template<class ...ms, int...Ns> struct p<u<ms,Ns>...>
{};



template<class ...units>
class unit_system
{
public:
    template<class u>
    static constexpr std::size_t index= index_of_this_type<Cs<units...>,u>::value;


    typedef p<u<units,0>...> scalar_unit;


};
typedef unit_system<m,s,kg,A> my_unit_sytem;

typedef typename my_unit_sytem::scalar_unit my_scalar_unit;




template<class...>struct simplify;
template<class...T>
using simplify_t=typename simplify<T...>::type;


template<class...us>
struct simplify<p<us...>>
{
    typedef transfer_t<
        decltype(std::tuple_cat(std::declval<typename wrap_if_not_unit<us>::type>() ...)),
        p> type;
};


template<class, int>struct mult_exponent;

template<class T, int N>
using mult_exponent_t=typename mult_exponent<T, N>::type;

template<class...ms, int...Ns, int N>
struct mult_exponent<p<u<ms,Ns>...>, N>
{
    typedef p<u<ms,Ns*N>...> type;
};

template<class, int>struct div_exponent;

template<class T, int N>
using div_exponent_t=typename div_exponent<T,N>::type;

template<class...ms, int...Ns, int N>
struct div_exponent<p<u<ms,Ns>...>, N>
{
    static_assert (((Ns % N == 0)&&... &&true), "some exponent is not divisible");
    typedef p<u<ms,Ns/N>...> type;
};




template<class...>struct add_exponent;

template<class...T>
using add_exponent_t=typename add_exponent<T...>::type;


template<class...us, class u>
struct add_exponent<p<us...>, u>
{
    typedef p<
        typename add_if_same_unit<us,u>::type ...
        > type;
};


template<class...us,class u, class... us2>
struct add_exponent<p<us...>, u, us2...>
{
    typedef add_exponent_t<add_exponent_t<p<us...>,u>,us2...> type;
};


template<class...>struct substr_exponent;

template<class...T>
using substr_exponent_t=typename substr_exponent<T...>::type;


template<class...us, class u>
struct substr_exponent<p<us...>, u>
{
    typedef p<
        typename substr_if_same_unit<us,u>::type ...
        > type;
};


template<class...us,class u, class... us2>
struct substr_exponent<p<us...>, u, us2...>
{
    typedef substr_exponent_t<substr_exponent_t<p<us...>,u>,us2...> type;
};




template<class ...us>
using p_t=simplify_t<add_exponent_t<
    typename my_unit_sytem::scalar_unit,us...>>
                                    ;

typedef simplify_t<my_scalar_unit> dimension_less;


typedef p_t<u<s,1>,u<m,-1>> ms_1;


template<class ...> struct v;
template<class ...> struct logv;




template<class T,class unit> class v<T,unit>
{
    T value_;
public:
    v(T&& x,unit): value_{std::move(x)}{}
    v(T&& x): value_{std::move(x)}{}
    const T& value()const &{return value_;}
    T value()&& {return value_;}


};

template<class T,class unit> class logv<T,unit>
{
    T value_;
    std::size_t n_;
public:
    logv(T&& x,unit, std::size_t n): value_{std::move(x)}, n_{n}{}
    logv(T&& x, std::size_t n): value_{std::move(x)}, n_{n}{}
    const T& value()const &{return value_;}
    T value()&& {return value_;}
    std::size_t size()const {return n_;}
};




template <class T>
v(T&& x)->v<T,dimension_less>;


template<class... units1, class ...units2>
auto operator*(p<units1...>,p<units2...>)
{ return  simplify_t<add_exponent_t<my_scalar_unit,units1...,units2...>>{};
}

template<class... units1, class ...units2>
auto operator/(p<units1...>,p<units2...>)
{ return  simplify_t<substr_exponent_t<add_exponent_t<my_scalar_unit,units1...>,units2...>>{};
}


template <class T, class U,class unit1, class unit2>
auto operator*( const v<T,unit1>& one, const  v<U,unit2>& other)
{
    return v(one.value()*other.value(),unit1{}*unit2{});
}


template <class T, class U,class unit1, class unit2>
auto operator/( const v<T,unit1>& one, const  v<U,unit2>& other)
{
    return v(one.value()/other.value(),unit1{}/unit2{});
}




template <class T, class U,class unit1>
auto operator+( const v<T,unit1>& one, const  v<U,unit1>& other)
{
    return v(one.value()+other.value(),unit1{});
}

template <class T, class U,class unit1>
auto operator+( const logv<T,unit1>& one, const  logv<U,unit1>& other)
{
    return logv(one.value()+other.value(),unit1{},one.size()+other.size());
}



template <class T, class U,class unit1>
auto operator-( const v<T,unit1>& one, const  v<U,unit1>& other)
{
    return v(one.value()-other.value(),unit1{});
}


template <class T, class U,class unit1>
auto operator-( const logv<T,unit1>& one, const  logv<U,unit1>& other)
{
    return logv(one.value()-other.value(),unit1{},one.size()-other.size());
}

template <class T,class unit1>
auto sqr(const v<T,unit1>& x) { return x * x; }




#include <cmath>

inline v<double,dimension_less> exp(const v<double,dimension_less>& x)
{
    return v<double,dimension_less>(std::exp(x.value()));
}

template <class unit1>
auto sqrt(const v<double,unit1>& x) { return v<double,div_exponent_t<unit1>>(std::sqrt(x.value())); }

template <class unit1>
auto log(const v<double,unit1>& x) { return logv<double,unit1>(std::log(x.value()),1); }




#endif // MYUNITSYSTEM_H
