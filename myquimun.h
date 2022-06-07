#ifndef MYQUIMUN_H
#define MYQUIMUN_H

#include "myunitsystem.h"
#include "myDistributions_units.h"
template <template<class...>class, class...> class Datum;
template <template<class...>class, class...> class Function;

template <template<class...>class, class...> class Measurement;
template <template<class...>class, class...> class Manipulation;


template<template<class...> class Tr,class Id> class Datum<Tr,Id>
{
public:

    typedef typename Tr<Id>::type T;
    typedef typename Tr<Id>::unit Unit;


private:
    v<T,Unit> value_;
public:
    Datum(T&& x):value_{std::move(x)}{}
    auto& value()const & {return value_;}
    v<T,Unit> value()&& {return value_;}
};

template<template<class...> class Tr,class F, class Id, class... Args> class Function<Tr,Id,F,Args...>
{
private:
    F f_;
public:
    Function(F&& x,Id,Args...):f_{std::move(x)}{}
    template<class Q>
    Datum<Tr,Id> operator()(const Q& q)const & {return std::invoke(f_,q(Args{})...);}
};





template<template<class...> class Tr,class Id> class Measurement<Tr,Id>: public Datum<Tr,Id>{
    using Datum<Tr,Id>::value;
    using Datum<Tr,Id>::Datum;

};

template<template<class...> class Tr,class Id> class Manipulation<Tr,Id>: public Datum<Tr,Id>{
    using Datum<Tr,Id>::value;
    using Datum<Tr,Id>::Datum;

};






template<class Id,class T,class unit> class EvenCoordinate
{
private:
    T start_;
    T step_;
    T end_;
    std::size_t nsteps_;
    static T get_end(T start,T step,std::size_t nsteps){ return start+nsteps*step;}
    static T get_step(T start, T end, std::size_t nsteps) { return (end-start)/nsteps;}
public:
    EvenCoordinate(v<T,unit> start, v<T,unit> end, std::size_t nsteps)
        :start_{start.value()},end_{end.value()},nsteps_{nsteps},step_{get_step(start.value(),end.value(),nsteps)} {}
    v<T,unit> start()const {return start_;}
    v<T,unit> step()const {return step_;}
    v<T,unit> end()const {return end_;}

    auto nsteps()const {return nsteps_;}
    v<T,unit> operator()(std::size_t i)const {return start_+i*step_;}

    std::size_t index(const v<T,unit>& x)
    {
        return (x.value()-start_)/step_;
    }

};



template<class Id,class T,class unit> class UnEvenCoordinate
{
private:private:

    std::map<T, std::size_t> m_;
    std::vector<T> v_;

    std::map<T, std::size_t> vector_to_map(const std::vector<T>& v)
    {
        std::map<T, std::size_t> out;
        for (std::size_t i=0; i<v.size(); ++i)
            out[v[i]]=i;
        return out;
    }
  public:
      UnEvenCoordinate(v<std::vector<T>,unit>&& c): m_{vector_to_map(c)},v_{std::move(c).value()}{}
    v<T,unit> start()const {return v_[0];}
    v<T,unit> end()const {return v_.back();}

    auto nsteps()const {return v_.size();}
    v<T,unit> operator()(std::size_t i)const {return v_[i];}

    std::size_t index(const v<T,unit>& x)
    {
        auto it= m_.lower_bound(x);
        return it->second;
    }

};


























#endif // MYQUIMUN_H
