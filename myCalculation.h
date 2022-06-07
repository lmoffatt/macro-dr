#ifndef MYCALCULATION_H
#define MYCALCULATION_H
#include "mytypetraits.h"
#include <functional>
template<class...>
class myCalculation;

template<class...>
class stepCalculation;


template <class Calculation, class Input, class... Parameters>
class stepCalculation<Calculation,Input,Cs<Parameters...>>
{
private:
    std::tuple<Parameters...> par_;
    Calculation c_;
public:


};


template<class... Calculator, class Output,class Input>
class myCalculation<Cs<Calculator...>,Input, Output>
{
public:
    static const std::size_t N=sizeof... (Calculator);
private:

    template<>
    auto calculate(const Input& input, double error_rate, std::size_t i, std::index_sequence<>)
    {
        return std::invoke(std::get<0>(c_),input,error_rate);

    }
    template<std::size_t I,std::size_t... Is>
    auto calculate(const Input& input, double error_rate, std::size_t i, std::index_sequence<I,Is...>)const
    {
        if (i==I)
            return std::invoke(std::get<I>(c_),input,error_rate);
        else
            return calculate(input,error_rate,i,std::index_sequence<Is...>());
    }
    auto calculate(const Input& input, double error_rate, std::size_t i)
    {
        return calculate(input, error_rate,i,std::index_sequence_for<Calculator...>());
    }

    std::tuple<Calculator...> c_;
    std::size_t best_one(const std::array<std::pair<double,double>,N>& times, double error_rate)
    {
        auto p=times[0];
        std::size_t n=0;
        for (std::size_t i=1;i<N; ++i )
        {
            auto& p1=times[i];
            if ((p.first-p1.first)/(p.second-p1.second)<error_rate)
            {
                p=p1;
                n=i;
            }
        }
        return n;
    }


public:
    Output operator()(const Input& input,double error_rate)
    {
        auto error_times=std::apply([&input,&error_rate](auto&... x){ return std::array<std::pair<double,double>,N>(x.estimate(input,error_rate)...);});

        auto i=best_one(error_times,error_rate);
        return calculate(input,error_rate,i);

    }



};


#endif //
