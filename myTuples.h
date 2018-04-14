#ifndef MYTUPLES_H
#define MYTUPLES_H
#include <tuple>
#include <utility>

// ------------- UTILITY---------------
template<int...> struct index_tuple{};

template<int I, typename IndexTuple, typename... Types>
struct make_indexes_impl;

template<int I, int... Indexes, typename T, typename ... Types>
struct make_indexes_impl<I, index_tuple<Indexes...>, T, Types...>
{
    typedef typename make_indexes_impl<I + 1, index_tuple<Indexes..., I>, Types...>::type type;
};

template<int I, int... Indexes>
struct make_indexes_impl<I, index_tuple<Indexes...> >
{
    typedef index_tuple<Indexes...> type;
};

template<typename ... Types>
struct make_indexes : make_indexes_impl<0, index_tuple<>, Types...>
{};

// ----------- FOR EACH -----------------
template<typename Func, typename Last>
void for_each_impl(Func&& f, Last&& last)
{
    f(last);
}

template<typename Func, typename First, typename ... Rest>
void for_each_impl(Func&& f, First&& first, Rest&&...rest)
{
    f(first);
    for_each_impl( std::forward<Func>(f), rest...);
}

template<typename Func, int ... Indexes, typename ... Args>
void for_each_helper( Func&& f, index_tuple<Indexes...>, std::tuple<Args...>&& tup)
{
    for_each_impl( std::forward<Func>(f), std::forward<Args>(std::get<Indexes>(tup))...);
}

template<typename Func, typename ... Args>
void for_each( std::tuple<Args...>& tup, Func&& f)
{
   for_each_helper(std::forward<Func>(f),
                   typename make_indexes<Args...>::type(),
                   std::forward<std::tuple<Args...>>(tup) );
}

template<typename Func, typename ... Args>
void for_each( std::tuple<Args...>&& tup, Func&& f)
{
   for_each_helper(std::forward<Func>(f),
                   typename make_indexes<Args...>::type(),
                   std::forward<std::tuple<Args...>>(tup) );
}




template<typename T>
class C
{

};

template<typename... T>
class Cs
{

};

template<template<typename...> class C>
class Co
{

};



template<typename...>
struct is_container
{
  static constexpr bool value=false;
};

template<template<typename T,typename> class Co,typename T, typename Alloc>
struct is_container<Co<T,Alloc>>
{
  static constexpr bool value=true;
};

template<typename T> void hash_combine(std::size_t & seed, T const& v)
{
  seed ^= std::hash<T>(v)() + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template<typename T, typename...Ts>
std::size_t get_hash_impl(std::size_t seed)
{
   return seed;
}

template<typename T, typename...Ts>
std::size_t get_hash_impl(std::size_t seed,T const & v,Ts const &... ts)
{
   hash_combine(seed,v);
   return get_hash_impl(seed,ts...);
}

template<typename T, typename...Ts>
std::size_t get_hash(T const & v,Ts const &... ts)
{
   std::size_t h=std::hash<T>(v)();
   return get_hash_impl(h,ts...);
}

namespace std
{
template<template<typename T,typename> class Co,typename T, typename Alloc>
struct hash<Co<T,Alloc>>
{
   std::size_t operator()(const Co<T,Alloc>& v)
   {
      std::size_t s=std::hash<std::size_t>()(0);
      for (auto e:v)
        hash_combine(s,e);
      return s;
   }
};

};







#endif // MYTUPLES_H
