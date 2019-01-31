#ifndef MYTUPLES_H
#define MYTUPLES_H
#include <tuple>
#include <utility>
#include <type_traits>
#include <optional>
#include <functional>


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
}

template<class F, class... Args, class... Args2, std::size_t ...Is>
auto
myApply(const F &f, std::tuple<Args...> &&t1, std::tuple<Args2...>&& t2, std::index_sequence<Is...>) {

  return std::invoke(f,std::make_pair(std::get<Is>(t1),std::get<Is>(t2))...);
}

  template <class F, class...Args, class...Args2>
auto myApply(const F& f, std::tuple<Args...>&& t1, std::tuple<Args2...> t2)
{
    static_assert(sizeof... (Args)==sizeof... (Args2) );
  return myApply(f,std::move(t1),std::move(t2),std::index_sequence_for<Args...>());

}







#endif // MYTUPLES_H
