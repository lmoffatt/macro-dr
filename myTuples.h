#ifndef MYTUPLES_H
#define MYTUPLES_H
#include <tuple>
#include <utility>
#include <type_traits>
#include <optional>


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



#endif // MYTUPLES_H
