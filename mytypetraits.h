#ifndef MYTYPETRAITS_H
#define MYTYPETRAITS_H
#include <type_traits>
#include <variant>
#include <set>
#include <map>
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

template <class...>
struct my_trait
{

};

template <>
struct my_trait<std::size_t>
{
    constexpr static const char* name="count";
};

template <class C>
struct my_trait<C>
{
    constexpr static const char* name=C::classname;
};



template <class T> struct is_container: public std::false_type{};

template<template<typename T,typename> class Co,typename T, typename Alloc>
struct is_container<Co<T,Alloc>>: public std::true_type{};

template <typename T, typename = void>
struct is_std_container : std::false_type { };

template <typename T>
struct is_std_container<T,
    std::void_t<decltype(std::declval<T&>().begin()),
           decltype(std::declval<T&>().end()),
           typename T::value_type
           >>
    : std::true_type { };

template <class T> struct is_map: public std::false_type{};

template <class K, class T, class Comp, class Alloc>
struct is_map<std::map<K, T, Comp, Alloc>> : public std::true_type {};

template <class K, class T, class Comp, class Alloc>
struct is_map<std::multimap<K, T, Comp, Alloc>> : public std::true_type {};


template <class T> struct is_set: public std::false_type{};

template <class T, class Comp, class Alloc>
struct is_set<std::set< T, Comp, Alloc>> : public std::true_type {};


template <class T> struct is_pair: public std::false_type{};

template <class T, class K>
struct is_pair<std::pair< T, K>> : public std::true_type {};



template <class> struct is_variant: public std::false_type{};

template<typename... T>
struct is_variant<std::variant<T...>>
        : std::true_type { };

template <typename T, typename = void>
struct is_field_Object : std::false_type { };

template <typename T>
struct is_field_Object<T,
    std::void_t<decltype(std::declval<T&>().get_constructor_fields())
           >>
    : std::true_type { };



template <typename>
struct is_tuple : std::false_type { };

template <typename... T>
struct is_tuple<std::tuple<T...>>
    : std::true_type { };



#endif // MYTYPETRAITS_H
