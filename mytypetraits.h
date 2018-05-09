#ifndef MYTYPETRAITS_H
#define MYTYPETRAITS_H
#include <type_traits>
#include <variant>
#include <set>
#include <map>
#include <vector>
#include <utility>
#include <string_view>


template<std::size_t N>
struct my_static_string
{
private:
   std::array<char,N+1> c_;
 public:


     constexpr my_static_string(const char (&c)[N]):
         c_{}
    {
        for (std::size_t i=0; i<N; ++i)
            c_[i]=c[i];
        c_[N]='\0';
    }

     constexpr char operator[](std::size_t i)const {
         return c_[i];}


    template<std::size_t N0>
    constexpr my_static_string(my_static_string<N0> one, my_static_string<N-N0> two): c_{}
    {
        for (std::size_t i=0; i<N0; ++i)
            c_[i]=one[i];
        for (std::size_t i=N0; i<N; ++i)
            c_[i]=two[i-N0];

        c_[N]='\0';

    }

    constexpr const char *c_str()const{
        return &c_[0];
    }




};

template <int N>
 my_static_string(const char (&lit)[N])   // <- match this
  -> my_static_string<N>;







template <std::size_t N1, std::size_t N2>
constexpr auto operator+(const my_static_string<N1>& s1,
                         const my_static_string<N2>& s2)
-> my_static_string<N1 + N2>
{
return my_static_string<N1 + N2>(s1, s2);
}




template<typename T>
class C
{

};

template<typename... Ts>
struct Cs
{
};

template<class...> struct has_this_type{};

template< template<typename...> class Cs, typename T, typename... Ts>
struct has_this_type<Cs<Ts...>,T>
{
  static constexpr bool value=(std::is_same_v<T,Ts >||...||false);

};





template<template<typename...> class C>
class Co
{

};

template <class C>
struct Constructor
{
    typedef C myClass;
};




template <class...>
struct my_trait
{

};

template <>
struct my_trait<std::size_t>
{
    constexpr static auto name=my_static_string("count");
};

template <>
struct my_trait<double>
{
    constexpr static auto name=my_static_string("real");
};


template <class C>
struct my_trait<C>
{
     constexpr static auto name=C::name;
};


template <typename T, typename = void>
struct has_value_type : std::false_type { };

template <typename T>
struct has_value_type<T,
    std::void_t<typename T::value_type>>
    : std::true_type { };

template <typename T> inline constexpr bool has_value_type_v=has_value_type<T>::value;


template <typename T, typename = void>
struct has_mapped_type : std::false_type { };

template <typename T>
struct has_mapped_type<T,
    std::void_t<typename T::key_type, typename T::mapped_type>>
    : std::true_type { };

template <typename T> inline constexpr bool has_mapped_type_v=has_mapped_type<T>::value;




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

template <class T> struct contains_constructor: public std::false_type{};
template <class C>
struct contains_constructor<Constructor<C>> : public std::true_type {};



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


template <typename T, typename = void>
struct is_arg_Command : std::false_type { };

template <typename T>
struct is_arg_Command<T,
    std::void_t<decltype(T::get_arguments())
           >>
    : std::true_type { };



template <typename>
struct is_tuple : std::false_type { };

template <typename... T>
struct is_tuple<std::tuple<T...>>
    : std::true_type { };






template<typename,class=void>
struct has_get_global: std::false_type{};

template<typename T>
struct has_get_global<T,
        std::void_t<decltype(get(std::declval<std::stringstream&>(),
                                 std::declval<T&>()))>> : std::true_type{};



template<typename,class=void>
struct has_get_method: std::false_type{};

template<typename T>
struct has_get_method<T,
        std::void_t<decltype(std::declval<T&>().get(std::declval<std::stringstream&>()))>> : std::true_type{};


template<typename,class=void>
struct has_global_extractor: std::false_type{};

template<typename T>
struct has_global_extractor<T,
        std::void_t<decltype(operator<<(std::declval<std::stringstream&>(),std::declval<T&>()))>> : std::true_type{};



template<typename...>
struct arg_types{};


template<typename...Ts>
using arg_types_t=typename arg_types<Ts...>::type;

struct elem_tag{};

struct value_tag{};

struct map_tag{};

struct pair_tag{};

struct object_tag{};

struct command_tag{};


template <class T>
struct my_tag
{
    typedef
    std::conditional_t<is_field_Object<T>::value, object_tag,
    std::conditional_t<is_arg_Command<T>::value, command_tag,
    std::conditional_t<is_pair<T>::value, pair_tag,
    std::conditional_t<has_value_type_v<T>,
           std::conditional_t<has_mapped_type_v<T>,map_tag, value_tag>,
     elem_tag>>>> type;
};

template <class T> using my_tag_t=typename my_tag<T>::type;










template <class T>
struct arg_types<T,elem_tag>
{
    typedef Cs<T> type;
};

template <class T>
struct arg_types<T,pair_tag>
{
    typedef Cs<arg_types_t<typename T::first_type>,arg_types_t<typename T::second_type>,T> type;
};


template <typename T>
struct arg_types<T,value_tag>
{
    typedef Cs<arg_types_t<typename T::value_type>,T> type;
};

template <typename T>
struct arg_types<T,map_tag>
{
    typedef Cs<arg_types_t<typename T::key_type>,arg_types_t<typename T::mapped_type> ,T> type;
};


//static_assert (has_key_type_v<std::vector<double>> );

template<typename T>
struct arg_types<T,object_tag>;

template<typename T>
struct arg_types<T,command_tag>;


template< typename T>
struct arg_types<T>
{
   typedef arg_types_t<T,my_tag_t<T>> type;
};


template<class ...>struct class_concatenate{};


template <class ...T>
using class_concatenate_t=typename class_concatenate<T...>::type;

template<template<class...>class Op, class... T0, template<class...> class Cs, class... T>
struct class_concatenate<Op<T0...>,Cs<T...>>
{
    typedef Op<T0...,T...> type;
};



template <class A, class B>
struct class_set_union{};

 template <class A, class B>
using class_set_union_t=typename class_set_union<A,B>::type;


template <class A, class B>
struct class_set_union_impl{};
template <class A, class B>
using class_set_union_impl_t=typename class_set_union_impl<A,B>::type;




template <template<class...> class Cs, class... Ts>
struct class_set_union_impl<Cs<Ts...>,Cs<>>
{
    typedef Cs<Ts...> type;
};


template <template<class...> class Cs, class T0,class... T, class... Ts>
struct class_set_union_impl<Cs<Ts...>,Cs<T0,T...>>
{
    typedef class_set_union_impl_t<class_set_union_impl_t<Cs<Ts...>,T0>,Cs<T...>> type;
};

template <template<class...> class Cs, class T, class... Ts>
struct class_set_union_impl<Cs<Ts...>,T>
{
    typedef std::conditional_t<(std::is_same_v<Ts,T > ||...||false),Cs<Ts...>,Cs<Ts...,T>> type;

};

template <template<class...> class Cs, class T, class... Ts>
struct class_set_union<Cs<Ts...>,T>
{
    typedef class_set_union_impl_t<Cs<Ts...>,arg_types_t<T>> type;
};

template <template<class...> class Cs, class T0,class... T, class... Ts>
struct class_set_union<Cs<Ts...>,Cs<T0,T...>>
{
    typedef class_set_union_t<class_set_union_t<Cs<Ts...>,T0>,Cs<T...>> type;
};



template<class...> struct remove_void{};

template<class...C>
using remove_void_t=typename remove_void<C...>::type;


template<template <class...>class Cs,typename ...Ts>
struct remove_void<Cs<Ts...>,Cs<>>
{
  typedef Cs<Ts...> type;
};



template<template<class...>class Cs,typename ...Ts, typename...T>
struct remove_void<Cs<Ts...>,Cs<void,T...>>
{
  typedef remove_void_t<Cs<Ts...>,Cs<T...>> type;
};


template<template<class...>class Cs,typename ...Ts, typename T0,typename...T>
struct remove_void<Cs<Ts...>,Cs<T0,T...>>
{
  typedef remove_void_t<Cs<Ts...,T0>,Cs<T...>> type;
};




template<template<class...>class Cs, typename...T>
struct remove_void<Cs<T...>>
{
  typedef remove_void_t<Cs<>,Cs<T...>> type;
};





template<class...> struct filter_pointer{};

template<class...C>
using filter_pointer_t=typename filter_pointer<C...>::type;


template <template<class...> class Cs,class...T0>
struct filter_pointer<Cs<T0...>>
{
    typedef remove_void_t<Cs<std::conditional_t<std::is_pointer_v<T0>,T0,void>...>>  type;
};






template<class...> struct filter_regular{};

template<class...C>
using filter_regular_t=typename filter_regular<C...>::type;


template <template<class...> class Cs,class...T0>
struct filter_regular<Cs<T0...>>
{
    typedef remove_void_t<Cs<std::conditional_t<!std::is_pointer_v<T0>,T0,void>...>>  type;
};




//class_set_union_t<Cs<double>,Cs<std::vector<std::map<double,int>>,std::vector<int>>>::kl a;




#endif // MYTYPETRAITS_H

