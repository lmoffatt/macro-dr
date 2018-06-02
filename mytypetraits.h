#ifndef MYTYPETRAITS_H
#define MYTYPETRAITS_H
#include <type_traits>
#include <variant>
#include <set>
#include <map>
#include <vector>
#include <utility>
#include <string_view>
#include <memory>


template<std::size_t N>
struct my_static_string
{
private:
    std::array<char,N> c_;
public:


    constexpr my_static_string(const char (&c)[N]):
        c_{}
    {
        for (std::size_t i=0; i<N; ++i)
            c_[i]=c[i];
    }

    constexpr char operator[](std::size_t i)const {
        return c_[i];}


    template<std::size_t N0>
    constexpr my_static_string(my_static_string<N0> one, my_static_string<N-N0+1> two): c_{}
    {
        for (std::size_t i=0; i<N0-1; ++i)
            c_[i]=one[i];
        for (std::size_t i=N0-1; i<N-1; ++i)
            c_[i]=two[i+1-N0];

        c_[N-1]='\0';

    }

    constexpr const char *c_str()const{
        return &c_[0];
    }


    std::string str()const { return c_str();}


};

template <int N>
my_static_string(const char (&lit)[N])   // <- match this
-> my_static_string<N>;







template <std::size_t N1, std::size_t N2>
constexpr auto operator+(const my_static_string<N1>& s1,
                         const my_static_string<N2>& s2)
-> my_static_string<N1 + N2 - 1>
{
    return my_static_string<N1 + N2 - 1>(s1, s2);
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

template<class ...Ts>
constexpr inline bool has_this_type_v=has_this_type<Ts...>::value;




template<template<typename...> class C>
class Co
{

};

template <class C>
struct Constructor
{
    typedef C myClass;
};

template <class C>
struct Loader
{
    typedef C myClass;
};

template <class C>
struct Valuer
{
    typedef C myClass;
};






template <class...>
struct my_trait
{

};

template <>
struct my_trait<std::string>
{
    constexpr static auto className=my_static_string("string");
};

template <>
struct my_trait<std::size_t>
{
    constexpr static auto className=my_static_string("count");
};

template <>
struct my_trait<int>
{
    constexpr static auto className=my_static_string("integer");
};

template <>
struct my_trait<bool>
{
    constexpr static auto className=my_static_string("boolean");
};


template <>
struct my_trait<char>
{
    constexpr static auto className=my_static_string("char");
};

template <>
struct my_trait<double>
{
    constexpr static auto className=my_static_string("real");
};


template <typename T, typename K>
struct my_trait<std::pair<T,K>>
{
    constexpr static auto className=my_static_string("pair_")+my_trait<T>::className+my_static_string("_")+my_trait<K>::className;
};

template <typename T, typename K>
struct my_trait<std::map<T,K>>
{
    constexpr static auto className=my_static_string("map_")+my_trait<T>::className+my_static_string("_")+my_trait<K>::className;
};

template <typename T>
struct my_trait<std::set<T>>
{
    constexpr static auto className=my_static_string("set_")+my_trait<T>::className;
};

template <typename T>
struct my_trait<std::vector<T>>
{
    constexpr static auto className=my_static_string("vector_")+my_trait<T>::className;
};

template <class C>
struct my_trait<const C&>
{
    constexpr static auto className=my_static_string("const_")+my_trait<std::decay_t<C>>::className+my_static_string("_ref");
};


template <class C>
struct my_trait<C>
{
    constexpr static auto className=C::className;
};

template<>
struct my_trait<Cs<>>
{
    constexpr static auto className=my_static_string("");
};

template <class T, class... Ts>
struct my_trait<Cs<T, Ts...>>
{
    constexpr static auto className=my_static_string("_")+my_trait<T>::className+my_trait<Cs<Ts...>>::className;
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

template <typename T> class myOptional;



template <class> struct is_optional: public std::false_type{};

template <typename T>
struct is_optional<myOptional<T>>: public std::true_type{};


template <typename T>
struct is_optional<std::optional<T>>: public std::true_type{};
template <typename T>
struct is_optional<std::optional<T>&>: public std::true_type{};

template <class> struct is_optional_ref: public std::false_type{};

template <typename T>
struct is_optional_ref<std::optional<std::reference_wrapper<T>>>: public std::true_type{};


template <typename T>
using is_optional_v=typename is_optional<std::decay_t<T>>::value;

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

template <class T> struct contains_loader: public std::false_type{};
template <class C>
struct contains_loader<Loader<C>> : public std::true_type {};

template <class T> struct contains_valuer: public std::false_type{};
template <class C>
struct contains_valuer<Valuer<C>> : public std::true_type {};



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
struct is_write_Object : std::false_type { };

template <typename T>
struct is_write_Object<T,
        std::void_t<decltype(std::declval<T&>().write(std::declval<std::ostream&>()))
>>
   : std::true_type { };

template <typename T, typename = void>
struct is_read_Object : std::false_type { };

template <typename T>
struct is_read_Object<T,
        std::void_t<decltype(std::declval<T&>().read(std::declval<std::istream&>()))
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

template <typename>
struct is_unique_ptr : std::false_type { };

template <typename T>
struct is_unique_ptr<std::unique_ptr<T>>
        : std::true_type { };






template <class> struct is_pointer_to_const: public std::false_type{};

template <typename T>
struct is_pointer_to_const<T const *>: public std::true_type{};

template<typename,class=void>
struct is_variable_ref: public std::false_type{};

template <typename T>
struct is_variable_ref<T,
        std::enable_if_t<!std::is_const_v<T>&&(std::is_lvalue_reference_v<T>||std::is_pointer_v<T>), int>>: public std::true_type{};

template< typename T>
inline constexpr bool is_variable_ref_v=is_variable_ref<T>::value;


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





template<typename,class=void>
struct has_push_back: std::false_type{};

template<typename T>
struct has_push_back<T,
        std::void_t<decltype(std::declval<T>().push_back(std::declval<typename T::value_type>()))>> : std::true_type{};

template<typename T> inline static constexpr bool has_push_back_v=has_push_back<T>::value;



template<typename,class=void>
struct has_insert: std::false_type{};

template<typename T>
struct has_insert<T,
        std::void_t<decltype(std::declval<T>().insert(std::declval<typename T::value_type>()))>> : std::true_type{};


template<typename T> inline static constexpr bool has_insert_v=has_push_back<T>::value;




static_assert(has_insert<std::map<double,double>>::value,"");


static_assert (true,"" );

template<typename...>
struct included_types{};


template<typename...Ts>
using included_types_t=typename included_types<Ts...>::type;

template<class...C> struct arg_types{
    typedef void type;
};

template<class...C> using arg_types_t =typename arg_types<C...>::type;



struct elem_tag{};

struct value_tag{};

struct map_tag{};

struct pair_tag{};

struct object_tag{};

struct command_tag{};

struct push_back_tag{};

struct insert_tag{};

struct set_tag{};



template <class T>
struct my_tag
{
    typedef std::decay_t<T> dT;
    //    typedef typename T::_in_my_tag dff;
    typedef
    std::conditional_t<is_field_Object<dT>::value, object_tag,
    std::conditional_t<is_arg_Command<dT>::value, command_tag,
    std::conditional_t<is_pair<dT>::value, pair_tag,
    std::conditional_t<has_value_type_v<dT>,
    std::conditional_t<has_mapped_type_v<dT>,map_tag, value_tag>,
    elem_tag>>>> type;
    //  typedef typename type::_in_my_tag dfff;

};



template <class T>
struct my_tag_arg
{
    typedef std::decay_t<T> dT;
    //    typedef typename T::_in_my_tag dff;
    typedef
    std::conditional_t<is_field_Object<dT>::value, std::pair<object_tag,arg_types_t<dT>>,
    std::conditional_t<is_arg_Command<dT>::value, command_tag,
    std::conditional_t<is_pair<dT>::value, pair_tag,
    std::conditional_t<has_push_back<dT>::value,push_back_tag,
    std::conditional_t<has_insert<dT>::value,
    std::conditional_t<has_mapped_type_v<dT>,map_tag,set_tag>,
    elem_tag>>>>> type;
    //  typedef typename type::_in_my_tag dfff;

};



template <class T> using my_tag_t=typename my_tag<T>::type;

template <class T> using my_tag_arg_t=typename my_tag_arg<T>::type;









template <class T>
struct included_types<T,elem_tag>
{
    typedef Cs<T> type;
    //typedef typename type::elem_type test_type;
};

template <class T>
struct included_types<T,pair_tag>
{
    typedef Cs<included_types_t<typename std::decay_t<T>::first_type>,included_types_t<typename std::decay_t<T>::second_type>,T> type;
    //  typedef typename type::pair_tag test_type;

};


template <typename T>
struct included_types<T,value_tag>
{
    typedef Cs<included_types_t<typename std::decay_t<T>::value_type>,T> type;
    // typedef typename type::value_tag test_type;

};

template <typename T>
struct included_types<T,map_tag>
{
    typedef typename std::decay_t<T>::key_type key_type;
    typedef typename std::decay_t<T>::mapped_type mapped_type;
    typedef std::pair<key_type,mapped_type> value_type;


    typedef Cs<included_types_t<value_type>,T> type;
    //typedef typename type::map_tag test_type;

};




//static_assert (has_key_type_v<std::vector<double>> );

template<typename T>
struct included_types<T,object_tag>;

template<typename T>
struct included_types<T,command_tag>;

template<typename T>
struct arg_types<T,object_tag>;

template<typename T>
struct arg_types<T,command_tag>;



template< typename T>
struct included_types<T>
{
    typedef included_types_t<T,my_tag_t<T>> type;
};

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
struct class_arg_set_union{};

template <class A, class B>
using class_arg_set_union_t=typename class_arg_set_union<A,B>::type;


template <class A, class B>
struct class_set_union_impl{};
template <class A, class B>
using class_set_union_impl_t=typename class_set_union_impl<A,B>::type;


template <class A, class B>
struct class_arg_set_union_impl{};
template <class A, class B>
using class_arg_set_union_impl_t=typename class_arg_set_union_impl<A,B>::type;




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



template<bool,class...>
struct class_set_do_union{};

template <template<class...> class Cs, class T, class... Ts>
struct class_set_do_union<false,Cs<Ts...>,T>
{
    typedef  Cs<Ts...,T> type;
};
template <template<class...> class Cs, class T, class... Ts>
struct class_set_do_union<true,Cs<Ts...>,T>
{
    typedef  Cs<Ts...> type;
};




template <template<class...> class Cs, class... Ts>
struct class_set_union<Cs<Ts...>,Cs<>>
{
    typedef Cs<Ts...> type;
};


template <template<class...> class Cs, class T, class... Ts>
struct class_set_union<Cs<Ts...>,T>
{
    typedef typename class_set_do_union<(std::is_same_v<Ts,T > ||...||false),Cs<Ts...>,T>::type type;

};




template <template<class...> class Cs, class T0,class... T, class... Ts>
struct class_set_union<Cs<Ts...>,Cs<T0,T...>>
{
    typedef class_set_union_t<class_set_union_t<Cs<Ts...>,T0>,Cs<T...>> type;
};


template <template<class...> class Cs, class... T>
struct class_arg_set_union<Cs<T...>,Cs<>>
{
    typedef Cs<T...> type;
};



template <template<class...> class Cs, class T, class... Ts>
struct class_arg_set_union<Cs<Ts...>,T>
{
    typedef class_set_union_t<Cs<Ts...>,included_types_t<T>> type;
};





template <template<class...> class Cs, class T0,class... T, class... Ts>
struct class_arg_set_union<Cs<Ts...>,Cs<T0,T...>>
{
    typedef class_arg_set_union_t<class_arg_set_union_t<Cs<Ts...>,T0>,Cs<T...>> type;
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


template<template<class>class,class ...> struct apply_Op_t_to_all{};




template<template<class>class Op_t, typename T>
struct apply_Op_t_to_all<Op_t,T>
{

    typedef Op_t<T> type;
};

template<template<class>class Op_t,template<class...>class Cs, typename...Ts>
struct apply_Op_t_to_all<Op_t,Cs<Ts...>>
{
    //static typename Cs<Ts...>::hurra k;

    typedef Cs<typename apply_Op_t_to_all<Op_t,Ts>::type...> type;
};






//class_set_union_t<Cs<double>,Cs<std::vector<std::map<double,int>>,std::vector<int>>>::kl a;




#endif // MYTYPETRAITS_H

