#ifndef MYFIELDS_H
#define MYFIELDS_H

#include "mytypetraits.h"

#include <string>
#include <optional>
#include <type_traits>
namespace grammar {



template <typename T>
struct argument
{
    typedef T  argument_type;
    typedef std::optional<std::decay_t<argument_type>> default_type;
    std::string idField;
    default_type default_value;
    argument(C<T>,const char* id):idField{id}, default_value{}{}
    argument(C<T>,const char* id, T val):idField{id}, default_value{std::move(val)}{}
};



}
template <auto>
struct function_trait{};





namespace grammar {


template <class Object,class Method>
struct field
{
    typedef
    std::invoke_result_t<Method,Object> return_type;

    typedef Method member_type;
    typedef Object object_type;
    typedef  std::decay_t<return_type> result_type;
    typedef std::optional<std::decay_t<return_type>> default_type;

    std::string idField;
    member_type access_method;
    default_type default_value;
    field(C<Object>,const char* id, member_type get):idField{id},access_method{get}, default_value{}{}
};



template <class Object,class... Method>
bool has_all(const std::tuple<field<Object,Method>...>& fs)
{
    return std::apply([](auto&...x){return (x.default_value.has_value()&&...&&true);},fs);
}

template <class... Field>
auto getIdFields(const std::tuple<Field...>& fs)
{
    return std::apply([](auto&...x){return std::set<std::string>{x.idField...};},fs);
}





} // namespace grammar


template <class Object,class... Method>
struct arg_types<std::tuple<grammar::field<Object,Method>...>>
{
    typedef Cs<std::invoke_result_t<Method,Object> ...>  type;
};



template <class... T>
struct arg_types<std::tuple<grammar::argument<T>...>>
{
    typedef Cs<T...> type;
};


template <class C>
struct arg_types<C,object_tag>
{
    typedef arg_types_t<decltype (C::get_constructor_fields())> type;
};




template <class Object,class... Method>
struct included_types<std::tuple<grammar::field<Object,Method>...>>
{
    typedef Cs<included_types_t<std::invoke_result_t<Method,Object>> ...>  type;
};


template <class... T>
struct included_types<std::tuple<grammar::argument<T>...>>
{
    typedef Cs<included_types_t<T>...> type;
};



template<typename...> struct extract_all_types_and_decay{};


template<class command>
struct result_of_command
{

    typedef  typename class_concatenate_t<
    std::invoke_result<decltype(&command::run)>,
    arg_types_t<std::invoke_result_t<decltype(command::get_arguments)>>
    >::type type;
  //  static typename type::chota d;

};

template<class command>
using result_of_command_t=typename result_of_command<command>::type;



template<typename T>
struct included_types<T,object_tag>
{
    typedef Cs<T,included_types_t<std::invoke_result_t<decltype(std::decay_t<T>::get_constructor_fields)> >> type;
    static std::pair<std::string,std::set<std::string>> getIdArgs()
    {
        return {my_trait<T>::className.c_str(),getIdFields(std::decay_t<T>::get_constructor_fields())};
    }
    //typedef typename type::object_tag test_type;


};



template<typename T>
struct included_types<T,command_tag>
{
    typedef included_types_t<std::invoke_result_t<decltype(T::get_arguments)> > type;
    static std::pair<std::string,std::set<std::string>> getIdArgs()
    {
        return {T::className.c_str(),getIdFields(T::get_arguments())};
    }
   // typedef typename type::command_tag test_type;

};





template <template <class...>class Cs,class ...Ts>
struct extract_all_types_and_decay<Cs<Ts...>>
{
    typedef class_set_union_t<Cs<>,typename apply_Op_t_to_all<std::decay_t,Cs<included_types_t<Ts>...>>::type> type;
};
template <class T>
using extract_all_types_and_decay_t=typename  extract_all_types_and_decay<T>::type;



template< class...>
struct extract_function_return_types{};

template< class...>
struct extract_types{};


template<template<class...> class Cs, class... types,class... commands>
struct extract_function_return_types<Cs<types...>,Cs<commands...>>
{
    typedef     class_set_union_t<class_set_union_t<Cs<>,remove_void_t<Cs<std::conditional_t<is_field_Object<types>::value,types,void>...>>>,Cs<result_of_command_t<commands>...>> type;

};




template<template<class...> class Cs,class... types, class... commands>
struct extract_types<Cs<types...>,Cs<commands...>>
{

    typedef  class_set_union_t<
    extract_all_types_and_decay_t<Cs<types...>>,extract_all_types_and_decay_t<Cs<commands...>>> type;

};






#endif // MYFIELDS_H
