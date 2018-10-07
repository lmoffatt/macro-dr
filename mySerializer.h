#ifndef MYSERIALIZER_H
#define MYSERIALIZER_H

#include <set>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

#include "myTuples.h"
#include "myfields.h"
#include "mytypetraits.h"

#include <cctype>
#include <variant>
#include <functional>
#include <type_traits>
#include <sstream>
#include <string>
//template<typename T>
//  std::ostream& write(std::ostream& os, const T& e);


///----------------- Declarations-----------------------------------------///

/*-     POD  */


namespace io {


std::string my_to_string(double x)
{
    std::stringstream strstr;
    if ((std::abs(x)<1e-3)||(std::abs(x)>1e3))
    {
        strstr << std::scientific << x;
    }
    else
        strstr  << x;

    return strstr.str();
}


template <typename T>
auto  ToString(const T& x)->decltype ((std::declval<std::stringstream&>()<<x), std::string())
{std::stringstream ss; ss<<x; return ss.str();}



template <typename T>
auto  ToString(const T& x)->decltype (x.put(std::declval<std::ostream&>()), std::string())
{std::stringstream ss; x.put(ss); return ss.str();}

template <>
std::string  ToString(const double& x)
{std::stringstream ss; ss<<std::fixed; ss<<x; return ss.str();}


template<typename T> std::ostream& write(std::ostream& s, const T& v);
template<typename T> std::istream& read(std::istream& s, T& v);


inline std::istream& read(std::istream& is, double& e);

inline std::ostream& write(std::ostream& os, const std::string& s);
inline std::istream& read(std::istream& is,  std::string& s);


template<typename T> std::istream& read(std::istream& is, T*& e);

template<typename T> std::istream& read_optional(std::istream& is, myOptional_t<T>& e);


inline std::istream &safeGetline(std::istream &is, std::string &t);


std::istream& read(std::istream& is, double& e);

template<typename ...Args> std::ostream& operator<<(std::ostream& os, const std::tuple<Args...>& tu);
template<typename ...Args> std::ostream& write_tuple(std::ostream& os, const std::tuple<Args...>& tu);
template<typename ...Args> std::ostream& write_tuple_i(std::ostream& os, const std::tuple<Args...>& tu, std::size_t i);

template<typename ...Args> std::istream& operator>>(std::istream& is,  std::tuple<Args...>& tu);
template<typename ...Args> std::istream& read_tuple(std::istream& is,  std::tuple<Args...>& tu);


template<typename T1, typename T2> std::ostream& operator<<(std::ostream& os,const std::pair<T1,T2>& other);
template<typename T1, typename T2> std::istream& operator>>(std::istream& is,const std::pair<T1,T2>& other);
template<typename T1, typename T2> std::ostream& write_pair(std::ostream& os,const std::pair<T1,T2>& other);
template<typename T1, typename T2> std::istream& read_pair(std::istream& is, std::pair<T1,T2>& other);


template<typename T> auto operator<<(std::ostream& s, const T& v)->std::enable_if_t<is_std_container<T>::value&&!std::is_same_v<T,std::string >,std::ostream&>;


template<typename T> auto operator>>(std::istream& s, T& v)
->std::enable_if_t<is_std_container<T>::value&&!std::is_same_v<T,std::string >, std::istream&>;
template<typename T> std::istream& read_vector(std::istream& s, std::vector<T>& v);

template<typename K,typename T> std::ostream& operator<<(std::ostream& s, const std::map<K,T>& v);

template<typename K,typename T> std::istream& operator>>(std::istream& s, std::map<K,T>& v);
template<typename K,typename T> std::ostream& write_map(std::ostream& s, const std::map<K,T>& v);

template<typename K,typename T> std::istream& read_map(std::istream& s, std::map<K,T>& v);


template<typename K,typename T> std::istream& operator>>(std::istream& s, std::multimap<K,T>& v);
template<typename K,typename T> std::ostream& operator<<(std::ostream& s, const std::multimap<K,T>& v);

template<typename K,typename T> std::istream& read_multimap(std::istream& s, std::multimap<K,T>& v);

template<typename K,typename T> std::ostream& write_multimap(std::ostream& s, const std::multimap<K,T>& v);
template<typename K> std::ostream& operator<<(std::ostream& s, const std::set<K>& v);
template<typename K> std::istream& operator>>(std::istream& s, const std::set<K>& v);

template<typename K> std::ostream& write_set(std::ostream& s, const std::set<K>& v);
template<typename K> std::istream& read_set(std::istream& s,  std::set<K>& v);

template<typename T> auto operator<<(std::ostream& s, T const& v)->std::enable_if_t<is_field_Object<T>::value,std::ostream&>;
template<typename T> auto operator<<(std::ostream& s, T const& v)->std::enable_if_t<is_write_Object<T>::value,std::ostream&>;

template<typename T> auto operator>>(std::istream& s, T & v)->std::enable_if_t<is_read_Object<T>::value,std::istream&>;

template<typename T> auto operator>>(std::istream& s, T*& v)
->std::enable_if_t<std::is_abstract<T>::value,std::istream&>;

template<typename T> auto operator>>(std::istream& s, T& v)->std::enable_if_t<is_field_Object<T>::value,std::istream&>;
template<typename ...Args> std::istream& operator>>(std::istream& is,  std::tuple<Args...>& tu);




template<typename Token,typename T> std::istream& read_variant(std::istream& is, std::variant<Token,T>& e);
template<typename Token,typename T> std::istream& input_operator_variant(std::istream& is, std::variant<Token,T>& e);



template<char... cs> struct token {
    static constexpr char value[]={cs...};
    typedef int this_is_token;
};
template <typename >
struct is_token : std::false_type { };

template <char...c>
struct is_token<token<c...>>
        : std::true_type { };

template<template<char...>class token,char... c>
std::ostream& operator<<(std::ostream& os, token<c...> );



template<template<char...>class token,char c, char... chs>
std::istream& operator>>(std::istream& is, token<c,chs...> );




template<template<char...>class token>
std::istream& operator>>(std::istream& is, token<> );

template<template<char...>class token,char c, char... chs>
std::istream& operator>>(std::istream& is, token<c,chs...> );

template<template<char...>class token,char c, char... chs>
bool get(std::stringstream& ss, token<c,chs...> t);


typedef token<'{'> start_of_Container;

typedef token<'}'> end_of_Container;

typedef token<'\n'> end_of_line;


typedef start_of_Container start_of_tuple;

typedef end_of_Container end_of_tuple;

typedef start_of_Container start_of_pair;

typedef end_of_Container end_of_pair;


typedef token<'\t'> separator;

typedef token<'='> equal;

typedef separator write_separator;

typedef io::token<'('> arg_start;
typedef io::token<','> arg_separator;

typedef io::token<')'> arg_end;


struct size_of_container
{
    std::size_t size;
};

std::ostream& operator<<(std::ostream& os, const size_of_container& n);

std::istream& operator>>(std::istream& is, size_of_container& n);

template<typename T> std::ostream& output_operator_on_element(std::ostream& os, const T& e);

template<typename T> std::ostream& write_on_element(std::ostream& os, const T& e);

template<typename T> std::ostream& output_operator_on_element(std::ostream& os, const T *&  e);
template<typename T> std::ostream& write_on_element(std::ostream& os, const T*& e);
template<typename T> std::istream& input_operator_on_element(std::istream& is, T& e);
template<typename T> std::istream& input_operator_on_element(std::istream& is, T*& e);


template<typename T>
std::istream& read_on_element(std::istream& is, T& e);
template<typename T>
std::istream& read_on_element(std::istream& is, T*& e);
template<typename Token,typename T>
std::istream& read_on_element(std::istream& is, std::variant<Token,T>& e);


inline
std::istream& extract_infinite(std::istream& ss, double& val, bool is_negative);
inline
std::istream& extract_nan(std::istream& ss, double& val, bool is_negative);
inline
std::istream& extract_finite(std::istream& ss, double& val, bool is_negative);

inline
std::istream& extract_double(std::istream& ss, double& val);

std::istream& input_operator_on_element(std::istream& is, double& e);
std::istream& input_operator_on_element(std::istream& is, double*& e);


template <class Object,class Method> std::ostream& write_on_Field(std::ostream& os,const grammar::field<Object,Method>& myField, const Object& myObject);
template <class Object,class Method> std::ostream& output_operator_on_Field(std::ostream& os,const grammar::field<Object,Method>& myField, const Object& myObject);


template<class FieldObject> auto write_on_Object(std::ostream& os,const FieldObject& myObject)->decltype (myObject.get_constructor_fields(),os)&;
template<class FieldObject> auto output_operator_on_Object(std::ostream& os,const FieldObject& myObject)->decltype (myObject.get_constructor_fields(),os)&;

template<class AbstractObject> auto output_operator_on_Abstract_Object(std::ostream& os,const AbstractObject& myObject)->decltype (myObject.get_constructor_fields(),os)&;



template <class Object,class Method>
bool read_on_if_field(std::istream& is,const std::string id,  grammar::field<Object,Method>& myField);

template <class Object,class Method>
bool input_operator_on_if_field(std::istream& is,const std::string id,  grammar::field<Object,Method>& myField);


template <class Object,class... Method>
std::istream& read_on_this_field(std::istream& is,const std::string id,  std::tuple<grammar::field<Object,Method>...>& myFields);
template <class Object,class... Method>
std::istream& input_operator_on_this_field(std::istream& is,const std::string id,  std::tuple<grammar::field<Object,Method>...>& myFields);



template<class FieldObject>
auto read_on_Object(std::istream& is, FieldObject& myObject)->decltype (myObject.get_constructor_fields(),is)&;

template<class FieldObject>
auto input_operator_on_Object(std::istream& is, FieldObject& myObject)->decltype (myObject.get_constructor_fields(),is)&;

template<class AbstractObject>
auto input_operator_on_Abstract_Object(std::istream& is, AbstractObject*& myObject)->decltype (std::declval<AbstractObject>().get_constructor_fields(),is)&;


template<typename T1, typename T2> std::ostream& operator<<(std::ostream& os,const std::pair<T1,T2>& other);
template<typename T1, typename T2> std::istream& operator>>(std::istream& is,std::pair<T1,T2>& other);


template<class Container> std::ostream& output_operator_on_container(std::ostream& os,const Container& myContainer);



template<class Container> std::ostream& write_on_container(std::ostream& os,const Container& myContainer);

template<template <class,class> class Container,typename T, class Allocator> std::istream& input_operator_on_container(std::istream& is, Container<T,Allocator> & myContainer);
template<template <class,class> class Container,typename T, class Allocator> std::istream& read_on_container(std::istream& is, Container<T,Allocator> & myContainer);



template<template <class,class, class...> class Map,typename K,typename T, class... Os> std::istream& input_operator_on_Map(std::istream& is, Map<K,T,Os...> & myMap);
template<template <class,class, class...> class Map,typename K,typename T, class... Os> std::istream& read_on_Map(std::istream& is, Map<K,T,Os...> & myMap);


template<template <class, class...> class Set,typename K,class... Os> std::istream& input_operator_on_set(std::istream& is, Set<K,Os...> & mySet);

template<template <class, class...> class Set,typename K,class... Os> std::istream& read_on_set(std::istream& is, Set<K,Os...> & mySet);


///-----------------Implementations------------------------------------------///


inline std::ostream& write(std::ostream& os, const std::string& s) {return os<<s;}

inline std::istream& read(std::istream& is,  std::string& s){return is>>s;}



template<typename T> std::istream& read(std::istream& is, T*& e)
{
    auto o=new T{};
    read(is,*o);
    e=o;
    return is;
}

template<typename T> std::istream& read_optional(std::istream& is, myOptional_t<T>& e)
{
    T x;
    if (read(is,x))
    {
        e=x;
    }
    else
    {
        is.clear();
        e.reset();
    }
    return is;
}

template<typename T> std::istream& read_optional(std::istream& is, std::optional<T>& e)
{
    T x;
    if (read(is,x))
    {
        e=x;
    }
    else
    {
        is.clear();
        e.reset();
    }
    return is;
}


template<typename Token,typename T> std::istream& read_variant(std::istream& is, std::variant<Token,T>& e)
{
    Token t;
    if (read(is,t))
    {
        e=t;
    }
    else
    {
        is.clear();
        T x;
        if (read(is,x))
            e=x;
    }
    return is;
}


template<typename Token,typename T> std::istream& input_operator_variant(std::istream& is, std::variant<Token,T>& e)
{
    Token t;
    if (is>>t)
    {
        e=t;
    }
    else
    {
        is.clear();
        T x;
        if (is>>x)
            e=x;
    }
    return is;
}



inline std::istream &safeGetline(std::istream &is, std::string &t)
{
    is.clear();
    std::getline(is,t);
    auto it=t.find('\r');
    if (it!=t.npos)
        t.erase(it);
    return is;
}



template<char... c>
std::ostream& operator<<(std::ostream& os, token<c...> )
{
    (os<<...<<c);
    return os;
}


template<char c, char... chs>
std::istream& operator>>(std::istream& is, token<c,chs...> );


std::istream& operator>>(std::istream& is, token<> )
{
    return is;
}


std::istream& set_fail_bit(std::istream& is)
{
    is.setstate(std::ios::failbit);
    return is;
}

template<char c, char... chs>
std::istream& operator>>(std::istream& is, token<c,chs...> )
{
    if constexpr(c=='\n')
    {
        char ch;
        while(is.get(ch)&& (ch==' '|| ch=='\t'|| ch=='\r')){};
        if (ch!=c)
        {
            is.putback(ch);
            return set_fail_bit(is);
        }
    }
    else if constexpr (c!=' '&& c!='\t')
    {
        char ch;
        is>>ch;
        if (ch!=c)
        {
            is.putback(ch);
            is.setstate(std::ios::failbit);
            return is;
        }
    }
    return is>>token<chs...>{};
}

template<char c, char... chs>
bool get(std::stringstream& ss, token<c,chs...> t)
{
    auto pos=ss.tellg();
    if (ss>>t)
    {
        return true;
    }
    else
    {
        ss.clear();
        ss.seekg(pos);
        return false;
    }
}



typedef token<'{'> start_of_Container;

typedef token<'}'> end_of_Container;

typedef token<'\n'> end_of_line;

typedef token<'\\','\\','-','-','-','\n'> end_of_Object;


typedef start_of_Container start_of_tuple;

typedef end_of_Container end_of_tuple;

typedef token<'['> start_of_title;

typedef token<']'> end_of_title;

typedef start_of_Container start_of_row;

typedef end_of_Container end_of_row;

typedef token<'#'> end_of_data_frame;


typedef start_of_Container start_of_pair;

typedef end_of_Container end_of_pair;


typedef token<'\t'> separator;

typedef token<'='> equal;

typedef separator write_separator;



std::ostream& operator<<(std::ostream& os, const size_of_container& n)
{
    os<<"size"<<equal{}<<n.size;
    return os;
}

std::istream& operator>>(std::istream& is, size_of_container& n)
{
    token<'s','i','z','e'> t;
    is>>t>>equal{}>>n.size;
    return is;
}


template<typename T>
std::ostream& output_operator_on_element(std::ostream& os, const T& e)
{
    return os<<e;
}

template<typename T>
std::ostream& write_on_element(std::ostream& os, const T& e)
{
    write(os,e);
    return os;
}

template<typename T>
std::ostream& output_operator_on_element(std::ostream& os, const T * & e)
{
    return os<<*e;
}

template<typename T>
std::ostream& write_on_element(std::ostream& os, const T*& e)
{
    write(os,*e);
    return os;
}

template<typename T>
std::istream& input_operator_on_element(std::istream& is, T& e)
{
    return  is>>e;
}

template<typename T>
std::istream& input_operator_on_element(std::istream& is, T*& e)
{
    e=new T;
    return is>>*e;
}



template<typename T>
std::istream& read_on_element(std::istream& is, T& e)
{
    read(is,e);
    return is;
}

template<typename T>
std::istream& read_on_element(std::istream& is, T*& e)
{
    e=new T;
    read(is,*e);
    return is;
}


template<typename T>
std::istream& read_one_element_on_vector(std::istream& is, std::vector<T>& v)
{
    T e;
    if (!read(is,e))
        return is;
    else
    {
        v.emplace_back(std::move(e));
    }
    return is;
}

template<typename T>
std::istream& read_one_element_on_vector(std::istream& is, std::vector<T*>& v)
{
    T* e=new T;
    if (!read(is,*e))
        return is;
    else
    {
        v.emplace_back(std::move(e));
    }
}




template<typename Token,typename T>
std::istream& read_on_element(std::istream& is, std::variant<Token,T>& e)
{
    Token t;
    T x;
    if (!read(is,t))
    {
        is.clear();
        if (!read(is,x)) return is;
        else
        {
            e=x;
            return is;
        }
    }
    else{
        e=t;
        return is;

    }
    return is;
}



inline
std::istream& extract_infinite(std::istream& ss, double& val, bool is_negative)
{
    std::string s;
    ss>>s;
    if ((s=="inf")||(s=="infinity")||(s=="INF")||(s=="INFINITY"))
    {
        if (is_negative)
            val=-std::numeric_limits<double>::infinity();
        else
            val=std::numeric_limits<double>::infinity();
        return ss;
    }
    else
    {
        ss.setstate(std::ios::failbit);
        return ss;
    }
}

inline
std::istream& extract_nan(std::istream& ss, double& val, bool is_negative)
{
    std::string s;
    ss>>s;
    if ((s=="nan")||(s=="NAN")||(s=="NaN"))
    {
        if (is_negative)
            val=-std::numeric_limits<double>::quiet_NaN();
        else
            val=std::numeric_limits<double>::quiet_NaN();
        return ss;
    }
    else
    {
        ss.setstate(std::ios::failbit);
        return ss;
    }
}

inline
std::istream& extract_finite(std::istream& ss, double& val, bool is_negative)
{
    ss>>val;
    if (is_negative)
        val=-val;
    return ss;
}


inline
std::istream& extract_double(std::istream& ss, double& val)
{
    char c;
    bool is_negative=false;
    ss>>c;
    if (c=='-')
    {
        is_negative=true;
        ss>>c;
    }
    else if (c=='+')
        ss>>c;
    switch(c)
    {
    case '0':
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
    {
        ss.putback(c);
        return extract_finite(ss,val,is_negative);
    }

    case 'i':
    case 'I':
    {
        ss.putback(c);
        return extract_infinite(ss,val,is_negative);
    }
    case 'n':
    case 'N':
    {
        ss.putback(c);
        return extract_nan(ss,val,is_negative);
    default:
        {
            ss.setstate(std::ios::failbit);
            return ss;
        }
        }

    }

}


std::istream& input_operator_on_element(std::istream& is, double& e)
{
    return extract_double(is,e);
}

std::istream& input_operator_on_element(std::istream& is, double*& e)
{
    e=new double;
    return extract_double(is,*e);
}



std::istream& read(std::istream& is, double& e)
{
    return io::extract_double(is,e);
}

template<typename ...Args> std::ostream& operator<<(std::ostream& os, const std::tuple<Args...>& tu)
{
    os<<io::start_of_Container{};
    std::apply([&os](const auto&... v){return ((os<<v<<io::separator{}),...,0);},tu);
    os<<io::end_of_Container{};
    return os;
}

template<typename ...Args> std::ostream& write_tuple(std::ostream& os, const std::tuple<Args...>& tu)
{
    std::apply([&os](const auto&... v){return ((os<<v<<io::separator{}),...,0);},tu);
    os<<io::end_of_line{};
    return os;
}

template<std::size_t I,typename ...Args >
void write_tuple_i(std::ostream& os, const std::tuple<Args...>& tu, std::size_t i)
{
    if (i==I)
       os<<std::get<I>(tu);
}

template<std::size_t...Is,typename ...Args>
std::ostream& write_tuple_i(std::ostream& os, const std::tuple<Args...>& tu, std::size_t i, std::index_sequence<Is...>)
{
    (write_tuple_i<Is>(os,tu, i),...);
    return os;
}
template<typename ...Args> std::ostream& write_tuple_i(std::ostream& os, const std::tuple<Args...>& tu, std::size_t i)
{
   return write_tuple_i(os,tu,i,std::index_sequence_for<Args...>());
}


template<typename ...Args> std::istream& read_tuple(std::istream& is,  std::tuple<Args...>& tu)
{
    std::apply([&is](auto&... v){return ((is>>io::separator{}>>v),...,0);},tu);
    is>>io::end_of_line{};
    return is;
}




template <class Object,class Method>
std::ostream& write_on_Field(std::ostream& os,const grammar::field<Object,Method>& myField, const Object& myObject){
    io::write_on_element (os,myField.idField);
    os<<io::end_of_line{};
    write(os,std::invoke(myField.access_method,myObject));
    os<<io::end_of_line{};
    return os;
}


template <class Object,class Method>
std::ostream& output_operator_on_Field(std::ostream& os,const grammar::field<Object,Method>& myField, const Object& myObject){
    io::output_operator_on_element (os,myField.idField);
    os<<io::end_of_line{};
    os<<std::invoke(myField.access_method,myObject)<<io::end_of_line{};
    return os;
}



template<class FieldObject>
auto write_on_Object(std::ostream& os,const FieldObject& myObject)->decltype (myObject.get_constructor_fields(),os)&
{
    auto fields=myObject.get_constructor_fields();
    std::apply([&os, &myObject](const auto&... v){return ((write_on_Field(os,v,myObject)),...,0);},fields);
    os<<io::end_of_Object{};
    return os;
}


template<class FieldObject>
auto output_operator_on_Object(std::ostream& os,const FieldObject& myObject)->decltype (myObject.get_constructor_fields(),os)&
{
    auto fields=myObject.get_constructor_fields();
    std::apply([&os, &myObject](const auto&... v){return ((output_operator_on_Field(os,v,myObject)),...,0);},fields);
    os<<io::end_of_Object{};
    return os;
}


template<class AbstractObject> auto output_operator_on_Abstract_Object(std::ostream& os,const AbstractObject& myObject)->decltype (myObject.get_constructor_fields(),os)&
{
    os<<myObject.myclass()<<io::arg_start();
    output_operator_on_Object(os,myObject);
    os<<io::arg_end();
    return os;
}

template<class AbstractObject>
auto& output_operator_on_Abstract_Object(std::ostream& os,std::enable_if_t<std::is_abstract_v<AbstractObject>,AbstractObject&> myObject,Cs<>)
{

    os<<myObject.myClass()<<" is not completely implemented ";
    return os;
}

template<class AbstractObject, class Derived, class... Deriveds>
auto& output_operator_on_Abstract_Object(std::ostream& os, std::enable_if_t<std::is_abstract_v<AbstractObject>,AbstractObject&> myObject,Cs<Derived,Deriveds...>)
{
    if (myObject.myClass()==Derived::className.str())
    {
        output_operator_on_Object(os,dynamic_cast<Derived&>(myObject));
        return os;

    }
    else return output_operator_on_Abstract_Object(os,myObject,C<Deriveds...>());
}

template<class AbstractObject>
auto& output_operator_on_Abstract_Object(std::ostream& os, std::enable_if_t<std::is_abstract_v<AbstractObject>,AbstractObject*&> myObject)
{
    std::string name;
    os<<myObject;
    os<<io::arg_start();
    output_operator_on_Abstract_Object(os,myObject, Derived_types_t< AbstractObject>());
    return os<<io::arg_end{};

}





template <class Object,class Method>
bool read_on_if_field(std::istream& is,const std::string id,  grammar::field<Object,Method>& myField)
{
    if (id!=myField.idField) return false;
    else
    {
        typename grammar::field<Object,Method>::result_type x;
        if (!read(is,x))
            return false;
        else
        {
            myField.default_value=x;
            return true;
        }
    }

}


template <class Object,class Method>
bool input_operator_on_if_field(std::istream& is,const std::string id,  grammar::field<Object,Method>& myField)
{
    if (id!=myField.idField) return false;
    else
    {
        typename grammar::field<Object,Method>::result_type x;
        if (!(is>>x))
        {
            return false;
        }
        else
        {
            myField.default_value=x;
            return true;
        }
    }

}



template <class Object,class... Method>
std::istream& read_on_this_field(std::istream& is,const std::string id,  std::tuple<grammar::field<Object,Method>...>& myFields){

    return std::apply([&](auto&...x)->auto& {if((read_on_if_field(is,id,x)||...)) return is;else {is.setstate(std::ios::failbit);return is;}},myFields);

}

std::istream& input_operator_on_this_field(std::istream& is,const std::string ,  std::tuple<>& ){
    return is;
}



template <class Object,class... Method>
std::istream& input_operator_on_this_field(std::istream& is,const std::string id,  std::tuple<grammar::field<Object,Method>...>& myFields){
    return std::apply([&](auto&...x)->auto& {if((input_operator_on_if_field(is,id,x)||...)) return is;else {is.setstate(std::ios::failbit);return is;}},myFields);
}



template<class FieldObject>
auto read_on_Object(std::istream& is, FieldObject& myObject)->decltype (myObject.get_constructor_fields(),is)&
{
    auto fields=myObject.get_constructor_fields();
    is>>io::start_of_Container{};
    is>>io::end_of_line{};
    std::string id;
    while (true)
    {
        std::variant<io::end_of_Object, std::string> id;
        static_assert(is_variant<decltype (id)>::value);
        if (!read_variant(is,id))
            return is;
        if (id.index()==1)
        {
            is>>io::end_of_line{};
            if (!read_on_this_field(is,std::get<1>(id),fields))
                return is;
            is>>io::end_of_line{};
        }
        else break;
    }

    if (grammar::has_all(fields))
    {
        myObject=std::apply([](auto& ...x){return FieldObject(x.default_value.value()...);},fields);
        return is;
    }
    else
    {
        is.setstate(std::ios::failbit); return is;
    }
}

template<class FieldObject>
auto input_operator_on_Object(std::istream& is, FieldObject& myObject)->decltype (myObject.get_constructor_fields(),is)&
{
    auto fields=myObject.get_constructor_fields();
    is>>io::start_of_Container{};
    is>>io::end_of_line{};
    std::string id;
    while (true)
    {
        std::variant<io::end_of_Object, std::string> id;
        static_assert(is_variant<decltype (id)>::value);
        if (!input_operator_variant(is,id))
            return is;
        if (id.index()==1)
        {
            is>>io::end_of_line{};
            if (!input_operator_on_this_field(is,std::get<1>(id),fields))
                return is;
            is>>io::end_of_line{};
        }
        else break;
    }

    if (grammar::has_all(fields))
    {
        myObject=std::apply([](auto& ...x){return FieldObject(std::move(x.default_value.value())...);},fields);
        return is;
    }
    else
    {
        is.setstate(std::ios::failbit); return is;
    }
}


template<class AbstractObject>
auto& input_operator_on_Abstract_Object(std::istream& is, AbstractObject*& , std::string ,Cs<>)
{

    return set_fail_bit(is);
}

template<class AbstractObject, class Derived, class... Deriveds>
std::istream& input_operator_on_Abstract_Object(std::istream& is, AbstractObject*& myObject, std::string name,Cs<Derived,Deriveds...>)
{
    if (name==Derived::className.str())
    {
        Derived* d=new Derived;
        input_operator_on_Object(is,*d);
        if (is.good())
        {
            myObject=d;
            return is;
        }
        else return is;

    }
    else return input_operator_on_Abstract_Object(is,myObject,name,Cs<Deriveds...>());
}

template<class AbstractObject>
auto& input_operator_on_Abstract_Object(std::istream& is, AbstractObject*& myObject)
{
    std::string name;
    is>>name;
    is>>io::arg_start();
    input_operator_on_Abstract_Object(is,myObject,name,Derived_types_t<AbstractObject>());
    return is>>io::arg_end{};

}



template<typename T1, typename T2> std::ostream& operator<<(std::ostream& os,const std::pair<T1,T2>& other)
{
    os<<io::start_of_tuple{};
    io::output_operator_on_element(os,other.first);
    os<<io::separator{};
    io::output_operator_on_element(os,other.second);
    os<<io::end_of_tuple{};
    return os;
}

template<typename T1, typename T2> std::istream& operator>>(std::istream& is,const std::pair<T1,T2>& other)
{
    is>>io::start_of_tuple{};
    io::input_operator_on_element(is,other.first);
    io::input_operator_on_element(is,other.second);
    is>>io::end_of_tuple{};
    return is;
}

template<typename T1, typename T2> std::ostream& write_pair(std::ostream& os,const std::pair<T1,T2>& other)
{
    io::write_on_element(os,other.first);
    os<<io::separator{};
    io::write_on_element(os,other.second);
    os<<io::end_of_line{};
    return os;
}

template<typename T1, typename T2> std::istream& read_pair(std::istream& is, std::pair<T1,T2>& other)
{
    io::read_on_element(is,other.first);
    is>>separator{};
    io::read_on_element(is,other.second);
    is>>io::end_of_line{};
    return is;
}



template<typename T1, typename T2>
std::istream& operator>>(std::istream& is,std::pair<T1,T2>& other)
{

    is>>io::start_of_tuple{};
    input_operator_on_element(is,other.first);
    input_operator_on_element(is,other.second);
    is>>io::end_of_tuple{};
    return is;
}



template<class Container>
std::ostream& output_operator_on_container(std::ostream& os,const Container& myContainer){

    os<<io::start_of_Container{};
    os<<io::size_of_container{myContainer.size()};
    for (auto it=myContainer.begin(); it!=myContainer.end();++it)
    {
        os<<io::separator{};
        io::output_operator_on_element(os,*it);
    }
    os<<io::end_of_Container{};
    return os;
}




template<class Container>
std::ostream& write_on_container(std::ostream& os,const Container& myContainer){

    os<<io::size_of_container{myContainer.size()}<<io::end_of_line{};
    for (auto it=myContainer.begin(); it!=myContainer.end();)
    {
        io::write_on_element(os,*it);
        if (++it!=myContainer.end())
            os<<io::separator{};
    }
    os<<io::end_of_line{};
    return os;
}


template<template <class,class> class Container,typename T, class Allocator>
std::istream& input_operator_on_container(std::istream& is, Container<T,Allocator> & myContainer)
{
    io::size_of_container s;
    is>>io::start_of_Container{}>>s;

    myContainer.reserve(s.size);
    for (std::size_t i=0; i<s.size; ++i)
    {
        T e{};
        input_operator_on_element(is,e);
        if (i+1<s.size)
            is>>io::separator{};
        myContainer.push_back(std::move(e));
    }
    is>>io::end_of_Container{};
    return is;
}

template<template <class,class> class Container,typename T, class Allocator>
std::istream& read_on_container(std::istream& is, Container<T,Allocator> & myContainer)
{
    std::optional<io::size_of_container> s;
    read_optional(is,s);
    if (s.has_value())
    {
        myContainer.reserve(s->size);
        for (std::size_t i=0; i<s->size; ++i)
        {
            T e{};
            read_on_element(is,e);
            if (i+1<s->size)
                is>>io::separator{};
            myContainer.push_back(std::move(e));
        }
        is>>io::end_of_line{};
        return is;
    }
    else
    {
        std::variant<io::end_of_line,T> x;
        while(read_on_element(is,x)&&x.index()==1)
        {
            myContainer.push_back(std::move(std::get<1>(x)));

        }
        return is;

    }

}




template<template <class,class, class...> class Map,typename K,typename T, class... Os>
std::istream& input_operator_on_Map(std::istream& is, Map<K,T,Os...> & myMap)
{
    io::size_of_container s;
    is>>io::start_of_Container{}>>s;
    for (std::size_t i=0; i<s.size; ++i)
    {
        std::pair<K,T> p{};
        if (!(is>>p))
            std::cerr<<"input error for std::vector of ";
        if (i+1<s.size)
            is>>io::separator{};
        myMap.insert(std::move(p));
    }

    is>>io::end_of_Container{};
    return is;
}

template<template <class,class, class...> class Map,typename K,typename T, class... Os>
std::istream& read_on_Map(std::istream& is, Map<K,T,Os...> & myMap)
{
    myOptional_t<io::size_of_container> s(false,"");
    read_optional(is,s);
    if (s.has_value())
    {
        for (std::size_t i=0; i<s.value().size; ++i)
        {
            std::pair<K,T> p{};
            read(is,p);
            if (i+1<s.value().size)
                is>>io::separator{};
            myMap.insert(std::move(p));
        }
        is>>io::end_of_line{};
        return is;
    }
    else
    {
        std::variant<io::end_of_line,std::pair<K,T>> x;
        static_assert (is_variant<decltype(x)>::value,"" );
        while(read(is,x)&&x.index()==1)
            myMap.insert(std::move(std::get<1>(x)));
        return is;
    }
}



template<template <class, class...> class Set,typename K,class... Os>
std::istream& input_operator_on_set(std::istream& is, Set<K,Os...> & mySet)
{
    io::size_of_container s;
    is>>io::start_of_Container{}>>s;
    for (std::size_t i=0; i<s.size; ++i)
    {
        K p{};
        if (!(input_operator_on_element(is,p)))
            std::cerr<<"input error for std::vector of ";
        if (i+1<s.size)
            is>>io::separator{};
        mySet.insert(std::move(p));
    }
    is>>io::end_of_Container{};
    return is;
}


template<template <class, class...> class Set,typename K,class... Os>
std::istream& read_on_set(std::istream& is, Set<K,Os...> & mySet)
{
    myOptional_t<io::size_of_container> s(false,"");
    read_optional(is,s);
    if (s.has_value())
    {
        for (std::size_t i=0; i<s.value().size; ++i)
        {
            K p{};
            read(is,p);
            mySet.insert(std::move(p));
            if (i+1<s.value().size)
                is>>io::separator{};
        }
        is>>io::end_of_line{};
        return is;
    }
    else
    {
        std::variant<io::end_of_line,K> x;
        while(read(is,x)&&x.index()==1)
            mySet.insert(std::move(std::get<1>(x)));
        return is;
    }


}





template<typename T> auto operator<<(std::ostream& s, const T& v)->std::enable_if_t<is_std_container<T>::value&&!std::is_same_v<T,std::string >,std::ostream&>
{
    return io::output_operator_on_container(s,v);
}


template<typename T> auto operator>>(std::istream& s, T& v)->std::enable_if_t<is_std_container<T>::value&&!std::is_same_v<T,std::string >, std::istream&>
{
    if constexpr(is_set<T>::value) {
        return io::input_operator_on_set(s,v);
    }
    else if constexpr(is_map<T>::value) {
        return io::input_operator_on_Map(s,v);
    }
    else
    return io::input_operator_on_container(s,v);
}

template<typename T> std::istream& read_on_vector(std::istream& s, std::vector<T>& v)
{
    return io::read_on_container(s,v);
}


template<typename K,typename T> std::ostream& operator<<(std::ostream& s, const std::map<K,T>& v)
{
    return io::output_operator_on_container(s,v);
}


template<typename K,typename T> std::istream& operator>>(std::istream& s, std::map<K,T>& v)
{
    return io::input_operator_on_Map(s,v);
}

template<typename K,typename T> std::ostream& write_map(std::ostream& s, const std::map<K,T>& v)
{
    return io::write_on_container(s,v);
}


template<typename K,typename T> std::istream& read_Map(std::istream& s, std::map<K,T>& v)
{
    return io::read_on_Map(s,v);
}



template<typename K,typename T> std::istream& operator>>(std::istream& s, std::multimap<K,T>& v)
{
    return io::input_operator_on_Map(s,v);
}

template<typename K,typename T> std::ostream& operator<<(std::ostream& s, const std::multimap<K,T>& v)
{
    return io::output_operator_on_container(s,v);
}


template<typename K,typename T> std::istream& read_multimap(std::istream& s, std::multimap<K,T>& v)
{
    return io::read_on_Map(s,v);
}

template<typename K,typename T> std::ostream& write_multimap(std::ostream& s, const std::multimap<K,T>& v)
{
    return io::write_on_container(s,v);
}

template<typename K> std::ostream& operator<<(std::ostream& s, const std::set<K>& v)
{
    return io::output_operator_on_container(s,v);
}

template<typename K> std::istream& operator>>(std::istream& s,  std::set<K>& v)
{
    return io::input_operator_on_set(s,v);
}


template<typename K> std::ostream& write_set(std::ostream& s, const std::set<K>& v)
{
    return io::write_on_container(s,v);
}

template<typename K> std::istream& read_set(std::istream& s,  std::set<K>& v)
{
    return io::read_on_set(s,v);
}






template<typename T> std::ostream& write(std::ostream& s, const T& v)
{
    if constexpr(std::is_arithmetic_v<T>)
            return s<<v;
    else if constexpr (is_std_container<T>::value)
            return io::write_on_container(s,v);
    else if constexpr(is_field_Object<T>::value)
            return io::write_on_Object(s,v);
    else if constexpr(is_write_Object<T>::value)
            return v.write(s);
    else if constexpr(is_tuple<T>::value)
            return write_tuple(s,v);
    else if constexpr(is_pair<T>::value)
            return write_pair(s,v);
    else return s<<v;
    // else statica_assert (false,"not managed" );
}

template<typename T> std::istream& read(std::istream& s, T& v)
{
    if constexpr(std::is_arithmetic_v<T>)
            return    s>>v;
    else if constexpr(is_set<T>::value) {
        return io::read_on_set(s,v);
    }
    else if constexpr(is_map<T>::value) {
        return io::read_on_Map(s,v);
    }
    else if constexpr (is_std_container<T>::value)
            return io::read_on_container(s,v);
    else if constexpr(is_field_Object<T>::value)
            return io::read_on_Object(s,v);
    else if constexpr(is_pair<T>::value)
            return read_pair(s,v);
    else if constexpr(is_tuple<T>::value)
            return read_tuple(s,v);

    else if constexpr(io::is_token<T>::value)
            return s>>v;
    else if constexpr(is_variant<T>::value)
            return read_variant(s,v);
    else return s>>v;
    // else static_assert (false,"not managed" );
}
template<typename T> auto operator<<(std::ostream& s, std::enable_if_t<std::is_abstract_v<T>,T const&> v)
{
    return io::output_operator_on_Abstract_Object(s,v);
}



template<typename T> auto operator<<(std::ostream& s, T const& v)->std::enable_if_t<is_field_Object<T>::value,std::ostream&>
{
    return io::output_operator_on_Object(s,v);
}

template<typename T> auto operator<<(std::ostream& s, T const& v)->std::enable_if_t<is_write_Object<T>::value,std::ostream&>
{
    return v.write(s);
}

template<typename T> auto operator>>(std::istream& s, T & v)->std::enable_if_t<is_read_Object<T>::value,std::istream&>
{
    return v.read(s);
}

template<typename T> auto operator>>(std::istream& s, T*& v)
->std::enable_if_t<std::is_abstract<T>::value,std::istream&>

{
    return io::input_operator_on_Abstract_Object(s,v);
}

template<typename T> auto operator>>(std::istream& s, T& v)->std::enable_if_t<is_field_Object<T>::value,std::istream&>
{
    if constexpr  (std::is_abstract_v<T>)
            return io::input_operator_on_Abstract_Object(s,v);
    else
    return io::input_operator_on_Object(s,v);
}
template<typename ...Args> std::istream& operator>>(std::istream& is,  std::tuple<Args...>& tu)
{
    using io::operator>>;
    is>>io::start_of_Container{};
    std::apply([&is](auto&... v){return ((is>>io::separator{}>>v),...,0);},tu);
    is>>io::end_of_Container{};
    return is;
}


template <typename T> std::string my_to_string(const T& x)
{
    std::ostringstream os;
    os<<x;
    return os.str();
}



} // namespace io



using namespace io;

#endif // MYSERIALIZER_H
