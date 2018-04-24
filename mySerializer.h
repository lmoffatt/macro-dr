#ifndef MYOUTPUTSERIALIZER_H
#define MYOUTPUTSERIALIZER_H

#include <set>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

#include "myTuples.h"

#include <cctype>
#include <variant>
template<typename T>
std::ostream& write(std::ostream& os, const T& e);

template<typename T>
std::istream& read(std::istream& os, T& e);


template<typename T>
std::istream& read(std::istream& is, T*& e)
{
    auto o=new T{};
    read(is,*o);
    e=o;
    return is;
}


template<typename T>
std::istream& read(std::istream& is, std::optional<T>& e)
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




template<typename Token,typename T>
std::istream& read(std::istream& is, std::variant<Token,T>& e)
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






inline
std::istream &safeGetline(std::istream &is, std::string &t)
{
    is.clear();
    std::getline(is,t);
    auto it=t.find('\r');
    if (it!=t.npos)
        t.erase(it);
    return is;
}




namespace io
{




template<char... cs> struct token {
    static constexpr char value[]={cs...};
};

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
            is.setstate(std::ios::failbit);
            return is;
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


typedef start_of_Container start_of_tuple;

typedef end_of_Container end_of_tuple;

typedef start_of_Container start_of_pair;

typedef end_of_Container end_of_pair;


typedef token<'\t'> separator;

typedef token<'='> equal;

typedef separator write_separator;


struct size_of_container
{
    std::size_t size;
};

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
    return os<<separator{}<<e;
}

template<typename T>
std::ostream& write_on_element(std::ostream& os, const T& e)
{
    return os<<io::separator{};
    write(os,e);
}

template<typename T>
std::ostream& output_operator_on_element(std::ostream& os, const T * const e)
{
    return os<<separator{}<<*e;
}

template<typename T>
std::ostream& write_on_element(std::ostream& os, const T*& e)
{
    return os<<io::separator{};
    write(os,*e);
}

template<typename T>
std::istream& input_operator_on_element(std::istream& is, T& e)
{
    return is>>separator{}>>e;
}

template<typename T>
std::istream& input_operator_on_element(std::istream& is, T*& e)
{
    e=new T;
    return is>>separator{}>>*e;
}



template<typename T>
std::istream& read_on_element(std::istream& is, T& e)
{
    return is>>io::separator{};
    read(is,e);
}

template<typename T>
std::istream& read_on_element(std::istream& is, T*& e)
{
    return is>>io::separator{};
    e=new T;
    read(is,*e);
}

template<typename Token,typename T>
std::istream& read_on_element(std::istream& is, std::variant<Token,T>& e)
{
    return is>>io::separator{};
    Token t;
    if (is<<t)
    {
        e=t;
        return is;
    }
    else
    {
        is.clear();
        read(is,e);

    }
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
    if ((s=="nan")||(s=="NAN"))
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
    is>>separator{};
    return extract_double(is,e);
}

std::istream& input_operator_on_element(std::istream& is, double*& e)
{
    is>>separator{};
    e=new double;
    return extract_double(is,*e);
}




}


std::istream& read(std::istream& is, double& e)
{
    return io::extract_double(is,e);
}



template<typename ...Args>
std::ostream& operator<<(std::ostream& os, const std::tuple<Args...>& tu)
{
    os<<io::start_of_Container{};
    std::apply([&os](const auto&... v){return ((os<<v<<io::separator{}),...,0);},tu);
    os<<io::end_of_Container{};
    return os;
}


template<typename ...Args>
std::ostream& write(std::ostream& os, const std::tuple<Args...>& tu)
{
    std::apply([&os](const auto&... v){return ((os<<v<<io::separator{}),...,0);},tu);
    os<<io::end_of_line{};
    return os;
}







template<typename ...Args>
std::istream& operator>>(std::istream& is,  std::tuple<Args...>& tu)
{
    is>>io::start_of_Container{};
    std::apply([&is](auto&... v){return ((is>>io::separator{}>>v),...,0);},tu);
    is>>io::end_of_Container{};
    return is;
}

template<typename ...Args>
std::istream& read(std::istream& is,  std::tuple<Args...>& tu)
{
    std::apply([&is](auto&... v){return ((is>>io::separator{}>>v),...,0);},tu);
    is>>io::end_of_line{};
    return is;
}


template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os,const std::pair<T1,T2>& other)
{
    os<<io::start_of_tuple{};
    io::output_operator_on_element(os,other.first);
    io::output_operator_on_element(os,other.second);
    os<<io::end_of_tuple{};
    return os;
}

template<typename T1, typename T2>
std::istream& operator>>(std::istream& is,const std::pair<T1,T2>& other)
{
    is>>io::start_of_tuple{};
    io::input_operator_on_element(is,other.first);
    io::input_operator_on_element(is,other.second);
    is>>io::end_of_tuple{};
    return is;
}


template<typename T1, typename T2>
std::ostream& write(std::ostream& os,const std::pair<T1,T2>& other)
{
    io::write_on_element(os,other.first);
    io::write_on_element(os,other.second);
    os<<io::end_of_line{};
    return os;
}

template<typename T1, typename T2>
std::istream& read(std::istream& is, std::pair<T1,T2>& other)
{
    io::read_on_element(is,other.first);
    io::read_on_element(is,other.second);
    is>>io::end_of_line{};
    return is;
}


namespace io
{
template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os,const std::pair<T1,T2>& other)
{
    os<<io::start_of_tuple{};
    output_operator_on_element(os,other.first);
    output_operator_on_element(os,other.second);
    os<<io::end_of_tuple{};
    return os;
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
    for (auto e:myContainer)
    {
        io::output_operator_on_element(os,e);
    }
    os<<io::end_of_Container{};
    return os;
}




template<class Container>
std::ostream& write_on_container(std::ostream& os,const Container& myContainer){

    os<<io::size_of_container{myContainer.size()}<<io::end_of_line{};
    for (auto e:myContainer)
    {
        io::write_on_element(os,e);
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
        input_operator_on_element(e);
        myContainer.push_back(std::move(e));
    }
    is>>io::end_of_Container{};
    return is;
}

template<template <class,class> class Container,typename T, class Allocator>
std::istream& read_on_container(std::istream& is, Container<T,Allocator> & myContainer)
{
    std::optional<io::size_of_container> s;
    read(is,s);
    if (s.has_value())
    {
        myContainer.reserve(s->size);
        for (std::size_t i=0; i<s->size; ++i)
        {
            T e{};
            read_on_element(is,e);
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
        if (!(is>>io::separator{}>>p))
            std::cerr<<"input error for std::vector of ";
        myMap.insert(std::move(p));
    }
    is>>io::end_of_Container{};
    return is;
}

template<template <class,class, class...> class Map,typename K,typename T, class... Os>
std::istream& read_on_Map(std::istream& is, Map<K,T,Os...> & myMap)
{
    std::optional<io::size_of_container> s;
    read(is,s);
    if (s.has_value())
    {
        for (std::size_t i=0; i<s->size; ++i)
        {
            std::pair<K,T> p{};
            read(is,p);
            myMap.insert(std::move(p));
        }
        is>>io::end_of_line{};
        return is;
    }
    else
    {
        std::variant<io::end_of_line,std::pair<K,T>> x;
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
        mySet.insert(std::move(p));
    }
    is>>io::end_of_Container{};
    return is;
}


template<template <class, class...> class Set,typename K,class... Os>
std::istream& read_on_set(std::istream& is, Set<K,Os...> & mySet)
{
    std::optional<io::size_of_container> s;
    read(is,s);
    if (s.has_value())
    {
        for (std::size_t i=0; i<s->size; ++i)
        {
            K p{};
            read(is,p);
            mySet.insert(std::move(p));
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

};

template<typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<T>& v)
{
    return io::output_operator_on_container(s,v);
}

template<typename T>
std::ostream& write(std::ostream& s, const std::vector<T>& v)
{
    return io::write_on_container(s,v);
}

template<typename T>
std::istream& operator>>(std::istream& s, std::vector<T>& v)
{
    return io::input_operator_on_container(s,v);
}

template<typename T>
std::istream& read(std::istream& s, std::vector<T>& v)
{
    return io::read_on_container(s,v);
}


template<typename K,typename T>
std::ostream& operator<<(std::ostream& s, const std::map<K,T>& v)
{
    return io::output_operator_on_container(s,v);
}


template<typename K,typename T>
std::istream& operator>>(std::istream& s, std::map<K,T>& v)
{
    return io::input_operator_on_Map(s,v);
}

template<typename K,typename T>
std::ostream& write(std::ostream& s, const std::map<K,T>& v)
{
    return io::write_on_container(s,v);
}


template<typename K,typename T>
std::istream& read(std::istream& s, std::map<K,T>& v)
{
    return io::read_on_Map(s,v);
}



template<typename K,typename T>
std::istream& operator>>(std::istream& s, std::multimap<K,T>& v)
{
    return io::input_operator_on_Map(s,v);
}

template<typename K,typename T>
std::ostream& operator<<(std::ostream& s, const std::multimap<K,T>& v)
{
    return io::output_operator_on_container(s,v);
}


template<typename K,typename T>
std::istream& read(std::istream& s, std::multimap<K,T>& v)
{
    return io::read_on_Map(s,v);
}

template<typename K,typename T>
std::ostream& write(std::ostream& s, const std::multimap<K,T>& v)
{
    return io::write_on_container(s,v);
}



template<typename K>
std::ostream& operator<<(std::ostream& s, const std::set<K>& v)
{
    return io::output_operator_on_container(s,v);
}

template<typename K>
std::istream& operator>>(std::istream& s, const std::set<K>& v)
{
    return io::input_operator_on_set(s,v);
}


template<typename K>
std::ostream& write(std::ostream& s, const std::set<K>& v)
{
    return io::write_on_container(s,v);
}

template<typename K>
std::istream& read(std::istream& s, const std::set<K>& v)
{
    return io::read_on_set(s,v);
}







#endif // MYOUTPUTSERIALIZER_H
