#ifndef MYOUTPUTSERIALIZER_H
#define MYOUTPUTSERIALIZER_H

#include <set>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

#include "myTuples.h"

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os,const std::pair<T1,T2>& other)
{
  os<<other.first<<","<<other.second;
  return os;
}







inline
std::ostream& operator<<(
    std::ostream& s,const std::vector< std::vector< double> >& matrix)
{
    s<<"\n";
    for (std::size_t i=0; i<matrix.size();++i)
    {
        for (std::size_t j=0; j<matrix[i].size();j++)
            s<<matrix[i][j]<<"\t";
        s<<"\n";
    }
    return s;
}

inline
std::ostream& operator<<(
    std::ostream& s,const std::vector<  double> & aVector){

  s<<"[";
  for (std::size_t j=0; j<aVector.size();j++)
        s<<aVector[j]<<"\t";
    s<<"]";
return s;
}



template<typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<T>& v)
{
  s<<"[";
  for (T x:v)
    s<<x<<"\t";
  s<<"]";
  return s;
}

template<typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<std::vector<T>>& m)
{
  s<<"[";
  for (const std::vector<T>& v:m)
    s<<v;
  s<<"]";
  return s;
}

template<typename K,typename T>
std::ostream& operator<<(std::ostream& s, const std::map<K,T>& v)
{
  s<<"{";
  for (auto it=v.begin(); it!=v.end(); ++it)
    s<<*it<<"\t";
  s<<"}";
  return s;
}

template<typename K,typename T>
std::ostream& operator<<(std::ostream& s, const std::multimap<K,T>& v)
{
  for (auto& it=v.begin(); it!=v.end(); ++it)
    s<<*it<<"\t";
  s<<"\n";
  return s;
}

template<typename T>
std::ostream& operator<<(std::ostream& s, const std::set<T>& v)
{
  for (auto it=v.begin(); it!=v.end(); ++it)
    s<<*it<<"\t";
  s<<"\n";
  return s;
}

template<typename T>
std::ostream& operator<<(std::ostream& s, const std::multiset<T>& v)
{
 s<<"{";
  for (auto it=v.begin(); it!=v.end(); ++it)
    s<<*it<<"\t";
  s<<"}";
  return s;
}




template<typename Last>
std::ostream& print_impl(std::ostream& os, const Last& last)
{
   os<<last;
   return os;
}

template<typename First, typename ... Rest>
std::ostream& print_impl(std::ostream& os, const First& first, const Rest&...rest)
{
    print_impl(os,first);
    os<<"\t";
    print_impl(os,rest...);
    return os;
}

template< int ... Indexes, typename ... Args>
std::ostream& print_helper(std::ostream& os, index_tuple<Indexes...>, const std::tuple<Args...>& tup)
{
    print_impl( os, std::get<Indexes>(tup)...);
    return os;
}


template<typename ...Args>
std::ostream& operator<<(std::ostream& os, const std::tuple<Args...> tu)
{
  return print_helper(os,typename make_indexes<Args...>::type(),
                      tu);
}

#endif // MYOUTPUTSERIALIZER_H
