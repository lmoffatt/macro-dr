#ifndef MYINPUTSERIALIZER_H
#define MYINPUTSERIALIZER_H

#include <set>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <limits>
#include "myTuples.h"

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


template<typename T>
struct value_wrapper
{
  value_wrapper(T& x):value(x){}


  T& value;
};

inline
std::istream& extract_infinite(std::istream& ss, value_wrapper<double>& val, bool is_negative)
{
  std::string s;
  ss>>s;
  if ((s=="inf")||(s=="infinity")||(s=="INF")||(s=="INFINITY"))
    {
      if (is_negative)
        val.value=-std::numeric_limits<double>::infinity();
      else
        val.value=std::numeric_limits<double>::infinity();
      return ss;
    }
  else
    {
      ss.setstate(std::ios::failbit);
      return ss;
    }
}

inline
std::istream& extract_nan(std::istream& ss, value_wrapper<double>& val, bool is_negative)
{
  std::string s;
  ss>>s;
  if ((s=="nan")||(s=="NAN"))
    {
      if (is_negative)
        val.value=-std::numeric_limits<double>::quiet_NaN();
      else
        val.value=std::numeric_limits<double>::quiet_NaN();
      return ss;
    }
  else
    {
      ss.setstate(std::ios::failbit);
      return ss;
    }
}

inline
std::istream& extract_finite(std::istream& ss, value_wrapper<double>& val, bool is_negative)
{
  ss>>val.value;
  if (is_negative)
    val.value=-val.value;
  return ss;

}


inline
std::istream& operator>>(std::istream& ss, value_wrapper<double>& val)
{
  char c;
  bool is_negative=false;
  ss>>c;
  while (std::isspace(c))
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

template <typename T>
std::istream& operator>>(std::istream& ss, value_wrapper<T>& val)
{
  ss>>val.value;
  return ss;
}

template<typename T>
auto operator>>(std::istream& is, T& v)
->decltype(v.read(std::declval<std::string&>(),is,std::declval<std::ostream&>()),is)
{
  std::string s;
  if (!v.read(s,is, std::cerr))
    {
      is.setstate(std::ios::failbit);

    }
  return is;

}




template<typename T>
auto operator>>(std::istream& is, T*& v)
->decltype(v->read(std::declval<std::string&>(),is,std::declval<std::ostream&>() ),is)
{
  std::string s;
  if (v!=nullptr)
    {
      v->read(s,is,std::cerr);
    }
  else
    {
      is.setstate(std::ios::failbit);
    }
  return is;

}




template<typename T
         ,typename std::enable_if<!std::is_pointer<T>::value,int>::type = 0>
std::istream& operator>>(std::istream& is, std::vector<T>& v)
{
  std::string line;
  char c;
  is>>c;
  if (c=='[')
    {
      v.clear();
      while (c!=']')
        {
          is>>c;
          if (c!=']')
            {
              is.putback(c);
              T e;
              value_wrapper<T> w(e);
              if(is>>w)
                v.push_back(e);
              else
                break;
            }
        }
      return is;
    }
  else
    {
      safeGetline(is,line);
      while (v.empty()  &&line.empty()&& is.good())
        safeGetline(is,line);

      T x;
      value_wrapper<T> w(x);
      std::stringstream ss(line);
      while (is>>w)
        v.push_back(x);
      return is;
    }
}
template<typename T>

std::istream& operator>>(std::istream& is, std::vector<T*>& v)
{
  std::string line;
  char c;
  is>>c;
  if (c=='[')
    {
      while (c!=']')
        {
          is>>c;
          if (c!=']')
            {
              T* e=new T{};
              value_wrapper<T> w(*e);
              is.putback(c);
              if(is>>w)
                v.push_back(e);
              else
                break;
            }
        }
      return is;
    }
  else
    {
      safeGetline(is,line);
      while (v.empty()  &&line.empty()&& is.good())
        safeGetline(is,line);

      while (true)
        {
          T* x=new T;
          value_wrapper<T> w(*x);
          std::stringstream ss(line);
          if (is>>w)
            {
              v.push_back(x);
            }
          else
            {
              delete x;
              break;
            }
        }
      return is;
    }
}


template<typename T>
std::istream& operator>>(std::istream& is, std::vector<std::vector<T>>& m)
{
  char c;
  is>>c;
  if (c=='[')
    {
      while (c!=']')
        {
          is>>c;
          if (c!=']')
            {
              is.putback(c);
              std::vector<T> e;
              if(is>>e)
                m.push_back(e);
              else
                break;
            }
        }
      return is;
    }
  else
    {
      is.putback(c);
      std::vector<T> v;
      while((is>>v)&& !v.empty())
        {
          m.push_back(v);
          v.clear();
        }
      return is;
    }
}




template<typename K,typename T>
std::istream& operator>>(std::istream& is, std::map<K,T>& v)
{ char c;
  is>>c;
  if (c=='{')
    {
      while (c!='}')
        {
          is>>c;
          while (std::isspace(c))
            is>>c;
          if (c!='}')
            {
              is.putback(c);
              std::pair<K,T> e;
              if(is>>e)
                v.insert(e);
              else
                break;
            }
        }
      return is;
    }
  else
    {
      is.setstate(std::ios::failbit);
      return is;
    }
}

template<typename T1, typename T2>
std::istream& operator>>(std::istream& os,std::pair<T1,T2>& other)
{
  char ch;
  value_wrapper<T1> w(other.first);
  value_wrapper<T2> w2(other.second);

  os>>w>>ch>>w2;
  return os;
}

template<typename K,typename T>
std::istream& operator>>(std::istream& is,  std::multimap<K,T>& v)
{
  char c;
  is>>c;
  if (c=='{')
    {
      while (c!='}')
        {
          is>>c;
          if (c!='}')
            {
              std::pair<K,T> e;
              is.putback(c);
              if(is>>e)
                v.insert(e);
              else
                break;
            }
        }
      return is;
    }
  else
    {
      is.setstate(std::ios::failbit);
    }
}
template<typename T>
std::istream& operator>>(std::istream& is, std::multiset<T>& v)
{
  char c;
  is>>c;
  if (c=='{')
    {
      while (c!='}')
        {
          is>>c;
          if (c!='}')
            {
              T e;
              is.putback(c);
              if(is>>value_wrapper<T>(e))
                v.insert(e);
              else
                break;
            }
        }
      return is;
    }
  else
    {
      is.setstate(std::ios::failbit);
      return is;
    }
}


template<typename Last>
std::istream& get_impl(std::istream& is, Last& last)
{
  is>>last;
  return is;
}

template<typename First, typename ... Rest>
std::istream& get_impl(std::istream& is,  First& first, Rest&...rest)
{
  get_impl(is,first);
  get_impl(is,rest...);
  return is;
}

template< int ... Indexes, typename ... Args>
std::istream& get_helper(std::istream& is, index_tuple<Indexes...>, std::tuple<Args...>& tup)
{
  get_impl( is, std::get<Indexes>(tup)...);
  return is;
}


template<typename ...Args>
std::istream& operator>>(std::istream& is,  std::tuple<Args...>& tu)
{
  return get_helper(is,typename make_indexes<Args...>::type(),
                    tu);
}



#endif // MYINPUTSERIALIZER_H
