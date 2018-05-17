#ifndef MYSMARTPOINTERSTOOLS_H
#define MYSMARTPOINTERSTOOLS_H

#include <memory>
#include <optional>


template<typename T, typename=void>
struct unique_if_ptr
{
    typedef T type;
};

template<typename T>
struct unique_if_ptr<T, std::enable_if<std::is_pointer_v<T>,int>>
{
    typedef typename std::unique_ptr<T> type;
};

template<typename T>
using unique_if_ptr_t=typename unique_if_ptr<T>::type;

template< typename T>
using optional_unique_t=std::optional<std::unique_ptr<T>>;




template <template<class...>class Vector, class T>
Vector<std::unique_ptr<T>> clone_vector(const Vector<std::unique_ptr<T>>& v){
   Vector<std::unique_ptr<T>> out;
   for (auto& e: v)
       out.emplace_back(e->clone());
   return out;
}

template <template<class...>class Map, class K, class T>
Map<K,std::unique_ptr<T>> clone_map(const Map<K,std::unique_ptr<T>>& v){
   Map<K,std::unique_ptr<T>> out;
   for (auto& e: v)
       out.emplace(e.first,e.second->clone());
   return out;
}

template<class T>
optional_unique_t<T> clone_optional(const std::optional<std::unique_ptr<T>> one)
{
    if (one.has_value())
    {
        return std::unique_ptr<T>(one.value()->clone());
    }
    else return{};
}





#endif // MYSMARTPOINTERSTOOLS_H
