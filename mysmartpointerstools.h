#ifndef MYSMARTPOINTERSTOOLS_H
#define MYSMARTPOINTERSTOOLS_H

#include <memory>

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



#endif // MYSMARTPOINTERSTOOLS_H
