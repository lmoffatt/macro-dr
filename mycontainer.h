#ifndef MYCONTAINER_H
#define MYCONTAINER_H


#include <set>
#include <map>
#include <iostream>

template <typename T>
bool contains(const std::set<T>& big, const std::set<T>& small)
{
    for (auto& e:small)
     if (big.find(e)==big.end()) return false;
    return true;
}

template <typename T, typename K>
bool contains(const std::map<T,K>& big, const std::map<T,K>& small)
{
    for (auto& e:small)
    {
      auto it=big.find(e.first);
      if (it==big.end())
      {
          return false;
      }
      if (!e.second.empty()&& (e.second!=it->second))
      {
          return false;
      }
    }
    return true;
}



#endif // MYCONTAINER_H
