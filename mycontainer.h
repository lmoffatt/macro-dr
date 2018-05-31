#ifndef MYCONTAINER_H
#define MYCONTAINER_H


#include <set>


template <typename T>
bool contains(const std::set<T>& big, const std::set<T>& small)
{
    for (auto& e:small)
     if (big.find(e)==big.end()) return false;
    return true;
}


#endif // MYCONTAINER_H
