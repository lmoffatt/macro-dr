#ifndef MYOPERATORS_H
#define MYOPERATORS_H
#include<vector>
#include<algorithm>
#include<numeric>
#include<Matrix.h>

M_Matrix<double> operator+(M_Matrix<double> total, const M_Matrix<double>& next)
{
    return Matrix_Binary_Transformations::operator+=(total,next);
}

namespace op {




template<typename T>
T sqr(const T& x){ return x*x;}

template<typename T>
M_Matrix<T> sqr(const M_Matrix<T>& x)
{
    return quadraticForm_XTX(x);
}


template<class T>
T sum(const std::vector<T>& v)
{
     return std::accumulate(v.begin(),v.end(),T());
}

template<class T>
T mean(const std::vector<T>& v)
{
     return sum(v)/v.size();
}





template<class C, class UnOp>
auto sum(const std::vector<C>& v, const UnOp& u)->std::decay_t<std::invoke_result_t<UnOp,C>>
{
    typedef std::decay_t<std::invoke_result_t<UnOp,C>> T;
    return std::accumulate(v.begin(),v.end(),T(),[&u](auto total,const C& next )
    {return total+ std::invoke(u,next);});
}
template<class T>
T sum_sqr(const std::vector<T>& v)
{
    return std::accumulate(v.begin(),v.end(),T(),[](auto total,const T& next ){return total+sqr(next);});
}

template<class C, class UnOp>
auto sum_sqr(const std::vector<C>& v, const UnOp& u)->std::invoke_result_t<UnOp,C>
{
    typedef std::invoke_result_t<UnOp,C> T;
    return std::accumulate(v.begin(),v.end(),T(),[&u](auto total,C next ){return total+sqr(std::invoke(u,next));});
}



template<class C, class UnOp>
auto mean(const std::vector<C>& v, const UnOp& u)->std::decay_t<std::invoke_result_t<UnOp,C>>
{
    return sum(v,u)/v.size();
}


template<class T>
T mean_sqr(const std::vector<T>& v)
{
    return sum_sqr(v)/v.size();
}

template<class C, class UnOp>
auto mean_sqr(const std::vector<C>& v, const UnOp& u)->std::invoke_result_t<UnOp,C>
{
    return sum_sqr(v,u)/v.size();
}

template<class T>
T variance(const std::vector<T>& v)
{
    return mean_sqr(v)-op::sqr(mean(v));
}

template<class C, class UnOp>
auto variance(const std::vector<C>& v, const UnOp& u)->std::decay_t<std::invoke_result_t<UnOp,C>>
{
    return mean_sqr(v,u)-op::sqr(mean(v,u));
}







} // namespace op






#endif // MYOPERATORS_H
