#ifndef MYDISTRIBUTIONS_H
#define MYDISTRIBUTIONS_H

#include "Matrix.h"
#include <random>



template <typename T, typename S,template <typename> class C>
C<T> emptyCopy(const C<S>& x);

template <typename T, typename S,template <typename> class C>
C<T> emptyCopy(const C<S>& x, const T& e);

template <typename T, typename S>
std::vector<T> emptyCopy(const std::vector<S>& x)
{
    return std::vector<T>(x.size());
}

template <typename T, typename S>
M_Matrix<T> emptyCopy(const M_Matrix<S>& x)
{
    return M_Matrix<T>(x.nrows(),x.ncols(),static_cast<typename M_Matrix<T>::TYPE>(x.type()));
}

template <typename T, typename S>
std::vector<T> emptyCopy(const std::vector<S>& x, const T& e)
{
    return std::vector<T>(x.size(),e);
}

template <typename T, typename S>
M_Matrix<T> emptyCopy(const M_Matrix<S>& x, const T& e)
{
    return M_Matrix<T>(x.nrows(),x.ncols(),e);
}






template<template <typename...> class  My_vec,typename _IntType>
class multinomial_distribution

{
    static_assert(std::is_integral<_IntType>::value,
                  "template argument not an integral type");

public:
    /** The type of the range of the distribution. */
    typedef My_vec<_IntType> result_type;
    /** Parameter type. */

    struct param_type
    {
        typedef multinomial_distribution<My_vec,_IntType> distribution_type;
        friend class multinomial_distribution<My_vec,_IntType>;

        explicit
        param_type(_IntType _N, My_vec<double> P )
            : N_(_N), P_(P)
        {
            _M_initialize();
        }

        _IntType
        N() const
        { return N_; }

        const My_vec<double>&
        P() const
        { return P_; }

        friend bool
        operator==(const param_type& __p1, const param_type& __p2)
        { return __p1.N_ == __p2.N_ && __p1.P_ == __p2.P_; }

    private:
        void
        _M_initialize()
        {
            std::size_t n=P_.size();
            rP_=P_-P_;
            auto s=rP_;
            s[n-1]=P_[n-1];
            for (std::size_t i=1; i<P_.size(); ++i)
            {
                s[n-1-i]=P_[n-1-i]+s[n-i];
            }
            for (std::size_t i=0; i<P_.size(); ++i)
            {
                rP_[i]=P_[i]/s[i];
                P_[i]=P_[i]/s[0];
            }
        }

        _IntType N_;
        My_vec<double> P_;
        My_vec<double> rP_;

    };



    // constructors and member function
    explicit
    multinomial_distribution(_IntType __t,
                             My_vec<double> __p )
        : _M_param(__t, __p)
    { }



    /**
     * @brief Returns the distribution @p t parameter.
     */
    _IntType
    N() const
    { return _M_param.N(); }

    /**
     * @brief Returns the distribution @p p parameter.
     */
    const My_vec<double>&
    P() const
    { return _M_param.P(); }

    /**
     * @brief Returns the parameter set of the distribution.
     */
    param_type
    param() const
    { return _M_param; }

    /**
     * @brief Sets the parameter set of the distribution.
     * @param __param The new parameter set of the distribution.
     */
    void
    param(const param_type& __param)
    { _M_param = __param; }

    /**
     * @brief Returns the greatest lower bound value of the distribution.
     */
    _IntType
    min() const
    { return 0; }

    /**
     * @brief Returns the least upper bound value of the distribution.
     */
    _IntType
    max() const
    { return _M_param.N(); }

    /**
     * @brief Generating functions.
     */
    template<typename _UniformRandomNumberGenerator>
    result_type
    operator()(_UniformRandomNumberGenerator& __urng)
    { return this->operator()(__urng, _M_param); }

    template<typename _UniformRandomNumberGenerator>
    result_type
    operator()(_UniformRandomNumberGenerator& __urng,
               const param_type& __p)
    {
        result_type out(emptyCopy<_IntType>(__p.P_));
        _IntType Nr=__p.N_;
        std::binomial_distribution<_IntType> Bi_;
        typedef typename std::binomial_distribution<_IntType>::param_type biPar;
        for (std::size_t i=0; i< out.size()-1; ++i)
        {
            _IntType ni=Bi_(__urng,biPar(Nr,__p.rP_[i]));
            Nr-=ni;
            out[i]=ni;
        }
        out[out.size()-1]=Nr;
        return out;
    }

    static
    double logP(const result_type& x,const param_type& __p)
    {
        assert(sum(x)==__p.N());
        assert(x.size()==__p.P().size());
        double out=lgamma(__p.N()+1);
        for (std::size_t i=0; i<x.size(); ++i)
            out+= x[i]*log(__p.P()[i])-lgamma(x[i]+1);
        return out;

    }
    double logP(const result_type& x)
    {
        return logP(x,_M_param);
    }
    double P(const result_type& x)
    {
        return std::exp(logP(x));
    }

    static
    double P(const result_type& x,const param_type& __p)
    {
        return std::exp(logP(x,__p));
    }



private:

    param_type _M_param;

};




template<typename T>
class markov_process

{

public:
    /*  * The type of the range of the distribution. */
    typedef M_Matrix<T> result_type;
    /** Parameter type. */

    struct param_type
    {
        typedef markov_process<T> distribution_type;
        friend class markov_process<T>;

        explicit
        param_type(const M_Matrix<T>& _N, const M_Matrix<double>& P )
            : N_(_N), P_(P),rP_(M_Matrix<double>(nrows(P),ncols(P)))

        {
            _M_initialize();
        }

        const M_Matrix<T>&
        N() const
        { return N_; }

        void set_P(M_Matrix<double>&&
                   _P)
        { P_=std::move(_P);
            _M_initialize();
        }

        void set_N(M_Matrix<T>&& _N)
        {  N_=std::move(_N); }

        const M_Matrix<double>&
        P() const
        { return P_; }

        friend bool
        operator==(const param_type& __p1, const param_type& __p2)
        { return __p1.N_ == __p2.N_ && __p1.P_ == __p2.P_; }

    private:
        void
        _M_initialize()
        {
            std::size_t n=ncols(P_);

            auto s=M_Matrix<double>(nrows(P_),ncols(P_));

            for (std::size_t i=0; i<nrows(P_); ++i)
            {
                s(i,n-1)=P_(i,n-1);
                for (std::size_t j=1; j<ncols(P_); ++j)
                {
                    s(i,n-1-j)=P_(i,n-1-j)+s(i,n-j);
                }
                for (std::size_t j=0; j<ncols(P_); ++j)
                {
                    rP_(i,j)=P_(i,j)/s(i,j);
                    P_(i,j)=P_(i,j)/s(i,0);
                }
            }
        }

        M_Matrix<T> N_;
        M_Matrix<double> P_;
        M_Matrix<double> rP_;
    };


    // constructors and member function
    explicit
    markov_process(const M_Matrix<T>& __t,
                   const M_Matrix<double>& __p)
        : _M_param(__t, __p)
    { }



    /**
       * @brief Returns the distribution @p t parameter.
       */
    const M_Matrix<T>&
    N() const
    { return _M_param.N(); }




    /**
       * @brief Returns the distribution @p p parameter.
       */
    const M_Matrix<double>&
    P() const
    { return _M_param.P(); }


    void set_P(M_Matrix<double>&& _P)
    { _M_param.set_P(std::move(_P));
    }

    void set_N(M_Matrix<T>&& _N)
    {
        _M_param.set_N(std::move(_N));
    }


    /**
       * @brief Returns the parameter set of the distribution.
       */
    param_type
    param() const
    { return _M_param; }

    /**
       * @brief Sets the parameter set of the distribution.
       * @param __param The new parameter set of the distribution.
       */
    void
    param(const param_type& __param)
    { _M_param = __param; }

    /**
       * @brief Returns the greatest lower bound value of the distribution.
       */
    T
    min() const
    { return 0; }

    /**
       * @brief Returns the least upper bound value of the distribution.
       */
    T
    max() const
    { return totalsum(_M_param.N()); }

    /**
       * @brief Generating functions.
       */
    template<typename _UniformRandomNumberGenerator>
    result_type
    operator()(_UniformRandomNumberGenerator& __urng)
    { return this->operator()(__urng, _M_param); }

    template<typename _UniformRandomNumberGenerator>
    result_type
    operator()(_UniformRandomNumberGenerator& __urng,
               const param_type& _p, bool sumcols=true)
    {
        if (sumcols)
        {
            std::size_t nc=ncols(_p.P_);
            M_Matrix<T> out(1,nc,0);
            for (std::size_t i=0; i< nrows(_p.P_); ++i)
            {
                T Nr=_p.N_[i];

                for (std::size_t j=1; j< nc-1; ++j)
                {
                    double p=_p.rP_(i,j);
                    auto bi=std::binomial_distribution<T>(Nr,p);
                    T n=bi(__urng);
                    Nr-=n;
                    out(0,j)+=n;
                }
                out(0,nc-1)+=Nr;
            }
            return out;
        }
        else
        {
            std::size_t nc=ncols(_p.P_);
            std::size_t nr=nrows(_p.P_);
            M_Matrix<T> out(nr,nc);
            for (std::size_t i=0; i< nrows(_p.P_); ++i)
            {
                T Nr=_p.N_[i];

                for (std::size_t j=1; j< nc-1; ++j)
                {
                    double p=_p.rP_(i,j);
                    auto bi=std::binomial_distribution<T>(Nr,p);
                    T n=bi(__urng);
                    Nr-=n;
                    out(i,j)=n;
                }
                out(i,nc-1)+=Nr;
            }
            return out;
        }

    }


private:

    param_type _M_param;

};



#endif // MYDISTRIBUTIONS_H
