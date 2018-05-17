#ifndef MYOPTIONAL_H
#define MYOPTIONAL_H

#include "mytypetraits.h"
#include <type_traits>
#include <optional>
#include <functional>



template<bool,typename T> struct optional_ref{};
template<typename T>
struct optional_ref<false,T>
{
   // typedef typename T::garca ga;
    typedef T value_type;
    typedef typename std::optional<T> type;
};

template<typename T>
struct optional_ref<true, T>
{
  //  typedef typename T::chino ga;
    typedef T value_type;
    typedef typename std::optional<std::reference_wrapper<std::remove_reference_t<T>>> type;
};

template<typename T>
using optional_ref_t=typename optional_ref<std::is_lvalue_reference_v<T>,T>::type;





template <class T> struct optional_decay {
    typedef T    type;
};

template <class T> struct optional_decay<std::optional<T>> {
    typedef T    type;
};

template <class T> struct optional_decay<std::optional<T>&> {
    typedef T    type;
};

template <class T> struct optional_decay<std::optional<std::reference_wrapper<T>>> {
    typedef T&    type;
};



template< class T> using optional_decay_t=typename optional_decay<T>::type;




template< class F, class... Args>
using invoke_optional_result_t=typename std::invoke_result<optional_decay_t<F>, optional_decay_t<Args>...>::type;


template <class T>
optional_decay_t<std::decay_t<T>>  optional_value(T&& x)
{
    //typename T::gik k;

    //typename optional_decay_t<std::decay_t<T>>::out gf;

    if constexpr(is_optional<std::decay_t<T>>::value)
    {        if constexpr(is_optional_ref<std::decay_t<T>>::value)
                return x.value().get();
             else
                 return x.value();
    }
    else
    return x;
}



template <class T>
auto optional_has_value(const T& ){return true;}

template <class T>
auto optional_has_value(const std::optional<T>& x){return x.has_value();}


template< class F, class... Args>
auto
invoke_optional(F&& f, Args&&... args)
{

    typename F::function q;
    typename Cs<Args...>::argumentos ar;
    if ((optional_has_value(args)&&...&&true)&&optional_has_value(f) )
        return std::invoke(optional_value(f),optional_value(std::forward<Args>(args))...);
    else
        return invoke_optional_result_t<F, Args...>{};
}

template< class F, class... Args>
auto
invoke_optional_functor(F&& f, Args&&... args)
{
    typedef std::decay_t<std::remove_pointer_t<optional_decay_t<F>>> object_type;

    if ((optional_has_value(std::forward<Args>(args))&&...)&&optional_has_value(f) )
    {

        return std::invoke(&object_type::operator(),optional_value(f),optional_value(std::forward<Args>(args))...);
    }
    else
    {
        return invoke_optional_result_t<decltype(&object_type::operator()),F, Args...>{};
    }
}



template< class F, class... Args>
auto apply_optional(F&& f, std::tuple<optional_ref_t<Args>...>&& args)
{
    auto res=std::apply([](auto&... x){ return ((x.has_value()&&...));},args);

    if constexpr (contains_constructor<F>::value)
    {

        if (res)
            return std::make_optional(std::apply([](auto&... x){ return typename F::myClass(std::forward<Args>(x.value())...);}, args));
        else return std::optional<typename F::myClass>{};
    }
    else
    {
        if (res)
            return std::make_optional(std::apply([&f](auto&... x){ return std::invoke<F,Args...>(std::forward<F>(f),std::forward<Args>(optional_value(x))...);}, args));
        else
            return std::optional<std::invoke_result_t<F, Args...>>();

    }
}

template< class F, class... Args>
auto apply_optional(const F& f, std::optional<Args>&& ...args)
{
    auto res=(args.has_value()&&...&&true);

    if constexpr (contains_constructor<F>::value)
    {

        if (res)
            return std::make_optional(typename F::myClass(args.value()...));
        else
            return std::optional<typename F::myClass>{};
    }
    else
    {
        if (res)
            return std::make_optional(f(args.value()...));
        else
            return std::optional<invoke_optional_result_t<F, Args...>>{};

    }
}




#endif // MYOPTIONAL_H
