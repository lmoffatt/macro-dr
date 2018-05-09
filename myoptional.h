#ifndef MYOPTIONAL_H
#define MYOPTIONAL_H

#include "mytypetraits.h"
#include <type_traits>
#include <optional>
#include <functional>

template <class> struct is_optional: public std::false_type{};

template <typename T>
struct is_optional<std::optional<T>>: public std::true_type{};
template <typename T>
struct is_optional<std::optional<T>&>: public std::true_type{};

template <typename T>
using is_optional_v=typename is_optional<std::decay_t<T>>::value;

template <class T> struct optional_decay {
    typedef T    type;
};

template <class T> struct optional_decay<std::optional<T>> {
    typedef T    type;
};

template <class T> struct optional_decay<std::optional<T>&> {
    typedef T    type;
};


template< class T> using optional_decay_t=typename optional_decay<T>::type;




template< class F, class... Args>
using invoke_optional_result_t=typename std::invoke_result<optional_decay_t<F>, optional_decay_t<Args>...>::type;


template <class T>
auto optional_value(T&& x)
{
    if constexpr(is_optional<std::decay_t<T>>::value)
       return x.value();
    else
        return std::forward<T>(x);
}


template <class T>
auto optional_has_value(T&& ){return true;}

template <class T>
auto optional_has_value(std::optional<T>&& x){return x.has_value();}


template< class F, class... Args>
auto
invoke_optional(F&& f, Args&&... args)
{
    if ((optional_has_value(std::forward<Args>(args))&&...)&&optional_has_value(f) )
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
auto apply_optional(const F& f, std::tuple<std::optional<Args>...>&& args)
{
    auto res=std::apply([](auto&... x){ return ((x.has_value()&&...));},args);

    if constexpr (contains_constructor<F>::value)
    {

        if (res)
            return std::make_optional(std::apply([](auto&... x){ return typename F::myClass(x.value()...);}, args));
        else return std::optional<typename F::myClass>{};
    }
    else
    {
    if (res)
        return std::make_optional(std::apply([&f](auto&... x){ return f(x.value()...);}, args));
    else
        return std::optional<std::invoke_result_t<F, Args...>>{};

    }
}

#endif // MYOPTIONAL_H