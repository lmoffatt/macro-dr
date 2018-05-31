#ifndef MYOPTIONAL_H
#define MYOPTIONAL_H

#include "mytypetraits.h"
#include <type_traits>
#include <optional>
#include <functional>
#include <fstream>



template<typename T>
class myOptional
{
  T x_;
  bool has_;
  std::string error_;

public :
  T& value(){return x_;}
  T const & value()const {return x_;}

  bool has_value()const { return has_;}
  std::string error()const {return error_;}

  //template<typename U=typename T::element_type,std::enable_if_t<std::is_same_v<T, std::unique_ptr<U>>, bool> =true>
  //myOptional(U* x):x_{x},has_{true},error_{}{}


  myOptional(const std::decay_t<T>& x): x_{x},has_{true},error_{}{}
  myOptional(std::decay_t<T>&& x): x_{std::move(x)},has_{true},error_{}{}

  myOptional(bool,const std::string& msg):x_{},has_{false},error_{msg}{}
  myOptional(bool,std::string&& msg):x_{},has_{false},error_{std::move(msg)}{}

  myOptional():x_{},has_{false},error_{}{}

  operator bool() { return has_;}
  template<class... Args>
  T& emplace(Args... a)
  {
      x_=T(a...);
      return x_;
  }

 std::remove_reference_t<T>* operator->() {
     return &x_;}

  T& operator*(){return x_;}

  void reset(){
      has_=false;
      error_.clear();
  }

};



template<typename T>
class myOptional<std::unique_ptr<T>>
{
  std::unique_ptr<T> x_;
  bool has_;
  std::string error_;

public :
  std::unique_ptr<T>& value(){return x_;}
  std::unique_ptr<T> const & value()const {return x_;}

  bool has_value()const { return has_;}
  std::string error()const {return error_;}

  //template<typename U=typename T::element_type,std::enable_if_t<std::is_same_v<T, std::unique_ptr<U>>, bool> =true>
  //myOptional(U* x):x_{x},has_{true},error_{}{}


  template<typename U>
  myOptional(std::unique_ptr<U>&& x):x_{x.release()},has_{true},error_{}{}

  template<typename U>
  myOptional(myOptional<U>&& x):x_{x.has_value()?x.value().release():nullptr},has_{x.has_value()},error_{x.error()}{}


  myOptional(T*&& x): x_{std::move(x)},has_{true},error_{}{}
  myOptional(std::unique_ptr<T>&& x): x_{std::move(x)},has_{true},error_{}{}

  myOptional(bool,const std::string& msg):x_{},has_{false},error_{msg}{}
  myOptional(bool,std::string&& msg):x_{},has_{false},error_{std::move(msg)}{}

  myOptional():x_{},has_{false},error_{}{}

  operator bool() { return has_;}
  template<class... Args>
  std::unique_ptr<T>& emplace(Args... a)
  {
      x_(new T(a...));
      return x_;
  }

 // std::enable_if_t<!std::is_reference_v<T>,T*> operator->() {return &x_;}

  std::unique_ptr<T>& operator*(){return x_;}
};



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
using optional_ref_t_old=typename optional_ref<std::is_lvalue_reference_v<T>,T>::type;

template<typename T>
using optional_ref_t=myOptional<T>;




template <class T> struct optional_decay {
    typedef T    type;
};

template <class T> struct optional_decay<myOptional<T>> {
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
    auto res=std::apply([](auto&... x){
        return ((
                    x.has_value()
                    &&...&&true));
    },
    args);

    if constexpr (contains_constructor<F>::value)
    {

        typedef   typename F::myClass T;
        if (res)
            return optional_ref_t<T>(std::apply([](auto&... x){ return T(std::forward<Args>(x.value())...);}, args));
        else return optional_ref_t<T>{};
    }
    else if constexpr (contains_loader<F>::value)
    {
        typedef typename F::myClass T;

        if (res)
        {

            auto f=std::get<optional_ref_t<std::string>>(args);
            if (f.has_value())
            {
                std::ifstream fe;
                fe.open(f.value().c_str());
                std::decay_t<T> x;
                x.read(fe);
                if (fe.good()||fe.eof())
                    return optional_ref_t<T>(std::move(x));
            }
        }
        return optional_ref_t<T>{};
    }
    else if constexpr (contains_valuer<F>::value)
    {
        typedef typename F::myClass T;
        if (res)
        {

            return std::get<optional_ref_t<T>>(args);
        }
        else return optional_ref_t<typename F::myClass>{};
    }
    else
    {
        typedef std::invoke_result_t<F, Args...> T;
        if (res)
            return optional_ref_t<T>(std::apply([&f](auto&... x){ return std::invoke<F,Args...>(std::forward<F>(f),std::forward<Args>(optional_value(x))...);}, args));
        else
            return optional_ref_t<T>();

    }
}

template< class F, class... Args>
auto apply_optional(const F& f, optional_ref_t<Args>&& ...args)
{
    auto res=(args.has_value()&&...&&true);

    if constexpr (contains_constructor<F>::value)
    {

        typedef typename F::myClass T;
        if (res)
            return optional_ref_t<T>(typename F::myClass(args.value()...));
        else
            return optional_ref_t<T>{};
    }
    else
    {
        typedef  invoke_optional_result_t<F, Args...> T;
        if (res)
            return optional_ref_t<T>(f(args.value()...));
        else
            return optional_ref_t<T>{};

    }
}




#endif // MYOPTIONAL_H
