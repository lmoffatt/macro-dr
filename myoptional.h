#ifndef MYOPTIONAL_H
#define MYOPTIONAL_H

#include "mytypetraits.h"
#include <type_traits>
#include <optional>
#include <functional>
#include <fstream>

struct reg_tag{};
struct lref_tag{};
struct const_ref_tag{};
struct base_ptr_tag{};
struct derived_ptr_tag{};
struct void_tag{};


template<typename...> class myOptional;


template<typename T>
struct optional_tag
{
    typedef  std::conditional_t<std::is_same_v<T,void>,void_tag,
             std::conditional_t< std::is_lvalue_reference_v<T>,
             std::conditional_t<std::is_const_v<std::remove_reference_t<T>>,const_ref_tag,lref_tag>,
             std::conditional_t<std::is_pointer_v<T>,
             std::conditional_t<has_base_type_v<std::remove_pointer_t<T> >, derived_ptr_tag, base_ptr_tag>,
              reg_tag>>> type;

};


template<typename T>
using optional_tag_t =typename optional_tag<T>::type;


template<typename T>
class myOptional<T, void_tag>;

template<typename T>
class myOptional<T, lref_tag>;
template<typename T>
class myOptional<T, const_ref_tag>;
template<typename T>
class myOptional<T, derived_ptr_tag>;
template<typename T>
class myOptional<T, base_ptr_tag>;
template<typename T>
class myOptional<T, reg_tag>;


template <typename T>
using myOptional_t=myOptional<T,typename optional_tag<T>::type>;

template<typename T>
class myOptional<T, reg_tag>
{
  T x_;
  bool has_;
  std::string error_;

public :
  //typename T::test test;
  T& value(){return x_;}
  T const & value()const {return x_;}

  bool has_value()const { return has_;}
  std::string error()const {return error_;}

  void setError(std::string&& e){has_=false; error_=std::move(e);}
  //template<typename U=typename T::element_type,std::enable_if_t<std::is_same_v<T, std::unique_ptr<U>>, bool> =true>
  //myOptional(U* x):x_{x},has_{true},error_{}{}


  myOptional(const T& x): x_{x},has_{true},error_{}{}
  myOptional(T&& x): x_{std::move(x)},has_{true},error_{}{}

  myOptional(bool,const std::string& msg):x_{},has_{false},error_{msg}{}
  myOptional(bool,std::string&& msg):x_{},has_{false},error_{std::move(msg)}{}


  operator bool() { return has_;}
  template<class... Args>
  static myOptional* make(Args... a)
  {
      return new myOptional(T(a...));
  }


  void reset(){
      x_={};
      has_=false;
      error_.clear();
  }

  T release(){ has_=false; error_.clear();
               return  std::move(x_);}

};



template<typename T>
class myOptional<T&, lref_tag>
{
  T empty_;
  T& x_;
  bool has_;
  std::string error_;

public :
  //typename T::test test;
  T& value(){return x_;}
  T const & value()const {return x_;}

  bool has_value()const { return has_;}
  std::string error()const {return error_;}

  void setError(std::string&& e){has_=false; error_=std::move(e);}
  //template<typename U=typename T::element_type,std::enable_if_t<std::is_same_v<T, std::unique_ptr<U>>, bool> =true>
  //myOptional(U* x):x_{x},has_{true},error_{}{}


  myOptional(T& x):empty_{}, x_{x},has_{true},error_{}{}

  myOptional(bool,const std::string& msg):empty_{},x_{empty_},has_{false},error_{msg}{}
  myOptional(bool,std::string&& msg):empty_{},x_{empty_},has_{false},error_{std::move(msg)}{}


  operator bool() { return has_;}

};

template<typename T>
class myOptional<const T&, const_ref_tag>
{
  const T& x_;
  bool has_;
  std::string error_;

public :
  //typename T::test test;
  T const & value()const {return x_;}

  bool has_value()const { return has_;}
  std::string error()const {return error_;}

  void setError(std::string&& e){has_=false; error_=std::move(e);}
  //template<typename U=typename T::element_type,std::enable_if_t<std::is_same_v<T, std::unique_ptr<U>>, bool> =true>
  //myOptional(U* x):x_{x},has_{true},error_{}{}


  myOptional(const T& x): x_{x},has_{true},error_{}{}

  myOptional(bool,const std::string& msg):x_{std::move(T{})},has_{false},error_{msg}{}
  myOptional(bool,std::string&& msg):x_{std::move(T{})},has_{false},error_{std::move(msg)}{}


  operator bool() { return has_;}

};



template<typename T>
class myOptional<T*, base_ptr_tag>
{
  std::unique_ptr<T> x_;
  bool has_;
  std::string error_;

public :
  T* value(){return x_.get();}
  const T *   value()const {return x_.get();}
  //typename T::test test;

  bool has_value()const { return has_;}
  std::string error()const {return error_;}

  void setError(std::string&& e){has_=false; error_=std::move(e);}
  //template<typename U=typename T::element_type,std::enable_if_t<std::is_same_v<T, std::unique_ptr<U>>, bool> =true>
  //myOptional(U* x):x_{x},has_{true},error_{}{}


  myOptional(T* x): x_{x},has_{true},error_{}{}

  myOptional(bool,const std::string& msg):x_{nullptr},has_{false},error_{msg}{}
  myOptional(bool,std::string&& msg):x_{nullptr},has_{false},error_{std::move(msg)}{}

  operator bool() { return has_;}

  template<class... Args>
  static myOptional* make(Args... a)
  {
      return new myOptional(new T(a...));
  }


  void reset(){
      x_.reset();
      has_=false;
      error_.clear();
  }

  T* release(){ has_=false; error_.clear(); return  x_.release();}

};

template<typename T>
class myOptional<T*, derived_ptr_tag> :public myOptional_t<typename T::base_type*>
{
  T* x_;
public :
  typedef myOptional_t<typename T::base_type*> base_type;
  typedef typename T::base_type base_element;

  T* value(){return x_;}
  T * /*const*/  value()const {return x_;}
  //typename T::test test;



  myOptional(T* x): base_type(x),x_{x}{}

  myOptional(bool,const std::string& msg)
      :base_type(false,msg),x_{nullptr}{}
  myOptional(bool,std::string&& msg)
      :base_type(false,msg),x_{nullptr}{}

  template<class... Args>
  static myOptional* make(Args... a)
  {
      return new myOptional(new T(a...));
  }


  void reset(){
      x_=nullptr;
      base_type::reset();
  }

  T* release(){ base_type::release(); return  x_;}

};



template<>
class myOptional<void, void_tag>
{
  bool has_;
  std::string error_;

public :

  bool has_value()const { return has_;}
  std::string error()const {return error_;}

  void setError(std::string&& e){has_=false; error_=std::move(e);}
  //template<typename U=typename T::element_type,std::enable_if_t<std::is_same_v<T, std::unique_ptr<U>>, bool> =true>
  //myOptional(U* x):x_{x},has_{true},error_{}{}



  myOptional(bool has,const std::string& msg):has_{has},error_{msg}{}
  myOptional(bool has, std::string&& msg):has_{has},error_{std::move(msg)}{}


  operator bool() { return has_;}

};











template <class T> struct optional_decay {
    typedef T    type;
};

template <class T, class tag> struct optional_decay<myOptional<T,tag>> {
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
auto apply_optional(F&& f, std::tuple<myOptional_t<Args>...>&& args)
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
            return myOptional_t<T>(std::apply([](auto&... x){ return T(std::forward<Args>(x.value())...);}, args));
        else return myOptional_t<T>{};
    }
    else if constexpr (contains_loader<F>::value)
    {
        typedef typename F::myClass T;

        if (res)
        {

            auto f=std::get<myOptional_t<std::string>>(args);
            if (f.has_value())
            {
                std::ifstream fe;
                fe.open(f.value().c_str());
                std::decay_t<T> x;
                x.read(fe);
                if (fe.good()||fe.eof())
                    return myOptional_t<T>(std::move(x));
            }
        }
        return myOptional_t<T>{};
    }
    else if constexpr (contains_valuer<F>::value)
    {
        typedef typename F::myClass T;
        if (res)
        {

            return std::get<myOptional_t<T>>(args);
        }
        else return myOptional_t<typename F::myClass>{};
    }
    else
    {
        typedef std::invoke_result_t<F, Args...> T;
        if (res)
            return myOptional_t<T>(std::apply([&f](auto&... x){ return std::invoke<F,Args...>(std::forward<F>(f),std::forward<Args>(optional_value(x))...);}, args));
        else
            return myOptional_t<T>();

    }
}

template< class F, class... Args>
auto apply_optional(const F& f, myOptional_t<Args>&& ...args)
{
    auto res=(args.has_value()&&...&&true);

    if constexpr (contains_constructor<F>::value)
    {

        typedef typename F::myClass T;
        if (res)
            return myOptional_t<T>(typename F::myClass(args.value()...));
        else
            return myOptional_t<T>{};
    }
    else
    {
        typedef  invoke_optional_result_t<F, Args...> T;
        if (res)
            return myOptional_t<T>(f(args.value()...));
        else
            return myOptional_t<T>{};

    }
}




#endif // MYOPTIONAL_H
