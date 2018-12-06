#ifndef MYOPTIONAL_H
#define MYOPTIONAL_H

#include "mytypetraits.h"
#include <type_traits>
#include <optional>
#include <functional>
#include <fstream>
#include <sstream>

struct reg_tag {};
struct lref_tag {};
struct const_ref_tag {};
struct base_ptr_tag {};
struct derived_ptr_tag {};
struct const_base_ptr_tag {};
struct const_derived_ptr_tag {};
struct void_tag {};

template <typename...> class myOptional;

template <typename T> struct optional_tag {
  typedef std::conditional_t<
      std::is_same_v<T, void>, void_tag,
      std::conditional_t<
          std::is_lvalue_reference_v<T>,
          std::conditional_t<std::is_const_v<std::remove_reference_t<T>>,
                             const_ref_tag, lref_tag>,
          std::conditional_t<
              std::is_pointer_v<T>,
              std::conditional_t<
                  std::is_same_v<T, std::add_pointer_t<std::add_const_t<
                                        std::remove_pointer_t<T>>>>,
                  std::conditional_t<has_base_type_v<std::remove_pointer_t<T>>,
                                     const_derived_ptr_tag, const_base_ptr_tag>,
                  std::conditional_t<has_base_type_v<std::remove_pointer_t<T>>,
                                     derived_ptr_tag, base_ptr_tag>>,
              reg_tag>>>
      type;
};

template <typename T> using optional_tag_t = typename optional_tag<T>::type;

template <typename T> class myOptional<T, void_tag>;

template <typename T> class myOptional<T, lref_tag>;
template <typename T> class myOptional<T, const_ref_tag>;
template <typename T> class myOptional<T, derived_ptr_tag>;
template <typename T> class myOptional<T, base_ptr_tag>;
template <typename T> class myOptional<T, const_derived_ptr_tag>;
template <typename T> class myOptional<T, const_base_ptr_tag>;
template <typename T> class myOptional<T, reg_tag>;

template< class T>
struct myOptional_impl
{
  typedef myOptional<T, typename optional_tag<T>::type> type;
};

template< class... Ts>
struct myOptional_impl<myOptional<Ts...>>
{
  typedef myOptional<Ts ...> type;
};



template <typename T>
using myOptional_t = myOptional<T, typename optional_tag<T>::type>;

template <typename T>
using myOptional_op = std::conditional_t<is_optional_v<T>,T,myOptional_t<T>>;



typedef myOptional_t<void> Op_void;


template <typename T> class myOptional<T, reg_tag> {
  T x_;
  bool has_;
  std::string error_;

public:
  typedef T value_type;
  // typename T::test test;
  T &value() & { return x_; }
  T const &value() const & { return x_; }
  T &&value() && { return std::move(x_); }
  //  T const && value()const&& {return std::move(x_);}

  bool has_value() const { return has_; }
  std::string error() const { return error_; }

  void setError(std::string &&e) {
    has_ = false;
    error_ = std::move(e);
  }
  // template<typename U=typename
  // T::element_type,std::enable_if_t<std::is_same_v<T, std::unique_ptr<U>>,
  // bool> =true> myOptional(U* x):x_{x},has_{true},error_{}{}

  myOptional(const T &x) : x_{x}, has_{true}, error_{} {}
  myOptional(T &&x) : x_{std::move(x)}, has_{true}, error_{} {}

  myOptional(bool, const std::string &msg) : x_{}, has_{false}, error_{msg} {}
  myOptional(bool, std::string &&msg)
      : x_{}, has_{false}, error_{std::move(msg)} {}

  operator bool() { return has_; }
  template <class... Args> static myOptional *make(Args... a) {
    return new myOptional(T(a...));
  }

  myOptional(const myOptional &other)
      : x_{other.x_}, has_{other.has_}, error_{other.error_} {}
  myOptional(myOptional &&other)
      : x_{std::move(other.x_)}, has_{std::move(other.has_)},
        error_{std::move(other.error_)} {}
  myOptional &operator=(const myOptional &other) {
    myOptional tmp(other);
    *this = std::move(tmp);
    return *this;
  }
  myOptional &operator=(myOptional &&other) {
    x_ = std::move(other.x_);
    has_ = std::move(other.has_);
    error_ = std::move(other.error_);
    return *this;
  }

  void reset() {
    x_ = {};
    has_ = false;
    error_.clear();
  }

  T release() {
    has_ = false;
    error_.clear();
    return std::move(x_);
  }
};

template <typename T> class myOptional<T &, lref_tag> {
  T empty_;
  T &x_;
  bool has_;
  std::string error_;

public:
  typedef T &value_type;

  // typename T::test test;
  T &value() { return x_; }
  T const &value() const { return x_; }

  bool has_value() const { return has_; }
  std::string error() const { return error_; }

  void setError(std::string &&e) {
    has_ = false;
    error_ = std::move(e);
  }
  // template<typename U=typename
  // T::element_type,std::enable_if_t<std::is_same_v<T, std::unique_ptr<U>>,
  // bool> =true> myOptional(U* x):x_{x},has_{true},error_{}{}

  myOptional(T &x) : empty_{}, x_{x}, has_{true}, error_{} {}

  myOptional(bool, const std::string &msg)
      : empty_{}, x_{empty_}, has_{false}, error_{msg} {}
  myOptional(bool, std::string &&msg)
      : empty_{}, x_{empty_}, has_{false}, error_{std::move(msg)} {}

  operator bool() { return has_; }
  myOptional(const myOptional &other)
      : x_{other.x_}, has_{other.has_}, error_{other.error_} {}
  myOptional(myOptional &&other)
      : x_{other.x_}, has_{std::move(other.has_)}, error_{std::move(
                                                       other.error_)} {}
  myOptional &operator=(const myOptional &other) {
    myOptional tmp(other);
    *this = std::move(tmp);
    return *this;
  }
  myOptional &operator=(myOptional &&other) {
    x_ = std::move(other.x_);
    has_ = std::move(other.has_);
    error_ = std::move(other.error_);
    return *this;
  }
};

template <typename T> class myOptional<const T &, const_ref_tag> {
  const T &x_;
  bool has_;
  std::string error_;

public:
  typedef T const &value_type;

  // typename T::test test;
  T const &value() const { return x_; }

  bool has_value() const { return has_; }
  std::string error() const { return error_; }

  void setError(std::string &&e) {
    has_ = false;
    error_ = std::move(e);
  }
  // template<typename U=typename
  // T::element_type,std::enable_if_t<std::is_same_v<T, std::unique_ptr<U>>,
  // bool> =true> myOptional(U* x):x_{x},has_{true},error_{}{}

  myOptional(const T &x) : x_{x}, has_{true}, error_{} {}

  myOptional(bool, const std::string &msg)
      : x_{std::move(T{})}, has_{false}, error_{msg} {}
  myOptional(bool, std::string &&msg)
      : x_{std::move(T{})}, has_{false}, error_{std::move(msg)} {}
  myOptional(const myOptional &other)
      : x_{other.x_}, has_{other.has_}, error_{other.error_} {}
  myOptional(myOptional &&other)
      : x_{std::move(other.x_)}, has_{std::move(other.has_)},
        error_{std::move(other.error_)} {}
  myOptional &operator=(const myOptional &other) {
    myOptional tmp(other);
    *this = std::move(tmp);
    return *this;
  }
  myOptional &operator=(myOptional &&other) {
    x_ = std::move(other.x_);
    has_ = std::move(other.has_);
    error_ = std::move(other.error_);
    return *this;
  }

  operator bool() { return has_; }
};

template <typename T> class myOptional<T *, base_ptr_tag> {
  std::unique_ptr<T> x_;
  bool has_;
  std::string error_;

public:
  virtual ~myOptional() {}
  typedef T *value_type;

  T *value() & { return x_.get(); }
  const T *value() const & { return x_.get(); }
  T *value() && { return release(); }
  const T *value() const && { return release(); }
  // typename T::test test;

  bool has_value() const { return has_; }
  std::string error() const { return error_; }

  void setError(std::string &&e) {
    has_ = false;
    error_ = std::move(e);
  }
  // template<typename U=typename
  // T::element_type,std::enable_if_t<std::is_same_v<T, std::unique_ptr<U>>,
  // bool> =true> myOptional(U* x):x_{x},has_{true},error_{}{}

  myOptional(T *x) : x_{x}, has_{true}, error_{} {}

  myOptional(bool, const std::string &msg)
      : x_{nullptr}, has_{false}, error_{msg} {}
  myOptional(bool, std::string &&msg)
      : x_{nullptr}, has_{false}, error_{std::move(msg)} {}

  myOptional(const myOptional &other)
      : x_(other.has_value() ? other.value()->clone() : nullptr),
        has_(other.has_value()), error_(other.error()) {}
  myOptional(myOptional &&other)
      : x_(std::move(other.x_)), has_(std::move(other.has_)),
        error_(std::move(other.error_)) {}
  myOptional &operator=(const myOptional &) = default;
  myOptional &operator=(myOptional &&other) {
    x_ = std::move(other.x_);
    has_ = std::move(other.has_);
    error_ = std::move(other.error_);
    return *this;
  };

  operator bool() { return has_; }

  template <class... Args> static myOptional *make(Args... a) {
    return new myOptional(new T(a...));
  }

  void reset() {
    x_.reset();
    has_ = false;
    error_.clear();
  }

  T *release() {
    has_ = false;
    error_.clear();
    return x_.release();
  }
};

template <typename T> class myOptional<T const *, const_base_ptr_tag> {
  T const * x_;
  bool has_;
  std::string error_;

public:
  virtual ~myOptional() {}
  typedef T *value_type;

  T const *value() & { return x_; }
  const T *value() const & { return x_; }
  T const *value() && { return release(); }
  const T *value() const && { return release(); }
  // typename T::test test;

  bool has_value() const { return has_; }
  std::string error() const { return error_; }

  void setError(std::string &&e) {
    has_ = false;
    error_ = std::move(e);
  }
  // template<typename U=typename
  // T::element_type,std::enable_if_t<std::is_same_v<T, std::unique_ptr<U>>,
  // bool> =true> myOptional(U* x):x_{x},has_{true},error_{}{}

  myOptional(T const  *x) : x_{x}, has_{true}, error_{} {}

  myOptional(bool, const std::string &msg)
      : x_{nullptr}, has_{false}, error_{msg} {}
  myOptional(bool, std::string &&msg)
      : x_{nullptr}, has_{false}, error_{std::move(msg)} {}

  myOptional(const myOptional &other)
      : x_(other.x_),
        has_(other.has_value()), error_(other.error()) {}
  myOptional(myOptional &&other)
      : x_(std::move(other.x_)), has_(std::move(other.has_)),
        error_(std::move(other.error_)) {}
  myOptional &operator=(const myOptional &) = default;
  myOptional &operator=(myOptional &&other) {
    x_ = std::move(other.x_);
    has_ = std::move(other.has_);
    error_ = std::move(other.error_);
    return *this;
  };

  operator bool() { return has_; }

  template <class... Args> static myOptional *make(Args... a) {
    return new myOptional(new T(a...));
  }

  void reset() {
    x_=nullptr;
    has_ = false;
    error_.clear();
  }

  T const *release() {
    has_ = false;
    error_.clear();
    auto out=x_;
    x_=nullptr;
    return out;
  }
};

template <typename T>
class myOptional<T *, derived_ptr_tag>
    : public myOptional_t<typename T::base_type *> {
  T *x_;

public:
  virtual ~myOptional() override {}

  typedef myOptional_t<typename T::base_type *> base_type;
  typedef typename T::base_type base_element;
  typedef T *value_type;

  T *value() & { return x_; }
  T * /*const*/ value() const & { return x_; }
  T *value() && { return release(); }
  T * /*const*/ value() const && { return release(); }
  // typename T::test test;

  myOptional(T *x) : base_type(x), x_{x} {}

  myOptional(bool, const std::string &msg)
      : base_type(false, msg), x_{nullptr} {}
  myOptional(bool, std::string &&msg) : base_type(false, msg), x_{nullptr} {}

  template <class... Args> static myOptional *make(Args... a) {
    return new myOptional(new T(a...));
  }

  void reset() {
    x_ = nullptr;
    base_type::reset();
  }

  T *release() {
    base_type::release();
    return x_;
  }
  myOptional(const myOptional &other) : base_type(other), x_{other.x_} {}
  myOptional(myOptional &&other) : base_type(std::move(other)), x_{other.x_} {}
  myOptional &operator=(const myOptional &other) {
    myOptional tmp(other);
    *this = std::move(tmp);
    return *this;
  }
  myOptional &operator=(myOptional &&other) {
    base_type::operator=(std::move(other));
    x_ = std::move(other.x_);
    return *this;
  }
};

template <typename T>
class myOptional<T const *, const_derived_ptr_tag>
    : public myOptional_t<typename T::base_type const *> {
  T const *x_;

public:
  virtual ~myOptional() override {}

  typedef myOptional_t<typename T::base_type *> base_type;
  typedef typename T::base_type base_element;
  typedef T *value_type;

  T const *  value() & { return x_; }
  T const * /*const*/ value() const & { return x_; }
  T const *value() && { return release(); }
  T const * /*const*/ value() const && { return release(); }
  // typename T::test test;

  myOptional(T const *x) : base_type(x), x_{x} {}

  myOptional(bool, const std::string &msg)
      : base_type(false, msg), x_{nullptr} {}
  myOptional(bool, std::string &&msg) : base_type(false, msg), x_{nullptr} {}

  template <class... Args> static myOptional *make(Args... a) {
    return new myOptional(new T(a...));
  }

  void reset() {
    x_ = nullptr;
    base_type::reset();
  }

  T *release() {
    base_type::release();
    return x_;
  }
  myOptional(const myOptional &other) : base_type(other), x_{other.x_} {}
  myOptional(myOptional &&other) : base_type(std::move(other)), x_{other.x_} {}
  myOptional &operator=(const myOptional &other) {
    myOptional tmp(other);
    *this = std::move(tmp);
    return *this;
  }
  myOptional &operator=(myOptional &&other) {
    base_type::operator=(std::move(other));
    x_ = std::move(other.x_);
    return *this;
  }
};

template <> class myOptional<void, void_tag> {
  bool has_;
  std::string error_;

public:
  bool has_value() const { return has_; }
  std::string error() const { return error_; }

  void setError(std::string &&e) {
    has_ = false;
    error_ = std::move(e);
  }
  // template<typename U=typename
  // T::element_type,std::enable_if_t<std::is_same_v<T, std::unique_ptr<U>>,
  // bool> =true> myOptional(U* x):x_{x},has_{true},error_{}{}

  myOptional(bool has, const std::string &msg) : has_{has}, error_{msg} {}
  myOptional(bool has, std::string &&msg) : has_{has}, error_{std::move(msg)} {}

  myOptional() = default;

  operator bool() { return has_; }

  template <class T, class arg> myOptional &operator=(const myOptional<T,arg> &other) {
    has_ = other.has_value();
    error_ = other.error();
    return *this;
  }

  template <class invariant, class... Ts>
  static myOptional test(const invariant &inv, Ts... t) {
    std::stringstream ss;
    if (inv.test(t..., ss))
      return myOptional(true, ss.str());
    else
      return myOptional(false, ss.str());
  }

  myOptional &operator+=(myOptional &&other) {
    has_ = has_ && other.has_value();
    if (!other.error().empty()) {
      if (error().empty())
        error_ = std::move(other.error_);
      else
        error_ = error() + "; " + std::move(other.error_);
    }
    return *this;
  }
};

Op_void operator+(Op_void &&one, Op_void &&other) {
  one += std::move(other);
  return one;
}

template <typename T, typename K>
myOptional_t<T> &operator<<(myOptional_t<T> &o, const K &x) {
  std::stringstream ss;
  ss << x;
  o.setError(o.error() + " " + ss.str());
  return o;
}

template <typename T, typename K>
myOptional_t<T> &operator<<(myOptional_t<T> &o, const myOptional_t<K> &x) {
  o.setError(o.error() + x.error());
  return o;
}

template <template <typename...> class V> Op_void consolidate(V<Op_void> &&v) {
  bool res = true;
  std::string error;
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (!v[i].has_value()) {
      res = false;
      if (error.empty())
        error += std::to_string(i) + ": " + v[i].error();
      else
        error += "; " + std::to_string(i) + ": " + v[i].error();
    }
  }
  return Op_void(res, error);
}

template <template <typename...> class V>
Op_void consolidate(V<std::pair<Op_void, Op_void>> &&v) {
  bool res = true;
  std::string error;
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (!v[i].first.has_value()) {
      res = false;
      if (error.empty())
        error += std::to_string(i) + ": " + v[i].error();
      else
        error += "; " + std::to_string(i) + ": " + v[i].error();
    }
    if (!v[i].second.has_value()) {
      res = false;
      if (error.empty())
        error += std::to_string(i) + ": " + v[i].error();
      else
        error += "; " + std::to_string(i) + ": " + v[i].error();
    }
  }
  return Op_void(res, error);
}



template <template <typename...> class V>
 std::string consolidate(const V<std::string>& errors)
{
  std::string out;
  for (auto& e: errors)
    out+=e;
  return out;
}



template <class T> struct remove_myOptional { typedef T type; };

template <class T, class tag> struct remove_myOptional<myOptional<T, tag>> {
  typedef T type;
};

template <class T>
using remove_myOptional_t = typename remove_myOptional<T>::type;

template <class T> struct optional_decay { typedef T type; };

template <class T, class tag> struct optional_decay<myOptional<T, tag>> {
  typedef T type;
};

template <class T> struct optional_decay<std::optional<T>> { typedef T type; };

template <class T> struct optional_decay<std::optional<T> &> {
  typedef T type;
};

template <class T>
struct optional_decay<std::optional<std::reference_wrapper<T>>> {
  typedef T &type;
};

template <class T> using optional_decay_t = typename optional_decay<T>::type;

template <class F, class... Args>
using invoke_optional_result_t =
    typename std::invoke_result<optional_decay_t<F>,
                                optional_decay_t<Args>...>::type;

template <class T> optional_decay_t<std::decay_t<T>> optional_value(T &&x) {
  // typename T::gik k;

  // typename optional_decay_t<std::decay_t<T>>::out gf;

  if constexpr (is_optional<std::decay_t<T>>::value) {
    if constexpr (is_optional_ref<std::decay_t<T>>::value)
      return x.value().get();
    else
      return x.value();
  } else
    return x;
}

template <class T> auto optional_has_value(const T &) { return true; }

template <class T> auto optional_has_value(const std::optional<T> &x) {
  return x.has_value();
}

template <class F, class... Args> auto invoke_optional(F &&f, Args &&... args) {

  typename F::function q;
  typename Cs<Args...>::argumentos ar;
  if ((optional_has_value(args) && ... && true) && optional_has_value(f))
    return std::invoke(optional_value(f),
                       optional_value(std::forward<Args>(args))...);
  else
    return invoke_optional_result_t<F, Args...>{};
}

template <class F, class... Args>
auto invoke_optional_functor(F &&f, Args &&... args) {
  typedef std::decay_t<std::remove_pointer_t<optional_decay_t<F>>> object_type;

  if ((optional_has_value(std::forward<Args>(args)) && ...) &&
      optional_has_value(f)) {

    return std::invoke(&object_type::operator(), optional_value(f),
                       optional_value(std::forward<Args>(args))...);
  } else {
    return invoke_optional_result_t<decltype(&object_type::operator()), F,
                                    Args...>{};
  }
}

#endif // MYOPTIONAL_H
