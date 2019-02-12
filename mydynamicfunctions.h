#ifndef MYDYNAMICFUNCTIONS_HPP
#define MYDYNAMICFUNCTIONS_HPP

#include "myparameters.h"
#include "myparameters_derivative.h"

template <typename T, class Model> class Base_Function {
public:
  constexpr static auto const className =
      my_static_string("Base_Function") + my_trait<T>::className;
  virtual std::string myClass() const = 0;

  virtual Base_Function *clone() const = 0;
  virtual T operator()(const Parameters_values<Model> &x) const = 0;
  virtual Derivative<T> operator()(const Derivative<Parameters_values<Model>> &x) const = 0;

  virtual void putme(std::ostream &os) const = 0;

  virtual ~Base_Function() {}
};

template <typename T, class Model>
std::ostream &operator<<(std::ostream &os, const Base_Function<T, Model> &f) {
  f.putme(os);
  return os;
}







struct compiler_grammar {
  inline static const char group_start = '(';
  inline static const char group_end = ')';
  inline static const char argument_separator = ',';

  inline static const char sum_operator = '+';
  inline static const char substraction_operator = '-';
  inline static const char prod_operator = '*';
  inline static const char power_operator = '^';
  inline static const char div_operator = '/';

  inline static const char plus_sign = '+';
  inline static const char minus_sign = '-';

  inline static const std::string space = " \t";
  inline static const std::string number = "0123456789";
  inline static const std::string realnumber = "0123456789.eE+-";

  inline static const std::string alfa =
      "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
  inline static const std::string alfanum = alfa + number + "_";
  inline static const std::string sign = "+-";

  inline static const std::string binary_operators = {
      sum_operator, substraction_operator, prod_operator, power_operator,
      div_operator};
  inline static const std::string unary_operators = {plus_sign, minus_sign};
};




template <typename T, class Model>
class ParValue : public Base_Function<T, Model> {
public:
  inline static const std::string legal_chars =
      Parameters_values<Model>::legal_chars;
  constexpr static auto const className =
      my_static_string("ParValue") + my_trait<double>::className;
  typedef typename Parameters_values<Model>::par_label par_label;
  virtual std::string myClass() const override { return className.str(); }

  virtual ParValue *clone() const override{ return new ParValue(*this); };
  virtual T operator()(const Parameters_values<Model> &x) const override {
    return x.at(p_);
  }
  virtual Derivative<T> operator ()(const Derivative<Parameters_values<Model> > &x) const override
  {
    return x.at(p_);

  }
  par_label getParameter() const { return p_; }

  virtual void putme(std::ostream &os) const override { os << getParameter(); }

  ParValue(const std::string &p) : p_{p} {}

   ~ParValue() override {}

private:
  par_label p_;


};

template <class Model> class ConstValue : public Base_Function<double, Model> {
public:
  constexpr static auto const className =
      my_static_string("ConstValue") + my_trait<double>::className;
  typedef typename Parameters_values<Model>::par_label par_label;
  virtual std::string myClass() const override { return className.str(); }

  virtual ConstValue *clone() const  override { return new ConstValue(*this); };
  virtual double operator()(const Parameters_values<Model> &) const override {
    return x_;
  }

  virtual Derivative<double> operator()(const Derivative<Parameters_values<Model>> &p) const override
      {
    return Derivative<double>(x_,p.x());
      };


  double getValue() const { return x_; }
  ConstValue(double x) : x_{x} {}
   ~ConstValue() override {}

  virtual void putme(std::ostream &os) const { os << getValue(); }

private:
  double x_;
};

template <class BiOp, typename T, class Model>
class BinaryOperator : public Base_Function<T, Model> {
public:
  constexpr static auto const className =
      my_trait<BiOp>::className + my_static_string("_BinaryOp");
  typedef typename Parameters_values<Model>::par_label par_label;
  virtual std::string myClass() const override { return className.str(); }
  virtual BinaryOperator *clone() const override{ return new BinaryOperator(*this); };
  virtual T operator()(const Parameters_values<Model> &x) const override {
    return BiOp()(first()(x), second()(x));
  }
  virtual Derivative<T> operator()(const Derivative<Parameters_values<Model>> &x) const override{
    return Derivative_t<BiOp>()(first()(x), second()(x));
  }

  BinaryOperator(Base_Function<T, Model> *one, Base_Function<T, Model> *two)
      : one_{std::move(one)}, two_{std::move(two)} {}
   ~BinaryOperator() override {}
  auto &first() const { return *one_; }
  auto &second() const { return *two_; }

  virtual void putme(std::ostream &os) const override{
    os << first() << BiOp::symbol << second();
  }

  BinaryOperator(const BinaryOperator& other):one_(other.one_->clone()),two_{other.two_->clone()}{}
  BinaryOperator(BinaryOperator&& other):one_(std::move(other.one_)),two_{std::move(other.two_)}{}
  BinaryOperator& operator=(const BinaryOperator& other)
  {
    BinaryOperator tmp(other);
    *this=std::move(tmp);
    return *this;
  }
  BinaryOperator &operator=(BinaryOperator &&other) {
    one_=std::move(other.one_);
    two_=std::move(other.two_);
    return *this;
  }

  private:
  std::unique_ptr<Base_Function<T, Model>> one_;
  std::unique_ptr<Base_Function<T, Model>> two_;
};

template <class UnOp, typename T, class Model>
class UnaryOperator : public Base_Function<T, Model> {
public:
  constexpr static auto const className =
      my_trait<UnOp>::className + my_static_string("_BinaryOp");
  typedef typename Parameters_values<Model>::par_label par_label;
  virtual std::string myClass() const { return className.str(); }
  virtual UnaryOperator *clone() const { return new UnaryOperator(*this); };
  virtual T operator()(const Parameters_values<Model> &x) const {
    return UnOp()(first()(x));
  }
  virtual Derivative<T> operator ()(const Derivative<Parameters_values<Model> > &x) const override
  {
    return Derivative_t<UnOp>()(first()(x));
  }


  UnaryOperator(Base_Function<T, Model> *one)
      : one_{std::move(one)} {}
   ~UnaryOperator() override {}
  auto &first() const { return *one_; }
  virtual void putme(std::ostream &os) const { os << UnOp::symbol << first(); }

  UnaryOperator(const UnaryOperator& other):one_(other.one_->clone()){}
  UnaryOperator(const UnaryOperator&& other):one_(std::move(other.one_)){}
  UnaryOperator& operator=(const UnaryOperator& other){UnaryOperator tmp(other); *this=std::move(tmp);}
  UnaryOperator& operator=(UnaryOperator&& other){one_=std::move(other.one_);}
private:
  std::unique_ptr<Base_Function<T, Model>> one_;

};

template <typename T, class Model>
class GroupOperator : public Base_Function<T, Model> {
public:
  inline constexpr static auto const start_symbol =
      compiler_grammar::group_start;
  inline constexpr static auto const end_symbol = compiler_grammar::group_start;

  constexpr static auto const className = my_static_string("Group_Op");
  typedef typename Parameters_values<Model>::par_label par_label;
  virtual std::string myClass() const { return className.str(); }
  virtual GroupOperator *clone() const { return new BinaryOperator(*this); }
  virtual T operator()(const Parameters_values<Model> &x) const {
    if (negative_)
      return -me_(x);
    else
      return me_(x);
  }
  GroupOperator(bool is_negative, Base_Function<T, Model> *one)
      : me_{std::move(one)}, negative_{is_negative} {}
   ~GroupOperator() override{}
  auto const &expression() const { return *me_.value(); }
  virtual void putme(std::ostream &os) const {
    os << start_symbol << expression() << end_symbol;
  }

private:
  std::unique_ptr<Base_Function<T, Model>> me_;
  bool negative_;
};

template <class F, typename Args, typename T, class Model> class Function;

template <typename T, class Model, class F, typename... Args>
class Function<F, Cs<Args...>, T, Model> : public Base_Function<T, Model> {
public:
  constexpr static auto const className =
      my_trait<F>::className + my_static_string("_f");
  typedef typename Parameters_values<Model>::par_label par_label;
  virtual std::string myClass() const override{ return className.str(); }
  virtual Function *clone() const override{ return new Function(*this); }
  virtual T operator()(const Parameters_values<Model> &x) const override {
    return std::apply([&x](auto &... t) { return F()((*t)(x)...); }, tu_);
  }

  virtual Derivative<T> operator()(const Derivative<Parameters_values<Model>> &x) const override {
    return std::apply([&x](auto &... t) { return Derivative_t<F>()((*t)(x)...); }, tu_);
  }

  Function(Base_Function<Args, Model> *... args) : tu_{args...} {}
   ~Function() override {}
  template <std::size_t I> auto &get_Arg() const { return *std::get<I>(tu_); }
  virtual void putme(std::ostream &os) const override {
    os << F::className.str() << compiler_grammar::group_start;
    std::apply([&os](auto &... a) { print_arg(os, a...); }, tu_);
    os << compiler_grammar::group_end;
  }

  Function(const Function& other): tu_{clone_tuple(other.tu_)}{}
  Function(Function&& other): tu_{std::move(other.tu_)}{}
  Function& operator=(const Function& other){Function tmp(other); *this=std::move(tmp); return *this;}
  Function& operator=(Function&& other){tu_=std::move(other.tu_); return *this;}


private:
  std::tuple<std::unique_ptr<Base_Function<Args, Model>>...> tu_;

  template <typename A, typename... As>
  static std::ostream &
  print_arg(std::ostream &os, const std::unique_ptr<Base_Function<A, Model>>& a,
            const std::unique_ptr<Base_Function<As, Model>>&... as) {
    if constexpr (sizeof...(As) > 0) {
      os << *a << compiler_grammar::argument_separator << " ";
      return print_arg(os, as...);
    } else {
      return os << *a;
    }
  }
};

template <typename T, class Model> struct Base_Compiler {
public:
  typedef myOptional_t<std::pair<Base_Function<T, Model> *, std::size_t>> Op;

  virtual Op
  compile(const std::string& s, std::size_t pos) const = 0;
  virtual ~Base_Compiler() {}
};

template <typename T, class Model> struct Base_Binary_Compiler {
public:
  typedef myOptional_t<std::pair<Base_Function<T, Model> *, std::size_t>> Op;

  virtual myOptional_t<std::pair<Base_Function<T, Model> *, std::size_t>>
  compile(Base_Function<T, Model> *first_term, const std::string s,
          std::size_t pos) const = 0;

  virtual std::size_t operator_precedence() const = 0;
  virtual ~Base_Binary_Compiler() {}
};

template <typename T, class Model> struct compile_expression;

template <typename T, class Model> struct compile_term;

template <typename T, class Model> struct compile_identifier;

template <typename T, class Model> struct compile_number;

template <typename T, class Model> struct compile_unary;

template <typename T, class Model> struct compile_binary;

template <class F, typename Args, typename T, class Model>
struct compile_function_typed;
template <class BinOp, typename T, class Model>
struct compile_binary_operator_typed;

template <class Model, typename... Ts, typename R, typename... Rs>
myOptional_t<std::pair<
    std::tuple<Base_Function<Ts, Model> *..., Base_Function<R, Model> *,
               Base_Function<Rs, Model> *...>,
    std::size_t>>
extract_argument_list(const std::string &s, std::size_t pos, Cs<R, Rs...>,
                      Base_Function<Ts, Model> *... tu) {
  typedef myOptional_t<std::pair<
      std::tuple<Base_Function<Ts, Model> *..., Base_Function<R, Model> *,
                 Base_Function<Rs, Model> *...>,
      std::size_t>>
      Op;
  auto res = compile_expression<R, Model>().compile(s, pos);
  if (!res)
    return Op(false, " Argument " + std::to_string(sizeof...(Ts) + 1) +
                         "th error: " + res.error());
  else {

    auto [arg,end] = std::move(res).value();
    if constexpr (sizeof...(Rs) == 0) {
      if (s[end] == compiler_grammar::group_end)
        return Op({std::tuple(tu...,arg), pos});
      else {
        return Op(false, std::string("missing final "
                                     "") +
                             compiler_grammar::group_end +
                             ""
                             " for function");
      }
    } else if (s[end] == compiler_grammar::argument_separator) {
      return extract_argument_list<Model>(s, res.value().second + 1,
                                          Cs<Rs...>(), tu...,arg);
    } else {
      return Op(false, std::string("missing "
                                   "") +
                           compiler_grammar::argument_separator +
                           ""
                           " after the " +
                           std::to_string(sizeof...(Ts) + 1) + " argument");
    }
  }
}

template <typename T, class Model, class F, typename... Args>
struct compile_function_typed<F, Cs<Args...>, T, Model>
    : public Base_Compiler<T, Model> {
  typedef Base_Compiler<T, Model> base_type;
  using typename base_type::Op;
  virtual Op compile(const std::string &s, std::size_t pos) const override {
    auto res = extract_argument_list<Model>(s, pos, Cs<Args...>());
    if (!res)
      return Op(false,
                "Error in position" + std::to_string(pos) + "of function " +
                        my_trait<Cs<F, Args...>>::className.str() + res.error());
    else {
      auto [tu,next] = std::move(res).value();
      auto f=std::apply([](auto...p){ return new Function<F, Cs<Args...>, T, Model>(p...);},tu);
      return Op(
          {f,next});
    }
  }
   ~compile_function_typed() override{}
};

template <typename T, class Model> struct get_compile_function {
  struct exp_derivative {
    constexpr static auto const className = my_static_string("exp_derivative");
    Derivative<double> operator()(const Derivative<double>& x) const {
      auto f=std::exp(x.f());
      auto dfdx=x.dfdx().apply([&f](auto const & x){return f*x;});
      return Derivative<double>(f,x.x(),std::move(dfdx));
    }
  };

  struct exp {
    constexpr static auto const className = my_static_string("exp");
    double operator()(double x) const { return std::exp(x); }
    typedef exp_derivative Derivative;
  };


  struct log_derivative {
    constexpr static auto const className = my_static_string("log_derivative");
    Derivative<double> operator()(const Derivative<double> &x) const {
      auto f = std::log(x.f());
      auto dfdx = x.dfdx().apply([&x](auto const &dx) { return (1.0/x.f()) * dx; });
      return Derivative<double>(f, x.x(), std::move(dfdx));
    }
  };

  struct log {
    constexpr static auto const className = my_static_string("log");
    double operator()(double x) const { return std::log(x); }
    typedef log_derivative Derivative;
  };

  struct log10_derivative {
    constexpr static auto const className = my_static_string("log10_derivative");
    Derivative<double> operator()(const Derivative<double> &x) const {
      auto f = std::log10(x.f());
      auto dfdx =
          x.dfdx().apply([&x](auto const &dx) { return (std::log(10.0) / x.f()) * dx; });
      return Derivative<double>(f, x.x(), std::move(dfdx));
    }
  };

  struct log10 {
    constexpr static auto const className = my_static_string("log10");
    double operator()(double x) const { return std::log(x); }
    typedef log10_derivative Derivative;
  };

inline static const std::map<std::string,
                               std::unique_ptr<Base_Compiler<T, Model>>>
      function_map = make_unique_map(std::map<std::string,Base_Compiler<T,Model>*>({
    {exp::className.str(),new compile_function_typed<exp, Cs<double>, T, Model>},
          {log::className.str(),new compile_function_typed<log, Cs<double>, T, Model>},
          {log10::className.str(),new compile_function_typed<log10, Cs<double>, T, Model>}

      }));

  myOptional_t<Base_Compiler<T, Model> const *>
  get_function(const std::string &id) const {
    typedef myOptional_t<Base_Compiler<T, Model> const *> Op;

    auto it = function_map.find(id);
    if (it != function_map.end()) {
      return Op(it->second.get());
    } else
      return Op(false, "function " + id + " not found");
  }
};

template <typename T, class Model>
struct compile_identifier : public Base_Compiler<T, Model> {
  typedef Base_Compiler<T, Model> base_type;
  using typename base_type::Op;
  virtual Op compile(const std::string &s, std::size_t pos) const override {
    assert((compiler_grammar::alfa.find(s[pos]) != std::string::npos));
    auto end = s.find_first_not_of(compiler_grammar::alfanum, pos);
    std::string id = s.substr(pos, end - pos);
    end = s.find_first_not_of(compiler_grammar::space, end);
    if ((end < s.size()) && (s[end] == compiler_grammar::group_start)) {
      auto f = get_compile_function<T, Model>().get_function(id);
      if (!f)
        return Op(false, f.error());
      else {
        return f.value()->compile(s, end);
      }
    } else {
      return Op({new ParValue<T, Model>(id), end});
    }
  }
   ~compile_identifier()override {}
};

template <typename T, class Model>
struct compile_number : public Base_Compiler<T, Model> {
  typedef Base_Compiler<T, Model> base_type;
  ~compile_number()override{}
  using typename base_type::Op;
  virtual Op compile(const std::string &s, std::size_t pos) const override {

    assert((compiler_grammar::number.find(s[pos]) != std::string::npos));
    auto end = s.find_first_not_of(compiler_grammar::realnumber, pos);
    std::stringstream ss(s.substr(pos, end - pos));
    double d;
    if (!(ss >> d))
      return Op(false, " invalid number :" + ss.str());
    else {
      return Op({new ConstValue<Model>(d), end});
    }
  }
};

template <typename T, class Model>
struct compile_group : public Base_Compiler<T, Model> {
  typedef Base_Compiler<T, Model> base_type;
  using typename base_type::Op;
  virtual myOptional_t<std::pair<Base_Function<T, Model> *, std::size_t>>
  compile(const std::string &s, std::size_t pos) const override {
    assert(s[pos] == compiler_grammar::group_start);
    auto res = compile_expression<T, Model>().compile(s, pos + 1);
    if (!res)
      return Op(false, res.error());
    else {
      auto [e, pos] = std::move(res).value();
      if (s[pos] != compiler_grammar::group_end)
        return Op(false, std::string("missing group end "
                                     "") +
                             compiler_grammar::group_end +
                             ""
                             " in pos=" +
                             std::to_string(pos));
      else
        return Op({e, pos + 1});
    }
  }
  ~compile_group()override{}

};


template <class Base_type,typename K>
std::map<K, std::unique_ptr<Base_type>>
make_unique_map(std::map<K,Base_type*>&& in)
{
  std::map<K, std::unique_ptr<Base_type>> out;
  for (auto& e: in)
    out[e.first]=std::unique_ptr<Base_type>(e.second);
  return out;
}




template <class UnOp, typename T, class Model>
struct compile_unary_operator_typed : public Base_Compiler<T, Model> {
  typedef Base_Compiler<T, Model> base_type;
  using typename base_type::Op;

  virtual Op compile(const std::string &s, std::size_t pos) const override {
    assert(s[pos] == UnOp::symbol);
    auto res = compile_term<T, Model>().compile(s, pos + 1);
    if (!res)
      return Op(false, res.error());
    else {
      auto [e, next] = std::move(res).value();
      return Op({new UnaryOperator<UnOp, T, Model>(e), next});
    }
  }

  virtual ~compile_unary_operator_typed() override{}
};

template <typename T, class Model> struct get_compile_unary_operator {

  struct minus_derivative {
    constexpr static auto const className = my_static_string("minus_op_derivative");
    inline constexpr static auto const symbol = compiler_grammar::minus_sign;
    Derivative<T> operator()(const Derivative<T>& x) const
    {
      return Derivative<T> (-x.f(),x.x(),-x.dfdx());
    }

  };

  struct minus {
    constexpr static auto const className = my_static_string("minus_op");
    inline constexpr static auto const symbol = compiler_grammar::minus_sign;
    T operator()(const T x) const { return -x; }
    typedef minus_derivative Derivative;
  };

  inline static const std::map<char, std::unique_ptr<Base_Compiler<T, Model>>>
      operator_map=make_unique_map<Base_Compiler<T,Model>,char>(std::map<char,Base_Compiler<T, Model>*>({{minus::symbol,
                                                new compile_unary_operator_typed<minus, T, Model>}}));

  myOptional_t<Base_Compiler<T, Model> const *>
  get_Operator(char symbol) const {
    typedef myOptional_t<Base_Compiler<T, Model> const *> Op;

    auto it = operator_map.find(symbol);
    if (it != operator_map.end()) {
      return Op(it->second.get());
    } else
      return Op(false, std::string("Operator "
                                   "") +
                           symbol +
                           ""
                           " is not a unary operator!!");
  }
};

template <typename T, class Model>
struct compile_unary : public Base_Compiler<T, Model> {
  typedef Base_Compiler<T, Model> base_type;
  using typename base_type::Op;
  ~compile_unary()override{}

  virtual Op compile(const std::string &s, std::size_t pos) const override {
    assert(compiler_grammar::unary_operators.find(s[pos]) != std::string::npos);
    auto symbol = s[pos];
    auto res = get_compile_unary_operator<T, Model>().get_Operator(symbol);
    if (!res)
      return Op(false, res.error());
    else {
      return res.value()->compile(s, pos);
    }
  }
};

template <typename T, class Model> struct get_compile_binary_operator {

  struct sum_derivative {
    constexpr static auto const className = my_static_string("sum_derivative");
    Derivative<double> operator()(const Derivative<double> &x,const Derivative<double> &y ) const {
      auto f = x.f()+y.f();
      auto dfdx = x.dfdx()+y.dfdx();
      return Derivative<double>(f, x.x(), std::move(dfdx));
    }
  };

  struct sum {
    constexpr static auto const className = my_static_string("sum_op");
    inline constexpr static auto const symbol = compiler_grammar::sum_operator;
    inline constexpr static const std::size_t precedence = 6;
    T operator()(const T &x, const T &y) const { return x + y; }
    typedef sum_derivative Derivative;
  };


  struct substraction_derivative {
    constexpr static auto const className = my_static_string("sum_derivative");
    Derivative<double> operator()(const Derivative<double> &x,
                                  const Derivative<double> &y) const {
      auto f = x.f() - y.f();
      auto dfdx = x.dfdx() - y.dfdx();
      return Derivative<double>(f, x.x(), std::move(dfdx));
    }
  };

  struct substraction {
    constexpr static auto const className = my_static_string("substraction_op");
    inline constexpr static auto const symbol = compiler_grammar::substraction_operator;
    inline constexpr static const std::size_t precedence = 6;
    T operator()(const T &x, const T &y) const { return x - y; }
    typedef substraction_derivative Derivative;
  };

  struct product_derivative {
    constexpr static auto const className = my_static_string("product_derivative");
    Derivative<double> operator()(const Derivative<double> &x,
                                  const Derivative<double> &y) const {
      auto f = x.f() * y.f();
      auto dfdx = x.dfdx().apply([&y](auto &m){ return m*y.f();}) + y.dfdx().apply([&x](auto &m){ return x.f()*m;});
      return Derivative<double>(f, x.x(), std::move(dfdx));
    }
  };

  struct product {
    constexpr static auto const className = my_static_string("product_op");
    inline constexpr static auto const symbol = compiler_grammar::prod_operator;
    inline constexpr static const std::size_t precedence = 6;
    T operator()(const T &x, const T &y) const { return x * y; }
    typedef product_derivative Derivative;
  };

  struct division_derivative {
    constexpr static auto const className =
        my_static_string("product_derivative");
    Derivative<double> operator()(const Derivative<double> &x,
                                  const Derivative<double> &y) const {
      auto f = x.f() / y.f();
      auto dfdx = x.dfdx().apply([&y](auto &dx) { return dx / y.f(); }) +
                  y.dfdx().apply([&f,&x,&y](auto &dy) { return f / y.f() * dy; });
      return Derivative<double>(f, x.x(), std::move(dfdx));
    }
  };

  struct division {
    constexpr static auto const className = my_static_string("div_op");
    inline constexpr static auto const symbol = compiler_grammar::div_operator;
    inline constexpr static const std::size_t precedence = 6;
    T operator()(const T &x, const T &y) const { return x / y; }
    typedef division_derivative Derivative;
  };

  struct power_derivative {
    constexpr static auto const className =
        my_static_string("product_derivative");
    Derivative<double> operator()(const Derivative<double> &x,
                                  const Derivative<double> &y) const {
      auto f = std::pow(x.f() ,y.f());
      auto dfdx = x.dfdx().apply([&f,&x,&y](auto &dx) {
        return y.f()*std::pow(x.f() ,y.f()-1.0)*dx;
      }) + y.dfdx().apply([&f, &x, &y](auto &dy) { return std::log(x.f())*f * dy; });
      return Derivative<double>(f, x.x(), std::move(dfdx));
    }
  };

  struct power {
    constexpr static auto const className = my_static_string("div_op");
    inline constexpr static auto const symbol = compiler_grammar::power_operator;
    inline constexpr static const std::size_t precedence = 6;
    T operator()(const T &x, const T &y) const { return std::pow(x ,y); }
    typedef division_derivative Derivative;
  };

  inline static const std::map<char, std::unique_ptr<Base_Binary_Compiler<T, Model>>>
      operator_map = make_unique_map(std::map<char,Base_Binary_Compiler<T, Model>*>({
          {sum::symbol,new compile_binary_operator_typed<sum, T, Model>},
          {substraction::symbol,new compile_binary_operator_typed<substraction, T, Model>},
          {product::symbol,new compile_binary_operator_typed<product, T, Model>},
          {division::symbol,new compile_binary_operator_typed<division, T, Model>},

      }));

  myOptional_t<Base_Binary_Compiler<T, Model> const *>
  get_Operator(char symbol) const {
    typedef myOptional_t<Base_Binary_Compiler<T, Model> const *> Op;

    auto it = operator_map.find(symbol);
    if (it != operator_map.end()) {
      return Op(it->second.get());
    } else
      return Op(false, std::string("Operator "
                                   "") +
                           symbol +
                           ""
                           " is not a binary operator!!");
  }
};

template <class BinOp, typename T, class Model>
struct compile_binary_operator_typed : public Base_Binary_Compiler<T, Model> {
  typedef Base_Binary_Compiler<T, Model> base_type;
  using typename base_type::Op;
  ~compile_binary_operator_typed()override{}

  virtual std::size_t operator_precedence() const override {
    return BinOp::precedence;
  }

  Op compile(Base_Function<T, Model> *first_term,
             Base_Function<T, Model> *second_term, const std::string &s,
             std::size_t pos) const {
    if (!(pos < s.size() && compiler_grammar::binary_operators.find(s[pos]))) {
      return Op(
          {new BinaryOperator<BinOp, T, Model>(first_term, second_term), pos});

    } else {
      auto Op_2_ = get_compile_binary_operator<T, Model>().get_Operator(s[pos]);
      if (!Op_2_)
        return Op(false, Op_2_.error());
      else {
        auto op_2 = std::move(Op_2_).value();
        if (operator_precedence() <= op_2->operator_precedence()) {
          auto new_first_term =
              new BinaryOperator<BinOp, T, Model>(first_term, second_term);
          return op_2->compile(new_first_term, s, pos);
        } else {
          auto new_second_term_ = op_2->compile(second_term, s, pos);
          if (!new_second_term_)
            return Op(false, new_second_term_.error());
          else {
            auto [new_second_term, next] = std::move(new_second_term_).value();
            return compile(first_term, new_second_term, s, next);
          }
        }
      }
    }
  }

  virtual Op compile(Base_Function<T, Model> *first_term, const std::string s,
                     std::size_t pos) const override {
    assert(s[pos] == BinOp::symbol);
    auto res = compile_term<T, Model>().compile(s, pos + 1);
    if (!res)
      return Op(false, std::string("Error in second term in ") + BinOp::symbol +
                           " operator :" + res.error());
    else {
      auto [second_term, next] = std::move(res).value();
      return compile(first_term, second_term, s, next);
    }
  }


  // Base_Binary_Compiler interface
};

template <typename T, class Model>
struct compile_term : public Base_Compiler<T, Model> {
  typedef Base_Compiler<T, Model> base_type;
  using typename base_type::Op;
  ~compile_term()override{}


  virtual Op compile(const std::string &s, std::size_t pos0) const override {
    std::size_t pos = s.find_first_not_of(compiler_grammar::space,pos0);
    if (compiler_grammar::alfa.find(s[pos]) != compiler_grammar::alfa.npos) {
      return compile_identifier<T, Model>().compile(s, pos);
    } else if ((compiler_grammar::number.find(s[pos]) != std::string::npos) ||
               ((compiler_grammar::sign.find(s[pos]) != compiler_grammar::sign.npos) &&
                (compiler_grammar::number.find(s[pos + 1]) != std::string::npos))) {
      return compile_number<T, Model>().compile(s, pos);
    } else if (s[pos] == compiler_grammar::group_start)
      return compile_group<T, Model>().compile(s, pos);
    else if (compiler_grammar::unary_operators.find(s[pos]) != compiler_grammar::unary_operators.npos)
      return compile_unary<T, Model>().compile(s, pos);
    else {
      return Op(false, std::string("invalid term : "
                                   "") +
                           s[pos] +
                           ""
                           " is an invalid term start, pos=" +
                           std::to_string(pos) + "in " + s);
    }
  }
};

template <typename T, class Model>
struct compile_expression : public Base_Compiler<T, Model> {
  typedef Base_Compiler<T, Model> base_type;
  using typename base_type::Op;

  ~compile_expression()override {}
  virtual Op compile(const std::string &s, std::size_t pos) const override {
    auto term0 = compile_term<T, Model>().compile(s, pos);
    if (!term0)
      return Op(false, term0.error());
    else {
      auto next = term0.value().second;
      if (!(next < s.size() && compiler_grammar::binary_operators.find(s[next]))) {
        return term0;
      } else {
        auto BinOp_ =
            get_compile_binary_operator<T, Model>().get_Operator(s[next]);
        if (!BinOp_)
          return Op(false, BinOp_.error());
        else {
          return BinOp_.value()->compile(term0.value().first, s, next);
        }
      }
    }
  }
};



template <typename T, class Model, typename K>
auto compile_map(const std::map<K,std::string>& m)
{
  typedef myOptional_t<std::map<K,Base_Function<T,Model>*>> Op;
  std::map<K,Base_Function<T,Model>*> out;
  compile_expression<T,Model> c;
  bool succeed=true;
  std::string errors;
  for (auto& e: m)
  {
    auto res=c.compile(e.second,0);
    if (!res)
    {
      succeed=false; errors+="\n error in"+ToString(e.first)+res.error();
    }
    else
      out[e.first]=std::move(res).value().first;
  }
  if (succeed)
    return Op(out);
  else return Op(false,errors);
 }

 template <typename T, class Model, typename K>
 auto text_map(const std::map<K,std::unique_ptr<Base_Function<T,Model>>>& m)
 {
   std::map<K,std::string> out;
   for (auto& e: m)
   {

     out[e.first]=ToString(*e.second);
   }
   return out;

 }
#endif // MYDYAMICFUNCTIONS_HPP
