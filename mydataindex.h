#ifndef MYDATAINDEX_H
#define MYDATAINDEX_H

template <class... Fields> struct CT_s {
  constexpr static auto const title =
      (my_static_string("") + ... + Fields::title);
};

template <std::size_t N, class SuperFunction, class SubFunction> struct CF_s {
private:
  template <class Object, std::size_t M, std::size_t... Supers,
            std::size_t... Subs>

  std::invoke_result_t<
      SubFunction,
      std::invoke_result_t<SuperFunction, Object, decltype(Supers)...>,
      decltype(Subs)...>

  op_impl(const Object &x, std::array<std::size_t, M> i,
          std::index_sequence<Supers...>, std::index_sequence<Subs...>) const {
    return std::invoke(sub_, std::invoke(f_, x, std::get<Supers>(i)...),
                       std::get<N + Subs>(i)...);
  }

public:
  CF_s(const SuperFunction &f, SubFunction &&sub)
      : f_(f), sub_{std::move(sub)} {}

  template <class Object, typename... Int>
  auto operator()(const Object &x, Int... is) const
      -> decltype(this->template op_impl(
          x, std::declval<std::array<std::size_t, sizeof...(Int)>>(),
          std::make_index_sequence<N>(),
          std::make_index_sequence<sizeof...(Int) - N>()))

  {
    constexpr auto M = sizeof...(Int);
    return op_impl(x, std::array<std::size_t, M>{is...},
                   std::make_index_sequence<N>(),
                   std::make_index_sequence<M - N>());
  }

  SuperFunction f_;
  SubFunction sub_;
};

template <class Field, class Fun> class F_s {
public:
  F_s(Field, Fun &&f) : f_{f} {}

  template <class Object, typename... Int>
  std::invoke_result_t<Fun, Object, Int...> operator()(const Object &x,
                                                       Int... i) const {
    return std::invoke(f_, x, i...);
  }

  static auto title() { return Field::title(); }

  auto &f() { return f_; }
  auto &f() const { return f_; }

private:
  Fun f_;
};

template <class Index, class Size, class Fun> class I_s {
public:
  I_s(Index, Size &&s, Fun &&f) : s_{std::move(s)}, f_{std::move(f)} {}
  I_s(Index, Size &&s, const Fun &f) : s_{std::move(s)}, f_{f} {}

  template <class Object, typename... Int>
  std::size_t size(const Object &x, Int... i) const {
    return std::invoke(s_, x, i...);
  }
  template <class Object, typename... Int>
  auto operator()(const Object &x, Int... i) const {
    return std::invoke(f_, x, i...);
  }

  static auto title() { return Index::title(); }

  auto &f() { return f_; }
  auto &f() const { return f_; }
  auto &s() { return s_; }
  auto &s() const { return s_; }

private:
  Size s_;
  Fun f_;
};

template <class Index, class Size> class I_s<Index, Size, int> {
public:
  I_s(Index, Size &&s, int = 0) : s_{s} {}

  template <class Object, typename... Int>
  std::size_t size(const Object &x, Int... i) const {
    return std::invoke(s_, x, i...);
  }
  template <class Object, typename... Int>
  std::size_t operator()(const Object &, Int... is) const {
    return std::array<std::size_t, sizeof...(Int)>{is...}.back();
  }

  static auto title() { return Index::title(); }

  auto &s() { return s_; }
  auto &s() const { return s_; }

private:
  Size s_;
};

template <class Index, class Size, class Fun, class Size2, class Fun2>
I_s<Index, Size, Fun> chooseWithFun(I_s<Index, Size, Fun> &&one,
                                    I_s<Index, Size2, Fun2> &&) {
  return std::move(one);
}

template <class Index, class Size, class Fun, class Size2, class Fun2>
auto chooseWithFun(I_s<Index, Size, int> &&, I_s<Index, Size2, Fun2> &&two)
    -> std::enable_if<!std::is_same_v<Fun2, int>, I_s<Index, Size2, Fun2>> {
  return std::move(two);
}

template <std::size_t Parent_I, class ParentField, class ParentFun, class Field,
          class Fun>
auto Compose_fun(F_s<ParentField, ParentFun> pfun, F_s<Field, Fun> &&f) {
  return F_s<CT_s<ParentField, Field>, CF_s<Parent_I, ParentFun, Fun>>(
      CT_s<ParentField, Field>{},
      CF_s<Parent_I, ParentFun, Fun>(pfun.f(), std::move(f.f())));
}

template <std::size_t Parent_I, class ParentFun, class Field, class Fun>
auto Compose_fun(ParentFun pfun, F_s<Field, Fun> &&f) {
  return F_s<Field, CF_s<Parent_I, ParentFun, Fun>>(
      Field{}, CF_s<Parent_I, ParentFun, Fun>(pfun, std::move(f.f())));
}

template <std::size_t Parent_I, class ParentFun, class Index, class Size>
auto Compose_fun_Index(const ParentFun &pfun, I_s<Index, Size, int> &&index) {
  return I_s<Index, CF_s<Parent_I, ParentFun, Size>, int>(
      Index{}, CF_s<Parent_I, ParentFun, Size>(pfun, std::move(index.s())), 0);
}

template <std::size_t Parent_I, class ParentFun, class Index, class Size,
          class Fun>
auto Compose_fun_Index(const ParentFun &pfun, I_s<Index, Size, Fun> &&index) {
  return I_s<Index, CF_s<Parent_I, ParentFun, Size>,
             CF_s<Parent_I, ParentFun, Fun>>(
      Index{}, CF_s<Parent_I, ParentFun, Size>(pfun, std::move(index.s())),
      CF_s<Parent_I, ParentFun, Fun>(pfun, std::move(index.f())));
}

template <class Indexes, class Fun> class Data_Index_static;

template <class... Index, class... Size, class... IndexFun, class... Field,
          class... Fun>
class Data_Index_static<Cs<I_s<Index, Size, IndexFun>...>,
                        Cs<F_s<Field, Fun>...>> {
public:
  static constexpr std::size_t N = sizeof...(Index);

  Data_Index_static(std::tuple<I_s<Index, Size, IndexFun>...> &&Ind,
                    std::tuple<F_s<Field, Fun>...> &&ftu)
      : Ind_{std::move(Ind)}, ftu_{std::move(ftu)} {}

  static std::ostream &put_title(std::ostream &os, const std::string &sep) {
    ((os << Index::title.str() << sep), ...);
    ((os << Field::title.str() << sep), ...);
    return os << "\n";
  }

  static bool print_title(const std::string &filename, const std::string &ext,
                          const std::string &sep) {
    std::string fname = filename + get_index_names() + ext;
    std::ofstream f(fname.c_str());
    put_title(f, sep);
    return true;
  }
  template <class Object>
  bool print_data(const Object &x, const std::string &filename,
                  const std::string &ext, const std::string &sep) const {
    std::stringstream ss;
    put_index_sample(ss, x, sep);
    std::string fname = filename + get_index_names() + ext;
    std::ofstream f(fname.c_str(), std::ios_base::app);
    f << ss.str();
    f.close();
    return true;
  }

  template <class Object>
  std::ostream &put_index_sample(std::ostream &os, const Object &x,
                                 const std::string &sep) const {
    return put_index_sample_impl(os, x, sep, std::tuple<>(),
                                 std::index_sequence_for<Index...>());
  }

  static std::string get_index_names() {
    return (Index::title.str() + ... + "");
  }
  auto &Indexes() { return Ind_; }
  auto &Indexes() const { return Ind_; }
  auto &Data() { return ftu_; }
  auto &Data() const { return ftu_; }

private:
  template <class Object, class... Index_Type, typename... Int>
  std::ostream &put_index_sample_impl(std::ostream &os, const Object &x,
                                      const std::string &sep,
                                      std::tuple<Index_Type...> out,
                                      std::index_sequence<>, Int... is) const {
    std::apply([&os, &sep](auto &... f) { ((os << f << sep), ...); }, out);
    std::apply([&x, &os, &sep,
                is...](auto &... f) { ((os << f(x, is...) << sep), ...); },
               ftu_);
    return os << "\n";
  }

  template <std::size_t I, std::size_t... Is, class Object, class... Index_Type,
            typename... Int>
  std::ostream &
  put_index_sample_impl(std::ostream &os, const Object &x,
                        const std::string &sep, std::tuple<Index_Type...> out,
                        std::index_sequence<I, Is...>, Int... is) const {
    for (std::size_t i = 0; i < Index_size<I>(x, is...); ++i)
      put_index_sample_impl(
          os, x, sep,
          std::tuple_cat(out, std::make_tuple(std::get<I>(Ind_)(x, is..., i))),
          std::index_sequence<Is...>(), is..., i);
    return os;
  }
  template <std::size_t I, class Object, typename... Int>
  std::size_t Index_size(const Object &x, Int... is) const {
    return std::get<I>(Ind_).size(x, is...);
  }
  std::tuple<I_s<Index, Size, IndexFun>...> Ind_;
  std::tuple<F_s<Field, Fun>...> ftu_;
};

template <class... Index, class... Size, class... FunIndex, class... Field,
          class... Fun>
Data_Index_static<Cs<I_s<Index, Size, FunIndex>...>, Cs<F_s<Field, Fun>...>>
make_data_static(std::tuple<I_s<Index, Size, FunIndex>...> &&Ind,
                 std::tuple<F_s<Field, Fun>...> &&ftu) {
  return Data_Index_static<Cs<I_s<Index, Size, FunIndex>...>,
                           Cs<F_s<Field, Fun>...>>(std::move(Ind),
                                                   std::move(ftu));
}

// template<class... NewIndex, class... NewSize,class ...FunIndex,class PF,
// class PFun, class I, class D> auto Compose_static(std::tuple<I_s<NewIndex,
// NewSize,FunIndex>...> s,
//                          F_s<PF,PFun> f,
//                    Data_Index_static<I,D> && d)
//{
//    constexpr auto N=sizeof...(NewIndex);
//    return
//        make_data_static(std::tuple_cat(std::move(s),
//                                         std::apply([&f](auto&&... ind)
//                                                    { return
//                                                    std::make_tuple(Compose_fun_Index<N>(f.f(),std::move(ind))...);},
//                                                    std::move(d.Indexes()))),
//                          std::apply([&f](auto&&...tu){return
//                          std::make_tuple(Compose_fun<N>(f,std::move(tu))...);},std::move(d.Data()))
//            );
//}

template <class... NewIndex, class... NewSize, class... FunIndex, class PFun,
          class I, class D>
auto Compose_static(std::tuple<I_s<NewIndex, NewSize, FunIndex>...> s, PFun f,
                    Data_Index_static<I, D> &&d) {
  constexpr auto N = sizeof...(NewIndex);
  return make_data_static(
      std::tuple_cat(std::move(s),
                     std::apply(
                         [&f](auto &&... ind) {
                           return std::make_tuple(
                               Compose_fun_Index<N>(f, std::move(ind))...);
                         },
                         std::move(d.Indexes()))),
      std::apply(
          [&f](auto &&... tu) {
            return std::make_tuple(Compose_fun<N>(f, std::move(tu))...);
          },
          std::move(d.Data())));
}

template <class PF, class PFun, class I, class D>
auto Compose_static(F_s<PF, PFun> f, Data_Index_static<I, D> &&d) {
  return make_data_static(
      std::apply(
          [&f](auto &&... ind) {
            return std::make_tuple(Compose_fun_Index<0>(f, std::move(ind))...);
          },
          std::move(d.Indexes())),
      std::apply(
          [&f](auto &&... tu) {
            return std::make_tuple(Compose_fun<0>(f, std::move(tu))...);
          },
          std::move(d.Data())));
}

template <class PFun, class I, class D>
auto Compose_static(PFun f, Data_Index_static<I, D> &&d) {
  return make_data_static(
      std::apply(
          [&f](auto &&... ind) {
            return std::make_tuple(Compose_fun_Index<0>(f, std::move(ind))...);
          },
          std::move(d.Indexes())),
      std::apply(
          [&f](auto &&... tu) {
            return std::make_tuple(Compose_fun<0>(f, std::move(tu))...);
          },
          std::move(d.Data())));
}

// template<class... NewIndex, class... NewSize,class... FunIndex,class PF,
// class PFun,class... Data_Index_static_R> auto
// Compose_static(std::tuple<I_s<NewIndex, NewSize,FunIndex>...> s,
//                    F_s<PF,PFun> f,
//                    std::tuple<Data_Index_static_R...> && d)
//{
//    return std::apply([&s,&f](auto&&...D){ return
//    std::make_tuple(Compose_static(s,f,std::move(D))...);},std::move(d));

//}

template <class... NewIndex, class... NewSize, class... FunIndex, class PFun,
          class... Data_Index_static_R>
auto Compose_static(std::tuple<I_s<NewIndex, NewSize, FunIndex>...> s, PFun f,
                    std::tuple<Data_Index_static_R...> &&d) {
  return std::apply(
      [&s, &f](auto &&... D) {
        return std::make_tuple(Compose_static(s, f, std::move(D))...);
      },
      std::move(d));
}

template <class PF, class PFun, class... Data_Index_static_R>
auto Compose_static(F_s<PF, PFun> f, std::tuple<Data_Index_static_R...> &&d) {
  return std::apply(
      [&f](auto &&... D) {
        return std::make_tuple(Compose_static(f, std::move(D))...);
      },
      std::move(d));
}

template <class PFun, class... Data_Index_static_R>
auto Compose_static(PFun f, std::tuple<Data_Index_static_R...> &&d) {
  return std::apply(
      [&f](auto &&... D) {
        return std::make_tuple(Compose_static(f, std::move(D))...);
      },
      std::move(d));
}

template <std::size_t... Is, class... Index, class... Size1, class... IFun1,
          class... Size2, class... IFun2>
auto chooseWithFun(std::tuple<I_s<Index, Size1, IFun1>...> &&one,
                   std::tuple<I_s<Index, Size2, IFun2>...> &&two,
                   std::index_sequence<Is...>) {
  return std::make_tuple(chooseWithFun(std::move(std::get<Is>(one)),
                                       std::move(std::get<Is>(two)))...);
}

template <class Index1, class Index2, class Size1, class IFun1, class Size2,
          class IFun2>
auto fuseWithFun_Is(I_s<Index1, Size1, IFun1> &&one,
                    I_s<Index2, Size2, IFun2> &&two)
    -> std::enable_if_t<
        !std::is_same_v<Index1, Index2> ||
            (std::is_same_v<IFun1, int> == std::is_same_v<IFun2, int>),
        std::pair<I_s<Index1, Size1, IFun1>, I_s<Index2, Size2, IFun2>>>

{

  return std::make_pair(std::move(one), std::move(two));
}

template <class Index, class Size1, class IFun, class Size2>
auto fuseWithFun_Is(I_s<Index, Size1, IFun> &&one, I_s<Index, Size2, int> &&two)
    -> std::enable_if_t<
        !std::is_same_v<IFun, int>,
        std::pair<I_s<Index, Size1, IFun>, I_s<Index, Size2, IFun>>> {

  auto newtwo = I_s(Index{}, std::move(two.s()), one.f());
  return std::make_pair(std::move(one), std::move(newtwo));
}

template <class Index, class Size1, class IFun, class Size2>
auto fuseWithFun_Is(I_s<Index, Size1, int> &&one, I_s<Index, Size2, IFun> &&two)
    -> std::enable_if_t<
        !std::is_same_v<IFun, int>,
        std::pair<I_s<Index, Size1, IFun>, I_s<Index, Size2, IFun>>> {

  // typedef typename IFun::df IFUN;
  auto newone = I_s(Index{}, std::move(one.s()), two.f());
  return std::make_pair(std::move(newone), std::move(two));
}

template <class... D01, class... D02, class... D11, class... D12,
          std::size_t... I1, std::size_t... I2, std::size_t I>
auto fuseWithFun_tuple(std::tuple<D01...> &&in1, std::tuple<D11...> &&out1,
                       std::index_sequence<I, I1...>,

                       std::tuple<D02...> &&in2, std::tuple<D12...> &&out2,
                       std::index_sequence<I, I2...>) {

  auto [ind1, ind2] =
      fuseWithFun_Is(std::move(std::get<I>(in1)), std::move(std::get<I>(in2)));
  return fuseWithFun_tuple(
      std::move(in1),
      std::tuple_cat(std::move(out1), std::make_tuple(std::move(ind1))),
      std::index_sequence<I1...>(), std::move(in2),
      std::tuple_cat(std::move(out2), std::make_tuple(std::move(ind2))),
      std::index_sequence<I2...>());
}

template <class... D01, class... D02, class... D11, class... D12,
          std::size_t... I1, std::size_t I>
auto fuseWithFun_tuple(std::tuple<D01...> &&in1, std::tuple<D11...> &&out1,
                       std::index_sequence<I, I1...>, std::tuple<D02...> &&,
                       std::tuple<D12...> &&out2, std::index_sequence<>) {
  return std::make_pair(
      std::tuple_cat(std::move(out1),
                     std::make_tuple(std::move(std::get<I>(in1)),
                                     std::move(std::get<I1>(in1))...)),
      std::move(out2));
}

template <class... D01, class... D02, class... D11, class... D12,
          std::size_t... I2>
auto fuseWithFun_tuple(std::tuple<D01...> &&, std::tuple<D11...> &&out1,
                       std::index_sequence<>, std::tuple<D02...> &&in2,
                       std::tuple<D12...> &&out2, std::index_sequence<I2...>) {
  return std::make_pair(
      std::move(out1),
      std::tuple_cat(std::move(out2),
                     std::make_tuple(std::move(std::get<I2>(in2))...)));
}

template <class... D1, class... D2>
auto fuseWithFun_tuple(std::tuple<D1...> &&tu1, std::tuple<D2...> &&tu2) {
  return fuseWithFun_tuple(std::move(tu1), std::tuple<>(),
                           std::index_sequence_for<D1...>(), std::move(tu2),
                           std::tuple<>(), std::index_sequence_for<D2...>());
}

template <bool TF, class... D0, class... D1, class... Index1, class... Size1,
          class... IFun1, class... Field1, class... Fun1, class... Index2,
          class... Size2, class... IFun2, class... Field2, class... Fun2,
          class... R>
auto Concatenate_tuple_static_impl(
    std::integral_constant<
        std::enable_if_t<(!std::is_same_v<Cs<Index1...>, Cs<Index2...>>), bool>,
        TF>,
    std::tuple<D0...> &&tu, std::tuple<D1...> &&tu1,
    Data_Index_static<Cs<I_s<Index1, Size1, IFun1>...>,
                      Cs<F_s<Field1, Fun1>...>> &&e,
    Data_Index_static<Cs<I_s<Index2, Size2, IFun2>...>,
                      Cs<F_s<Field2, Fun2>...>> &&eTu,
    R &&... r) {
  if constexpr (true) {

    auto t =
        fuseWithFun_tuple(std::move(e.Indexes()), std::move(eTu.Indexes()));

    auto new_e = make_data_static(std::move(t.first), std::move(e.Data()));
    auto new_eTu = make_data_static(std::move(t.second), std::move(eTu.Data()));

    //  typedef typename decltype(e)::e  old_e;
    //    typedef typename decltype(new_e)::newe  new_e_e;

    if constexpr (TF) {
      return Concatenate_tuple_static_impl(
          std::true_type{}, std::move(tu),
          std::tuple_cat(std::move(tu1), std::make_tuple(std::move(new_eTu))),
          std::move(new_e), std::move(r)...);
    } else {
      static_assert(sizeof...(D1) == 0);
      return Concatenate_tuple_static_impl(
          std::false_type{},
          std::tuple_cat(std::move(tu), std::make_tuple(std::move(new_eTu))),
          std::move(tu1), std::move(new_e), std::move(r)...);
    }
  } else {
    if constexpr (TF) {
      return Concatenate_tuple_static_impl(
          std::true_type{}, std::move(tu),
          std::tuple_cat(std::move(tu1), std::make_tuple(std::move(eTu))),
          std::move(e), std::move(r)...);
    } else {
      static_assert(sizeof...(D1) == 0);
      return Concatenate_tuple_static_impl(
          std::false_type{},
          std::tuple_cat(std::move(tu), std::make_tuple(std::move(eTu))),
          std::move(tu1), std::move(e), std::move(r)...);
    }
  }
}

template <class... D, class... Index, class... Size1, class... IFun1,
          class... Field1, class... Fun1, class... Size2, class... IFun2,
          class... Field2, class... Fun2, class... R>
auto Concatenate_tuple_static_impl(
    std::false_type, std::tuple<D...> &&tu, std::tuple<>,
    Data_Index_static<Cs<I_s<Index, Size2, IFun2>...>, Cs<F_s<Field2, Fun2>...>>
        &&e,
    Data_Index_static<Cs<I_s<Index, Size1, IFun1>...>, Cs<F_s<Field1, Fun1>...>>
        &&eTu,
    R &&... r) {

  auto new_e = make_data_static(
      chooseWithFun(std::move(eTu.Indexes()), std::move(e.Indexes()),
                    std::index_sequence_for<Index...>()),
      std::tuple_cat(std::move(eTu.Data()), std::move(e.Data())));

  return Concatenate_tuple_static_impl(std::true_type{}, std::move(tu),
                                       std::tuple<>(), std::move(new_e),
                                       std::move(r)...);
}

template <class... D, class... Index, class... Size, class... IFun,
          class... Field, class... Fun>
auto Concatenate_tuple_static_impl(
    std::false_type, std::tuple<D...> &&tu, std::tuple<>,
    Data_Index_static<Cs<I_s<Index, Size, IFun>...>, Cs<F_s<Field, Fun>...>>
        &&e) {
  return std::tuple_cat(std::move(tu), std::make_tuple(std::move(e)));
}
template <class... D0, class... D1, class... Index, class... Size,
          class... IFun, class... Field, class... Fun>
auto Concatenate_tuple_static_impl(
    std::true_type, std::tuple<D0...> &&tu, std::tuple<D1...> &&tu1,
    Data_Index_static<Cs<I_s<Index, Size, IFun>...>, Cs<F_s<Field, Fun>...>>
        &&e) {
  return std::tuple_cat(std::move(tu), std::make_tuple(std::move(e)),
                        std::move(tu1));
}

template <class D, class... Datas, class... Index, class... Size, class... IFun,
          class... Field, class... Fun>

auto Concatenate_tuple_static(std::tuple<D, Datas...> &&tu,
                              Data_Index_static<Cs<I_s<Index, Size, IFun>...>,
                                                Cs<F_s<Field, Fun>...>> &&e) {
  auto out = std::apply(
      [&e](auto &&... t) {
        return Concatenate_tuple_static_impl(std::false_type{}, std::tuple<>(),
                                             std::tuple<>(), std::move(e),
                                             std::move(t)...);
      },
      std::move(tu));

  // typedef typename decltype(out)::one one;

  return out;
}

template <class... Ds, class D2, class D3, class... R>
auto Concatenate_tuple_static(std::tuple<Ds...> &&one, D2 &&r, D3 &&s,
                              R &&... rest) {
  auto next = Concatenate_tuple_static(std::move(one), std::move(r));
  return Concatenate_tuple_static(std::move(next), std::move(s),
                                  std::move(rest)...);
}

template <class... Datas, class... Datas2>
auto Concatenate_tuple_static(std::tuple<Datas...> &&one,
                              std::tuple<Datas2...> &&two) {
  return std::apply(
      [&one](auto &&... d) {
        return Concatenate_tuple_static(std::move(one), std::move(d)...);
      },
      std::move(two));
}
template <class... Datas>
auto Concatenate_tuple_static(std::tuple<Datas...> &&one) {
  return std::move(one);
}

class Data_Index_point {
  std::map<std::string, std::pair<std::size_t, std::size_t>> pos_;

public:
  void push_back(std::string &&name, std::size_t n) {
    pos_[std::move(name)] = std::pair(0ul, n);
  }

  Data_Index_point begin() {
    Data_Index_point out(*this);
    for (auto &e : out.pos_)
      e.second.first = 0;
    return out;
  }

  bool isEnd() const {
    return !(pos_.rbegin()->second.first < pos_.rbegin()->second.second);
  }

  Data_Index_point &operator++() {
    auto it = pos_.begin();
    while (true) {
      ++(it->second.first);
      if (it->second.first < it->second.second)
        return *this;
      else {
        auto it2 = it;
        it2++;
        if ((it2) != pos_.end()) {
          it->second.first = 0;
          ++it;
        } else
          return *this;
      }
    }
  }
};

class Data_Index_scheme {
  std::map<std::set<std::string>, std::vector<std::string>> d_;

  static std::set<std::string> intersects(const std::set<std::string> &one,
                                          std::set<std::string> &two) {
    std::set<std::string> out;
    std::set_intersection(one.begin(), one.end(), two.begin(), two.end(),
                          std::inserter(out, out.begin()));
    return out;
  }

  static std::set<std::string> set_difference(const std::set<std::string> &one,
                                              std::set<std::string> &two) {
    std::set<std::string> out;
    std::set_difference(one.begin(), one.end(), two.begin(), two.end(),
                        std::inserter(out, out.begin()));
    return out;
  }

  template <class inputIter>
  static std::set<std::string> intersects(inputIter it, inputIter end) {
    auto init = it->first;
    ++it;
    for (; it != end; ++it) {
      init = intersects(init, it->first);
    }
    return init;
  }

public:
  auto &self_map() const { return d_; }

  std::set<std::string> first_index() const { return d_.begin()->first; }

  typedef Data_Index_scheme self_type;

  void push_back(std::string &&name, std::set<std::string> &&new_indexes) {
    d_[new_indexes].push_back(std::move(name));
  }

  friend Data_Index_scheme concatenate_impl(Data_Index_scheme &&one,
                                            const Data_Index_scheme &two) {
    for (auto &e : two.d_) {
      auto &p = one.d_[e.first];
      p.insert(p.end(), e.second.begin(), e.second.end());
    }
    return std::move(one);
  }

  friend Data_Index_scheme concatenate(Data_Index_scheme &&one) {
    return std::move(one);
  }
  template <class... Data_Index>
  friend Data_Index_scheme concatenate(Data_Index_scheme &&one,
                                       const Data_Index_scheme &two,
                                       const Data_Index &... three) {
    return concatenate(concatenate_impl(std::move(one), two), three...);
  }

  template <class... T> static Data_Index_scheme data_index();

  std::size_t size() const { return d_.size(); }

  std::vector<std::string> set_titles() const {
    std::vector<std::string> out(d_.size());

    std::size_t i = 0;
    for (auto &e : d_) {

      for (auto &in : e.first)
        out[i] += "_" + in;
      ++i;
    }
    return out;
  }

  std::vector<std::vector<std::string>> index_titles() const {
    std::vector<std::vector<std::string>> out(size());
    std::size_t i = 0;
    for (auto it = d_.begin(); it != d_.end(); ++it, ++i) {
      std::vector<std::string> e;
      e.insert(e.end(), it->first.begin(), it->first.end());
      out[i] = e;
    }
    return out;
  }

  std::vector<std::vector<std::string>> names_titles() const {
    std::vector<std::vector<std::string>> out(size());
    std::size_t i = 0;
    for (auto it = d_.begin(); it != d_.end(); ++it, ++i) {
      out[i] = it->second;
    }
    return out;
  }

 inline void insert_index(std::string &&new_index_name) {
    std::map<std::set<std::string>, std::vector<std::string>> new_d;
    for (auto &e : d_) {
      auto i = e.first;
      i.insert(new_index_name);
      new_d[std::move(i)] = std::move(e.second);
    }
    d_ = std::move(new_d);
  }
inline  Data_Index_scheme &
  insert_indexes(const std::vector<std::string> &new_index_name) {
    std::map<std::set<std::string>, std::vector<std::string>> new_d;
    for (auto &e : d_) {
      auto i = e.first;
      i.insert(new_index_name.begin(), new_index_name.end());
      new_d[std::move(i)] = std::move(e.second);
    }
    d_ = std::move(new_d);
    return *this;
  }

  Data_Index_scheme
  operator*(std::vector<std::string> &&new_index_names) const & {
    Data_Index_scheme out(*this);
    return out.insert_indexes(std::move(new_index_names));
  }

  Data_Index_scheme operator*(std::vector<std::string> &&new_index_names) && {
    insert_indexes(std::move(new_index_names));
    return *this;
  }

  /*
  static auto get_constructor_fields()
  {
      return std::make_tuple(
          grammar::field(C<self_type>{},"value_name",&self_type::value_name),
          grammar::field(C<self_type>{},"index_names",&self_type::index_names));
  }


  Data_Index_scheme(const std::string& name, const std::vector<std::string>&
  indexes) :value_name_{name},index_names_{indexes}{}

  Data_Index_scheme( std::string&& name,  std::vector<std::string>&& indexes)
      :value_name_{std::move(name)},index_names_{std::move(indexes)}{}

  Data_Index_scheme()=default;

*/
};

inline Data_Index_scheme Insert_Index(Data_Index_scheme &&v,
                               std::vector<std::string> &&indexes) {
  v.insert_indexes(std::move(indexes));
  return std::move(v);
}

template <class...> struct myData_Index;

template <class T> struct myData_Index<T> {
  static Data_Index_scheme data_index() { return T::data_index(); }

  template <class... Indexes> static auto get_data_index_static(Indexes... i) {
    return T::get_data_index_static(i...);
  }
};

template <class... T> struct myData_Index<std::tuple<T...>> {
  typedef std::tuple<T...> self_type;

  template <std::size_t... Is>
  static auto get_data_index_static_impl(std::index_sequence<Is...>) {
    return Concatenate_tuple_static(Compose_static(
        [](const self_type &self) { return std::get<Is>(self); },
        myData_Index<
            std::tuple_element_t<Is, self_type>>::get_data_index_static())...);
  }

  template <std::size_t... Is, class index>
  static auto get_data_index_static_impl(std::index_sequence<Is...>,
                                         index ind) {
    auto out = Concatenate_tuple_static(
        Compose_static([](const self_type &self) -> auto const & {
          return std::get<Is>(self);
        },
                       myData_Index<std::tuple_element_t<Is, self_type>>::
                           get_data_index_static(ind))...);

    // typedef typename decltype(out)::check test;
    return out;
  }

  template <class... indexes> static auto get_data_index_static(indexes... i) {
    return get_data_index_static_impl(std::index_sequence_for<T...>(), i...);
  }
};

template <class... T> Data_Index_scheme Data_Index_scheme::data_index() {
  return concatenate(myData_Index<T>::data_index()...);
}

#endif // MYDATAINDEX_H
