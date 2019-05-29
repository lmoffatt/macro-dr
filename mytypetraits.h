#ifndef MYTYPETRAITS_H
#define MYTYPETRAITS_H

#include <type_traits>
#include <variant>
#include <set>
#include <map>
#include <vector>
#include <utility>
#include <string_view>
#include <memory>

template <std::size_t N> struct my_static_string {
private:
  std::array<char, N> c_;

public:
  constexpr my_static_string(const char (&c)[N]) : c_{} {
    for (std::size_t i = 0; i < N; ++i)
      c_[i] = c[i];
  }

  constexpr char operator[](std::size_t i) const { return c_[i]; }

  template <std::size_t N0>
  constexpr my_static_string(my_static_string<N0> one,
                             my_static_string<N - N0 + 1> two)
      : c_{} {
    for (std::size_t i = 0; i < N0 - 1; ++i)
      c_[i] = one[i];
    for (std::size_t i = N0 - 1; i < N - 1; ++i)
      c_[i] = two[i + 1 - N0];

    c_[N - 1] = '\0';
  }

  constexpr const char *c_str() const { return &c_[0]; }

  std::string str() const { return c_str(); }
};

template <int N>
my_static_string(const char (&lit)[N]) // <- match this
    ->my_static_string<N>;

template <std::size_t N1, std::size_t N2>
constexpr auto operator+(const my_static_string<N1> &s1,
                         const my_static_string<N2> &s2)
    -> my_static_string<N1 + N2 - 1> {
  return my_static_string<N1 + N2 - 1>(s1, s2);
}

template <typename T> struct C { typedef T type; };

template <typename... Ts> struct Cs {};


template <typename T>
using Id_t=typename C<T>::type;




template <class> struct Derived_types {

  typedef Cs<> type;
  constexpr static bool value = false;
};

template <class T>
using Derived_types_t = typename Derived_types<std::decay_t<T>>::type;
template <class T>
constexpr bool static is_derived_type_v = Derived_types<T>::value;

template <class T> struct Derived_types<T *> {

  typedef Derived_types_t<T> type;
  constexpr static bool value = is_derived_type_v<T>;
};

template <class T, class = void> struct Base_type {

  constexpr static bool value = false;
};

template <class T> struct Base_type<T, std::void_t<typename T::base_type>> {
  typedef typename T::base_type type;
  constexpr static bool value = true;
};

template <class T, class = void> struct Enclosed_type {

  constexpr static bool value = false;
};

template <class T>
struct Enclosed_type<T, std::void_t<typename T::enclosed_type>> {
  typedef typename T::enclosed_type type;
  constexpr static bool value = true;
};

template <template <typename...> class... Ts> struct CCs {};

template <template <typename...> class Function,
          template <typename...> class useThisType>
struct templatePair {};

template <class...> struct has_this_type {};

template <template <typename...> class Cs, typename... Ts, typename T>
struct has_this_type<Cs<Ts...>, T> {
  static constexpr bool value = (std::is_same_v<T, Ts> || ... || false);
};

template <class... Ts>
constexpr inline bool has_this_type_v = has_this_type<Ts...>::value;

template <template <typename...> class C> class Co {};

template <class C> struct Constructor { typedef C myClass; };

template <class Base, class Derived> struct DerivedConstructor {
  typedef Base base_type;
  typedef Derived derived_type;
};

template <class C> struct Loader { typedef C myClass; };

template <class C> struct Valuer { typedef C myClass; };

template <class...> struct my_trait {};

template <> struct my_trait<std::string> {
  constexpr static auto className = my_static_string("string");
};

template <> struct my_trait<std::size_t> {
  constexpr static auto className = my_static_string("count");
};

template <> struct my_trait<int> {
  constexpr static auto className = my_static_string("integer");
};

template <> struct my_trait<bool> {
  constexpr static auto className = my_static_string("boolean");
};

template <> struct my_trait<char> {
  constexpr static auto className = my_static_string("char");
};

template <> struct my_trait<double> {
  constexpr static auto className = my_static_string("real");
};

template <typename T, typename K> struct my_trait<std::pair<T, K>> {
  constexpr static auto className =
      my_static_string("pair_") + my_trait<T>::className +
      my_static_string("_") + my_trait<K>::className;

};

template <typename T, typename K, typename... R>
struct my_trait<std::map<T, K, R...>> {
  constexpr static auto className =
      my_static_string("map_") + my_trait<T>::className +
      my_static_string("_") + my_trait<K>::className;
};

template <typename T> struct my_trait<std::set<T>> {
  constexpr static auto className =
      my_static_string("set_") + my_trait<T>::className;
};

template <typename T> struct my_trait<std::vector<T>> {
  constexpr static auto className =
      my_static_string("vector_") + my_trait<T>::className;
};

template <typename T> struct my_trait<std::unique_ptr<T>> {
  constexpr static auto className =
      my_static_string("upt_") + my_trait<T>::className;
};

template <class...> struct my_trait_tuple {};
template <typename T> struct my_trait_tuple<T> {
  constexpr static auto className = my_trait<T>::className;
};

template <typename T, typename K, typename... Ts>
struct my_trait_tuple<T, K, Ts...> {
  constexpr static auto className = my_trait<T>::className +
                                    my_static_string("_") +
                                    my_trait_tuple<K, Ts...>::className;
};
template <typename... T> struct my_trait<std::tuple<T...>> {
  constexpr static auto className =
      my_static_string("tuple_") + my_trait_tuple<T...>::className;
};

template <class C> struct my_trait<const C &> {
  constexpr static auto className = my_static_string("const_") +
                                    my_trait<std::decay_t<C>>::className +
                                    my_static_string("_ref");
};

template <class C> struct my_trait<C &> {
  constexpr static auto className =
      my_trait<std::decay_t<C>>::className + my_static_string("_ref");
};

template <class C> struct my_trait<C *> {
  constexpr static auto className =
      my_trait<std::decay_t<C>>::className + my_static_string("_ptr");
};

template <class C> struct my_trait<C> {
  constexpr static auto className = C::className;
};

template <> struct my_trait<Cs<>> {
  constexpr static auto className = my_static_string("");
};

template <class T, class... Ts> struct my_trait<Cs<T, Ts...>> {
  constexpr static auto className = my_static_string("_") +
                                    my_trait<T>::className +
                                    my_trait<Cs<Ts...>>::className;
};


template <template <class...> class > struct my_template_trait {};

template<>
struct my_template_trait<C>
{
    constexpr static auto className = my_static_string("");

};



template <typename T, typename = void>
struct has_value_type : public std::false_type {};

template <typename T>
struct has_value_type<T, std::void_t<typename T::value_type>> : public std::true_type {
};

template <typename T>
inline constexpr bool has_value_type_v = has_value_type<T>::value;

template <typename T, typename = void>
struct has_mapped_type : public std::false_type {};

template <typename T>
struct has_mapped_type<
    T, std::void_t<typename T::key_type, typename T::mapped_type>>
    : public std::true_type {};

template <typename T>
inline constexpr bool has_mapped_type_v = has_mapped_type<T>::value;

template <typename T, typename = void>
struct has_base_type  : public std::false_type {};

template <typename T>
struct has_base_type<T, std::void_t<typename T::base_type>> : public std::true_type {
    using std::true_type::value;
};


template <typename T>
inline constexpr bool has_base_type_v = has_base_type<T>();

template <typename...> class myOptional;

template <class> struct is_optional : public std::false_type {};

template <typename T, typename tag>
struct is_optional<myOptional<T, tag>> : public std::true_type {};

template <typename T>
struct is_optional<std::optional<T>> : public std::true_type {};
template <typename T>
struct is_optional<std::optional<T> &> : public std::true_type {};


template <class> struct is_myoptional : public std::false_type {};

template <typename T, typename tag>
struct is_myoptional<myOptional<T, tag>> : public std::true_type {};



template <class> struct is_optional_ref : public std::false_type {};

template <typename T>
struct is_optional_ref<std::optional<std::reference_wrapper<T>>>
    : public std::true_type {};

template <typename T>
static constexpr inline bool is_optional_v = is_myoptional<std::decay_t<T>>::value;

template <class T> struct is_container : public std::false_type {};

template <template <typename T, typename> class Co, typename T, typename Alloc>
struct is_container<Co<T, Alloc>> : public std::true_type {};

template <typename T, typename = void>
struct is_std_container :public  std::false_type {};

template <typename T>
struct is_std_container<
    T, std::void_t<decltype(std::declval<T &>().begin()),
                   decltype(std::declval<T &>().end()), typename T::value_type>>
    : public std::true_type {};

template <class T> struct contains_constructor : public std::false_type {};
template <class C>
struct contains_constructor<Constructor<C>> : public std::true_type {};

template <class T>
struct contains_derived_constructor : public std::false_type {};
template <class B, class C>
struct contains_derived_constructor<DerivedConstructor<B, C>>
    : public std::true_type {};

template <class T> struct contains_loader : public std::false_type {};
template <class C> struct contains_loader<Loader<C>> : public std::true_type {};

template <class T> struct contains_valuer : public std::false_type {};
template <class C> struct contains_valuer<Valuer<C>> : public std::true_type {};

template <class T> struct is_map : public std::false_type {};

template <class K, class T, class Comp, class Alloc>
struct is_map<std::map<K, T, Comp, Alloc>> : public std::true_type {};

template <class K, class T, class Comp, class Alloc>
struct is_map<std::multimap<K, T, Comp, Alloc>> : public std::true_type {};

template <class T> struct is_set : public std::false_type {};

template <class T, class Comp, class Alloc>
struct is_set<std::set<T, Comp, Alloc>> : public std::true_type {};

template <class T> struct is_pair : public std::false_type {};

template <class T, class K>
struct is_pair<std::pair<T, K>> : public std::true_type {};

template <class> struct is_variant : public std::false_type {};

template <typename... T>
struct is_variant<std::variant<T...>> : std::true_type {};

template <typename T, typename = void>
struct is_field_Object : public std::false_type {};

template <typename T>
struct is_field_Object<
    T, std::void_t<decltype(T::get_constructor_fields())>>
    : public std::true_type {};

template <typename T, typename = void>
struct is_label : public std::false_type {};

template <typename T>
struct is_label<
    T, std::void_t<decltype(std::declval<T &>().name())>>
    :public  std::true_type {};


template <typename T, typename = void>
struct is_read_Object : public std::false_type {};

template <typename T>
struct is_read_Object<T, std::void_t<decltype(std::declval<T &>().read(
                             std::declval<std::istream &>()))>>
    : public std::true_type {};

template <typename T, typename = void>
struct is_write_Object : public std::false_type {};

template <typename T>
struct is_write_Object<T, std::void_t<decltype(std::declval<T &>().read(
                              std::declval<std::istream &>()))>>
    :public  std::true_type {};

template <typename T, typename = void>
struct is_arg_Command : public std::false_type {};

template <typename T>
struct is_arg_Command<T, std::void_t<decltype(T::get_arguments())>>
    :public  std::true_type {};

template <typename> struct is_tuple :public  std::false_type {};

template <typename... T> struct is_tuple<std::tuple<T...>> : public std::true_type {};

template <typename> struct is_unique_ptr : public std::false_type {};

template <typename T>
struct is_unique_ptr<std::unique_ptr<T>> : public std::true_type {};

template <class> struct is_pointer_to_const : public std::false_type {};

template <typename T>
struct is_pointer_to_const<T const *> : public std::true_type {};

template <typename T>
constexpr inline bool is_pointer_to_const_v = is_pointer_to_const<T>::value;

template <class> struct is_const_ref : public std::false_type {};

template <typename T> struct is_const_ref<T const &> : public std::true_type {};

template <typename T>
constexpr inline bool is_const_ref_v = is_const_ref<T>::value;

template <typename, class = void>
struct is_variable_ref : public std::false_type {};

template <typename T>
struct is_variable_ref<
    T,
    std::enable_if_t<(!is_const_ref_v<T>)&&(!is_pointer_to_const_v<T>)&&(
                         std::is_lvalue_reference_v<T> || std::is_pointer_v<T>),
                     int>> : public std::true_type {};

template <typename T>
inline constexpr bool is_variable_ref_v = is_variable_ref<T>::value;

template <typename, class = void> struct has_get_global :public  std::false_type {};

template <typename T>
struct has_get_global<
    T, std::void_t<decltype(get(std::declval<std::stringstream &>(),
                                std::declval<T &>()))>> : std::true_type {};

template <typename, class = void>
struct has_get_global_optional : public std::false_type {};

template <typename T>
struct has_get_global_optional<
    T, std::void_t<decltype(
           get(std::declval<C<T>>(), std::declval<std::stringstream &>()))>>
    :public  std::true_type {};

template <typename, class = void> struct has_get_method : public std::false_type {};

template <typename T>
struct has_get_method<T, std::void_t<decltype(std::declval<T &>().get(
                             std::declval<std::stringstream &>()))>>
    :public  std::true_type {};

template <typename, class = void>
struct has_global_extractor :public  std::false_type {};

template <typename T>
struct has_global_extractor<
    T, std::void_t<decltype(std::declval<std::ostream &>()
                            << std::declval<std::decay_t<T> const &>())>>
    : public std::true_type {};

template <typename, class = void> struct has_push_back : public std::false_type {};

template <typename T>
struct has_push_back<T, std::void_t<decltype(std::declval<T>().push_back(
                            std::declval<typename T::value_type>()))>>
    : public std::true_type {};

template <typename T>
inline static constexpr bool has_push_back_v = has_push_back<T>::value;

template <typename, class = void> struct has_insert :public  std::false_type {};

template <typename T>
struct has_insert<T, std::void_t<decltype(std::declval<T>().insert(
                         std::declval<typename T::value_type>()))>>
    : public std::true_type {};

template <typename T>
inline static constexpr bool has_insert_v = has_push_back<T>::value;

static_assert(has_insert<std::map<double, double>>::value, "");

static_assert(true, "");

template <typename...> struct included_types {};

template <typename... Ts>
using included_types_t = typename included_types<Ts...>::type;

template <class... C> struct arg_types { typedef void type; };

template <class... C> using arg_types_t = typename arg_types<C...>::type;

struct elem_tag {};

struct value_tag {};

struct map_tag {};

struct tuple_tag {};

struct pair_tag {};

struct object_tag {};

struct base_tag {};

struct command_tag {};

struct push_back_tag {};

struct insert_tag {};

struct set_tag {};

struct myOptional_tag {};

template <class T> struct my_tag {
  typedef std::decay_t<std::remove_pointer_t<T>> dT;
  //    typedef typename T::_in_my_tag dff;
  typedef std::conditional_t<
      is_tuple<dT>::value, tuple_tag,
      std::conditional_t<
          is_field_Object<dT>::value, object_tag,
          std::conditional_t<
              is_derived_type_v<dT>, base_tag,
              std::conditional_t<
                  is_arg_Command<dT>::value, command_tag,
                  std::conditional_t<
                      is_pair<dT>::value, pair_tag,
                      std::conditional_t<
                          has_value_type_v<dT>,
                          std::conditional_t<has_mapped_type_v<dT>, map_tag,
                                             value_tag>,
                          elem_tag>>>>>>
      type;
  // typedef typename type::_in_my_tag dfff;
};

template <class T> struct my_tag_arg {
  typedef std::decay_t<T> dT;
  //    typedef typename T::_in_my_tag dff;
  typedef std::conditional_t<
      is_field_Object<dT>::value, std::pair<object_tag, arg_types_t<dT>>,
      std::conditional_t<
          is_tuple<dT>::value, std::pair<object_tag, arg_types_t<dT>>,
          std::conditional_t<
              is_arg_Command<dT>::value, command_tag,
              std::conditional_t<
                  is_pair<dT>::value, pair_tag,
                  std::conditional_t<
                      has_push_back<dT>::value, push_back_tag,
                      std::conditional_t<
                          has_insert<dT>::value,
                          std::conditional_t<has_mapped_type_v<dT>, map_tag,
                                             set_tag>,
                          elem_tag>>>>>>
      type;
  //  typedef typename type::_in_my_tag dfff;
};

template <class T> using my_tag_t = typename my_tag<T>::type;

template <class T> using my_tag_arg_t = typename my_tag_arg<T>::type;

template <typename T, class C, bool> struct Index_imp {};

template <typename T, template <typename...> class Cs, typename... Ts>
struct Index_imp<T, Cs<Ts...>, true> {
  // typedef typename T::same dgfd;
  static constexpr std::size_t value = 0;
};

template <typename T, typename U, template <typename...> class Cs,
          typename... Ts>
struct Index_imp<T, Cs<U, Ts...>, false> {
  // typedef typename  Cs<T,U>::different dg;
  static constexpr std::size_t value =
      1 + Index_imp<T, Cs<Ts...>, std::is_same<T, U>::value>::value;
};

template <typename T, class C> struct Index {};

template <typename T, typename U, template <typename...> class Cs,
          typename... Ts>
struct Index<T, Cs<U, Ts...>> {
  static constexpr std::size_t value =
      Index_imp<T, Cs<Ts...>, std::is_same<T, U>::value>::value;
};

template<typename T, typename CsT>
inline constexpr auto Index_v=Index<T,CsT>::value;


// lets try a constexpr version of Index<T,Cs<Ts...>>....


//template<class T, class... Ts>
//constexpr std::size_t findIndex(Cs<T>, Cs<Ts...>)
//{

//    for (std::size_t I=0; I<sizeof...(Ts); ++I)
//    {
//        constexpr std::size_t II=I;
//        if (std::is_same_v<T,std::tuple_element_t<II,std::tuple<Ts...>>)
//            return II;
//    }
//    return sizeof...(Ts);
//}
//Does not work!!!!





template <class T> struct included_types<T, elem_tag> {
  typedef Cs<T> type;
  // typedef typename type::elem_type test_type;
};

template <class T> struct included_types<T, pair_tag> {
  typedef Cs<included_types_t<typename std::decay_t<T>::first_type>,
             included_types_t<typename std::decay_t<T>::second_type>, T>
      type;
  //  typedef typename type::pair_tag test_type;
};

template <class...> struct included_types_tuple;

template <class... T> struct included_types_tuple<std::tuple<T...>> {
  typedef Cs<included_types_t<std::decay_t<T>>...> type;
};

template <class T> struct included_types<T, tuple_tag> {
  typedef Cs<typename included_types_tuple<std::decay_t<T>>::type, T> type;
  // typedef typename type::tuple_tag test_type;
};

template <typename T> struct included_types<T, value_tag> {
  typedef Cs<included_types_t<typename std::decay_t<T>::value_type>, T> type;
  // typedef typename type::value_tag test_type;
};

template <typename T> struct included_types<T, map_tag> {
  typedef typename std::decay_t<T>::key_type key_type;
  typedef typename std::decay_t<T>::mapped_type mapped_type;
  typedef std::pair<key_type, mapped_type> value_type;

  typedef Cs<included_types_t<value_type>, T> type;
  //  typedef typename type::here dgfd;
  // typedef typename type::map_tag test_type;
};

// static_assert (has_key_type_v<std::vector<double>> );

template <typename T> struct included_types<T, object_tag>;

template <typename T> struct included_types<T, command_tag>;

template <typename T> struct included_types<T, base_tag>;

template <typename T> struct arg_types<T, object_tag>;

template <typename T> struct arg_types<T, command_tag>;

template <typename T> struct arg_types<T, tuple_tag>;

template <typename T> struct included_types<T> {
  typedef included_types_t<T, my_tag_t<T>> type;
};

template <typename T> struct arg_types<T> {
  typedef arg_types_t<T, my_tag_t<T>> type;
};

template <class...> struct class_concatenate {};

template <class... T>
using class_concatenate_t = typename class_concatenate<T...>::type;

template <template <class...> class Op, class... T0,
          template <class...> class Cs, class... T>
struct class_concatenate<Op<T0...>, Cs<T...>> {
  typedef Op<T0..., T...> type;
};

template <template <class...> class Op, class... T0,
          template <class...> class Cs, class... T, class... Rs>
struct class_concatenate<Op<T0...>, Cs<T...>, Rs...> {
  typedef class_concatenate_t<class_concatenate_t<Op<T0...>, Cs<T...>>, Rs...>
      type;
};

template <class A, class B> struct class_set_union {};

template <class A, class B>
using class_set_union_t = typename class_set_union<A, B>::type;

template <class A, class B> struct class_arg_set_union {};

template <class A, class B>
using class_arg_set_union_t = typename class_arg_set_union<A, B>::type;

template <class A, class B> struct class_set_union_impl {};
template <class A, class B>
using class_set_union_impl_t = typename class_set_union_impl<A, B>::type;

template <class A, class B> struct class_arg_set_union_impl {};
template <class A, class B>
using class_arg_set_union_impl_t =
    typename class_arg_set_union_impl<A, B>::type;

template <template <class...> class Cs, class... Ts>
struct class_set_union_impl<Cs<Ts...>, Cs<>> {
  typedef Cs<Ts...> type;
};

template <template <class...> class Cs, class T0, class... T, class... Ts>
struct class_set_union_impl<Cs<Ts...>, Cs<T0, T...>> {
  typedef class_set_union_impl_t<class_set_union_impl_t<Cs<Ts...>, T0>,
                                 Cs<T...>>
      type;
};

template <bool, class...> struct class_set_do_union {};

template <template <class...> class Cs, class T, class... Ts>
struct class_set_do_union<false, Cs<Ts...>, T> {
  typedef Cs<Ts..., T> type;
};
template <template <class...> class Cs, class T, class... Ts>
struct class_set_do_union<true, Cs<Ts...>, T> {
  typedef Cs<Ts...> type;
};

template <template <class...> class Cs, class... Ts>
struct class_set_union<Cs<Ts...>, Cs<>> {
  typedef Cs<Ts...> type;
};

template <template <class...> class Cs, class T, class... Ts>
struct class_set_union<Cs<Ts...>, T> {
  typedef typename class_set_do_union<(std::is_same_v<Ts, T> || ... || false),
                                      Cs<Ts...>, T>::type type;
};

template <template <class...> class Cs, class T0, class... T, class... Ts>
struct class_set_union<Cs<Ts...>, Cs<T0, T...>> {
  typedef class_set_union_t<class_set_union_t<Cs<Ts...>, T0>, Cs<T...>> type;
};

template <template <class...> class Cs, class... T>
struct class_arg_set_union<Cs<T...>, Cs<>> {
  typedef Cs<T...> type;
};

template <template <class...> class Cs, class T, class... Ts>
struct class_arg_set_union<Cs<Ts...>, T> {
  typedef class_set_union_t<Cs<Ts...>, included_types_t<T>> type;
};

template <template <class...> class Cs, class T0, class... T, class... Ts>
struct class_arg_set_union<Cs<Ts...>, Cs<T0, T...>> {
  typedef class_arg_set_union_t<class_arg_set_union_t<Cs<Ts...>, T0>, Cs<T...>>
      type;
};

template <class...> struct remove_void {};

template <class... C> using remove_void_t = typename remove_void<C...>::type;

template <template <class...> class Cs, typename... Ts>
struct remove_void<Cs<Ts...>, Cs<>> {
  typedef Cs<Ts...> type;
};

template <template <class...> class Cs, typename... Ts, typename... T>
struct remove_void<Cs<Ts...>, Cs<void, T...>> {
  typedef remove_void_t<Cs<Ts...>, Cs<T...>> type;
};

template <template <class...> class Cs, typename... Ts, typename T0,
          typename... T>
struct remove_void<Cs<Ts...>, Cs<T0, T...>> {
  typedef remove_void_t<Cs<Ts..., T0>, Cs<T...>> type;
};

template <template <class...> class Cs, typename... T>
struct remove_void<Cs<T...>> {
  typedef remove_void_t<Cs<>, Cs<T...>> type;
};

template <class...> struct filter_pointer {};

template <class... C>
using filter_pointer_t = typename filter_pointer<C...>::type;

template <template <class...> class Cs, class... T0>
struct filter_pointer<Cs<T0...>> {
  typedef remove_void_t<
      Cs<std::conditional_t<std::is_pointer_v<T0>, T0, void>...>>
      type;
};

template <class...> struct filter_regular {};

template <class... C>
using filter_regular_t = typename filter_regular<C...>::type;

template <template <class...> class Cs, class... T0>
struct filter_regular<Cs<T0...>> {
  typedef remove_void_t<
      Cs<std::conditional_t<!std::is_pointer_v<T0>, T0, void>...>>
      type;
};

template <template <class> class, class...> struct apply_Op_t_to_all {};

template <template <class> class Op_t, typename T>
struct apply_Op_t_to_all<Op_t, T> {

  typedef Op_t<T> type;
};

template <template <class> class Op_t, typename... Ts>
struct apply_Op_t_to_all<Op_t, Cs<Ts...>> {
  // static typename Cs<Ts...>::hurra k;

  typedef Cs<typename apply_Op_t_to_all<Op_t, Ts>::type...> type;
};


template<class, template<class...>class> struct transfer{};


template<template<class...>class Co,template<class...>class D, class...T>
struct transfer<Co<T...>,D>
{
    typedef D<T...> type;
};

template<class S, template<class...>class D>
using transfer_t=typename transfer<S,D>::type;



// class_set_union_t<Cs<double>,Cs<std::vector<std::map<double,int>>,std::vector<int>>>::kl
// a;

#endif // MYTYPETRAITS_H
