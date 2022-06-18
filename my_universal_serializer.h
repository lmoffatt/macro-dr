#ifndef MY_UNIVERSAL_SERIALIZER_H
#define MY_UNIVERSAL_SERIALIZER_H

#include <array>
#include <functional>
#include <tuple>
#include <variant>
#include <string>
#include <set>
namespace us {

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

template <std::size_t N>
my_static_string(const char (&lit)[N]) // <- match this
    ->my_static_string<N>;

template <std::size_t N1, std::size_t N2>
constexpr auto operator+(const my_static_string<N1> &s1,
                         const my_static_string<N2> &s2)
    -> my_static_string<N1 + N2 - 1> {
    return my_static_string<N1 + N2 - 1>(s1, s2);
}

template <typename E>
constexpr auto to_underlying(E e) noexcept
{
    return static_cast<std::underlying_type_t<E>>(e);
}


template <std::size_t N1,std::size_t N2 >
constexpr auto concat(const char (&one)[N1], const char (&two)[N2])
{
    return (my_static_string(one)+my_static_string(two)).c_str();
}




class error_message
{
    std::string m_;
public:
    error_message(std::string error): m_{error}{}

     auto operator()()const {return m_;}
};

template<class T>
class Maybe_error: public std::variant<T,error_message>
{
public:
    using std::variant<T,error_message>::variant;
    constexpr explicit operator bool() const noexcept
    {
        return this->index()==0;
    }

    constexpr const T* operator->() const noexcept
    {
        return std::get_if<0>(*this);
    }

    constexpr T* operator->() noexcept
    {
        return std::get_if<0>(*this);
    }

    constexpr const T& operator*() const& noexcept
    {
        return std::get<0>(*this);
    }
    constexpr T& operator*() & noexcept
    {
        return std::get<0>(*this);
    }
    constexpr const T&& operator*() const&& noexcept
    {
        return std::get<0>(std::move(*this));
    }
    constexpr T&& operator*() && noexcept
    {
        return std::get<0>(std::move(*this));

    }
    Maybe_error(T t):std::variant<T,error_message>(t){}
    Maybe_error(error_message m):std::variant<T,error_message>(m){}
    const error_message&  error()const {return std::get<1>(*this);}
    error_message&  error() {return std::get<1>(*this);}
};




template<char c, char...cs>
struct token{
    constexpr bool operator()(char ch)const
    {
        return ((ch==cs)||...||(ch==c));
    }

    template<char...o>
    friend auto operator+(token,token<o...> ){return token<c,cs...,o...>{};}

   static  constexpr auto to_string(){
        return std::array<char,sizeof...(cs)+1> {c,cs...};
       }
   auto operator()()const{return std::string({c,cs...});}
};

namespace impl{
template<char, class...>
struct token_t;

template<char begin, char...cs>
struct token_t<begin, std::integer_sequence<char,cs...>>
{
    using type=token<char(begin+cs)...>;
};

};

template<char min, char max>
using token_r=typename impl::token_t<min,std::make_integer_sequence<char,1+max-min>>::type;





using sign_token=token<'+','-'>;
using exponent_token=token<'e','E'>;
using decimal_point_token=token<'.'>;
using digit_token=token_r<'0','9'>;
using lowercase_token=token_r<'a','z'>;
using uppercase_token=token_r<'A','Z'>;
using alphabetic_token=std::decay_t<decltype(lowercase_token{}+uppercase_token{})>;
using alphanumeric_token=std::decay_t<decltype(alphabetic_token{}+digit_token{})>;


template<class state,class token, class...more>
constexpr state token_to_state(char c, token t, state s, more...m )
{
    if (t(c))
        return s;
    else
        return token_to_state(c,m...);
}

template<class state>
constexpr state token_to_state(char ,state s)
{
    return s;
}


template<class state, class token,class... more>
struct token_map: public std::tuple<token,state, more...>
{
    using base_type=std::tuple<token,state, more...>;
    constexpr state operator()(char c)const
    {
        return std::apply([c](token t,state s, more...m){return token_to_state(c,t,s,m...);},*this);
    }

    constexpr token_map(token tm, state s,more...m): base_type{tm,s,m...}{}

    template<class new_token>
    constexpr auto operator()(new_token nt, state ns)const
    {
        return std::apply([nt,ns](token t,state s, more...m){
            return token_map<state,token,more...,new_token,state>(t,s,m...,nt,ns);},
                          *this);


    }
    constexpr auto operator()( state final_state)const
    {
        return std::apply([final_state](token t,state s, more...m){
            return token_map<state,token,more...,state>(t,s,m...,final_state);},
                          *this);


    }
};



template <class F, class Init, class Begin, class End >
auto map_reduce(F&& f, Init&& init, Begin be, End e)
{
    if (be!=e)
        return map_reduce(std::forward<F>(f),
                          std::invoke(std::forward<F>(f),std::forward<Init>(init),*be),
                          ++be,e);
    else
        return std::forward<Init>(init);
}

template <class F, class Init, class Begin, class End , class Cond>
auto map_reduce_while(F&& f, Init&& init, Begin be, End e, Cond cond)
{
    if (cond(init) && (be!=e))
        return map_reduce(std::forward<F>(f),
                          std::invoke(std::forward<F>(f),std::forward<Init>(init),*be),
                          ++be,e);
    else
        return std::pair(std::forward<Init>(init),be);
}


template <class F, class Init, class Container,class Cond>
constexpr std::pair<Init,std::size_t> map_reduce_while(F&& f, Init&& init,const Container& s, std::size_t ipos,std::size_t end_pos, Cond&& cond)
{
     auto c=s[ipos];
    if (cond(init) && (ipos<end_pos))
        return map_reduce_while(std::forward<F>(f),
                                f(std::forward<Init>(init),c),
                                s,
                          ++ipos,end_pos,std::forward<Cond>(cond));
    else
        return std::pair(std::forward<Init>(init),ipos);
}

struct real_serializer
{
    enum  class state{begin, end, error,initial_sign, integral_digit,decimal_point,decimal_digit,exponent_symbol,exponent_sign,exponent_digit};

    static constexpr const char *state_names []={"begin", "end", "error","initial_sign", "integral_digit","decimal_point","decimal_digit","exponent_symbol","exponent_sign","exponent_digit"};
    static constexpr bool valid_final_state(state s){return (s==state::end)|| (s==state::integral_digit)||(s==state::decimal_digit)|| (s==state::exponent_digit);}





    static std::string to_string(state s)
    {
        return state_names[to_underlying(s)];
    }


    constexpr static state next(state s,char c)
    {
        switch (s){
        case state::end:
            return state::end;

        case state::error:
            return state::error;

        case state::begin:
            if (sign_token{}(c))
                return state::initial_sign;
            else if (digit_token{}(c))
                return state::integral_digit;
            else
                return state::error;
        case state::initial_sign:
            if (digit_token{}(c))
                return state::integral_digit;
            else
                return state::error;
        case state::integral_digit:
            if (digit_token{}(c))
                return state::integral_digit;
            else if (decimal_point_token{}(c))
                return state::decimal_point;
            else if (exponent_token{}(c))
                return state::exponent_symbol;
            else
                return state::end;

        case state::decimal_point:
            if (digit_token{}(c))
                return state::decimal_digit;
            else
                return state::error;

        case state::exponent_symbol:
            if (sign_token{}(c))
                return state::exponent_sign;
            else if (digit_token{}(c))
                return state::exponent_digit;
            else
                return state::error;

        case state::exponent_sign:
             if (digit_token{}(c))
                return state::exponent_digit;
            else
                return state::error;

        case state::decimal_digit:
            if (digit_token{}(c))
                return state::decimal_digit;
            else if (decimal_point_token{}(c))
                return state::error;
            else if (exponent_token{}(c))
                return state::exponent_symbol;
            else
                return state::end;

        case state::exponent_digit:
            if (digit_token{}(c))
                return state::exponent_digit;
            else if (decimal_point_token{}(c))
                return state::error;
            else if (exponent_token{}(c))
                return state::exponent_digit;
            else
                return state::end;
        default:
            return state::error;



        }

    }


    template<class stringview>
    constexpr static auto check(stringview const& s, std::size_t pos, std::size_t maxPos)
    {
        return map_reduce_while([](state st,char c){return next(st,c);},
                                state::begin,
                                s,
                                pos,
                                maxPos,
                                [](state s){ return (s!=state::end)&&(s!=state::error);});
    }


     static Maybe_error<double> extract(std::string_view s, std::size_t pos)
    {
        auto [final_state,final_pos]=check(s,pos,s.size());
        if (valid_final_state(final_state))
            return std::strtod(s.substr(pos,final_pos-pos+1).data(),0);

        else
        {
            return std::string("Error: ")+to_string(final_state) +" "+ std::string(s.substr(pos,final_pos-pos))+ "\\/"+std::string(s.substr(final_pos,s.size()-final_pos));
        }
    }



};






} // namespace us

#endif // MY_UNIVERSAL_SERIALIZER_H
