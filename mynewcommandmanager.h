#ifndef MYNEWCOMMANDMANAGER_H
#define MYNEWCOMMANDMANAGER_H
#include "mygrammar.h"
#include "mycompilation.h"



namespace grammar {


template <class...> class CommandManager;
template <class...> class DataManager;


template <class self_type,class...Ts,  class... allT,  class... fn_results,class... Bop_terms, class ...Uop_terms>
class DataManager<self_type,Cs<Ts...>, Cs<allT...>,  Cs<fn_results...>, Cs<Bop_terms...>, Cs<Uop_terms...>>
{
public:

 //   typedef typename Cs<Ts...>::TS_porota g1;
   // typedef typename Cs<allT...>::all_porota g2;
   //typedef typename Cs<fn_results...>::porota g;

    typedef Cs<Ts...> myTypes;


    friend self_type;

private:
    template<class T>
    using data_map=std::map<std::string,unique_if_ptr_t<T>>;

    template<class T>
    using def_map=std::map<std::string,std::unique_ptr<Compiled_Expression<self_type,T>>>;

    template< class T>
    using fun_map=std::map<std::string,std::unique_ptr<Function_Compiler_Typed<self_type,T>>>;

    template<class T>
    using unary_map=std::map<std::string,std::unique_ptr<UnaryOperatorTyped<T>>>;

    template<class T>
    using binary_map=std::map<std::string,std::unique_ptr<BinaryOperatorTyped<T>>>;


    std::tuple<data_map<Ts>...> data_;

    //std::tuple<std::map<std::string,Ts>...> data_;
    //std::tuple<std::map<std::string,std::unique_ptr<Tptr>>...> dataPtr_;

    std::tuple<def_map<allT>...> fe_;
    //std::tuple<std::map<std::string,std::unique_ptr<Compiled_Expression<self_type,allT>>>...> fe_;
    //std::tuple<std::map<std::string,std::unique_ptr<Compiled_Expression<self_type,allTptr*>>>...> feptr_;

    std::tuple<fun_map<fn_results>...> f_;

    //std::tuple<std::map<std::string,std::unique_ptr<Compiled_Function_Typed<self_type,fn_results>>>...> f_;
    //std::tuple<std::map<std::string,std::unique_ptr<Compiled_Function_Typed<self_type,fn_resultsPtr*>>>...> fptr_;


    std::tuple<unary_map<Uop_terms>...> Uop_;

    //std::tuple<std::map<std::string,std::unique_ptr<UnaryOperatorTyped<Uop_terms>>>...> Uop_;

    std::tuple<binary_map<Bop_terms>...> Bop_;


    //std::tuple<std::map<std::string,std::unique_ptr<BinaryOperatorTyped<Bop_terms>>>...> Bop_;


    //std::tuple<std::map<std::string,std::unique_ptr<AssignmentOperatorTyped<self_type,Aop_terms>>>...> Aop_;

};




template <class...types, class...commands>
class CommandManager<Cs<types...>, Cs<commands...>>
{

public:
    typedef CommandManager<Cs<types...>, Cs<commands...>> self_type;

    typedef typename extract_types<Cs<types...>,Cs<commands...>>::type myTypes;

    typedef typename extract_function_return_types<Cs<types...>,Cs<commands...>>::type myFuncReturns;

    typedef class_set_union_t<myTypes,myFuncReturns> allTypes;


    typedef  Cs<> UopTypes;

    typedef  Cs<> BopTypes;

    //static typename myTypes::fdgse d;

    //static typename myPtrTypes::pointersfdgse g;


   typedef DataManager<self_type,myTypes,allTypes,myFuncReturns,UopTypes,BopTypes> Dm;
    Dm d_;




public:

    CommandManager():d_{}
    {
        auto c=(insert_constructor<types>(),...,0);
        auto d=(insert_command<commands>(),...,0);
    }


    void execute(const std::string& line, std::ostream& /*logstream*/)
    {
        auto sta=grammar::string_to_Statement(line);
        if (sta)
        {
            auto c=compile(this,sta.value());
            if (c)
                (*c)->execute(this);
        }

    }

    template<typename T>
    optional_ref_t<T> get(C<T>,std::enable_if_t<std::is_const_v<T>||!std::is_lvalue_reference_v<T>,const std::string&> id)const
    {
            auto &m=std::get<typename Dm::template data_map<std::decay_t<T>>>(d_.data_);
            auto it=m.find(id);
            if (it!=m.end())
                return it->second;
            else
                return {};
    }
    template<typename T>
    optional_ref_t<T> get(C<T>,const std::string& id)
    {
        auto &m=std::get<typename Dm::template data_map<std::decay_t<T>>>(d_.data_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second;
        else
            return {};
    }

    template<typename T>
    std::optional<T const *> get(Cs<T const *>,const std::string& id)const
    {
        auto &m=std::get<Dm::template data_map<T>>(d_.data_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second->get();
        else
            return {};
    }

    template<typename T>
    std::optional<T *> get(Cs<T*>,const std::string& id)
    {
        auto &m=std::get<Dm::template data_map<T>>(d_.data_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second->get();
        else
            return {};
    }


    ///------------------
    ///

    template<typename T>
    bool has(C<T>,std::enable_if_t<std::is_const_v<T>||!std::is_lvalue_reference_v<T>,const std::string&> id)const
    {
            auto &m=std::get<typename Dm::template data_map<std::decay_t<T>>>(d_.data_);
            return  m.find(id)!=m.end();
    }
    template<typename T>
    bool has(C<T>,const std::string& id)
    {
   //     typename  decltype (d_.data_)::joya k;
        auto &m=std::get<typename Dm::template data_map<std::decay_t<T>>>(d_.data_);
        return  m.find(id)!=m.end();
    }

    template<typename T>
    bool has(Cs<T const *>,const std::string& id)const
    {
        auto &m=std::get<typename Dm::template data_map<T>>(d_.data_);
        return  m.find(id)!=m.end();
    }

    template<typename T>
    bool has(Cs<T*>,const std::string& id)
    {
        auto &m=std::get<typename Dm::template data_map<T>>(d_.data_);
        return  m.find(id)!=m.end();
    }





///----------------------

    template<typename T>
    void define(C<T>,const std::string& id,Compiled_Expression<self_type,T>* c)
    {
        if (c!=nullptr)
        {
            auto&  m=  std::get<typename Dm::template def_map<T>>(d_.fe_);
            m.emplace(id,c);
        }
    }


    template <class T, class F,class ...Args>
    void insert_function(C<T>,std::string&& id,F&& f, std::tuple<grammar::argument<Args>...>&& args)
    {
        //typename T::yy a;
        auto&  m=  std::get<typename Dm::template fun_map<T>>(d_.f_);
        m.emplace(id,make_compiled_function<self_type,T>(this,C<T>{},std::forward<F>(f),
                                            std::forward<std::tuple<grammar::argument<Args>...>>(args)));
    }
    template <class T, class ...Args>
    void insert_function(C<T>,std::string&& id,Constructor<T>, std::tuple<grammar::field<T,Args>...>&& args)
    {
        //typename decltype (d_.f_)::ojo_aca d;
        auto&  m=  std::get<typename Dm::template fun_map<T>>(d_.f_);
        m.emplace(id,make_compiled_constructor(this,C<T>{},Constructor<T>{}, std::forward<std::tuple<grammar::field<T,Args>...>>(args)));
    }


    template<class type>
    void insert_constructor()
    {
        insert_function(C<type>{},my_trait<type>::name.c_str(),Constructor<type>{},type::get_constructor_fields());
    }

    template<class type>
    void insert_command()
    {
        typedef result_of_command_t<type> result_type;

        insert_function (C<result_type>{},type::name.c_str(), &type::run ,type::get_arguments());
    }






    template<typename T>
    Compiled_Expression<self_type,T>const * get_defined(C<T>,const std::string& id) const
    {

        auto&  m=  std::get<typename Dm::template def_map<T>>(d_.fe_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second.get();
        else return nullptr;
    }




    template<typename T>
    Function_Compiler_Typed<self_type,T> const* get_Function_Typed(C<T>,const std::string& id)const
    {
        auto&  m=  std::get<typename Dm:: template fun_map<T>>(d_.f_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second.get();
        else
            return nullptr;
    }




    template<typename T>
    void set(C<T>,const std::string& id, T&& o)
        {
        auto &m=std::get<std::map<std::string,T>>(d_.data_);
        m.emplace(id,std::forward<T>(o));
    }

    template<typename T>
    void set(C<T>,const std::string& id, const T& o)
        {
        auto &m=std::get<std::map<std::string,T>>(d_.data_);
        m.emplace(id,o);
    }



    template<typename T>
    BinaryOperatorTyped<T> const * get_Binary_Operator(C<T>,const std::string& id) const
    {
        auto&  m=  std::get<std::map<std::string,std::unique_ptr<BinaryOperatorTyped<T>>>>(d_.Bop_);
        auto it= m.find(id);
        if (it!=m.end())
        {
            return it->second.get();
        }
        else {
            return nullptr;
        }
    }

    template<typename T>
    UnaryOperatorTyped<T>* const get_Unary_Operator(C<T>,const std::string& id) const
    {
        auto&  m=  std::get<std::map<std::string,std::unique_ptr<UnaryOperatorTyped<T>>>>(d_.Uop_);
        auto it= m.find(id);
        if (it!=m.end())
        {
            return it->second.get();
        }
        else {
            return nullptr;
        }
    }

};


} // namespace grammar


#endif // MYNEWCOMMANDMANAGER_H
