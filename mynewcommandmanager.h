#ifndef MYNEWCOMMANDMANAGER_H
#define MYNEWCOMMANDMANAGER_H
#include "mygrammar.h"
#include "mycompilation.h"
#include "mycontainer.h"


namespace grammar {


template <class...> class CommandManager;
template <class...> class DataManager;


template <class self_type,class...Ts,  class... allT,  class... Bop_terms, class ...Uop_terms>
class DataManager<self_type,Cs<Ts...>, Cs<allT...>, Cs<Bop_terms...>, Cs<Uop_terms...>>
{
public:
    
    //   typedef typename Cs<Ts...>::TS_daaMap data_map_types;
    //  typedef typename Cs<allT...>::all_fe def_map_fe;

    typedef Cs<Ts...> myTypes;
    
    
    friend self_type;
    
private:
    template<class T>
    using data_map=std::map<std::string,unique_if_ptr_t<T>>;
    
    template<class T>
    using def_map=std::map<std::string,std::unique_ptr<Compiled_Expression<self_type,T>>>;
    
    template< class T>
    using fun_map=std::map<std::pair<std::string,std::map<std::string,std::string>>,std::unique_ptr<Function_Compiler_Typed<self_type,typename C<T>::type>>>;
    
    template<class T>
    using unary_map=std::map<std::string,std::unique_ptr<UnaryOperatorTyped<T>>>;
    
    template<class T>
    using binary_map=std::map<std::string,std::unique_ptr<BinaryOperatorTyped<T>>>;
    
    
    std::tuple<data_map<Ts>...> data_;
    std::map<std::string,std::unique_ptr<Identifier_Compiler<self_type>>> dataId_;
    
    std::map<std::string,std::unique_ptr<Literal_Compiler<self_type>>> literalId_;


    //std::tuple<std::map<std::string,Ts>...> data_;
    //std::tuple<std::map<std::string,std::unique_ptr<Tptr>>...> dataPtr_;
    
    std::tuple<def_map<allT>...> fe_;
    std::map<std::string,Compiled_Statement<self_type>*> sta_;   // naked pointer, shared with fe_
    
    //std::tuple<std::map<std::string,std::unique_ptr<Compiled_Expression<self_type,allT>>>...> fe_;
    //std::tuple<std::map<std::string,std::unique_ptr<Compiled_Expression<self_type,allTptr*>>>...> feptr_;
    
    std::tuple<fun_map<typename C<allT>::type>...> f_;
    //typedef  typename std::tuple<fun_map<allT>...>::test test2;


    std::map<std::pair<std::string,std::map<std::string, std::string>>,Function_Compiler<self_type>*> fgen_;  //naked pointer shared with f_:
    
    //std::tuple<std::map<std::string,std::unique_ptr<Compiled_Function_Typed<self_type,fn_results>>>...> f_;
    //std::tuple<std::map<std::string,std::unique_ptr<Compiled_Function_Typed<self_type,fn_resultsPtr*>>>...> fptr_;
    
    
    std::tuple<unary_map<Uop_terms>...> Uop_;
    
    //std::tuple<std::map<std::string,std::unique_ptr<UnaryOperatorTyped<Uop_terms>>>...> Uop_;
    
    std::tuple<binary_map<Bop_terms>...> Bop_;
    
    
    //std::tuple<std::map<std::string,std::unique_ptr<BinaryOperatorTyped<Bop_terms>>>...> Bop_;
    
    
    //std::tuple<std::map<std::string,std::unique_ptr<AssignmentOperatorTyped<self_type,Aop_terms>>>...> Aop_;
    
};




template <class...types, class...commands, template <typename...>class... template_commands>
class CommandManager<Cs<types...>, Cs<commands...>, CCs<template_commands...>>
{
    
public:
    typedef CommandManager<Cs<types...>, Cs<commands...>, CCs<template_commands...>> self_type;
    
    typedef typename extract_types<Cs<types...>,Cs<commands...>,CCs<template_commands...> >::type myTypes;
    
    typedef class_concatenate_t<Cs<commands...>,build_all_commands_t<CCs<template_commands...>,extract_types_t<Cs<types...>,Cs<commands...>> > > myCommands;

    typedef typename extract_function_return_types<myTypes,Cs<commands...>>::type myFuncReturns;

    typedef class_set_union_t<myTypes,myFuncReturns> allTypes;


    //typedef typename myTypes::f jaja;
    //
    //typedef typename myFuncReturns::f jajaja;

    //    typedef typename allTypes::f jaj;

    //   typedef typename myCommands::f jajtest;

    typedef  Cs<> UopTypes;
    
    typedef  Cs<> BopTypes;
    
    //static typename myTypes::fdgse d;
    
    //static typename myPtrTypes::pointersfdgse g;
    
    
    typedef DataManager<self_type,myTypes,allTypes,UopTypes,BopTypes> Dm;
    Dm d_;
    
    
    template< typename type>
    void insert_Literal()
    {
        d_.literalId_.emplace(my_trait<type>::className.str(),new Literal_Compiler_Typed<self_type,type>());
    }


    
public:
    
    template<typename... Ts>
    void insert_Literal(Cs<Ts...>)
    {
        (insert_Literal<Ts>(),...);
    }


    CommandManager():d_{}
    {
        insert_constructor(myTypes());
        insert_valuer(myTypes());
        insert_loader(myTypes());
        insert_Literal(myTypes());
        insert_command(myCommands());
    }
    
    
    void execute(const std::string& line, std::ostream& logstream)
    {
        if (!line.empty())
        {
            auto sta=grammar::string_to_Statement(line);
            if (sta)
            {
                auto c=compile(this,sta.value());
                if (c==nullptr) logstream<<"compiling error ";
                else if (*c)
                {
                    auto res=c->value()->execute(this);

                    logstream<<res<<std::endl;
                }
                else
                    logstream<<"compiling error "<<c->error()<<std::endl;
            }
            else
                logstream<<"parsing error "<<sta.error()<<std::endl;;
        }
    }
    
    //std::map<std::string,std::unique_ptr<Identifier_Compiler<self_type>>> dataId_;
    
    
    Identifier_Compiler<self_type> const * get_Identifier(const std::string& id)const
    {
        auto it=d_.dataId_.find(id);
        if (it!=d_.dataId_.end())
            return it->second.get();
        else return nullptr;
    }
    
    Literal_Compiler<self_type> const * get_Literal(const std::string& id)const
    {
        auto it=d_.literalId_.find(id);
        if (it!=d_.literalId_.end())
            return it->second.get();
        else return nullptr;
    }



    
    template<typename T>
    myOptional_t<T> get(C<T>,std::enable_if_t<!is_variable_ref_v<T>,const std::string&> id)const
    {
        auto constexpr I=Index<T,myTypes>::value;
        auto &m=std::get<I>(d_.data_);
        auto it=m.find(id);
        if (it!=m.end())
            return myOptional_t<T>(it->second);
        else
            return {false,"identifier"+id+" not found for "+my_trait<T>::className.str()};
    }
    template<typename T>
    myOptional_t<T> get(C<T>,const std::string& id)
    {
        auto &m=std::get<Index<std::decay_t<T>,allTypes>::value>(d_.data_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second;
        else
            return {false,"identifier"+id+" not found for "+my_trait<T>::className.str()};
    }
    
    template<typename T>
    myOptional_t<T const *> get(Cs<T const *>,const std::string& id)const
    {
        auto &m=std::get<Index<std::decay_t<T>,allTypes>::value>(d_.data_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second->get();
        else
            return {false,"identifier"+id+" not found for "+my_trait<T>::className.str()};
    }
    
    template<typename T>
    myOptional_t<T *> get(Cs<T*>,const std::string& id)
    {
        auto &m=std::get<Index<std::decay_t<T>,allTypes>::value>(d_.data_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second->get();
        else
            return {false,"identifier"+id+" not found for "+my_trait<T*>::className.str()};
    }
    
    
    ///------------------
    ///
    // <T> get(C<T>,std::enable_if_t<std::is_const_v<T>||!std::is_lvalue_reference_v<T>,const std::string&> id)const

    template<typename T>
    bool has(C<T>,std::enable_if_t<!is_variable_ref_v<T>,const std::string&> id)const
    {
        auto &m=std::get<Index<std::decay_t<T>,allTypes>::value>(d_.data_);
        return  m.find(id)!=m.end();
    }
    template<typename T>
    bool has(C<T>,const std::string& id)
    {
        //     typename  decltype (d_.data_)::joya k;
        auto &m=std::get<Index<std::decay_t<T>,allTypes>::value>(d_.data_);
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
    void define(C<T>,std::string id,Compiled_Expression<self_type,T>* c)
    {
        if (c!=nullptr)
        {
            auto&  m=  std::get<typename Dm::template def_map<T>>(d_.fe_);
            m.emplace(id,c);
            d_.sta_.emplace(id,c);
        }
    }
    
    
    template <class T, class F,class ...Args>
    void insert_function(C<T>,std::pair<std::string,std::map<std::string, std::string>>&& id_args,F&& f, std::tuple<grammar::argument<Args>...>&& args)
    {
        //typename T::yy a;
        auto&  m=  std::get<typename Dm::template fun_map<T>>(d_.f_);
        m.emplace(id_args,make_compiled_function<self_type,T>(this,C<T>{},std::forward<F>(f),
                                                              std::forward<std::tuple<grammar::argument<Args>...>>(args)));
        auto p=m.at(id_args).get();
        d_.fgen_.emplace(std::move(id_args),p);
    }


    template <class T, class ...Args>
    void insert_function_constructor(C<T>,std::pair<std::string,std::map<std::string, std::string>>&& id_args,Constructor<T>, std::tuple<grammar::field<T,Args>...>&& args)
    {
        //typename decltype (d_.f_)::ojo_aca d;
        auto&  m=  std::get<typename Dm::template fun_map<T>>(d_.f_);
        m.emplace(id_args,make_compiled_constructor(this,C<T>{},Constructor<T>{}, std::forward<std::tuple<grammar::field<T,Args>...>>(args)));

        auto p=m.at(id_args).get();
        d_.fgen_.emplace(std::move(id_args),p);
    }
    
    template <class T, class Derived,class ...Args>
    void insert_function_derived_constructor(C<T*>,std::pair<std::string,std::map<std::string, std::string>>&& id_args,DerivedConstructor<T,Derived>, std::tuple<grammar::field<Derived,Args>...>&& args)
    {
        //typename decltype (d_.f_)::ojo_aca d;
        auto&  m=  std::get<typename Dm::template fun_map<T*>>(d_.f_);
        m.emplace(id_args,make_compiled_derived_constructor(this,C<T*>{},DerivedConstructor<T,Derived>{}, std::forward<std::tuple<grammar::field<Derived,Args>...>>(args)));

    }


    
    template<class type>
    void insert_constructor()
    {
        if constexpr (is_field_Object<type>::value)
        {
            insert_function_constructor(C<type>{},included_types<type,my_tag_t<type>>::getIdArgs(),Constructor<type>{},type::get_constructor_fields());
            if constexpr(Base_type<type>::value)
                    if constexpr (has_this_type<myFuncReturns,typename type::base_type*>::value)
            {
                typedef typename type::base_type base;
                insert_function_derived_constructor(C<base*>{},included_types<type,my_tag_t<type>>::getIdArgs(),DerivedConstructor<base,type>{}
                                            ,type::get_constructor_fields());
            }
        }
    }
    template<class type>
    void insert_loader()
    {
        if constexpr (is_read_Object<type>::value)
        {
            std::pair<std::string,std::map<std::string, std::string>> id_args
            {my_trait<type>::className.c_str(),{{"filename",my_trait<type>::className.str()}}};
            auto&  m=  std::get<typename Dm::template fun_map<type>>(d_.f_);
            m.emplace(id_args,make_compiled_loader(this,C<type>{}));
            d_.fgen_.emplace(id_args,m.at(id_args).get());
        }
    }
    template<class type>
    void insert_valuer()
    {
        //   typedef typename type::test test;
        const std::pair<std::string,std::map<std::string,std::string>> id_args(std::string(my_trait<type>::className.c_str()),std::map<std::string, std::string>({{"value",my_trait<type>::className.str()}}));
        auto&  m=  std::get<Index<type,myTypes>::value>(d_.f_);
        typedef typename std::decay_t<decltype (m)>::value_type value_type;
        typedef typename value_type::first_type first_type;
        typedef typename value_type::second_type second_type;
        typedef typename second_type::element_type::element_type element_type;



        first_type f(id_args);
        second_type second(std::unique_ptr<Function_Compiler_Typed<self_type,element_type>>(make_compiled_valuer(this,C<element_type>{})));
        value_type p(std::move(f),std::move(second));
        m.insert(std::move(p));
        d_.fgen_.emplace(id_args,m.at(id_args).get());
    }


    template<class... type>
    void insert_constructor(Cs<type...>)
    {
        (insert_constructor<type>(),...);
    }
    
    template<class... type>
    void insert_valuer(Cs<type...>)
    {
        (insert_valuer<type>(),...);
    }

    template<class... type>
    void insert_loader(Cs<type...>)
    {
        (insert_loader<type>(),...);
    }



    template<class type>
    void insert_command()
    {
        typedef result_of_command_t<type> result_type;
        
        insert_function (C<result_type>{},std::pair(type::className.c_str(),getIdFields(type::get_arguments())), &type::run ,type::get_arguments());
    }
    
    template<class... type>
    void insert_command(Cs<type...>)
    {
        (insert_command<type>(),...);
    }


    
    Compiled_Statement<self_type>const * get_defined(const std::string& id) const
    {
        
        auto&  m=  d_.sta_;
        auto it=m.find(id);
        if (it!=m.end())
            return it->second.get();
        else return nullptr;
    }
    
    
    
    
    template<typename T>
    Compiled_Expression<self_type,T>const * get_defined_typed(C<T>,const std::string& id) const
    {
        
        auto&  m=  std::get<typename Dm::template def_map<T>>(d_.fe_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second.get();
        else return nullptr;
    }
    
    std::set<Function_Compiler<self_type> const*> get_Function(const std::pair<std::string,std::map<std::string, std::string>>& id)const
    {
        std::set<Function_Compiler<self_type> const*> out;

        auto&  m=  d_.fgen_;
        auto it=m.find(id);
        if (it!=m.end()) return {it->second};
        else
        {
            it =m.lower_bound(id);
            std::set<Function_Compiler<self_type> const*> out;
            it =m.lower_bound(id);
            while ((it!=m.end())&& (it->first.first==id.first))
            {
                if(contains(it->first.second,id.second))
                    out.insert(it->second);
                ++it;
            }
            return out;

        }
    }
    
    
    
    template<typename T>
    std::set<Function_Compiler_Typed<self_type,T> const*> get_Function_Typed(C<T>,const std::pair<std::string,std::map<std::string, std::string>>& id)const
    {
        auto&  m=  std::get<typename Dm:: template fun_map<T>>(d_.f_);
        auto it=m.find(id);
        if (it!=m.end())
            return {it->second.get()};
        else
        {
            std::set<Function_Compiler_Typed<self_type,T> const*> out;
            it =m.lower_bound(id);
            while ((it!=m.end())&& (it->first.first==id.first))
            {
                if(contains(it->first.second,id.second))
                    out.insert(it->second.get());
                ++it;
            }
            return out;
        }
    }
    
    
    
    
    
    
    
    template<typename T>
    void set(C<T>,const std::string& id, T&& o)
    {
        auto &m=std::get<std::map<std::string,T>>(d_.data_);
        m.emplace(id,std::forward<T>(o));
        d_.dataId_.emplace(id,new Identifier_Compiler_Typed<self_type,T>);
    }
    
    template<typename T>
    void set(C<T>,const std::string& id, const T& o)
    {
        auto &m=std::get<Index<std::decay_t<T>,allTypes>::value>(d_.data_);
        m.emplace(id,o);
        d_.dataId_.emplace(id,new Identifier_Compiler_Typed<self_type,T>);
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
    UnaryOperatorTyped<T> const * get_Unary_Operator(C<T>,const std::string& id) const
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
