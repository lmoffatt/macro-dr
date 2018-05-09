#ifndef MYNEWCOMMANDMANAGER_H
#define MYNEWCOMMANDMANAGER_H
#include "mygrammar.h"
#include "mycompilation.h"



namespace grammar {


template <class...> class CommandManager;
template <class...> class DataManager;


template <class self_type,class...Ts, class...Tptr>
class DataManager<self_type,Cs<Ts...>, Cs<Tptr...>>
{
public:

typedef Cs<Ts...> myTypes;
typedef Cs<Tptr...> myPtrTypes;


friend self_type;

private:
std::tuple<std::map<std::string,Ts>...> data_;
std::tuple<std::map<std::string,std::unique_ptr<Tptr>>...> dataPtr_;

std::tuple<std::map<std::string,std::unique_ptr<Compiled_Expression<self_type,Ts>>>...> fe_;
std::tuple<std::map<std::string,std::unique_ptr<Compiled_Expression<self_type,Tptr*>>>...> feptr_;


std::tuple<std::map<std::string,std::unique_ptr<Compiled_Function_Typed<self_type,Ts>>>...> f_;
std::tuple<std::map<std::string,std::unique_ptr<Compiled_Function_Typed<self_type,Tptr*>>>...> fptr_;


std::tuple<std::map<std::string,std::unique_ptr<UnaryOperatorTyped<Ts>>>...> Uop_;

std::tuple<std::map<std::string,std::unique_ptr<BinaryOperatorTyped<Ts>>>...> Bop_;

std::tuple<std::map<std::string,std::unique_ptr<AssignmentOperatorTyped<self_type,Ts>>>...> Aop_;

};




template <class...types, class...commands>
class CommandManager<Cs<types...>, Cs<commands...>>
{

public:
    typedef CommandManager<Cs<types...>, Cs<commands...>> self_type;

    typedef typename extract_regular_types<Cs<types...>,Cs<commands...>>::type myTypes;

    typedef typename extract_pointer_types<Cs<types...>,Cs<commands...>>::type myPtrTypes;


    static typename myTypes::fdgse d;

    DataManager<self_type,myTypes,myPtrTypes> d_;


public:
    void execute(const std::string& line, std::ostream& /*logstream*/)
    {
        auto sta=grammar::string_to_Statement(line);
        if (sta)
        {
            auto c=compile(this,sta);
            (*c)->execute(this);
        }

    }


    template<typename T>
    std::optional<T> get(C<T>,const std::string& id)const
    {
        auto &m=std::get<std::map<std::string,T>>(d_.data_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second;
        else
            return {};
    }
    template<typename T>
    std::optional<T const *> get(Cs<T*>,const std::string& id)const
    {
        auto &m=std::get<std::map<std::string,T>>(d_.dataPtr_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second->get();
        else
            return {};
    }

    template<typename T>
    std::optional<T *> get(Cs<T*>,const std::string& id)
    {
        auto &m=std::get<std::map<std::string,T>>(d_.dataPtr_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second->get();
        else
            return {};
    }




    template<typename T>
    void define(C<T>,std::string&& id,std::optional<Compiled_Expression<self_type,T>*> c)
    {
        if (c){
            //typename myTypes::j  f;

            //if constexpr (!myTypes::template has<T>)  typename T::noexta k;
            //if constexpr (myTypes::template has<T>)  typename T::esta k;
           // static_assert (myTypes::template has<T>,"class not in list");
            auto&  m=  std::get<std::map<std::string,std::unique_ptr<Compiled_Expression<self_type,T>>>>(d_.fe_);
            m.emplace(id,*c);}
    }
    template<typename T>
    void define(C<T*>,std::string&& id,std::optional<Compiled_Expression<self_type,T*>*> c)
    {
        if (c){
            auto&  m=  std::get<std::map<std::string,std::unique_ptr<Compiled_Expression<self_type,T*>>>>(d_.feptr_);
            m.emplace(id,c);}
    }


    template <class T, class F,class ...Args>
    void insert_function(C<T>,std::string&& id,F&& f, std::tuple<grammar::argument<Args>...>&& args)
    {
            auto&  m=  std::get<std::map<std::string,std::unique_ptr<Compiled_Function_Typed<self_type,T>>>>(d_.f_);
            m.emplace(id,make_compiled_function(C<T>{},std::forward<F>(f), std::forward(args)));
    }
    template <class T, class ...Args>
    void insert_function(C<T>,std::string&& id,Constructor<T>, std::tuple<grammar::field<T,Args>...>&& args)
    {
            auto&  m=  std::get<std::map<std::string,std::unique_ptr<Compiled_Function_Typed<self_type,T>>>>(d_.f_);
            m.emplace(id,make_compiled_constructor(this,C<T>{},Constructor<T>{}, std::forward<std::tuple<grammar::field<T,Args>...>>(args)));
    }


    template<class type>
    void insert_constructor()
    {
        insert_function(C<type>{},type::name.c_str(),Constructor<type>{},type::get_constructor_fields());
    }

    template<class type>
    void insert_command()
    {
        typedef typename result_of_command<type>::type result_type;
        insert_function(C<result_type>{},type::name.c_str(),type::operator(),type::get_arguments());
    }






    template<typename T>
    std::optional<Compiled_Expression<self_type,T>*> get_defined(C<T>,const std::string& id)
    {

        auto&  m=  std::get<std::map<std::string,std::unique_ptr<Compiled_Expression<self_type,T>>>>(d_.fe_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second.get();
        else return {};
    }


    template<typename T>
    std::optional<Compiled_Expression<self_type,T*>*> get_defined(C<T*>,const std::string& id)
    {

        auto&  m=  std::get<std::map<std::string,std::unique_ptr<Compiled_Expression<self_type,T*>>>>(d_.feptr_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second.get();
        else return {};
    }


    template<typename T>
    std::optional<Compiled_Function_Typed<self_type,T>*> get_Function_Typed(C<T>,const std::string& id)
    {
        auto&  m=  std::get<std::map<std::string,std::unique_ptr<Compiled_Function_Typed<self_type,T>>>>(d_.f_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second.get();
        else
            return {};
    }


    template<typename T>
    std::optional<Compiled_Function_Typed<self_type,T*>*> get_Function_Typed(C<T*>,const std::string& id)
    {
        auto&  m=  std::get<std::map<std::string,std::unique_ptr<Compiled_Function_Typed<self_type,T*>>>>(d_.fptr_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second.get();
        else
            return {};
    }


    template<typename T>
    void set(C<T>,std::string&& id, T&& o)
    {
        auto &m=std::get<std::map<std::string,T>>(d_.data_);
        m.emplace(id,o);
    }

    template<typename T>
    void set(C<T>,std::string&& id, T* o)
    {
        auto &m=std::get<std::map<std::string,T>>(d_.dataPtr_);
        m.emplace(id,o);
    }


    template<typename T>
    std::optional<BinaryOperatorTyped<T> const *> get_Binary_Operator(C<T>,const std::string& id)
    {
        auto&  m=  std::get<std::map<std::string,std::unique_ptr<BinaryOperatorTyped<T>>>>(d_.Bop_);
        auto it= m.find(id);
        if (it!=m.end())
        {
            return it->second.get();
        }
        else {
            return {};
        }
    }

    template<typename T>
    std::optional<UnaryOperatorTyped<T> const *> get_Unary_Operator(C<T>,const std::string& id)
    {
        auto&  m=  std::get<std::map<std::string,std::unique_ptr<UnaryOperatorTyped<T>>>>(d_.Uop_);
        auto it= m.find(id);
        if (it!=m.end())
        {
            return it->second.get();
        }
        else {
            return {};
        }
    }
    template<typename T>
    std::optional<AssignmentOperatorTyped<self_type,T> const *> get_Assignment_Operator(C<T>,const std::string& id)
    {
        auto&  m=  std::get<std::map<std::string,std::unique_ptr<AssignmentOperatorTyped<self_type,T>>>>(d_.Aop_);
        auto it= m.find(id);
        if (it!=m.end())
        {
            return it->second.get();
        }
        else {
            return {};
        }
    }

};


} // namespace grammar


#endif // MYNEWCOMMANDMANAGER_H
