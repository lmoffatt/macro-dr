#ifndef MYNEWCOMMANDMANAGER_H
#define MYNEWCOMMANDMANAGER_H
#include "mygrammar.h"
#include "mycompilation.h"






namespace grammar {


template <class...> class CommandManager;

template <class...Ts, class...Tptr>
class CommandManager<Cs<Ts...>, Cs<Tptr...>>
{
public:
    typedef CommandManager<Cs<Ts...>, Cs<Tptr*...>> self_type;

    typedef Cs<Ts...> myTypes;
    typedef Cs<Tptr...> myPtrTypes;


private:
    std::tuple<std::map<std::string,Ts>...> data_;
    std::tuple<std::map<std::string,Tptr*>...> dataPtr_;

    std::tuple<std::map<std::string,Compiled_Expression<self_type,Ts>*>...> fe_;
    std::tuple<std::map<std::string,Compiled_Expression<self_type,Tptr*>*>...> feptr_;


    std::tuple<std::map<std::string,Compiled_Function_Typed<self_type,Ts>*>...> f_;
    std::tuple<std::map<std::string,Compiled_Function_Typed<self_type,Tptr*>*>...> fptr_;


    std::tuple<std::map<std::string,UnaryOperatorTyped<Ts>*>...> Uop_;

    std::tuple<std::map<std::string,BinaryOperatorTyped<Ts>*>...> Bop_;

    std::tuple<std::map<std::string,AssignmentOperatorTyped<self_type,Ts>*>...> Aop_;

public:
    void execute(const std::string& line, std::ostream& logstream)
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
        auto &m=std::get<std::map<std::string,T>>(data_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second;
        else
            return {};
    }
    template<typename T>
    std::optional<T const *> get(Cs<T*>,const std::string& id)const
    {
        auto &m=std::get<std::map<std::string,T>>(dataPtr_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second;
        else
            return {};
    }

    template<typename T>
    std::optional<T *> get(Cs<T*>,const std::string& id)
    {
        auto &m=std::get<std::map<std::string,T>>(dataPtr_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second;
        else
            return {};
    }








    template<typename T>
    void define(C<T>,std::string&& id,std::optional<Compiled_Expression<self_type,T>*> c)
    {
        if (c){
            auto&  m=  std::get<std::map<std::string,Compiled_Expression<self_type,T>*>>(fe_);
            m.emplace(id,*c);}
    }
    template<typename T>
    void define(C<T*>,std::string&& id,std::optional<Compiled_Expression<self_type,T*>*> c)
    {
        if (c){
            auto&  m=  std::get<std::map<std::string,Compiled_Expression<self_type,T*>*>>(feptr_);
            m.emplace(id,c);}
    }

    template<typename T>
    std::optional<Compiled_Expression<self_type,T>*> get_defined(C<T>,const std::string& id)
    {

        auto&  m=  std::get<std::map<std::string,Compiled_Expression<self_type,T>*>>(fe_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second;
        else return {};
    }


    template<typename T>
    std::optional<Compiled_Expression<self_type,T*>*> get_defined(C<T*>,const std::string& id)
    {

        auto&  m=  std::get<std::map<std::string,Compiled_Expression<self_type,T>*>>(fe_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second;
        else return {};
    }


    template<typename T>
    std::optional<Compiled_Function_Typed<self_type,T>*> get_Function_Typed(C<T>,const std::string& id)
    {
        auto&  m=  std::get<std::map<std::string,Compiled_Function_Typed<self_type,T>*>>(f_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second;
        else
            return {};
    }


    template<typename T>
    std::optional<Compiled_Function_Typed<self_type,T*>*> get_Function_Typed(C<T*>,const std::string& id)
    {
        auto&  m=  std::get<std::map<std::string,Compiled_Function_Typed<self_type,T*>*>>(f_);
        auto it=m.find(id);
        if (it!=m.end())
            return it->second;
        else
            return {};
    }


    template<typename T>
    void set(C<T>,std::string&& id, T&& o)
    {
        auto &m=std::get<std::map<std::string,T>>(data_);
        m.emplace(id,o);
    }

    template<typename T>
    void set(C<T>,std::string&& id, T* o)
    {
        auto &m=std::get<std::map<std::string,T>>(dataPtr_);
        m.emplace(id,o);
    }


    template<typename T>
    std::optional<BinaryOperatorTyped<T> const *> get_Binary_Operator(C<T>,const std::string& id)
    {
        auto&  m=  std::get<std::map<std::string,BinaryOperatorTyped<T>*>>(Bop_);
        auto it= m.find(id);
        if (it!=m.end())
        {
            return it->second;
        }
        else {
            return {};
        }
    }

    template<typename T>
    std::optional<UnaryOperatorTyped<T> const *> get_Unary_Operator(C<T>,const std::string& id)
    {
        auto&  m=  std::get<std::map<std::string,UnaryOperatorTyped<T>*>>(Uop_);
        auto it= m.find(id);
        if (it!=m.end())
        {
            return it->second;
        }
        else {
            return {};
        }
    }
    template<typename T>
    std::optional<AssignmentOperatorTyped<self_type,T> const *> get_Assignment_Operator(C<T>,const std::string& id)
    {
        auto&  m=  std::get<std::map<std::string,AssignmentOperatorTyped<self_type,T>*>>(Aop_);
        auto it= m.find(id);
        if (it!=m.end())
        {
            return it->second;
        }
        else {
            return {};
        }
    }

};


} // namespace grammar


#endif // MYNEWCOMMANDMANAGER_H
