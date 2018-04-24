#ifndef MYCOMPILATION_H
#define MYCOMPILATION_H

#include "mygrammar.h"
#include <optional>
#include <functional>
#include "myTuples.h"
#include "myoptional.h"





namespace grammar {



template<class Cm>
class Compiled_Statement
{
public:
    virtual void execute(Cm *) const=0;
    ~Compiled_Statement(){}
};



template <class Cm,class T>
class Compiled_Expression: public Compiled_Statement<Cm>
{
public:
    typedef T result_type;
    typedef std::optional<T> oT;
    virtual oT run(Cm*) const=0;
    void execute(Cm *) const{};
    ~Compiled_Expression(){}
};


template <class Cm,class T>
class Compiled_General_Assignment: public Compiled_Expression<Cm,T>
{
public:
    typedef T result_type;
    typedef std::optional<T> oT;
    oT run(Cm*) const=0;
    void execute(Cm *) const=0;
    ~Compiled_General_Assignment(){}
};





template <class Cm,class T>
class Compiled_Identifier: public Compiled_Expression<Cm,T>
{
private:
    std::unique_ptr<Identifier const> id;
public:
    typedef std::optional<T> oT;
    virtual oT run(Cm* cm) const override
    {
        auto def=cm->get_defined(C<T>{},id->value());
        if (def.has_value())
            return (**def).run(cm);
        else
            return cm->get(C<T>{},id->value());
    }

    Compiled_Identifier(Identifier const* i):id{i}{}
};

template <class Cm,class T>
std::optional<Compiled_Expression<Cm,T>*> compile(C<T>,Cm* cm,const Expression* id);

template <class Cm>
Compiled_Statement<Cm>* compile(Cm* cm,const Statement* id);


template <class Cm,class T>
class Compiled_Literal: public Compiled_Expression<Cm,T>
{
private:
    std::unique_ptr<Literal<T> const> x_;
public:
    typedef std::optional<T> oT;
    virtual oT run(Cm* ) const override
    {
        return x_->getValue();
    }

    Compiled_Literal(Literal<T> const* x):x_{x}{}

};


template <class Cm,class T>
class Compiled_Definition: public Compiled_General_Assignment<Cm,T>
{
private:
    std::unique_ptr<Definition const> x_;
public:
    typedef std::optional<T> oT;
    virtual oT run(Cm* cm) const override
    {
        std::string id=x_->id()->value();
        std::unique_ptr<Compiled_Expression<Cm,T>> c{*compile(C<T>{},cm,x_->arg(0))};
        return c->run(cm);
    }

    Compiled_Definition(Definition const* x):x_{x}{}
    void execute(Cm *cm) const
    {
        std::string id=x_->id()->value();
        cm->define(C<T>{},std::move(id),compile(C<T>{},cm,x_->arg(0)));
    }

};




template <class Cm,class T>
class Compiled_Function_Typed: public Compiled_Expression<Cm,T>
{

public:
    virtual void set_arguments(Function const *f)=0;

    ~Compiled_Function_Typed(){}
};



template <class Cm,class T, class F,class ...Args>
class Compiled_Function_Arg: public Compiled_Function_Typed<Cm,T>
{
private:
    typedef std::tuple<std::pair<std::optional<std::unique_ptr<Compiled_Expression<Cm,Args>>>,std::string>...> args_types;
    F f_;
    args_types a_;
    std::unique_ptr<Function> fn;
    template<class Arg>
    static
    std::optional<Arg>
    run_argument(Compiled_Function_Arg* self,Cm* cm,const std::pair<std::optional<std::unique_ptr<Compiled_Expression<Cm,Arg>>>, std::string>& dargs,
                 const Function::arg_map& m)
    {
        auto it=m.find(dargs.second);
        if(it!=m.end())
        {

            auto e=std::make_unique(compile<Cm,Arg>(cm,it->second.get()));
            if (auto x=dynamic_cast<Compiled_Definition<Cm,Arg>*>(e.get()); x!=nullptr)
            {
                x->execute(self);
                return x->run(cm);
            }
            else
                return e->run(cm);
        }
        else if (dargs.first.has_value())
            return (**dargs.first).run(cm);
    }
    typedef std::optional<T> oT;

    template<class Arg1, class Arg2>
    static
    void try_This_arg(const std::string& id,Compiled_Expression<Cm,Arg1>* e
                      , std::pair<std::optional<std::unique_ptr<Compiled_Expression<Cm,Arg2>>>, std::string>& argn)
    {
        if constexpr(std::is_same_v<Arg1,Arg2 >)
                if (id==argn.second)
                argn.first=std::make_unique(e);
    }

public:

    template<class Arg>
    void define(const std::string& id,Compiled_Expression<Cm,Arg>* e)
    {
        std::apply([id,e](auto& ...x)
        {(...,try_This_arg(id,e,x));},a_);
    }


    virtual oT run(Cm* cm) const override
    {
        auto& argsMap=fn->getArgs();
        auto a=std::apply([this,cm,&argsMap](auto ...x){std::make_tuple(run_argument(this,cm,x,argsMap)...);} , a_);
        return apply_optional(f_,a);
    }



    void set_arguments(Function *f){fn= std::unique_ptr<Function>(f);}

    Compiled_Function_Arg(F f, args_types a):f_{f},fn{f}{}
};







template <class Cm,class T>
class Compiled_UnaryOperation: public Compiled_Expression<Cm,T>
{
private:
    std::unique_ptr<UnaryOperation const> x_;
public:
    typedef std::optional<T> oT;
    virtual oT run(Cm* cm) const override
    {
        std::unique_ptr<Compiled_Expression<Cm,T>> c(*compile(C<T>{},cm,x_->arg(0)));
        auto Op=cm->get_Unary_Operator(C<T>{},x_->value());
        return invoke_optional_functor(Op,c->run(cm));
    }
    Compiled_UnaryOperation(UnaryOperation const* x):x_{x}{}
};

template <class Cm,class T>
class Compiled_BinaryOperation: public Compiled_Expression<Cm,T>
{
private:
    std::unique_ptr<BinaryOperation const> x_;
public:
    typedef std::optional<T> oT;
    virtual oT run(Cm* cm) const override
    {
        std::unique_ptr<Compiled_Expression<Cm,T>> c0(*compile(C<T>{},cm,x_->arg(0)));
        std::unique_ptr<Compiled_Expression<Cm,T>>  c1(*compile(C<T>{},cm,x_->arg(1)));
        auto Op=cm->get_Binary_Operator(C<T>{},x_->value());

        static_assert (is_optional<decltype (Op)>::value,"e" );
        return invoke_optional_functor(Op,c0->run(cm),c1->run(cm));
    }

    Compiled_BinaryOperation(BinaryOperation const* x):x_{x}{}
};


template <class Cm,class T>
class Compiled_Assignment: public Compiled_General_Assignment<Cm,T>
{
private:
    std::unique_ptr<Assignment const> x_;
public:
    typedef std::optional<T> oT;
    virtual oT run(Cm* cm) const override
    {

        auto Op=cm->get_Assignment_Operator(C<T>{},x_->value());
        std::unique_ptr<Compiled_Expression<Cm,T>> c{*compile(C<T>{},cm,x_->arg(0))};
        auto x=c->run(cm);
        if (Op.has_value()&&x.has_value())
            return (**Op)(cm,x_->id()->value(),x.value());
        else return {};
    }

    Compiled_Assignment(Assignment const* x):x_{x}{}
    void execute(Cm *cm) const
    {
        std::string id=x_->id()->value();
        oT o=run(cm);
        if (o)
            cm->set(C<T>{},std::move(id),std::move(o.value()));
    }

};



namespace comp
{
template<class...T> struct compiler_imp{};

template <class Cm, typename...Ts, typename...Tptr >
struct compiler_imp<Cm,Cs<Ts...>,Cs<Tptr...>>
{

private:

    static std::optional<Compiled_Statement<Cm>*> get_Unary_Operation_impl(Cm*,Cs<>,UnaryOperation const *)
    {
        return {};
    }

    template<typename T,typename... Ks>
    static std::optional<Compiled_Statement<Cm>*> get_Unary_Operation_impl(Cm* cm, Cs<T,Ks...>,UnaryOperation const * a)
    {

        if (auto x=compile(C<T>{},cm,a->arg(0)); x)
            return new Compiled_UnaryOperation<Cm,T>(a);
        else return get_Unary_Operation_impl(cm,Cs<Ks...>{},a);
    }

public:


    static std::optional<Compiled_Statement<Cm>*> get_Unary_Operation(Cm* cm,UnaryOperation const * id)
    {
        auto out= get_Unary_Operation_impl(cm,Cs<Ts...>{},id);
        if (out) return out;
        else return get_Unary_Operation_impl(cm,Cs<Tptr...>{}, id);

    }

private:

    static std::optional<Compiled_Statement<Cm>*> get_Binary_Operation_impl(Cm*,Cs<>,BinaryOperation const * )
    {
        return {};
    }

    template<typename T,typename... Ks>
    static std::optional<Compiled_Statement<Cm>*> get_Binary_Operation_impl(Cm* cm, Cs<T,Ks...>,BinaryOperation const * a)
    {

        if (auto x=compile(C<T>{},cm,a->arg(0)); x)
            return new Compiled_BinaryOperation<Cm,T>(a);
        else return get_Binary_Operation_impl(cm,Cs<Ks...>{},a);
    }

public:


    static std::optional<Compiled_Statement<Cm>*> get_Binary_Operation(Cm* cm,BinaryOperation const * id)
    {
        auto out= get_Binary_Operation_impl(cm,Cs<Ts...>{},id);
        if (out) return out;
        else return get_Binary_Operation_impl(cm,Cs<Tptr...>{}, id);

    }

private:

    static std::optional<Compiled_Statement<Cm>*> get_Definition_impl(Cm*,Cs<>,Definition const * )
    {
        return {};
    }

    template<typename T,typename... Ks>
    static std::optional<Compiled_Statement<Cm>*> get_Definition_impl(Cm* cm, Cs<T,Ks...>,Definition const * a)
    {

        if (auto x=compile<Cm,T>(C<T>{},cm,a->arg(0)); x)
            return new Compiled_Definition<Cm,T>(a);
        else return get_Definition_impl(cm,Cs<Ks...>{},a);
    }

public:


    static std::optional<Compiled_Statement<Cm>*> get_Definition(Cm* cm,Definition const * id)
    {
        auto out= get_Definition_impl(cm,Cs<Ts...>{},id);
        if (out) return out;
        else return get_Definition_impl(cm,Cs<Tptr...>{}, id);

    }

private:

    static std::optional<Compiled_Statement<Cm>*> get_Assignment_impl(Cm*,Cs<>,Assignment const* )
    {
        return {};
    }

    template<typename T,typename... Ks>
    static std::optional<Compiled_Statement<Cm>*> get_Assignment_impl(Cm* cm, Cs<T,Ks...>,Assignment const* a)
    {

        if (auto x=compile<Cm,T>(C<T>{},cm,a->arg(0)); x)
            return new Compiled_Assignment<Cm,T>(a);
        else return get_Assignment_impl(cm,Cs<Ks...>{},a);
    }

public:


    static std::optional<Compiled_Statement<Cm>*> get_Assignment(Cm* cm,Assignment const* id)
    {
        auto out= get_Assignment_impl(cm,Cs<Ts...>{},id);
        if (out) return out;
        else return get_Assignment_impl(cm,Cs<Tptr...>{}, id);

    }
private:

    static std::optional<Compiled_Statement<Cm>*> get_Literal_impl(Cm*,Cs<>,LiteralGeneric const *  )
    {
        return {};
    }

    template<typename T,typename... Ks>
    static std::optional<Compiled_Statement<Cm>*> get_Literal_impl(Cm* cm, Cs<T,Ks...>,LiteralGeneric const *  id)
    {
        if (auto x=dynamic_cast<Literal<T> const*>(id); x!=nullptr)
            return new Compiled_Literal<Cm,T>(x);
        else return get_Literal_impl(cm,Cs<Ks...>{},id);
    }

public:


    static std::optional<Compiled_Statement<Cm>*> get_Literal(Cm* cm,LiteralGeneric const *  id)
    {
        auto out= get_Literal_impl(cm,Cs<Ts...>{},id);
        if (out) return out;
        else return get_Literal_impl(cm,Cs<Tptr...>{}, id);

    }

private:

    static std::optional<Compiled_Statement<Cm>*> get_Function_impl(Cm*,Cs<>,Function const  *)
    {
        return {};
    }

    template<typename T,typename... Ks>
    static std::optional<Compiled_Statement<Cm>*> get_Function_impl(Cm* cm, Cs<T,Ks...>,Function const * f)
    {
        auto out=cm->get_Function_Typed(C<T>{},f->value());
        if (out)
            return out;
        else return get_Function_impl(cm,Cs<Ks...>{},f);
    }

public:
    static std::optional<Compiled_Statement<Cm>*> get_Function(Cm* cm,Function const* f)
    {

        auto out= get_Function_impl(cm,Cs<Ts...>{},f);
        if (out) return out;
        else return get_Function_impl(cm,Cs<Tptr...>{}, f);

    }


private:

    static std::optional<Compiled_Statement<Cm>*> get_Identifier_impl(Cm*,Cs<>,Identifier const * )
    {
        return {};
    }

    template<typename T,typename... Ks>
    static std::optional<Compiled_Statement<Cm>*> get_Identifier_impl(Cm* cm, Cs<T,Ks...>,Identifier const * id)
    {
        auto out=cm->get(C<T>{},id->value());
        if (out)
        {
            auto c=new Compiled_Identifier<Cm,T>(id);
            Compiled_Statement<Cm>* t=c;
            return t;

        }
        else return get_Identifier_impl(cm,Cs<Ks...>{},id);
    }

public:


    static std::optional<Compiled_Statement<Cm>*> get_Identifier(Cm* cm,Identifier const * id)
    {
        auto out= get_Identifier_impl(cm,Cs<Ts...>{},id);
        if (out) return out;
        else return get_Identifier_impl(cm,Cs<Tptr...>{}, id);

    }






};



template <class Cm>
using compiler=compiler_imp<Cm, typename Cm::myTypes,typename Cm::myPtrTypes>;


}// namespace comp























template <class Cm,class T>
std::optional<Compiled_Expression<Cm,T>*> compile(C<T>,Cm* cm,Term const* id)
{
    if(auto x=dynamic_cast<Identifier const *>(id); x!=nullptr)
    {
        return new Compiled_Identifier<Cm,T>(x);
    }
    else if(auto x=dynamic_cast<Literal<T> const *>(id); x!=nullptr)
    {
        return new Compiled_Literal<Cm,T>(x);
    }
    else   if(auto x=dynamic_cast<Function const *>(id); x!=nullptr)
    {
        auto f=cm->template get_Function_Typed(C<T>{},x->value());
        if (f.has_value())
        {
            (**f).set_arguments(x);
            return f;
        }
        else return {};
    }
    else return {};
}



template <class Cm,class T>
std::optional<Compiled_Expression<Cm,T>*> compile(C<T>,Cm* cm,Expression const* id)
{
    if(auto x=dynamic_cast<Identifier const *>(id); x!=nullptr)
    {
        return new Compiled_Identifier<Cm,T>(x);
    }
    else if(auto x=dynamic_cast<Literal<T> const *>(id); x!=nullptr)
    {
        return new Compiled_Literal<Cm,T>(x);
    }
    else   if(auto x=dynamic_cast<Function const *>(id); x!=nullptr)
    {
        auto f=cm->template get_Function_Typed(C<T>{},x->value());
        if (f.has_value())
        {
            (**f).set_arguments(x);
            return f;
        }
        else return {};
    }
    else if(auto x=dynamic_cast<UnaryOperation const *>(id); x!=nullptr)
    {
        return new Compiled_UnaryOperation<Cm,T>(x);
    }
    else if(auto x=dynamic_cast<BinaryOperation const *>(id); x!=nullptr)
    {
        return new Compiled_BinaryOperation<Cm,T>(x);
    }
    else if(auto x=dynamic_cast<Assignment const *>(id); x!=nullptr)
    {
        return new Compiled_Assignment<Cm,T>(x);
    }
    else if(auto x=dynamic_cast<Definition const *>(id); x!=nullptr)
    {
        return new Compiled_Definition<Cm,T>(x);
    }
    else return {};

}

template <class Cm>
std::optional<Compiled_Statement<Cm>*> compile(Cm* cm,Function const* fn)
{

    auto f=cm->template get_Function(fn);
    if (f.has_value())
    {
        f->set_arguments(fn);
        return f;
    }
    else return {};
}

template <class Cm>
std::optional<Compiled_Statement<Cm>*> compile(Cm* cm,std::optional<const Statement*> sta)
{
    if (sta)
    {
        auto id=sta.value();

        if(auto x=dynamic_cast<const Identifier*>(id); x!=nullptr)
        {
            return comp::compiler<Cm>::get_Identifier(cm,x);
        }
        else if(auto x=dynamic_cast<LiteralGeneric const *>(id); x!=nullptr)
        {
            return comp::compiler<Cm>::get_Literal(cm,x);
        }
        else if(auto x=dynamic_cast<Function const *>(id); x!=nullptr)
        {
            return comp::compiler<Cm>::get_Function(cm,x);
        }
        else if(auto x=dynamic_cast<UnaryOperation const *>(id); x!=nullptr)
        {
            return comp::compiler<Cm>::get_Unary_Operation(cm,x);
        }
        else if(auto x=dynamic_cast<BinaryOperation const *>(id); x!=nullptr)
        {
            return comp::compiler<Cm>::get_Binary_Operation(cm,x);
        }
        else if(auto x=dynamic_cast<Assignment const *>(id); x!=nullptr)
        {
            return comp::compiler<Cm>::get_Assignment(cm,x);
        }
        else if(auto x=dynamic_cast<Definition const *>(id); x!=nullptr)
        {
            return comp::compiler<Cm>::get_Definition(cm,x);
            //get_Definition(x);
        }
        else return {};
    }
    else return {};
}








} // namespace grammar


#endif // MYCOMPILATION_H
