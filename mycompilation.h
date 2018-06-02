#ifndef MYCOMPILATION_H
#define MYCOMPILATION_H

#include "mygrammar.h"
#include <optional>
#include <functional>
#include "myTuples.h"
#include "myoptional.h"
#include "myfields.h"




namespace grammar {

template<class Cm>
class Compiled_Statement
{
public:
    virtual std::string execute(Cm * cm) const=0;
    virtual bool valid_execute(Cm * cm) const=0;

    virtual Compiled_Statement* clone()const =0;


    virtual Compiled_Statement* compile_Assignment(std::string&& id)=0;
    ~Compiled_Statement(){}
};


template <class Cm,class T>
class Compiled_Expression: public Compiled_Statement<Cm>
{
public:
    //static_assert (std::is_const_v<Cm>,"hehehe" );
    typedef T result_type;
    typedef optional_ref_t<T> oT;
    virtual oT run(Cm const* cm) const=0;
    virtual oT run(Cm* cm) const=0;
    virtual bool valid_run(Cm* cm) const=0;
    virtual bool valid_run(const Cm* cm) const=0;
    virtual std::string execute(Cm * ) const override { return {};}
    virtual bool valid_execute(Cm * ) const {return true;};

    virtual Compiled_Expression* clone()const =0;

    virtual Compiled_Expression* compile_Assignment(std::string&& id) override;



    ~Compiled_Expression(){}
};


template <class Cm,class T>
class Compiled_General_Assignment: public Compiled_Expression<Cm,T>
{
public:
    typedef T result_type;
    typedef optional_ref_t<T> oT;
    virtual bool valid_execute(Cm * cm) const override=0;
    virtual Compiled_General_Assignment* clone()const =0;

    ~Compiled_General_Assignment(){}
};



template <class Cm,class T>
class Compiled_Identifier: public Compiled_Expression<Cm,T>
{
private:
    std::string id;
public:
    typedef optional_ref_t<T> oT;
    virtual oT run(Cm* cm) const override
    {
        if constexpr (has_this_type_v<typename Cm::allTypes,T>)
        {
            auto def=cm->get_defined_typed(C<T>{},id);
            if (def!=nullptr)
                return def->run(cm);
            else
                return cm->get(C<T>{},id);
        }
        else
        return cm->get(C<T>{},id);


    }
    virtual oT run(const Cm * cm[[maybe_unused]]) const override
    {
        if constexpr(std::is_const_v<T>||(!std::is_lvalue_reference_v<T>&&!std::is_pointer_v<T>))
        {
            auto def=cm->get_defined_typed(C<T>{},id);
            if (def!=nullptr)
                return def->run(cm);
            else
                return cm->get(C<T>{},id);
        }
        else return {};
    }

    virtual bool valid_run(Cm *cm) const override
    {
        if constexpr (has_this_type_v<typename Cm::allTypes,T>)
        {
            auto def=cm->get_defined_typed(C<T>{},id);
            if (def!=nullptr)
                return true;
            else
                return cm->has(C<T>{},id);
        }
        else
        return cm->has(C<T>{},id);

    }
    virtual bool valid_run(const Cm *cm[[maybe_unused]]) const override
    {
        if constexpr(std::is_const_v<T>||(!std::is_lvalue_reference_v<T>&&!std::is_pointer_v<T>))
        {
            auto def=cm->get_defined_typed(C<T>{},id);
            if (def!=nullptr)
                return true;
            else
                return cm->has(C<T>{},id);
        }
        else
        return false;
    }
    virtual Compiled_Identifier* clone()const
    {
        return new Compiled_Identifier(*this);
    }

    Compiled_Identifier(std::string&& i):id{i}{}

    Compiled_Identifier(const std::string& i):id{i}{}

    static Compiled_Identifier* create(const Identifier* x)
    {
        return new Compiled_Identifier(x->value());
    }

};

template <class Cm,class T>
optional_unique_t<Compiled_Expression<Cm,T>> compile(const Cm* cm,C<T>,const Expression* id);

template <class Cm>
optional_unique_t<Compiled_Statement<Cm>> compile(const Cm* cm,const Statement* id);


template <class Cm,class T>
class Compiled_Literal: public Compiled_Expression<Cm,T>
{
private:
    std::decay_t<T> x_;
public:
    typedef optional_ref_t<T> oT;
    virtual oT run(const Cm* ) const override
    {
        return oT(x_);
    }

    Compiled_Literal(T x):x_{x}{}


    virtual oT run(Cm *) const  override {
        return oT(x_);
    }

    virtual bool valid_run(Cm *) const override
    {
        return true;
    }
    virtual bool valid_run(const Cm *) const override
    {
        return true;
    }
    virtual Compiled_Literal* clone()const override
    {
        return new Compiled_Literal(*this);
    }


    static Compiled_Literal* create(Literal<T> const * l)
    {
        return new Compiled_Literal(l->getValue());
    }
    static Compiled_Literal* create(LiteralGeneric const * l)
    {
        std::decay_t<T> x;
        std::stringstream ss(l->value());
        ss>>x;
        return new Compiled_Literal(x);
    }

};

template<class...> class Compiled_Array{};


template <class Cm, class T>
using Compiled_Array_t=
Compiled_Array<Cm,T,my_tag_arg_t<T>>;

template <class Cm,class T>
class Compiled_Array<Cm,T,push_back_tag>
        : public Compiled_Expression<Cm,T>
{
public:
    typedef typename std::decay_t<T>::value_type value_type;

private:

    std::vector<std::unique_ptr<Compiled_Expression<Cm,value_type>>> c_;
public:
    typedef optional_ref_t<T> oT;
    typedef std::decay_t<T> dT;
    virtual oT run(const Cm* cm ) const override
    {
        dT out;
        for (auto& e: c_)
        {
            auto c=e->run(cm);
            if (! c) return {};
            else out.push_back(c.value());
        }
        return oT(out);
    }


    virtual oT run( Cm* cm ) const override
    {
        dT out;
        for (auto& e: c_)
        {
            auto c=e->run(cm);
            if (! c) return {};
            else out.push_back(c.value());
        }
        return oT(out);
    }
    Compiled_Array(const Compiled_Array& other):c_{clone_vector(other.c_)}{}

    Compiled_Array(std::vector<std::unique_ptr<Compiled_Expression<Cm,value_type>>> x):c_{std::move(x)}{}

    virtual bool valid_run(Cm * cm) const override
    {
        for (auto& e: c_)
        {
            if (!e->valid_run(cm))
                return false;
        }
        return true;
    }
    virtual bool valid_run(const Cm * cm) const override
    {
        for (auto& e: c_)
        {
            if (!e->valid_run(cm))
                return false;
        }
        return true;
    }
    virtual Compiled_Array* clone()const override
    {
        return new Compiled_Array(*this);
    }


    static Compiled_Array* create(const Cm* cm,ArrayOperation const * l)
    {

        std::vector<std::unique_ptr<Compiled_Expression<Cm,value_type>>> v;
        for (std::size_t i=0; i<l->nArgin(); ++i)
        {
            auto c=compile<Cm,value_type>(cm,C<value_type>(),l->arg(i));
            if (c)
                v.push_back(std::move(c.value()));
            else
                return nullptr;
        }
        return new Compiled_Array(std::move(v));
    }

};


template <class Cm,class T>
class Compiled_Array<Cm,T,set_tag>
        : public Compiled_Expression<Cm,T>
{
public:
    typedef typename std::decay_t<T>::value_type value_type;

private:

    std::vector<std::unique_ptr<Compiled_Expression<Cm,value_type>>> c_;
public:
    typedef optional_ref_t<T> oT;
    typedef std::decay_t<T> dT;
    virtual oT run(const Cm* cm ) const override
    {
        dT out;
        for (auto& e: c_)
        {
            auto c=e->run(cm);
            if (! c) return {};
            else out.insert(c.value());
        }
        return oT(out);
    }


    virtual oT run( Cm* cm ) const override
    {
        dT out;
        for (auto& e: c_)
        {
            auto c=e->run(cm);
            if (! c) return {};
            else out.insert(c.value());
        }
        return oT(out);
    }
    Compiled_Array(const Compiled_Array& other):c_{clone_vector(other.c_)}{}

    Compiled_Array(std::vector<std::unique_ptr<Compiled_Expression<Cm,value_type>>> x):c_{std::move(x)}{}

    virtual bool valid_run(Cm * cm) const override
    {
        for (auto& e: c_)
        {
            if (!e->valid_run(cm))
                return false;
        }
        return true;
    }
    virtual bool valid_run(const Cm * cm) const override
    {
        for (auto& e: c_)
        {
            if (!e->valid_run(cm))
                return false;
        }
        return true;
    }
    virtual Compiled_Array* clone()const override
    {
        return new Compiled_Array(*this);
    }


    static Compiled_Array* create(const Cm* cm,ArrayOperation const * l)
    {

        std::vector<std::unique_ptr<Compiled_Expression<Cm,value_type>>> v;
        for (std::size_t i=0; i<l->nArgin(); ++i)
        {
            auto c=compile<Cm,value_type>(cm,C<value_type>(),l->arg(i));
            if (c)
                v.push_back(std::move(c.value()));
            else
                return nullptr;
        }
        return new Compiled_Array(std::move(v));
    }

};



template <class Cm,class T>
class Compiled_Array<Cm,T,map_tag>
        : public Compiled_Expression<Cm,T>
{
public:
    typedef typename std::decay_t<T>::key_type key_type;
    typedef typename std::decay_t<T>::mapped_type mapped_type;
    typedef std::pair<key_type,mapped_type> value_type;



private:

    std::vector<std::unique_ptr<Compiled_Expression<Cm,value_type>>> c_;
public:
    typedef optional_ref_t<T> oT;
    typedef std::decay_t<T> dT;
    virtual oT run(const Cm* cm ) const override
    {
        dT out;

        for (auto& e: c_)
        {
            auto c=e->run(cm);
            if (! c) return {};
            else out.insert(c.value());
        }
        return oT(out);
    }


    virtual oT run( Cm* cm ) const override
    {
        dT out;

        for (auto& e: c_)
        {
            auto c=e->run(cm);
            if (! c)
                return {};
            else out.insert(c.value());
        }
        return oT(out);
    }

    Compiled_Array(const Compiled_Array& other):c_{clone_vector(other.c_)}{}

    Compiled_Array(std::vector<std::unique_ptr<Compiled_Expression<Cm,value_type>>> x):c_{std::move(x)}{}

    virtual bool valid_run(Cm * cm) const override
    {
        for (auto& e: c_)
        {
            if (!e->valid_run(cm))
                return false;
        }
        return true;
    }
    virtual bool valid_run(const Cm * cm) const override
    {
        for (auto& e: c_)
        {
            if (!e->valid_run(cm))
                return false;
        }
        return true;
    }
    virtual Compiled_Array* clone()const override
    {
        return new Compiled_Array(*this);
    }


    static Compiled_Array* create(const Cm* cm,ArrayOperation const * l)
    {

        std::vector<std::unique_ptr<Compiled_Expression<Cm,value_type>>> v;
        for (std::size_t i=0; i<l->nArgin(); ++i)
        {
            auto c=compile<Cm,value_type>(cm,C<value_type>(),l->arg(i));
            if (c)
                v.push_back(std::move(c.value()));
            else
                return nullptr;
        }
        return new Compiled_Array(std::move(v));
    }

};



template <class Cm,class T>
class Compiled_Array<Cm,T,pair_tag>
        : public Compiled_Expression<Cm,T>
{
public:
    typedef typename std::decay_t<T>::first_type first_type;
    typedef typename std::decay_t<T>::second_type second_type;

private:

    std::unique_ptr<Compiled_Expression<Cm,first_type>> f_;
    std::unique_ptr<Compiled_Expression<Cm,second_type>> s_;
public:
    typedef optional_ref_t<T> oT;
    typedef std::decay_t<T> dT;
    virtual oT run(const Cm* cm ) const override
    {
        auto f=f_->run(cm);
        auto s=s_->run(cm);
        if (f.has_value()&&s.has_value())
            return oT(std::pair(f.value(),s.value()));
        else return {};
    }


    virtual oT run( Cm* cm ) const override
    {
        auto f=f_->run(cm);
        auto s=s_->run(cm);
        if (f.has_value()&&s.has_value())
            return oT(std::pair(f.value(),s.value()));
        else return {};
    }
    Compiled_Array(const Compiled_Array& other):f_{other.f_->clone()}, s_{other.s_->clone()}{}

    Compiled_Array(std::unique_ptr<Compiled_Expression<Cm,first_type>>&& f,
                   std::unique_ptr<Compiled_Expression<Cm,second_type>>&& s):f_{std::move(f)}, s_{std::move(s)}{}

    virtual bool valid_run(Cm * cm) const override
    {
        return f_->valid_run(cm)&& s_->valid_run(cm);
    }
    virtual bool valid_run(const Cm * cm) const override
    {
        return f_->valid_run(cm)&& s_->valid_run(cm);
    }
    virtual Compiled_Array* clone()const override
    {
        return new Compiled_Array(*this);
    }


    static Compiled_Array* create(const Cm* cm,ArrayOperation const * l)
    {

        auto f=compile<Cm,first_type>(cm,C<first_type>(),l->arg(0));
        auto s=compile<Cm,second_type>(cm,C<second_type>(),l->arg(1));
        if (f.has_value()&&s.has_value())
            return new Compiled_Array(std::move(f.value()),std::move(s.value()));
        else
            return nullptr;
    }


};


template <class Cm,class T, typename...args>
class Compiled_Array<Cm,T,std::pair<object_tag,Cs<args...>>>
        : public Compiled_Expression<Cm,T>
{
public:

private:

    //typedef typename Cs<args...>::Compiled_array  object_args;
    std::tuple<std::unique_ptr<Compiled_Expression<Cm,args>>...> c_;
public:
    typedef std::decay_t<T> dT;
    typedef optional_ref_t<T> oT;
    virtual oT run(const Cm* cm ) const override
    {

        if (std::apply([&cm](auto& ...x)
        {return (x->valid_run(cm)&&...&&true);},c_))
        {
            auto out=std::apply([&cm](auto& ...x)
            {
                return dT(x->run(cm).value()...);
            },c_);
            return oT(out);
        }
        else return oT{};
    }
    virtual oT run( Cm* cm ) const override
    {

        if (std::apply([&cm](auto& ...x)
        {return (x->valid_run(cm)&&...&&true);},c_))
        {
            return std::apply([&cm](auto& ...x)
            {
                auto o=dT(x->run(cm).value()...);
                return oT(o);
            },c_);
        }
        else return oT{};
    }


    Compiled_Array(const Compiled_Array& other):c_{clone_tuple(other.c_)}{}

    Compiled_Array(std::tuple<std::unique_ptr<Compiled_Expression<Cm,args>>...>&& x):c_{std::move(x)}{}

    virtual bool valid_run(Cm * cm) const override
    {
        return  std::apply([&cm](auto& ...x)
        {return (x->valid_run(cm)&&...&&true);},c_);
    };
    virtual bool valid_run(const Cm * cm) const override
    {
        return  std::apply([&cm](auto& ...x)
        {return (x->valid_run(cm)&&...&&true);},c_);
    }
    virtual Compiled_Array* clone()const override
    {
        return new Compiled_Array(*this);
    }

    template<std::size_t ...I>
    static Compiled_Array* create(const Cm* cm,ArrayOperation const * l, std::index_sequence<I...>)
    {
        std::tuple<std::unique_ptr<Compiled_Expression<Cm,args>>...> tu=
                std::make_tuple(std::move(compile<Cm,args>(cm,C<args>(),l->arg(I)).value())...);
        if (std::apply([](auto& ...x){ return ((x.get()!=nullptr)&&...&&true);}, tu))
            return new Compiled_Array(std::move(tu));
        else
            return nullptr;

    }


    static Compiled_Array* create(const Cm* cm,ArrayOperation const * l)
    {
        return create(cm,l,std::index_sequence_for<args...>());
    }

};


template <class Cm,class T>
class Compiled_Array<Cm,T,elem_tag>
        : public Compiled_Expression<Cm,T>
{

private:

    std::unique_ptr<Compiled_Expression<Cm,T>> c_;
public:
    typedef optional_ref_t<T> oT;
    typedef std::decay_t<T> dT;
    virtual oT run(const Cm* cm ) const override
    {
        return c_->run(cm);
    }

    virtual oT run( Cm* cm ) const override
    {
        return c_->run(cm);
    }
    Compiled_Array(const Compiled_Array& other):c_{other.c_->clone()}{}

    Compiled_Array(std::unique_ptr<Compiled_Expression<Cm,T>>&& x):c_{std::move(x)}{}

    virtual bool valid_run(Cm * cm) const override
    {
        return c_->valid_run(cm);
    }
    virtual bool valid_run(const Cm * cm) const override
    {
        return c_->valid_run(cm);
    }
    virtual Compiled_Array* clone()const override
    {
        return new Compiled_Array(*this);
    }


    static Compiled_Array* create(const Cm* cm,ArrayOperation const * l)
    {


        auto c=compile<Cm,T>(cm,C<T>(),l->arg(0));
        if (c)
            return new Compiled_Array(std::move(c.value()));
        else
            return nullptr;
    }

};




template <class Cm,class T>
class Compiled_Assignment: public Compiled_General_Assignment<Cm,T>
{
private:
    std::string id_;
    std::unique_ptr<Compiled_Expression<Cm,T>> c_;
public:
    typedef optional_ref_t<T> oT;
    virtual oT run(Cm* cm) const override
    {
        return c_->run(cm);
    }

    virtual oT run(const Cm* cm) const override
    {
        return c_->run(cm);
    }


    Compiled_Assignment(const Compiled_Assignment& other): id_{other.id_}, c_{other.c_->clone()}{}

    virtual Compiled_Assignment* clone()const
    {
        return new Compiled_Assignment(*this);
    }


    Compiled_Assignment(std::string&& id,Compiled_Expression<Cm,T>* c):
        id_{std::move(id)},c_{c}{}

    Compiled_Assignment(std::string&& id,std::unique_ptr<Compiled_Expression<Cm,T>>&& c):
        id_{std::move(id)},c_{std::move(c)}{}

    static Compiled_Assignment* create(const Cm* cm,Assignment const* o )
    {
        auto c=compile<Cm,T>(cm,C<T>{},o->expr());
        if (!o->value().empty()&&c.has_value())
            return new Compiled_Assignment(o->value(),std::move(c.value()));
        else
            return nullptr;
    }



    std::string execute(Cm *cm) const
    {
        oT o=run(cm);
        if constexpr (false && std::is_lvalue_reference_v<T>)
        {
            if (o)
            {
                cm->set(C<std::decay_t<T>>{},id_,o.value().get());
                return {};
            }
            else return o.error();
        }
        else
        {
            if (o)
            {
                cm->set(C<std::decay_t<T>>{},id_,o.value());
                return {};
            }
            else return o.error();
        }
    }


    // Compiled_Statement interface
public:
    virtual bool valid_execute(Cm *cm) const override
    {
        return valid_run(cm);
    }

    // Compiled_Expression interface
public:
    virtual bool valid_run(Cm *cm) const override
    {
        return c_->valid_run(cm);
    }
    virtual bool valid_run(const Cm *cm) const override
    {
        return c_->valid_run(cm);
    }
};


template<class Cm, class T>
Compiled_Expression<Cm,T> *Compiled_Expression<Cm,T>::compile_Assignment(std::string&& id)
{
    return new Compiled_Assignment<Cm,T>(std::move(id), this);
}


template <class Cm,class T>
class Compiled_Definition: public Compiled_General_Assignment<Cm,T>
{
private:
    std::string id_;
    std::unique_ptr<Compiled_Expression<Cm,T>> c_;
public:
    typedef optional_ref_t<T> oT;
    virtual oT run(Cm* cm) const override{ return c_->run(cm);}
    Compiled_Definition(const Compiled_Definition& other): id_{other.id_}, c_{other.c_->clone()}{}

    virtual Compiled_Definition* clone()const override
    {
        return new Compiled_Definition(*this);
    }

    Compiled_Definition(std::string && id, std::unique_ptr<Compiled_Expression<Cm,T>>&& c):id_{std::move(id)}, c_{std::move(c)}{}
    std::string execute(Cm *cm) const
    {
        if constexpr(has_this_type_v<typename Cm::allTypes,T>)
        {
            cm->define(C<T>{},id_,c_->clone()); return {};
        }
        else
        {
            cm=cm;
            return std::string("try to define an unknown type: ")+ my_trait<T>::className.c_str();
        }

    }
    static Compiled_Definition* create(const Cm* cm,Definition const* o )
    {
        auto c=compile<Cm,T>(cm,C<T>{},o->expr());
        if (!o->value().empty()&&c.has_value())
            return new Compiled_Definition(o->value(),std::move(c.value()));
        else
            return nullptr;
    }

    virtual oT run(const Cm *cm) const override
    {
        if constexpr (is_variable_ref_v<T>)  return {};
        else  return c_->run(cm);
    }
    virtual bool valid_run(Cm *cm) const override { return c_->valid_run(cm);}
    virtual bool valid_run(const Cm * cm) const override
    {
        if constexpr (is_variable_ref_v<T>)  return false;
        else
        return c_->valid_run(cm);
    }

    virtual bool valid_execute(Cm *cm) const override
    {
        return c_->valid_run(cm);
    }
};


template <class Cm,class T>
class Compiled_UnaryOperation: public Compiled_Expression<Cm,T>
{
private:

    UnaryOperatorTyped<T> const * op_;
    std::unique_ptr<Compiled_Expression<Cm,T>> c_;
public:
    typedef optional_ref_t<T> oT;
    virtual oT run(Cm* cm) const override
    {
        return invoke_optional_functor(op_,c_->run(cm));
    }
    virtual oT run(Cm const* cm) const override
    {
        return invoke_optional_functor(op_,c_->run(cm));
    }

    Compiled_UnaryOperation(const Compiled_UnaryOperation& other): op_{other.op_}, c_{other.c_->clone()}{}

    virtual Compiled_UnaryOperation* clone()const override
    {
        return Compiled_UnaryOperation(*this);
    }

    Compiled_UnaryOperation(UnaryOperatorTyped<T> const * op,
                            std::unique_ptr<Compiled_Expression<Cm,T>>&& c):op_{op},c_{std::move(c)}{}

    virtual bool valid_run(Cm *cm) const override{
        return c_->valid_run(cm);
    }
    virtual bool valid_run(const Cm *cm) const override{
        return c_->valid_run(cm);
    }

    static Compiled_UnaryOperation* create(const Cm* cm,UnaryOperation const* o )
    {
        UnaryOperatorTyped<T> const * op=cm->get_Unary_Operator(C<T>{},o->value());
        if (op!=nullptr)
        {
            auto c=compile<Cm,T>(cm,C<T>{},o->arg(0));
            if (c)
                return new Compiled_UnaryOperation(op,std::move(c.value()));
            else return nullptr;
        }
        else return nullptr;
    }
};

template <class Cm,class T>
class Compiled_BinaryOperation: public Compiled_Expression<Cm,T>
{
private:
    BinaryOperatorTyped<T> const * op_;
    std::unique_ptr<Compiled_Expression<Cm,T>> c0_;
    std::unique_ptr<Compiled_Expression<Cm,T>>  c1_;
public:

    Compiled_BinaryOperation(const Compiled_BinaryOperation& other): op_{other.op_}, c0_{other.c0_->clone()},c1_{other.c1_->clone()}{}

    virtual Compiled_BinaryOperation* clone()const override
    {
        return Compiled_BinaryOperation(*this);
    }

    typedef optional_ref_t<T> oT;
    virtual oT run(Cm* cm) const override
    {
        return invoke_optional_functor(op_,c0_->run(cm),c1_->run(cm));
    }

    virtual oT run(const Cm* cm) const override
    {
        return invoke_optional_functor(op_,c0_->run(cm),c1_->run(cm));
    }

    virtual bool valid_run(Cm *cm) const override{
        return c0_->valid_run(cm)&&c1_->valid_run(cm);
    }
    virtual bool valid_run(const Cm *cm) const override{
        return c0_->valid_run(cm)&&c1_->valid_run(cm);
    }



    Compiled_BinaryOperation(BinaryOperatorTyped<T> const * op,
                             std::unique_ptr<Compiled_Expression<Cm,T>>&& c0,
                             std::unique_ptr<Compiled_Expression<Cm,T>>&&  c1):op_{op},c0_{std::move(c0)}, c1_{std::move(c1)}{}

    static Compiled_BinaryOperation* create(const Cm* cm,BinaryOperation const* o )
    {
        BinaryOperatorTyped<T> const * op=cm->get_Binary_Operator(C<T>{},o->value());
        if (op!=nullptr)
        {
            auto c=compile<Cm,T>(cm,C<T>{},o->arg(0));
            auto d=compile<Cm,T>(cm,C<T>{},o->arg(1));
            if ((c)&&(d))
                return new Compiled_BinaryOperation(op,std::move(c.value()),std::move(d.value()));
            else return nullptr;
        }
        else
            return nullptr;
    }


};

template <class Cm>
class Identifier_Compiler
{

public:

    virtual Compiled_Statement<Cm>* compile_this(const Cm * cm,const Identifier *f) const=0;
    virtual ~Identifier_Compiler(){}
};

template <class Cm,class T>
class Identifier_Compiler_Typed: public Identifier_Compiler<Cm>
{
public:
    virtual Compiled_Identifier<Cm,T>* compile_this(const Cm * ,const Identifier *f) const override
    {
        return new Compiled_Identifier<Cm,T>(f->value());
    }
    virtual ~Identifier_Compiler_Typed(){}
};


template <class Cm>
class Literal_Compiler
{

public:

    virtual Compiled_Statement<Cm>* compile_this(const Cm * cm,const LiteralGeneric *f) const=0;
    virtual ~Literal_Compiler(){}
};

template <class Cm,class T>
class Literal_Compiler_Typed: public Literal_Compiler<Cm>
{
public:
    virtual Compiled_Literal<Cm,T>* compile_this(const Cm * ,const LiteralGeneric *f) const override
    {
        auto x=dynamic_cast<const Literal<T>*>(f);
        if (x!=nullptr)
            return new Compiled_Literal<Cm,T>(x->getValue());
        else return nullptr;
    }
    virtual ~Literal_Compiler_Typed(){}
};





template <class Cm,class T>
class Compiled_Function_Typed: public Compiled_Expression<Cm,T>
{

public:

    virtual Compiled_Function_Typed* clone()const override=0;

    ~Compiled_Function_Typed(){}
};

template <class Cm>
class Function_Compiler
{

public:

    virtual Compiled_Statement<Cm>* compile_this(const Cm * cm,const Function *f) const=0;
    virtual ~Function_Compiler(){}
};

template <class Cm,class T>
class Function_Compiler_Typed: public Function_Compiler<Cm>
{

public:

    virtual Compiled_Function_Typed<Cm,T>* compile_this(const Cm * cm,const Function *f) const override=0;
    virtual ~Function_Compiler_Typed(){}
};




template <class Cm,class T, class F,class ...Args>
class Compiled_Function_Arg: public Compiled_Function_Typed<Cm,T>
{
public:
    typedef std::tuple<std::unique_ptr<Compiled_Expression<Cm,Args>>...> args_types;

private:

    F f_;
    args_types a_;


    static args_types clone(const args_types& other)
    {
        return std::apply([](auto&...x){return args_types(std::unique_ptr<Compiled_Expression<Cm,Args>>(x->clone())...);},other);
    }


    static constexpr bool has_variable_ref=(is_variable_ref<Args>::value&&...&&false);


public:
    virtual optional_ref_t<T> run(Cm* cm) const override
    {
        F f(f_);
        return apply_optional<F,Args...>(std::move(f),std::apply([&cm](auto&...x){
            return std::make_tuple(
                        x->run(cm)
                        ...
                        );},a_));
    }


    virtual optional_ref_t<T> run(const Cm* cm) const override
    {
        if constexpr(!is_variable_ref<T>::value&& !has_variable_ref)
        {
            auto c=std::apply([&cm](auto&...x){return std::make_tuple(x->run(cm)...);},a_);
            F f(f_);
            return apply_optional<F,Args...>(std::move(f),std::move(c));
        }
        else
        return {};
    }

    virtual bool valid_run(Cm* cm) const
    {
        return std::apply([&cm](auto&...x){return (x->valid_run(cm)&&...&&true);},a_);
    }

    virtual bool valid_run(const Cm* cm) const
    {
        if constexpr(!is_variable_ref<T>::value&& !has_variable_ref)
        {
            return std::apply([&cm](auto&...x){return (x->valid_run(cm)&&...&&true);},a_);
        } else
        return false;
    }

    Compiled_Function_Arg(const F& f,args_types&& a):f_{f},a_{std::move(a)}{}


    Compiled_Function_Arg(const Compiled_Function_Arg& other): f_{other.f_}, a_{clone(other.a_)}{}

    virtual Compiled_Function_Arg* clone()const override
    {
        return new Compiled_Function_Arg(*this);
    }


};




template <class Cm,class T, class F,class ...Args>
class Function_Arg: public Function_Compiler_Typed<Cm,T>
{
private:
    typedef std::tuple<std::pair<std::string, optional_unique_t<Compiled_Expression<Cm,Args>>>...> args_types;
    F f_;
    args_types a_;
private:
    template<class Arg>
    static
    optional_unique_t<Compiled_Expression<Cm,Arg>>
    compile_argument(const Cm* cm,
                     const std::pair< std::string,optional_unique_t<Compiled_Expression<Cm,Arg>>> & dargs,
                     const Function::arg_map& m)
    {
        auto it=m.find(dargs.first);
        if(it!=m.end())
        {
            return compile<Cm,Arg>(cm,C<Arg>{},it->second.get());
        }
        else if (dargs.second.has_value())
            return optional_unique_t<Compiled_Expression<Cm,Arg>>(dargs.second.value()->clone());
        else return {};
    }

public :


    virtual Compiled_Function_Typed<Cm,T>*
    compile_this(const Cm* cm,const Function *f) const override
    {
        auto& argsMap=f->getArgs();
        auto opt_args=std::apply([&cm,&argsMap](auto& ...x)
        {return std::make_tuple(compile_argument(cm,x,argsMap)...);},a_);

        if (std::apply([](auto&...x){return (x.has_value()&&...&&true);},opt_args))
        {
            auto args=std::apply([](auto&...x){return std::make_tuple(std::move(x.value())...);},opt_args);
            return  new Compiled_Function_Arg<Cm,T,F,Args...>(f_,std::move(args));
        }
        else
            return nullptr;
    }

    Function_Arg(C<T>,F f, args_types&& a):f_{f},a_{std::move(a)}{}
};




template <class Cm,class Arg>
std::pair<std::string, optional_unique_t<Compiled_Expression<Cm,Arg>>>
make_compiled_argument(const Cm* ,grammar::argument<Arg> a)
{
    if (a.default_value.has_value())
    {
        if constexpr (!std::is_lvalue_reference_v<Arg>)
                return std::pair(a.idField,optional_unique_t<Compiled_Expression<Cm,Arg>>
                                 (std::unique_ptr<Compiled_Expression<Cm,Arg>>(new Compiled_Literal<Cm,Arg>(a.default_value.value()))));
    }
    return {a.idField,optional_ref_t<std::unique_ptr<Compiled_Expression<Cm,Arg>>>()};

}

template <class Cm,class ...Args>
auto compile_all_arguments(Cm const* cm,std::tuple<grammar::argument<Args>...>&& args)
{
    return std::apply([&cm](auto& ...x){return std::make_tuple(make_compiled_argument(cm,x)...);},args);

}



template <class Cm,class T, class F,class ...Args>
auto make_compiled_function(Cm const * cm,C<T>,F&& f, std::tuple<grammar::argument<Args>...>&& args)
{
    return new Function_Arg<Cm,T,F,Args...>(C<T>{},f,compile_all_arguments<Cm,Args...>(cm,std::forward<std::tuple<grammar::argument<Args>...>>(std::move(args))));
}





template <class Cm,class C,class Arg>
auto
make_compiled_field(const Cm* /*cm*/,const grammar::field<C,Arg>& a)
{
    typedef  typename grammar::field<C,Arg>::result_type T;
    typedef  optional_unique_t<Compiled_Expression<Cm,T>> ot;
    if (a.default_value.has_value())
    {
        return std::pair(a.idField,ot(std::unique_ptr<Compiled_Expression<Cm,T>>(new Compiled_Literal<Cm,T>(a.default_value.value()))));

    }
    else
    {
        ot o;
        return std::pair<std::string,ot>(a.idField,std::move(o));
    }
}

template <class Cm,class C, class ...Args>
auto compile_all_fields(const Cm* cm,const std::tuple<grammar::field<C,Args>...>& args)
{
    return std::apply([&cm](auto& ...x){return std::make_tuple(make_compiled_field(cm,x)...);},args);

}



template <class Cm,class T, class ...Args>
auto make_compiled_constructor(const Cm *cm,C<T>,Constructor<T>,const std::tuple<grammar::field<T,Args>...>& args)
{

    auto a=compile_all_fields(cm,args);
    return new Function_Arg<Cm,T,Constructor<T>,typename grammar::field<T,Args>::result_type...>(C<T>{},Constructor<T>{},std::move(a));
}

template <class Cm,class T>
auto make_compiled_loader(const Cm *,C<T>)
{
    auto args=std::make_tuple(std::pair(std::string("filename"),optional_unique_t<Compiled_Expression<Cm,std::string>>()));
    return new Function_Arg<Cm,T,Loader<T>,std::string>(C<T>{},Loader<T>{},std::move(args));
}

template <class Cm,class T>
auto make_compiled_valuer(const Cm *,C<T>)
{
    auto args=std::make_tuple(std::pair(std::string("value"),optional_unique_t<Compiled_Expression<Cm,T>>()));
    return new Function_Arg<Cm,T,Valuer<T>,T>(C<T>{},Valuer<T>{},std::move(args));
}











namespace comp
{
template<class...T> struct compiler_imp{};

template <class Cm,class...Ts,  class... allT,  class... fn_results,class... Bop_terms, class ...Uop_terms>
struct compiler_imp<Cm,Cs<Ts...>, Cs<allT...>,  Cs<fn_results...>, Cs<Bop_terms...>, Cs<Uop_terms...>>
{

private:

    static optional_unique_t<Compiled_Statement<Cm>> get_Unary_Operation_impl(const Cm*,Cs<>,UnaryOperation const *)
    {
        return {};
    }

    template<typename T,typename... Ks>
    static optional_unique_t<Compiled_Statement<Cm>> get_Unary_Operation_impl(const Cm* cm, Cs<T,Ks...>,UnaryOperation const * a)
    {

        if (auto x=compile(cm,C<T>{},a->arg(0)); x)
            return new Compiled_UnaryOperation<Cm,T>(a->op()->clone(),x);
        else return get_Unary_Operation_impl(cm,Cs<Ks...>{},a);
    }

public:


    static optional_unique_t<Compiled_Statement<Cm>> get_Unary_Operation(const Cm* cm,UnaryOperation const * id)
    {
        auto out= get_Unary_Operation_impl(cm,Cs<Uop_terms...>{},id);
        if (out) return out;
        else return {};

    }

private:

    static optional_unique_t<Compiled_Statement<Cm>> get_Binary_Operation_impl(const Cm*,Cs<>,BinaryOperation const * )
    {
        return {};
    }

    template<typename T,typename... Ks>
    static optional_unique_t<Compiled_Statement<Cm>> get_Binary_Operation_impl(const Cm* cm, Cs<T,Ks...>,BinaryOperation const * a)
    {

        if (auto x=compile(cm,C<T>{},a->arg(0)); x)
        {
            if (auto y=compile(cm,C<T>{},a->arg(1)); y)
                return new Compiled_BinaryOperation<Cm,T>(a->op()->clone(),std::move(x),std::move(y));
        }
        else return get_Binary_Operation_impl(cm,Cs<Ks...>{},a);
    }

public:


    static optional_unique_t<Compiled_Statement<Cm>> get_Binary_Operation(const Cm* cm,BinaryOperation const * id)
    {
        return get_Binary_Operation_impl(cm,Cs<Bop_terms...>{},id);

    }

private:

    static optional_unique_t<Compiled_Statement<Cm>> get_Definition_impl(const Cm*,Cs<>,Definition const * )
    {
        return {};
    }

    template<typename T,typename... Ks>
    static optional_unique_t<Compiled_Statement<Cm>> get_Definition_impl(const Cm* cm, Cs<T,Ks...>,Definition const * a)
    {

        if (auto x=compile<Cm,T>(cm,C<T>{},a->arg(0)); x)
            return optional_unique_t<Compiled_Statement<Cm>>(new Compiled_Definition<Cm,T>(a->id()->value(),std::move(x.value())));
        else return get_Definition_impl(cm,Cs<Ks...>{},a);
    }

public:


    static optional_unique_t<Compiled_Statement<Cm>> get_Definition(const Cm* cm,Definition const * id)
    {
        return get_Definition_impl(cm,Cs<allT...>{},id);
    }

private:

    static optional_unique_t<Compiled_Statement<Cm>> get_Assignment_impl(const Cm*,Cs<>,Assignment const* )
    {
        return {};
    }

    template<typename T,typename... Ks>
    static optional_unique_t<Compiled_Statement<Cm>> get_Assignment_impl(const Cm* cm, Cs<T,Ks...>,Assignment const* a)
    {

        if (auto x=compile<Cm,T>(cm,C<T>{},a->arg(0)); x)
            return optional_unique_t<Compiled_Statement<Cm>>(new Compiled_Assignment<Cm,T>(a->id()->value(),std::move(x.value())));
        else return get_Assignment_impl(cm,Cs<Ks...>{},a);
    }

public:


    static optional_unique_t<Compiled_Statement<Cm>> get_Assignment(const Cm* cm,Assignment const* id)
    {
        return get_Assignment_impl(cm,Cs<allT...>{},id);
    }
private:

    static optional_unique_t<Compiled_Statement<Cm>> get_Literal_impl(const Cm*,Cs<>,LiteralGeneric const *  )
    {
        return {};
    }

    template<typename T,typename... Ks>
    static optional_unique_t<Compiled_Statement<Cm>> get_Literal_impl(const Cm* cm, Cs<T,Ks...>,LiteralGeneric const *  id)
    {
        if (auto x=dynamic_cast<Literal<T> const*>(id); x!=nullptr)
            return optional_unique_t<Compiled_Statement<Cm>>(new Compiled_Literal<Cm,T>(x->getValue()));
        else return get_Literal_impl(cm,Cs<Ks...>{},id);
    }

public:


    static optional_unique_t<Compiled_Statement<Cm>> get_Literal(const Cm* cm,LiteralGeneric const *  id)
    {
        return get_Literal_impl(cm,Cs<Ts...>{},id);
    }

private:

    static optional_unique_t<Compiled_Statement<Cm>> get_Function_impl(const Cm*,Cs<>,Function const  *)
    {
        return {};
    }

    template<typename T,typename... Ks>
    static optional_unique_t<Compiled_Statement<Cm>> get_Function_impl(const Cm* cm, Cs<T,Ks...>,Function const * f)
    {
        auto fun=cm->get_Function_Typed(C<T>{},f->getIdArg());
        if(!fun.empty())
        {
            for  (auto &e:fun)
            {
                auto c= e->compile_this(cm,f);
                if (c!=nullptr)
                {
                    return optional_unique_t<Compiled_Statement<Cm>>(c);

                }
            }
        }
        return get_Function_impl(cm,Cs<Ks...>{},f);
    }

public:
    static optional_unique_t<Compiled_Statement<Cm>> get_Function(const Cm* cm,Function const* f)
    {
        return get_Function_impl(cm,Cs<allT...>{},f);
    }


private:

    static optional_unique_t<Compiled_Statement<Cm>> get_Identifier_impl(const Cm*,Cs<>,Identifier const * )
    {
        return {};
    }

    template<typename T,typename... Ks>
    static optional_unique_t<Compiled_Statement<Cm>> get_Identifier_impl(const Cm* cm, Cs<T,Ks...>,Identifier const * id)
    {
        if (cm->has(C<T>{},id->value()))
        {
            auto c=new Compiled_Identifier<Cm,T>(id->value());
            return optional_unique_t<Compiled_Statement<Cm>>(std::move(c));
        }
        else
            return get_Identifier_impl(cm,Cs<Ks...>{},id);
    }

public:


    static optional_unique_t<Compiled_Statement<Cm>> get_Identifier(const Cm* cm,Identifier const * id)
    {
        return get_Identifier_impl(cm,Cs<allT...>{},id);

    }






};



template <class Cm>
using compiler=compiler_imp<Cm, typename Cm::myTypes, typename Cm::allTypes,typename Cm::myFuncReturns, typename Cm::UopTypes, typename Cm::BopTypes>;


}// namespace comp





template <class Cm,class T>
optional_unique_t<Compiled_Expression<Cm,T>> compile(const Cm* cm,C<T>,Term const* id)
{
    typedef   optional_unique_t<Compiled_Expression<Cm,T>> res_type;
    if(auto x=dynamic_cast<Identifier const *>(id); x!=nullptr)
    {
        return res_type(Compiled_Identifier<Cm,T>::create(x));
    }
    else if(auto x=dynamic_cast<Literal<T> const *>(id); x!=nullptr)
    {
        return res_type(Compiled_Literal<Cm,T>::create(x));
    }
    else if(auto x=dynamic_cast<Function const *>(id); x!=nullptr)
    {
        if constexpr (has_this_type_v<typename Cm::myFuncReturns,T>)
        {
            auto f=cm->template get_Function_Typed(C<T>{},x->value());
            if (f!=nullptr)
            {
                return f->compile_this(cm,x);
            }
            else return {};
        }
        else return {};
    }
    else return {};
}



template <class Cm,class T>
optional_unique_t<Compiled_Expression<Cm,T>> compile(const Cm* cm,C<T>,Expression const* id)
{
    typedef   optional_unique_t<Compiled_Expression<Cm,T>> res_type;
    if(auto x=dynamic_cast<Identifier const *>(id); x!=nullptr)
    {
        return res_type(Compiled_Identifier<Cm,T>::create(x));
    }
    else if(auto x=dynamic_cast<Literal<T> const *>(id); x!=nullptr)
    {
        return res_type(Compiled_Literal<Cm,T>::create(x));
    }
    else if(auto x=dynamic_cast<LiteralGeneric const *>(id); x!=nullptr)
    {
        return res_type(Compiled_Literal<Cm,T>::create(x));
    }
    else if constexpr (std::is_reference_v<T>)
    {
        if (auto x=dynamic_cast<LiteralGeneric const *>(id); x!=nullptr)
            return res_type(Compiled_Literal<Cm,T>::create(x));
    }
    if(auto x=dynamic_cast<Function const *>(id); x!=nullptr)
    {
        if constexpr (has_this_type_v<typename Cm::allTypes,T>)
        {
            auto f=cm->template get_Function_Typed(C<T>{},x->getIdArg());
            for (auto &e:f)
            {
                auto c=e->compile_this(cm,x);
                if (c!=nullptr)
                    return optional_unique_t<Compiled_Expression<Cm,T>>(c);
            }

            return {};
        }
        else return {};
    }
    else if(auto x=dynamic_cast<UnaryOperation const *>(id); x!=nullptr)
    {
        if constexpr (has_this_type_v<typename Cm::UopTypes,T>)
        {
            return res_type(Compiled_UnaryOperation<Cm,T>::create(cm,x));
        }
        else return {};
    }
    else if(auto x=dynamic_cast<BinaryOperation const *>(id); x!=nullptr)
    {
        if constexpr (has_this_type_v<typename Cm::BopTypes,T>)
        {
            return res_type(Compiled_BinaryOperation<Cm,T>::create(cm,x));
        }
        else return {};
    }
    else if(auto x=dynamic_cast<Assignment const *>(id); x!=nullptr)
    {
        return res_type(Compiled_Assignment<Cm,T>::create(cm,x));
    }
    else if(auto x=dynamic_cast<Definition const *>(id); x!=nullptr)
    {
        return res_type(Compiled_Definition<Cm,T>::create(cm,x));
    }
    else if(auto x=dynamic_cast<ArrayOperation const *>(id); x!=nullptr)
    {
        //  typedef typename res_type::res_type aqui;
        // typedef typename T::T aqui2;

        return res_type(Compiled_Array_t<Cm,T>::create(cm,x));
    }
    else return {};

}

template <class Cm>
optional_unique_t<Compiled_Statement<Cm>> compile(const Cm* cm,Function const* fn)
{

    auto f=cm->template get_Function(fn->getIdArg());
    if (!f.empty()) for  (auto &e:f)
    {
        auto c= e->compile_this(cm,fn);
        if (c!=nullptr)       return optional_unique_t<Compiled_Statement<Cm>>(std::move(c));
    }
    return {};
}

template <class Cm>
optional_unique_t<Compiled_Statement<Cm>> compile(const Cm* cm,const Assignment* a)
{

    if (a!=nullptr)
    {
        auto c=compile(cm,a->arg(0));
        if (c.has_value())
            return c.value().release()->compile_Assignment(a->id()->value());
        else return {};
    }
    else return  {};
}


template <class Cm>
optional_unique_t<Compiled_Statement<Cm>> compile(const Cm* cm,const Statement* id)

{
    if(auto x=dynamic_cast<const Identifier*>(id); x!=nullptr)
    {
        return cm->get_Identifier(x->value())->compile_this(cm,x);
    }
    else if(auto x=dynamic_cast<LiteralGeneric const *>(id); x!=nullptr)
    {
        return cm->get_Literal(x->myType())->compile_this(cm,x);
    }
    else if(auto x=dynamic_cast<Function const *>(id); x!=nullptr)
    {
        return compile(cm,x);
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
        return compile(cm,x);
    }
    else if(auto x=dynamic_cast<Definition const *>(id); x!=nullptr)
    {
        return comp::compiler<Cm>::get_Definition(cm,x);
        //get_Definition(x);
    }
    else return {};
}




template <class Cm,class T>
optional_unique_t<Compiled_General_Assignment<Cm,T>> compile(const Cm* cm ,C<T>,const Generic_Assignment* id)
{
    typedef  optional_unique_t<Compiled_General_Assignment<Cm,T>> res_type;

    if(auto x=dynamic_cast<Assignment const *>(id); x!=nullptr)
    {
        return res_type(Compiled_Assignment<Cm,T>::create(cm,x));
    }
    else if(auto x=dynamic_cast<Definition const *>(id); x!=nullptr)
    {
        return res_type(Compiled_Definition<Cm,T>::create(cm,x));
    }
    else return {};
}

















} // namespace grammar


#endif // MYCOMPILATION_H
