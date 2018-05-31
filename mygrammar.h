#ifndef MYGRAMMAR_H
#define MYGRAMMAR_H

#include "mySerializer.h"
#include "mysmartpointerstools.h"
#include <string>
#include <map>
#include <vector>
#include <strstream>

#include <type_traits>
#include <memory>



namespace grammar {


typedef io::token<'('> arg_start;
typedef io::token<','> arg_separator;

typedef io::token<')'> arg_end;

typedef io::token<';'> statement_end;

typedef io::token<'='> assigment_symbol;
typedef io::token<':','='> definition_symbol;


typedef  io::token<'{'> array_start;
typedef  io::token<'}'> array_end;
typedef  io::token<'\t'> array_separator;

typedef  io::token<'('> group_start;
typedef  io::token<')'> group_end;



typedef  io::token<'"'> label_start;
typedef  io::token<'"'> label_end;


class Statement
{
public:
    virtual ~Statement()=default;
   virtual std::string value()const=0;

    virtual std::size_t nArgin()const=0;
    virtual Statement const * arg(std::size_t i)const=0;
    virtual std::ostream& put(std::ostream& os)const=0;
    virtual bool get(std::stringstream& ss)=0;

    virtual Statement* clone() const=0;

};


class Expression : public Statement
{
public:
    virtual ~Expression(){}
     virtual Expression* clone() const override=0;
};

class Term : public Expression
{
public:
    virtual ~Term(){}
    virtual Term* clone() const override=0;
 };


inline bool get(std::stringstream& ss, Term*& e);
inline bool get(std::stringstream& ss, Expression*& e);
inline bool get(std::stringstream& ss, Statement*& e);


inline
myOptional<Statement*> string_to_Statement(const std::string& s)
{
    std::stringstream ss(s);
    Statement* sta;
    if (get(ss,sta))
        return sta;
    else return {};
}

class Identifier: public Term
{
    std::string id_;
    // Expression interface
public:
    static bool isValidIdentfier(const std::string& s)
    {
        if (!std::isalpha(s[0])) return false;
        if (s.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_0123456789")!=s.npos) return false;
        return true;
    }
    virtual std::string value() const override
    {
        return id_;
    }
    virtual std::size_t nArgin() const override
    {
        return 0;
    }
    virtual const Expression *arg(std::size_t ) const override
    {
        return nullptr;
    }

    virtual Identifier* clone() const override
    {
        return new Identifier(*this);
    }

    Identifier()=default;
    Identifier(const Identifier&)=default;
    Identifier(Identifier&&)=default;
    Identifier& operator=(const Identifier&)=default;
    Identifier& operator=(Identifier&&)=default;

    Identifier(std::string&& s): id_(s){}
    Identifier(const std::string& s): id_(s){}

    virtual ~Identifier()override{}
    // Expression interface
public:
    virtual std::ostream &put(std::ostream &os) const override
    {
        os<<value();
        return os;
    }



    virtual bool get(std::stringstream& ss)  override
    {
        auto pos=ss.tellg();
        std::string idCandidate;
        auto line=ss.str();
        ss>>idCandidate;
        if (ss.good()&&isValidIdentfier(idCandidate))
        {
            *this=Identifier(std::move(idCandidate));
            return true;
        }
        else
        {
            ss.clear();
            ss.seekg(pos);
            return false;
        }



    }
};





namespace get_variant_impl {





template<class...Ts>
bool get(std::stringstream& , std::variant<Ts...>& ,Cs<Ts...>, Cs<>)
{
    return false;
}






template<class... T0s, class T,class...Ts>
bool get(std::stringstream& ss, std::variant<T0s...,T,Ts...>& v, Cs<T0s...>,Cs<T,Ts...>)
{
    if constexpr (has_get_global<T>::value)
    {
        T e;
        if (get(ss,e))
        {
            v=std::move(e);
            return true;
        }
    }
    else if constexpr(std::is_pointer_v<T>)
    {
        std::unique_ptr<std::remove_pointer_t<T>> e(new  std::remove_pointer_t<T>{});
        if (e->get(ss))
        {
            v=e.release();
            return true;
        }
    }
    else
    {
        T e;
        if (e.get(ss))
        {
            v=std::move(e);
            return true;
        }
    }
    return get(ss,v,Cs<T0s...,T>{},Cs<Ts...>{});
}
} // namespace get_variant_impl

template<class...Ts>
bool get_variant(std::stringstream& ss, std::variant<Ts...>& v)
{
    return get_variant_impl::get(ss,v,Cs<>{},Cs<Ts...>{});
}


template<class Base>
bool get_Derived(std::stringstream&,  Base*&, Cs<>)
{
    return false;
}

template<class Base,class D,class...Ds>
bool get_Derived(std::stringstream& ss,Base*& x, Cs<D,Ds...>)
{
    if constexpr (has_get_global<D>::value&& ! std::is_same_v<Base,D >)
            //                                                |    |
            //                                      prevent infinite recursion when Base is concrete and it is its own child
    {
        D* e;
        if (get(ss,e))
            return true;
        else
            return get_Derived(ss,x,Cs<Ds...>{});

    }
    else
    {
        std::unique_ptr<D> e(new D{});
        if (e->get(ss))
        {
            x=e.release();
            return true;
        }
        else
        {
            return get_Derived(ss,x,Cs<Ds...>{});
        }
    }
}


class Operator
{
public:
    virtual std::string  name()const =0;
    virtual std::ostream &put(std::ostream &os) const=0;
    virtual bool get(std::stringstream& ss)=0;
    virtual std::size_t order()const=0;
    virtual Operator* clone()const=0;
    virtual ~Operator(){}
};



class AssignmentGenericOperator: public Operator{
public:
    virtual AssignmentGenericOperator* clone()const=0;
    virtual ~AssignmentGenericOperator(){}


};
class DefinitionGenericOperator: public AssignmentGenericOperator{
public:
    virtual DefinitionGenericOperator* clone()const=0;
    virtual ~DefinitionGenericOperator(){}


};

class Generic_Assignment: public Statement
{
public:
    ~Generic_Assignment(){}
    virtual Identifier const* id()const =0;
    virtual AssignmentGenericOperator const * op()const =0;
    virtual Expression const * expr()const=0;
    virtual Generic_Assignment* clone()const override=0;
};


inline bool get(std::stringstream& ss, Generic_Assignment*& e);


class Function: public Term
{
 public :

    typedef std::map<std::string,std::unique_ptr<Generic_Assignment>> arg_map;
private:
    Identifier idfn_;
    arg_map m_;

public:
    typedef arg_start start_symbol;
    typedef arg_end end_symbol;
    typedef arg_separator separator;

    std::pair<std::string, std::set<std::string>> getIdArg()const
    {
        std::pair <std::string,std::set<std::string>> out(idfn_.value(),{});
        for (auto& e: m_)
            out.second.insert(e.first);
        return out;
    }

    arg_map const & getArgs()const {return m_;}

    Function(Identifier id,arg_map a):idfn_{std::move(id)},m_{std::move(a)}
    {


    }
    Function(){}
    Function(const Function& other):
        idfn_{other.idfn_}, m_{clone_map(other.m_)}

    {

    }
    Function& operator=(const Function& other)
    {
        if (this!=&other)
        {
            Function temp(other);
            *this=std::move(temp);
        }
        return  *this;
    }
    Function(Function&&)=default;
    Function& operator=(Function&&)=default;
    virtual Function* clone()const override
    {
        return new Function(*this);
    }
    virtual ~Function()override{}

    std::string value()const override
    {
        return idfn_.value();
    }
    std::size_t nArgin()const override
    {
        return m_.size();
    }
    Generic_Assignment const * arg(std::size_t i)const override{
        auto it=m_.begin();
        for (std::size_t j=0; j<i; ++j) ++it;
        return it->second.get();
    }

    virtual std::ostream &put(std::ostream &os) const override;
    static bool getArguments(std::stringstream& ss, arg_map &out);

    virtual bool get(std::stringstream &is) override;

};



template <class Optype,typename symbol>
class Operator_symbol: public Optype
{
public:
    virtual std::ostream &put(std::ostream &os) const override
    {
        os<<symbol{};
        return os;
    }

    virtual bool get(std::stringstream& ss)override
    {
        auto pos=ss.tellg();
        if (ss>>symbol{})
            return true;
        else
        {
            ss.clear();
            ss.seekg(pos);
            return false;
        }
    }
};

class UnaryOperator: public Operator{

public:
    virtual UnaryOperator* clone()const=0;
    virtual ~UnaryOperator(){}


};
class BinaryOperator: public Operator{

public:
    virtual BinaryOperator* clone()const=0;
    virtual ~BinaryOperator(){}


};

template <class T>
class UnaryOperatorTyped: public UnaryOperator{
public:
    virtual T operator()(T A)const =0;
    virtual UnaryOperatorTyped* clone()const=0;
    virtual ~UnaryOperatorTyped(){}

};

template <class T>
struct BinaryOperatorTyped: public BinaryOperator{

    virtual T operator()(T A, T B)const =0;
    virtual BinaryOperatorTyped* clone()const=0;
    virtual ~BinaryOperatorTyped(){}

};

template <class Cm,class T>
struct AssignmentOperatorTyped: public Operator{

  virtual  T operator()(Cm* cm, std::string id, T x)const =0;
    virtual AssignmentOperatorTyped* clone()const=0;
    virtual ~AssignmentOperatorTyped(){}

};





inline bool get(std::stringstream& ss, BinaryOperator*&);
inline bool get(std::stringstream& ss, UnaryOperator*&);
inline bool get(std::stringstream& ss, AssignmentGenericOperator*&);
inline bool get(std::stringstream& ss, DefinitionGenericOperator*&);
inline bool get(std::stringstream& ss, Identifier*&);




class AssignmentOperator: public Operator_symbol<AssignmentGenericOperator,assigment_symbol>
{
public:
    static constexpr auto className(){return  my_static_string("Assignment_operator");}
    virtual std::string name() const override {return className().c_str();}
    virtual std::size_t order()const override {return 14;}
    virtual AssignmentOperator* clone()const override{return new AssignmentOperator;}
    virtual ~AssignmentOperator(){}
    AssignmentOperator()=default;



};
class DefinitionOperator: public Operator_symbol<DefinitionGenericOperator,definition_symbol>
{
public:
    static constexpr auto className(){return  my_static_string("Definition_operator");}
    virtual std::string name() const override {return className().c_str();}
    virtual std::size_t order()const override {return 14;}
    virtual DefinitionOperator* clone()const override{return new DefinitionOperator;}
    virtual ~DefinitionOperator(){}
    DefinitionOperator()=default;
};




class UnaryOperation: public Term
{
private:
    std::unique_ptr<UnaryOperator> op_;
    std::unique_ptr<Term> a_;

public:
    explicit UnaryOperation(UnaryOperator* op,Term* e): op_{op},a_{e}{}
    UnaryOperation(){}
    UnaryOperator const * op()const
    {return  op_.get();}
    virtual std::string value() const override;
    virtual std::size_t nArgin() const override
    {
        return 1;
    }
    virtual const Term *arg(std::size_t ) const override
    {
        return a_.get();
    }

    // Expression interface
public:
    virtual std::ostream &put(std::ostream &os) const override;
    virtual bool get(std::stringstream &ss)  override;
    virtual UnaryOperation* clone()const override
    {return new UnaryOperation(*this);}
    virtual ~UnaryOperation(){}
    UnaryOperation(const UnaryOperation& other)
        :op_{other.op_->clone()},a_{other.a_->clone()}{}
    UnaryOperation& operator=(const UnaryOperation& other)
    {
        if(this!=&other)
        {
            UnaryOperation tmp(other);
            *this=std::move(tmp);

        }
        return *this;
    }
    UnaryOperation(UnaryOperation&& other)=default;
    UnaryOperation& operator=(UnaryOperation&& other)=default;

};

class BinaryOperation: public Expression
{
    std::unique_ptr<BinaryOperator> op_;
    std::unique_ptr<Expression> la_;
    std::unique_ptr<Expression> ra_;


public:
    explicit BinaryOperation(Expression* e1,BinaryOperator* op, Expression* e2): op_{op},la_{e1},ra_{e2}{}
    BinaryOperation(){}
    BinaryOperator const * op()const{return op_.get();}
    Expression const * larg()const {return la_.get();}
    Expression const * rarg()const {return ra_.get();}


    // Expression interface
public:
    virtual std::string value() const override;
    virtual std::size_t nArgin() const override
    {
        return 2;
    }
    virtual const Expression *arg(std::size_t i) const override
    {
        switch (i) {
        case 0: return larg();
            break;
        case 1: return rarg();
        default:
            return nullptr;

        }
    }
    virtual std::ostream &put(std::ostream &os) const override;
    virtual bool get(std::stringstream& ss) override;
    virtual BinaryOperation* clone()const override
    {return new BinaryOperation(*this);}
    virtual ~BinaryOperation(){}
    BinaryOperation(const BinaryOperation& other)
        :op_{other.op_->clone()},la_{other.la_->clone()},ra_{other.ra_->clone()}{}
    BinaryOperation& operator=(const BinaryOperation& other)
    {
        if(this!=&other)
        {
            BinaryOperation tmp(other);
            *this=std::move(tmp);

        }
        return *this;
    }
    BinaryOperation(BinaryOperation&& other)=default;
    BinaryOperation& operator=(BinaryOperation&& other)=default;


};


class ArrayOperation: public Term
{
    std::vector<std::unique_ptr<Expression>> v_;
    // Expression interface
public:
    ArrayOperation(){}
    ArrayOperation(std::vector<std::unique_ptr<Expression>>&& v):v_{std::move(v)}{}
    static constexpr const char* className(){ return "ArrayOperation";}

    virtual std::string value() const override
    {
        return className();
    }
    virtual std::size_t nArgin() const override
    {
        return v_.size();
    }
    virtual const Expression *arg(std::size_t i) const override
    {
        return v_[i].get();
    }
    virtual std::ostream &put(std::ostream &os) const override
    {
        os<<array_start{};
        for (std::size_t i=0; i<v_.size(); ++i)
        {
            v_[i]->put(os);
            if (i+1<v_.size())
                os<<array_separator{};
        }
        os<<array_end{};
        return os;
    }


    // Expression interface
public:
    virtual bool get(std::stringstream &ss) override
    {
        auto pos=ss.tellg();
        auto line=ss.str();
        if (ss>>array_start{})
        {
            std::vector<std::unique_ptr<Expression>> out;
            std::variant<array_end,Expression*> v;
            std::string s;
            bool succeed;
            while (succeed=get_variant(ss,v),succeed&&(v.index()==1))
            {
                out.emplace_back(std::get<Expression*>(v));
            }
            if (succeed)
            {
                *this=ArrayOperation(std::move(out));
                return  true;
            }
        }
        ss.clear();
        ss.seekg(pos);
        return false;
    }

    virtual ArrayOperation* clone()const override
    {return new ArrayOperation(*this);}
    virtual ~ArrayOperation(){}
    ArrayOperation(const ArrayOperation& other)
        :v_{clone_vector(other.v_)}{}
    ArrayOperation& operator=(const ArrayOperation& other)
    {
        if(this!=&other)
        {
            ArrayOperation tmp(other);
            *this=std::move(tmp);

        }
        return *this;
    }
    ArrayOperation(ArrayOperation&& other)=default;
    ArrayOperation& operator=(ArrayOperation&& other)=default;



};





class GroupOperation: public Term
{
    std::unique_ptr<Expression> e_;
    // Expression interface
public:
    typedef  group_start start_symbol;
    typedef  group_end end_symbol;

    GroupOperation(){}
    GroupOperation(Expression* e):e_{e}{}
    static constexpr const char* className(){ return "GroupOperation";}

    virtual std::string value() const override
    {
        return className();
    }
    virtual std::size_t nArgin() const override
    {
        return 1;
    }
    virtual const Expression *arg(std::size_t) const override
    {
        return e_.get();
    }
    virtual std::ostream &put(std::ostream &os) const override
    {
        os<<start_symbol{}<<e_.get()<<end_symbol{};
        return os;
    }


    // Expression interface
public:
    virtual bool get(std::stringstream &ss) override
    {
        Expression* e;
        auto pos=ss.tellg();
        if ((ss>>group_start{})&&grammar::get(ss,e)&&(ss>>group_end{}))
        {
            *this=GroupOperation(e);
            return  true;
        }
        ss.clear();
        ss.seekg(pos);
        return false;
    }

    virtual GroupOperation* clone()const override
    {return new GroupOperation(*this);}
    virtual ~GroupOperation(){}
    GroupOperation(const GroupOperation& other)
        :e_{other.e_->clone()}{}
    GroupOperation& operator=(const GroupOperation& other)
    {
        if(this!=&other)
        {
            GroupOperation tmp(other);
            *this=std::move(tmp);

        }
        return *this;
    }
    GroupOperation(GroupOperation&& other)=default;
    GroupOperation& operator=(GroupOperation&& other)=default;

};




template<bool B,class AssignmentGenericOperator>
class basic_Assignment: public Generic_Assignment
{
private:
    std::unique_ptr<Identifier> id_;
    std::unique_ptr<AssignmentGenericOperator> op_;
    std::unique_ptr<Expression> exp_;
public:
    constexpr static bool delayed=B;
    basic_Assignment(){}
    basic_Assignment(Identifier* id, AssignmentGenericOperator* op, Expression* exp):
        id_{id},op_{op},exp_{exp}{}
    Identifier const* id()const override {return id_.get();}
    AssignmentGenericOperator const * op()const override {return op_.get();}
    Expression const * expr()const override {return exp_.get();}


    virtual std::ostream &put(std::ostream &os) const override;
    virtual bool get(std::stringstream &ss) override;

    // Expression interface
public:
    virtual std::string value() const override;
    virtual std::size_t nArgin() const override
    {
        return 1;
    }
    virtual const Expression *arg(std::size_t ) const override
    {
        return expr();
    }
    virtual basic_Assignment* clone()const override
    {return new basic_Assignment(*this);}
    virtual ~basic_Assignment(){}
    basic_Assignment(const basic_Assignment& other)
        :id_{other.id_->clone()},op_{other.op_->clone()},exp_{other.exp_->clone()}{}
    basic_Assignment& operator=(const basic_Assignment& other)
    {
        if(this!=&other)
        {
            basic_Assignment tmp(other);
            *this=std::move(tmp);

        }
        return *this;
    }
    basic_Assignment(basic_Assignment&& other)=default;
    basic_Assignment& operator=(basic_Assignment&& other)=default;




};

typedef basic_Assignment<false,AssignmentGenericOperator> Assignment;

typedef basic_Assignment<true,DefinitionGenericOperator> Definition;



class LiteralGeneric: public Term
{
public:
    virtual ~LiteralGeneric();
    virtual LiteralGeneric* clone()const override=0;

};
LiteralGeneric::~LiteralGeneric()=default;

template<typename T>
class Literal: public LiteralGeneric
{
    T value_;
public:
    //Literal(T&& value):value_{std::move(value)}{}
    Literal(const T& value):value_{value}{}
    Literal(){}
    T getValue()const {return value_;}
    virtual std::string value() const override
    {
        return my_to_string(value_);
    }
    virtual std::size_t nArgin() const override
    {
        return 0;
    }
    virtual const Assignment *arg(std::size_t ) const override
    {
        return nullptr;
    }

    virtual std::ostream &put(std::ostream &os) const override
    {
        os<<value_;
        return  os;
    }
    virtual bool get(std::stringstream &ss) override
    {
        auto pos=ss.tellg();
        T value;
        if (ss>>value)
        {
            value_=std::move(value);
            return true;
        }
        else
        {
            ss.clear();
            ss.seekg(pos);
            return false;
        }
    }

    virtual Literal* clone()const override { return new Literal(*this);}

    Literal(const Literal&)=default;
    Literal(Literal&&)=default;
    Literal& operator=(const Literal&)=default;
    Literal& operator=(Literal&&)=default;

};


template<>
class Literal<std::string>: public LiteralGeneric
{
    std::string value_;
public:
    explicit Literal(std::string&& out): value_{std::move(out)}{}
    Literal(){}
    std::string getValue()const {return value_;}
    virtual std::string value() const override
    {
        return value_;
    }
    virtual std::size_t nArgin() const override
    {
        return 0;
    }
    virtual const Assignment *arg(std::size_t ) const override
    {
        return nullptr;
    }

    virtual std::ostream &put(std::ostream &os) const override
    {
        return os<<label_start{}<<value_<<label_end{};
    }
    virtual bool get(std::stringstream &ss) override
    {
        auto line=ss.str();
        auto pos=ss.tellg();
        if (ss>>label_start{})
        {
            std::string out;
            std::string end=label_end{}.value;
            auto n=end.size();
            char c;
            ss.get(c);
            out.push_back(c);
            auto p=out.find(end,std::max(n,out.size())-n);
            while(p==out.npos)
            {

                if (!(ss.get(c))) break;
                out.push_back(c);
                p=out.find(end,std::max(n,out.size())-n);
            }
            if(ss.good())
            {
                out.erase(p);
                *this=Literal(std::move(out));
                return  true;
            }
        }
        ss.clear();
        ss.seekg(pos);
        return false;
    }
    virtual Literal* clone()const override { return new Literal(*this);}

    Literal(const Literal&)=default;
    Literal(Literal&&)=default;
    Literal& operator=(const Literal&)=default;
    Literal& operator=(Literal&&)=default;

};









std::string UnaryOperation::value() const
{
    return op()->name();

}

std::ostream &UnaryOperation::put(std::ostream &os) const
{
    op()->put(os);
    arg(0)->put(os);
    return os;
}

bool UnaryOperation::get(std::stringstream &ss)
{
    auto pos=ss.tellg();
    UnaryOperator* op;
    Term* e;
    if (grammar::get(ss,op)&&grammar::get(ss,e))
    {
        *this=UnaryOperation(op,e);
        return true;
    }
    else
    {
        ss.clear();
        ss.seekg(pos);
        return false;
    }
}




std::ostream & Function::put(std::ostream &os) const
{
    os<<value()<<arg_start{};
    for (std::size_t i=0; i<nArgin(); ++i )
    {
        arg(i)->put(os);
        if (i<nArgin()) os<<arg_separator{};
    }
    os<<arg_end{};
    return  os;

}


bool Function::get(std::stringstream &ss)
{
    auto pos=ss.tellg();
    Identifier id;
    arg_map out;

    if (id.get(ss)&&(ss>>arg_start{})&&getArguments(ss,out))
    {
        *this=std::move(Function(std::move(id),std::move(out)));
        return  true;
    }
    else{
        ss.clear();
        ss.seekg(pos);
        return false;
    }
}
bool Function::getArguments(std::stringstream &ss, arg_map& out)
{
    std::variant<arg_end,Generic_Assignment*> v;
    bool succeed;
    while (succeed=get_variant(ss,v), succeed&&(v.index()==1))
    {
        auto a=std::get<Generic_Assignment*>(v);
        out.emplace(a->id()->value(),a);
    }
    if (succeed)
    {
        return  true;
    }
    return false;

}


template<bool B,class AssignmentGenericOperator>
std::ostream &basic_Assignment<B,AssignmentGenericOperator>::put(std::ostream &os) const
{
    id()->put(os);
    op()->put(os);
    expr()->put(os);
    return os;

}


template<bool B,class AssignmentGenericOperator>
bool basic_Assignment<B,AssignmentGenericOperator>::get(std::stringstream &ss)
{
    Identifier* id;
    AssignmentGenericOperator* a;
    Expression* e;
    auto pos=ss.tellg();
    if (grammar::get(ss,id)&& grammar::get(ss,a)
            &&grammar::get(ss,e))
    {
        *this=basic_Assignment(id,a,e);
        return true;
    }
    else
    {
        ss.clear();
        ss.seekg(pos);
        return false;
    }

}



template<bool B,class AssignmentGenericOperator>
std::string basic_Assignment<B,AssignmentGenericOperator>
 ::value() const
{
    return op()->name();
}

typedef std::variant<LiteralGeneric*,GroupOperation*,ArrayOperation*,UnaryOperator*,Identifier*> valid_start;
typedef  std::variant< BinaryOperator*, AssignmentGenericOperator*, DefinitionGenericOperator*> after_idenfier;

inline
bool resolve_operator_order(std::stringstream& ss, BinaryOperation*& b,Expression* e0, BinaryOperator* op0, Expression* e1)
{
    BinaryOperator* op1;
    if(get(ss,op1))
    {
        Term *e2;
        if (!get(ss,e2))
            return false;
        if (op0->order()<=op1->order())
        {
            auto e0_op0_e1=new BinaryOperation(e0,op0,e1);
            return resolve_operator_order(ss,b,e0_op0_e1,op1,e2);
        }
        else
        {
            BinaryOperation* e1_op1_e2;
            if (resolve_operator_order(ss,e1_op1_e2,e1,op1,e2))
            {
                b=new BinaryOperation(e0,op0,e1_op1_e2);
                return true;
            }
            else return false;
        }
    }
    b=new BinaryOperation(e0,op0,e1);
    return true;
}



inline
bool get_term_from_identifier(std::stringstream& ss, Term*& t, Identifier* id)
{
    if (get(ss,Function::start_symbol{}))
    {
        Function::arg_map out;
        if (!Function::getArguments(ss,out))
            return false;
        else {
            t=new Function(std::move(*id), std::move(out));
            return true;
        }
    }
    else
    {
        t=id;
        return true;
    }

}

inline
bool get_term_from_the_rest(std::stringstream& ss, Term*& e, valid_start& first)
{
    if (std::holds_alternative<UnaryOperator*>(first))
    {
        auto unary=std::get<UnaryOperator*>(first);
        Term *e2;
        if (!get(ss,e2))  {return false;}
        else
        {
            e=new UnaryOperation(unary,e2); return true;
        }
    }
    else if (std::holds_alternative<LiteralGeneric*>(first))
    {e=std::get<LiteralGeneric*>(first); return true;}
    else if (std::holds_alternative<GroupOperation*>(first))
    { e=std::get<GroupOperation*>(first); return true;}
    else if (std::holds_alternative<ArrayOperation*>(first))
    { e=std::get<ArrayOperation*>(first); return true;}
    else
        // something very wrong
    {
        return false;
    }


}



inline
bool get(std::stringstream& ss, Term*& e)
{
    auto pos=ss.tellg();
    auto line=ss.str();
    valid_start first;
    if (!get_variant(ss,first))
    { ss.clear(); ss.seekg(pos); return false; }
    else
    {
        if (std::holds_alternative<Identifier*>(first))
        {
            auto id=std::get<Identifier*>(first);
            if (get_term_from_identifier(ss,e,id))
                return true;
            else {ss.clear(); ss.seekg(pos); return false;}
        }
        else if (get_term_from_the_rest(ss, e,first))
            return true;
        else {ss.clear(); ss.seekg(pos); return false;}
    }
}


inline
bool get_expression_from_term(std::stringstream& ss, Expression*& e, Term* t)
{
    BinaryOperator* op;
    if(!get(ss,op))
    {
        e=t;
        return true;
    }
    else
    {
        Term *t2;
        if (!get(ss,t2))
            return false;
        else
        {
            BinaryOperation* b;
            if (resolve_operator_order(ss,b,t,op,t2))
            {
                e=b;
                return true;
            }
            else
                return false;
        }
    }
}



inline
bool get(std::stringstream& ss, Expression*& e)
{
    auto pos=ss.tellg();
    Term* t;
    if (get(ss,t)&&get_expression_from_term(ss,e,t))
        return true;
    else
    { ss.clear(); ss.seekg(pos); return false; }
}

inline
bool get(std::stringstream& ss, Statement*& sta)
{
    auto pos=ss.tellg();
    Term *t;
    Expression *e;
    valid_start first;
    if (!get_variant(ss,first))
    { ss.clear(); ss.seekg(pos); return false; }
    else
    {
        if (std::holds_alternative<Identifier*>(first))
        {
            auto id=std::get<Identifier*>(first);
            AssignmentGenericOperator* opa;
            DefinitionGenericOperator* opd;
            if (get(ss,opa))
            {
                Expression *e2;
                if (!get(ss,e2))
                {ss.clear(); ss.seekg(pos); return false;}
                else {
                    sta=new Assignment(id, opa,e2);
                    return true;
                }
            }
            else if (get(ss,opd))
            {
                Expression *e2;
                if (!get(ss,e2))
                {ss.clear(); ss.seekg(pos); return false;}
                else {
                    sta=new Definition(id, opd,e2);
                    return true;
                }
            }
            else if (get_term_from_identifier(ss,t,id)&&get_expression_from_term(ss,e,t))
            {
                sta=e;
                return true;
            }
            else {ss.clear(); ss.seekg(pos); return false;}
        }
        else if (get_term_from_the_rest(ss, t,first)&&get_expression_from_term(ss,e,t))
        {
            sta=e;
            return true;
        }
        else {ss.clear(); ss.seekg(pos); return false;}

    }

}



inline
bool get(std::stringstream& ss, Generic_Assignment*& a)
{
    auto pos=ss.tellg();
    Identifier* id;
    Expression* e;
    std::variant<AssignmentGenericOperator*,DefinitionGenericOperator*> op;
    if (get(ss,id)&&get_variant(ss,op)&&get(ss,e))
    {
        if (std::holds_alternative<AssignmentGenericOperator*>(op))
        {
            auto opa=std::get<AssignmentGenericOperator*>(op);
            a= new Assignment(id,opa,e);
            return true;
        }
        else if (std::holds_alternative<DefinitionGenericOperator*>(op))
        {
            auto opd=std::get<DefinitionGenericOperator*>(op);
            a= new Definition(id,opd,e);
            return true;

        }

    }
    ss.clear();
    ss.seekg(pos);
    return false;
}
bool get(std::stringstream& ss, Identifier*& e)
{
    return get_Derived(ss,e,Cs<Identifier>{});
}

inline
bool get(std::stringstream& ss, AssignmentGenericOperator*& e)
{
    return get_Derived(ss,e,Cs<AssignmentOperator>{});
}

inline
bool get(std::stringstream& ss, DefinitionGenericOperator*& e)
{
    return get_Derived(ss,e,Cs<DefinitionOperator>{});
}

inline bool get(std::stringstream& ss, UnaryOperator*& e)
{
    return get_Derived(ss,e,Cs<>{});
}
inline bool get(std::stringstream& ss, BinaryOperator*& e)
{
    return get_Derived(ss,e,Cs<>{});
}




inline bool get(std::stringstream& ss, LiteralGeneric*& e)
{
    return get_Derived(ss,e,Cs<Literal<std::string>, Literal<double>, Literal<std::size_t>, Literal<int>, Literal<bool>>{});
}



std::string BinaryOperation::value() const
{
    return op()->name();
}
std::ostream &BinaryOperation::put(std::ostream &os) const
{
    arg(0)->put(os);
    op()->put(os);
    arg(1)->put(os);
    return os;

}

bool BinaryOperation::get(std::stringstream &ss)
{
    auto pos=ss.tellg();
    BinaryOperator* op;
    Term* t1;
    Term* t2;

    if (grammar::get(ss,t1)&&grammar::get(ss,op)&&grammar::get(ss,t2))
    {
        BinaryOperation* b;
        if (resolve_operator_order(ss,b,t1,op,t2))
        {
            *this=std::move(*b);
            delete b;
            return true;
        }
        else
        {
            ss.clear();
            ss.seekg(pos);
            return false;
        }
    }
    else
    {
        ss.clear();
        ss.seekg(pos);
        return false;
    }
}





} // namespace grammar



#endif // MYGRAMMAR_H

