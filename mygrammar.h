#ifndef MYGRAMMAR_H
#define MYGRAMMAR_H

#include "mySerializer.h"
#include "mysmartpointerstools.h"
#include <string>
#include <map>
#include <vector>

#include <type_traits>
#include <memory>
#include <random>


namespace grammar {



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


template <typename Pos>
std::string  print_pos(const Pos& pos) { return "@ pos="+std::to_string(pos)+" :";}

class Statement
{
public:
  virtual ~Statement(){}
    virtual std::string value()const=0;

    virtual std::size_t nArgin()const=0;
    virtual Statement const * arg(std::size_t i)const=0;
    virtual std::ostream& put(std::ostream& os)const=0;
    virtual myOptional_t<void> get(std::stringstream& ss)=0;

    virtual Statement* clone() const=0;

};


class Expression : public Statement
{
public:
    typedef Statement base_type;
    virtual ~Expression()override{}
    virtual Expression* clone() const override=0;
};

class Term : public Expression
{
public:
    typedef Expression base_type;

    virtual ~Term()override{}
    virtual Term* clone() const override=0;
};

template<typename T>
myOptional_t<T> get(C<T>, std::stringstream& ss);

inline  myOptional_t<Term*> get(C<Term*>,std::stringstream& ss);
inline  myOptional_t<Expression*> get(C<Expression*>,std::stringstream& ss);
inline  myOptional_t<Statement*> get(C<Statement*>,std::stringstream& ss);

template<char... cs>
inline  myOptional_t<io::token<cs...>> get(C<io::token<cs...>>,std::stringstream& ss)
{
    std::string line=ss.str();

    io::token<cs...> t;
    if (!ss>>t)
        return {false," invalid token"};
    else return t;
}



inline
myOptional_t<Statement*> string_to_Statement(const std::string& s)
{
    std::stringstream ss(s);
    return get(C<Statement*>{},ss);
}

class Identifier: public Term
{
    std::string id_;
    // Expression interface
public:
    constexpr static char valid_chars[]="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_0123456789";

    static  bool valid_char(char ch) { return std::string(valid_chars).find_first_of(ch)!=std::string::npos; }

    typedef Term base_type;
    static bool isValidIdentfier(const std::string& s)
    {
        if (!std::isalpha(s[0])) return false;
        if (s.find_first_not_of(valid_chars)!=s.npos) return false;
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

    Identifier(std::string&& s): id_(std::move(s)){}
    Identifier(const std::string& s): id_(s){}

    virtual ~Identifier()override{}
    // Expression interface
public:
    virtual std::ostream &put(std::ostream &os) const override
    {
        os<<value();
        return os;
    }



    virtual myOptional_t<void> get(std::stringstream& ss)  override
    {
        auto pos=ss.tellg();
        std::string idCandidate;
        auto line=ss.str();
        char ch;
        while (ss.get(ch))
        {
            if (valid_char(ch)) idCandidate.push_back(ch);
            else if (!std::isspace(ch)||(idCandidate.size()>0))
                    {ss.putback(ch); break;}
        }
        if ((!ss.good())||(!isValidIdentfier(idCandidate)))
        {
            ss.clear();
            ss.seekg(pos);
            if (!ss.good())
                return {false,print_pos(pos)+"not a string"};
            else return {false,print_pos(pos)+"not a valid candidate"};
        }
        else
        {
            *this=Identifier(std::move(idCandidate));
            return {true,""};
        }

    }
};







namespace get_variant_templ_impl {





template<class...Ts>
myOptional_t<std::variant<Ts...>> get(Cs<Ts...>,Cs<>,std::stringstream& ss, std::string error);

template<class... T0s, class T,class...Ts>
myOptional_t<std::variant<T0s...,T,Ts...>> get(Cs<T0s...>,Cs<T,Ts...>,std::stringstream& ss, std::string error);

} // namespace get_variant_impl




template<class...Ts>
myOptional_t<std::variant<Ts...>> get_variant(C<std::variant<Ts...>>,std::stringstream& ss)
{
    auto pos=ss.tellg();
    std::string line=ss.str().substr(pos);
    return get_variant_templ_impl::get(Cs<>{},Cs<Ts...>{},ss,"");

}

template<class Base>
myOptional_t<Base*>  get_Derived(Cs<Base>,std::stringstream&, std::string error)
{
    return {false,error};
}

template<class Base,class D,class...Ds>
myOptional_t<Base*> get_Derived(Cs<Base,D,Ds...>,std::stringstream& ss, std::string error="")
{
    if constexpr (has_get_global<D>::value&& ! std::is_same_v<Base,D >)
            //                                                |    |
            //                                      prevent infinite recursion when Base is concrete and it is its own child
    {
        D* e;
        auto o=get(C<D*>{},ss);
        if (o)
            return o.release();
        else
            return get_Derived(Cs<Base,Ds...>{},ss, error+", not : "+o.error());

    }
    else
    {
        std::unique_ptr<D> e(new D{});
        auto o=e->get(ss);
        if (o)
        {
            return e.release();
        }
        else
        {
            return get_Derived(Cs<Base,Ds...>{},ss, error+", not : "+o.error());
        }
    }
}







class Operator
{
public:
    virtual std::string  name()const =0;
    virtual std::ostream &put(std::ostream &os) const=0;
    virtual myOptional_t<void> get(std::stringstream& ss)=0;
    virtual std::size_t order()const=0;
    virtual Operator* clone()const=0;
    virtual ~Operator(){}
};



class AssignmentGenericOperator: public Operator{
public:
    typedef Operator base_type;

    virtual AssignmentGenericOperator* clone() const override=0;
    virtual ~AssignmentGenericOperator()override{}


};
class DefinitionGenericOperator: public AssignmentGenericOperator{
public:
    typedef AssignmentGenericOperator base_type;

    virtual DefinitionGenericOperator* clone()const override=0;
    virtual ~DefinitionGenericOperator()override{}


};

class Generic_Assignment: public Statement
{
public:
    typedef Statement base_type;

    ~Generic_Assignment(){}
    virtual Identifier const* id()const =0;
    virtual AssignmentGenericOperator const * op()const =0;
    virtual Expression const * expr()const=0;
    virtual Generic_Assignment* clone()const override=0;
};



inline myOptional_t<Generic_Assignment*> get(C<Generic_Assignment>,std::stringstream& ss);


class Function: public Term
{
public :
    typedef Term base_type;

    typedef std::map<std::string,std::unique_ptr<Generic_Assignment>> arg_map;
private:
    Identifier idfn_;
    arg_map m_;

public:
    typedef arg_start start_symbol;
    typedef arg_end end_symbol;
    typedef arg_separator separator;

    std::pair<std::string, std::map<std::string, std::string>> getIdArg()const
    {
        std::pair <std::string,std::map<std::string, std::string>> out(idfn_.value(),{});
        for (auto& e: m_)
            out.second.emplace(e.first,"");
        return out;
    }

    arg_map const & getArgs()const {return m_;}

    Function(Identifier&& id,arg_map&& a):idfn_{std::move(id)},m_{std::move(a)}
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

    static myOptional_t<arg_map> getArguments(C<arg_map>,std::stringstream& ss);
    virtual myOptional_t<void> get(std::stringstream &is) override;

};



template <class Optype,typename symbol>
class Operator_symbol: public Optype
{
public:
    typedef Optype base_type;

    virtual ~Operator_symbol()override{}

    virtual std::ostream &put(std::ostream &os) const override
    {
        os<<symbol{};
        return os;
    }

    virtual myOptional_t<void> get(std::stringstream& ss)override
    {
        auto pos=ss.tellg();
        std::string line=ss.str().substr(pos);
        if (ss>>symbol{})
            return {true,""};
        else
        {
            ss.clear();
            ss.seekg(pos);
            return {false, print_pos(pos)+"operator symbol not found"};
        }
    }
};

class UnaryOperator: public Operator{

public:
    typedef Operator base_type;

    virtual UnaryOperator* clone()const override=0;
    virtual ~UnaryOperator()override{}


};
class BinaryOperator: public Operator{

public:
    typedef Operator base_type;

    virtual BinaryOperator* clone()const override=0;
    virtual ~BinaryOperator()override{}


};

template <class T>
class UnaryOperatorTyped: public UnaryOperator{
public:
    typedef UnaryOperator base_type;

    virtual T operator()(T A)const =0;
    virtual UnaryOperatorTyped* clone()const=0;
    virtual ~UnaryOperatorTyped()override{}

};

template <class T>
struct BinaryOperatorTyped: public BinaryOperator{
    typedef BinaryOperator base_type;

    virtual T operator()(T A, T B)const =0;
    virtual BinaryOperatorTyped* clone()const=0;
    virtual ~BinaryOperatorTyped()override{}

};

template <class Cm,class T>
struct AssignmentOperatorTyped: public Operator{
    typedef Operator base_type;

    virtual  T operator()(Cm* cm, std::string id, T x)const =0;
    virtual AssignmentOperatorTyped* clone()const=0;
    virtual ~AssignmentOperatorTyped()override{}

};





inline myOptional_t<BinaryOperator*> get(C<BinaryOperator*>,std::stringstream& ss);
inline myOptional_t<UnaryOperator*> get(C<UnaryOperator*>,std::stringstream& ss);
inline myOptional_t<AssignmentGenericOperator*> get(C<AssignmentGenericOperator*>,std::stringstream& ss);
inline myOptional_t<DefinitionGenericOperator*> get(C<DefinitionGenericOperator*>,std::stringstream& ss);
inline myOptional_t<Identifier*> get(C<Identifier*>,std::stringstream& ss);




class AssignmentOperator: public Operator_symbol<AssignmentGenericOperator,assigment_symbol>
{
public:
    typedef Operator_symbol<AssignmentGenericOperator,assigment_symbol> base_type;

    static constexpr auto className(){return  my_static_string("Assignment_operator");}
    virtual std::string name() const override {return className().c_str();}
    virtual std::size_t order()const override {return 14;}
    virtual AssignmentOperator* clone()const override{return new AssignmentOperator;}
    virtual ~AssignmentOperator()override{}
    AssignmentOperator()=default;



};
class DefinitionOperator: public Operator_symbol<DefinitionGenericOperator,definition_symbol>
{
public:
    typedef Operator_symbol<DefinitionGenericOperator,definition_symbol> base_type;

    static constexpr auto className(){return  my_static_string("Definition_operator");}
    virtual std::string name() const override {return className().c_str();}
    virtual std::size_t order()const override {return 14;}
    virtual DefinitionOperator* clone()const override{return new DefinitionOperator;}
    virtual ~DefinitionOperator()override{}
    DefinitionOperator()=default;
};




class UnaryOperation: public Term
{
private:
    std::unique_ptr<UnaryOperator> op_;
    std::unique_ptr<Term> a_;

public:

    typedef Term base_type;
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
    virtual myOptional_t<void> get(std::stringstream &ss)  override;
    virtual UnaryOperation* clone()const override
    {return new UnaryOperation(*this);}
    virtual ~UnaryOperation()override{}
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

    typedef Expression base_type;
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
    virtual myOptional_t<void> get(std::stringstream& ss) override;
    virtual BinaryOperation* clone()const override
    {return new BinaryOperation(*this);}
    virtual ~BinaryOperation()override{}
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
    typedef Term base_type;

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
    virtual myOptional_t<void> get(std::stringstream &ss) override
    {
        auto pos=ss.tellg();
        auto line=ss.str();
        if (!(ss>>array_start{}))
        {
            ss.clear();
            ss.seekg(pos);
            return {false,print_pos(pos)+"lacks an array start"};
        }
        else
        {
            std::vector<std::unique_ptr<Expression>> out;
            //std::variant<array_end,Expression*> v;
            std::string s;
            while(true)
            {
                auto v=get_variant(C<std::variant<array_end,Expression*>>{},ss);
                if (! v)
                {
                    ss.clear();
                    ss.seekg(pos);
                    return {false, print_pos(pos)+"invalid input at the "+std::to_string(out.size())+"th place of the array: "+v.error()};
                }
                else if (v.value().index()==1)
                    out.emplace_back(std::get<Expression*>(v.value()));
                else
                {
                    *this=ArrayOperation(std::move(out));
                    return  {true,""};

                }
            }
        }
    }

    virtual ArrayOperation* clone()const override
    {return new ArrayOperation(*this);}
    virtual ~ArrayOperation()override{}
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
    typedef Term base_type;

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
    virtual myOptional_t<void> get(std::stringstream &ss) override
    {
        auto pos=ss.tellg();
        if (!(ss>>group_start{}))
        {
            ss.clear(); ss.seekg(pos); return {false,print_pos(pos)+"lacks group start"};
        }
        else
        {
            auto o=grammar::get(C<Expression*>{},ss);
            if  (!o)
            {
                ss.clear(); ss.seekg(pos); return {false,print_pos(pos)+"error inside group: "+o.error()};
            }
            else if (!(ss>>group_end{}))
            {
                ss.clear(); ss.seekg(pos); return {false,print_pos(pos)+"group end is lacking"};
            }
            else
            {
                *this=GroupOperation(o.release());
                return  {true,""};
            }
        }
    }
    virtual GroupOperation* clone()const override
    {return new GroupOperation(*this);}
    virtual ~GroupOperation()override{}
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

inline myOptional_t<GroupOperation*> get(C<GroupOperation*>, std::stringstream& ss)
{
    auto pos=ss.tellg();
    std::string line=ss.str().substr(pos);
    auto g=std::unique_ptr<GroupOperation>(new GroupOperation);
    auto res=g->get(ss);
    if (res) return g.release();
    else return {false,res.error()};
}

inline myOptional_t<ArrayOperation*> get(C<ArrayOperation*>, std::stringstream& ss)
{
    std::string line=ss.str();
    auto g=std::unique_ptr<ArrayOperation>(new ArrayOperation);
    auto res=g->get(ss);
    if (res) return g.release();
    else return {false,res.error()};
}


template<bool B,class AssignmentGenericOperator>
class basic_Assignment: public Generic_Assignment
{
private:
    std::unique_ptr<Identifier> id_;
    std::unique_ptr<AssignmentGenericOperator> op_;
    std::unique_ptr<Expression> exp_;
public:
    typedef Generic_Assignment base_type;

    constexpr static bool delayed=B;
    virtual ~basic_Assignment()override{}
    basic_Assignment()=default;
    basic_Assignment(Identifier* id, AssignmentGenericOperator* op, Expression* exp):
        id_{id},op_{op},exp_{exp}{}
    Identifier const* id()const override {return id_.get();}
    AssignmentGenericOperator const * op()const override {return op_.get();}
    Expression const * expr()const override {return exp_.get();}


    virtual std::ostream &put(std::ostream &os) const override;
    virtual myOptional_t<void> get(std::stringstream &ss) override;

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
    virtual ~LiteralGeneric()override=default;
    virtual LiteralGeneric* clone()const override=0;

    virtual std::string myType()const =0;
    typedef Term base_type;

};

template<typename T>
class Literal: public LiteralGeneric
{
    T value_;
public:
    //Literal(T&& value):value_{std::move(value)}{}

    typedef LiteralGeneric base_type;

    virtual std::string myType()const override
    {
        return my_trait<T>::className.str();
    }
    Literal(const T& value):value_{value}{}
    Literal(){}
    virtual ~Literal()override{}

    T getValue()const {return value_;}
    virtual std::string value() const override
    {
        return ToString(value_);
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
    virtual myOptional_t<void> get(std::stringstream &ss) override
    {
        auto pos=ss.tellg();
        T value;
        if (ss>>value)
        {
            value_=std::move(value);
            return {true,""};
        }
        else
        {
            ss.clear();
            ss.seekg(pos);
            return {false,print_pos(pos)+ "not a"+my_trait<T>::className.str()};
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
    virtual std::string myType()const override
    {
        return my_trait<std::string>::className.str();
    }
    typedef LiteralGeneric base_type;

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
    virtual myOptional_t<void> get(std::stringstream &ss) override
    {
        auto pos=ss.tellg();
        std::string line=ss.str().substr(pos);
        if (!(ss>>label_start{}))
        {
            ss.clear();
            ss.seekg(pos);
            return {false, print_pos(pos)+"lacks label start"};

        }
        else
        {
            std::string out;
            std::string end=label_end{}.str();
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
                return { true,""};
            }
            else
            {
                ss.clear();
                ss.seekg(pos);
                return {false, print_pos(pos)+"some eror in label that I do not currently understand"};

            }
        }
    }
    virtual Literal* clone()const override { return new Literal(*this);}

    Literal(const Literal&)=default;
    Literal(Literal&&)=default;
    Literal& operator=(const Literal&)=default;
    Literal& operator=(Literal&&)=default;
    ~Literal()override=default;

};









inline std::string UnaryOperation::value() const
{
    return op()->name();

}

inline std::ostream &UnaryOperation::put(std::ostream &os) const
{
    op()->put(os);
    arg(0)->put(os);
    return os;
}

inline std::ostream & Function::put(std::ostream &os) const
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


inline myOptional_t<void> Function::get(std::stringstream &ss)
{
    auto pos=ss.tellg();
    Identifier id;

    auto oid=id.get(ss);
    if (!oid)
    {ss.clear(); ss.seekg(pos); return {false,print_pos(pos)+"Lacks identifier: "+oid.error()};}
    else if (!(ss>>arg_start{}))
    {ss.clear(); ss.seekg(pos); return {false,print_pos(pos)+"Lacks argument starter: "};}
    else{
        auto oarg=getArguments(C<typename Function::arg_map>{},ss);
        if (!oarg)
        {ss.clear(); ss.seekg(pos); return {false,print_pos(pos)+"error in function arguments: "+oarg.error()};}
        else
        {
            *this=Function(std::move(id),oarg.release());
            return  {true,""};
        }
    }
}


inline myOptional_t<typename Function::arg_map> Function::getArguments(C<typename Function::arg_map>,std::stringstream& ss){
    typename Function::arg_map out;
    typedef std::variant<arg_end,Generic_Assignment*> init_type;
    auto pos=ss.tellg();
    std::string line=ss.str().substr(pos);
    while (true)
    {
        auto v=get_variant(C<init_type>{},ss);
        if (!v.has_value())
        {
            return {false,"unexpected end of function arguments: "+v.error()};
        }
        else if (v.value().index()==1)
        {
            auto a=std::get<Generic_Assignment*>(v.value());
            out.emplace(a->id()->value(),a);
        }
        else return myOptional_t<typename Function::arg_map>(std::move(out));

    }

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
myOptional_t<void> basic_Assignment<B,AssignmentGenericOperator>::get(std::stringstream &ss)
{
    auto pos=ss.tellg();
    auto oid=grammar::get(C<Identifier*>{},ss);
    if (!oid) {ss.clear(); ss.seekg(pos); return {false,print_pos(pos)+"Assignment without Identifier: "+oid.error()};}
    else
    {
        auto oa=grammar::get(C<AssignmentGenericOperator*>{},ss);
        if (!oa ) {ss.clear(); ss.seekg(pos); return {false,print_pos(pos)+"Assignment without assignment operator: "+oa.error()};}
        else {
            auto oe=grammar::get(C<Expression*>{},ss);
            if (!oe ) {ss.clear(); ss.seekg(pos); return {false,print_pos(pos)+"Assignment with an error in Expression: "+oe.error()};}
            else {
                *this=basic_Assignment(oid.release(),oa.release(),oe.release());
                return {true,""};
            }
        }
    }

}



template<bool B,class AssignmentGenericOperator>
std::string basic_Assignment<B,AssignmentGenericOperator>
::value() const
{
    return op()->name();
}

typedef std::variant<LiteralGeneric*,Identifier*,GroupOperation*,ArrayOperation*,UnaryOperator*> valid_start;
typedef  std::variant< BinaryOperator*, AssignmentGenericOperator*, DefinitionGenericOperator*> after_idenfier;




inline
myOptional_t<BinaryOperation*> resolve_operator_order(C<BinaryOperation*>,std::stringstream& ss, std::unique_ptr<Expression>&& e0, std::unique_ptr<BinaryOperator>&& op0, std::unique_ptr<Expression>&& e1)
{
    auto oop1=get(C<BinaryOperator*>{},ss);
    if (oop1.has_value())
    {
        auto op1=std::unique_ptr<BinaryOperator>(oop1.release());
        auto ee2=get(C<Term*>{},ss);
        if (!ee2.has_value())
            return {false,"second term missing :"+ee2.error()};
        auto e2=std::unique_ptr<Expression>(ee2.release());
        if (op0->order()<=op1->order())
        {
            auto e0_op0_e1=std::unique_ptr<Expression>(new BinaryOperation(e0.release(),op0.release(),e1.release()));
            return resolve_operator_order(C<BinaryOperation*>{},ss,std::move(e0_op0_e1),std::move(op1),std::move(e2));
        }
        else
        {
            auto e1_op1_e2=resolve_operator_order(C<BinaryOperation*>{},ss,std::move(e1),std::move(op1),std::move(e2));
            if (e1_op1_e2.has_value())
            {
                return new BinaryOperation(e0.release(),op0.release(),e1_op1_e2.release());
            }
            else
                return {false,e1_op1_e2.error()};
        }
    }
    return new BinaryOperation(e0.release(),op0.release(),e1.release());
}


inline
myOptional_t<Term*> get_term_from_identifier(C<Term*>,std::stringstream& ss, std::unique_ptr<Identifier>&& id)
{
    auto pos=ss.tellg();
    std::string line=ss.str().substr(pos);
    if (get(ss,Function::start_symbol{}))
    {
        auto out=Function::getArguments(C<typename Function::arg_map>{},ss);
        if (!out){ return {false,"error on function arguments: "+out.error()};}
        else {
          return new Function(std::move(*id.release()), out.release());
        }
    }
    else
    {
        return id.release();
    }

}





inline
myOptional_t<Term*> get_term_from_the_rest(C<Term*>,std::stringstream& ss,  valid_start& first)
{
    if (std::holds_alternative<UnaryOperator*>(first))
    {
        auto unary=std::get<UnaryOperator*>(first);
        auto ot=grammar::get(C<Term*>{},ss);
        if (!ot){return {false,"has no term after UnaryOperator: "+ot.error()};}
        else
        {
            return new UnaryOperation(unary,ot.release());
        }
    }
    else if (std::holds_alternative<LiteralGeneric*>(first))
    {return std::get<LiteralGeneric*>(first); }
    else if (std::holds_alternative<GroupOperation*>(first))
    { return std::get<GroupOperation*>(first); }
    else if (std::holds_alternative<ArrayOperation*>(first))
    { return std::get<ArrayOperation*>(first);}
    else
        // something very wrong
    {
        return {false, "something fishy in get_term_from_the_rest"};
    }


}



inline
myOptional_t<Term*> get(C<Term*>,std::stringstream& ss)
{
    auto pos=ss.tellg();
    std::string line=ss.str().substr(pos);
    auto ofirst=get_variant(C<valid_start>{},ss);
    if (!ofirst) { ss.clear(); ss.seekg(pos); return {false,print_pos(pos)+ "invalid start of term:"+ofirst.error()}; }
    else
    {
        if (std::holds_alternative<Identifier*>(ofirst.value()))
        {
            auto id=std::get<Identifier*>(ofirst.value());
            auto ot=grammar::get_term_from_identifier(C<Term*>{},ss,std::unique_ptr<Identifier>(id));
            if (!ot) {ss.clear(); ss.seekg(pos); return {false,print_pos(pos)+"error after Id: "+ot.error()};}
            else return ot;
        }
        else return get_term_from_the_rest(C<Term*>{},ss, ofirst.value());

    }
}




inline
myOptional_t<Expression*> get_expression_from_term(C<Expression*>,std::stringstream& ss, Term* t)
{
    auto oop=grammar::get(C<BinaryOperator*>{},ss);
    if (!oop)
    {
        return t;
    }
    else
    {
        auto ot2=grammar::get(C<Term*>{},ss);
        if (!ot2) return {false,"term expected after BinaryOperator, found error: "+ot2.error()};
        else
        {
            return resolve_operator_order(C<BinaryOperation*>{},ss,std::unique_ptr<Expression>(t),
                                          std::unique_ptr<BinaryOperator>(oop.release()),std::unique_ptr<Expression>(ot2.release()));
        }
    }
}




inline
myOptional_t<Expression*> get(C<Expression*>,std::stringstream& ss)
{
    auto pos=ss.tellg();
    std::string line=ss.str().substr(pos);
    auto ot=get(C<Term*>{},ss);
    if (!ot)
    { ss.clear(); ss.seekg(pos); return {false, print_pos(pos)+"error in Term :"+ot.error()}; }
    return get_expression_from_term(C<Expression*>{},ss,ot.release());
}







inline myOptional_t<Statement*> get(C<Statement*>,std::stringstream& ss)
{
    auto pos=ss.tellg();
    std::string line=ss.str().substr(pos);

    auto ofirst=get_variant(C<valid_start>{},ss);
    if (!ofirst)
    { ss.clear(); ss.seekg(pos); return {false,print_pos(pos)+"Statement failed start: "+ofirst.error()}; }
    else
    {
        if (std::holds_alternative<Identifier*>(ofirst.value()))
        {
            auto id=std::get<Identifier*>(ofirst.value());
            auto oopa=get(C<AssignmentGenericOperator*>{},ss);
            if (oopa)
            {
                auto oe2=grammar::get(C<Expression*>{},ss);
                if (!oe2)
                {ss.clear(); ss.seekg(pos); return {false,print_pos(pos)+"failed expression after asignment :"+ oe2.error()};}
                else {
                    return new Assignment(id, oopa.release(),oe2.release());
                }
            }
            else
            {
                auto oopd=get(C<DefinitionGenericOperator*>{},ss);
                if (oopd)
                {
                    auto oe2=grammar::get(C<Expression*>{},ss);
                    if (!oe2)
                    {ss.clear(); ss.seekg(pos); return {false,print_pos(pos)+"failed expression after definition :"+ oe2.error()};}
                    else {
                        return new Definition(id, oopd.release(),oe2.release());
                    }
                }
                else
                {
                    auto ot=get_term_from_identifier(C<Term*>{},ss,std::unique_ptr<Identifier>(id));
                    if (!ot) {ss.clear(); ss.seekg(pos); return {false, print_pos(pos)+"error after identifier : "+ot.error()};}
                    else return get_expression_from_term(C<Expression*>{},ss,ot.release());

                }
            }
        }
        else
        {
            auto ot=get_term_from_the_rest(C<Term*>{},ss,ofirst.value());
            if (!ot)
            {ss.clear(); ss.seekg(pos); return {false, print_pos(pos)+"error after rest: "+ot.error()};}
            else return get_expression_from_term(C<Expression*>{},ss,ot.release());

        }

    }

}



inline
myOptional_t<Generic_Assignment*> get(C<Generic_Assignment*>,std::stringstream& ss)
{
    auto pos=ss.tellg();
    std::string line=ss.str().substr(pos);
    auto oid=get(C<Identifier*>{},ss);
    if (!oid) {ss.clear(); ss.seekg(pos);return {false,print_pos(pos)+"assignment with invalid identifier: "+oid.error()};}
    auto oop=get_variant(C<std::variant<AssignmentGenericOperator*,DefinitionGenericOperator*>>{},ss);
    if (!oop)
    {ss.clear(); ss.seekg(pos);return {false,print_pos(pos)+"assignment with invalid operator : "+oop.error()};}
    auto oe=grammar::get(C<Expression*>{},ss);
    if (!oe)
    {ss.clear(); ss.seekg(pos);return {false,print_pos(pos)+"assignment with invalid expression : "+oe.error()};}
    else
    {
        if (std::holds_alternative<AssignmentGenericOperator*>(oop.value()))
        {
            auto opa=std::get<AssignmentGenericOperator*>(oop.value());
            return new Assignment(oid.release(),opa,oe.release());
        }
        else if (std::holds_alternative<DefinitionGenericOperator*>(oop.value()))
        {
            auto opd=std::get<DefinitionGenericOperator*>(oop.value());
            return new Definition(oid.release(),opd,oe.release());
        }
        else return {false,"not an Assignment nor a Definition, error"};

    }
}




inline
myOptional_t<Identifier*> get(C<Identifier*>,std::stringstream& ss)
{
    auto pos=ss.tellg();
    std::string line=ss.str().substr(pos);
    return get_Derived(Cs<Identifier, Identifier>{},ss,"");
}

inline
myOptional_t<AssignmentGenericOperator*> get(C<AssignmentGenericOperator*>,std::stringstream& ss)
{
    auto pos=ss.tellg();
    std::string line=ss.str().substr(pos);
    return get_Derived(Cs<AssignmentGenericOperator,AssignmentOperator>{},ss,"");
}





inline
myOptional_t<DefinitionGenericOperator*> get(C<DefinitionGenericOperator*>,std::stringstream& ss)
{
    auto pos=ss.tellg();
    std::string line=ss.str().substr(pos);
    return get_Derived(Cs<DefinitionGenericOperator>{},ss,"");
}
inline
myOptional_t<UnaryOperator*> get(C<UnaryOperator*>,std::stringstream& ss)
{
    auto pos=ss.tellg();
    std::string line=ss.str().substr(pos);
    return get_Derived(Cs<UnaryOperator>{},ss,"");
}
inline
myOptional_t<BinaryOperator*> get(C<BinaryOperator*>,std::stringstream& ss)
{
    auto pos=ss.tellg();
    std::string line=ss.str().substr(pos);
    return get_Derived(Cs<BinaryOperator>{},ss,"");
}


inline myOptional_t<LiteralGeneric*> get(C<LiteralGeneric*>,std::stringstream& ss)
{
    auto pos=ss.tellg();
    std::string line=ss.str().substr(pos);
    return get_Derived(Cs<LiteralGeneric,Literal<std::string>, Literal<double>, Literal<std::size_t>,Literal<std::mt19937_64::result_type>, Literal<int>, Literal<bool>>{},ss,"");
}





inline std::string BinaryOperation::value() const
{
    return op()->name();
}
inline std::ostream &BinaryOperation::put(std::ostream &os) const
{
    arg(0)->put(os);
    op()->put(os);
    arg(1)->put(os);
    return os;

}

inline myOptional_t<void> BinaryOperation::get(std::stringstream &ss)
{
    auto pos=ss.tellg();
    auto ot1=grammar::get(C<Term*>{},ss);
    if (!ot1)
    {ss.clear();ss.seekg(pos); return {false, print_pos(pos)+"binary operation failed on first term: "+ot1.error()};}
    else
    {
        auto oop=grammar::get(C<BinaryOperator*>{},ss);
        if (!oop)
        {ss.clear();ss.seekg(pos); return {false,print_pos(pos)+ "binary operation failed on operator: "+oop.error()};}
        else
        {
            auto ot2=grammar::get(C<Term*>{},ss);
            if (!ot2)
            {ss.clear();ss.seekg(pos); return {false,print_pos(pos)+ "binary operation failed on second term: "+ot2.error()};}
            else
            {
                auto ob= resolve_operator_order(C<BinaryOperation*>{},ss,std::unique_ptr<Term>(ot1.release()),std::unique_ptr<BinaryOperator>(oop.release())
                                                ,std::unique_ptr<Term>(ot2.release()));
                if (ob) {*this=std::move(*ob.release());
                    return {true,""};}
                else return {false,ob.error()};
            }
        }
    }
}


inline myOptional_t<void> UnaryOperation::get(std::stringstream &ss)
{
    auto pos=ss.tellg();
    auto oop= grammar::get(C<UnaryOperator*>{},ss);
    if (!oop)
    {
        ss.clear();
        ss.seekg(pos);
        return {false,print_pos(pos)+"Unary Operator not recognized:"+oop.error()};
    }
    else
    {
        auto oe=grammar::get(C<Term*>{},ss);
        if (!oe)
        {
            ss.clear(); ss.seekg(pos);return {false,print_pos(pos)+"Term error after Unary Operator:"+oe.error()};

        }    else
        {
            *this=UnaryOperation(oop.release(),oe.release());
            return {true,""};
        }
    }
}

namespace get_variant_templ_impl {





template<class...Ts>
myOptional_t<std::variant<Ts...>> get(Cs<Ts...>,Cs<>,std::stringstream& , std::string error)
{
    return {false,error};
}

template<class... T0s, class T,class...Ts>
myOptional_t<std::variant<T0s...,T,Ts...>> get(Cs<T0s...>,Cs<T,Ts...>,std::stringstream& ss, std::string error)
{
    if constexpr (has_get_global_optional<T>::value)
    {
        auto o= grammar::get(C<T>{},ss);
        if (o.has_value())
        {
            return std::variant<T0s...,T,Ts...>(o.release());
        }
        else error+=o.error()+" ,";
    }
    else if constexpr (has_get_global<T>::value)
    {
        T e;

        if (get(ss,e))
        {
            return std::variant<T0s...,T,Ts...>(e);
        }
        else error+" not valid input ,";
    }

    else if constexpr(std::is_pointer_v<T>)
    {
        std::unique_ptr<std::remove_pointer_t<T>> e(new  std::remove_pointer_t<T>{});
        auto o= e->get(ss);
        if (o.has_value())
        {
            return std::variant<T0s...,T,Ts...>(e.release());
        }
        else error+=o.error()+" ,";

    }
    else
    {
        T e;
        auto o=e.get(ss);
        if (o)
        {
            return std::move(e);
        }
        else error+=o.error()+" ,";
    }
    return get(Cs<T0s...,T>{},Cs<Ts...>{},ss, error);
}


} // namespace get_variant_impl


} // namespace grammar



#endif // MYGRAMMAR_H

