#ifndef MYPARAMETERS_H
#define MYPARAMETERS_H

#include "mytypetraits.h"
#include "mygrammar.h"
#include "myDistributions.h"
#include <vector>
#include <map>
#include <set>
#include <cmath>


template< typename T>
class Label
{
    std::string name_;
public :
     constexpr static const char legal_chars[]="abcdefghijklmopqrstuvwxz_0123456789ABCDEFGHIJKLMOPQRSTUVWXYZ";

      typedef T referred_type;
    std::string name()const {return name_;} // unique at all levels

    Label(std::string&& label): name_{std::move(label)}{}
    Label(const std::string& label): name_{label}{}

    Label()=default;

    operator std::string()const {return name_;}

    bool operator < (const Label& other)const
    {
        return name()<other.name();
    }
    bool operator==(const Label& other){return name()==other.name();}

    std::ostream& write(std::ostream& os)const { return os<<name_;}
    std::istream& read(std::istream& is) { return is>>name_;}

};





template<bool output,typename X>
struct are_Equal<output,X,std::enable_if_t<is_label<X>::value, void>>
{

public:
  are_Equal(double, double)  {}
  are_Equal(){}
  bool test(const X& x, const X& y, std::ostream &os ) const
  {
    if (x.name()!=y.name())
    {
      os<<"Differ in Label :"<<x.name()<<" vs "<<y.name();
      return false;
    }
    else return true;
  }


};



template< typename T>
bool operator==(const Label<T>& one,const Label<T>& other){return one.name()==other.name();}


template<class T>
std::ostream& operator<<(std::ostream& os, const Label<T>& x)
{
    return os<<x.name();
}
template<class T>
std::istream& operator>>(std::istream& is,  Label<T>& x)
{
    return x.read(is);
}


template <class label>
struct LabelMap
{
    //typedef typename label::referred_type referred_type;
    typedef   std::map<std::string,double > type;

};

template <class label>
using LabelMap_t=typename LabelMap<label>::type;


class Parameter_label
        : public Label<double> {public: using Label::Label;  using Label::legal_chars;
    using Label::referred_type;
};

template<typename Model, typename aParameter_label >
class Parameters_labels
{
    std::vector<aParameter_label const&> v_;
public:

    Parameters_labels()=default;

    Parameters_labels(std::vector<aParameter_label const&> labels):v_{std::move(labels)}{}



};






template<typename Model>
class Parameters_values
{
public:
    typedef typename Model::myParameter_label par_label;

  private:
    //LabelMap_t<par_label> d_;
      std::map<std::string,double> d_;

  public:

  //  constexpr static const char legal_chars[]=par_label::legal_chars;
    std::size_t size()const { return d_.size();}
    Parameters_values(std::map<std::string,double> values): d_{values}{};
    Parameters_values()=default;
    double& at(const par_label& label)
    {
        return d_.at(label);
    }

    double const& at(const par_label& label)const
    {
        return d_.at(label);
    }

    auto& getParameterMap()const { return d_;}
    typedef  Parameters_values<Model> self_type;
    constexpr static auto  className=my_trait<Model>::className+my_static_string("_Parameters");
    static auto get_constructor_fields()
    {
        return std::make_tuple(
                    grammar::field(C<self_type>{},"values",&self_type::getParameterMap));
    }

    Parameters_values& operator<<(Parameters_values&& other)
    {
        d_.merge(std::move(other.getParameterMap()));
        return *this;
    }
    Parameters_values& operator<<(const Parameters_values& other)
    {
        for (auto e:other.getParameterMap())
            d_[e.first]=e.second;
        return *this;
    }

};

template<template<class...> class Tr,typename Model>
class Parameters_values_new_;
template<typename Model>
using Parameters_values_new=Parameters_values_new_<C,Model>;


template<template<class...> class Tr,typename Model>
class Parameters_values_new_
{
public:

    template<class T> using Tr_t=typename Tr<T>::type;
    typedef typename Model::myParameter_label par_label;

private:
   // Tr_t<LabelMap_t<par_label>> d_;
    Tr_t<std::map<std::string,double>> d_;

public:

    //  constexpr static const char legal_chars[]=par_label::legal_chars;
    std::size_t size()const { return d_.size();}
    Parameters_values_new_(Tr_t<std::map<std::string,double>> values): d_{values}{}
    Parameters_values_new_()=default;
    Tr_t<double>& at(const par_label& label)
    {
        return d_.at(label);
    }

    Tr_t<double> const& at(const par_label& label)const
    {
        return d_.at(label);
    }

    auto& getParameterMap()const { return d_;}
    typedef  Parameters_values_new_ self_type;
    constexpr static auto  className=my_template_trait<Tr>::className+my_trait<Model>::className+my_static_string("_Parameters_new");
    static auto get_constructor_fields()
    {
        return std::make_tuple(
            grammar::field(C<self_type>{},"values",&self_type::getParameterMap));
    }

    Parameters_values_new_& operator<<(Parameters_values_new_&& other)
    {
        d_.merge(std::move(other.getParameterMap()));
        return *this;
    }
    Parameters_values_new_& operator<<(const Parameters_values_new_& other)
    {
        for (auto e:other.getParameterMap())
            d_[e.first]=e.second;
        return *this;
    }

};





class Range
{
public:
  constexpr static auto const className = my_static_string("Range");
  typedef Range self_type;
  static auto get_constructor_fields() {
    return std::make_tuple(
        grammar::field(C<self_type>{}, "min", &self_type::min),
        grammar::field(C<self_type>{}, "max", &self_type::max));
  }

  Range (double minValue, double  maxValue ): min_{minValue}, max_{maxValue}{}

  Range()=default;

  bool operator()(double x)const {return x>min_&&x<max_;}

 double min()const {return min_;}
  double max()const {return max_;}

private:
  double min_=-std::numeric_limits<double>::infinity();
  double max_=+std::numeric_limits<double>::infinity();
};




template<typename Model>
class Parameters_distribution: public Base_Distribution<M_Matrix<double>>
{
public:
    typedef Base_Distribution<M_Matrix<double>> base_type;
    virtual Parameters_distribution* clone()const override{ return new Parameters_distribution(*this);}
    std::string myClass()const override { return className.str();}

    typedef typename Model::myParameter_label par_label;
    typedef std::tuple<std::string,std::unique_ptr<Base_Transformation<double>>, std::unique_ptr<Base_Distribution<double>>,Range> tuple;
    typedef std::tuple<std::string,Base_Transformation<double>*, Base_Distribution<double>*,Range> tuple_ptr;

    // Base_Distribution interface
public:
    virtual M_Matrix<double> sample(std::mt19937_64 &mt) const override
    {
        M_Matrix<double> out(1,size());
        for (std::size_t i=0; i<size(); ++i)
        {
            out[i]=dist(i)->sample(mt);
        }
        return out;
    }
    virtual double p(const M_Matrix<double> &x) const override
    {
        return std::exp(logP(x));
    }
    virtual double logP(const M_Matrix<double> &x) const override
    {
        double sum=0;
        for (std::size_t i=0; i<size(); ++i)
        {
            sum+=dist(i)->logP(x[i]);
        }
        return sum;
    }

    virtual double vlogP(const M_Matrix<double> &x)const override {
      double sum = 0;
      for (std::size_t i = 0; i < size(); ++i) {
        sum += dist(i)->vlogP(x[i]);
      }
      return sum;
    }


    virtual double expected_logP()const override
    {
        double sum=0;
        for (std::size_t i=0; i<size(); ++i)
        {
            sum+=dist(i)->expected_logP();
        }
        return sum;
    }

    virtual double variance_logP()const override
    {
        double sum=0;
        for (std::size_t i=0; i<size(); ++i)
        {
            sum+=dist(i)->variance_logP();
        }
        return sum;
    }


    virtual  M_Matrix<M_Matrix<double>> param() const override
    {
        M_Matrix<M_Matrix<double>> out(1,size());
        for (std::size_t i=0; i<size(); ++i)
        {
            out[i]=dist(i)->param();
        }
        return out;
    }
    virtual M_Matrix<M_Matrix<double>> Fisher_Information() const override
    {
        M_Matrix<M_Matrix<double>> out(1,size());
        for (std::size_t i=0; i<size(); ++i)
        {
            out[i]=dist(i)->Fisher_Information();
        }
        return out;

    }
    virtual M_Matrix<double> dlogL_dx(const M_Matrix<double> &x) const override
    {
        M_Matrix<double> out(1,size());
        for (std::size_t i=0; i<size(); ++i)
        {
            out[i]=dist(i)->dlogL_dx(x[i]);
        }
        return out;
    }
    virtual M_Matrix<double> dlogL_dx2(const M_Matrix<double> &x) const override
    {
        M_Matrix<double> out(size(),size(),Matrix_TYPE::DIAGONAL);
        for (std::size_t i=0; i<size(); ++i)
        {
            out(i,i)=dist(i)->dlogL_dx2(x[i]);
        }
        return out;

    }
    virtual M_Matrix<double> mean() const override
    {
        M_Matrix<double> out(1,size());
        for (std::size_t i=0; i<size(); ++i)
        {
            out[i]=dist(i)->mean();
        }
        return out;
    }
    virtual M_Matrix<double> stddev() const override
    {
        M_Matrix<double> out(1,size());
        for (std::size_t i=0; i<size(); ++i)
        {
            out[i]=dist(i)->stddev();
        }
        return out;

    }




    virtual myOptional_t<Parameters_values<Model>> tr_to_Parameter(const M_Matrix<double>& val)const
    {
      typedef myOptional_t<Parameters_values<Model>> Op;
        assert(val.size()==tu_.size());
        LabelMap_t<par_label> m;
        std::vector<Op_void> res;
        for ( std::size_t i=0; i<val.size(); ++i)
        {
          res.emplace_back(is_in_range(i,val[i]));
          if (res[i])
          {
            m[name(i)]=tr(i)->apply_inv(val[i]);
          }
        }
        auto r=consolidate(std::move(res));
        if (r)
           return Op(Parameters_values<Model>(m));
        else return Op(false,r.error());
    }
    virtual double tr_to_Parameter(double val, std::size_t ipar)const
    {
        assert(ipar<size());
        return tr(ipar)->apply_inv(val);
    }




    M_Matrix<double> Parameter_to_tr(const Parameters_values<Model> & val) const
    {
        assert(val.size()==tu_.size());
        M_Matrix<double> out(1,val.size());
        for ( std::size_t i=0; i<val.size(); ++i)
        {
            out[i]=tr(i)->apply(val.at(name(i)));
        }
        return out;
    }


    
    std::string name(std::size_t i)const { return  std::get<0>(tu_[i]);}

    Base_Transformation<double> * tr(std::size_t i)const{return std::get<1>(tu_[i]).get();}

    Base_Distribution<double> * dist(std::size_t i)const {return std::get<2>(tu_[i]).get();}

    Range const & range(std::size_t i)const {return std::get<3>(tu_[i]);}

    Op_void is_in_range(std::size_t i, double x)const { if (range(i)(x)) return Op_void(true,"");
      else return Op_void(false,name(i)+"("+ToString(i)+"th parameter)="+ToString(x)+ "is outside "+ToString(range(i)));
    }

    std::size_t size()const { return tu_.size();}
    
    auto getParameterDistributionMap()const {
      std::vector<tuple_ptr> out(size());
        for (std::size_t i=0; i<size(); ++i)
          out[i]={name(i),tr(i),dist(i),range(i)};
        return out;}
    typedef  Parameters_distribution<Model> self_type;
    constexpr static auto  className=my_trait<Model>::className+my_static_string("_Parameters_Distribution");
    static auto get_constructor_fields()
    {
        return std::make_tuple(
                    grammar::field(C<self_type>{},"values",&self_type::getParameterDistributionMap));
    }
    virtual std::ostream& put(std::ostream& os)const override
    {
      return io::output_operator_on_Object (os,*this);
    }

    Parameters_distribution(const std::vector<tuple_ptr>& in)
        :tu_{tu2tp(in)}{}

    Parameters_distribution()=default;

    Parameters_distribution(const Parameters_distribution& other): tu_{tucopy(other.tu_)}{};
    Parameters_distribution( Parameters_distribution&&)=default;
    Parameters_distribution& operator=(const Parameters_distribution& other)
    {if (this!=&other){
            Parameters_distribution tmp(other); *this=std::move(tmp);} return *this;
    }
    Parameters_distribution& operator=( Parameters_distribution&&)=default;



  protected:

    std::vector<tuple> tu_;


    static
    std::vector<tuple>
        tu2tp(const std::vector<tuple_ptr>& in)
    {
        std::vector<tuple> out;
        for (auto& e:in)
          out.emplace_back(std::get<0>(e),std::get<1>(e), std::get<2>(e),std::get<3>(e));
        return out;

    }

    std::vector<tuple>
    tucopy(const std::vector<tuple>& in
           )
    {
        std::vector<tuple> out;
        for (auto& e:in)
            out.emplace_back(std::get<0>(e),
                             std::unique_ptr<Base_Transformation<double>>(std::get<1>(e)->clone()),
                           std::unique_ptr<Base_Distribution<double>>(std::get<2>(e)->clone()),
                           std::get<3>(e));
        return out;

    }

};


template<template<class...> class Tr,typename Model>
class Parameters_distribution_new_;
template<typename Model>
using Parameters_distribution_new=Parameters_distribution_new_<C,Model>;

//template <class... > class Der;
template<template<class...> class Tr,typename Model>
myOptional_t<typename Tr<Parameters_values_new<Model>>::type>
tr_to_Parameter_impl(const Parameters_distribution_new_<Tr,Model>* self,const M_Matrix<double>& val);

template<template<class...> class Tr,typename Model>
class Parameters_distribution_new_: public Base_Distribution<M_Matrix<double>>
{
public:
    template<class T> using Tr_t=typename Tr<T>::type;
    typedef Base_Distribution<M_Matrix<double>> base_type;
    virtual Parameters_distribution_new_* clone()const override{ return new Parameters_distribution_new_(*this);}
    std::string myClass()const override { return className.str();}

    typedef typename Model::myParameter_label par_label;
    typedef std::tuple<std::string,std::unique_ptr<Base_Transformation<double>>, std::unique_ptr<Base_Distribution<double>>,Range> tuple;
    typedef std::tuple<std::string,Base_Transformation<double>*, Base_Distribution<double>*,Range> tuple_ptr;

    // Base_Distribution interface
public:
    virtual M_Matrix<double> sample(std::mt19937_64 &mt) const override
    {
        M_Matrix<double> out(1,size());
        for (std::size_t i=0; i<size(); ++i)
        {
            out[i]=dist(i)->sample(mt);
        }
        return out;
    }
    virtual double p(const M_Matrix<double> &x) const override
    {
        return std::exp(logP(x));
    }
    virtual double logP(const M_Matrix<double> &x) const override
    {
        double sum=0;
        for (std::size_t i=0; i<size(); ++i)
        {
            sum+=dist(i)->logP(x[i]);
        }
        return sum;
    }

    virtual double vlogP(const M_Matrix<double> &x)const override {
        double sum = 0;
        for (std::size_t i = 0; i < size(); ++i) {
            sum += dist(i)->vlogP(x[i]);
        }
        return sum;
    }


    virtual double expected_logP()const override
    {
        double sum=0;
        for (std::size_t i=0; i<size(); ++i)
        {
            sum+=dist(i)->expected_logP();
        }
        return sum;
    }

    virtual double variance_logP()const override
    {
        double sum=0;
        for (std::size_t i=0; i<size(); ++i)
        {
            sum+=dist(i)->variance_logP();
        }
        return sum;
    }


    virtual  M_Matrix<M_Matrix<double>> param() const override
    {
        M_Matrix<M_Matrix<double>> out(1,size());
        for (std::size_t i=0; i<size(); ++i)
        {
            out[i]=dist(i)->param();
        }
        return out;
    }
    virtual M_Matrix<M_Matrix<double>> Fisher_Information() const override
    {
        M_Matrix<M_Matrix<double>> out(1,size());
        for (std::size_t i=0; i<size(); ++i)
        {
            out[i]=dist(i)->Fisher_Information();
        }
        return out;

    }
    virtual M_Matrix<double> dlogL_dx(const M_Matrix<double> &x) const override
    {
        M_Matrix<double> out(1,size());
        for (std::size_t i=0; i<size(); ++i)
        {
            out[i]=dist(i)->dlogL_dx(x[i]);
        }
        return out;
    }
    virtual M_Matrix<double> dlogL_dx2(const M_Matrix<double> &x) const override
    {
        M_Matrix<double> out(size(),size(),Matrix_TYPE::DIAGONAL);
        for (std::size_t i=0; i<size(); ++i)
        {
            out(i,i)=dist(i)->dlogL_dx2(x[i]);
        }
        return out;

    }
    virtual M_Matrix<double> mean() const override
    {
        M_Matrix<double> out(1,size());
        for (std::size_t i=0; i<size(); ++i)
        {
            out[i]=dist(i)->mean();
        }
        return out;
    }
    virtual M_Matrix<double> stddev() const override
    {
        M_Matrix<double> out(1,size());
        for (std::size_t i=0; i<size(); ++i)
        {
            out[i]=dist(i)->stddev();
        }
        return out;

    }




    auto tr_to_Parameter(const M_Matrix<double>& val)const
    {
        return tr_to_Parameter_impl(this,val);
    }


    virtual double tr_to_Parameter(double val, std::size_t ipar)const
    {
        assert(ipar<size());
        return tr(ipar)->apply_inv(val);
    }




    M_Matrix<double> Parameter_to_tr(const Parameters_values_new<Model> & val) const
    {
        assert(val.size()==tu_.size());
        M_Matrix<double> out(1,val.size());
        for ( std::size_t i=0; i<val.size(); ++i)
        {
            out[i]=tr(i)->apply(val.at(name(i)));
        }
        return out;
    }





    std::string name(std::size_t i)const { return  std::get<0>(tu_[i]);}

    Base_Transformation<double> * tr(std::size_t i)const{return std::get<1>(tu_[i]).get();}

    Base_Distribution<double> * dist(std::size_t i)const {return std::get<2>(tu_[i]).get();}

    Range const & range(std::size_t i)const {return std::get<3>(tu_[i]);}

    Op_void is_in_range(std::size_t i, double x)const { if (range(i)(x)) return Op_void(true,"");
        else return Op_void(false,name(i)+"("+ToString(i)+"th parameter)="+ToString(x)+ "is outside "+ToString(range(i)));
    }

    std::size_t size()const { return tu_.size();}

    auto getParameterDistributionMap()const {
        std::vector<tuple_ptr> out(size());
        for (std::size_t i=0; i<size(); ++i)
            out[i]={name(i),tr(i),dist(i),range(i)};
        return out;}
    typedef  Parameters_distribution_new_ self_type;
    constexpr static auto  className=my_trait<Model>::className+my_static_string("_Parameters_Distribution_new");
    static auto get_constructor_fields()
    {
        return std::make_tuple(
            grammar::field(C<self_type>{},"values",&self_type::getParameterDistributionMap));
    }
    virtual std::ostream& put(std::ostream& os)const override
    {
        return io::output_operator_on_Object (os,*this);
    }

    Parameters_distribution_new_(const std::vector<tuple_ptr>& in)
        :tu_{tu2tp(in)}{}

    Parameters_distribution_new_()=default;

    template<template<class...> class Tr_other>
    Parameters_distribution_new_(const Parameters_distribution_new_<Tr_other,Model>& other): tu_{tucopy(other.get_data())}{};

    Parameters_distribution_new_(const Parameters_distribution_new_& other): tu_{tucopy(other.tu_)}{};
    Parameters_distribution_new_( Parameters_distribution_new_&&)=default;
    Parameters_distribution_new_& operator=(const Parameters_distribution_new_& other)
    {if (this!=&other){
            Parameters_distribution_new_ tmp(other); *this=std::move(tmp);} return *this;
    }
    Parameters_distribution_new_& operator=( Parameters_distribution_new_&&)=default;



    auto& get_data()const {return tu_;}

protected:

    std::vector<tuple> tu_;


    static
        std::vector<tuple>
        tu2tp(const std::vector<tuple_ptr>& in)
    {
        std::vector<tuple> out;
        for (auto& e:in)
            out.emplace_back(std::get<0>(e),std::get<1>(e), std::get<2>(e),std::get<3>(e));
        return out;

    }

    std::vector<tuple>
    tucopy(const std::vector<tuple>& in
           )
    {
        std::vector<tuple> out;
        for (auto& e:in)
            out.emplace_back(std::get<0>(e),
                             std::unique_ptr<Base_Transformation<double>>(std::get<1>(e)->clone()),
                             std::unique_ptr<Base_Distribution<double>>(std::get<2>(e)->clone()),
                             std::get<3>(e));
        return out;

    }

};



template<template<class...> class Tr,typename Model>
myOptional_t<typename Tr<Parameters_values_new<Model>>::type>
tr_to_Parameter_impl(const Parameters_distribution_new_<Tr,Model>* self,const M_Matrix<double>& val)
{
    typedef myOptional_t<typename Tr<Parameters_values<Model>>::type> Op;
    typename Tr<LabelMap_t<typename Model::par_label>>::type m;
    std::vector<Op_void> res;
    for ( std::size_t i=0; i<val.size(); ++i)
    {
        res.emplace_back(self->is_in_range(i,val[i]));
        if (res[i])
        {
            m[self->name(i)]=self->tr(i)->apply_inv(val[i]);
        }
    }
    auto r=consolidate(std::move(res));
    if (r)
        return Op(Parameters_values<Model>(m));
    else return Op(false,r.error());
}



template<class Model>
class Parameters_partial_distribution: public Parameters_distribution<Model>
{


    // Base_Distribution interface
public:
    typedef  Parameters_partial_distribution<Model> self_type;

    typedef Parameters_distribution<Model> base_type;
    constexpr static auto  className=my_trait<Model>::className+my_static_string("_Parameters_partial_Distribution");

    virtual Parameters_partial_distribution *clone() const override
    {return new Parameters_partial_distribution(*this);}
    virtual std::string myClass() const override
    {
        return className.str();
    }

    virtual myOptional_t<Parameters_values<Model>> tr_to_Parameter(const M_Matrix<double>& val)const override
    {
      typedef myOptional_t<Parameters_values<Model>> Op;
        assert(val.size()==base_type::size());

        auto out=fixed_values_;
        auto var=base_type::tr_to_Parameter(val);
        if (var)
        {
          auto v=std::move(var).value();
          return Op(out<<v);
        }
        else
          return var;
    }

    Parameters_distribution<Model> const & variable_parameters()const { return *this;}

    Parameters_values<Model> const & fixed_parameters()const { return fixed_values_;}

  Parameters_partial_distribution(Parameters_distribution<Model>&& variable,
                                  Parameters_values<Model>&& fixed):
      Parameters_distribution<Model>{std::move(variable)}, fixed_values_{std::move(fixed)}{}

  Parameters_partial_distribution(const Parameters_distribution<Model>& variable,
                                  const Parameters_values<Model>& fixed):
      Parameters_distribution<Model>{std::move(variable)}, fixed_values_{std::move(fixed)}{}

  static auto get_constructor_fields()
  {
      return std::make_tuple(
                  grammar::field(C<self_type>{},"variable_parameters",&self_type::variable_parameters),
                  grammar::field(C<self_type>{},"fixed_parameters",&self_type::fixed_parameters));
  }

   Parameters_partial_distribution()=default;

private:

      Parameters_values<Model> fixed_values_;
};



template<class Model>
struct Derived_types<Parameters_distribution<Model>>{
  typedef Cs<Parameters_distribution<Model>,Parameters_partial_distribution<Model>> type;
    constexpr bool static value=true;

};





#endif // MYPARAMETERS_H
