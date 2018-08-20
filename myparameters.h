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
    typedef T referred_type;
    std::string name()const {return name_;} // unique at all levels

    Label(std::string&& label): name_{label}{}
    Label(const std::string& label): name_{label}{}

    Label()=default;

    operator std::string()const {return name_;}

    bool operator < (const Label& other)const
    {
        return name()<other.name();
    }


    std::ostream& write(std::ostream& os)const { return os<<name_;}
    std::istream& read(std::istream& is) { return is>>name_;}

};
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
    typedef typename label::referred_type referred_type;
    typedef   std::map<std::string,referred_type > type;
};

template <class label>
using LabelMap_t=typename LabelMap<label>::type;


class Parameter_label
        : public Label<double> {public: using Label::Label;  };

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
    typedef typename Model::myParameter_label par_label;
    LabelMap_t<par_label> d_;
public:

    std::size_t size()const { return d_.size();}
    Parameters_values(LabelMap_t<par_label> values): d_{values}{};
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


};



template<typename Model>
class Parameters_distribution: public Base_Distribution<M_Matrix<double>>
{
public:
    virtual Parameters_distribution* clone()const override{ return new Parameters_distribution(*this);};
    std::string myClass()const override { return className.str();}

    typedef typename Model::myParameter_label par_label;

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
        M_Matrix<double> out(size(),size(),M_Matrix<double>::DIAGONAL);
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

    Parameters_values<Model> tr_to_Parameter(const M_Matrix<double>& val)const
    {
        assert(val.size()==tu_.size());
        LabelMap_t<par_label> m;
        for ( std::size_t i=0; i<val.size(); ++i)
        {
            m[name(i)]=tr(i)->apply_inv(val[i]);
        }
        return Parameters_values<Model>(m);
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

    std::size_t size()const { return tu_.size();}
    
    auto getParameterDistributionMap()const {
        std::vector<std::tuple<std::string,Base_Transformation<double>*, Base_Distribution<double>*>> out;
        for (std::size_t i=0; i<size(); ++i)
            out[i]={name(i),tr(i),dist(i)};
        return out;}
    typedef  Parameters_distribution<Model> self_type;
    constexpr static auto  className=my_trait<Model>::className+my_static_string("_Parameters_Distribution");
    static auto get_constructor_fields()
    {
        return std::make_tuple(
                    grammar::field(C<self_type>{},"values",&self_type::getParameterDistributionMap));
    }

    Parameters_distribution(const std::vector<std::tuple<std::string,Base_Transformation<double>*, Base_Distribution<double>*>>& in)
        :tu_{tu2tp(in)}{}

    Parameters_distribution()=default;

    Parameters_distribution(const Parameters_distribution& other): tu_{tucopy(other.tu_)}{};
    Parameters_distribution( Parameters_distribution&&)=default;
    Parameters_distribution& operator=(const Parameters_distribution& other)
    {if (this!=&other){
            Parameters_distribution tmp(other); *this=std::move(tmp);} return *this;
    }
    Parameters_distribution& operator=( Parameters_distribution&&)=default;


private:
    std::vector<std::tuple<std::string,std::unique_ptr<Base_Transformation<double>>, std::unique_ptr<Base_Distribution<double>>>> tu_;

    static
    std::vector<std::tuple<std::string,std::unique_ptr<Base_Transformation<double>>, std::unique_ptr<Base_Distribution<double>>>>
    tu2tp(const std::vector<std::tuple<std::string,Base_Transformation<double>*, Base_Distribution<double>*>>& in)
    {
        std::vector<std::tuple<std::string,std::unique_ptr<Base_Transformation<double>>, std::unique_ptr<Base_Distribution<double>>>> out;
        for (auto& e:in)
            out.emplace_back(std::get<0>(e),std::get<1>(e), std::get<2>(e));
        return out;

    }

    std::vector<std::tuple<std::string,std::unique_ptr<Base_Transformation<double>>, std::unique_ptr<Base_Distribution<double>>>>
    tucopy(const std::vector<std::tuple<std::string,std::unique_ptr<Base_Transformation<double>>, std::unique_ptr<Base_Distribution<double>>>>& in
           )
    {
        std::vector<std::tuple<std::string,std::unique_ptr<Base_Transformation<double>>, std::unique_ptr<Base_Distribution<double>>>> out;
        for (auto& e:in)
            out.emplace_back(std::get<0>(e),
                             std::unique_ptr<Base_Transformation<double>>(std::get<1>(e)->clone()),
                             std::unique_ptr<Base_Distribution<double>>(std::get<2>(e)->clone()));
        return out;

    }

};







#endif // MYPARAMETERS_H
