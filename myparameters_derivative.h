#ifndef MYPARAMETERS_DERIVATIVE_H
#define MYPARAMETERS_DERIVATIVE_H
#include "myparameters.h"
#include "matrixderivative.h"





template<typename Model>
class Derivative<Parameters_values<Model>>
{
    typedef typename Model::myParameter_label par_label;
    Derivative_t<LabelMap_t<par_label>> d_;
public:
  auto& x()const { return d_.begin()->second.x();}

    std::size_t size()const { return d_.size();}
    Derivative(Derivative_t<LabelMap_t<par_label>> values): d_{values}{};
    Derivative()=default;
    auto& at(const par_label& label)
    {
        return d_.at(label);
    }

    auto& at(const par_label& label)const
    {
        return d_.at(label);
    }

    auto& getParameterMap()const { return d_;}
    typedef  Derivative self_type;
    constexpr static auto  className=my_trait<Model>::className+my_static_string("_Parameters_derivative");
    static auto get_constructor_fields()
    {
        return std::make_tuple(
                    grammar::field(C<self_type>{},"values",&self_type::getParameterMap));
    }

    Derivative& operator<<(Derivative&& other)
    {
        d_.merge(std::move(other.getParameterMap()));
        return *this;
    }
    Derivative& operator<<(const Derivative& other)
    {
        for (auto e:other.getParameterMap())
            d_[e.first]=e.second;
        return *this;
    }

};


template<typename Model>
class Derivative<Parameters_distribution<Model>>: public Parameters_distribution<Model>
{
public:
  typedef Parameters_distribution<Model> base_type;

     base_type const& f()const { return static_cast<base_type const &>(*this);}
     Derivative (const base_type& p): base_type(p){}
     Derivative<Parameters_values<Model>> tr_to_Parameter_derivative(const M_Matrix<double>& val)const
    {
        assert(val.size()==base_type::tu_.size());
        Derivative_t<LabelMap_t<typename base_type::par_label>> m;
        for ( std::size_t i=0; i<val.size(); ++i)
        {
            M_Matrix<double> d(val.nrows(),val.ncols(),val.type(),0.0);
            d[i]=base_type::tr(i)->dapply_inv(val[i]);

            m[base_type::name(i)]=Derivative<double>(base_type::tr(i)->apply_inv(val[i]),val,
                                                     std::move(d));
        }
        return Derivative<Parameters_values<Model>>(m);
    }
    using base_type::tr_to_Parameter;
    Derivative()=default;
};



#endif // MYPARAMETERS_DERIVATIVE_H
