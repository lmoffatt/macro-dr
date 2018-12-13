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

    Parameters_values<Model> f()const
    {
      LabelMap_t<par_label> m;
      for (auto & e: d_)
        m.insert(e.first,e.second.f());
      return Parameters_values<Model>(m);
    }
    Parameters_values<Model> f_dfdx(std::size_t i,double dx)const
    {
      LabelMap_t<par_label> m;
      for (auto &e : d_)
        m.insert(e.first, e.second.f()+e.second.dfdx()[i]*dx);
      return Parameters_values<Model>(m);
    }

    Derivative Directional_Derivative(const M_Matrix<double> &newx, std::size_t i,
                                                    double eps)const
    {
      Derivative_t<LabelMap_t<par_label>> m;
      for (auto &e : d_)
        m.emplace(e.first, Derivative<double>(e.second.f()+e.second.dfdx()[i]*eps,newx,M_Matrix<double>(1,1,e.second.dfdx()[i])));
      return Derivative(m);
    }


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
     myOptional_t<Derivative<Parameters_values<Model>>> tr_to_Parameter_derivative(const M_Matrix<double>& val)const
    {
        assert(val.size()==base_type::tu_.size());

        typedef myOptional_t<Derivative<Parameters_values<Model>>> Op;
        Derivative_t<LabelMap_t<typename base_type::par_label>> m;
        std::vector<Op_void> res;
        for ( std::size_t i=0; i<val.size(); ++i)
        {
          res.emplace_back(base_type::is_in_range(i,val[i]));
          if (res[i])
          {
            M_Matrix<double> d(val.nrows(),val.ncols(),val.type(),0.0);
            d[i]=base_type::tr(i)->dapply_inv(val[i]);

            m[base_type::name(i)]=Derivative<double>(base_type::tr(i)->apply_inv(val[i]),val,
                                                     std::move(d));
          }
        }
        auto r=consolidate(std::move(res));
        if (r)

          return Op(Derivative<Parameters_values<Model>>(m));
        else {
          return Op(false,r.error());
        }
    }
    using base_type::tr_to_Parameter;
    Derivative()=default;
};

template<typename Model>
inline Parameters_values<Model> Taylor_first(const Derivative<Parameters_values<Model>> &dx, std::size_t i,
                           double eps) {
  return dx.f_dfdx(i,eps);
}



template<typename Model>
inline const Derivative<Parameters_values<Model>>  Directional_Derivative(const Derivative<Parameters_values<Model>> &dx,
                                                 const M_Matrix<double> &new_x,
                                                 std::size_t i, double eps) {

  return dx.Directional_Derivative(new_x,i,eps);
}





#endif // MYPARAMETERS_DERIVATIVE_H
