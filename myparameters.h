#ifndef MYPARAMETERS_H
#define MYPARAMETERS_H

#include "mytypetraits.h"
#include "mygrammar.h"
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



#endif // MYPARAMETERS_H
