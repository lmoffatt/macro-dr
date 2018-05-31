#ifndef QMODEL_H
#define QMODEL_H

#include "Matrix.h"
#include "myfields.h"
#include "mymath.h"
#include "mytypetraits.h"
#include "myparameters.h"
#include <vector>
#include <map>
#include <set>
#include <cmath>


class Model_Parameter_label
        : public Parameter_label {
public: using Parameter_label::Parameter_label;
    constexpr static auto className=my_static_string("Model_Parameter_label");

};

class rate_Parameter_label
        : public Model_Parameter_label {
public: using Model_Parameter_label::Model_Parameter_label;
    constexpr static auto className=my_static_string("rate_Parameter_label");
};

class Equilibrium_Parameter_label
        : public Model_Parameter_label {public: using Model_Parameter_label::Model_Parameter_label;
                                        constexpr static auto className=my_static_string("Equilibrium_Parameter_label");
                                       };

class Time_Constant_Parameter_label
        : public Model_Parameter_label {public: using Model_Parameter_label::Model_Parameter_label;
                                        constexpr static auto className=my_static_string("Time_Constant_Parameter_label");
                                       };

class Coupling_factor_Parameter_label
        : public Model_Parameter_label {public: using Model_Parameter_label::Model_Parameter_label;
                                        constexpr static auto className=my_static_string("Coupling_factor_Parameter_label");
                                       };

class Coupling_coefficient_Parameter_label
        : public Model_Parameter_label {public: using Model_Parameter_label::Model_Parameter_label;
                                        constexpr static auto className=my_static_string("Coupling_coefficient_Parameter_label");
                                       };

class Conductance_Parameter_label
        : public Model_Parameter_label {public: using Model_Parameter_label::Model_Parameter_label;
                                        constexpr static auto className=my_static_string("Conductance_Parameter_label");
                                       };

typedef rate_Parameter_label::referred_type hiho;

inline static rate_Parameter_label p={"gsds"};

class Conformational_change;
class Conformational_change_label
        : public Label<Conformational_change> {
public: using Label::Label;
    constexpr static auto className=my_static_string("Conformational_change_label");

};

class Conformational_change
{
    int change_in_agonist_;
    int change_in_conductance_;
    Conformational_change_label label_;
    rate_Parameter_label par_on_;
    rate_Parameter_label par_off_;
    /*   Equilibrium_Parameter_label par_Eq_;
    Time_Constant_Parameter_label par_tc_; */
public:
    typedef  Conformational_change self_type;
    Conformational_change_label label()const {return label_;}
    const rate_Parameter_label& par_on()const {return par_on_;}
    const rate_Parameter_label& par_off()const {return par_off_;}
    /*  const Equilibrium_Parameter_label&  par_Eq()const {return par_Eq_;};
    const Time_Constant_Parameter_label& par_tc()const {return par_tc_;}*/


    int change_in_agonist()const {return change_in_agonist_;}
    int change_in_conductance()const { return change_in_conductance_;}
    Conformational_change()=default;
    Conformational_change(int _change_in_agonist,
                          int _change_in_conductance,
                          Conformational_change_label _label,
                          rate_Parameter_label _par_on,
                          rate_Parameter_label _par_off /*,
                                                                                            Equilibrium_Parameter_label&& _par_Eq,
                                                                                            Time_Constant_Parameter_label&& _par_tc*/):
        change_in_agonist_{_change_in_agonist},change_in_conductance_{_change_in_conductance},
        label_{std::move(_label)},par_on_{std::move(_par_on)},par_off_{std::move(_par_off)}/*,par_Eq_{std::move(_par_Eq)},par_tc_{std::move(_par_tc)}*/{}


    constexpr static auto className=my_static_string("Conformational_change");

    static auto get_constructor_fields()
    {
        return std::make_tuple(grammar::field(C<self_type>{},"change_in_agonist", &self_type::change_in_agonist),
                               grammar::field(C<self_type>{},"change_in_conductance", &self_type::change_in_conductance),
                               grammar::field(C<self_type>{},"label",&self_type::label),
                               grammar::field(C<self_type>{},"par_on", &self_type::par_on),
                               grammar::field(C<self_type>{},"par_off", &self_type::par_off)
                               /*,
                                                                                                    grammar::field(C<self_type>{},"par_Eq", &self_type::par_Eq),
                                                                                                    grammar::field(C<self_type>{},"par_tc", &self_type::par_tc)*/);
    }
};
class Conformational_interaction;
class Conformational_interaction_label
        : public Label<Conformational_interaction> {public: using Label::Label;};



class Conformational_interaction
{
    std::vector<Conformational_change_label> conf_changes_;
    Coupling_factor_Parameter_label factor_;
    std::vector<Coupling_coefficient_Parameter_label> par_inter_f_;

public:
    std::vector<Conformational_change_label> const & interacting_conformational_changes()const {return conf_changes_;}
    Coupling_factor_Parameter_label const& factor_label()const {return factor_;};
    std::vector<Coupling_coefficient_Parameter_label> const &coefficient_labels()const {return par_inter_f_;}

    typedef  Conformational_interaction self_type;

    Conformational_interaction()=default;
    Conformational_interaction(std::vector<Conformational_change_label> _conf_changes,
                               Coupling_factor_Parameter_label _factor,
                               std::vector<Coupling_coefficient_Parameter_label> _par_inter_f):
        conf_changes_{std::move(_conf_changes)},factor_{std::move(_factor)},par_inter_f_{std::move(_par_inter_f)}
    {}
    static auto get_constructor_fields()
    {
        return std::make_tuple(
                    grammar::field(C<self_type>{},"interacting_conformational_changes",&self_type::interacting_conformational_changes),
                    grammar::field(C<self_type>{},"factor_label", &self_type::factor_label),
                    grammar::field(C<self_type>{},"coefficient_labels", &self_type::coefficient_labels));
    }
    constexpr static auto className=my_static_string("Conformational_interaction");

    bool operator<(const Conformational_interaction& other) const
    {
        return factor_label()<other.factor_label();
    }

};


/**
 * @brief The Allosteric_Model class build a kinetic rate model starting with a set of conformational changes and their interactions.
 *
 *
 */
class Allosteric_Model
{
public:
    typedef  Allosteric_Model self_type;
    typedef Model_Parameter_label myParameter_label;
    struct transitions{
        bool on;
        bool agonist;
        std::string conformation;
        std::map<std::vector<std::pair<std::string, std::string>>, std::size_t> coupling;

    };
    struct model_definition
    {
        std::vector<std::string> conformational_changes;
        std::map<std::pair<std::string,bool>,std::string> conformational_changes_names;
        std::set<std::size_t> agonist_changes;
        std::set<std::size_t> conductance_changes;
        std::map<std::size_t, std::string> conductance_names;
        std::multimap<std::size_t, std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>> conformational_inter_unit_cell;
    };

    struct new_model_definition
    {
        std::size_t number_of_units;
        std::map<Conformational_change_label,Conformational_change> conformational_changes;
        std::vector<Conformational_change_label> unit_of_conformational_changes;
        std::set<Conformational_interaction> conformational_interactions;
        std::map<std::size_t, Conductance_Parameter_label> conductance_names;

    };

    model_definition new_to_old(const new_model_definition& m)
    {
        model_definition out;
        for (auto& e:m.conductance_names)
            out.conductance_names[e.first]=e.second;
        auto k=m.unit_of_conformational_changes.size();
        std::size_t n=m.number_of_units*k;
        out.conformational_changes.resize(n);
        for (std::size_t i=0; i<m.number_of_units; ++i)
            for (std::size_t j=0; j<m.unit_of_conformational_changes.size(); ++j)
            {
                auto ii=i*k+j;
                auto& cf=m.conformational_changes.at(m.unit_of_conformational_changes[j].name());
                out.conformational_changes[ii]=cf.label();
                if (cf.change_in_agonist()!=0)
                    out.agonist_changes.insert(ii);
                if (cf.change_in_conductance()!=0)
                    out.conductance_changes.insert(ii);
            }
        for (auto& e:m.conformational_changes)
        {
            out.conformational_changes_names[std::pair(e.second.label(), true)]=e.second.par_on();
            out.conformational_changes_names[std::pair(e.second.label(), false)]=e.second.par_off();
        }
        for (auto &e: m.conformational_interactions)
        {
            std::vector<std::size_t> cc;
            std::size_t current_i=0;
            for (std::size_t j=0; j<e.interacting_conformational_changes().size(); ++j)
            {
                auto j_n=m.conformational_changes.at(e.interacting_conformational_changes()[j]);
                std::size_t i=current_i;
                while (j_n.label().name()!=out.conformational_changes[i]) ++i;
                cc.push_back(i);
                ++current_i;
            }

            for (std::size_t i=0; i<cc.size(); ++i)
            {
                auto x=cc[i];
                auto ix=x%k;
                auto nx=x/k;
                auto shift=(n-nx)*k;
                std::set<std::size_t> s;
                for (std::size_t j=0; j<cc.size(); ++j)
                    if (j!=i) s.insert(cc[j]);
                s=rotate(s,n*k,shift);
                out.conformational_inter_unit_cell.emplace(ix,std::pair(s,std::pair(e.factor_label(),e.coefficient_labels()[i])));
            }

        }
        return out;

    }



    constexpr static auto  className=my_static_string("Allosteric_Model");



private:
    new_model_definition new_d_;
    std::vector<std::vector<std::string>> conformer_;
    std::vector<std::map<std::size_t,transitions>> transitions_;
    std::set<std::string> paramNames_;
    std::vector<std::string> conductances_;
    model_definition d_;

    /**
     * @brief getParameterNamesFrom
     * @param conformational_changes_names
     * @param conductance_names
     * @param conformational_inter
     */
    static auto getParameterNamesFrom(const std::map<std::pair<std::string,bool>,std::string> conformational_changes_names,
                                      const std::map<std::size_t, std::string>& conductance_names,
                                      const std::multimap<std::size_t, std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>>& conformational_inter)
    {
        std::set<std::string> out;
        for (auto &e: conformational_changes_names)
            out.insert(e.second);
        for (auto& e: conductance_names)
            out.insert(e.second);
        for (auto &e : conformational_inter)
        {
            out.insert(e.second.second.first);
            out.insert(e.second.second.second);
        }
        return out;
    }



    static auto getConformers(const std::vector<std::string>& conformational_changes,std::size_t p)
    {
        std::vector<std::vector<std::string>> conformer;
        std::map<std::vector<std::string>,std::size_t> state_to_conformer;
        for (std::size_t i=0; i< (1u<<conformational_changes.size()); ++i)
        {
            auto c=index_to_conformational_change(conformational_changes,i);
            if (state_to_conformer.find(c)==state_to_conformer.end())
            {
                state_to_conformer[c]=conformer.size();
                std::size_t n=1;
                auto permute=rotate(c,p*n);
                while(c!=permute)
                {
                    state_to_conformer[permute]=conformer.size();
                    ++n;
                    permute=rotate(c,p*n);
                }
                conformer.push_back(c);

            }

        }
        return std::make_tuple(conformer, state_to_conformer);
    }


    static
    auto getTransitions(const std::vector<std::string>& conformational_changes,
                        const std::vector<std::vector<std::string>>& conformer,
                        const  std::map<std::vector<std::string>,std::size_t>& state_to_conformer,
                        const std::set<std::size_t>& agonist_changes,
                        const std::map<std::pair<std::string,bool>,std::string> conformational_changes_names,
                        std::multimap<std::size_t, std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>> conformational_interactions)
    {
        std::vector<std::map<std::size_t,transitions>> transitions_out;
        for (std::size_t i=0; i<conformer.size(); ++i)
        {
            std::map<std::size_t,transitions> myTransition;
            auto c=conformer[i];
            for (std::size_t k=0; k<c.size(); ++k)
            {
                auto change=c;
                if (c[k]=="")
                    change[k]=conformational_changes[k];
                else
                    change[k]="";
                auto j=state_to_conformer.at(change);
                if (myTransition.find(j)==myTransition.end())
                {

                    myTransition[j].agonist=agonist_changes.find(k)!=agonist_changes.end();
                    bool on=change[k]==conformational_changes[k];
                    myTransition[j].on=on;
                    auto name=conformational_changes_names.at(std::pair{conformational_changes[k],on});
                    myTransition[j].conformation=name;
                }
                auto m=conformational_interactions.equal_range(k);
                std::vector<std::pair<std::string, std::string>> coupling;
                for (auto it=m.first; it!=m.second; ++it)
                {
                    std::set<std::size_t> s=it->second.first;
                    bool all=true;
                    for (auto e:s)
                    {
                        if (c[e].empty()) {
                            all=false; break;
                        }
                    }
                    if (all)
                        coupling.push_back(it->second.second);


                }
                ++myTransition[j].coupling[coupling];

            }
            transitions_out.push_back(std::move(myTransition));
        }
        return transitions_out;

    }



    static

    std::string g_name_of_conformer(const std::vector<std::string>& c, const std::set<std::size_t>& conductance_changes,const std::map<std::size_t, std::string>& conductance_names)
    {
        std::size_t i=0;
        for (auto e:conductance_changes)
            if (!c.at(e).empty()) ++i;
        return conductance_names.at(i);
    }


    static     std::vector<std::string> get_g_names(const std::vector<std::vector<std::string>>& conformers, const std::set<std::size_t>& conductance_changes,const std::map<std::size_t, std::string>& conductance_names)
    {
        std::vector<std::string> out(conformers.size());
        for (std::size_t i=0; i<conformers.size(); ++i)
            out[i]=g_name_of_conformer(conformers[i],conductance_changes, conductance_names);
        return out;

    }

public:


    static auto get_constructor_fields_old()
    {
        return std::make_tuple(
                    grammar::field(C<self_type>{},"conformational_changes",&self_type::get_conformational_changes),
                    grammar::field(C<self_type>{},"conformational_changes_names", &self_type::get_conformational_changes_names),
                    grammar::field(C<self_type>{},"agonist_changes", &self_type::get_agonist_changes),
                    grammar::field(C<self_type>{},"conductance_changes", &self_type::get_conductance_changes),
                    grammar::field(C<self_type>{},"conductance_names", &self_type::get_conductance_names),
                    grammar::field(C<self_type>{},"conformational_inter_unit_cell", &self_type::get_conformational_inter_unit_cell));
    }


    const std::vector<std::string>& get_conformational_changes()const {return d_.conformational_changes;}
    const std::map<std::pair<std::string,bool>,std::string>& get_conformational_changes_names()const{return d_.conformational_changes_names;}
    const std::set<std::size_t>& get_agonist_changes()const{return d_.agonist_changes;}
    const std::set<std::size_t>& get_conductance_changes()const {return d_.conductance_changes;}
    const std::map<std::size_t, std::string>& get_conductance_names()const{ return d_.conductance_names;}
    const std::multimap<std::size_t, std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>>& get_conformational_inter_unit_cell()const{return  d_.conformational_inter_unit_cell;}


    std::size_t number_of_units()const {return new_d_.number_of_units;}
    std::map<Conformational_change_label,Conformational_change> const &conformational_changes()const {return new_d_.conformational_changes;}
    std::vector<Conformational_change_label> const &unit_of_conformational_changes()const { return new_d_.unit_of_conformational_changes;}
    std::set<Conformational_interaction> const & conformational_interactions() const { return new_d_.conformational_interactions; }
    std::map<std::size_t, Conductance_Parameter_label> const &conductance_names()const { return  new_d_.conductance_names;}


    static auto get_constructor_fields()
    {
        return std::make_tuple(
                    grammar::field(C<self_type>{},"number_of_units",&self_type::number_of_units),
                    grammar::field(C<self_type>{},"conformational_changes", &self_type::conformational_changes),
                    grammar::field(C<self_type>{},"unit_of_conformational_changes", &self_type::unit_of_conformational_changes),
                    grammar::field(C<self_type>{},"conformational_interactions", &self_type::conformational_interactions),
                    grammar::field(C<self_type>{},"conductance_names", &self_type::conductance_names));
    }



    Allosteric_Model()=default;



    Allosteric_Model(std::size_t number_of_units,
                     std::map<Conformational_change_label,Conformational_change> conformational_changes,
                     std::vector<Conformational_change_label> unit_of_conformational_changes,
                     std::set<Conformational_interaction> conformational_interactions,
                     std::map<std::size_t, Conductance_Parameter_label> conductance_names):
        new_d_{number_of_units,conformational_changes,unit_of_conformational_changes,conformational_interactions,conductance_names},d_{new_to_old(new_d_)}
    {
        init();
    }



    Allosteric_Model(const std::vector<std::string>& conformational_changes,
                     const std::map<std::pair<std::string,bool>,std::string> conformational_changes_names,
                     const std::set<std::size_t>& agonist_changes,
                     const std::set<std::size_t>& conductance_changes,
                     const std::map<std::size_t, std::string>& conductance_names,
                     const std::multimap<std::size_t, std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>>& conformational_inter_unit_cell)
        :d_{conformational_changes,conformational_changes_names,agonist_changes,conductance_changes,conductance_names,conformational_inter_unit_cell}
    {
        init();
    }




    template<class Parameters>
    auto Qs(const Parameters& p) const
    {

        M_Matrix<double> Q0(conformer_.size(),conformer_.size());
        M_Matrix<double> Qa(conformer_.size(),conformer_.size());
        for (std::size_t i=0; i<Q0.nrows(); ++i)
        {
            Q0(i,i)=0;
            Qa(i,i)=0;
            for (auto it=transitions_[i].begin(); it!=transitions_[i].end(); ++it)
            {
                std::size_t j=it->first;
                if ((it->second.agonist)&&(it->second.on))
                {   Qa(i,j)=rate(it->second,p);
                    Qa(i,i)-=Qa(i,j);
                }else{
                    Q0(i,j)=rate(it->second,p);
                    Q0(i,i)-=Q0(i,j);
                }}
        }
        auto r=Q0*ones<double>(Q0.ncols(),1);
        auto q=Qa*ones<double>(Q0.ncols(),1);

        return std::pair(Q0,Qa);
    }


    template<class Parameters>
    auto g(const Parameters& p) const
    {
        M_Matrix<double> out(conformer_.size(),1);
        for (std::size_t i=0; i<conformer_.size(); ++i)
            out(i,0)=p.at(conductances_.at(i));
        return out;
    }




    auto getParameterNames()const
    {
        return paramNames_;
    }


private:
    void init()
    {
        paramNames_=getParameterNamesFrom(d_.conformational_changes_names,d_.conductance_names,d_.conformational_inter_unit_cell);
        std::cout<<paramNames_;
        auto p=periodicity(d_.conformational_changes);

        auto conformational_interactions=fill_conformational_interactions(d_.conformational_inter_unit_cell,d_.conformational_changes.size(),p);

        std::map<std::vector<std::string>,std::size_t> state_to_conformer;

        std::tie(conformer_, state_to_conformer)=getConformers(d_.conformational_changes,p);

        transitions_=getTransitions(d_.conformational_changes,conformer_,state_to_conformer,d_.agonist_changes,d_.conformational_changes_names,conformational_interactions);
        conductances_=get_g_names(conformer_,d_.conductance_changes,d_.conductance_names);

    }


    static
    std::vector<std::string> index_to_conformational_change(const std::vector<std::string>& cc, std::size_t index)
    {
        std::vector<std::string> out(cc.size(),"");
        for (std::size_t i=0; i<cc.size(); ++i)
            if ( ((index>>i) & 1)==1)   out[i]=cc[i];
        return out;
    }

    static std::vector<std::string> rotate(const std::vector<std::string>& x, std::size_t n)
    {
        std::vector<std::string> out(x.size());
        for (std::size_t i=0; i<out.size(); ++i)
            out[(i+n)%out.size()]=x[i];
        return out;
    }




    static std::size_t periodicity(const std::vector<std::string> & x)
    {
        std::size_t p=1;
        auto t=rotate(x,p);
        while(t!=x && p<x.size())
        {
            ++p;
            t=rotate(x,p);
        }
        return p;
    }

    static std::set<std::size_t> rotate(const std::set<std::size_t>& c , std::size_t n, std::size_t i)
    {
        std::set<std::size_t> out;
        for (auto e:c)
            out.insert((e+i)%n);
        return out;
    }


    static
    std::multimap<std::size_t, std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>>
    fill_conformational_interactions(const std::multimap<std::size_t, std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>>& conformational_interactions,
                                     std::size_t n, std::size_t p)
    {
        std::multimap<std::size_t, std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>> out;
        for (std::size_t i=0; i<n/p; ++i)
        {
            for (auto& e:conformational_interactions)
            {
                out.insert({(e.first+i*p)%n,{rotate(e.second.first,n,i*p),e.second.second}});
            }
        }
        return out;
    }



    template <class P>
    static double rate(const transitions& tr,const P& p)
    {
        double out=0;
        if (tr.on)
        {
            for (auto& e:tr.coupling)
            {
                double b=e.second;
                for (auto& e2:e.first)
                    b*=std::pow(p.at(e2.first),p.at(e2.second));
                out+=b;
            }
        }
        else
        {
            for (auto& e:tr.coupling)
            {
                double b=e.second;
                for (auto& e2:e.first)
                    b*=std::pow(p.at(e2.first),p.at(e2.second)-1.0);
                out+=b;
            }
        }
        return out*p.at(tr.conformation);
    }

};




class Markov_Transition_rate
{

    M_Matrix<double> Qrun_;// transition rate matrix at time zero
    M_Matrix<double> g_;

    M_Matrix<double> V_;// eigenvector of Qrun
    M_Matrix<double> W_;// eigenvector of Qrun
    M_Matrix<double> landa_; // eigenvalues

    M_Matrix<double> Wg_;
    M_Matrix<double> WgV_;

    void clean_landa(M_Matrix<double>& la)
    {
        for (std::size_t i=1; i<la.size(); ++i)
            if (la[i]>0) {la[i]=0; }
    }

    void init()
    {
        std::tie(W_,landa_,V_)=Matrix_Decompositions::EigenSystem_full_real_eigenvalue(Qrun_);

        clean_landa(landa_);
        Wg_=W_*g_;
        WgV_=W_*Matrix_Unary_Transformations::diag(g_)*V_;
    }

public:
    Markov_Transition_rate()=default;
    ///virtual copy constructors


    Markov_Transition_rate(M_Matrix<double>&& _Qrun,
                           const M_Matrix<double>& _g):Qrun_{std::move(_Qrun)},g_{_g}{ init();}


    const M_Matrix<double>&  Qrun()const { return Qrun_;}// transition rate matrix at time zero
    const M_Matrix<double>&  V()const{ return V_;}// eigenvector of Qrun
    const M_Matrix<double>&  W()const{ return W_;}// eigenvector of Qrun
    const M_Matrix<double>&  landa()const{ return landa_;} // eigenvalues

    const M_Matrix<double>&  g()const{ return g_;}
    const M_Matrix<double>&  Wg()const{ return Wg_;}
    const M_Matrix<double>&  WgV()const{ return WgV_;}

    M_Matrix<double>  Peq()const
    {
        M_Matrix<double> p0(1,Qrun().nrows(),1.0/Qrun().nrows());

        auto OO=Matrix_Generators::ones<double>(Qrun_.nrows(),Qrun_.ncols());
        auto SS=inv(Qrun()+OO);
        if (!SS.second.empty())

        {

            return p0*V()*exp(landa()*1e8)*W();
        }
        else
            return p0*SS.first;

    }
};



class Markov_Transition_step
{
    M_Matrix<double> P_;// transition matrix

    M_Matrix<double> ladt_; // exp(la dt)

    M_Matrix<double> exp_ladt_; // exp(la dt)

    M_Matrix<double> g_; // conductance matrix


    double dt_;

    std::size_t nsamples_;
public:

    M_Matrix<double> const&  P() const {return P_;}// transition matrix

    double dt()const {return dt_;}

    M_Matrix<double> const & ladt(){return ladt_;}; // exp(la dt)

    M_Matrix<double> const& exp_ladt(){return exp_ladt_;} // exp(la dt)


    std::size_t nsamples(){ return nsamples_;}

    /// mean conductance for each starting state i
    M_Matrix<double> const&  gmean_i() const {return g_;} // conductance matrix

    Markov_Transition_step(const Markov_Transition_rate& Qx, std::size_t nsamples, double fs)
        :P_{},ladt_{},exp_ladt_{},g_{Qx.g()},dt_{1.0/fs*nsamples},nsamples_{nsamples}{
        init(Qx);
    }
    Markov_Transition_step()=default;

private:
    void init(const Markov_Transition_rate& Qx)
    {
        ladt_=Qx.landa()*dt();
        exp_ladt_=exp(ladt_);
        P_=Qx.V()*exp_ladt_*Qx.W();
        g_=Qx.g();
    }

};


class Markov_Transition_step_double: public Markov_Transition_step
{
    M_Matrix<double> gmean_i_; // conductance matrix
    M_Matrix<double> gtotal_ij_; // conductance matrix
    M_Matrix<double> gmean_ij_; // conductance matrix

    M_Matrix<double> gtotal_sqr_ij_; // conductance matrix

    M_Matrix<double> gsqr_i_; // conductance matrix

    M_Matrix<double> gvar_i_; //variance of the conductance matrix
    M_Matrix<double> gtotal_var_ij_; //variance of the conductance matrix

    M_Matrix<double> gtotal_var_i_; //variance of the conductance matrix

    M_Matrix<double> gvar_ij_; //variance of the conductance matrix

public:


    /// mean conductance for each starting state i
    M_Matrix<double> const&  gmean_i() const {return gmean_i_;} // conductance matrix
    ///total conductance for each starting state i and ending state j
    M_Matrix<double> const&  gtotal_ij() const {return gtotal_ij_;} // conductance matrix
    ///mean conductance for each starting state i and ending state j
    M_Matrix<double> const&  gmean_ij() const {return gmean_ij_;} // conductance matrix

    /// squared mean conductance for each starting state i and ending state j
    M_Matrix<double> const&  gtotal_sqr_ij() const {return gtotal_sqr_ij_;} // conductance matrix

    /// squared mean conductance for each starting state i
    M_Matrix<double> const&  gsqr_i() const {return gsqr_i_ ;} // conductance matrix

    /// variance of the mean conductance for each starting state i
    M_Matrix<double> const&  gvar_i() const {return gvar_i_;} //variance of the conductance matrix
    /// variance of the mean conductance for each starting state i contributed by the ones ending at state j
    M_Matrix<double> const&  gtotal_var_ij() const {return gtotal_var_ij_;} //variance of the conductance matrix

    /// variance of the mean conductance for each starting state i summed all over  j
    M_Matrix<double> const&  gtotal_var_i() const {return gtotal_var_i_;} //variance of the conductance matrix

    /// variance of the mean conductance for each starting state i and ending state j
    M_Matrix<double> const&  gvar_ij() const {return gvar_ij_;} //variance of the conductance matrix

    Markov_Transition_step_double(const Markov_Transition_rate& Qx,std::size_t nsamples, double fs):
        Markov_Transition_step(Qx,nsamples,fs){init(Qx);}



private:

    static double E1(double x){
        if (std::abs(x)<std::numeric_limits<double>::epsilon()*100) return 1.0;
        else if (std::abs(x)<1e-2) return std::expm1(x)/x;
        else return (std::exp(x)-1.0)/x;
    }

    static double E2(double x, double y)
    {
        const double eps= std::numeric_limits<double>::epsilon();
        if (x*x<eps)
        {   if (y*y<eps)
                return 0.5;
            else
                return (E1(y)-1.0)/y;
        }
        else if (y*y<eps)
            return (E1(x)-1.0)/x;
        else if ((y-x)*(y-x)<eps)
            return  (std::exp(x)-E1(x))/x;
        else
            return    (E1(y)-E1(x))/(y-x);
    }

    static double Ee(double x, double y, double exp_x, double exp_y)
    {
        const double eps= std::numeric_limits<double>::epsilon();
        if (sqr(x-y)<eps)
            return exp_x;
        else
            return (exp_x-exp_y)/(x-y);
    };

    static double EX_111(double x, double y, double z, double exp_x)
    {
        return exp_x/((x-y)*(x-z));
    }

    static double E111(double x, double y, double z, double exp_x, double exp_y, double exp_z)
    {
        return EX_111(x,y,z,exp_x)+EX_111(y,x,z,exp_y)+EX_111(z,y,x, exp_z);
    }
    static double E12(double x,double y, double exp_x, double exp_y)
    {
        return EX_111(x,y,y,exp_x)+exp_y/(y-x)*(1.0-1.0/(y-x));
    }

    static double E3(double x, double y, double z, double exp_x, double exp_y, double exp_z)
    {
        const double eps= std::numeric_limits<double>::epsilon();
        if (sqr(x-y)<eps)   // x==y
        {
            if (sqr(y-z)<eps)   // y==z
                return exp_x/2.0;  // x==y==z
            else
                return E12(z,x,exp_z,exp_x); // x==y!=z
        }
        else if (sqr(y-z)<eps)   // x!=y==z
        {
            return E12(x,y,exp_x,exp_y);
        }
        else if(sqr(x-z)>eps)   // y!=z==x!=y
        {
            return E12(y,x,exp_y,exp_x);
        }
        else return E111(x,y,z,exp_x,exp_y,exp_z); // x!=y!=z!=x
    }


    void init(const Markov_Transition_rate& Qx)
    {
        //const double eps=std::numeric_limits<double>::epsilon();
        double dt=Markov_Transition_step::dt();

        std::size_t N=P().ncols();
        std::vector<double> Wg_E0(N);

        for (std::size_t k0=0; k0<N; k0++)
        {
            double rladt=Qx.landa()[k0]*dt;
            if (rladt<std::numeric_limits<double>::epsilon())
                Wg_E0[k0]=Qx.Wg()[k0];
            else
                Wg_E0[k0]=Qx.Wg()[k0]*(exp(rladt)-1.0)/rladt;
        };


        gmean_i_=M_Matrix<double>(N,1,0);
        for (std::size_t k0=0; k0<N; k0++)
        {
            for (std::size_t j=0; j<N; j++)
                gmean_i_[k0]+=Qx.V()(k0,j)*Wg_E0[j];
        }

        auto ladt=Qx.landa()*dt;
        auto exp_ladt=ladt.apply([](double x){return std::exp(x);});
        auto Wg_E0_=ladt.apply(&E1)*Qx.Wg();
        gmean_i_=Qx.V()*Wg_E0_;

        auto k_u=N;

        //build E2
        M_Matrix<double> WgV_Wg_E2(k_u,k_u);


        for (std::size_t i=0; i<N; ++i)
            for (std::size_t j=0; j<N; ++j)
                WgV_Wg_E2(i,j)=Qx.WgV()(i,j)*Qx.Wg()[j]*E2(ladt[i], ladt[j]);


        // vmean(i)=A(k,i,j) g(j) A(k2,j,j2) g(j2) E(k,k)
        // =V(i,k) W(k,j) g(j) V(j,k2) W(k2,j2) g(j2) E(k,k2)
        //         ___________________ ______________
        //         WgV(k,k2)           Wg(k2)
        //___________________________________________________

        // G_(i)=V(i,n1), W(n1,k)g(k) V(k,n2) W(n2,j)       E(n1,n2)
        //                ___________________ ______________
        //                WgV(n1,n2)          W(n2)
        //___________________________________________________


        M_Matrix<double> u_col(k_u,1,1);

        gsqr_i_=Qx.V()*(WgV_Wg_E2*u_col);
        gvar_i_=gsqr_i_-elemMult(gmean_i_,gmean_i_);



        M_Matrix<double> EE2(k_u,k_u);
        for (std::size_t k0=0; k0<k_u; k0++)
            for (std::size_t k2=0; k2<k_u; k2++)
                EE2(k0,k2)=Ee(ladt[k0],ladt[k2],exp_ladt[k0],exp_ladt[k2]);

        gtotal_ij_=Qx.V()*elemMult(Qx.WgV(),EE2)*Qx.W();



        M_Matrix<double> WgV_E3=Matrix_Generators::zeros<double>(k_u,k_u);
        for (std::size_t n1=0; n1<k_u; n1++)
            for (std::size_t n3=0; n3<k_u; n3++)
                for (std::size_t n2=0; n2<k_u; n2++)
                    WgV_E3(n1,n3)+=Qx.WgV()(n1,n2)*Qx.WgV()(n2,n3)*E3(ladt[n1],ladt[n2],ladt[n3],exp_ladt[n1],exp_ladt[n2],exp_ladt[n3]);
        gtotal_sqr_ij_=(Qx.V()*(WgV_E3*2.0)*Qx.W());

        gmean_ij_=M_Matrix<double>(k_u,k_u);
        gvar_ij_=M_Matrix<double>(k_u,k_u);
        for (std::size_t i=0; i<k_u; i++)
            for (std::size_t j=0; j<k_u; j++)
                if (P()(i,j)>1e-9)
                {
                    gmean_ij_(i,j)=gtotal_ij_(i,j)/P()(i,j);
                    gvar_ij_(i,j)=gtotal_sqr_ij_(i,j)/P()(i,j)-gmean_ij_(i,j)*gmean_ij_(i,j);
                }
                else
                {
                    gmean_ij_(i,j)=0;
                    gvar_ij_(i,j)=0;
                }
        M_Matrix<double> u=Matrix_Generators::ones<double>(k_u,1);
        gtotal_var_ij_=elemMult(gvar_ij_,P());
        gtotal_var_i_=gtotal_var_ij_*u;
    }

    void init_old(const Markov_Transition_rate& Qx, std::size_t nsamples, double fs)
    {
        const double eps=std::numeric_limits<double>::epsilon();

        double dt=1.0/fs*nsamples;

        std::size_t N=P().ncols();
        std::vector<double> Wg_E0(N);


        for (std::size_t k0=0; k0<N; k0++)
        {
            double rladt=Qx.landa()[k0]*dt;
            if (rladt<std::numeric_limits<double>::epsilon())
                Wg_E0[k0]=Qx.Wg()[k0];
            else
                Wg_E0[k0]=Qx.Wg()[k0]*(exp(rladt)-1.0)/rladt;
        };


        gmean_i_=M_Matrix<double>(N,1,0);
        for (std::size_t k0=0; k0<N; k0++)
        {
            for (std::size_t j=0; j<N; j++)
                gmean_i_[k0]+=Qx.V()(k0,j)*Wg_E0[j];
        }


        auto k_u=N;

        //build E2
        M_Matrix<double> WgV_Wg_E2(k_u,k_u);
        for (std::size_t k0=0; k0<k_u; k0++)
        {
            double rladt=Qx.landa()[k0]*dt;
            if (rladt*rladt<eps)
            {
                for (std::size_t k2=0; k2<k_u; k2++)
                {
                    double rla2dt=Qx.landa()[k2]*dt;
                    if (std::abs(rla2dt*rla2dt)<eps)
                        WgV_Wg_E2(k0,k2)=Qx.WgV()(k0,k2)*
                                Qx.Wg()[k2]*0.5;
                    else
                        WgV_Wg_E2(k0,k2)=Qx.WgV()(k0,k2)*Qx.Wg()[k2]*
                                (exp(rla2dt)-rla2dt-1.0)/rla2dt/rla2dt;
                }
            }
            else
            {
                for (std::size_t k2=0; k2<k_u; k2++)
                {
                    double rla2dt=Qx.landa()[k2]*dt;
                    if (rla2dt*rla2dt<eps)
                    {
                        WgV_Wg_E2(k0,k2)=Qx.WgV()(k0,k2)*Qx.Wg()[k2]*
                                (exp(rladt)-rladt-1.0)/rladt/rladt;
                    }
                    else if ((rla2dt-rladt)*(rla2dt-rladt)<eps)   //comparing squared difference
                    {
                        WgV_Wg_E2(k0,k2)=Qx.WgV()(k0,k2)*Qx.Wg()[k2]*
                                (1.0-exp(rladt)*(1.0-rladt))/rladt/rladt;
                    }
                    else
                    {
                        WgV_Wg_E2(k0,k2)=Qx.WgV()(k0,k2)*Qx.Wg()[k2]*
                                (1.0/rladt/rla2dt+
                                 exp(rla2dt)/rla2dt/(rla2dt-rladt)+
                                 exp(rladt)/rladt/(rladt-rla2dt));
                    }
                }
            }
        }


        // vmean(i)=A(k,i,j) g(j) A(k2,j,j2) g(j2) E(k,k)
        // =V(i,k) W(k,j) g(j) V(j,k2) W(k2,j2) g(j2) E(k,k2)
        //         ___________________ ______________
        //         WgV(k,k2)           Wg(k2)
        //___________________________________________________

        // G_(i)=V(i,n1), W(n1,k)g(k) V(k,n2) W(n2,j)       E(n1,n2)
        //                ___________________ ______________
        //                WgV(n1,n2)          W(n2)
        //___________________________________________________

        gsqr_i_=Matrix_Generators::zeros<double>(k_u,1);

        for (std::size_t i=0; i<k_u; i++)
        {
            for (std::size_t k0=0; k0<k_u; k0++)
                for (std::size_t k2=0; k2<k_u; k2++)
                    gsqr_i_[i]+=2*Qx.V()(i,k0)*WgV_Wg_E2(k0,k2);

        }

        gvar_i_=M_Matrix<double>(k_u,1);
        for (std::size_t i=0; i<k_u; i++)
            gvar_i_[i]=gsqr_i_[i]-gmean_i_[i]*gmean_i_[i];



        M_Matrix<double> E2(k_u,k_u);
        for (std::size_t k0=0; k0<k_u; k0++)
        {
            double rladt=Qx.landa()[k0]*dt;
            for (std::size_t k2=0; k2<k_u; k2++)
            {
                double rladt2=Qx.landa()[k2]*dt;
                if (sqr(rladt-rladt2)<eps)
                    E2(k0,k2)=exp(rladt2);
                else
                    E2(k0,k2)=(exp(rladt)-exp(rladt2))/(rladt-rladt2);
            };
        };
        //std::cerr<<"\nEE0\n"<<EE0;
        gtotal_ij_=Qx.V()*elemMult(Qx.WgV(),E2)*Qx.W();



        M_Matrix<double> WgV_E3=Matrix_Generators::zeros<double>(k_u,k_u);
        for (std::size_t n1=0; n1<k_u; n1++)
        {
            double x1=Qx.landa()[n1]*dt;
            double expx1=exp(x1);
            for (std::size_t n3=0; n3<k_u; n3++)
            {
                double x3=Qx.landa()[n3]*dt;
                if (sqr(x1-x3)<eps)
                {
                    for (std::size_t n2=0; n2<k_u; n2++)
                    {
                        double x2=Qx.landa()[n2]*dt;
                        if (sqr(x2-x1)<eps)
                        {
                            WgV_E3(n1,n3)+=Qx.WgV()(n1,n2)*Qx.WgV()(n2,n3)*expx1/2;
                        }
                        else
                        {
                            double expx2=exp(x2);
                            WgV_E3(n1,n3)+=Qx.WgV()(n1,n2)*Qx.WgV()(n2,n3)*
                                    (expx2+expx1*(-x2+x1-1))/(x1-x2)/(x1-x2);
                        };
                    }
                }
                else
                {
                    double expx3=exp(x3);
                    for (std::size_t n2=0; n2<k_u; n2++)
                    {
                        double x2=Qx.landa()[n2]*dt;
                        if (sqr(x2-x1)<eps)
                        {
                            WgV_E3(n1,n3)+=Qx.WgV()(n1,n2)*Qx.WgV()(n2,n3)*
                                    (expx3+expx1*(-x3+x1-1))/(x1-x3)/(x1-x3);
                        }
                        else if (sqr(x2-x3)<eps)
                        {
                            WgV_E3(n1,n3)+=Qx.WgV()(n1,n2)*Qx.WgV()(n2,n3)*
                                    (expx1+expx3*(x3-x1-1))/(x1-x3)/(x1-x3);
                        }
                        else
                        {
                            double expx2=exp(x2);
                            WgV_E3(n1,n3)+=Qx.WgV()(n1,n2)*Qx.WgV()(n2,n3)*
                                    (expx3*(x1-x2)+expx1*(x2-x3)+
                                     expx2*(x3-x1))/((x1-x2)*(x1-x3)*(x2-x3));

                        };
                    };
                };
            };
        };
        //  std::cerr<<"WgV_E3"<<WgV_E3;
        gtotal_sqr_ij_=(Qx.V()*(WgV_E3*2.0)*Qx.W());

        gmean_ij_=M_Matrix<double>(k_u,k_u);
        gvar_ij_=M_Matrix<double>(k_u,k_u);
        for (std::size_t i=0; i<k_u; i++)
            for (std::size_t j=0; j<k_u; j++)
                if (P()(i,j)>1e-9)
                {
                    gmean_ij_(i,j)=gtotal_ij_(i,j)/P()(i,j);
                    gvar_ij_(i,j)=gtotal_sqr_ij_(i,j)/P()(i,j)-gmean_ij_(i,j)*gmean_ij_(i,j);
                }
                else
                {
                    gmean_ij_(i,j)=0;
                    gvar_ij_(i,j)=0;
                }

        M_Matrix<double> u=Matrix_Generators::ones<double>(k_u,1);
        gtotal_var_ij_=elemMult(gvar_ij_,P());
        gtotal_var_i_=gtotal_var_ij_*u;
    }

};






class SingleLigandModel
{
public:
    auto Q(double x)const { return Q0_+Qa_*x;}

    auto& g(double)const { return g_;}
    auto Qx(double x) const { return Markov_Transition_rate(Q(x),g(x));}

    auto P(const Markov_Transition_rate& aQx, std::size_t n, double fs) const {return Markov_Transition_step(aQx,n,fs); }

    auto Pdt(const Markov_Transition_rate& x, std::size_t n, double fs) const {return Markov_Transition_step_double(x,n,fs); }

    std::size_t size()const {return Q0_.nrows();}

    SingleLigandModel(std::pair<M_Matrix<double>,M_Matrix<double>>&& Qs, M_Matrix<double>&& g)
        :Q0_{std::move(Qs.first)},Qa_{std::move(Qs.second)}, g_{std::move(g)}{}
private:
    M_Matrix<double> Q0_;
    M_Matrix<double> Qa_;
    M_Matrix<double> g_;

};








template <class M, class Experiment,class X>
class Markov_Model_calculations
{
    void schedule_P(X x, std::size_t nsamples)
    {
        ++map_[x][nsamples];
    }


    std::map<std::size_t, M_Matrix<double>> calc_P_map(const X& x, const std::map<std::size_t, std::size_t>& m, std::size_t min_rep) const
    {
        std::map<std::size_t, M_Matrix<double>> out;
        auto nmax=m.rbegin()->first;
        auto P=calc_P(m_,x,1,fs_);
        for (std::size_t i=1; i<nmax>>1; i*=2)
        {
            out[i]=P;
            P=P*P;
        }
        for (auto& e: m)
        {
            if (e.second>=min_rep)
            {
                out[e.first]=get_P_given_x(x,out,e.first);
            }
        }
        return out;
    }






public:
    auto Peq(const X& x)const
    {
        auto Qx=m_.Qx(x);
        return Qx.Peq();
    };

    std::size_t N(const X&)const
    {
        return N_;
    }

    auto g(const X& x)const { return m_.g(x);}

    Markov_Transition_step const & calc_P_given_x(const X& x,const std::map<std::size_t,Markov_Transition_step>& m,std::size_t nsamples) const
    {
        auto it=m.find(nsamples);
        if (it!=m.end())
            return  it->second;
        else
        {
            auto Qx=rate_map_.at(x);


        }




    }


    Markov_Transition_step const & get_P_given_x(const X& x,const std::map<std::size_t,Markov_Transition_step>& m,std::size_t nsamples) const
    {
        auto it=m.find(nsamples);
        if (it!=m.end())
            return  it->second;
        else {
            std::size_t n0=base2_floor(nsamples);
            if (n0==nsamples)
            {
                n0>>=1;
                auto out0=get_P_given_x(x,m,n0);
                return out0*out0;

            }
            else
            {
                std::size_t n1=nsamples-n0;

                auto out0=get_P_given_x(x,m,n0);
                auto out1=get_P_given_x(x,m,n1);
                return out0*out1;
            }
        }
    }

    Markov_Transition_rate calc_Qx(X x)const
    {
        return m_.Qx(x);
    }

    Markov_Transition_step calc_P(const Markov_Transition_rate& Qx,X /*x*/,std::size_t nsamples)const
    {
        return m_.P(Qx,nsamples,fs());
    }

    Markov_Transition_step calc_P(X x,std::size_t nsamples)const
    {

        return m_.P(calc_Qx(x),nsamples,fs());
    }


    Markov_Transition_step const& get_P(X x,std::size_t nsamples) const
    {
        if ((x==current_x_)&&(nsamples==current_nsamples_))
            return current_P_;

        else
        {
            auto it=step_map_.find(x);
            if (it==step_map_.end())
            {
                current_x_=x;
                current_nsamples_=nsamples;
                current_P_= calc_P(x,nsamples);
                return current_P_;
            }
            else
            {
                auto it2=it->second.find(nsamples);
                if (it2!=it->second.end())
                {
                    return it2->second;
                }
                else
                {
                    Markov_Transition_rate Qx=rate_map_.at(x);
                    current_x_=x;
                    current_nsamples_=nsamples;
                    current_P_= calc_P(Qx,x,nsamples);
                    return current_P_;

                }

            }
        }
    }



    double fs()const {return fs_;}

    Markov_Model_calculations(const M& m, std::size_t N, const Experiment& e)
        :    m_{m},N_{N},fs_{e.frequency_of_sampling()},e_{e},rate_map_{}, map_{},step_map_{}{}

private:
    const M& m_;
    size_t N_;
    double fs_;

    const Experiment& e_;


    mutable X current_x_;
    mutable std::size_t current_nsamples_;
    mutable Markov_Transition_step current_P_;

    std::map<X, Markov_Transition_rate> rate_map_;

    std::map<X,std::map<std::size_t, std::size_t>> map_;

    std::map<X,std::map<std::size_t, Markov_Transition_step>> step_map_;



};




#endif // QMODEL_H
