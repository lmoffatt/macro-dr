#ifndef QMODEL_H
#define QMODEL_H

#include "Matrix.h"
#include "myfields.h"
#include <vector>
#include <map>
#include <set>
#include <cmath>






class Allosteric_Model
{
public:
    typedef  Allosteric_Model self_type;
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



private:
    std::vector<std::vector<std::string>> conformer_;
    std::vector<std::map<std::size_t,transitions>> transitions_;
    std::set<std::string> paramNames_;
    std::vector<std::string> conductances_;
    model_definition d_;
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


    static auto get_constructor_fields()
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


    Allosteric_Model()=default;

    Allosteric_Model(const std::vector<std::string>& conformational_changes,
                     const std::map<std::pair<std::string,bool>,std::string> conformational_changes_names,
                     const std::set<std::size_t>& agonist_changes,
                     const std::set<std::size_t>& conductance_changes,
                     const std::map<std::size_t, std::string>& conductance_names,
                     const std::multimap<std::size_t, std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>>& conformational_inter_unit_cell)
        :d_{conformational_changes,conformational_changes_names,agonist_changes,conductance_changes,conductance_names,conformational_inter_unit_cell}
    {
        paramNames_=getParameterNamesFrom(conformational_changes_names,conductance_names,conformational_inter_unit_cell);
        std::cout<<paramNames_;
        auto p=periodicity(conformational_changes);

        auto conformational_interactions=fill_conformational_interactions(conformational_inter_unit_cell,conformational_changes.size(),p);

        std::map<std::vector<std::string>,std::size_t> state_to_conformer;

        std::tie(conformer_, state_to_conformer)=getConformers(conformational_changes,p);

        transitions_=getTransitions(conformational_changes,conformer_,state_to_conformer,agonist_changes,conformational_changes_names,conformational_interactions);
        conductances_=get_g_names(conformer_,conductance_changes,conductance_names);

    }


    template<class Parameters>
    auto Qs(const Parameters& p)
    {

        M_Matrix<double> Q0(conformer_.size(),conformer_.size());
        M_Matrix<double> Qa(conformer_.size(),conformer_.size());
        for (std::size_t i=0; i<Q0.nrows(); ++i)
        {

            for (auto it=transitions_[i].begin(); it!=transitions_[i].end(); ++it)
            {
                std::size_t j=it->first;
                if (it->second.agonist)
                    Qa(i,j)=rate(it->second,p);
                else
                    Q0(i,j)=rate(it->second,p);
            }
        }
        return std::pair(Q0,Qa);
    }


    template<class Parameters>
    auto g(const Parameters& p)
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

#endif // QMODEL_H
