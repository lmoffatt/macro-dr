#ifndef MARKOV_H
#define MARKOV_H



#include "Matrix.h"
#include "myDistributions.h"
#include "Experiment.h"




#include <type_traits>


template<class T,class Operation,class State,class Model, class IntervalIterator, template<class>class Vector>
auto apply(std::mt19937_64& mt,const Operation& op,const State& s,const Model& m, IntervalIterator first, const IntervalIterator& last, Vector<std::invoke_result_t<Operation>>& out)
{
    State r=s;
    for (std::size_t i=0; first!=last; ++first, ++i)
    {

        std::tie(out[i],r)=op(mt,r,m,*first);
    }
    return out;
}



double getCurrent(const markov_process<std::size_t>& mp,const Model& m)
{
    return mp.N()*m.g();
}


template<class F,class State,class Model, class Input>
auto measure_point(const F& f,std::mt19937_64& mt,markov_process<std::size_t>& mp,State& s,const Model& m, const Input& x, double ddt, std::size_t n)
{
    bool change=false;
    if(s.x!=x)
    {
        change=true;
        s.Q=m.Q(x);
        s.x=x;
    }
    if (s.dt!=ddt)
    {
        change=true;
        s.dt=ddt;
    }
    if (change)
    {
        mp.set_P(expm(s.Q*ddt));
    }

    mp.set_N(mp(mt));
    auto y=f(mp, m);
    for (auto i=1u; i<n; ++i)
    {
        mp.set_N(mp(mt));
        y+=f(mp,m);
    }
    y/=n;
    return y;
}


template<class State,class Model, class Input>
void  skip_point(std::mt19937_64& mt,markov_process<std::size_t>& mp,State& s,const Model& m, const Input& x, double dt)
{
    bool change=false;
    if(s.x!=x)
    {
        change=true;
        s.Q=m.Q(x);
        s.x=x;
    }
    if (s.dt!=dt)
    {
        change=true;
        s.dt=dt;
    }
    if (change)
    {
        mp.set_P(expm(s.Q*dt));
    }

    mp.set_N(mp(mt));
 }





template<class F,class State,class Model, class PointIterator, class Point>
void measure_step(const F& f,std::mt19937_64& mt,std::vector<Point>& p, markov_process<std::size_t>& mp,State& s,const Model& m, PointIterator first, PointIterator last, double dt_max)
{
    auto x=first->x();
    double dt=first->dt();
    std::size_t n=std::round(dt/dt_max);
    double ddt=dt/n;

    auto y=measure_point(f,mt,mp,s,m,x,ddt,n);
    p.emplace_back(point{dt,x,std::move(y)});
    ++first;
    if (first!=last)
    {
        for (;first!=last; ++first)
        {

            x=first->x();
            double dt=first->dt();
            std::size_t n=std::round(dt/dt_max);
            ddt=dt/n;
            auto y=meanCurrent_point(mt,mp,s,m,x,ddt,n);
            p.emplace_back(point{dt,x,y});
        }
    }
}

template<class State,class Model, class PointIterator, class Point>
void skip_step(std::mt19937_64& mt,std::vector<Point>& p, markov_process<std::size_t>& mp,State& s,const Model& m, PointIterator first, PointIterator last)
{
    auto x=first->x();
    double dt=first->dt();
    skip_point(mt,mp,s,m,x,ddt);
    p.emplace_back(Point{dt,x,std::numeric_limits<double>::quiet_NaN()});
    ++first;
    if (first!=last)
    {
        for (;first!=last; ++first)
        {
            x=first->x();
            double dt=first->dt();
            skip_point(mt,mp,s,m,x,dt);
            p.emplace_back(Point{dt,x,std::numeric_limits<double>::quiet_NaN()});
        }
    }
}



template<class F, class State,class Model, class StepIterator, class Point>
void meansure_trace(const F& f,std::mt19937_64& mt,std::vector<Point>& p, markov_process<std::size_t>& mp,State& s,const Model& m, StepIterator first, StepIterator last, std::size_t n)
{
    for (;first!=last; ++first)
    {
        double dt=first->dt();
        double dtmin=dt/n;
        measure_step(f,mt,p,mp,s,m,first->begin(), first->end(),dtmin);
    }
}


template<class F,class Model, class Experiment>
auto measure_experiment(const F& f, std::mt19937_64& mt,const Model& m,Experiment e, std::size_t n)
{
    auto it_step=e.pre_begin();

    auto x=it_step->begin()->x;

    struct state{
        decltype(m.Q) Q;
        double dt;
        decltype(x)  x;
    };

    state s;
    double dt=s.x=it_step->begin()->x;
    s.Q=m.Q(s.x);

    auto Peq=m.Peq(s.x);
    std::size_t N=m.N(s.x);

    auto Nini=multinomial_distribution<M_Matrix,std::size_t>(N,Peq)(mt);
    auto P=expm(s.Q*dt);
    markov_process<std::size_t> mp(Nini, P);

    typedef std::invoke_result_t<F,markov_process<std::size_t>,Model> R;

    typedef std::vector<experiment::point<decltype(s.x),R>> V;

    V p;
    for (;it_step!=e.begin()->begin(); ++it_step)
    {
       skip_step(mt,p,mp,s,m,it_step->begin(), it_step->end());
    }
    for (auto it_trace=e.begin(); it_trace!=e.end(); ++it_trace)
    {
       meansure_trace(f,mt,p,mp,s,m,it_trace->begin(), it_trace->end(),n);
       auto it_next=it_trace+1;
       if (it_next!=e.end())
       for(it_step=it_trace->end(); it_step!=it_next->begin(); it_step++)
       {
           skip_step(mt,p,mp,s,m,it_step->begin(), it_step->end());
       }

   }
    return p;
}














#endif // MARKOV_H
