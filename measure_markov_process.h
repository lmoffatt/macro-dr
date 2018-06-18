#ifndef MEASURE_MARKOV_PROCESS_H
#define MEASURE_MARKOV_PROCESS_H




#include "Matrix.h"
#include "myDistributions.h"
#include "Experiment.h"




#include <type_traits>

namespace markov {
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


template<class Model>
double getCurrent(const markov_process<std::size_t>& mp,const Model& m)
{
    return mp.N()*m.g();
}

namespace old {


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

}

template <typename X>
struct measure_algorithm_state
{
    X x;
    std::size_t nsubsamples;
    std::size_t sum_samples;
};


template<class Model,  class measure_algorithm_state,typename X>
void actualize_mp(measure_algorithm_state& current,markov_process<std::size_t>& mp, Model& m,X new_x,std::size_t actual_nsubsamples)
{
if ((current.nsubsamples!=actual_nsubsamples)||(current.x!=new_x))
{
    current.x=new_x;
    current.nsubsamples=actual_nsubsamples;
    auto P=m.get_P(current.x,current.nsubsamples).P();
    mp.set_P(std::move(P));
}
}


template<class Model,  typename X>
auto init_mp(std::mt19937_64& mt, Model& m, X x, std::size_t nsamples)
{
    measure_algorithm_state<X> current;
    auto Peq=m.Peq(x);
    std::size_t N=m.N(x);
    auto M=multinomial_distribution<M_Matrix,std::size_t>(N,Peq);
    auto Nini=M(mt);
    current.x=x;
    current.sum_samples=nsamples;
    current.sum_samples=0;
    auto P=m.get_P(x,nsamples);
    return std::make_tuple(markov_process<std::size_t>(Nini, P.P()),current);
}



template<class F,class Model, template <typename, typename>class Point, class measure_algorithm_state,typename X, typename Y>
auto measure_point(measure_algorithm_state& current,const F& f,std::mt19937_64& mt,markov_process<std::size_t>& mp, Model& m, const Point<X,Y>& in_point,  std::size_t desired_nsubsamples)
{

    auto n_samples=in_point.nsamples();
    auto new_x=in_point.x();
    auto actual_nsubsamples=std::min(desired_nsubsamples,n_samples);
    auto n_times=n_samples/actual_nsubsamples;
    auto remainder_samples=n_samples-n_times*actual_nsubsamples;
    actualize_mp(current,mp,m,new_x,actual_nsubsamples);
    mp.set_N(mp(mt));
    auto y=f(mp, m,new_x);
    for (auto i=1u; i<n_times; ++i)
    {
        mp.set_N(mp(mt));
        y+=f(mp,m,new_x)*actual_nsubsamples;
    }
    if (remainder_samples>0)
    {
        actualize_mp(current,mp,m,current.x,remainder_samples);
        mp.set_N(mp(mt));
        y+=f(mp, m,new_x)*remainder_samples;
    }
    y=y/n_samples;
    current.nsubsamples+=n_samples;
    return Point<X,decltype (y)> (in_point,y);
}






template<class F,class Model, class measure_algorithm_state, class Step, class Point>
void measure_step(std::vector<Point>& out,measure_algorithm_state& current,const F& f,std::mt19937_64& mt,markov_process<std::size_t>& mp, const Step& step,   Model& m,std::size_t n_substeps, double /*fs*/)
{
    auto nsamples=step.nsamples();
    auto desired_subsamples=std::max(nsamples/n_substeps,1ul);
    current.sum_samples=0;
    for (auto it=step.begin(); it!=step.end(); ++it)
    {
        auto y=measure_point(current,f,mt,mp,m,*it,desired_subsamples);
        out.push_back(y);
    }
}
template<class F,class Model, class measure_algorithm_state, class Step, class Point>
void skip_step(std::vector<Point>& out,measure_algorithm_state& current,const F& /*f*/,std::mt19937_64& mt,markov_process<std::size_t>& mp, const Step& step,
               Model& m, double/* fs*/)
{

    current.sum_samples=0;
    for (auto it=step.begin(); it!=step.end(); ++it)
    {
        actualize_mp(current,mp,m,it->x(),it->nsamples());
        mp.set_N(mp(mt));
        out.push_back(Point(*it,std::numeric_limits<double>::quiet_NaN()));
    }
}



template<class F, class Model,class measure_algorithm_state, class trace, class Point>
void meansure_trace(std::vector<Point>& out,measure_algorithm_state& current,const F& f,std::mt19937_64& mt,markov_process<std::size_t>& mp, Model& m,const trace& t, std::size_t n_substeps,double fs)
{

    for (auto it=t.begin();it!=t.end(); ++it)
    {
        measure_step(out,current,f,mt,mp,*it,m,n_substeps,fs);
    }
    skip_step(out,current,f,mt,mp,*t.end(),m,fs);
}


template<class F,class Model,template<class, class > class Experiment, class Point, class measure>
auto measure_experiment(const F& f, std::mt19937_64& mt, Model& m,const Experiment<Point,measure>& e, std::size_t n_substeps)
{
    double fs=e.frequency_of_sampling();
    auto first_point=*e.begin_begin_begin();
    auto x=first_point.x();
    std::size_t n_sub_samples=first_point.nsamples()/n_substeps;
    auto [mp, current]=init_mp(mt,m,x,n_sub_samples);

    auto y=measure_point(current,f,mt,mp,m,first_point,n_sub_samples);

    std::vector<decltype (y)> points;

    for (auto it=e.begin();it!=e.end(); ++it)
    {
       meansure_trace(points,current,f,mt,mp,m,*it,n_substeps,fs);
    }
    return Experiment<Point, measure>(std::move(points),fs);
 }





} // namespace markov





#endif // MEASURE_MARKOV_PROCESS_H
