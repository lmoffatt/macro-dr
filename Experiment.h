#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <vector>
#include "Markov.h"
#include "mydataframe.h"
namespace experiment {


struct experiment_tag;

template<typename X, typename Y>
struct myPulsesPoints
{
    typedef myPulsesPoints<X,Y> self_type;
    double t;
    std::size_t ns;
    double dt;
    X x;
    Y y;

    std::ostream& operator<<(std::ostream& os)
    {
        return   os<<t<<io::separator{}<<ns<<io::separator{}<<dt<<io::separator{}<<x<<io::separator{}<<y<<io::separator{};
    }

    std::istream& operator>>(std::istream& is)
    {
        return   is>>t>>io::separator{}>>ns>>io::separator{}>>dt>>io::separator{}>>x>>io::separator{}>>y>>io::separator{};
    }


};





template<typename X=double, typename Y=double>
struct point{

    constexpr static auto const name=my_static_string("point_of_")+my_trait<X>::name+ my_static_string("_")+my_trait<Y>::name;
    double t_;
    std::size_t nsamples_;
    X x_;
    Y y_;

    double t()const {return t_;}
    std::size_t nsamples()const {return nsamples_;}
    X x()const {return x_;}
    Y y() const {return y_;}
    std::ostream& operator<<(std::ostream& os)
    {
        return   os<<t<<io::separator{}<<t()<<io::separator{}<<nsamples()<<io::separator{}<<x()<<io::separator{}<<y()<<io::separator{};
    }

    std::istream& operator>>(std::istream& is)
    {
        return   is>>t>>io::separator{}>>t()>>io::separator{}>>nsamples()>>io::separator{}>>x()>>io::separator{}>>y()>>io::separator{};
    }



    point(double _t, std::size_t _nsamples,X _x,Y _y):t_{_t},nsamples_{_nsamples}, x_{_x},y_{_y}{}

    template<typename Z>
    point(const point<X,Z>& p, Y new_y):t_{p.t()},nsamples_{p.nsamples()},
        x_{p.x()},y_{new_y}{}

};

template<class...>class basic_Experiment;
using namespace std::string_view_literals;


template<template<typename, typename>class point, typename X, typename Y>
class basic_Experiment<point<X,Y>>
{
public:


    typedef basic_Experiment<point<X,Y>>  self_type;
    typedef point<X,Y> point_type;
    auto cbegin_begin_begin()const {return points_.cbegin();}
    auto cbegin_begin()const {return steps_.cbegin();}
    auto cbegin() const {return traces_.cbegin();}

    auto begin_begin_begin()const {return points_.begin();}
    auto begin_begin()const {return steps_.begin();}
    auto begin() const {return traces_.begin();}

    auto begin_begin_begin() {return points_.begin();}
    auto begin_begin(){return steps_.begin();}

    auto begin(){return traces_.begin();}


    auto end()const {return traces_.end();}
    auto end_end()const {return steps_.end();}
    auto end_end_end()const {return points_.end();}
    auto cend()const {return traces_.cend();}
    auto cend_end()const {return steps_.cend();}
    auto cend_end_end()const {return points_.cend();}
    class trace;
    class step{
        basic_Experiment& e_;
        std::size_t myIndex_;
        std::size_t index_of_start_point_;
        std::size_t nsamples_;
        Y y_; //invariant: average of point::y


    public:
        static constexpr auto name=my_static_string("Experiment_on_")+point<X,Y>::name;

        void calc()
        {
            if (++begin()==end())
            {
                nsamples_=(*begin()).nsamples();
                y_=(*begin()).y();
            }
            else
            {
                nsamples_=0;
                y_={};
                for (auto it=begin(); it!=end(); ++it)
                {
                    nsamples_+=(*it).nsamples();
                    y_+=(*it).y() * (*it).nsamples();
                }
                y_/=nsamples_;
            }
        }
        auto nsamples()const {return nsamples_;}
        Y y()const {return y_;}
        auto begin(){return e_.begin_begin_begin()+index_of_start_point_;}
        auto begin()const {return e_.begin_begin_begin()+index_of_start_point_;}
        auto cbegin()const {return e_.cbegin_begin_begin()+index_of_start_point_;}
        auto end()const
        {
            if (myIndex_+1<e_.steps_.size())
                return e_.steps_[myIndex_+1].cbegin();
            else
                return e_.cend_end_end();
        }
        step(basic_Experiment& e,std::size_t myIndex, std::size_t index_of_start_point)
            :e_{e}, myIndex_{myIndex},index_of_start_point_{index_of_start_point},nsamples_{},y_{}
        {
        }

    };
    class trace
    {

        basic_Experiment& e_;
        std::size_t index_of_start_step_;
        std::size_t index_of_end_step_;
    public:
        auto begin(){return e_.begin_begin()+index_of_start_step_;}
        auto begin()const {return e_.begin_begin()+index_of_start_step_;}
        auto end()const{return e_.begin_begin()+index_of_end_step_;}
        trace(basic_Experiment& e, std::size_t index_of_start_step, std::size_t index_of_end_step )
            :
              e_{e},index_of_start_step_{index_of_start_step}, index_of_end_step_{index_of_end_step}
        {
        }

    };


    std::vector<std::size_t> get_step_start_points()const
    {
        std::vector<std::size_t> out(steps_.size());
        for (std::size_t i=0; i<steps_.size(); ++i)
            out[i]=steps_[i].index_of_start_point_;
        return out;
    }

    std::vector<std::pair<std::size_t, std::size_t>> get_step_start_trace()const
    {
        std::vector<std::pair<std::size_t, std::size_t>> out(traces_.size());
        for (std::size_t i=0; i<traces_.size(); ++i)
            out[i]={traces_[i].index_of_start_step_, traces_[i].index_of_end_step_};
        return out;
    }



    //----  observers of Experiment initial


    auto const & get_Points()const
    {
        return  points_;
    }

    static auto get_constructor_fields()
    {
        return std::make_tuple(
                    grammar::field(C<self_type>{},"points",&self_type::get_Points));

    }

    basic_Experiment(std::vector<point<X,Y>>&& points,const std::vector<std::size_t>& start_points, const std::vector<std::pair<std::size_t, std::size_t>>& step_start_trace)
        : points_{std::move(points)}, steps_{},traces_{}

    {
        for (std::size_t i=0; i<start_points.size(); ++i)
            steps_.emplace_back(step(*this,i,start_points[i]));
        for (auto &e: steps_) e.calc();
        for (auto e:step_start_trace)
            traces_.emplace_back(trace(*this,e.first,e.second));
    }

    basic_Experiment(std::vector<point<X,Y>>&& points,const std::vector<std::size_t>& start_points)
        : points_{std::move(points)}, steps_{},traces_{}

    {
        for (std::size_t i=0; i<start_points.size(); ++i)
            steps_.emplace_back(step(*this,i,start_points[i]));
        extract_traces_from_Nan();
    }


    basic_Experiment(std::vector<point<X,Y>>&& points, double fs)
        : frequency_of_sampling_{fs},points_{std::move(points)}, steps_{},traces_{}

    {
        for (std::size_t i=0; i<points_.size(); ++i)
            steps_.emplace_back(*this,i,i);
        for (auto &e:steps_) e.calc();

        extract_traces_from_Nan();
    }




    basic_Experiment limit_dt(double min_dt)const
    {
        auto points=points_;
        return basic_Experiment(std::move(points),getIntervalsForMinDt(points_,min_dt));
    }

    template<class Z>
    static basic_Experiment set_Points(std::vector<point<X,Y>>& newPoints,const basic_Experiment<point<X,Z>>& old)
    {
        return basic_Experiment(newPoints, old.get_step_start_points(),old.get_step_start_trace());
    }


    static    std::vector<std::size_t> getIntervalsForMinDt(const std::vector<point<X,Y>>& points, double min_dt)
    {
        double sumdt=0;
        std::vector<std::size_t> out;
        for (std::size_t i=0; i<points.size(); ++i)
        {
            sumdt+=points[i].dt;
            if (sumdt>=min_dt)
            {
                out.push_back(i);
                sumdt=0;
            }
        }
        return out;

    }

    double frequency_of_sampling()const {return  frequency_of_sampling_;}

private:
    void extract_traces_from_Nan()
    {
        std::size_t iStart=0;
        for (std::size_t i=0; i<steps_.size(); ++i)
            if (std::isnan(steps_[i].y()))
            {
                if (i>iStart)
                {
                    traces_.emplace_back(trace(*this,iStart,i));
                }
                iStart=i+1;
            }

    }
    double frequency_of_sampling_;
    std::vector<point<X,Y>> points_;
    std::vector<step> steps_;
    std::vector<trace> traces_;
}; // namespace experiment







template <typename...Ts>
basic_Experiment<point<double,double>>
DataFrame_to_Experiment(const io::myDataFrame<Ts...>& d
                        ,const std::string& colname_time,
                        const std::string& colname_nsample,
                        const std::string& colname_x,
                        const std::string& colname_y,
                        double frequency_of_sampling)
{

    std::vector<point<double,double>> out;
    for (std::size_t i=0; i<d.nrows(); ++i)
    {
        out.push_back(point<double,double>(*d.template get<double>(colname_time,i),
                                           *d.template get<std::size_t>(colname_nsample,i),
                                           *d.template get<double>(colname_x,i),
                                           *d.template get<double>(colname_y,i)));
    }
    return basic_Experiment<point<double,double>>(std::move(out),frequency_of_sampling);
}





auto Experiment_to_DataFrame(basic_Experiment<point<double,double>> e,
                        const std::string& colname_time="time",
                        const std::string& colname_nsample="nsample",
                        const std::string& colname_x="x",
                        const std::string& colname_y="y")
{
    io::myDataFrame<double,std::size_t> d(std::pair(colname_time,C<double>{}),
                              std::pair(colname_nsample,C<std::size_t>{}),
                              std::pair(colname_x,C<double>{}),
                              std::pair(colname_y,C<double>{}));

    for (auto it=e.begin_begin_begin(); it!=e.end_end_end(); ++it)
    {
        d.push_back(it->t(), it->nsamples(), it->x(), it->y());
    }
    return d;

}



}








#endif // EXPERIMENT_H

