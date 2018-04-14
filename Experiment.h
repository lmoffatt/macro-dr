#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <vector>
#include "Markov.h"

namespace experiment {


template<typename X=double, typename Y=double>
struct point{
    double dt;
    X x;
    Y y;
    template<typename Z>
    point(const point<X,Z>& p, Y new_y):
        x{p.x},y{new_y}{}

};


template<typename X=double, typename Y=double>
class Experiment
{
public:
    auto pre_begin() {return steps_.begin();}
    auto begin(){return traces_.begin();}
    auto end()const {return traces_.cend();}
    

    class step{
        Experiment& e_;
        std::size_t myIndex_;
        std::size_t index_of_start_point_;
        double dt_;
        Y y_; //invariant: average of point::y
        void calc()
        {
            if (++begin()==end())
            {
                dt_=(*begin()).dt;
                y_=(*begin()).y;
            }
            else
            {
                dt_=0;
                y_={};
                for (auto it=begin(); it!=end(); ++it)
                {
                    dt_+=(*it).dt;
                    y_+=(*it).y * (*it).dt;
                }
                y_/=dt_;
            }
        }


    public:
        double dt()const {return dt_;}
        Y y()const {return y_;}
        auto begin(){return e_.points_.begin()+index_of_start_point_;}
        auto begin()const {return e_.points_.begin()+index_of_start_point_;}
        auto end()const
        {
            if (myIndex_+1<e_.steps_.size())
                return e_.points_.begin()+e_.steps_[myIndex_+1].index_of_start_point_;
            else
                return e_.points_.end();
        }
        step(trace& tr,std::size_t myIndex, std::size_t index_of_start_point)
            :e_{tr}, myIndex_{myIndex},index_of_start_point_{index_of_start_point},dt_{},y_{}
        {
            calc();
        }

    };

    class trace
    {

        Experiment& e_;
        std::size_t index_of_start_step_;
        std::size_t index_of_end_step_;
    public:
        auto begin(){return e_.steps_.begin()+index_of_start_step_;}
        auto end()const{return e_.steps_.begin()+index_of_end_step_;}
        trace(Experiment& e, std::size_t index_of_start_step, std::size_t index_of_end_step )
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


    Experiment(std::vector<point>&& points,const std::vector<std::size_t>& start_points, const std::vector<std::pair<std::size_t, std::size_t>>& step_start_trace)
        : points_{std::move(points)}, steps_{},traces_{}

    {
        for (std::size_t i=0; i<start_points.size(); ++i)
            steps_.emplace_back(step(*this,i,start_points[i]));
        for (auto e:step_start_trace)
            traces_.emplace_back(trace(*this,e.first,e.second));
    }

    Experiment(std::vector<point>&& points,const std::vector<std::size_t>& start_points)
        : points_{std::move(points)}, steps_{},traces_{}

    {
        for (std::size_t i=0; i<start_points.size(); ++i)
            steps_.emplace_back(step(*this,i,start_points[i]));
        extract_traces_from_Nan();
    }


    Experiment(std::vector<point>&& points)
        : points_{std::move(points)}, steps_{},traces_{}

    {
        for (std::size_t i=0; i<points_.size(); ++i)
            steps_.emplace_back(step(*this,i,i));
        extract_traces_from_Nan();
    }


    Experiment limit_dt(double min_dt)const
    {
        auto points=points_;
        return Experiment(std::move(points),getIntervalsForMinDt(points_,min_dt));
    }

    template<class Z>
    static Experiment set_Points(std::vector<point>& newPoints,const Experiment<X,Z>& old)
    {
         return Experiment(newPoints, old.get_step_start_points(),old.get_step_start_trace());
    }


static    std::vector<std::size_t> getIntervalsForMinDt(const std::vector<point>& points, double min_dt)
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
    friend class step;
    friend class point;
    std::vector<point> points_;
    std::vector<step> steps_;
    std::vector<trace> traces_;
}; // namespace experiment











}

#endif // EXPERIMENT_H
