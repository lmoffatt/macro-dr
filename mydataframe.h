#ifndef MYDATAFRAME_H
#define MYDATAFRAME_H

#include "mytypetraits.h"
#include "mySerializer.h"
#include "myoptional.h"
#include "Matrix.h"
#include <vector>
#include <string>
#include <map>

namespace io {

template <typename T>
std::vector<T> create_empty_vector(const std::vector<T>& )
{
    return std::vector<T>();
}


template <typename ...Ts>
class myDataFrame
{

public:
    struct col
    {
        template<typename T>
        col(std::pair<std::string,C<T>> x, std::size_t nrows)
            :title{x.first}, data{std::vector<T>(nrows)}{}
        template<typename T>
        col(std::pair<std::string,C<T>> x)
            :title{x.first}, data{std::vector<T>()}{}
        template<typename T>
        col(const std::string& name,std::vector<T>&& x)
            :title{name}, data{std::move(x)}{}
        col()=default;
        std::string title;
        std::variant<std::vector<Ts>...> data;
        typedef  std::variant<std::vector<Ts>...> data_type;

        friend bool compatible(const col& one, const col& two)
        {
            return ((one.title==two.title) && (one.data.index()==two.data.index()));
        }
        col emptyCopy()const
        {
            col out;
            out.title=title;
             out.data=std::visit([](auto& e){ return data_type(create_empty_vector(e));},data);
            return out;
        }

        col& concatenate(const col& other)
        {
            assert((title==other.title)&&(data.index()==other.data.index()));
            std::visit([&other](auto& me){
                auto& next=std::get<std::decay_t<decltype(me)>>(other.data);
                me.insert(me.begin(),next.begin(),next.end());},data);
            return *this;
        }

    };

    struct row
    {
        std::vector<std::variant<Ts...>> data;
    };


    template<typename... Ks>
    bool same_types_i(const std::tuple<Ks...>&,std::index_sequence<>)
    {
        return true;
    }

    template<std::size_t I,std::size_t... Is,typename... Ks>
    bool same_types_i(std::tuple<Ks...> t, std::index_sequence<I,Is...>)
    {

        if(!std::holds_alternative<std::vector<std::decay_t<std::tuple_element_t<I,std::tuple<Ks...>>>>>(data_[I].data))
            return false;
        else return same_types_i(t, std::index_sequence<Is...>());
    }



    template<typename... Ks>
    bool same_types(Ks&&... row_of_data)
    {
        constexpr auto N=sizeof...(Ks);
        if (N!=data_.size())
            return false;
        else
            return same_types_i(std::forward_as_tuple(std::forward<Ks>(row_of_data)...), std::make_index_sequence<N>());
    }

    template<typename...K0>
    bool push_back_i(std::tuple<K0...> ,Cs<K0...>,Cs<>,std::index_sequence<>){
        return true;
    }

    template<std::size_t I,std::size_t... Is,typename...K0, typename K, typename ...K1>
    bool push_back_i(std::tuple<K0...,K,K1...> t,Cs<K0...>,Cs<K,K1...>,std::index_sequence<I,Is...>)
    {
        std::get<std::vector<std::decay_t<K>>>(data_[I].data).push_back(std::move(std::get<I>(t)));
        return push_back_i(t,Cs<K0...,K>{},Cs<K1...>{},std::index_sequence<Is...>{});
    }



    template<typename... Ks>
    bool push_back_t(const std::tuple<Ks...>& data)
    {
        return std::apply([this](auto... x){return push_back(x...);},data);
    }

    template<typename... Ks>
    bool push_back(Ks&&... row_of_data)
    {
        assert (same_types(row_of_data...));
            return push_back_i(std::forward_as_tuple<Ks...>(std::forward<Ks>(row_of_data)...),Cs<>{},Cs<Ks...>{},std::index_sequence_for<Ks...>{});
    }




    template<typename... Ks>
    std::ostream& write_back(std::ostream& os,Ks&&... row_of_data)
    {
        assert(same_types(row_of_data...));
        bool sucess= ((io::write_on_element(os,row_of_data))&&...&&true);
        return os;
    }


    std::size_t ncols()const { return data_.size();}

    std::size_t nrows() const {
        if (ncols()>0)
            return std::visit([](auto const& arg){return arg.size();},data_[0].data);
        else return 0;
    }

    template<class T>
    bool insert_column(std::string&& title, C<T>)
    {
        if (nrows()>0)
            return  false;
        else {
            data_.push_back(col(std::pair(std::move(title),C<T>{})));
            return true;
        }
    }

    myDataFrame& concatenate(const myDataFrame& other )
    {
        assert(compatible(*this,other));
        for (std::size_t j=0; j<ncols(); ++j)
            data_[j].concatenate(other.data_[j]);
        return *this;

    }

    static myDataFrame concatenate(const std::vector<myDataFrame>& v )
    {
        if (v.size()==0) return myDataFrame();
        else
        {
            myDataFrame out=v[0].emptyCopy();
            auto n=rows_size(v);
            out.reserve(n);
            for (auto &e: v)
                out.concatenate(e);
            return out;
        }

    }

    static myDataFrame merge(myDataFrame&& one, myDataFrame&& other)
    {
        assert(one.nrows()==other.nrows());
        std::vector<col> out(std::move(one.data_));
        out.insert(out.end(),std::make_move_iterator(other.data_.begin()), std::make_move_iterator(other.data_.end()));
        return myDataFrame(std::move(out));
    }




    static std::size_t rows_size(const std::vector<myDataFrame>& v){ std::size_t n=0; for (auto& e:v) n+=e.nrows(); return n;}

    static myDataFrame consolidate(const std::vector<myDataFrame>& v, const std::string& col_title)
    {
        std::vector<std::size_t> index(rows_size(v));
        std::size_t n=0;
        for (std::size_t i=0; i<v.size(); ++i)
        {
            for (std::size_t j=0; j<v[i].nrows(); ++j) {index[n]=i; ++n;}
        }
       return merge(myDataFrame(std::pair(col_title,index)),concatenate(v));
    }


    template<typename T>
    myOptional_t<T> get(const std::string& id, std::size_t j)const
    {
        auto it=map_.find(id);
        if (it==map_.end()) return{false,"field not found"};
        else
            return std::visit([j](auto const& a)
            {return a[j];}, data_[it->second].data);
    }

    std::ostream& write_title(std::ostream& os)const
    {
        for (std::size_t i=0; i<ncols(); ++i)
            if (i+1<ncols())
                os<<data_[i].title<<io::separator{};
            else
                os<<data_[i].title;
        os<<io::end_of_line{};
        return os;
    }

    std::ostream& write_row(std::ostream& os, std::size_t j)const
    {
        for (std::size_t i=0; i<ncols(); ++i)
        {
            if (i+1<ncols())
                std::visit([j, &os](auto const & a){io::write_on_element(os,a[j]);os<<io::separator{};},data_[i].data);
            else
                std::visit([j, &os](auto const & a){io::write_on_element(os,a[j]);},data_[i].data);

        }
        os<<io::end_of_line{};
        return os;

    }

    std::ostream& write(std::ostream& os)const
    {
        write_title(os);
        for (std::size_t j=0; j<nrows(); ++j)
            write_row(os,j);
        return os;
    }


    std::ostream& operator<<(std::ostream& os)const
    {
        return write(os);
    }



    struct lambda
    {
        std::stringstream& ss;
        template<typename T>
        auto operator()(T& a){
            io::read_one_element_on_vector(ss,a);}
        lambda(std::stringstream& s):ss{s}{}
    };
    std::istream& read(std::istream& is)
    {
        std::vector<std::string> colnames;
        read_on_vector(is,colnames);
        std::vector<col> d(colnames.size());
        for (std::size_t i=0; i<d.size(); ++i )
            d[i].title=std::move(colnames[i]);

        std::string line;
        safeGetline(is, line);
        while (!line.empty())
        {
            std::stringstream ss(line);
            for (std::size_t i=0; i<d.size(); ++i )
            {
                std::visit(lambda(ss),d[i].data);
                if (i+1<d.size()) ss>>io::separator{};
            }safeGetline(is,line);
        }
        *this=myDataFrame(std::move(d));
        return is;

    }

    std::istream& operator>>(std::istream& is)
    {
        return read(is);

    }





    static std::map<std::string, std::size_t> getMap(const std::vector<col>& d)
    {
        std::map<std::string, std::size_t> m;
        for (std::size_t i=0; i<d.size(); ++i)
            m[d[i].title]=i;
        return m;
    }

    template<typename... Ks>
    myDataFrame(std::pair<std::string,C<Ks>>&&... titles):
        data_{{col(titles)...}},map_{getMap(data_)}{}

    template<typename... Ks>
    myDataFrame(std::pair<std::string,std::vector<Ks>>&&... titles):
        data_{{col(std::move(titles.first),std::move(titles.second))...}},map_{getMap(data_)}{}



    template<typename... Ks>
    myDataFrame(std::pair<std::string,C<Ks>>&&... titles, std::size_t nrows):
        data_{{col(titles,nrows)...}},map_{getMap(data_)}{}


    myDataFrame(std::vector<col>&& data):data_{std::move(data)}, map_{getMap(data_)}{}
    myDataFrame()=default;

    friend bool compatible(const myDataFrame& one, const myDataFrame& two)
    {
        if (one.ncols()!=two.ncols()) return false;
        else for (std::size_t j=0; j<one.ncols(); ++j)
            if (!compatible(one.data_[j], two.data_[j])) return false;
        return true;
    }

  myDataFrame emptyCopy()const
  {
     std::vector<col> d;
     for (auto& e:data_)
         d.push_back(e.emptyCopy());
     return myDataFrame(std::move(d));
  }


    constexpr static auto className=my_static_string("DataFrame")+my_trait<Cs<Ts...>>::className;


private:
    std::vector<col> data_;
    std::map<std::string, std::size_t> map_;

    void reserve(std::size_t n)
    {
        for (auto& e:data_)
            std::visit([n](auto& v){v.reserve(n);},e.data);
    }

};

template <typename ...Ts>
std::ostream& write(std::ostream& os, const myDataFrame<Ts...>& d)
{
    return d.write(os);
}

template <typename ...Ts>
std::istream& read(std::istream& is, myDataFrame<Ts...>& d)
{
    return d.read(is);
}

template <typename ...Ts>
std::ostream& operator<<(std::ostream& os, const myDataFrame<Ts...>& d)
{
    return d.write(os);
}

template <typename ...Ts>
std::istream& operator>>(std::istream& is, myDataFrame<Ts...>& d)
{
    return d.read(is);
}



} // namespace io
#endif // MYDATAFRAME_H
