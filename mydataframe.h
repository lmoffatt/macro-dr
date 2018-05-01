#ifndef MYDATAFRAME_H
#define MYDATAFRAME_H

#include "mytypetraits.h"
#include "mySerializer.h"
#include "Matrix.h"
#include <vector>
#include <string>
#include <map>

namespace io {


template <typename ...Ts>
class myDataFrame
{

public:
   struct col
   {
       template<typename T>
       col(std::pair<std::string,C<T>> x)
           :title{x.first}, data{std::vector<T>()}{}
       col()=default;
       std::string title;
       std::variant<std::vector<Ts>...> data;
       typedef  std::variant<std::vector<Ts>...> data_type;
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
   void push_back_i(std::tuple<K0...> ,Cs<K0...>,Cs<>,std::index_sequence<>){}

   template<std::size_t I,std::size_t... Is,typename...K0, typename K, typename ...K1>
   void push_back_i(std::tuple<K0...,K,K1...> t,Cs<K0...>,Cs<K,K1...>,std::index_sequence<I,Is...>)
   {
       std::get<std::vector<K>>(data_[I].data).push_back(std::move(std::get<I>(t)));
       push_back_i(t,Cs<K0...,K>{},Cs<K1...>{},std::index_sequence<Is...>{});
   }


   template<typename... Ks>
   bool push_back(Ks&&... row_of_data)
   {
       if (!same_types(row_of_data...))
           return false;
       else
           push_back_i(std::forward_as_tuple<Ks...>(std::forward<Ks>(row_of_data)...),Cs<>{},Cs<Ks...>{},std::index_sequence_for<Ks...>{});
   }


   std::size_t ncols()const { return data_.size();}

   std::size_t nrows() const {
       if (ncols()>0)
       return std::visit([](auto const& arg){return arg.size();},data_[0].data);
   else return 0;
   }


   template<typename T>
   std::optional<T> get(const std::string& id, std::size_t j)const
   {
       auto it=map_.find(id);
       if (it==map_.end()) return{};
       else
       return std::visit([j](auto const& a)
       {return a[j];}, data_[it->second].data);
   }


   std::ostream& write(std::ostream& os)const
   {
      for (std::size_t i=0; i<ncols(); ++i)
       io::write_on_element(os,data_[i].title);
      os<<io::end_of_line{};
      for (std::size_t j=0; j<nrows(); ++j)
      {
          for (std::size_t i=0; i<ncols(); ++i)
          {
              std::visit([j, &os](auto const & a){io::write_on_element(os,a[j]);},data_[i].data);
          }
          os<<io::end_of_line{};
      }
      return os;
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
            std::visit(lambda(ss),d[i].data);
          safeGetline(is,line);
      }
      *this=myDataFrame(std::move(d));
      return is;

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

   myDataFrame(std::vector<col>&& data):data_{std::move(data)}, map_{getMap(data_)}{}
  myDataFrame()=default;
private:
   std::vector<col> data_;
   std::map<std::string, std::size_t> map_;

};

template <typename ...Ts>
std::ostream& write(std::ostream& os, const myDataFrame<Ts...>& d)
{
    return d.write(os);
}

template <typename ...Ts>
std::istream& read(std::ostream& is, myDataFrame<Ts...>& d)
{
    return d.read(is);
}




} // namespace io
#endif // MYDATAFRAME_H
