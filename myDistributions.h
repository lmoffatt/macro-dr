#ifndef MYDISTRIBUTIONS_H
#define MYDISTRIBUTIONS_H

#include "Matrix.h"
#include "mymath.h"
# include "mytests.h"
#include "myoperators.h"
#include "mySerializer.h"
#include <random>
#include <cmath>
#include <algorithm>
#include<type_traits>



template<class Object, class Fun,class... Sizes >
class Data_Index
{
public:

    Data_Index(Cs<Object>,const std::string& fieldname,const Fun& data ,const std::pair<std::string,Sizes>&... indexes)
        :fieldname_{fieldname},indexes_{indexes...},datafun_{data}{}

    auto getFieldName()const { return fieldname_;}
    auto& dataFun()const { return datafun_;}
    auto& indexes()const {return indexes_;}

    template<class... Ints>
    auto get_Index_Size(const Object& x, Ints... is)const {
        static_assert((true&&...&&std::is_arithmetic_v<Ints>),"sizes!!");
        return std::invoke(std::get<sizeof...(Sizes)>(indexes()).second,x,is...);
    }

    std::vector<std::string> indexed_by()const {
        return std::apply([](auto&... t){ return  std::vector<std::string> {t.first...};},indexes_);
    }

    static constexpr std::size_t IndexNumber=sizeof...(Sizes);



    template<std::size_t... Is>
    bool sameIndexes(const std::vector<std::string>& ind, std::index_sequence<Is...>)const
    {
        if (ind.size()!=IndexNumber) return false;
        else
            return (bool((ind[Is])==(std::get<Is>(indexes()).first))&&...&&true);
    }

    bool isAnIndex(const std::vector<std::string>& ind)const
    {
        return std::find(ind.begin(),ind.end(),getFieldName())!=ind.end();
    }

    std::ostream& put(std::ostream& os, const Object& x,const  std::array<std::size_t, IndexNumber>& ind)const
    {

        return std::apply([&os,&x,this](auto&... i)->std::ostream&
                          {
            return os<<std::invoke(dataFun(),x,i...);
        }, ind);

    }

    std::vector<std::size_t> start(const Object& )const{
        return std::vector<std::size_t>(IndexNumber,0ul);
    }

        std::vector<std::size_t> end(const Object& x)const{
        std::vector<std::size_t> out(IndexNumber,0ul);
        out[0]=std::invoke(std::get<0>(indexes()).second,x);
        return out;
    }


    template<std::size_t I>
    std::size_t size(const Object& x,const std::vector<std::size_t>& is)const
    {
        if constexpr (I>0)
        return std::apply([this,&x](auto&&... i){return std::invoke(std::get<I>(indexes()).second,x,i...);}, to_subarray<I>(is));
        else return std::invoke(std::get<I>(indexes()).second,x);

    }

    template<std::size_t I>
    std::vector<std::size_t>& next(const Object& x,std::vector<std::size_t>& i)const
    {
        auto n=size<I>(x,i);
        i[I]++;
        if(i[I]<n)
            return i;
        else
        {
            if constexpr (I>0)
            {
                i[I]=0;
                return next<I-1>(x,i);
            }
            else
                return i;
        }

    }
    std::vector<std::size_t>& next(const Object& x,std::vector<std::size_t>& i)const
    {
        return next<IndexNumber-1>(x,i);
    }

    std::ostream& put(std::ostream& os,const Object& x, const std::vector<std::string>& index_names , const std::vector<std::size_t>& is, const std::string& sep)const
    {
        if (isAnIndex(index_names))
            return put(os,x,to_array(index_names,is,std::make_index_sequence<IndexNumber>()))<<sep;
        else if(sameIndexes(index_names, std::make_index_sequence<IndexNumber>()))
            return put(os,x,to_array(is))<<sep;
        else
            return os;
    }

    std::ostream& put_title(std::ostream& os, const std::vector<std::string>& index_names , const std::string& sep)const
    {
        if (isAnIndex(index_names))
            return os<<getFieldName()<<sep;
        else if(sameIndexes(index_names, std::make_index_sequence<IndexNumber>()))
            return os<<getFieldName()<<sep;
        else
            return os;
    }


    bool next_if_same(const Object& x, const std::vector<std::string>& indexes , std::vector<std::size_t>& is)const
    {
        if (sameIndexes(indexes, std::make_index_sequence<IndexNumber>()))
        {
            next(x,is);
            return true;
        }
        else
            return false;
    }
    bool start_if_same(const Object& x, const std::vector<std::string>& indexes , std::vector<std::size_t>& is)const
    {
        if (sameIndexes(indexes, std::make_index_sequence<IndexNumber>()))
        {
            is=start(x);
            return true;
        }
        else
            return false;
    }
    bool end_if_same(const Object& x, const std::vector<std::string>& indexes , std::vector<std::size_t>& is)const
    {
        if (sameIndexes(indexes, std::make_index_sequence<IndexNumber>()))
        {
            is=end(x);
            return true;
        }
        else
            return false;
    }


private:
    template<std::size_t I>
    static std::array<std::size_t,I> to_subarray(const std::vector<std::size_t>& is)
    {
        assert(is.size()==IndexNumber);
        std::array<std::size_t,I> ind;
        for (std::size_t i=0; i<I; ++i) ind[i]=is[i];
        return ind;
    }

    static std::array<std::size_t,IndexNumber> to_array(const std::vector<std::size_t>& is)
    {
        assert(is.size()==IndexNumber);
        std::array<std::size_t,IndexNumber> ind;
        std::copy(is.begin(),is.end(),ind.begin());
        return ind;
    }

    static std::size_t name_to_i(const std::string& name,const std::vector<std::string>& names,const std::vector<std::size_t>& is)
    {
        std::size_t i=std::find(names.begin(),names.end(),name)-names.begin();
        assert(i<is.size());
        return is[i];
    }

    template<std::size_t... Is>
    std::array<std::size_t,IndexNumber> to_array(const std::vector<std::string>& names,const std::vector<std::size_t>& is, std::index_sequence<Is...>)const
    {
        return std::array<std::size_t,IndexNumber>{name_to_i(std::get<Is>(indexes()).first,names,is)...};
    }

    std::string fieldname_;
    std::tuple<std::pair<std::string,Sizes>...> indexes_;
    Fun datafun_;
};


template<class... Fields>
struct CT_s
{
    constexpr static auto const title =  (my_static_string("")+...+Fields::title);
};

template<std::size_t N, class SuperFunction, class SubFunction>
struct CF_s
{
private:




    template< class Object, std::size_t M,std::size_t ...Supers, std::size_t ...Subs>

    std::invoke_result_t<SubFunction,
                         std::invoke_result_t<SuperFunction,Object,decltype(Supers)...>,
                         decltype (Subs)...>

    op_impl(const Object& x, std::array< std::size_t,M> i, std::index_sequence<Supers...>, std::index_sequence<Subs...>)const
    {
        return std::invoke(sub_,std::invoke(f_,x,std::get<Supers>(i)...),std::get<N+Subs>(i)...);
    }

public:
    CF_s(const SuperFunction& f, SubFunction&& sub):f_(f),sub_{std::move(sub)}{}

      template<class Object, typename... Int>
        auto operator ()(const Object& x, Int... is) const->
        decltype (this->template op_impl(x,std::declval<std::array<std::size_t,sizeof...(Int)>>(),                                                                         std::make_index_sequence<N>(), std::make_index_sequence<sizeof...(Int)-N>()))

    {
          constexpr auto M=sizeof...(Int);
          return op_impl(x,std::array<std::size_t,M>{is...},
                         std::make_index_sequence<N>(),
                         std::make_index_sequence<M-N>());
    }

    SuperFunction f_;
    SubFunction sub_;
};


template<class Field, class Fun>
class F_s
{
public:
    F_s(Field, Fun&& f):f_{f}{}


    template<class Object, typename... Int>
    std::invoke_result_t<Fun,Object,Int...>
    operator ()(const Object& x, Int... i)const {
        return std::invoke(f_,x,i...);
    }

    static auto title() { return Field::title();}

    auto & f(){return f_;}
    auto & f()const {return f_;}

private:
    Fun f_;
};

template<class Index, class Size, class Fun>
class I_s
{
public:
    I_s(Index, Size&& s, Fun&& f):s_{std::move(s)},f_{std::move(f)}{}
    I_s(Index, Size&& s, const Fun& f):s_{std::move(s)},f_{f}{}

    template<class Object, typename... Int>
    std::size_t size (const Object& x, Int... i)const {
        return std::invoke(s_,x,i...);
    }
    template<class Object, typename... Int>
    auto operator () (const Object& x, Int... i)const {
        return std::invoke(f_,x,i...);
    }

    static auto title() { return Index::title();}

    auto & f(){return f_;}
    auto & f()const {return f_;}
    auto & s(){return s_;}
    auto & s()const {return s_;}


private:
    Size s_;
    Fun f_;
};


template<class Index, class Size>
class I_s<Index,Size,int>
{
public:
    I_s(Index, Size&& s,int=0):s_{s}{}

    template<class Object, typename... Int>
    std::size_t size (const Object& x, Int... i)const {
        return std::invoke(s_,x,i...);
    }
    template<class Object, typename... Int>
    std::size_t operator ()(const Object& , Int... is) const
    {
        return std::array<std::size_t, sizeof...(Int)>{is...}.back();
    }

    static auto title() { return Index::title();}

    auto & s(){return s_;}
    auto & s()const {return s_;}

private:
    Size s_;
};



template<class Index, class Size,class Fun, class Size2,class Fun2>
I_s<Index,Size,Fun> chooseWithFun(I_s<Index,Size,Fun>&& one,I_s<Index,Size2,Fun2>&& )
{
    return std::move(one);
}

template<class Index, class Size,class Fun, class Size2,class Fun2>
    auto chooseWithFun(I_s<Index,Size,int>&& ,I_s<Index,Size2,Fun2>&& two)->
    std::enable_if<!std::is_same_v<Fun2,int>,I_s<Index,Size2,Fun2>>
{
    return std::move(two);
}



template<std::size_t Parent_I,class ParentField,class  ParentFun, class Field, class Fun>
auto Compose_fun(F_s<ParentField,ParentFun> pfun,F_s<Field,Fun>&& f)
{
    return F_s<CT_s<ParentField,Field>,CF_s<Parent_I,ParentFun,Fun>>(CT_s<ParentField,Field>{},CF_s<Parent_I,ParentFun,Fun>(pfun.f(),std::move(f.f())));
}

template<std::size_t Parent_I,class  ParentFun, class Field, class Fun>
auto Compose_fun( ParentFun pfun,F_s<Field,Fun>&& f)
{
    return F_s<Field,CF_s<Parent_I,ParentFun,Fun>>(Field{},CF_s<Parent_I,ParentFun,Fun>(pfun,std::move(f.f())));
}


template<std::size_t Parent_I,class  ParentFun, class Index, class Size>
auto Compose_fun_Index(const ParentFun& pfun,I_s<Index,Size,int>&& index)
{
    return I_s<Index,CF_s<Parent_I,ParentFun,Size>,int>(Index{},CF_s<Parent_I,ParentFun,Size>(pfun,std::move(index.s())),0);
}

template<std::size_t Parent_I,class  ParentFun, class Index, class Size, class Fun>
auto Compose_fun_Index(const ParentFun& pfun,I_s<Index,Size,Fun>&& index)
{
    return I_s<Index,CF_s<Parent_I,ParentFun,Size>,CF_s<Parent_I,ParentFun,Fun>>
        (Index{},CF_s<Parent_I,ParentFun,Size>(pfun,std::move(index.s())),
         CF_s<Parent_I,ParentFun,Fun>(pfun,std::move(index.f())));
}





template<class Indexes,class Fun >
class Data_Index_static;

template<class... Index, class... Size,class...IndexFun,class... Field, class... Fun>
class Data_Index_static<Cs<I_s<Index,Size,IndexFun>...>,Cs<F_s<Field,Fun>...>>
{
public:
    static constexpr std::size_t N=sizeof...(Index);

   Data_Index_static(std::tuple<I_s<Index,Size,IndexFun>...>&& Ind,
                      std::tuple<F_s<Field,Fun>...>&& ftu)
        :Ind_{std::move(Ind)},ftu_{std::move(ftu)}{}

   static std::ostream& put_title(std::ostream& os, const std::string& sep)
    {
       ((os<<Index::title.str()<<sep),...);
       ((os<<Field::title.str()<<sep),...);
       return os<<"\n";
    }

    static bool print_title(const std::string& filename ,const std::string& ext, const std::string& sep)
    {
        std::string fname=filename+get_index_names()+ext;
        std::ofstream f(fname.c_str());
        put_title(f,sep);
        return true;
    }
    template<class Object>
    bool print_data(const Object& x,const std::string& filename ,const std::string& ext, const std::string& sep)const
    {
        std::stringstream ss;
        put_index_sample(ss,x,sep);
        std::string fname=filename+get_index_names()+ext;
        std::ofstream f(fname.c_str(), std::ios_base::app);
        f<<ss.str();
        f.close();
        return true;
    }


    template<class Object>
    std::ostream& put_index_sample(std::ostream& os, const Object& x, const std::string& sep)const
    {
        return put_index_sample_impl(os,x,sep,std::tuple<>(),std::index_sequence_for<Index...>());
    }

    static std::string get_index_names() {
        return (Index::title.str()+...+"");
    }
    auto& Indexes(){return Ind_;}
    auto& Indexes()const {return Ind_;}
    auto& Data(){return ftu_;}
    auto& Data()const{return ftu_;}


private:
    template<class Object, class... Index_Type,typename... Int>
    std::ostream& put_index_sample_impl(std::ostream& os, const Object& x,const std::string& sep,std::tuple<Index_Type...> out, std::index_sequence<>,Int...is)const
    {
        std::apply([this,&x,&os,&sep](auto&...f){
            ((os<<f<<sep),...);},out);
         std::apply([this,&x,&os,&sep,is...](auto&...f){
            ((os<<f(x,is...)<<sep),...);},ftu_);
         return os<<"\n";
    }

    template<std::size_t I,std::size_t... Is, class Object, class...Index_Type,typename... Int>
    std::ostream& put_index_sample_impl(std::ostream& os, const Object& x, const std::string& sep,std::tuple<Index_Type...> out,std::index_sequence<I,Is...>,Int...is)const
    {
        for (std::size_t i=0; i<Index_size<I>(x,is...); ++i)
            put_index_sample_impl(os,x,sep,std::tuple_cat(out,std::make_tuple(std::get<I>(Ind_)(x,is...,i))),
                                  std::index_sequence<Is...>(),is...,i);
        return os;
    }
    template<std::size_t I, class Object,typename...Int>
    std::size_t Index_size(const Object& x, Int... is)const
    {
        return std::get<I>(Ind_).size(x,is...);
    }
    std::tuple<I_s<Index,Size,IndexFun>...> Ind_;
    std::tuple<F_s<Field,Fun>...> ftu_;
};



template<class... Index, class... Size,class...FunIndex,class... Field, class... Fun>
Data_Index_static<Cs<I_s<Index,Size,FunIndex>...>,Cs<F_s<Field,Fun>...>>
make_data_static(std::tuple<I_s<Index,Size, FunIndex>...>&& Ind,
                 std::tuple<F_s<Field,Fun>...>&& ftu)
{
    return Data_Index_static<Cs<I_s<Index,Size, FunIndex>...>,Cs<F_s<Field,Fun>...>>(std::move(Ind), std::move(ftu));
}

//template<class... NewIndex, class... NewSize,class ...FunIndex,class PF, class PFun, class I, class D>
//auto Compose_static(std::tuple<I_s<NewIndex, NewSize,FunIndex>...> s,
//                          F_s<PF,PFun> f,
//                    Data_Index_static<I,D> && d)
//{
//    constexpr auto N=sizeof...(NewIndex);
//    return
//        make_data_static(std::tuple_cat(std::move(s),
//                                         std::apply([&f](auto&&... ind)
//                                                    { return std::make_tuple(Compose_fun_Index<N>(f.f(),std::move(ind))...);},
//                                                    std::move(d.Indexes()))),
//                          std::apply([&f](auto&&...tu){return std::make_tuple(Compose_fun<N>(f,std::move(tu))...);},std::move(d.Data()))
//            );
//}

template<class... NewIndex, class... NewSize,class ...FunIndex,class PFun,class I, class D>
auto Compose_static(std::tuple<I_s<NewIndex, NewSize,FunIndex>...> s,
                    PFun f,
                    Data_Index_static<I,D> && d)
{
    constexpr auto N=sizeof...(NewIndex);
    return
        make_data_static(std::tuple_cat(std::move(s),std::apply([&f](auto&&... ind)
                                                    { return std::make_tuple(Compose_fun_Index<N>(f,std::move(ind))...);},
                                                    std::move(d.Indexes()))),
                          std::apply([&f](auto&&...tu){return std::make_tuple(Compose_fun<N>(f,std::move(tu))...);},std::move(d.Data()))
                                          );
}


template<class PF, class PFun,class I, class D>
auto Compose_static(F_s<PF,PFun> f,
                    Data_Index_static<I,D> && d)
{
    return
        make_data_static(std::apply([&f](auto&&... ind)
                                                    { return std::make_tuple(Compose_fun_Index<0>(f,std::move(ind))...);},
                                                    std::move(d.Indexes())),
                          std::apply([&f](auto&&...tu)
                                    {return std::make_tuple(Compose_fun<0>(f,std::move(tu))...);},std::move(d.Data()))
                                          );
}


template<class PFun,class I, class D>
auto Compose_static(PFun f,
                    Data_Index_static<I,D> && d)
{
    return
        make_data_static(std::apply([&f](auto&&... ind)
                                                    { return std::make_tuple(Compose_fun_Index<0>(f,std::move(ind))...);},
                                                    std::move(d.Indexes())),
                          std::apply([&f](auto&&...tu)
                                    {return std::make_tuple(Compose_fun<0>(f,std::move(tu))...);},
                                     std::move(d.Data()))
                                          );
}


//template<class... NewIndex, class... NewSize,class... FunIndex,class PF, class PFun,class... Data_Index_static_R>
//auto Compose_static(std::tuple<I_s<NewIndex, NewSize,FunIndex>...> s,
//                    F_s<PF,PFun> f,
//                    std::tuple<Data_Index_static_R...> && d)
//{
//    return std::apply([&s,&f](auto&&...D){ return std::make_tuple(Compose_static(s,f,std::move(D))...);},std::move(d));

//}

template<class... NewIndex, class... NewSize,class... FunIndex,class PFun,class... Data_Index_static_R>
auto Compose_static(std::tuple<I_s<NewIndex, NewSize,FunIndex>...> s,
                    PFun f,
                    std::tuple<Data_Index_static_R...> && d)
{
    return std::apply([&s,&f](auto&&...D){ return std::make_tuple(Compose_static(s,f,std::move(D))...);},std::move(d));

}



template<class PF, class PFun,class... Data_Index_static_R>
auto Compose_static(F_s<PF,PFun> f,
                    std::tuple<Data_Index_static_R...> && d)
{
    return std::apply([&f](auto&&...D){ return std::make_tuple(Compose_static(f,std::move(D))...);},std::move(d));

}


template<class PFun,class... Data_Index_static_R>
auto Compose_static(PFun f,
                    std::tuple<Data_Index_static_R...> && d)
{
    return std::apply([&f](auto&&...D){ return std::make_tuple(Compose_static(f,std::move(D))...);},std::move(d));
}





template<std::size_t...Is,class... Index, class... Size1,class ...IFun1,class... Size2,class ...IFun2>
auto chooseWithFun(std::tuple<I_s<Index,Size1,IFun1>...>&&one,
                   std::tuple<I_s<Index,Size2,IFun2>...>&&two,
                   std::index_sequence<Is...>)
{
    return std::make_tuple(chooseWithFun(std::move(std::get<Is>(one)),std::move(std::get<Is>(two)))...);
}


template<class Index1,class Index2,class Size1,class IFun1,class Size2,class IFun2>
auto fuseWithFun_Is(I_s<Index1,Size1,IFun1> &&one,
                        I_s<Index2,Size2,IFun2> &&two)->
    std::enable_if_t<
    !std::is_same_v<Index1,Index2> || (std::is_same_v<IFun1,int> == std::is_same_v<IFun2,int>),
    std::pair<I_s<Index1,Size1,IFun1>,I_s<Index2,Size2,IFun2>>
        >

{

    return std::make_pair(std::move(one), std::move(two));
}


template<class Index,class Size1,class IFun,class Size2>
auto fuseWithFun_Is(I_s<Index,Size1,IFun> &&one,
               I_s<Index,Size2,int> &&two)->
    std::enable_if_t<
    ! std::is_same_v<IFun,int>,
    std::pair<I_s<Index,Size1,IFun>,I_s<Index,Size2,IFun>>
        >
{

    auto newtwo=I_s(Index{},std::move(two.s()),one.f());
    return std::make_pair(std::move(one), std::move(newtwo));
}

template<class Index,class Size1,class IFun,class Size2>
auto fuseWithFun_Is(I_s<Index,Size1,int> &&one,
                   I_s<Index,Size2,IFun> &&two)->
    std::enable_if_t<
        ! std::is_same_v<IFun,int>,
    std::pair<I_s<Index,Size1,IFun>,I_s<Index,Size2,IFun>>
    >
{

    //typedef typename IFun::df IFUN;
    auto newone=I_s(Index{},std::move(one.s()),two.f());
    return std::make_pair(std::move(newone), std::move(two));
}


template<class...D01, class ...D02,class...D11, class ...D12, std::size_t...I1,std::size_t...I2,std::size_t I>
auto fuseWithFun_tuple(std::tuple<D01...>&& in1,
                       std::tuple<D11...>&& out1,
                       std::index_sequence<I,I1...>,

                       std::tuple<D02...>&& in2,
                       std::tuple<D12...>&& out2,
                       std::index_sequence<I,I2...>)
{

    auto [ind1,ind2]=fuseWithFun_Is(std::move(std::get<I>(in1)),std::move(std::get<I>(in2)));
    return fuseWithFun_tuple(std::move(in1),
                             std::tuple_cat(std::move(out1), std::make_tuple(std::move(ind1))),
                             std::index_sequence<I1...>(),
                             std::move(in2),
                             std::tuple_cat(std::move(out2), std::make_tuple(std::move(ind2))),
                             std::index_sequence<I2...>());

}

template<class...D01, class ...D02,class...D11, class ...D12, std::size_t...I1,std::size_t I>
auto fuseWithFun_tuple(std::tuple<D01...>&& in1,std::tuple<D11...>&& out1,std::index_sequence<I,I1...>,
                       std::tuple<D02...>&& ,std::tuple<D12...>&& out2,std::index_sequence<>)
{
    return std::make_pair(std::tuple_cat(std::move(out1), std::make_tuple(std::move(std::get<I>(in1)),std::move(std::get<I1>(in1))...)),
                          std::move(out2));

}


template<class...D01, class ...D02,class...D11, class ...D12, std::size_t...I2>
auto fuseWithFun_tuple(std::tuple<D01...>&& ,std::tuple<D11...>&& out1,std::index_sequence<>,
                       std::tuple<D02...>&& in2,std::tuple<D12...>&& out2,std::index_sequence<I2...>)
{
    return std::make_pair(std::move(out1),std::tuple_cat(std::move(out2), std::make_tuple(std::move(std::get<I2>(in2))...)));

}



template<class ...D1, class... D2>
auto fuseWithFun_tuple(std::tuple<D1...>&& tu1, std::tuple<D2...>&& tu2)
{
    return fuseWithFun_tuple(std::move(tu1),std::tuple<>(),std::index_sequence_for<D1...>(),
                       std::move(tu2),std::tuple<>(),std::index_sequence_for<D2...>());

}





template<bool TF,class... D0,class... D1,class... Index1, class... Size1,
          class...IFun1,class... Field1, class... Fun1,class... Index2,class... Size2,class...IFun2,class... Field2, class... Fun2, class...R>
auto Concatenate_tuple_static_impl(std::integral_constant<
                                       std::enable_if_t<(!std::is_same_v<Cs<Index1...>,Cs<Index2...> >),bool>,
                                       TF>,
                                   std::tuple<D0...>&& tu,
                                   std::tuple<D1...>&& tu1,
                                   Data_Index_static<Cs<I_s<Index1,Size1,IFun1>...>,Cs<F_s<Field1,Fun1>...>>&& e,
                                   Data_Index_static<Cs<I_s<Index2,Size2,IFun2>...>,Cs<F_s<Field2,Fun2>...>> &&   eTu,
                                   R&&...r)
{
    if constexpr (true)
    {


        auto t=fuseWithFun_tuple(std::move(e.Indexes()),std::move(eTu.Indexes()));

        auto  new_e=make_data_static(std::move(t.first), std::move(e.Data()));
        auto  new_eTu=make_data_static(std::move(t.second), std::move(eTu.Data()));


      //  typedef typename decltype(e)::e  old_e;
    //    typedef typename decltype(new_e)::newe  new_e_e;




        if constexpr (TF)
        {
            return Concatenate_tuple_static_impl(std::true_type{},
                                                 std::move(tu),
                                                 std::tuple_cat(std::move(tu1),std::make_tuple(std::move(new_eTu))),
                                                 std::move(new_e),std::move(r)...);
        }
        else {
            static_assert(sizeof...(D1)==0);
            return Concatenate_tuple_static_impl(std::false_type{},
                                                 std::tuple_cat(std::move(tu),std::make_tuple(std::move(new_eTu))),
                                                 std::move(tu1),
                                                 std::move(new_e),std::move(r)...);

        }
    }
    else {
        if constexpr (TF)
        {
            return Concatenate_tuple_static_impl(std::true_type{},
                                                 std::move(tu),
                                                 std::tuple_cat(std::move(tu1),std::make_tuple(std::move(eTu))),
                                                 std::move(e),std::move(r)...);
        }
        else {
            static_assert(sizeof...(D1)==0);
            return Concatenate_tuple_static_impl(std::false_type{},
                                                 std::tuple_cat(std::move(tu),std::make_tuple(std::move(eTu))),
                                                 std::move(tu1),
                                                 std::move(e),std::move(r)...);

        }


    }

  }

  template<class...D,class... Index, class... Size1,class...IFun1,class... Field1, class... Fun1,class... Size2,class...IFun2,class... Field2, class... Fun2, class...R>
  auto Concatenate_tuple_static_impl(std::false_type,
                                     std::tuple<D...>&& tu,
                                     std::tuple<>,
                                     Data_Index_static<Cs<I_s<Index,Size2,IFun2>...>,Cs<F_s<Field2,Fun2>...>>&& e,
                                     Data_Index_static<Cs<I_s<Index,Size1,IFun1>...>,Cs<F_s<Field1,Fun1>...>>  &&eTu,
                                     R&&...r)
  {

        auto new_e=make_data_static(
              chooseWithFun(std::move(eTu.Indexes()),std::move(e.Indexes()),std::index_sequence_for<Index...>()),
              std::tuple_cat(std::move(eTu.Data()), std::move(e.Data())));


        return Concatenate_tuple_static_impl(std::true_type{},
                                            std::move(tu),std::tuple<>(),
                                           std::move(new_e),std::move(r)...);

  }

  template<class...D,class... Index, class... Size,class...IFun,class... Field, class... Fun>
  auto Concatenate_tuple_static_impl(std::false_type,std::tuple<D...>&& tu,std::tuple<>,
                                     Data_Index_static<Cs<I_s<Index,Size,IFun>...>,Cs<F_s<Field,Fun>...>>&& e)
  {
      return std::tuple_cat(std::move(tu),std::make_tuple(std::move(e)));
  }
  template<class...D0,class...D1,class... Index, class... Size,class...IFun,class... Field, class... Fun>
  auto Concatenate_tuple_static_impl(std::true_type,std::tuple<D0...>&& tu,std::tuple<D1...>&&tu1,
                                     Data_Index_static<Cs<I_s<Index,Size,IFun>...>,Cs<F_s<Field,Fun>...>>&& e)
  {
      return std::tuple_cat(std::move(tu),std::make_tuple(std::move(e)), std::move(tu1));
  }



template<class D,class...Datas,class... Index, class... Size,class...IFun,class... Field, class... Fun>

auto Concatenate_tuple_static(std::tuple<D,Datas...>&& tu, Data_Index_static<Cs<I_s<Index,Size,IFun>...>,Cs<F_s<Field,Fun>...>>&& e)
{
    auto out= std::apply([&e](auto&&...t){
        return Concatenate_tuple_static_impl(std::false_type{},
                                             std::tuple<>(),
                                             std::tuple<>(),
                                             std::move(e),
                                             std::move(t)...);},std::move(tu));


    //typedef typename decltype(out)::one one;

    return out;
}

template<class...Ds,class D2,class D3,class... R>
auto Concatenate_tuple_static(std::tuple<Ds...>&& one, D2&& r,D3&& s,R&&... rest)
{
    auto next=Concatenate_tuple_static(std::move(one),std::move(r));
    return Concatenate_tuple_static(std::move(next),std::move(s),std::move(rest)...);
}

template<class...Datas,class... Datas2>
auto Concatenate_tuple_static(std::tuple<Datas...>&& one,std::tuple<Datas2...>&& two)
{
    return std::apply([&one](auto&&...d){return Concatenate_tuple_static(std::move(one),std::move(d)...);},std::move(two));
}





template<class Parent, class getObject, class F, class... Size_t>
struct Inner_Composition{
    Inner_Composition(Cs<Parent>,const getObject& get, Cs<Size_t...>, const F& childFunc): g_{get},data_{childFunc} {}
    template<class... Ts>
    auto operator()(const Parent& p,Size_t...is, Ts...x)const
    {
        return std::invoke(data_,std::invoke(g_,p,is...),x...);
    }
    getObject g_;
    F data_;


};
template<class Parent, class getObject>
auto Inner_Composition_tuple(Cs<Parent>,const getObject& ,Cs<>,const std::tuple<>& )
{
    return std::tuple<>();
}

template<class Parent, class getObject, class... Size_t, class... Sizes >
auto Inner_Composition_tuple(Cs<Parent>,const getObject& get,Cs<Size_t...>,const std::tuple<Sizes...>& indexes)
{
    return std::apply([&get](auto& ...p){
        return std::make_tuple(std::make_pair(p.first,
                                              Inner_Composition(Cs<Parent>(),get,Cs<Size_t...>(),p.second ))...);
    },indexes);
}


template<class Object, class getSub,class DataIndex >
auto Insert(Cs<Object>,const getSub& sub,const DataIndex& v)
{

    auto sizes=Inner_Composition_tuple(Cs<Object>(),sub,Cs<>(),v.indexes());
    if constexpr (std::tuple_size_v<decltype(sizes)> >0)
    {

        return std::apply([&sub,&v](auto&...x){
            return Data_Index(Cs<Object>(),v.getFieldName(),
                              Inner_Composition(Cs<Object>(),sub,Cs<>(),v.dataFun()),x...);},
                          sizes);
    }
    else {
        return Data_Index(Cs<Object>(),v.getFieldName(),
                          Inner_Composition(Cs<Object>(),sub,Cs<>(),v.dataFun()));
    }
}

template<class Object, class getSub,class... Size_t,class DataIndex,class OuterSizes >
auto Insert(Cs<Object>,const getSub& sub,Cs<Size_t...>,const DataIndex& v,const OuterSizes  indexes)
{
    auto tu=Inner_Composition_tuple(Cs<Object>(),sub,Cs<Size_t...>(),v.indexes());
    return std::apply([&sub,&v,&indexes](auto&...x){
        return Data_Index(Cs<Object>(),v.getFieldName(),
                          Inner_Composition(Cs<Object>(),sub,Cs<Size_t...>(),v.dataFun()),indexes,x...);},
                      tu);
}


template<class Object, class getSub,class... Size_t,class DataIndex,class OuterS,class... OuterSizes >
auto Insert(Cs<Object>,const getSub& sub,Cs<Size_t...>,const DataIndex& v,const OuterS& ind,const OuterSizes&... indexes)
{
    return std::apply([&sub,&v,&ind,&indexes...](auto&&...x){
        return Data_Index(Cs<Object>(),v.getFieldName(),
                          Inner_Composition(Cs<Object>(),sub,Cs<Size_t...>(),v.dataFun()),ind,indexes...,x...);},
                      Inner_Composition_tuple(Cs<Object>(),sub,Cs<Size_t...>(),v.indexes()));
}


template<class Object, class getSub,class... Data_Indexes >
auto Insert_tuple(Cs<Object>,const getSub& sub,const std::tuple<Data_Indexes...>& v)
{
    return std::apply([&sub](auto&...d){return std::tuple(Insert(Cs<Object>(),sub,d)...);},v);
}



template<class Object, class getSub,class... S,class Data_Tuple, class OuterSizes>
auto Insert_tuple(Cs<Object>,const getSub& sub,Cs<S...>,const Data_Tuple& v,const OuterSizes& indexes)
{
    return std::apply([&sub, &indexes](auto&...d){return std::make_tuple(Insert(Cs<Object>{},sub,Cs<S...>{},d,indexes)...);},v);
}

template<class Object, class getSub,class... S,class Data_Tuple, class OuterS, class... OuterSizes>
auto Insert_tuple(Cs<Object>,const getSub& sub,Cs<S...>,const Data_Tuple& v,const OuterS ind,const OuterSizes&... indexes)
{
    return std::apply([&sub,&ind, &indexes...](auto&...d){return std::make_tuple(Insert(Cs<Object>(),sub,Cs<S...>(),d,ind,indexes...)...);},v);
}



template< std::size_t I,class Object,typename... Data_Indexes >
std::vector<std::size_t> start_index(const Object& x,const std::tuple<Data_Indexes...>& tu,  const std::vector<std::string>& index)
{
    std::vector<std::size_t> i;
    if constexpr (I<sizeof... (Data_Indexes))
    {
        if (std::get<I>(tu).start_if_same(x,index,i))
            return i;
        else
            return start_index<I+1>(x,tu,index);
    }
    else
    {
        assert(false);
        return i;
    }
 }

 template< std::size_t I,class Object,typename... Data_Indexes >
 std::vector<std::size_t> end_index(const Object& x,const std::tuple<Data_Indexes...>& tu,  const std::vector<std::string>& index)
 {
     std::vector<std::size_t> i;
     if constexpr (I<sizeof... (Data_Indexes))
     {
         if (std::get<I>(tu).end_if_same(x,index,i))
             return i;
         else
             return end_index<I+1>(x,tu,index);
     }
     else
     {
         assert(false);
         return i;
     }
 }


 template< std::size_t I,class Object,typename... Data_Indexes >
 std::vector<std::size_t> next_index(const Object& x,const std::tuple<Data_Indexes...>& tu,  const std::vector<std::string>& index, std::vector<std::size_t>& i)
 {
     if constexpr (I<sizeof... (Data_Indexes))
     {
         if (std::get<I>(tu).next_if_same(x,index,i))
             return i;
         else
             return next_index<I+1>(x,tu,index,i);
     }
     else
     {
         assert(false);
         return i;
     }
 }

 template< std::size_t I,class Object,typename... Data_Indexes >
 std::ostream& put_index(std::ostream& os,const Object& x,const std::tuple<Data_Indexes...>& tu,  const std::vector<std::string>& index, const std::vector<std::size_t>& i,const std::string& sep)
 {
     if constexpr (I<sizeof... (Data_Indexes))
     {
         std::get<I>(tu).put(os,x,index,i,sep);
         return put_index<I+1>(os,x,tu,index,i,sep);
     }
     else
     {
         return os;
     }
 }

 template< std::size_t I,typename... Data_Indexes >
 std::ostream& put_index_title(std::ostream& os,const std::tuple<Data_Indexes...>& tu,  const std::vector<std::string>& index,const std::string& sep)
 {
     if constexpr (I<sizeof... (Data_Indexes))
     {
         std::get<I>(tu).put_title(os,index,sep);
         return put_index_title<I+1>(os,tu,index,sep);
     }
     else
     {
         return os;
     }
 }




 template<typename... Data_Indexes >
std::map<std::vector<std::string>, std::vector<std::string>>&
get_index_map(const std::tuple<Data_Indexes...>& ,  std::map<std::vector<std::string>, std::vector<std::string>>& map, std::index_sequence<>)
{
    return map;
}


template<std::size_t I,std::size_t...Is,typename... Data_Indexes >
std::map<std::vector<std::string>, std::vector<std::string>>&
get_index_map(const std::tuple<Data_Indexes...>& tu,  std::map<std::vector<std::string>, std::vector<std::string>>& map, std::index_sequence<I,Is...>)
{
    map[std::get<I>(tu).indexed_by()].push_back(std::get<I>(tu).getFieldName());
    return get_index_map(tu, map, std::index_sequence<Is...>());
}


template<typename... Data_Indexes>
auto get_index_titles(const std::tuple<Data_Indexes...>& tu)
{
    std::map<std::vector<std::string>, std::vector<std::string>> map;
    map=get_index_map(tu,map,std::index_sequence_for<Data_Indexes...>());



    std::vector<std::string> set_titles(map.size());

    std::size_t i=0;
    for (auto& e: map)
    {

        for (auto& in:e.first)
            set_titles[i]+="_"+in;
        ++i;

    }

    std::vector<std::vector<std::string>> names_titles(map.size());
    i=0;
    for (auto it=map.begin(); it!=map.end(); ++it,++i)
    {
        names_titles[i]=it->second;
    }
    return std::make_pair(set_titles,names_titles);

}

template<typename... Data_Indexes>
auto get_index_files(const std::tuple<Data_Indexes...>& tu)
{
    std::map<std::vector<std::string>, std::vector<std::string>> map;
    map=get_index_map(tu,map,std::index_sequence_for<Data_Indexes...>());



    std::vector<std::string> set_titles(map.size());

    std::size_t i=0;
    for (auto& e: map)
    {

        for (auto& in:e.first)
            set_titles[i]+="_"+in;
        ++i;

    }

    std::vector<std::vector<std::string>> index_names(map.size());
    i=0;
    for (auto it=map.begin(); it!=map.end(); ++it,++i)
    {
        index_names[i]=it->first;
    }
    return std::make_pair(set_titles,index_names);

}





template<class Object,typename... Data_Indexes>
std::ostream& put_index_sample(std::ostream& os, const Object& x,const std::tuple<Data_Indexes...>& tu,const std::vector<std::string>& indexes, const std::string& sep)
{
        std::vector<std::size_t> is=start_index<0>(x,tu,indexes);
        std::vector<std::size_t> end=end_index<0>(x,tu,indexes);
        while (is<end)
        {

            put_index<0>(os,x,tu,indexes,is,sep)<<"\n";
            next_index<0>(x,tu,indexes,is);
        }
        return os;
}

template<typename... Data_Indexes>
std::ostream& put_index_title(std::ostream& os, const std::tuple<Data_Indexes...>& tu,const std::vector<std::string>& indexes, const std::string& sep)
{
    put_index_title<0>(os,tu,indexes,sep)<<"\n";
    return os;
}







class Data_Index_point
{
    std::map<std::string, std::pair<std::size_t, std::size_t>> pos_;

public:

    void push_back(std::string&& name, std::size_t n)
    {
        pos_[std::move(name)]=std::pair(0ul,n);
    }

    Data_Index_point begin(){
        Data_Index_point out(*this);
        for (auto& e:out.pos_) e.second.first=0;
        return out;
    }

    bool isEnd()const
    {
        return ! (pos_.rbegin()->second.first<pos_.rbegin()->second.second);
    }

    Data_Index_point& operator++(){
        auto it=pos_.begin();
        while (true)
        {
            ++(it->second.first);
            if ( it->second.first<it->second.second)
                return *this;
            else {
                auto it2=it;
                it2++;
                if ((it2) !=pos_.end()){
                    it->second.first=0;
                    ++it;
                }
                else return *this;
            }
        }
    }



};


class Data_Index_scheme
{
    std::map<std::set<std::string>, std::vector<std::string> > d_;


    static std::set<std::string> intersects(const std::set<std::string>& one, std::set<std::string>& two)
    {
        std::set<std::string> out;
        std::set_intersection(one.begin(),one.end(),two.begin(),two.end(),std::inserter(out,out.begin()));
        return out;
    }

    static std::set<std::string> set_difference(const std::set<std::string>& one, std::set<std::string>& two)
    {
        std::set<std::string> out;
        std::set_difference(one.begin(),one.end(),two.begin(),two.end(),std::inserter(out,out.begin()));
        return out;
    }


    template<class inputIter>
    static std::set<std::string> intersects(inputIter it, inputIter end)
    {
        auto init=it->first;
        ++it;
        for (; it!=end; ++it)
        {
            init=intersects(init,it->first);
        }
        return init;
    }

public:

    auto & self_map()const {return d_;}

    std::set<std::string> first_index()const {return d_.begin()->first;}



    typedef   Data_Index_scheme self_type;

    void push_back(std::string&& name, std::set<std::string>&& new_indexes)
    {
        d_[new_indexes].push_back(std::move(name));
    }




    friend Data_Index_scheme concatenate_impl( Data_Index_scheme&& one, const Data_Index_scheme& two)
    {
        for (auto& e:two.d_)
        {
            auto& p=one.d_[e.first];
            p.insert(p.end(),e.second.begin(),e.second.end());
        }
        return std::move(one);
    }

    friend Data_Index_scheme concatenate( Data_Index_scheme&& one){return std::move(one);}
    template<class... Data_Index>
    friend Data_Index_scheme concatenate( Data_Index_scheme&& one, const Data_Index_scheme&  two, const Data_Index&... three)
    {
        return concatenate(concatenate_impl(std::move(one),two),three...);
    }





    template<class ...T>
    static Data_Index_scheme data_index();



    std::size_t size()const { return d_.size();}


    std::vector<std::string> set_titles()const
    {
        std::vector<std::string> out(d_.size());

        std::size_t i=0;
        for (auto& e: d_)
        {

            for (auto& in:e.first)
                out[i]+="_"+in;
            ++i;

        }
        return out;
    }

    std::vector<std::vector<std::string>> index_titles()const
    {
        std::vector<std::vector<std::string>> out(size());
        std::size_t i=0;
        for (auto it=d_.begin(); it!=d_.end(); ++it,++i)
        {
            std::vector<std::string> e;
            e.insert(e.end(),it->first.begin(),it->first.end());
            out[i]=e;
        }
        return out;

    }

    std::vector<std::vector<std::string>> names_titles()const
    {
        std::vector<std::vector<std::string>> out(size());
        std::size_t i=0;
        for (auto it=d_.begin(); it!=d_.end(); ++it,++i)
        {
            out[i]=it->second;
        }
        return out;

    }



    void insert_index(std::string&& new_index_name)
    {
        std::map<std::set<std::string>, std::vector<std::string> > new_d;
        for (auto& e:d_)
        {
            auto i=e.first;
            i.insert(new_index_name);
            new_d[std::move(i)]=std::move(e.second);
        }
        d_=std::move(new_d);
    }
    Data_Index_scheme& insert_indexes(const std::vector<std::string>& new_index_name)
    {
        std::map<std::set<std::string>, std::vector<std::string> > new_d;
        for (auto& e:d_)
        {
            auto i=e.first;
            i.insert(new_index_name.begin(), new_index_name.end());
            new_d[std::move(i)]=std::move(e.second);
        }
        d_=std::move(new_d);
        return *this;
    }

    Data_Index_scheme operator*(std::vector<std::string>&& new_index_names)const &
    {
        Data_Index_scheme out(*this);
        return out.insert_indexes(std::move(new_index_names));
    }

    Data_Index_scheme operator*(std::vector<std::string>&& new_index_names)&&
    {
        insert_indexes(std::move(new_index_names));
        return *this;
    }






    /*
    static auto get_constructor_fields()
    {
        return std::make_tuple(
            grammar::field(C<self_type>{},"value_name",&self_type::value_name),
            grammar::field(C<self_type>{},"index_names",&self_type::index_names));
    }


    Data_Index_scheme(const std::string& name, const std::vector<std::string>& indexes)
        :value_name_{name},index_names_{indexes}{}

    Data_Index_scheme( std::string&& name,  std::vector<std::string>&& indexes)
        :value_name_{std::move(name)},index_names_{std::move(indexes)}{}

    Data_Index_scheme()=default;

*/


};




Data_Index_scheme Insert_Index(Data_Index_scheme&& v, std::vector<std::string>&& indexes)
{
    v.insert_indexes(std::move(indexes));
    return std::move(v);
}


template<class...> struct myData_Index;




template<class T>
struct myData_Index<T>
{
    static Data_Index_scheme data_index(){ return T::data_index();}

    static auto get_data_index(){return T::get_data_index();}
    static auto get_data_index(const std::string& name){return T::get_data_index(name);}
    static auto get_data_index_static(){return T::get_data_index_static();}

};

template<class... T>
struct myData_Index<std::tuple<T...>>
{
    typedef std::tuple<T...> self_type;

    template<std::size_t... Is>
    static auto get_data_index_impl(const std::string& name, std::index_sequence<Is...>)
    {
        return std::tuple_cat(
            Insert_tuple(Cs<self_type>(),
                         [](const self_type& self) {return std::get<Is>(self);},
                         myData_Index<std::tuple_element_t<Is,self_type> >::get_data_index(name))...);
    }

    static auto get_data_index(const std::string& name){
        return get_data_index_impl(name,std::index_sequence_for<T...>());
    }
        template<std::size_t... Is>
        static auto get_data_index_static_impl(std::index_sequence<Is...>)
        {
            return std::tuple_cat(
                Compose_static([](const self_type& self) {return std::get<Is>(self);},
                             myData_Index<std::tuple_element_t<Is,self_type> >::get_data_index_static())...);
        }

        static auto get_data_index_static(){
            return get_data_index_static_impl(std::index_sequence_for<T...>());

        }
};


template<class ...T>
Data_Index_scheme Data_Index_scheme::data_index()
{
    return concatenate(myData_Index<T>::data_index()...);
}


template <typename T, typename S,template <typename> class C>
C<T> emptyCopy(const C<S>& x);


template <typename T, typename S,template <typename> class C>
C<T> emptyCopy(const C<S>& x, const T& e);

template <typename T, typename S>
std::vector<T> emptyCopy(const std::vector<S>& x)
{
    return std::vector<T>(x.size());
}


template <typename T, typename S>
M_Matrix<T> emptyCopy(const M_Matrix<S>& x)
{
    return M_Matrix<T>(x.nrows(),x.ncols(),static_cast<typename M_Matrix<T>::TYPE>(x.type()));
}

template <typename T, typename S>
std::vector<T> emptyCopy(const std::vector<S>& x, const T& e)
{
    return std::vector<T>(x.size(),e);
}

template <typename T, typename S>
M_Matrix<T> emptyCopy(const M_Matrix<S>& x, const T& e)
{
    return M_Matrix<T>(x.nrows(),x.ncols(),e);
}

std::mt19937_64
init_mt(
    std::mt19937_64::result_type initseed,
    std::ostream& os)
{
    if (initseed==0)
    {
        std::random_device rd;
        initseed=rd();
        os<<"\n ## init seed generated ## \n initseed="<<initseed<<"\n";
    }
    else
        os<<"\n ## provided seed  ## \n initseed="<<initseed<<"\n";
    return std::mt19937_64(initseed);
}






template<typename... >
class moments;


enum MOMENTS {COUNT,SUM,SUM_SQR,MEAN,VARIANCE,CHI_VALUE,T_VALUE,T2_VALUE};

const std::string MOMENTS_string[]={"count","sum","sum_sqr","mean", "variance","chi_value","t_value","T2_value"};

constexpr MOMENTS MOMENTS_VALUES[]={COUNT,SUM,SUM_SQR,MEAN,VARIANCE,CHI_VALUE,T_VALUE,T2_VALUE};

template<>
class moments<double>
{

public:
    typedef  moments self_type;
    constexpr static auto  className=my_static_string("moments_")+my_trait<double>::className;
    static auto get_constructor_fields()
    {
        return std::make_tuple(
            grammar::field(C<self_type>{},"count",&self_type::count),
            grammar::field(C<self_type>{},"sum",&self_type::sum),
            grammar::field(C<self_type>{},"sum_sqr",&self_type::sum_sqr)
                                                                                                                    );
    }

    moments(std::size_t count, double sum, double sum_sqr):count_(count),sum_{sum},sum_sqr_{sum_sqr}{};



    moments(const std::vector<double>& data):
                                               count_(data.size()),sum_(op::sum(data)),sum_sqr_(op::sum_sqr(data)){}

    template<class C,class F>
    moments(const std::vector<C>& data, const F& f):
                                                      count_(data.size()),sum_(op::sum(data,f)),sum_sqr_(op::sum_sqr(data,f)){}


    std::size_t data_size()const { return 1;}


    std::size_t data_dim(MOMENTS )const { return 0;}


    std::tuple<double>
    data_row(MOMENTS m)const {
        switch(m){
        case COUNT: return count();
        case SUM: return sum();
        case SUM_SQR: return  sum_sqr();
        case MEAN: return mean();
        case VARIANCE:  return variance();
        case CHI_VALUE: return chi_value();
        case T_VALUE: return t_value();
        case T2_VALUE:
        default:
            return T2_value();
        }
    }

    std::size_t count()const {return count_;}
    double sum()const { return sum_;}
    double sum_sqr()const {return sum_sqr_;}
    double mean()const {return sum()*(1.0/count());}
    double variance()const { return sum_sqr()*(1.0/count())-sqr(mean());}
    double chi_value()const { return mean()/sqrt(variance());}
    double t_value()const { return mean()/sqrt(variance()/count());}
    double T2_value()const { return sqr(t_value());}


    moments()=default;

    moments& operator+=(const moments& other){count_+=other.count(); sum_+=other.sum(); sum_sqr_+=other.sum_sqr(); return *this;}


    void push_back(double val)
    {
        count_++;
        sum_+=val;
        sum_sqr_+=sqr(val);
    }
    moments<double> operator*(double t)
    {
        double newcount=count()*t;
        double newsum=sum()*t;
        double newsumsqr=sum_sqr()*t;
        return moments<double> (newcount,newsum,newsumsqr);
    }

private:
    double count_;
    double sum_;
    double sum_sqr_;
};


template<typename T>
class moments<M_Matrix<T>>
{

public:
    typedef  moments self_type;
    constexpr static auto  className=my_static_string("moments_")+my_trait<M_Matrix<T>>::className;
    static auto get_constructor_fields()
    {
        return std::make_tuple(
            grammar::field(C<self_type>{},"count",&self_type::count),
            grammar::field(C<self_type>{},"sum",&self_type::sum),
            grammar::field(C<self_type>{},"sum_sqr",&self_type::sum_sqr)
                                                                                                                    );
    }

    moments(std::size_t count, M_Matrix<T> sum, M_Matrix<T> sum_sqr):count_(count),sum_{sum},sum_sqr_{sum_sqr},
                                                                       mean_{calc_mean()},variance_{calc_variance()},chi_value_{calc_chi_value()},t_value_{calc_t_value()}{};

    moments(const std::vector<M_Matrix<T>>& data):
                                                    count_(data.size()),sum_(op::sum(data)),sum_sqr_(op::sum_sqr(data)),
                                                    mean_{calc_mean()},variance_{calc_variance()},chi_value_{calc_chi_value()},t_value_{calc_t_value()}{}

    template<class C,class F>
    moments(const std::vector<C>& data, const F& f):
                                                      count_(data.size()),sum_(op::sum(data,f)),sum_sqr_(op::sum_sqr(data,f)),
                                                      mean_{calc_mean()},variance_{calc_variance()},chi_value_{calc_chi_value()},t_value_{calc_t_value()}{}


    std::size_t data_size()const { return size();}


    std::size_t data_dim(MOMENTS m)const {
        switch(m){
        case COUNT: return 0;
        case SUM: return 1;
        case SUM_SQR: return  2;
        case MEAN: return 1;
        case VARIANCE:return 2;
        case CHI_VALUE: return 1;
        case T_VALUE: return 1;
        case T2_VALUE:
        default:
            return 0;

        }
    }
    std::tuple<double>
    data_row(MOMENTS m, std::size_t i, std::size_t j)const {
        switch(m){
        case COUNT: return count();
        case SUM: return sum()[i];
        case SUM_SQR: return  sum_sqr()(i,j);
        case MEAN: return mean()[i];
        case VARIANCE:return variance()(i,j);
        case CHI_VALUE: return chi_value()[i];
        case T_VALUE: return t_value()[i];
        case T2_VALUE:
        default:
            return T2_value();

        }
    }

    auto count()const {return count_;}
    auto& sum()const { return sum_;}
    auto& sum_sqr()const {return sum_sqr_;}
    auto& mean()const {return mean_;}
    auto& variance()const { return variance_;}

    auto& chi_value()const {return chi_value_;}

    auto& t_value()const {return t_value_;}

    double T2_value()const {
        auto covinv=inv(variance()/count());
        if (!covinv) return std::numeric_limits<double>::quiet_NaN();
        else
            return xTSigmaX(mean(),covinv.value());
    }


    moments()=default;

    std::size_t size()const {return sum_.size();}
    moments& operator+=(const moments& other){count_+=other.count(); sum_+=other.sum(); sum_sqr_+=other.sum_sqr();
        mean_=calc_mean();variance_=calc_variance();chi_value_=calc_chi_value();t_value_=calc_t_value();}

private:
    std::size_t count_;
    M_Matrix<T> sum_;
    M_Matrix<T> sum_sqr_;
    M_Matrix<T> mean_;
    M_Matrix<T> variance_;
    M_Matrix<T> chi_value_;
    M_Matrix<T> t_value_;

    M_Matrix<T> calc_variance(){ return sum_sqr()*(1.0/count())-op::sqr(mean());}
    M_Matrix<T> calc_mean(){return sum()*(1.0/count());}
    M_Matrix<T> calc_chi_value(){ return mean()*(diag(variance())).apply([](double x){return 1.0/std::sqrt(x);});}
    M_Matrix<T> calc_t_value(){ return mean()*(diag(variance())/count()).apply([](double x){return 1.0/std::sqrt(x);});}
};

template<typename T>
class moments_matrix
{

public:
    typedef  moments_matrix self_type;
    constexpr static auto  className=my_static_string("moments_matrix_")+my_trait<M_Matrix<T>>::className;
    static auto get_constructor_fields()
    {
        return std::make_tuple(
            grammar::field(C<self_type>{},"count",&self_type::count),
            grammar::field(C<self_type>{},"sum",&self_type::sum),
            grammar::field(C<self_type>{},"sum_sqr",&self_type::sum_elem_sqr)
                                                                                                                    );
    }

    moments_matrix(std::size_t count, M_Matrix<T> sum, M_Matrix<T> sum_elem_sqr):count_(count),sum_{sum},sum_elem_sqr_{sum_elem_sqr},
                                                                                   mean_{calc_mean()},variance_{calc_variance()},chi_value_{calc_chi_value()},t_value_{calc_t_value()}
                                                                                   {};

    moments_matrix(const std::vector<M_Matrix<T>>& data):
                                                           count_(data.size()),sum_(op::sum(data)),sum_elem_sqr_(op::sum_sqr_elem(data)),
                                                           mean_{calc_mean()},variance_{calc_variance()},chi_value_{calc_chi_value()},t_value_{calc_t_value()}
    {}

    template<class C,class F>
    moments_matrix(const std::vector<C>& data, const F& f):
                                                             count_(data.size()),sum_(op::sum(data,f)),sum_elem_sqr_(op::sum_sqr_elem(data,f)),        mean_{calc_mean()},variance_{calc_variance()},chi_value_{calc_chi_value()},t_value_{calc_t_value()}
    {}


    std::size_t data_size()const { return size();}


    std::size_t data_dim(MOMENTS m)const {
        switch(m){
        case COUNT: return 0;
        case SUM: return 2;
        case SUM_SQR: return  2;
        case MEAN: return 2;
        case VARIANCE:return 2;
        case CHI_VALUE: return 2;
        case T_VALUE: return 2;
        case T2_VALUE:
        default:
            return 0;

        }
    }
    std::tuple<double>
    data_row(MOMENTS m, std::size_t i, std::size_t j)const {
        switch(m){
        case COUNT: return count();
        case SUM: return sum()(i,j);
        case SUM_SQR: return  sum_elem_sqr()(i,j);
        case MEAN: return mean()(i,j);
        case VARIANCE:return variance()(i,j);
        case CHI_VALUE: return chi_value()(i,j);
        case T_VALUE: return t_value()(i,j);
        case T2_VALUE:
        default:
            return T2_value();

        }
    }

    auto count()const {return count_;}
    auto& sum()const { return sum_;}
    auto& sum_elem_sqr()const {return sum_elem_sqr_;}
    auto& mean()const {return mean_;}
    auto& variance()const { return variance_;}

    auto& chi_value()const {return chi_value_;}

    auto& t_value()const {return t_value_;}

    double T2_value()const { return fullSum(elemMult(t_value(),t_value()));}


    moments_matrix()=default;

    std::size_t size()const {return sum_.size();}
    moments_matrix& operator+=(const moments_matrix& other){count_+=other.count(); sum_+=other.sum(); sum_elem_sqr_+=other.sum_sqr();
        mean_=calc_mean();variance_=calc_variance();chi_value_=calc_chi_value();t_value_=calc_t_value();}


private:
    std::size_t count_;
    M_Matrix<T> sum_;
    M_Matrix<T> sum_elem_sqr_;
    M_Matrix<T> mean_;
    M_Matrix<T> variance_;
    M_Matrix<T> chi_value_;
    M_Matrix<T> t_value_;

    M_Matrix<T> calc_variance(){ return sum_elem_sqr()*(1.0/count())-elemMult(mean(),mean());}
    M_Matrix<T> calc_mean(){return sum()*(1.0/count());}
    M_Matrix<T> calc_chi_value(){ return elemDiv(mean(),(variance()).apply([](auto x){return std::sqrt(x);}));}
    M_Matrix<T> calc_t_value(){ return elemDiv(mean(),(variance()/count()).apply([](auto x){return std::sqrt(x);}));}
};





template <class T,class =void >
struct getMoment{

    typedef moments<T> type;
};

template <class T>
struct getMoment<T, std::void_t<typename T::moment_type>>
{
    typedef typename T::moment_type type;
};

template <class C>
using moments_t=typename getMoment<C>::type;


template<class X,class F, typename... args>
auto moments_by_row(const std::vector<X>& x, const F& f, args...t)
{
    typedef std::invoke_result_t<F,X,args...> C;
    typedef typename std::decay_t<C>::value_type T;
    auto& c=std::invoke(f,x[0],t...);

    std::vector<moments_t<T>> out(c.size());
    for (std::size_t i=0; i<out.size(); ++i)
    {
        out[i]=moments_t<T>(x,[&f,&i,&t...](auto& e){return std::invoke(f,e,t...)[i];});
    }
    return out;
}

template<class X,class F, class G>
auto moments_by_row_2(const std::vector<X>& x, const F& f, const G& g)
{
    typedef std::invoke_result_t<F,X> C;
    typedef typename std::decay_t<C>::value_type T;
    typedef std::decay_t<std::invoke_result_t<G,T>> K;
    auto& c=std::invoke(f,x[0]);

    std::vector<moments_t<K>> out(c.size());
    for (std::size_t i=0; i<out.size(); ++i)
    {
        out[i]=moments_t<K>(x,[&f,&g,&i](auto& e){return std::invoke(g,std::invoke(f,e)[i]);});
    }
    return out;
}





template<typename T>
class Normal_Distribution;

class Beta_Distribution;



class stretch_move_Distribution;



template<typename T>
class Base_Distribution
{
public:
    virtual Base_Distribution* clone()const=0;
    virtual std::ostream& put(std::ostream& os)const=0;

    constexpr static auto const className=my_static_string("Base_Distribution_")+my_trait<T>::className;
    virtual std::string myClass()const=0;

    virtual T sample(std::mt19937_64& mt)const =0;

    virtual double p(const T& x)const =0;

    virtual double logP(const T& x)const =0;

    virtual double vlogP(const T& x)const =0;


    virtual M_Matrix<T> param()const =0;

    // virtual M_Matrix<T> Score(const T& x)const =0;    // V(x)

    virtual M_Matrix<T> Fisher_Information()const =0;    //I()


    virtual T dlogL_dx(const T& x)const =0;

    virtual T dlogL_dx2(const T& x)const =0;


    virtual T mean()const =0;

    virtual T stddev()const =0;

    virtual double expected_logP()const { return 0;}
    virtual double variance_logP()const { return 0.5;}


    virtual ~Base_Distribution(){}
};


template<bool output,typename X>
struct are_Equal<output,X,std::enable_if_t<std::is_base_of_v<Base_Distribution<double>,X >,void>>
{
private :
    double absolute_;
    double relative_;

public:
    are_Equal(double abs, double rel) : absolute_{abs},relative_{rel}{}
    are_Equal(){}
    bool test(const X& x, const X& y, std::ostream &os ) const
    {
        if (x.myClass()!=y.myClass())
        {
            os<<"Differ in Transformation :"<<x.myClass()<<" vs "<<y.myClass();
            return false;
        }
        else return are_Equal_v(x.param(),y.param(),absolute_,relative_,os);
    }
    bool test(const X* x, const X* y, std::ostream &os ) const
    {
        if (x->myClass()!=y->myClass())
        {
            os<<"Differ in Transformation :"<<x->myClass()<<" vs "<<y->myClass();
            return false;
        }
        else return are_Equal_v(x->param(),y->param(),absolute_,relative_,os);
    }


};



template<bool output,typename X>
struct are_Equal<output,X,std::enable_if_t<!is_field_Object<X>::value&&std::is_base_of_v<Base_Distribution<M_Matrix<double>>,X >,void>>
{
private :
    double absolute_;
    double relative_;

public:
    are_Equal(double abs, double rel) : absolute_{abs},relative_{rel}{}
    are_Equal(){}
    bool test(const X& x, const X& y, std::ostream &os ) const
    {
        if (x.myClass()!=y.myClass())
        {
            os<<"Differ in Transformation :"<<x.myClass()<<" vs "<<y.myClass();
            return false;
        }
        else return are_Equal_v(x.param(),y.param(),absolute_,relative_,os);
    }
    bool test(const X* x, const X* y, std::ostream &os ) const
    {
        if (x->myClass()!=y->myClass())
        {
            os<<"Differ in Transformation :"<<x->myClass()<<" vs "<<y->myClass();
            return false;
        }
        else return are_Equal_v(x->param(),y->param(),absolute_,relative_,os);
    }


};



template <>
struct Derived_types<Base_Distribution<double>>
{
    typedef   Cs<Normal_Distribution<double>,Beta_Distribution, stretch_move_Distribution> type;
    constexpr bool static value=true;
};
inline double logit(double x){return std::log(x/(1.0-x));}



inline std::pair<double,double> logit(double x,double sd){
    return {std::log(x/(1.0-x)),sd/(x*(1.0-x))};}


inline double logistic(double x){return 1.0/(1.0+std::exp(-x));}
inline double log10it(double x){return std::log10(x/(1.0-x));}

inline double dlog10it(double x){return 1.0/(std::log(10.0)*(x-sqr(x)));}


inline std::pair<double,double> log10it(double x,double sd){
    return {std::log10(x/(1.0-x)),sd/(x*(1.0-x))};}


inline double log10istic(double x){return 1.0/(1.0+std::pow(10.0,-x));}
inline double dlog10istic(double x){return -std::log(10)*std::pow(10,-x)/sqr(1.0+std::pow(10.0,-x));}


template <class T> class Identity_Transformation;
template <class T> class  Logarithm_Transformation;
template <class T> class  Logit_Transformation;


template<typename T>
class Base_Transformation
{
public:
    constexpr static auto const className=my_static_string("Base_Transformation_")+my_trait<T>::className;
    virtual std::string myClass()const=0;

    virtual Base_Transformation* clone()const =0;
    virtual std::ostream& put(std::ostream& os)const=0;

    virtual T apply(const T& x)const =0;
    virtual T apply_inv(const T& x)const =0;

    virtual T dapply(const T& x)const =0;
    virtual T dapply_inv(const T& x)const =0;


    virtual ~Base_Transformation(){}
};

template<bool output,typename X>
struct are_Equal<output,X,std::enable_if_t<std::is_base_of_v<Base_Transformation<double>,X >,void>>
{
public:
    are_Equal(double,double){}
    are_Equal(){}
    bool test(const X& x, const X& y, std::ostream &os ) const
    {
        if (x.myClass()!=y.myClass())
        {
            os<<"Differ in Transformation :"<<x.myClass()<<" vs "<<y.myClass();
            return false;
        }
        else return true;
    }
    bool test(const X* x, const X* y, std::ostream &os ) const
    {
        if (x->myClass()!=y->myClass())
        {
            os<<"Differ in Transformation :"<<x->myClass()<<" vs "<<y->myClass();
            return false;
        }
        else return true;
    }


};


struct s_logDetCov{    constexpr static auto const title = my_static_string("logDetCov");};
struct s_mean{    constexpr static auto const title = my_static_string("mean");};
struct s_variance{    constexpr static auto const title = my_static_string("variance");};
struct s_stddev{    constexpr static auto const title = my_static_string("stddev");};
struct s_num_samples{    constexpr static auto const title = my_static_string("num_samples");};

struct s_Cov{    constexpr static auto const title = my_static_string("Cov");};
struct s_Transpose{    constexpr static auto const title = my_static_string("_T");};

struct s_alpha{    constexpr static auto const title = my_static_string("alpha");};
struct s_beta{    constexpr static auto const title = my_static_string("beta");};
struct s_betaMap{    constexpr static auto const title = my_static_string("betaMap");};
struct s_probability{ constexpr static auto const title = my_static_string("probability");};
struct s_use{ constexpr static auto const title = my_static_string("use");};
struct s_obs{ constexpr static auto const title = my_static_string("obs");};
struct s_pred{ constexpr static auto const title = my_static_string("pred");};
struct s_logEvidence{ constexpr static auto const title = my_static_string("logEvidence");};
struct s_logL{ constexpr static auto const title = my_static_string("logL");};
struct s_elogL{ constexpr static auto const title = my_static_string("elogL");};
struct s_vlogL{    constexpr static auto const title = my_static_string("vlogL");};
struct s_evlogL{    constexpr static auto const title = my_static_string("evlogL");};



template<>
struct Derived_types<Base_Transformation<double>>{
    typedef Cs<Identity_Transformation<double>,Logarithm_Transformation<double>, Logit_Transformation<double>> type;
    constexpr bool static value=true;

};

template<typename T>
class Identity_Transformation: public Base_Transformation<T>
{
public:
    typedef  Base_Transformation<T> base_type;
    typedef T value_type;
    static std::tuple<> get_constructor_fields() {return std::tuple<>();}
    virtual std::string myClass()const override{return className.str();}
    virtual std::ostream& put(std::ostream& os)const override
    {
        return io::output_operator_on_Object (os,*this);
    }


    constexpr static auto const className=my_static_string("Identity_Transformation")+my_trait<T>::className;

    virtual Identity_Transformation* clone()const override{ return new Identity_Transformation(*this);};

    virtual T apply(const T& x)const override {return x;}
    virtual T apply_inv(const T& x)const override {return x;}
    virtual T dapply(const T& )const override {return T(1);}
    virtual T dapply_inv(const T& )const override {return T(1);}

    virtual ~Identity_Transformation(){}

    Identity_Transformation(){}
};



template<typename T>
class Logarithm_Transformation;

template<>
class Logarithm_Transformation<double>: public Base_Transformation<double>
{
public:
    typedef  Base_Transformation<double> base_type;
    static std::tuple<> get_constructor_fields() {return std::tuple<>();}
    virtual std::string myClass()const override{return className.str();}
    virtual std::ostream& put(std::ostream& os)const override
    {
        return io::output_operator_on_Object (os,*this);
    }

    constexpr static auto const className=my_static_string("Logarithm_Transformation");

    virtual Logarithm_Transformation* clone()const override{ return new Logarithm_Transformation(*this);};
    virtual double apply(const double& x)const override  {return std::log10(x);}
    virtual double apply_inv(const double& x)const override {return std::pow(10.0,x);}
    virtual double dapply(const double& x)const override {return 1.0/(std::log(10.0)*x);}
    virtual double dapply_inv(const double& x)const override {return std::log(10)*std::pow(10.0,x);}

    virtual ~Logarithm_Transformation(){}

    Logarithm_Transformation() {}

};

template<>
class Logit_Transformation<double>: public Base_Transformation<double>
{
public:
    typedef  Base_Transformation<double> base_type;
    static std::tuple<> get_constructor_fields() {return std::tuple<>();}
    virtual std::string myClass()const override{return className.str();}
    virtual std::ostream& put(std::ostream& os)const override
    {
        return io::output_operator_on_Object (os,*this);
    }

    constexpr static auto const className=my_static_string("Logit_Transformation");

    virtual Logit_Transformation* clone()const override{ return new Logit_Transformation(*this);};
    virtual double apply(const double& x)const override  {return log10it(x);}
    virtual double apply_inv(const double& x)const override {return log10istic(x);}
    virtual double dapply(const double& x)const override  {return dlog10it(x);}
    virtual double dapply_inv(const double& x)const override {return dlog10istic(x);}

    virtual ~Logit_Transformation(){}

    Logit_Transformation() {}

};

template<template <typename...> class  My_vec,typename _IntType>
class multinomial_distribution // not a child of  Base_Distribution, modeled after STL random library

{
    static_assert(std::is_integral<_IntType>::value,
                  "template argument not an integral type");

public:
    /** The type of the range of the distribution. */
    typedef My_vec<_IntType> result_type;
    /** Parameter type. */

    struct param_type
    {
        typedef multinomial_distribution<My_vec,_IntType> distribution_type;
        friend class multinomial_distribution<My_vec,_IntType>;

        explicit
            param_type(_IntType _N, My_vec<double> P )
            : N_(_N), P_(P)
        {
            _M_initialize();
        }

        _IntType
        N() const
        { return N_; }

        const My_vec<double>&
        P() const
        { return P_; }

        friend bool
        operator==(const param_type& __p1, const param_type& __p2)
        { return __p1.N_ == __p2.N_ && __p1.P_ == __p2.P_; }

    private:
        void
        _M_initialize()
        {
            std::size_t n=P_.size();
            rP_=P_-P_;
            auto s=rP_;
            s[n-1]=P_[n-1];
            for (std::size_t i=1; i<P_.size(); ++i)
            {
                s[n-1-i]=P_[n-1-i]+s[n-i];
            }
            for (std::size_t i=0; i<P_.size(); ++i)
            {
                rP_[i]=P_[i]/s[i];
                P_[i]=P_[i]/s[0];
            }
        }

        _IntType N_;
        My_vec<double> P_;
        My_vec<double> rP_;

    };



    // constructors and member function
    explicit
        multinomial_distribution(_IntType __t,
                                 My_vec<double> __p )
        : _M_param(__t, __p)
    { }



    /**
     * @brief Returns the distribution @p t parameter.
     */
    _IntType
    N() const
    { return _M_param.N(); }

    /**
     * @brief Returns the distribution @p p parameter.
     */
    const My_vec<double>&
    P() const
    { return _M_param.P(); }

    /**
     * @brief Returns the parameter set of the distribution.
     */
    param_type
    param() const
    { return _M_param; }

    /**
     * @brief Sets the parameter set of the distribution.
     * @param __param The new parameter set of the distribution.
     */
    void
    param(const param_type& __param)
    { _M_param = __param; }

    /**
     * @brief Returns the greatest lower bound value of the distribution.
     */
    _IntType
    min() const
    { return 0; }

    /**
     * @brief Returns the least upper bound value of the distribution.
     */
    _IntType
    max() const
    { return _M_param.N(); }

    /**
     * @brief Generating functions.
     */
    template<typename _UniformRandomNumberGenerator>
    result_type
    operator()(_UniformRandomNumberGenerator& __urng)
    { return this->operator()(__urng, _M_param); }

    template<typename _UniformRandomNumberGenerator>
    result_type
    operator()(_UniformRandomNumberGenerator& __urng,
               const param_type& __p)
    {
        result_type out(emptyCopy<_IntType>(__p.P_));
        _IntType Nr=__p.N_;
        std::binomial_distribution<_IntType> Bi_;
        typedef typename std::binomial_distribution<_IntType>::param_type biPar;
        for (std::size_t i=0; i< out.size()-1; ++i)
        {
            _IntType ni=Bi_(__urng,biPar(Nr,__p.rP_[i]));
            Nr-=ni;
            out[i]=ni;
        }
        out[out.size()-1]=Nr;
        return out;
    }

    static
        double logP(const result_type& x,const param_type& __p)
    {
        assert(sum(x)==__p.N());
        assert(x.size()==__p.P().size());
        double out=lgamma(__p.N()+1);
        for (std::size_t i=0; i<x.size(); ++i)
            out+= x[i]*log(__p.P()[i])-lgamma(x[i]+1);
        return out;

    }
    double logP(const result_type& x)
    {
        return logP(x,_M_param);
    }
    double P(const result_type& x)
    {
        return std::exp(logP(x));
    }

    static
        double P(const result_type& x,const param_type& __p)
    {
        return std::exp(logP(x,__p));
    }



private:

    param_type _M_param;

};

template <typename T,typename _UniformRandomNumberGenerator>
T sample_rev_map(const std::map<double,T>& reverse_prior,_UniformRandomNumberGenerator& mt);

template<typename T>
class markov_process;



template<typename T>
class markov_process<M_Matrix<T>>

{

public:
    /*  * The type of the range of the distribution. */
    typedef M_Matrix<T> result_type;
    /** Parameter type. */

    markov_process()=default;
    struct param_type
    {
        typedef markov_process<M_Matrix<T>> distribution_type;
        friend class markov_process<M_Matrix<T>>;

        param_type()=default;
        explicit
            param_type(const M_Matrix<T>& _N, const M_Matrix<double>& P )
            : N_(_N), P_(P),cP_(M_Matrix<double>(nrows(P),ncols(P)))

        {
            _M_initialize();
        }

        const M_Matrix<T>&
        N() const
        { return N_; }

        void set_P(M_Matrix<double>&&
                       _P)
        {
            P_=std::move(_P);
            _M_initialize();
        }

        void set_N(M_Matrix<T>&& _N)
        {  N_=std::move(_N); }

        const M_Matrix<double>&
        P() const
        { return P_; }

        friend bool
        operator==(const param_type& __p1, const param_type& __p2)
        { return __p1.N_ == __p2.N_ && __p1.P_ == __p2.P_; }

    private:
        void
        _M_initialize()
        {
            std::size_t n=ncols(P_);

            auto s=M_Matrix<double>(nrows(P_),ncols(P_));

            for (std::size_t i=0; i<nrows(P_); ++i)
            {
                s(i,n-1)=P_(i,n-1);
                for (std::size_t j=1; j<ncols(P_); ++j)
                {
                    s(i,n-1-j)=P_(i,n-1-j)+s(i,n-j);
                }
                for (std::size_t j=0; j<ncols(P_); ++j)
                {
                    if (s(i,j)>0)
                        cP_(i,j)=P_(i,j)/s(i,j);
                    else cP_(i,j)=0;
                    P_(i,j)=P_(i,j)/s(i,0);
                }
            }
        }

        M_Matrix<T> N_;
        M_Matrix<double> P_;
        M_Matrix<double> cP_;
    };


    // constructors and member function
    explicit
        markov_process(const M_Matrix<T>& __t,
                       const M_Matrix<double>& __p)
        : _M_param(__t, __p)
    { }



    /**
       * @brief Returns the distribution @p t parameter.
       */
    const M_Matrix<T>&
    N() const
    { return _M_param.N(); }




    /**
       * @brief Returns the distribution @p p parameter.
       */
    const M_Matrix<double>&
    P() const
    { return _M_param.P(); }


    void set_P(M_Matrix<double>&& _P)
    { _M_param.set_P(std::move(_P));
    }

    void set_N(M_Matrix<T>&& _N)
    {
        _M_param.set_N(std::move(_N));
    }


    /**
       * @brief Returns the parameter set of the distribution.
       */
    param_type
    param() const
    { return _M_param; }

    /**
       * @brief Sets the parameter set of the distribution.
       * @param __param The new parameter set of the distribution.
       */
    void
    param(const param_type& __param)
    { _M_param = __param; }

    /**
       * @brief Returns the greatest lower bound value of the distribution.
       */
    T
    min() const
    { return 0; }

    /**
       * @brief Returns the least upper bound value of the distribution.
       */
    T
    max() const
    { return totalsum(_M_param.N()); }

    /**
       * @brief Generating functions.
       */
    template<typename _UniformRandomNumberGenerator>
    result_type
    operator()(_UniformRandomNumberGenerator& __urng)
    { return this->operator()(__urng, _M_param); }

    template<typename _UniformRandomNumberGenerator>
    result_type
    operator()(_UniformRandomNumberGenerator& __urng,
               const param_type& _p, bool sumcols=true)
    {
        if (sumcols)
        {
            std::size_t nc=ncols(_p.P());
            M_Matrix<T> out(1,nc,0);
            for (std::size_t i=0; i< nrows(_p.P()); ++i)
            {
                T Nr=_p.N_[i];

                for (std::size_t j=0; j< nc-1; ++j)
                {
                    double p=_p.cP_(i,j);
                    auto bi=std::binomial_distribution<T>(Nr,p);
                    T n=bi(__urng);
                    Nr-=n;
                    out(0,j)+=n;
                }
                out(0,nc-1)+=Nr;
            }
            return out;
        }
        else
        {
            std::size_t nc=ncols(_p.P_);
            std::size_t nr=nrows(_p.P_);
            M_Matrix<T> out(nr,nc);
            for (std::size_t i=0; i< nrows(_p.P_); ++i)
            {
                T Nr=_p.N_[i];

                for (std::size_t j=1; j< nc-1; ++j)
                {
                    double p=_p.cP_(i,j);
                    auto bi=std::binomial_distribution<T>(Nr,p);
                    T n=bi(__urng);
                    Nr-=n;
                    out(i,j)=n;
                }
                out(i,nc-1)+=Nr;
            }
            return out;
        }

    }


private:

    param_type _M_param;

    };


    class white_noise
    {
    public:
        double operator()(std::mt19937_64& mt, double dt)const
        {
            return std::normal_distribution<double>(0,std::sqrt(noise_at_1_Hz_/dt))(mt);
        }
        double variance(double dt)const
        {
            return noise_at_1_Hz_/dt;
        };

        white_noise(double noise_at_1Hz):noise_at_1_Hz_{noise_at_1Hz}{}
    private:
        double noise_at_1_Hz_;
    };



    M_Matrix<double> normalize_to_prob(M_Matrix<double>&& P)
    {
        if (P.nrows()==1)
        {
            double sum=0;
            for (std::size_t i=0; i<P.ncols(); ++i)
                sum+=std::abs(P(0,i));
            for (std::size_t i=0; i<P.ncols(); ++i)
                P(0,i)=std::abs(P(0,i))/sum;
        }
        else if (P.ncols()==1)
        {      double sum=0;
            for (std::size_t i=0; i<P.nrows(); ++i)
                sum+=std::abs(P(i,0));
            for (std::size_t i=0; i<P.nrows(); ++i)
                P(i,0)=std::abs(P(i,0))/sum;
        }
        else
        {
            for (std::size_t j=0; j<P.ncols();++j)
            {
                double sum=0;
                for (std::size_t i=0; i<P.nrows(); ++i)
                    sum+=std::abs(P(i,j));
                for (std::size_t i=0; i<P.nrows(); ++i)
                    P(i,j)=std::abs(P(i,j))/sum;

            }
        }
        return std::move(P);
    }



    template<typename T>
    class Normal_Distribution;


    template<>
    class Normal_Distribution<double>: public Base_Distribution<double>
    {
    public:
        typedef  Base_Distribution<double> base_type;

        constexpr static auto const className=my_static_string("Normal_Distribution");
        std::string myClass()const override { return className.str();}

        typedef   Normal_Distribution<double> self_type ;
        static auto get_constructor_fields()
        {
            return std::make_tuple(
                grammar::field(C<self_type>{},"mean",&self_type::mean),
                grammar::field(C<self_type>{},"variance",&self_type::variance)
                                                                                                                                            );
        }

        virtual std::ostream& put(std::ostream& os)const override
        {
            return io::output_operator_on_Object (os,*this);
        }

        virtual Normal_Distribution<double>* clone()const override{ return new Normal_Distribution<double>(*this);};

        virtual double sample(std::mt19937_64& mt) const override{ return std::normal_distribution<double> {param_[0],param_[1]}(mt);}

        virtual double p(const double& x)const override { return 1.0/(std::sqrt(2*PI)*stddev())*std::exp(-0.5*sqr((x-mean())/stddev()));}

        virtual double logP(const double& x)const override {return -0.5*std::log(2*PI)-std::log(stddev())-0.5*sqr((x-mean())/stddev());}
        virtual double vlogP(const double& x)const override {return sqr(logP(x)-expected_logP());}


        /*
    virtual M_Matrix<double> Score(const double& x) const override
    { return M_Matrix<double>(1,2,std::vector<double>{dlogL_dmean(x),dlogL_dstddev(x)});};
*/
        virtual  M_Matrix<double> Fisher_Information()const override{
            return M_Matrix<double>(2,2,Matrix_TYPE::DIAGONAL,std::vector<double>{d2logL_dmean2(),d2logL_dvariance2()});};

        virtual double mean()const override {return param_[0];}

        virtual double stddev()const override  {return std::sqrt(param_[1]);};

        virtual double variance()const {return param_[1];}

        virtual M_Matrix<double> param()const override{return param_;}

        Normal_Distribution<double>(double mean, double variance)
            :param_(1,2,std::vector<double>{mean,variance}){}

        Normal_Distribution()=default;
        virtual ~Normal_Distribution(){}
        virtual double dlogL_dx(const double& x)const override
        {
            return (mean()-x)/variance();
        }

        virtual double dlogL_dx2(const double& )const override
        {
            return -1.0/variance();
        }


        virtual double expected_logP()const override
        {
            return -0.5*(1.0+log(2*PI*variance()));
        }

        virtual double variance_logP()const override
        {
            return 0.5;
        }

        Normal_Distribution(const Normal_Distribution&)=default;
        Normal_Distribution(Normal_Distribution&&)=default;
        Normal_Distribution&operator=(const Normal_Distribution&)=default;
        Normal_Distribution&operator=(Normal_Distribution&&)=default;


        static auto
        get_data_index_static()
        {
            return std::make_tuple(
                make_data_static(std::tuple<>(),
                                 std::make_tuple(
                                         F_s(s_mean{},&self_type::mean),
                                         F_s(s_stddev{},&self_type::stddev),
                                         F_s(s_variance{},&self_type::variance))
));
        }




    protected:
        M_Matrix<double> param_;


        double d2logL_dmean2()const {return 1.0/variance();}

        double d2logL_dvariance2()const {return  0.5/sqr(variance());}

    };



    template<typename E>
    class Normal_Distribution<M_Matrix<E>>: public Base_Distribution<M_Matrix<E>>
    {
    public:
        constexpr static auto const className=my_static_string("Normal_Distribution_")+my_trait<M_Matrix<E>>::className;
        std::string myClass()const override { return className.str();}
        virtual Normal_Distribution* clone()const override{ return new Normal_Distribution(*this);}

        typedef   Normal_Distribution<M_Matrix<E>> self_type ;
        static auto get_constructor_fields()
        {
            return std::make_tuple(
                grammar::field(C<self_type>{},"mean",&self_type::mean),
                grammar::field(C<self_type>{},"Cov",&self_type::Cov),
                grammar::field(C<self_type>{},"CovInv",&self_type::CovInv),
                grammar::field(C<self_type>{},"cholesky",&self_type::Chol)
                                                                                                                                            );
        }
        virtual std::ostream& put(std::ostream& os)const override
        {
            return io::output_operator_on_Object (os,*this);
        }


        M_Matrix<E> sample(std::mt19937_64& mt)const override
        {
            M_Matrix<E> r;
            std::normal_distribution<> normal;
            if (param_.size()>0)
            {
                auto z=Rand(mean(),normal,mt);
                r=mean()+multTransp(z,Chol());
            }
            return r;
        }


        virtual  M_Matrix<M_Matrix<E>> Fisher_Information()const override{
            return M_Matrix<M_Matrix<E>>(2,2,Matrix_TYPE::DIAGONAL,
                                         std::vector<M_Matrix<E>>{d2logL_dmean2(),d2logL_dCov2()});}


        double logP(const M_Matrix<E> &x) const override
        {
            if (param_.size()>0)
                return -0.5*(mean().size()*log(PI)+logDetCov()+chi2(x));
            else
                return std::numeric_limits<double>::quiet_NaN();
        }
        double vlogP(const M_Matrix<E> &x) const override {
            return sqr(logP(x)-expected_logP());
        }

        double p(const M_Matrix<E>& x)const override
        {
            return std::exp(logP(x));
        }

        virtual M_Matrix<M_Matrix<E>>  param()const override {return param_;}

        const M_Matrix<E>& CovInv()const  {   return covinv_;  }

        const M_Matrix<E>& Chol()const{return cho_cov_;}

        double logDetCov()const {return logDetCov_;}


        virtual ~Normal_Distribution(){}

        M_Matrix<E> mean()const override {return param_[0];}
        M_Matrix<E> Cov()const {return param_[1];}
        M_Matrix<E> stddev()const override{return cho_cov_;}

        double chi2(const M_Matrix<E> &x) const
        {
            if (!param_.empty())
                return xTSigmaX(x-mean(),covinv_);
            else return std::numeric_limits<double>::quiet_NaN();
        }


        static    myOptional_t<Normal_Distribution> make(const M_Matrix<E> &mean,
                                                      const M_Matrix<E>& cov)
        {
            typedef myOptional_t<Normal_Distribution> Op;
            auto covinv=inv(cov);
            auto cho_cov=chol(cov,"lower");
            if (covinv.has_value()&&cho_cov.has_value())
                return Op(Normal_Distribution(mean,cov,std::move(covinv).value(),std::move(cho_cov).value()));
            else
                return Op(false,"singular matrix: inverse error:"+covinv.error()+ " and cholesky error: "+cho_cov.error());
        }

        Normal_Distribution(const M_Matrix<E> &mean,
                            const M_Matrix<E>& cov,
                            const M_Matrix<E>& covInv,
                            const M_Matrix<E>& cho):
                                                      param_(M_Matrix<M_Matrix<E>>(1,2,std::vector<M_Matrix<E>>{mean,cov})),
                                                      covinv_(covInv),
                                                      cho_cov_(cho),
                                                      logDetCov_(logDiagProduct(cho_cov_))
        {
            assert(cov.isSymmetric());
            assert(Cov().isSymmetric());
            assert(covinv_.isSymmetric());

            assert(!cho_cov_.empty());

        }


        Normal_Distribution(M_Matrix<E> &&mean,
                            M_Matrix<E>&& cov,
                            M_Matrix<E>&& covInv,
                            M_Matrix<E>&& cho):
                                                 param_(M_Matrix<M_Matrix<E>>(1,2,std::vector<M_Matrix<E>>{std::move(mean),std::move(cov)})),
                                                 covinv_(std::move(covInv)),
                                                 cho_cov_(std::move(cho)),
                                                 logDetCov_(logDiagProduct(cho_cov_))
        {
            assert(cov.isSymmetric());
            assert(Cov().isSymmetric());
            assert(covinv_.isSymmetric());

            assert(!cho_cov_.empty());

        }

        Normal_Distribution()=default;
        const M_Matrix<E>& d2logL_dmean2()const
        {
            return CovInv();
        }

        M_Matrix<E> d2logL_dCov2()const
        {
            std::size_t n=CovInv().size();
            M_Matrix<E> out(n,n,Matrix_TYPE::DIAGONAL);
            for (std::size_t i=0; i<n; ++i)
                for (std::size_t j=0; j<i+1; ++j)
                {
                    out(i,j)=0.5*CovInv()[i]*CovInv()[j];
                }

            return out;

        }

        void autoTest(std::mt19937_64& mt,std::size_t n)const
        {
            //  std::cerr<<"\n chi test n="<<mean().size()<<" chis\n";
            double chisum=0;
            double chisqr=0;
            for (std::size_t i=0; i<n; ++i)
            {
                auto s=sample(mt);
                auto chi=chi2(s);
                //  std::cerr<<chi<<" ";
                chisum+=chi;
                chisqr+=chi*chi;
            }
            chisum/=n;
            chisqr-=n*chisum*chisum;
            chisqr/=(n-1);
            std::cerr<<"\n chi test n="<<mean().size()<<" chimean="<<chisum<<" chisqr="<<chisqr<<"\n";
        }
        bool isValid()const
        {
            return !cho_cov_.empty();
        }

        virtual M_Matrix<E> dlogL_dx(const M_Matrix<E>& x)const override
        {
            return CovInv()*(mean()-x);
        }

        virtual M_Matrix<E> dlogL_dx2(const M_Matrix<E>& )const override
        {
            return -CovInv();
        }
        virtual double expected_logP()const override
        {
            return -0.5*(mean().size()*(log(PI)+1.0)+logDetCov());
        }

        virtual double variance_logP()const override
        {
            return 0.5*mean().size();
        }

        template<class i_index>
        static auto
        get_data_index_static(i_index)
        {
                return std::make_tuple(
                make_data_static(std::tuple<>(),
                                 std::make_tuple(
                                     F_s(s_logDetCov{},&self_type::logDetCov))
                                                                              ),

                    make_data_static(std::make_tuple(
                        I_s(i_index{},[](const self_type& s){return s.mean().size();},0)),
                                     std::make_tuple(
                                         F_s(s_mean{},[](const self_type& s, std::size_t i){return s.mean()[i];}))
                                         ),
                    make_data_static(std::make_tuple(
                                         I_s(i_index{},[](const self_type& s){return s.Cov().nrows();},0),
                                     I_s(CT_s<i_index,s_Transpose>{},[](const self_type& s){return s.Cov().ncols();},0)),
                                     std::make_tuple(
                                         F_s(s_Cov{},[](const self_type& s, std::size_t i, std::size_t j)
                                                {return s.Cov()(i,j);}))
                                     ));

        }



    protected:
        M_Matrix<M_Matrix<E>> param_;
        M_Matrix<E> covinv_;
        M_Matrix<E> cho_cov_;
        double logDetCov_=std::numeric_limits<double>::quiet_NaN();
    };


    inline double BetaDistribution(double p, std::size_t success, std::size_t failures)
    {
        return std::pow(p,success)*std::pow(1.0-p,failures)/std::exp(log_beta_f(1.0+success,1.0+failures));
    }



    class Beta_Distribution:public Base_Distribution<double>
    {
    public:
        typedef  Base_Distribution<double> base_type;
        constexpr static auto const className=my_static_string("Beta_Distribution");
        std::string myClass()const override { return className.str();}
        typedef   Beta_Distribution self_type ;
        static auto get_constructor_fields()
        {
            return std::make_tuple(
                grammar::field(C<self_type>{},"alfa",&self_type::alfa),
                grammar::field(C<self_type>{},"beta",&self_type::beta)
                                                                                                                                            );
        }

        static Data_Index_scheme data_index(const std::string variable_name)
        {
            Data_Index_scheme out;
            out.push_back(variable_name+"_alfa", {});
            out.push_back(variable_name+"_beta", {});

            return out;

        }

        static auto
        get_data_index(const std::string variable_name)
        {
            using namespace std::literals::string_literals;
            return std::make_tuple(
                Data_Index(Cs<self_type>(),variable_name+"_alfa"s,&self_type::alfa),
                Data_Index(Cs<self_type>(),variable_name+"_beta"s,&self_type::beta)
                                                                        );
        }

        static auto
        get_data_index_static()
        {
            return std::make_tuple(
                make_data_static(std::tuple<>(), std::make_tuple(
                                                   F_s(s_alpha(),&self_type::alfa),
                                                   F_s(s_beta(),&self_type::beta),
                                                   F_s(s_mean(),&self_type::mean))));
        }

        virtual Base_Distribution* clone()const override{ return new Beta_Distribution(*this);}

        Beta_Distribution(double alfa, double beta):
                                                      a_(1,2,std::vector<double>{alfa,beta}){}

        Beta_Distribution():Beta_Distribution(0.5,0.5){}


        double count()const {return alfa()+beta();}

        static Beta_Distribution UniformPrior()
        {
            return Beta_Distribution(1.0,1.0);
        }
        static Beta_Distribution UnInformativePrior()
        {
            return Beta_Distribution(0.5,0.5);
        }
        virtual std::ostream& put(std::ostream& os)const override
        {
            return os<<*this;
        }


        double p()const {return mean();}

        void push_accept()
        {
            ++a_[0];
        }
        void push_reject()
        {
            ++a_[1];
        }


        double sample(std::mt19937_64& mt)const override
        {
            std::gamma_distribution<double> ga(alfa(),2.0);
            std::gamma_distribution<double> gb(beta(),2.0);

            double a=ga(mt);
            double b=gb(mt);
            return a/(a+b);
        }


    private:
        M_Matrix<double> a_;


        // Base_Distribution interface
    public:
        virtual double p(const double &x) const override{
            return std::exp(logP(x));
        }
        virtual double logP(const double &x) const override
        {
            return alfa()*std::log(x)+beta()*std::log(1.0-x)-log_beta_f(1.0+alfa(),1.0+beta());
        }
        virtual double vlogP(const double &x) const override
        {
            return sqr(logP(x)-expected_logP());
        }

        double alfa()const{return a_[0];}
        double beta()const{return a_[1];}

        virtual  M_Matrix<double> param() const override{ return a_;}
        virtual M_Matrix<double> Fisher_Information() const override
        { return M_Matrix<double>(2,2,std::vector<double>{d2logLik_dalfa2(),d2logLik_dalfadbeta(),d2logLik_dalfadbeta(),d2logLik_dbeta2()});}

        virtual double mean() const override { return alfa()/(alfa()+beta());}
        virtual double stddev() const override { return std::sqrt(variance());}
        double variance() const {return alfa()*beta()/(sqr(alfa()+beta())*(alfa()+beta()+1));}

        double d2logLik_dalfadbeta()const { return -digamma(alfa()+beta());}

        double d2logLik_dalfa2()const { return digamma(alfa())+d2logLik_dalfadbeta();}
        double d2logLik_dbeta2()const { return digamma(beta())+d2logLik_dalfadbeta();}

        virtual double dlogL_dx(const double& x)const override
        {
            return alfa()/x-beta()/(1.0-x);

        }

        virtual double dlogL_dx2(const double& x)const override
        {
            return -alfa()/sqr(x)+ beta()/sqr(1-x);

        }


    };



    class stretch_move_Distribution: public Base_Distribution<double>
    {


        // Base_Distribution interface
    public:
        virtual stretch_move_Distribution* clone()const override{ return new stretch_move_Distribution(*this);};
        constexpr static auto const className=my_static_string("stretch_move_Distribution");
        std::string myClass()const override { return className.str();}
        typedef   stretch_move_Distribution self_type ;
        static auto get_constructor_fields()
        {
            return std::make_tuple(
                grammar::field(C<self_type>{},"alpha_stretch",&self_type::alpha)
                                                                                                                    );
        }
        virtual std::ostream& put(std::ostream& os)const override
        {
            return os<<*this;
        }

        static Data_Index_scheme data_index()
        {
            Data_Index_scheme out;
            out.push_back("alpha_stretch", {});
            return out;
        }
        static auto
        get_data_index()
        {
            using namespace std::literals::string_literals;
            return std::make_tuple(
                Data_Index(Cs<self_type>(),"alpha_stretch"s,&self_type::alpha)
                                                                                    );
        }
        static auto
        get_data_index_static()
        {
            using namespace std::literals::string_literals;
            return std::make_tuple(
                make_data_static(std::tuple<>(),std::make_tuple(F_s(s_alpha{},&self_type::alpha))
                                     ));
        }

        stretch_move_Distribution()=default;
        virtual double sample(std::mt19937_64 &mt) const override
        {
            std::uniform_real_distribution<double> U(std::sqrt(1.0/a_[0]),std::sqrt(a_[0]));
            double z=sqr(U(mt));
            return z;
        }
        virtual double p(const double &x) const override{ return std::pow(x,-0.5)/Z_;}
        virtual double logP(const double &x) const override { return -0.5*std::log(x)-log(Z_);}
        virtual double vlogP(const double& x)const override {return sqr(logP(x)-expected_logP());}

        virtual  M_Matrix<double> param() const override {return a_;}
        virtual M_Matrix<double> Fisher_Information() const override {
            assert(false);
            return {};}
        virtual double mean() const override {assert(false); return{};};
        virtual double stddev() const override {assert(false); return{};};
        stretch_move_Distribution(double a): a_{1,1,std::vector<double>{a}}, Z_(2*(std::sqrt(a)-std::pow(a,-0.5))) {}
        double alpha()const {return a_[0];}

    private:
        /// parameter alpha between 1 an infinity, 2 is good.
        M_Matrix<double> a_;
        double Z_;

        // Base_Distribution interface
    public:
        virtual double dlogL_dx(const double &x) const override
        {
            return -0.5/x;
        }
        virtual double dlogL_dx2(const double &x) const override{
            return 0.5/sqr(x);
        }
    };


    template<class P>
    struct complement_prob
    {
        complement_prob(const P& p):p_{p}{}
        template<typename... Ts>
        double operator()(Ts... xs)const
        {
            return 1.0-p_(xs...);
        }
    private:
        const P& p_;
    };

    template<class P>
    complement_prob<P>
    Complement_prob(const P& p){return complement_prob<P>(p);}


    template<class P>
    struct log_of
    {
        log_of(const P& p):p_{p}{}
        template<typename... Ts>
        double operator()(Ts... xs)const
        {
            return std::log(p_(xs...));
        }
    private:
        const P& p_;
    };


    template<class P>
    log_of<P> Log_of(const P& p){return log_of<P>(p);}

    template<class P>
    struct exp_of
    {
        exp_of(const P& p):p_{p}{}
        template<typename... Ts>
        double operator()(Ts... xs)const
        {
            return std::exp(p_(xs...));
        }
    private:
        const P& p_;
    };



    template <class T>
    std::pair<std::map<T,double>,double>
    logLik_to_p(const std::map<T,double>& logLikelihoods)
    {
        std::map<T,double> out(logLikelihoods);
        double Evidence=0;
        double maxlog=out.begin()->second;
        for (auto& e:out)
        {
            if (e.second>maxlog) maxlog=e.second;
        }
        for (auto& e:out)
        {
            e.second=std::exp(e.second-maxlog);
            Evidence+=e.second;
        }

        for (auto& e:out) e.second/=Evidence;

        return {out,Evidence};
    }

    template <typename T,typename _UniformRandomNumberGenerator>
    T sample_rev_map(const std::map<double,T>& reverse_prior,_UniformRandomNumberGenerator& mt)
    {
        std::uniform_real_distribution<> u;
        double r= u(mt);
        auto it=reverse_prior.lower_bound(r);
        return it->second;

    }


    template< typename T, class F, class Likelihood,class P_map>
    double Expectance(const F& f, const Likelihood& lik,const P_map& pm, const T& landa  )
    {
        auto p=pm;
        double sum=0;
        for (auto& e:p)
            sum+=e.second*lik(e.first,landa)*f(landa);
        return sum;
    }
    struct Probability: public invariant{};


    struct variable_value: public invariant
    {
        template<bool output>
        static Op_void test(double value)
        {
            if (std::isfinite(value))
                return Op_void(true,"");
            else if constexpr (output)
                return Op_void(false," not finite value: "+std::to_string(value));
            else
                return Op_void(true,"");
        }
    };
    struct variance_value: public invariant
    {
        template<bool output>
        static Op_void test(double value, double min_variance,double tolerance)
        {
            if (!std::isfinite(value))
                return Op_void(false," not finite value: "+std::to_string(value));
            else if (value+tolerance>min_variance)
                return Op_void(true,"");
            else if (value<0)
                return Op_void(false,"negative variance!! value="+std::to_string(value));
            else if (value==0)
                return Op_void(false,"zero variance!! value="+std::to_string(value));
            else    return Op_void(false,"variance too small!! value="+std::to_string(value)+ "min variance="+std::to_string(min_variance)+ "tolerance="+std::to_string(tolerance));

        }


        static double adjust(double value, double min_variance)
        {
            return std::max(value,min_variance);
        }

    };

    struct logLikelihood_value: public invariant
    {
        template<bool output>
        static Op_void test(double value)
        {
            if (std::isfinite(value))
                return Op_void(true,"");
            else
                return Op_void(false," not finite value: "+std::to_string(value));
        }
    };

    struct Probability_value: public Probability
    {
        template <bool output,class C>
        static Op_void test(const C& e, double tolerance)
        {
            if (!std::isfinite(e))
            {    if constexpr(output)
                    return  Op_void(false," not finite value="+ToString(e)+"\n");
                else return Op_void(false,"");
            }
            else if (e+tolerance<0)
            {
                if constexpr(output)
                    return  Op_void(false," negative prob="+ToString(e)+"\n");
                else return Op_void(false,"");
            }
            else if (e-tolerance>1)
            {
                if constexpr(output)
                    return  Op_void(false,"  prob greater than one" +ToString(e)+ " 1- prob="+ToString(1-e)+"\n");
                else return Op_void(false,"");

            }
            else return Op_void(true,"");
        }
    };



    struct Probability_distribution: public Probability
    {
        template<bool output,class P>
        static myOptional_t<void> test(const P& p,double tolerance)
        {
            double sum=0;
            std::stringstream ss;
            for (std::size_t i=0; i<p.size(); ++i)
                if (auto ch=Probability_value::test<output>(p[i],tolerance);!ch.has_value())
                {
                    if constexpr(output)
                        return Op_void(false,"p=\n"+ToString(p)+" value test on p[i] i="+std::to_string(i)+ch.error());
                    else
                        return Op_void(false,"");
                }
                else sum+=p[i];
            if  (std::abs(sum-1.0)<tolerance)
                return Op_void(true,"");
            else
            {
                if constexpr(output)
                    return Op_void(false," p=\n"+ToString(p)+" sum test sum="+ToString(sum)+"\n");
                return Op_void(false,"");
            }

        }

        template<class P>
        P static normalize(P&& p, double min_p)
        {
            double sum=0;
            for (std::size_t i=0; i<p.size(); ++i)
            {
                if (p[i]<min_p)
                {
                    p[i]=0;
                }
                else
                {
                    if (p[i]+min_p>1)
                    {
                        p[i]=1;
                    }
                    sum+=p[i];
                }
            }
            for (std::size_t i=0; i<p.size(); ++i)
                p[i]=p[i]/sum;
            return p;
        }
    };



    struct Probability_transition: public Probability
    {

        template<bool output>
        static bool test(const M_Matrix<double>& p,double tolerance)
        {
            for (std::size_t i=0; i<p.nrows(); ++i)
            {
                double sum=0;
                for (std::size_t j=0; j<p.ncols(); ++j)
                    if (!Probability_value::test<output>(p(i,j),tolerance))
                    {
                        if constexpr(output)
                            std::cerr<<" p=\n"<<p<<"\n tolerance="<<tolerance<<" i="<<i<<"j="<<j<<" p(i,j)="<<p(i,j)<<" \n";
                        return false;

                    }
                    else sum+=p(i,j);
                if (std::abs(sum-1.0)>tolerance)
                {
                    if constexpr(output)
                        std::cerr<<" p=\n"<<p<<"\n i="<<i<<" sum ="<<sum<<" 1-sum"<<1-sum<<" tolerance="<<tolerance<<"\n";
                    return false;
                }
            }
            return true;

        }

        M_Matrix<double> static normalize(M_Matrix<double>&& p, double min_p)
        {
            for (std::size_t i=0; i<p.nrows(); ++i)
            {
                if (p(i,i)>1-min_p)
                {
                    for (std::size_t j=0; j<p.ncols(); ++j)
                        p(i,j)=0;
                    p(i,i)=1;
                }
                else
                {
                    double sum=0;
                    for (std::size_t j=0; j<p.ncols(); ++j)
                    {
                        if (p(i,j)<min_p)
                            p(i,j)=0;
                        else
                            sum+=p(i,j);
                    }
                    if (sum!=1)
                        for (std::size_t j=0; j<p.ncols(); ++j) p(i,j)=p(i,j)/sum;
                }
            }
            return std::move(p);
        }



    };

    struct Probability_distribution_covariance: public Probability
    {

        static Op_void test_correlation(const M_Matrix<double>& p, std::size_t i, std::size_t j, double eps=std::sqrt(std::numeric_limits<double>::epsilon())*1000)
        {
            double corr=sqr(p(i,j))/p(i,i)/p(j,j);
            if ((corr+eps<0)||(corr-eps>1))
            {
                std::stringstream ss;
                ss<<"p: \n"<<p<<"\n p(ij)"<<p(i,j)<<"\n p(ii)"<<p(i,i)<<"\n p(jj)"<<p(j,j)<<"\n corr"<<corr<<" i="<<i<<" j="<<j;
                return Op_void(false,ss.str());
            }
            else
                return Op_void(true,"");

        }

        static M_Matrix<double> normalize(M_Matrix<double>&& p,double min_p)
        {
            std::set<std::size_t> non_zero_i;

            for (std::size_t i=0; i<p.nrows(); ++i) {
                if(p(i,i)<min_p)
                {
                    for (std::size_t j=0; j<p.ncols(); ++j)
                        p(i,j)=0;
                }
                else if(p(i,i)+min_p>1.0)
                {
                    for (std::size_t n=0; n<p.size(); ++n) p[n]=0;
                    p(i,i)=1;
                    return std::move(p);
                }
                else
                    non_zero_i.insert(i);
            }
            for (auto i:non_zero_i)
            {

                double sum=0;
                for (auto j: non_zero_i)
                    if (j!=i)
                        sum+=p(i,j);

                if (-sum!=p(i,i))
                {
                    auto sum_new=(sum-p(i,i))*0.5;
                    double f=sum_new/sum;
                    p(i,i)=-sum_new;
                    for (auto j: non_zero_i)
                        if (i!=j)
                            p(i,j)*=f;

                }
            }

            return std::move(p);
        }


        template<bool output,class P>
        static Op_void test(const P& p,double tolerance)
        {

            for (std::size_t i=0; i<p.nrows(); ++i)
            {
                if (!Probability_value::test<output>(p(i,i),tolerance))
                {
                    if constexpr(output)
                    {
                        std::stringstream ss;
                        ss<<"\ntolerance="<<tolerance<<"\n";
                        ss<<" pcov=\n"<<p<<"\n i="<<i<<" pcov(i,i)="<<p(i,i)<<"\n";
                        return Op_void(false,ss.str());
                    }
                    return Op_void(false,"");

                }double sum=0;
                for (std::size_t j=0; j<p.ncols(); ++j)
                {

                    if (i!=j)
                    {
                        if ((p(i,i)*p(j,j)-sqr(p(i,j))+tolerance<0)
                            &&(p(i,i)>tolerance*tolerance)&&(p(j,j)>tolerance*tolerance))
                        {
                            if constexpr(output){
                                double corr=sqr(p(i,j))/p(i,i)/p(j,j);
                                std::stringstream ss;
                                ss<<"tolerance="<<tolerance<<"\n";
                                ss<<" pcov=\n"<<p<<"\n i="<<i<<"j="<<j<<" pcov(i,j)="<<p(i,j)<<" corr="<<corr<<" pcov(i,i)"<<p(i,i)<<" pcov(j,j)"<<p(j,j)<<"\n";
                                return  Op_void(false, ss.str());
                            }
                            return Op_void(false,"");
                        }
                        else
                            sum+=p(i,j);
                    }
                }
                if (std::abs(p(i,i)+sum)>tolerance)
                {
                    if constexpr(output)
                    {
                        std::stringstream ss;
                        ss<<"tolerance="<<tolerance<<"\n";
                        ss<<" p=\n"<<p<<"\n i="<<i<<" p(i,j)="<<p(i,i)<<" sum="<<sum<<"\n";
                        return  Op_void(false, ss.str());

                    }
                    return Op_void(false,"");
                }
            }
            return Op_void(true,"");
        }

    };









    template<typename T>
    class Base_Probability_map
    {
    public:
        virtual T sample(std::mt19937_64& mt)const=0;

        virtual const std::map<T,double>& p() const=0;
        virtual std::ostream& put(std::ostream& os)const=0;

        virtual void reduce(double nmax)=0;
        virtual double nsamples()const =0;
        virtual ~Base_Probability_map(){}
    };

    template <class T>
    std::map<double,T>
    cumulative_reverse_map(const std::map<T,double>& normalized_map)
    {
        std::map<double,T> out;
        double sump=0;
        for (auto& e:normalized_map)
        {
            sump+=e.second;
            out[sump]=e.first;
        }
        return out;
    }

    template <class T>
    std::pair<std::map<T,double>,double>
    normalize_map(const std::map<T,double>& unnormalized_map)
    {
        std::map<T,double> out(unnormalized_map);
        double Evidence=0;
        for (auto& e:out)
        {
            Evidence+=e.second;
        }
        for (auto& e:out) e.second/=Evidence;

        return {out,Evidence};
    }


    template<typename T>
    class Probability_map: public Base_Probability_map<T>
    {
        static std::vector<T> to_vec(const std::map<T,double>& m)
        {
            std::vector<T> out(m.size());
            std::size_t i=0;
            for (auto& e: m) {out[i]=e.first; ++i;}
            return out;
        }

    public:
        typedef  Probability_map self_type;
        typedef  Cs<T> template_types;
        constexpr static auto const className=my_static_string("Probability_map")+my_trait<template_types>::className;
        ~Probability_map()override{}

        T sample(std::mt19937_64& mt)const override
        {
            return sample_rev_map(rev_,mt);
        }

        auto & x()const { return x_;}

        const std::map<T,double>& p() const override
        {
            return p_;
        }

        Probability_map(const std::map<T,double>& myNormalized_map, double nsamples)
            :
              p_{myNormalized_map},rev_{cumulative_reverse_map(p_)},x_{to_vec(myNormalized_map)}, nsamples_(nsamples)
        {

        }

        template<template<typename...>class V, typename...As>
        Probability_map(const V<T,As...>& x):p_(Uniform(x)),rev_(cumulative_reverse_map(p_)), nsamples_(0)
        {}
        Probability_map()=default;


        template<template<typename...>class V>
        static
            std::map<T,double> Uniform(const V<T>& x)
        {
            std::map<T,double> out;
            std::size_t n=x.size();
            double p=1.0/n;
            for (std::size_t i=0; i<n;++i )
                out[x[i]]+=p;
            return out;
        }

        template<class K>
        static
            std::map<T,double> Uniform(const std::map<T,K>& x)
        {
            std::map<T,double> out;
            std::size_t n=x.size();
            double p=1.0/n;
            for (auto& e:x )
                out[e.first]+=p;
            return out;
        }




        void reduce(double nmax) override
        {
            double f=nsamples_/nmax;
            if (f<1.0)
            {
                auto o=p_;
                for (auto& e:o) e.second=std::pow(e.second,f);
                *this=normalize(o,nmax).first;
            }

        }
        double nsamples()const override{return nsamples_;}

        static std::pair<Probability_map,double> normalize(const std::map<T,double>& myposterior , double nsamples)
        {
            auto out= normalize_map(myposterior);
            return {Probability_map(out.first,nsamples),out.second};
        }

        static auto get_constructor_fields()
        {
            return std::make_tuple(
                grammar::field(C<self_type>{},"probability_map",&self_type::p),
                grammar::field(C<self_type>{},"reverse_map",&self_type::reverse),
                grammar::field(C<self_type>{},"nsamples",&self_type::nsamples));
        }
        virtual std::ostream& put(std::ostream& os)const override
        {
            return io::output_operator_on_Object (os,*this);
        }

        std::map<double,T> const & reverse()const {return rev_;}

        Probability_map(const std::map<T,double>& myNormalized_map, const std::map<double,T>& reverseMap, double nsamples)
            :
              p_{myNormalized_map},rev_{reverseMap}, x_{to_vec(myNormalized_map)},nsamples_(nsamples)
        {

        }

        static Data_Index_scheme data_index(std::string&& index_name, const std::string& variable_name)
        {
            Data_Index_scheme out;
            out.push_back("prob_"+variable_name,{index_name});
            out.push_back("rev_prob_"+variable_name,{index_name});
            out.push_back("nsamples_"+variable_name,{});

            return out;

        }

        static auto
        get_data_index(const std::string& variable_name, const std::string& index_name)
        {
            using namespace std::literals::string_literals;


            return std::make_tuple(
                Data_Index(Cs<self_type>(),"prob_"+variable_name,
                           [](const self_type& s, std::size_t i_x){return s.p().at(s.x()[i_x]);},
                           std::make_pair(index_name,[](const self_type& s){return s.x().size();})
                                                                                             ),
                Data_Index(Cs<self_type>(),"nsamples"s+variable_name,&self_type::nsamples)
                                                            );
        }

        template<class i_index>
        static auto
        get_data_index_static()
        {
            return std::make_tuple(
                make_data_static
                (std::tuple<>(),
                 std::make_tuple(
                     F_s(s_num_samples{},&self_type::nsamples))),
                make_data_static
                (std::make_tuple(
                     I_s(i_index{},
                         [](const self_type& s){return s.x().size();},
                         [](const self_type& s, std::size_t i) {return s.x()[i];})),
                 std::make_tuple(
                     F_s(CT_s<s_use,s_probability>{},[](const self_type& s, std::size_t i_x){return s.p().at(s.x()[i_x]);}))));
        }


    private:
        std::map<T,double> p_;
        std::map<double,T> rev_;
        std::vector<T> x_;
        double nsamples_;
    };
    template <typename T, typename K>
    std::vector<T> to_vec(const std::map<T,K>& m)
    {
        std::vector<T> out(m.size());
        std::size_t i=0;
        for (auto& e: m) {out[i]=e.first; ++i;}
        return out;
    }

    template<typename T>
    class logLikelihood_map: public Base_Probability_map<T>
    {

    public:
        typedef  logLikelihood_map self_type;
        typedef  Cs<T> template_types;
        constexpr static auto const className=my_static_string("logLikelihood_map")+my_trait<template_types>::className;

        virtual std::ostream& put(std::ostream& os)const override
        {
            return io::output_operator_on_Object (os,*this);
        }

        T sample(std::mt19937_64& mt) const override
        {
            return sample_rev_map(rev_,mt);
        }

        const std::map<T,double>& logLik() const
        {
            return logLik_;
        }

        std::map<T,double>const & p()const override
        {
            return p_;
        }

        logLikelihood_map(const std::map<T,double>& mylogLikelihood_Map, double nsamples)
            :
              logLik_{mylogLikelihood_Map}
        {
            auto p=logLik_to_p(mylogLikelihood_Map);
            p_=std::move(p.first);
            Evidence_=p.second;
            rev_=cumulative_reverse_map(p_);
            x_=to_vec(p_);
            nsamples_=nsamples;
        }

        logLikelihood_map()=default;
        void reduce(double nmax) override
        {
            double f=nmax/nsamples_;
            if (f<1.0)
            {
                auto o=logLik_;
                for (auto& e:o) e.second*=f;
                *this=logLikelihood_map(o,nmax);
            }
        }


        double nsamples()const override {return nsamples_;}
        static auto get_constructor_fields()
        {
            return std::make_tuple(
                grammar::field(C<self_type>{},"log_probability_map",&self_type::logLik),
                grammar::field(C<self_type>{},"probability_map",&self_type::p),
                grammar::field(C<self_type>{},"reverse_map",&self_type::reverse),
                grammar::field(C<self_type>{},"evidence",&self_type::Evidence),
                grammar::field(C<self_type>{},"nsamples",&self_type::nsamples));
        }
        std::map<double,T> const & reverse()const {return rev_;}

        double Evidence()const {return Evidence_;}

        logLikelihood_map(const std::map<T,double>& mylogLikelihood_Map, const std::map<T,double>& myNormalized_map, const std::map<double,T>& reverseMap, double evidence, double nsamples)
            : logLik_{mylogLikelihood_Map},
              p_{myNormalized_map},rev_{reverseMap}, x_{to_vec(p_)},Evidence_{evidence},nsamples_(nsamples)
        {

        }

        static Data_Index_scheme data_index(std::string&& index_name, const std::string& variable_name)
        {
            Data_Index_scheme out;
            out.push_back("logLik_"+variable_name,{index_name});
            out.push_back("prob_"+variable_name,{index_name});
            out.push_back("rev_prob_"+variable_name,{index_name});
            out.push_back("Evidence_"+variable_name,{});
            out.push_back("nsamples_"+variable_name,{});

            return out;

        }
        static auto data_index(const std::string& variable_name)->decltype (T::data_index())
        {
            Data_Index_scheme out=myData_Index<T>::data_index();
            out.push_back("logLik_"+variable_name,out.first_index());
            out.push_back("prob_"+variable_name,out.first_index());
            out.push_back("rev_prob_"+variable_name,out.first_index());
            out.push_back("Evidence_"+variable_name,{});
            out.push_back("nsamples_"+variable_name,{});

            return out;

        }

        static auto
        get_data_index(const std::string& variable_name, const std::string& index_name)
        {
            using namespace std::literals::string_literals;
                return std::make_tuple(
                    Data_Index(Cs<self_type>(),"logLik_"+variable_name,
                               [](const self_type& s, std::size_t i_x){return s.logLik().at(s.x()[i_x]);},
                               std::make_pair(index_name,[](const self_type& s){return s.logLik().size();})
                                                                                                                                                                   ),
                    Data_Index(Cs<self_type>(),"prob_"+variable_name,
                               [](const self_type& s, std::size_t i_x){return s.p().at(s.x()[i_x]);},
                               std::make_pair(index_name,[](const self_type& s){return s.x().size();})
                                                                                                                                                                                                  ),
                    Data_Index(Cs<self_type>(),"Evidence"s+variable_name,&self_type::Evidence),
                    Data_Index(Cs<self_type>(),"nsamples"s+variable_name,&self_type::nsamples)

                                                                                                                                                        );

        }

        template<class i_index>
        static auto
        get_data_index_static(i_index)
        {

            return std::make_tuple(
                Compose_static(
                    std::make_tuple(I_s(i_index{},[](const self_type& s){return s.x().size();},0)),
                               [](const self_type& s, std::size_t i){return s.x()[i];},
                               myData_Index<T>::get_data_index_static()
                               ),
                make_data_static
                (std::tuple<>(),
                 std::make_tuple(
                     F_s(s_logEvidence{},&self_type::Evidence),
                     F_s(s_num_samples{},&self_type::nsamples))),
                make_data_static
                (std::make_tuple(
                     I_s(i_index{},
                         [](const self_type& s){return s.logLik().size();},0)
                         ),
                 std::make_tuple(
                     F_s(CT_s<s_obs,s_logL>{},[](const self_type& s, std::size_t i_x){return s.logLik().at(s.x()[i_x]);}),
                     F_s(CT_s<s_use,s_probability>{},[](const self_type& s, std::size_t i_x){return s.p().at(s.x()[i_x]);}))));
        }


        auto & x()const {return x_;}

    private:
        std::map<T,double> logLik_;
        std::map<T,double>p_;
        std::map<double,T> rev_;
        std::vector<T> x_;
        double Evidence_;
        double nsamples_;

    };


    template<class Likelihood,  typename T, class...Data>
    std::pair<Probability_map<T>,double>
    Bayes_rule(const Probability_map<T>& prior,const Likelihood& lik, const Data&... data)
    {
        auto p=prior.p();
        for (auto& e:p)
        {
            double l=lik(e.first,data...);
            e.second*=l;
        }
        double nsamples=prior.nsamples()+1;
        return Probability_map<T>::normalize(p,nsamples);
    }


    template<typename T>
    class Beta_map //: public Base_Probability_map<T>
    {

    public:
        typedef   Beta_map self_type;
        Beta_map(const std::map<T,Beta_Distribution>& a):a_(a), x_{to_vec(a)}{}


        auto &get_Map()const{return a_;}

        auto & x() const { return x_;}

        static auto get_constructor_fields()
        {
            return std::make_tuple(
                grammar::field(C<self_type>{},"probability_map",&self_type::get_Map));
        }

        void reduce(double nmax)
        {
            double f=nmax/count();
            if (f<1.0)
            {
                for (auto& e:a_)
                {
                    e=Beta_Distribution(e.alfa()*f,e.beta()*f);

                }
            }
        }



        static Beta_map UniformPrior(const std::map<T,double>& a)
        {
            std::map<T,Beta_Distribution> o;
            for (auto& e:a)
                o[e.first]=Beta_Distribution::UniformPrior();
            return Beta_map(o);
        }


        static Beta_map UnInformativePrior(const std::map<T,double>& a)
        {
            std::map<T,Beta_Distribution> o;
            for (auto& e:a)
                o[e.first]=Beta_Distribution::UnInformativePrior();
            return Beta_map(o);
        }

        static Beta_map UniformPrior(const std::vector<T>& a)
        {
            std::map<T,Beta_Distribution> o;
            for (auto& e:a)
                o[e]=Beta_Distribution::UniformPrior();
            return Beta_map(o);
        }


        static Beta_map UnInformativePrior(const std::vector<T>& a)
        {
            std::map<T,Beta_Distribution> o;
            for (auto& e:a)
                o[e]=Beta_Distribution::UnInformativePrior();
            return Beta_map(o);
        }




        Beta_map()=default;

        std::size_t size()const {return a_.size();}

        double count()const {
            double sum=0;
            for (auto& e:a_) sum+=e.second.count();
            return sum;
        }


        std::map<T,double> sample(std::mt19937_64& mt)
        {
            std::map<T,double> out;
            double sum=0;

            for (auto it=a_.begin(); it!=a_.end(); ++it)
            {
                std::gamma_distribution<double> g(it->second);
                out[it->first]=g(mt);
                sum+=out[it->first];
            }
            for (auto& o:out)
                o.second/=sum;
            return out;
        }

        Beta_map& operator+=(const Beta_map& other)
        {
            for (auto& e:a_)
            {
                auto it=other.a_.find(e.first);
                if (it!=other.a_.end())
                    e.second.Parameters()+=it->second.Parameters();
            }
            return *this;
        }

        Beta_map operator+(const Beta_map& other)const
        {
            Beta_map out(a_);
            out+=other;
            return out;
        }

        std::map<T,double> p()const
        {
            std::map<T,double> out;
            for (auto& e:a_)
                out[e.first]=e.second.p();
            return out;
        }

        Beta_Distribution& operator[](const T& x)
        {
            return a_[x];
        }

        Beta_Distribution operator[](const T& x)const
        {
            auto it=a_.find(x);
            if (it!=a_.end())
                return it.second;
            else
                return {};
        }

        Probability_map<T> Distribute_on_p(double desired_probability)const
        {
            auto prior= Probability_map<T>(a_);
            auto const & d=a_;
            return Bayes_rule(prior,[&d](const T& x,double target){return d.at(x).p(target);},desired_probability).first;
        }

        template<class F>
        Probability_map<T> Distribute_on_f(const F& f)const
        {
            auto prior= Probability_map<T>(a_);
            auto const & d=a_;
            return Bayes_rule(prior,[&d,&f](const T& x){return f(x)*d.at(x).p();}).first;
        }


        static Data_Index_scheme data_index(std::string&& index_name, const std::string& variable_name)
        {
            Data_Index_scheme out=Beta_Distribution::data_index(variable_name);
            out.insert_index(std::move(index_name));
            return out;

        }
        static auto
        get_data_index(const std::string& variable_name, const std::string& index_name)
        {
            using namespace std::literals::string_literals;
            return std::tuple_cat(

                Insert_tuple(Cs<self_type>(),
                             [](const self_type& s, std::size_t i_x){return s.get_Map().at(s.x()[i_x]);},Cs<std::size_t>(),
                             Beta_Distribution::get_data_index(variable_name),
                             std::make_pair(index_name,[](const self_type& s){return s.get_Map().size();})
                                                                  )

                                        );
        }
        template<class i_index>
        static auto
        get_data_index_static(i_index)
        {
            return Compose_static(
                std::make_tuple(I_s(i_index(),[](const self_type& s){return s.get_Map().size();},0)),
                [](const self_type& s, std::size_t i_x){return s.get_Map().at(s.x()[i_x]);},
                Beta_Distribution::get_data_index_static());
        }

    private:
        std::map<T,Beta_Distribution> a_;
        std::vector<T> x_;
    };


    template<class logLikelihood, class Data, typename T>
    logLikelihood_map<T>
    logBayes_rule(const logLikelihood& loglik, const Data& data, const logLikelihood_map<T>& prior)
    {
        auto logP=prior.logLik();
        for (auto& e:logP)
        {
            double logL=loglik(e.first,data);
            e.second+=logL;
        }
        double nsamples=prior.nsamples()+1;
        return logLikelihood_map<T>(logP,nsamples);
    }







    struct TargetProb
    {
        double operator()(const std::pair<std::size_t, std::size_t>& p)const
        {
            return BetaDistribution(p_target_,p.first,p.second);
        }
        TargetProb(double p_target):p_target_(p_target){}
        TargetProb(){}
    private:
        double p_target_;

    };

    struct PascalProb
    {
        double operator()(const std::pair<std::size_t, std::size_t>& p)const
        {
            return (1.0+p.first)/(2.0+p.first+p.second);
        }
    };



    template<class AP>
    struct One
    {
        double operator()(const AP& )const {return 1.0;}
    };



    template<typename T>
    class markov_process
    {

    public:
        markov_process(){};

        /*  * The type of the range of the distribution. */
        typedef T result_type;

        /** Parameter type. */
        struct param_type
        {
            typedef markov_process<T> distribution_type;
            friend class markov_process<T>;

            param_type()=default;
            explicit
                param_type(T _N, const M_Matrix<double>& P )
                : N_(_N), P_(P),rm_(P.nrows())

            {
                _M_initialize();
            }

            T
            N() const
            { return N_; }

            void set_P(M_Matrix<double>&&
                           _P)
            {
                P_=std::move(_P);
                _M_initialize();
            }

            void set_N(T _N)
            {  N_=_N; }

            const M_Matrix<double>&
            P() const
            { return P_; }


            friend bool
            operator==(const param_type& __p1, const param_type& __p2)
            { return __p1.N_ == __p2.N_ && __p1.P_ == __p2.P_; }

        private:
            void
            _M_initialize()
            {

                std::size_t n=P_.nrows();
                for (std::size_t i=0; i<n; ++i)
                {
                    double sum=0;
                    std::vector<std::pair<double, std::size_t>> o(n);
                    for (std::size_t j=0; j<n; ++j)
                    {o[j].second=j; o[j].first=P_(i,j);}
                    std::sort(o.begin(),o.end());
                    //  std::cerr<<" \n--- o---\n"<<o;
                    for (std::size_t j=0; j<n; ++j)
                    {
                        sum+=o[j].first;
                        rm_[i][sum]=o[j].second;
                    }
                }
                //     std::cerr<<"\n mapa inverso -----------------\n"<<rm_;
            }

            T N_;
            M_Matrix<double> P_;
            std::vector<std::map<double,T>> rm_;
        };



        // constructors and member function

        markov_process(T __t,
                       const M_Matrix<double>& __p)
            : _M_param(__t, __p)
        { }



        /**
       * @brief Returns the distribution @p t parameter.
       */
        T
        N() const
        { return _M_param.N(); }




        /**
       * @brief Returns the distribution @p p parameter.
       */
        const M_Matrix<double>&
        P() const
        { return _M_param.P(); }


        void set_P(M_Matrix<double>&& _P)
        { _M_param.set_P(std::move(_P));
        }

        void set_N(T _N)
        {
            _M_param.set_N(_N);
        }


        /**
       * @brief Returns the parameter set of the distribution.
       */
        param_type
        param() const
        { return _M_param; }

        /**
       * @brief Sets the parameter set of the distribution.
       * @param __param The new parameter set of the distribution.
       */
        void
        param(const param_type& __param)
        { _M_param = __param; }

        /**
       * @brief Returns the greatest lower bound value of the distribution.
       */
        T
        min() const
        { return 0; }

        /**
       * @brief Returns the least upper bound value of the distribution.
       */
        T
        max() const
        { return _M_param.P().nrows(); }

        /**
       * @brief Generating functions.
       */
        template<typename _UniformRandomNumberGenerator>
        result_type
        operator()(_UniformRandomNumberGenerator& __urng)
        { return this->operator()(__urng, _M_param); }

        template<typename _UniformRandomNumberGenerator>
        result_type
        operator()(_UniformRandomNumberGenerator& __urng,
                   const param_type& _p)const
        {
            return sample_rev_map(_p.rm_[_p.N()],__urng);
        }


    private:

        param_type _M_param;

    };




#endif // MYDISTRIBUTIONS_H
