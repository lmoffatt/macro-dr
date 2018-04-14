#ifndef MYCOMMANDMANAGEMENT_H
#define MYCOMMANDMANAGEMENT_H

#include <string>
#include <map>

#include <vector>
#include <initializer_list>
#include <sstream>
#include <fstream>

#include <chrono>

#include <type_traits>
#include "myTuples.h"
#include "myInputSerializer.h"






template<template<typename> class Cls,template<typename, typename> class Co, typename T, typename Alloc,typename std::enable_if<!std::is_pointer<T>::value, int>::type=0>
bool load_from_file(C<Co<T,Alloc>>, const std::string& file_name, Co<T,Alloc>& v, std::ostream& logstream)
{
    std::vector<std::string> names;
    std::stringstream ss(file_name);
    if (ss>>names)
    {
        for (std::string filename:names)
        {

            T x;
            if (load_from_file<Cls>(C<T>(),filename,x,logstream))
                v.push_back(x);
            else
                return false;
        }
        return true;
    }
    else
        return false;
}


template<template<typename> class Cls,template<typename, typename> class Co, typename T, typename Alloc>
bool load_from_file(C<Co<T*,Alloc>>,const std::string& file_name, Co<T*,Alloc>& v, std::ostream& logstream)
{
    std::vector<std::string> names;
    std::stringstream ss(file_name);
    if (ss>>names)
    {
        for (std::string filename:names)
        {

            T* x=new T;
            if (load_from_file<Cls>(C<T*>(),filename,x,logstream))
                v.push_back(x);
            else
            {
                delete x;
                return false;
            }
        }
        return true;
    }
    else
        return false;
}


template<template<typename> class Cls,typename T,typename std::enable_if<!is_container<T>::value, int>::type=0>
bool load_from_file(C<T>,const std::string& file_name, T& x, std::ostream& logstream)
{
    std::string filename=file_name;
    std::ifstream fi;
    fi.open(filename.c_str(),std::ios_base::in);
    if (!fi)
    {
        fi.close();
        filename=file_name+".txt";
        fi.open(filename.c_str(),std::ios_base::in);
        if (!fi)
        {
            fi.close();
            filename=Cls<T>::name()+"_"+file_name+".txt";
            fi.open(filename.c_str(),std::ios_base::in);
            if (!fi)
            {
                fi.close();
                logstream<<"cannot load "<<file_name<<" file "<<filename<<" not found\n";
                return false;
            }
            else
                logstream<<"file "<<filename<<" opened successfully\n";
        }
        else
            logstream<<"file "<<filename<<" opened successfully\n";

    }
    else
        logstream<<"file "<<filename<<" opened successfully\n";
    if (Cls<T>::read(fi,x,logstream))
    {
        logstream<<Cls<T>::name()<<" "<<file_name<<" loaded successfully \n";
        fi.close();
        return true;
    }
    else
    {
        logstream<<Cls<T>::name()<<" "<<file_name<<" not loaded successfully \n";
        fi.close();
        return false;
    }

}



template<typename Cm>
class myBaseCommand
{
public:
    virtual bool run(Cm* cm,
                     const std::string& idResult,
                     const std::map<std::string,std::string>& args,
                     std::ostream& logstream)=0;
};


namespace impl_run {


template <class Cm,template<typename>class Cls,class F, typename R,typename... Argsout>
std::pair<R,bool> run_R_impl_1(Cm* /*cm*/,
                               Co<Cls>,
                               Cs<R>,
                               const F& f,
                               std::ostream& /*log*/,
                               std::map<std::string, std::string> /*in*/,
                               Argsout...argsOut)
{
    return {f(argsOut...),true};

}


template <class Cm,template<typename>class Cls,class F,typename R,typename T,typename... ArgsIn,typename... Argsout>
std::pair<R,bool> run_R_impl_1(Cm* cm,
                               Co<Cls>,
                               Cs<R>,
                               const F& f,
                               std::ostream& log,
                               const std::map<std::string, std::string>& in,
                               std::pair<T, std::string> a,
                               std::pair<ArgsIn,std::string>... argsIn,
                               Argsout...argsOut)
{
    auto it=in.find(a.second);
    if (it==in.end())
    {
        log<<a.second<<"=  [not listed, use default], ";
    }
    else
    {
        std::stringstream ss(it->second);
        if (Cls<T>::read(ss,a.first,log))
        {
            log<<it->first<<"="+it->second<<", ";
        }
        else if (cm->get(C<T>(),it->second,a.first))
        {
            log<<it->first<<"="<<it->second<<" [variable], ";
        }
        else  if (load_from_file<Cls>(C<T>(),it->second,a.first,log))
        {
            log<<it->first<<"="<<it->second<<" [in file], ";
        }
        else
        {
            return {{},false};


        }
    }
    return run_R_impl_1(cm,Co<Cls>(),Cs<R>(),f,log,in,argsIn...,argsOut...,a.first );
}


template <class Cm,template<typename>class Cls,class F,class R,typename T,typename... ArgsIn>
std::pair<R,bool> run_R_impl_1(Cm* cm,
                               Co<Cls>,
                               Cs<R>,
                               const F& f,
                               std::ostream& log,
                               const std::map<std::string, std::string>& in,
                               std::pair<T, std::string> a,
                               std::pair<ArgsIn,std::string>... argsIn)
{
    auto it=in.find(a.second);
    if (it==in.end())
    {
        log<<a.second<<"=  [not listed, use default], ";
    }
    else
    {
        if (cm->template get(C<T>(),it->second,a.first))
        {
            log<<it->first<<"="<<it->second<<" [variable], ";
        }
        else  if (load_from_file<Cls,T>(it->second,a.first,log))
        {
            log<<it->first<<"="<<it->second<<" [in file], ";
        }
        else
        {
            std::stringstream ss(it->second);
            if (ss>>a.second)
            {
                log<<it->first<<"="+it->second<<", ";
            }
            else
            {
                return {{},false};

            }
        }
    }
    return run_R_impl_1(cm,Co<Cls>(),Cs<R>(),f,log,in,argsIn...,a.first );
}



template <class Cm,template <typename T>class Cls,class F, class R,std::size_t... Is,typename... ArgsIn>
std::pair<R,bool> run_R_impl_0
(Cm* cm,
 Co<Cls>,
 Cs<R>,
 const F& f,
 std::ostream& ls,
 const std::map<std::string, std::string>& in,
 std::index_sequence<Is...>,
 std::tuple<std::pair<ArgsIn,std::string>...> argsIn
 )
{
    return run_R_impl_1(cm,Co<Cls>(),Cs<R>(),f,ls,in,std::get<Is>(argsIn)...);
}



} // namespace run_impl




/*!
 * \brief wrapper class around a function
 *  \tparam Cm Command Manager class
 * \tparam F function class
 * \tparam R return type
 * \tparam Args ordered list of function arguments
 */
template <class Cm,template<typename>class Cls,class F,class R,typename... Args>
class myCommand: public myBaseCommand<Cm>
{
public:

    /*!
   * \brief run run command on line input
   * \param cm_ the commandManager
   * \param line arguments: arg1=value1, arg2=value2 value could be text value or id
   * \param logstream logs the arguments actually used with their origin
   *
   * \return object of running the command and log text
   *
   * the arguments could be in any order. If an argument is absent, its default value is used instead. The log strean receives all the arguments actually used with their origin
   */
    std::pair<R,bool> runit(Cm* cm_,const std::map<std::string, std::string>& m, std::ostream& logstream)
    {
        return  impl_run::run_R_impl_0(cm_,
                             Co<Cls>(),
                             Cs<R>(),
                             f_,
                             logstream,
                             m,
                             std::index_sequence_for<Args...>(),
                             args_);
    }

    myCommand(const F& f,std::pair<Args,std::string>... args):
        f_(f),
        args_{args...}{}

private:
    const F& f_;
    std::tuple<std::pair<Args,std::string>...> args_;


    // myBaseCommand interface
public:
    virtual bool  run(Cm* cm,
                      const std::string& idResult,
                      const std::map<std::string,std::string>& args,
                      std::ostream& logstream) override
    {

        std::pair<R,bool> o=runit(cm,args,logstream);
        if (o.second)
            cm->template push_back<R>(idResult,o.first);
        return o.second;

    }
};




namespace impl_run_void {


template <class Cm,template<typename>class Cls,class F, typename... Argsout>
bool  run_void_impl_1(Cm* /*cm*/,
                      Co<Cls>,
                      const F& f,
                      std::ostream& /*log*/,
                      std::map<std::string, std::string> /*in*/,
                      Argsout...argsOut)
{
    std::cerr<<"f(argsOut..)\n";
    f(argsOut...);
    return true;

}


template <class Cm,template<typename>class Cls,class F,typename T,typename... ArgsIn,typename... Argsout>
bool run_void_impl_1(Cm* cm,
                     Co<Cls>,
                     const F& f,
                     std::ostream& log,
                     const std::map<std::string, std::string>& in,
                     std::pair<T, std::string> a,
                     std::pair<ArgsIn,std::string>... argsIn,
                     Argsout...argsOut)
{
    auto it=in.find(a.second);
    if (it==in.end())
    {
        log<<a.second<<"=  [not listed, use default], ";
    }
    else
    {
        std::stringstream ss(it->second);
        if (Cls<T>::read(ss,a.first,log))
        {
            log<<it->first<<"="+it->second<<", ";
        }
        else if (cm->template get(C<T>(),it->second,a.first))
        {
            log<<it->first<<"="<<it->second<<" [variable], ";
        }
        else  if (load_from_file<Cls>(C<T>(),it->second,a.first,log))
        {
            log<<it->first<<"="<<it->second<<" [in file], ";
        }
        else
        {
            return false;
        }
    }
    return run_void_impl_1(cm,Co<Cls>(),f,log,in,argsIn...,argsOut...,a.first );
}


template <class Cm,template<typename>class Cls,class F,typename T,typename... ArgsIn>
bool run_void_impl_1(Cm* cm,
                     Co<Cls>,
                     const F& f,
                     std::ostream& log,
                     const std::map<std::string, std::string>& in,
                     std::pair<T, std::string> a,
                     std::pair<ArgsIn,std::string>... argsIn)
{
    auto it=in.find(a.second);
    if (it==in.end())
    {
        log<<a.second<<"=  [not listed, use default], ";
    }
    else
    {
        if (cm->template get(C<T>(),it->second,a.first))
        {
            log<<it->first<<"="<<it->second<<" [variable], ";
        }
        else  if (load_from_file<Cls>(C<T>(),it->second,a.first,log))
        {
            log<<it->first<<"="<<it->second<<" [in file], ";
        }
        else
        {
            std::stringstream ss(it->second);
            if (ss>>a.second)
            {
                log<<it->first<<"="+it->second<<", ";
            }
            else
            {
                return false;

            }
        }
    }
    return run_void_impl_1(cm,Co<Cls>(),f,log,in,argsIn...,a.first );
}



template <class Cm,template <typename T>class Cls,class F, std::size_t... Is,typename... ArgsIn>
bool run_void_impl_0
(Cm* cm,
 Co<Cls>,
 const F& f,
 std::ostream& ls,
 const std::map<std::string, std::string>& in,
 std::index_sequence<Is...>,
 std::tuple<std::pair<ArgsIn,std::string>...> argsIn
 )
{
    return run_void_impl_1(cm,Co<Cls>(),f,ls,in,std::get<Is>(argsIn)...);
}




} // namespace run_void_impl



/*!
 * \brief wrapper class around a function
 *  \tparam Cm Command Manager class
 * \tparam F function class
 * \tparam R return type
 * \tparam Args ordered list of function arguments
 */
template <class Cm,template<typename>class Cls,class F,typename... Args>
class myCommand<Cm,Cls,F,void,Args...>: public myBaseCommand<Cm>
{
public:

    /*!
   * \brief run run command on line input
   * \param cm_ the commandManager
   * \param line arguments: arg1=value1, arg2=value2 value could be text value or id
   * \param logstream logs the arguments actually used with their origin
   *
   * \return object of running the command and log text
   *
   * the arguments could be in any order. If an argument is absent, its default value is used instead. The log strean receives all the arguments actually used with their origin
   */
    bool runit(Cm* cm_,const std::map<std::string, std::string>& m, std::ostream& logstream)
    {
        return impl_run_void::run_void_impl_0(cm_,
                               Co<Cls>(),
                               f_,
                               logstream,
                               m,
                               std::index_sequence_for<Args...>(),
                               args_);
    }

    myCommand(const F& f,std::pair<Args,std::string>... args):
        f_(f),
        args_{args...}{}

private:
    const F& f_;
    std::tuple<std::pair<Args,std::string>...> args_;


    // myBaseCommand interface
public:
    virtual bool  run(Cm* cm,
                      const std::string& ,
                      const std::map<std::string,std::string>& args,
                      std::ostream& logstream) override
    {

        return runit(cm,args,logstream);
    }
};




template <class Cm,template <typename>class Cls,class R,class F, typename... Args>
myCommand<Cm,Cls,F,R,Args...>* make_Command
(C<Cm>,C<R>,Co<Cls>,const F& f,std::pair<Args,std::string>... args)
{
    return new myCommand<Cm,Cls,F,R,Args...>(f,args...);
}


namespace impl_get {


template<typename T, typename...Ts>
bool get_impl_1(const std::string& /*id*/, T& /*x*/)
{
    return false;
}


template<typename T, typename K, typename...Ts>
bool get_impl_1(const std::string& id, T& x,
                std::map<std::string,K>,
                std::map<std::string,Ts>... m)
{
    return get_impl_1(id,x,m...);
}

template<typename T, typename...Ts>
bool get_impl_1(const std::string& id, T& x,
                std::map<std::string,T> mymap,
                std::map<std::string,Ts>... )
{
    auto it=mymap.find(id);
    if (it!=mymap.end())
    {
        x=it->second;
        return true;
    }
    else return false;
}

template<typename T, std::size_t...Is, typename...Ts>
bool get_impl_0(const std::string& id, T& x,
                std::index_sequence<Is...>,
                std::tuple<std::map<std::string,Ts>...> m)
{
    return get_impl_1(id,x,std::get<Is>(m)...);
}

} // namespace impl_get

template<typename T, typename...Ts>
bool get_map(const std::string& id, T& x, std::tuple<std::map<std::string,Ts>...> m)
{
    return impl_get::get_impl_0(id,x,std::index_sequence_for<Ts...>(),m);
}

template<template<typename> class Cls,typename...Ts>
class myCommandManager;


namespace impl_apply {

template<template<typename...> class F, class Cm,typename... Args>
bool apply_impl(const Co<F>&,Cm* ,const std::string& , Cs<>, Args... )
{
    return false;
}



template<template<typename...> class F, class Cm,typename... Args, typename T,typename... Ts>
bool apply_impl(const Co<F>& f,Cm* cm,const std::string& id, Cs<T,Ts...>, Args... args)
{
    T x;
    if (cm->get(C<T>(),id,x))
    {
        F<T>::apply(x,args...);
        return true;
    }
    else return apply_impl(f,cm,id,Cs<Ts...>(),args...);
}


} // namespace impl_apply



/*!
 * \tparam Cls
 * Cls<T>::name()
 * Cls<T>::read(std::istream&,T& x,std::ostream& logstream)
 */
template<template<typename> class Cls,typename...Ts, class... Tptrs>
class myCommandManager<Cls,Cs<Ts...>,Cs<Tptrs...>>
{
public:
    template<typename T>
    using myCls=Cls<T>;

    static std::string ClassName(){return "CommandManager";}

    typedef myBaseCommand<myCommandManager> Command;
    std::tuple<std::string,std::string, std::map<std::string,std::string>>
    line_to_CommandName_Arg_Result(const std::string& line)
    {
        auto n0=line.find_first_of('(');
        auto n1=line.find_first_of(')');
        auto n2=line.find_first_of('>',n1);
        std::string commandName=line.substr(0,n0);
        std::string argsList=line.substr(n0+1,n1-n0-1);
        std::string resultName;
        if (n2!=line.npos)
            resultName=line.substr(n2+1);
        std::map<std::string, std::string> out;
        // name0=value0, name1=value1,

        std::stringstream ss(argsList);
        std::string s;
        while (std::getline(ss,s,','))
        {
            auto n0=s.find_first_not_of(" ");
            auto npos=s.find('=');
            out[s.substr(n0,npos-n0)]=s.substr(npos+1);
        }


        return {commandName,resultName,out};
    }


    void execute(std::string line, std::ostream& logstream)
    {
        if (!line.empty())
        {
            std::stringstream ss(line);
            auto o=line_to_CommandName_Arg_Result(line);
            if (!std::get<0>(o).empty())
            {
                auto cv=command(std::get<0>(o));
                if (!cv.empty())
                {
                    std::size_t i;
                    for (i=0; i<cv.size();++i)
                    {
                        if(cv[i]->run(this,std::get<1>(o),std::get<2>(o),logstream))
                        {
                            logstream<<line<<"\n";
                            break;
                        }
                    }
                    if (i==cv.size())
                        logstream<<"error in"<<line<<" "<<std::get<0>(o)<<" substitution error";
                }
                else
                    logstream<<"error in"<<line<<" "<<std::get<0>(o)<<" is not a command";
            }
        }
    }


    //CommandBase* Command(std::string commandName);

    template<typename T, typename std::enable_if<!is_container<T>::value,int>::type=0>
    bool get(C<T>,const std::string& id, T& x)
    {
        return get_map(id,x,data_);
    }
    template<typename T, typename std::enable_if<!is_container<T>::value,int>::type=0>
    bool get(C<T*>,const std::string& id, T*& x)
    {
        return get_map(id,x,dataPtr_);
    }

    template< typename T, typename Alloc,
              typename std::enable_if<!std::is_pointer<T>::value,int>::type = 0>
    bool get(C<std::vector<T, Alloc>>,const std::string& id, std::vector<T,Alloc>& v)
    {
        if (get_map(id,v,data_))
            return true;
        else
        {
            std::vector<std::string> ids;
            std::stringstream ss(id);
            if (ss>>ids)
            {
                for (std::string name:ids)
                {
                    T x;
                    if (!get_map(name,x,data_))
                        return false;
                    else
                        v.push_back(x);
                }
                return true;
            }
            else return false;
        }
    }

    template< typename T, typename Alloc>
    bool get(C<std::vector<T*, Alloc>>, const std::string& id, std::vector<T*,Alloc>& v)
    {
        if (get_map(id,v,data_))
            return true;
        else
        {
            std::vector<std::string> ids;
            std::stringstream ss(id);
            if (ss>>ids)
            {
                for (std::string name:ids)
                {
                    T *x=new T;
                    if (!get_map(name,x,dataPtr_))
                    {
                        delete x;
                        return false;
                    }
                    else
                        v.push_back(x);
                }
                return true;
            }
            else return false;
        }
    }



    template<template<typename...>class F, typename... Args>
    bool apply(const Co<F>& f,const std::string& id, Args... args)
    {
        if (impl_apply::apply_impl(f,this,id,Cs<Tptrs*...>(),std::forward<Args>(args)...))
            return true;
        else return impl_apply::apply_impl(f,this,id,Cs<Ts...>(),args...);
    }



    template <typename T>
    void push_back(const std::string& id, T* x)
    {
        auto& m=std::get<typename std::map<std::string,T*>>(dataPtr_);
        m[id]=x;
    }

    template <typename T>
    void push_back(const std::string& id, T x)
    {
        auto m=std::get<typename std::map<std::string,T>>(data_);
        m[id]=x;
    }


    void push_command(const std::string& id, Command * cmd)
    {
        cmds_.insert({id,cmd});
    }

    std::vector<Command*> command(const std::string& id)
    {
        auto p=cmds_.equal_range(id);
        std::vector<Command*> o;
        for (auto it=p.first; it!=p.second; ++it)
            o.push_back(it->second);
        return o;
    }

    myCommandManager():
        data_{},
        dataPtr_{},
        cmds_{}{}



private:
    std::tuple<std::map <std::string, Ts>...> data_;
    std::tuple<std::map <std::string, Tptrs*>...> dataPtr_;

    std::multimap<std::string,myBaseCommand<myCommandManager>*> cmds_;


};




#endif // MYCOMMANDMANAGEMENT_H

