#ifndef MYOPTIMIZATION_H
#define MYOPTIMIZATION_H

#include "Matrix.h"
#include "mymath.h"
#include "mySerializer.h"

#include <string>
#include <random>
#include <chrono>
#include <iomanip>
namespace opt {

// ------------- UTILITY---------------
template<int...> struct index_tuple{};

template<int I, typename IndexTuple, typename... Types>
struct make_indexes_impl;

template<int I, int... Indexes, typename T, typename ... Types>
struct make_indexes_impl<I, index_tuple<Indexes...>, T, Types...>
{
    typedef typename make_indexes_impl<I + 1, index_tuple<Indexes..., I>, Types...>::type type;
};

template<int I, int... Indexes>
struct make_indexes_impl<I, index_tuple<Indexes...> >
{
    typedef index_tuple<Indexes...> type;
};

template<typename ... Types>
struct make_indexes : make_indexes_impl<0, index_tuple<>, Types...>
{};




inline
std::string leadingZero(int i)
{
  if (i==0)
    return "00";
  else if (i<10)
    return "0"+std::to_string(i);
  else return std::to_string(i);
}

inline
std::string leadingZeroZero(int i)
{
  if (i==0)
    return "000";
  else if (i<10)
    return "00"+std::to_string(i);
  else if (i<100)
    return "0"+std::to_string(i);
  else return std::to_string(i);
}

inline
int rename_done(const std::string& f_name)
{
  std::string f_name_done=f_name+".done";
  int res= std::rename(f_name.c_str(),f_name_done.c_str());
  if (res!=0)
    std::cerr<<"could not rename file"<<f_name<<"\n";
  return res;
}


inline
std::string time_now()
{
  auto tc=std::chrono::system_clock::now();
  std::time_t rawtime=std::chrono::system_clock::to_time_t(tc);

  auto tcount=(tc.time_since_epoch().count()/1000)%1000000;


  struct std::tm * t;
  time (&rawtime);
  t = localtime (&rawtime);
  return leadingZero(t->tm_hour)+leadingZero(t->tm_min)+leadingZero(t->tm_sec)+"s"+std::to_string(tcount);

}

inline
std::pair<std::size_t, std::string>
extract_Seed(const std::string& s)
{
   auto last=s.find("_state");
   auto first=s.find_last_of('_',last-1);
   auto v= s.substr(first+1,last-first-1);
   auto val=std::stoull(v);
   auto eviName=s.substr(0,last);
   return {val,eviName};
}


template<typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<T>& v)
{
  s<<"[";
  for (T x:v)
    s<<x<<"\t";
  s<<"]";
  return s;
}

template<typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<std::vector<T>>& m)
{
  s<<"[";
  for (const std::vector<T>& v:m)
    s<<v;
  s<<"]";
  return s;
}




template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os,const std::pair<T1,T2>& other)
{
  os<<other.first<<","<<other.second;
  return os;
}



template<typename K,typename T>
std::ostream& operator<<(std::ostream& s, const std::map<K,T>& v)
{
  s<<"{";
  for (auto it=v.begin(); it!=v.end(); ++it)
    s<<*it<<"\t";
  s<<"}";
  return s;
}

template<typename K,typename T>
std::ostream& operator<<(std::ostream& s, const std::multimap<K,T>& v)
{
  for (auto& it=v.begin(); it!=v.end(); ++it)
    s<<*it<<"\t";
  s<<"\n";
  return s;
}

template<typename T>
std::ostream& operator<<(std::ostream& s, const std::set<T>& v)
{
  for (auto& it=v.begin(); it!=v.end(); ++it)
    s<<*it<<"\t";
  s<<"\n";
  return s;
}

template<typename T>
std::ostream& operator<<(std::ostream& s, const std::multiset<T>& v)
{
 s<<"{";
  for (auto it=v.begin(); it!=v.end(); ++it)
    s<<*it<<"\t";
  s<<"}";
  return s;
}


inline
std::istream &safeGetline(std::istream &is, std::string &t)
{
  is.clear();
  std::getline(is,t);
  auto it=t.find('\r');
  if (it!=t.npos)
    t.erase(it);
  return is;
}


template<typename T>
struct value_wrapper
{
  value_wrapper(T& x):value(x){}


  T& value;
};

inline
std::istream& extract_infinite(std::istream& ss, value_wrapper<double>& val, bool is_negative)
{
  std::string s;
  ss>>s;
  if ((s=="inf")||(s=="infinity")||(s=="INF")||(s=="INFINITY"))
    {
      if (is_negative)
        val.value=-std::numeric_limits<double>::infinity();
      else
        val.value=std::numeric_limits<double>::infinity();
      return ss;
    }
  else
    {
      ss.setstate(std::ios::failbit);
      return ss;
    }
}

inline
std::istream& extract_nan(std::istream& ss, value_wrapper<double>& val, bool is_negative)
{
  std::string s;
  ss>>s;
  if ((s=="nan")||(s=="NAN"))
    {
      if (is_negative)
        val.value=-std::numeric_limits<double>::quiet_NaN();
      else
        val.value=std::numeric_limits<double>::quiet_NaN();
      return ss;
    }
  else
    {
      ss.setstate(std::ios::failbit);
      return ss;
    }
}

inline
std::istream& extract_finite(std::istream& ss, value_wrapper<double>& val, bool is_negative)
{
  ss>>val.value;
  if (is_negative)
    val.value=-val.value;
  return ss;

}


inline
std::istream& operator>>(std::istream& ss, value_wrapper<double>& val)
{
  char c;
  bool is_negative=false;
  ss>>c;
  while (std::isspace(c))
    ss>>c;
  if (c=='-')
    {
      is_negative=true;
      ss>>c;
    }
  else if (c=='+')
    ss>>c;
  switch(c)
    {
    case '0':
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
      {
        ss.putback(c);
        return extract_finite(ss,val,is_negative);
      }

    case 'i':
    case 'I':
      {
        ss.putback(c);
        return extract_infinite(ss,val,is_negative);
      }
    case 'n':
    case 'N':
      {
        ss.putback(c);
        return extract_nan(ss,val,is_negative);
    default:
          {
            ss.setstate(std::ios::failbit);
            return ss;
          }
      }

    }

}

template <typename T>
std::istream& operator>>(std::istream& ss, value_wrapper<T>& val)
{
  ss>>val.value;
  return ss;
}

template<typename T>
auto operator>>(std::istream& is, T& v)
->decltype(v.read(std::declval<std::string&>(),is,std::declval<std::ostream&>()),is)
{
  std::string s;
  if (!v.read(s,is, std::cerr))
    {
      is.setstate(std::ios::failbit);

    }
  return is;

}




template<typename T>
auto operator>>(std::istream& is, T*& v)
->decltype(v->read(std::declval<std::string&>(),is,std::declval<std::ostream&>() ),is)
{
  std::string s;
  if (v!=nullptr)
    {
      v->read(s,is,std::cerr);
    }
  else
    {
      is.setstate(std::ios::failbit);
    }
  return is;

}




template<typename T
         ,typename std::enable_if<!std::is_pointer<T>::value,int>::type = 0>
std::istream& operator>>(std::istream& is, std::vector<T>& v)
{
  std::string line;
  char c;
  is>>c;
  if (c=='[')
    {
      v.clear();
      while (c!=']')
        {
          is>>c;
          if (c!=']')
            {
              is.putback(c);
              T e;
              value_wrapper<T> w(e);
              if(is>>w)
                v.push_back(e);
              else
                break;
            }
        }
      return is;
    }
  else
    {
      safeGetline(is,line);
      while (v.empty()  &&line.empty()&& is.good())
        safeGetline(is,line);

      T x;
      value_wrapper<T> w(x);
      std::stringstream ss(line);
      while (is>>w)
        v.push_back(x);
      return is;
    }
}
template<typename T>

std::istream& operator>>(std::istream& is, std::vector<T*>& v)
{
  std::string line;
  char c;
  is>>c;
  if (c=='[')
    {
      while (c!=']')
        {
          is>>c;
          if (c!=']')
            {
              T* e=new T{};
              value_wrapper<T> w(*e);
              is.putback(c);
              if(is>>w)
                v.push_back(e);
              else
                break;
            }
        }
      return is;
    }
  else
    {
      safeGetline(is,line);
      while (v.empty()  &&line.empty()&& is.good())
        safeGetline(is,line);

      while (true)
        {
          T* x=new T;
          value_wrapper<T> w(*x);
          std::stringstream ss(line);
          if (is>>w)
            {
              v.push_back(x);
            }
          else
            {
              delete x;
              break;
            }
        }
      return is;
    }
}


template<typename T>
std::istream& operator>>(std::istream& is, std::vector<std::vector<T>>& m)
{
  char c;
  is>>c;
  if (c=='[')
    {
      while (c!=']')
        {
          is>>c;
          if (c!=']')
            {
              is.putback(c);
              std::vector<T> e;
              if(is>>e)
                m.push_back(e);
              else
                break;
            }
        }
      return is;
    }
  else
    {
      is.putback(c);
      std::vector<T> v;
      while((is>>v)&& !v.empty())
        {
          m.push_back(v);
          v.clear();
        }
      return is;
    }
}

template<typename T1, typename T2>
std::istream& operator>>(std::istream& os,std::pair<T1,T2>& other)
{
  char ch;
  value_wrapper<T1> w(other.first);
  value_wrapper<T2> w2(other.second);

  os>>w>>ch>>w2;
  return os;
}



template<typename K,typename T>
std::istream& operator>>(std::istream& is, std::map<K,T>& v)
{ char c;
  is>>c;
  if (c=='{')
    {
      while (c!='}')
        {
          is>>c;
          while (std::isspace(c))
            is>>c;
          if (c!='}')
            {
              is.putback(c);
              std::pair<K,T> e;
              if(is>>e)
                v.insert(e);
              else
                break;
            }
        }
      return is;
    }
  else
    {
      is.setstate(std::ios::failbit);
      return is;
    }
}


template<typename K,typename T>
std::istream& operator>>(std::istream& is,  std::multimap<K,T>& v)
{
  char c;
  is>>c;
  if (c=='{')
    {
      while (c!='}')
        {
          is>>c;
          if (c!='}')
            {
              std::pair<K,T> e;
              is.putback(c);
              if(is>>e)
                v.insert(e);
              else
                break;
            }
        }
      return is;
    }
  else
    {
      is.setstate(std::ios::failbit);
    }
}
template<typename T>
std::istream& operator>>(std::istream& is, std::multiset<T>& v)
{
  char c;
  is>>c;
  if (c=='{')
    {
      while (c!='}')
        {
          is>>c;
          if (c!='}')
            {
              T e;
              is.putback(c);
              if(is>>value_wrapper<T>(e))
                v.insert(e);
              else
                break;
            }
        }
      return is;
    }
  else
    {
      is.setstate(std::ios::failbit);
      return is;
    }
}


template<typename Last>
std::istream& get_impl(std::istream& is, Last& last)
{
  is>>last;
  return is;
}

template<typename First, typename ... Rest>
std::istream& get_impl(std::istream& is,  First& first, Rest&...rest)
{
  get_impl(is,first);
  get_impl(is,rest...);
  return is;
}

template< int ... Indexes, typename ... Args>
std::istream& get_helper(std::istream& is, index_tuple<Indexes...>, std::tuple<Args...>& tup)
{
  get_impl( is, std::get<Indexes>(tup)...);
  return is;
}


template<typename ...Args>
std::istream& operator>>(std::istream& is,  std::tuple<Args...>& tu)
{
  return get_helper(is,typename make_indexes<Args...>::type(),
                    tu);
}





using namespace io;

inline double log10_guard(double x)
{
  if (x==0) return 0;
  else return log10(x);

}




template<typename E>
struct D_logL
{
  M_Matrix<E> G;
  M_Matrix<E> H;
};



template<class E>
inline
std::ostream& operator<<(std::ostream& os,const D_logL<E>& d)
{
  os<<"G\n"<<d.G<<"\nH\n"<<d.H;
  return os;
}

template<class E>
inline
std::istream& operator>>(std::istream& is,D_logL<E>& d)
{
  std::string line;
  std::getline(is,line);
  is>>d.G;
  std::getline(is,line);
  std::getline(is,line);
  is>>d.H;
  return is;
}

template<class P>
struct mcmc_prior
{
  mcmc_prior();
  M_Matrix<P> param;
  double logPrior()const {return fullSum(logPriorLikelihood);}
  D_logL<P> D_prior;
  P logPriorLikelihood;
};

template <>
inline
mcmc_prior<double>::mcmc_prior():param(),D_prior(),logPriorLikelihood(std::numeric_limits<double>::quiet_NaN()){}


template <>
inline
mcmc_prior<M_Matrix<M_Matrix<double>>>::mcmc_prior():param(),D_prior(),logPriorLikelihood(){}



inline
mcmc_prior<double>
getElement(const mcmc_prior<M_Matrix<M_Matrix<double>>>& p, std::size_t j)
{
  mcmc_prior<double> out;
  out.logPriorLikelihood=p.logPriorLikelihood[1][j];
  out.param=p.param[1][j];
  out.D_prior.G=p.D_prior.G[1][j];
  out.D_prior.H=p.D_prior.H[1][j];
  return out;
}

inline
void setElement(mcmc_prior<M_Matrix<M_Matrix<double>>>& p,mcmc_prior<double> e, std::size_t j)
{
  p.param[1][j]=e.param;
  p.D_prior.G[1][j]=e.D_prior.G;
  p.D_prior.H[1][j]=e.D_prior.H;
  p.logPriorLikelihood[1][j]=e.logPriorLikelihood;
}


template<class E>
inline
std::ostream& operator<<(std::ostream& os,const mcmc_prior<E>& x)
{
  os<<"param\n"<<x.param<<"\nlogPrior\n"<<x.logPriorLikelihood<<"\nD_prior\n"<<x.D_prior;
  return os;
}


template<class E>
inline
std::istream& operator>>(std::istream& is, mcmc_prior<E>& x)
{
  std::string line;
  std::getline(is,line);
  is>>x.param;
  std::getline(is,line);
  std::getline(is,line);
  is>>x.logPrior;
  std::getline(is,line);
  std::getline(is,line);
  is>>x.D_prior;
  return is;
}

template<class T >struct myDts{};
template<>struct myDts<double>{
    typedef std::pair<std::vector<double>, std::vector<std::size_t>> type; };
template<>struct myDts<M_Matrix<M_Matrix<double>>>{
  typedef M_Matrix<std::pair<std::vector<double>, std::vector<std::size_t>>> type;};

template<class T> struct myL{};
template<>struct myL<double>{
  typedef double type;};
template<>struct myL<M_Matrix<M_Matrix<double>>>{
  typedef M_Matrix<double> type;};

template<class T>struct myTuple{};
template<>struct myTuple<double>{
  typedef std::tuple<double,std::size_t,double> type;};
template<>struct myTuple<M_Matrix<M_Matrix<double>>>{
  typedef M_Matrix<std::tuple<double,std::size_t,double>> type;};


template<class E, class L=typename myL<E>::type,class Dts=typename myDts<E>::type,
         class Tuple=typename myTuple<E>::type>
struct mcmc_post: public mcmc_prior<E>
{
  mcmc_post(const mcmc_prior<E>& p): mcmc_prior<E>(p), isValid(false), f(),logLikelihood(),vlogLikelihood(){}
  mcmc_post():mcmc_prior<E>(),isValid(false), f(),
    logLikelihood(),vlogLikelihood(){}
  bool isValid=false;
  M_Matrix<L> f;
  Dts dts;
  Tuple dtmin_Npoints_dtmax;
  double logLik()const {return fullSum(logLikelihood);}
  double slogLik()const {return std::sqrt(fullSum(vlogLikelihood));}


  L logLikelihood;
  L vlogLikelihood;
  double mlogbPL(double beta)const
  {
    return mcmc_prior<E>::logPrior()+logLik()*beta;
  }
  double mlogbPLikelihood(double beta, std::size_t i)const
  {
    return mcmc_prior<E>::logPriorLikelihood[i]+logLikelihood[i]*beta;
  }

  double logL(std::mt19937_64& mt)const
  {
    std::normal_distribution<double> logL_d(logLik(),slogLik());
    return logL_d(mt);

  }

  double logbPL(double beta, std::mt19937_64& mt)const {
    return mcmc_prior<E>::logPrior()+logL(mt)*beta;
  }
};


template<typename E,typename Dist>
struct mcmc_step;

template<typename Dist>
struct mcmc_step<double,Dist>;
template<class E>
struct mcmc_Dpost: public mcmc_post<E>
{

  mcmc_Dpost():mcmc_post<E>(),D_lik(){}
  mcmc_Dpost(mcmc_post<E> p): mcmc_post<E>(p), D_lik(){}
  D_logL<E> D_lik;

  double d_logLik_dBeta(double beta)const
  {
    auto Hbinv=inv(mcmc_post<E>::D_prior.H+beta*D_lik.H).first;
    auto Gb=mcmc_post<E>::D_prior.G+beta*D_lik.G;
    auto db=Gb*Hbinv;
    auto GdH=D_lik.G-db*D_lik.H;
    double s=xTSigmaX(GdH,Hbinv);
    //  for (std::size_t i=0; i<D_lik.H.nrows(); ++i)
    //     s+=sqr(D_lik.H(i,i));
    return s;
  }
};


template<typename Dist>
struct mcmc_step<double,Dist>: public mcmc_Dpost<double>
{
  mcmc_step(){}
  mcmc_step(mcmc_Dpost<double> p, double beta_, std::size_t iscout_): mcmc_Dpost<double>(p),beta(beta_),iscout(iscout_),proposed(){}
  double beta;
  std::size_t iscout;
  std::size_t dts_size()const {return dts.second.size();}


  double mlogbPL()const {return mcmc_post<double>::mlogbPL(beta);}
  double logbPL(std::mt19937_64& mt)const {return mcmc_post<double>::logbPL(beta,mt);}
  double logbPLb(double mybeta,std::mt19937_64& mt )const {return mcmc_post<double>::logbPL(mybeta,mt);}
  Dist proposed;


  static
  std::ostream& writelogLHeaderDataFrame(std::ostream& os)
  {
    os<<"beta\t";
    os<<"nscout\t";
    os<<"logPrior\t";
    os<<"logLik\t";
    os<<"slogLik\t";
    os<<"n_dts\t";
    os<<"logP";
    return os;
  }
  template<class M>
  static
  std::ostream& writeParamHeaderDataFrame(std::ostream& os, const M& model )
  {
    auto param=model.getPrior();
    param.writeHeaderDataFrame(os);
    return os;
  }

  template<class D,class M>
  std::ostream& writeSimulationHeaderDataFrame
  (std::ostream& os, const D& data,const M& model )const
  {

    auto sim=model.getSimulation(data.myExperiment(),param,dts);
    sim.writeHeaderDataFrame(os);
    return os;
  }
  template<class D,class M>
  std::ostream& writeFitHeaderDataFrame
  (std::ostream& os, const D& data,const M& model )const
  {

    auto sim=model.getSimulation(data.myExperiment(),param,dts);
    auto& cl=model.getLikelihood();
    cl.writeYfitHeaderDataFrame(data.myExperiment(),os,sim);
    return os;
  }


  std::ostream& writelogLRowDataFrame(std::ostream& os, const std::string& ss)
  {
    os<<ss;
    os<<beta<<"\t";
    os<<iscout<<"\t";
    os<<mcmc_Dpost<double>::logPriorLikelihood<<"\t";
    os<<mcmc_Dpost<double>::logLikelihood<<"\t";
    os<<mcmc_Dpost<double>::vlogLikelihood<<"\t";
    os<<mcmc_Dpost<double>::dts.first.size()<<"\t";
    os<<mlogbPL();
    return os;
  }

  std::ostream& writeYfitRowDataFrame(std::ostream& os, const std::string& s)
  {

    writelogLRowDataFrame(os,s);
    os<<"\t";
    for (std::size_t i=0; i+1< f.size(); ++i)
      os<< f[i]<<"\t";
    if (f.size()>0)
      os<< f[f.size()-1];
    return os;
  }
  template<class M>
  std::ostream& writeParamRowDataFrame(std::ostream& os, const M& model,const std::string& s)
  {
    writelogLRowDataFrame(os,s);
    os<<"\t";
    auto p=model.getParameter(param);
    p.writeRowDataFrame(os);
    return os;
  }

  template<class D,class M>
  std::ostream& writeSimulationRowDataFrame(std::ostream& os, const D& data,const M& model,const std::string& s)
  {
    writelogLRowDataFrame(os,s);
    os<<"\t";
    auto sim=model.getSimulation(data.myExperiment(),param,dts);
    sim.writeRowDataFrame(os);
    return os;
  }


};





template<class E>
class MultivariateGaussian
{
public:

  MultivariateGaussian(const M_Matrix<E> &mean,
                       const M_Matrix<E>& cov):
    mean_(mean),
    cov_(cov),
    covinv_(inv(cov).first),
    cho_cov_(chol(cov,"lower")),
    logDetCov_(logDiagProduct(cho_cov_))
  {
    assert(!cho_cov_.empty());

  }

  MultivariateGaussian(const M_Matrix<E>&mean
                       , const M_Matrix<E> &cov
                       , const M_Matrix<E> &covInv):
    mean_(mean),
    cov_(cov),
    covinv_(covInv),
    cho_cov_(chol(cov,"lower").first),
    logDetCov_(logDiagProduct(cho_cov_))
  {
    assert(!cho_cov_.empty());
  }






  MultivariateGaussian()=default;

  double logP(const M_Matrix<E> &x) const
  {
    if (mean_.size()>0)

      return -0.5*size()*log(PI)-logDetCov()-chi2(x);
    else return std::numeric_limits<double>::quiet_NaN();
  }
  double P(const M_Matrix<E>& x)const
  {
    return exp(logP(x));
  }

  void autoTest(std::mt19937_64& mt,std::size_t n)const
  {
    std::cerr<<"chi test n="<<size()<<" chis\n";
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
    std::cerr<<"\n chimean="<<chisum<<" chisqr="<<chisqr;
  }

  double chi2(const M_Matrix<E> &x) const
  {
    if (!mean_.empty())
      return 0.5*xTSigmaX(x-mean_,covinv_);
    else return std::numeric_limits<double>::quiet_NaN();
  }

  double operator()(const M_Matrix<E>& x)const
  {
    return logP(x);
  }
  M_Matrix<E> sample(std::mt19937_64& mt)const
  {
    M_Matrix<E> r;
    std::normal_distribution<> normal;
    if (this->size()>0)
      {
        auto z=Rand(mean_,normal,mt);
        r=mean_+multTransp(z,cho_cov_);
      }
    return r;
  }
  M_Matrix<E> operator()(std::mt19937_64& mt)const
  {
    return sample(mt);
  }
  friend
  std::istream& operator>>(std::istream& is, MultivariateGaussian& x)
  {
    std::string line;
    std::getline(is,line);
    std::getline(is,line);
    is>>x.mean_;
    std::getline(is,line);
    std::getline(is,line);
    is>>x.cov_;
    return is;
  }


  const M_Matrix<E>& Mean()const
  {
    return mean_;
  }
  M_Matrix<E> Cov()const
  {
    return inv(covinv_).first;
  }

  const M_Matrix<E>& CovInv()const
  {
    return covinv_;
  }




  const M_Matrix<E>& Chol()const{return cho_cov_;}

  double logDetCov()const {return logDetCov_;}

  virtual M_Matrix<E> logPGradient(const M_Matrix<E>& x)const
  {
    return M_Matrix<E>(1u,x.size(),(x-Mean())*Cov());
  }
  virtual M_Matrix<E> logPHessian(const M_Matrix<E>& )const
  {
    return Cov();
  }



  std::size_t size()const
  {
    return mean_.size();
  }

  bool isValid()const
  {
    return !cho_cov_.empty();
  }


  MultivariateGaussian(const MultivariateGaussian& other)=default;
  MultivariateGaussian& operator=(const MultivariateGaussian& other)=default;


  ~MultivariateGaussian(){};

private:
  M_Matrix<E> mean_;
  M_Matrix<E> cov_;
  M_Matrix<E> covinv_;
  M_Matrix<E> cho_cov_;
  double logDetCov_;

  // Distribution interface
public:
  virtual std::__cxx11::string myClass() const
  {
    return "MultivariateGaussian";
  }
};



template<class E>
std::ostream& operator<<(std::ostream& os,const MultivariateGaussian<E>& x)
{
  os<<"\nMean\n"<<x.Mean();
  os<<"\nCov\n"<<x.Cov();
  return os;
}

struct trust_region
{
  static trust_region min(){return {1E-6};}

  static trust_region max(){return {1.0};}

  static double logFactor(){return std::log(100.0);}

  double r_;

  double getValue()const {return r_;}
  void setValue(double r) { r_=r;}


  bool operator()(double expected, double found, std::size_t k)const
  {
    if (!std::isnan(found)&&(sqr(expected-found)<2.0*r_*k))
      return true;
    else
      return false;
  }

};


class Landa
{
public:
  typedef std::pair<Landa,double> myParameter;


  static std::string ParName(std::size_t i)
  {
    switch (i)
      {
      case 0: return "Landa_50";
      case 1: return "Hill_Coeff";
      default: return "";
      }
  }


  static std::string ClassName(){return "Landa";}
  static trust_region min(){return {0.0};}

  static trust_region max(){return {1E9};}

  static double logFactor(){return std::log(10.0);}



  double landa_;

  double getValue()const {return landa_;}
  void setValue(double r) { landa_=r;}


  struct myAcceptProb
  {
    double operator()(const Landa& landa,const myParameter& param)const
    {
      double landa50=param.first.getValue();
      double h=param.second;
      return 1.0/(1.0+std::pow(landa50/(landa.getValue()+1.0),h));
    }
    double operator()(const myParameter& param, const Landa& landa)const
    {
      return operator()(landa,param);
    }
  };
  typedef myAcceptProb AcceptanceProbability;

  struct myExpectVelocity
  {
    double operator()(const Landa& landa)const
    {
      return (1.0/(1.0+landa.getValue()));
    }
  };

  typedef myExpectVelocity ExpectedVelocity;

  static std::map<myParameter, double> uniform_parameter_prior(const  std::vector<std::vector<double>>& v, double p=-1)
  {
    std::map<myParameter, double> out;
    if (p==-1)
      p=1.0/(v[0].size()*v[1].size());
    for (std::size_t i=0; i<v[0].size(); ++i)
      for (std::size_t j=0; j<v[1].size(); ++j)
        out[{Landa{v[0][i]},v[1][j]}]+=p;
  return out;
}



};


inline bool operator<(const Landa& one,const Landa& two)
{ return one.getValue()<two.getValue();}


struct Likelihood_Record{
  struct Particle
  {
    bool isValid=false;
    double logLikInit;
    double logLikNext;
  };

  std::size_t size()const {return particle.size();}
  std::vector<std::vector<Particle>> particle;
  void newParticles()
  {
    std::vector<Particle> p(desc_beta.size());
    particle.push_back(p);
  }

  friend std::ostream& operator<<(std::ostream& os, const Particle& me)
  {
    os<<me.isValid<<"\t"<<me.logLikInit<<"\t"<<me.logLikNext<<"\t";
    return os;
  }
  friend std::istream& operator>>(std::istream& is,  Particle& me)
  {
    is>>me.isValid>>me.logLikInit>>me.logLikNext;
    return is;
  }


  std::vector<double> desc_beta;

  friend std::ostream& operator<<(std::ostream& os, const Likelihood_Record& me)
  {
    os<<me.particle;
    os<<me.desc_beta;
    return os;
  }
  friend std::istream& operator>>(std::istream& is, Likelihood_Record& me)
  {
    is>>me.particle;
    is>>me.desc_beta;
    return is;
  }


};


class Beta

{
public:
  static std::string ClassName(){return "Beta";}


  Beta(std::size_t n, double med_beta, std::size_t n2=4, double min_beta=1e-7):asc_beta_(n+n2)
  {
    double f=std::pow(med_beta,1.0/(n-1));
    double f2=std::pow(med_beta/min_beta,1.0/n2);

    for (std::size_t i=0; i<n; ++i)
      asc_beta_[i]=(std::pow(f,i)+(n-i-1)*med_beta)/(1.0+(n-i-1)*med_beta);
    for (std::size_t i=0; i<n2; ++i)
      asc_beta_[n+i]=med_beta/(std::pow(f2,i+1));


  }

  Beta():asc_beta_(){}

  Beta(std::vector<double> beta):asc_beta_(beta){}


  std::vector<double>const & getValue()const {return asc_beta_;}

  std::vector<double>& getValue() {return asc_beta_;}

  std::size_t size()const {return asc_beta_.size();}

  static
  std::vector<double> resize(std::size_t newsize, const std::vector<double>& beta)
  {
    std::vector<double> newBeta(newsize);
    double f=1.0*beta.size()/newsize;
    for (std::size_t i=0; i<newsize; ++i)
      newBeta[i]=interpolate(i*f,beta);
    return newBeta;
  }

private:
  std::vector<double> asc_beta_;

  static double interpolate(double n, const std::vector<double>& beta)
  {
    std::size_t i=std::floor(n);
    double f=n-i;
    if (i+1<beta.size())
      return beta[i]+f*beta[i+1];
    else return beta[i];


  }
};

inline
std::ostream& operator<<(std::ostream& os, const Beta& t)
{
  os<<t.getValue();
  return os;
}

inline
std::istream& operator>>(std::istream& is, Beta& t)
{
  std::vector<double> betas;
  is>>betas;
  t.getValue()=betas;
  return is;
}




class Master_Adaptive_Beta_New
{
public:




  void push_step()
  {
    data_.back().newParticles();
  }

  template<class mcmc>
  void new_acceptance(std::size_t i, const mcmc& sDist,const mcmc& cDist)
  {
    this->data_.back().particle.back()[i].isValid=true;
    this->data_.back().particle.back()[i].logLikInit=sDist.logLik();
    this->data_.back().particle.back()[i].logLikNext=cDist.logLik();
  }
  template<class mcmc>
  void new_rjection(std::size_t i, const mcmc& sDist)
  {
    this->data_.back().particle.back()[i].isValid=false;
    this->data_.back().particle.back()[i].logLikInit=sDist.logLik();
  }


  void push_acceptance(double betaless, double betamore)
  {
    ++accepts_[{betaless,betamore}].first;

  }

  void push_rejection(double betaless, double betamore)
  {
    ++accepts_[{betaless,betamore}].second;
  }




  void reset()
  {
    data_.clear();
  }

  Master_Adaptive_Beta_New(std::size_t Ninitial,  double beta_min, std::size_t N_2, double beta_infimo): desc_beta_{Ninitial,beta_min, N_2, beta_infimo}, data_{}
  {
    Likelihood_Record r;
    r.desc_beta=desc_beta_.getValue();
    data_.push_back(r);
  }




  std::size_t size()const {return desc_beta_.size();}

  Beta const& getBeta()const {return desc_beta_;}

  friend
  std::ostream& operator<<(std::ostream& os, const Master_Adaptive_Beta_New& me)
  {
    os<<"Beta\n"<<me.desc_beta_;
    os<<"accepts\n";
    os<<me.accepts_;
    os<<"data\n";
    os<<me.data_;
    return os;

  }

  friend
  std::istream& operator>>(std::istream& os, Master_Adaptive_Beta_New& me)
  {
    std::string line;
    std::getline(os,line);
    os>>me.desc_beta_;
    std::getline(os,line);
    os>>me.accepts_;
    std::getline(os,line);
    os>>me.data_;
    return os;

  }




private:

  Beta desc_beta_;
  std::map<std::pair<double,double>,std::pair<std::size_t,std::size_t>> accepts_;

  std::vector<Likelihood_Record>  data_;

};




template<typename E>
struct LM_MultivariateGaussian: public MultivariateGaussian<E>
{
  LM_MultivariateGaussian(MultivariateGaussian<E> m
                          , Landa mylanda
                          ,double my_exp_delta_next_logL):
    MultivariateGaussian<E>(m),
    landa(mylanda),
    exp_delta_next_logL(my_exp_delta_next_logL)
  {}
  LM_MultivariateGaussian():landa(),exp_delta_next_logL(){}
  Landa landa;
  double exp_delta_next_logL;
};

template<class E>
std::ostream& operator<<(std::ostream& os,const LM_MultivariateGaussian<E>& x)
{
  const MultivariateGaussian<E>& b=x;
  os<<b;
  os<<"\nlanda\n"<<x.landa;
  os<<"\nexp_delta_next_logL\n"<<x.exp_delta_next_logL<<"\n";
  return os;
}


template<class D,  class M>
class Poisson_Likelihood
{
public:



  static double logLikelihood(double landa,std::size_t k)
  {
    return k*log(landa)-landa-lgamma(k+1);
  }

  static
  M_Matrix<double> sample(const M& model,const D& ,std::mt19937_64& mt)
  {

    return model.sample(mt);
  }




  static double logL(const D& data,const M_Matrix<double>& landa)
  {
    M_Matrix<std::size_t> k=data();
    if (landa.empty()) return std::numeric_limits<double>::quiet_NaN();
    double sumLogL=0;
    for (std::size_t i=0; i<k.nrows(); ++i)
      {
        for (std::size_t j=0; j<k.ncols(); ++j)
          {
            if (std::isnan(landa(i,j)))
              return landa(i,j);
            else
              if (landa(i,j)!=0)
                sumLogL+=logLikelihood(landa(i,j),k(i,j));
              else if(k(i,j)!=0)
                return std::numeric_limits<double>::quiet_NaN();
          }
      }
    return sumLogL;
  }

  static mcmc_post<double> get_mcmc_Post(const M& model, const D& data, M_Matrix<double> param)
  {


    mcmc_prior<double> p=model.prior(data,param);
    mcmc_post<double> out(p);
    out.f=model.f(data,param,out.dts);
    if (out.f.size()==0)
      out.isValid=false;
    else
      {

        out.logLikelihood=logL(data,out.f);
        if (std::isnan(out.logLikelihood))
          out.isValid=false;
        else out.isValid=true;
      }
    return out;
  }

  static mcmc_post<double> get_mcmc_Post(const M& model,
                                         const D& data,
                                         M_Matrix<double> param,
                                         double slogL_max,
                                         std::size_t ndts_max_1)
  {
    std::size_t ndts_max_0=ndts_max_1/2;
    double vlogL_max=sqr(slogL_max);
    mcmc_prior<double> p=model.prior(data,param);
    mcmc_post<double> out(p);
    double dtmin_0=0, dtmin_1;
    std::size_t n_per10_0, n_per10_1;
    double dtmax_0, dtmax_1;
    std::pair<std::vector<double>, std::vector<std::size_t>> dts_0, dts_1;
    auto f0=model.f(data,param,dtmin_0,n_per10_0,dtmax_0, ndts_max_0,dts_0);
    double logLik0=logL(data,f0);
    bool ishope=dts_0.first.size()>0&&dts_0.first.size()<ndts_max_0;
    bool firstValid=false;
    std::size_t iloop=0;
    std::size_t maxloop=10;
    while (ishope&&!firstValid&& iloop<maxloop)
      {
        ++iloop;
        dtmin_0/=2;
        n_per10_0*=2;
        dtmax_0/=2;
        f0=model.f(data,param,dtmin_0,n_per10_0,dtmax_0, ndts_max_0,dts_0);
        logLik0=logL(data,f0);
        ishope=dts_0.first.size()>0&&dts_0.first.size()<ndts_max_0;
        firstValid=std::isfinite(logLik0);
      }

    if (!ishope )
      {
        out.isValid=false;
        out.logLikelihood=logLik0;
        out.dts=dts_0;
        out.f=f0;

      }
    else
      {
        dtmin_1=dtmin_0/2;
        n_per10_1=n_per10_0*2;
        dtmax_1=dtmax_0/2;
        auto f1=model.f(data,param,dtmin_1,n_per10_1,dtmax_1,ndts_max_1,dts_1);
        double logLik1=logL(data,f1);
        double vlogLik=sqr(logLik0-logLik1);
        bool hav_logL1=std::isfinite(logLik1);
        bool have_slogL=std::isfinite(vlogLik);
        bool good_slogL=vlogLik<vlogL_max;
        bool exceed_ndts=dts_1.first.size()>ndts_max_0;
        iloop=0;
        while (hav_logL1&&have_slogL&&!good_slogL&&!exceed_ndts&&iloop<maxloop)
          {
            ++iloop;
            f0=std::move(f1);
            dts_0=std::move(dts_1);

            logLik0=logLik1;
            dtmin_0=dtmin_1;
            n_per10_0=n_per10_1;

            dtmin_1=dtmin_0/2;
            n_per10_1=n_per10_0*2;
            dtmax_1=dtmax_0/2;

            f1=model.f(data,param,dtmin_1,n_per10_1,dtmax_1,ndts_max_1,dts_1);
            logLik1=logL(data,f1);
            hav_logL1=std::isfinite(logLik1);
            if (hav_logL1)
              vlogLik=sqr(logLik0-logLik1);
            else
              vlogLik/=4;
            have_slogL=std::isfinite(vlogLik);
            good_slogL=vlogLik<vlogL_max;
            exceed_ndts=dts_1.first.size()>ndts_max_0;

          }

        out.vlogLikelihood=vlogLik*4;
        out.logLikelihood=logLik0;
        out.f=f0;
        out.dtmin_Npoints_dtmax={dtmin_0,n_per10_0,dtmax_0};
        out.dts=dts_0;
        if (!std::isfinite(out.logLikelihood)
            ||out.vlogLikelihood>std::sqrt(std::abs(out.logLikelihood)))
          out.isValid=false;
        else out.isValid=true;
      }


    return out;
  }


};

template<class D,  class M>
class Poisson_DLikelihood: public Poisson_Likelihood<D,M>
{
public:
  typedef double  E;


  static mcmc_Dpost<double> get_mcmc_Dpost(const M& model, const D& data, const M_Matrix<double>& param,double slogL_max,std::size_t ndts_max)
  {
    return get_mcmc_Dpost(model,data,param,get_mcmc_Post(model,data,param,slogL_max,ndts_max));
  }
  static mcmc_post<double> get_mcmc_Post(const M& model, const D& data, M_Matrix<double> param,double slogL_max,std::size_t ndts_max)
  {
    return Poisson_Likelihood<D,M>::get_mcmc_Post(model,data,param,slogL_max,ndts_max);
  }


  static mcmc_Dpost<double> get_mcmc_Dpost(const M& model, const D& data, const M_Matrix<double>& param, mcmc_post<double> p)
  {
    mcmc_Dpost<double> out(p);
    if (out.isValid)
      {
        M_Matrix<std::size_t> k=data();
        M_Matrix<double> logLanda_0=out.f.apply([](double x)
        {return log10_guard(x);});
        if (isnan(logLanda_0))
          out.isValid=false;
        M_Matrix<double> J=get_J(model,  data, param,logLanda_0, out.dts  );
        if (J.size()==0)
          out.isValid=false;
        else
          {
            out.D_lik.G=get_G(out.f,k,J);
            out.D_lik.H=get_H(out.f,J);
          }
      }
    return out;
  }

private:
  static
  M_Matrix<double> get_G(const M_Matrix<double>& landa, const M_Matrix<std::size_t>& k, const M_Matrix<double>& J)
  {
    M_Matrix<double> out(1,J.ncols(),0.0);
    for (std::size_t j=0; j<J.ncols(); ++j)
      for (std::size_t i=0; i<landa.size(); ++i)
        out[j]+=(landa[i]-k[i])*J(i,j);
    return out;
  }

  static
  M_Matrix<double> get_H(const M_Matrix<double>& landa, const M_Matrix<double>& J)
  {
    std::size_t n=landa.size();
    std::size_t npar=J.ncols();
    M_Matrix<double> out(npar,npar,M_Matrix<double>::SYMMETRIC,0.0);
    for (std::size_t j=0; j<npar; ++j)
      for (std::size_t j2=j; j2<npar; ++j2)
        for (std::size_t i=0; i<n; ++i)
          out(j,j2)+=landa[i]*J(i,j)*J(i,j2);

    //auto test=out-Transpose(J)*diag(landa.toVector_of_Rows())*J;
    return out;

  }

  static
  M_Matrix<double>
  get_J(const M& model, const D& data, const M_Matrix<double>& param,
        const M_Matrix<double>& logLanda_0 , std::pair<std::vector<double>,std::vector<std::size_t>>& dts,
        double delta=1e-5, double delta_div=10, double deltamin=1e-7)
  {
    M_Matrix<double> k=data();
    std::size_t n=k.size();
    std::size_t npar=param.size();
    M_Matrix<double> out(n,npar,0.0);

    M_Matrix<double> logLanda2=model.logLanda(data,param,dts);
    auto logLandadif=logLanda_0-logLanda2;

    double maxdif=maxAbs(logLandadif);
    std::cerr<<"max dif="<<maxdif<<"\n";

    for (std::size_t j=0; j<npar; ++j)
      {
        double deltarun=delta;
        M_Matrix<double> logLanda_0_run=logLanda_0;
        auto dts_run=dts;

        M_Matrix<double> p(param);
        p[j]+=deltarun;
        M_Matrix<double> logLanda_i=model.logLanda(data,p,dts);
        while ((isnan(logLanda_i)||(logLanda_i.empty()))&&deltarun>deltamin)
          {
            deltarun=deltarun/delta_div;
            p=param;
            p[j]+=deltarun;
            logLanda_i=model.logLanda(data,p,dts);
          }
        if (isnan(logLanda_i)||logLanda_i.empty())
          {
            return {};
          }

        for (std::size_t i=0; i<n; ++i)
          out(i,j)=(logLanda_i[i]-logLanda_0[i])/deltarun;

      }
    return out;
  }

};




template<class D,  class M>
class Normal_Likelihood
{
public:



  static double logLikelihood(double landa,double k,double se2)
  {
    if (se2!=0)
      return -(sqr(landa-k)/(landa+se2)+std::log(2*PI*se2))/2.0;
    else
      return k*log(landa)-landa-lgamma(k+1);

  }

  static
  M_Matrix<double> sample(const M& model,const D& ,std::mt19937_64& mt)
  {

    return model.sample(mt);
  }




  static double logL(const D& data,const M_Matrix<double>& landa)
  {
    M_Matrix<double> k=data();
    M_Matrix<double> se2=data.se2();
    if (landa.empty()) return std::numeric_limits<double>::quiet_NaN();
    double sumLogL=0;
    for (std::size_t i=0; i<k.nrows(); ++i)
      {
        for (std::size_t j=0; j<k.ncols(); ++j)
          {
            if (std::isnan(landa(i,j)))
              return landa(i,j);
            else
              if ((landa(i,j)!=0)&&std::isfinite(se2(i,j)))
                sumLogL+=logLikelihood(landa(i,j),k(i,j), se2(i,j));
              else if ((k(i,j)!=0)&& std::isfinite(se2(i,j)))
                return std::numeric_limits<double>::quiet_NaN();
          }
      }
    return sumLogL;
  }

  static mcmc_post<double> get_mcmc_Post(const M& model, const D& data, M_Matrix<double> param)
  {


    mcmc_prior<double> p=model.prior(data,param);
    mcmc_post<double> out(p);
    out.f=model.f(data,param,out.dts);
    if (out.f.size()==0)
      out.isValid=false;
    else
      {

        out.logLikelihood=logL(data,out.f);
        if (std::isnan(out.logLikelihood))
          out.isValid=false;
        else out.isValid=true;
      }
    return out;
  }

  static mcmc_post<double> get_mcmc_Post(const M& model,
                                         const D& data,
                                         M_Matrix<double> param,
                                         double slogL_max,
                                         std::size_t ndts_max_1)
  {
    std::size_t ndts_max_0=ndts_max_1/2;
    double vlogL_max=sqr(slogL_max);
    mcmc_prior<double> p=model.prior(data,param);
    mcmc_post<double> out(p);
    double dtmin_0=0, dtmin_1;
    std::size_t n_per10_0, n_per10_1;
    double dtmax_0, dtmax_1;
    std::pair<std::vector<double>, std::vector<std::size_t>> dts_0, dts_1;
    auto f0=model.f(data,param,dtmin_0,n_per10_0,dtmax_0, ndts_max_0,dts_0);
    double logLik0=logL(data,f0);
    bool ishope=dts_0.first.size()>0&&dts_0.first.size()<ndts_max_0;
    bool firstValid=false;
    std::size_t iloop=0;
    std::size_t maxloop=10;
    while (ishope&&!firstValid&& iloop<maxloop)
      {
        ++iloop;
        dtmin_0/=2;
        n_per10_0*=2;
        dtmax_0/=2;
        f0=model.f(data,param,dtmin_0,n_per10_0,dtmax_0, ndts_max_0,dts_0);
        logLik0=logL(data,f0);
        ishope=dts_0.first.size()>0&&dts_0.first.size()<ndts_max_0;
        firstValid=std::isfinite(logLik0);
      }

    if (!ishope )
      {
        out.isValid=false;
        out.logLikelihood=logLik0;
        out.dts=dts_0;
        out.f=f0;

      }
    else
      {
        dtmin_1=dtmin_0/2;
        n_per10_1=n_per10_0*2;
        dtmax_1=dtmax_0/2;
        auto f1=model.f(data,param,dtmin_1,n_per10_1,dtmax_1,ndts_max_1,dts_1);
        double logLik1=logL(data,f1);
        double vlogLik=sqr(logLik0-logLik1);
        bool hav_logL1=std::isfinite(logLik1);
        bool have_slogL=std::isfinite(vlogLik);
        bool good_slogL=vlogLik<vlogL_max;
        bool exceed_ndts=dts_1.first.size()>ndts_max_0;
        iloop=0;
        while (hav_logL1&&have_slogL&&!good_slogL&&!exceed_ndts&&iloop<maxloop)
          {
            ++iloop;
            f0=std::move(f1);
            dts_0=std::move(dts_1);

            logLik0=logLik1;
            dtmin_0=dtmin_1;
            n_per10_0=n_per10_1;

            dtmin_1=dtmin_0/2;
            n_per10_1=n_per10_0*2;
            dtmax_1=dtmax_0/2;

            f1=model.f(data,param,dtmin_1,n_per10_1,dtmax_1,ndts_max_1,dts_1);
            logLik1=logL(data,f1);
            hav_logL1=std::isfinite(logLik1);
            if (hav_logL1)
              vlogLik=sqr(logLik0-logLik1);
            else
              vlogLik/=4;
            have_slogL=std::isfinite(vlogLik);
            good_slogL=vlogLik<vlogL_max;
            exceed_ndts=dts_1.first.size()>ndts_max_0;

          }

        out.vlogLikelihood=vlogLik*4;
        out.logLikelihood=logLik0;
        out.f=f0;
        out.dtmin_Npoints_dtmax={dtmin_0,n_per10_0,dtmax_0};
        out.dts=dts_0;
        if (!std::isfinite(out.logLikelihood))
          out.isValid=false;
        else out.isValid=true;
      }


    return out;
  }


};




template<
    class D
    ,  class M
    , class D_Lik=Poisson_DLikelihood<D,M>
    , class myPropDistribution=LM_MultivariateGaussian<double>
    , class AP=trust_region
    >
class LevenbergMarquardt_step
{
public:

  double landa0_;
  double v_;
  std::size_t maxLoop_;

  //  double optimal_acceptance_rate_=0.574; // Handbook of Markov Chain Montecarlo p.100
  // optimal acceptance rate for Metropolis-Adjusted Langevin Algorithm

  double optimal_acceptance_rate_=0.234; // Handbook of Markov Chain Montecarlo p.96
  // optimal acceptance rate for Metropolis Algotithm


public:

  double optimal_acceptance_rate()const {return optimal_acceptance_rate_;}

  LevenbergMarquardt_step( double landa0,
                           double v,
                           std::size_t maxLoop)
    :landa0_{landa0},v_(v),maxLoop_(maxLoop)
  {
  }





  template<class E>
  struct LM_logL:public D_logL<E>
  {
    double logL;
    M_Matrix<E> Hinv;
    M_Matrix<E> d;
    double exp_delta_next_logL;
  };


  template<class E>
  static LM_logL<E> update_landa(const mcmc_Dpost<E>& postL,double landa,double beta)
  {
    LM_logL<E> out;
    out.G=postL.D_prior.G+postL.D_lik.G*beta;
    out.H=postL.D_prior.H+postL.D_lik.H*beta;
    out.logL=postL.mlogbPL(beta);
    for (std::size_t i=0; i<out.H.nrows(); ++i)
      //    out.H(i,i)=out.H(i,i)+postL.D_lik.H(i,i)*beta*landa;
      // this alternative does not work
      out.H(i,i)=out.H(i,i)*(1+landa);
    out.Hinv=invSafe(out.H);
    if (!out.Hinv.empty())
      {
        out.d=-(out.G*out.Hinv);
        out.exp_delta_next_logL=-0.5*multTransp(out.d,out.G)[0];
      }
    return out;
  }

  template<class E>
  mcmc_step<E,myPropDistribution> get_mcmc_step(const D_Lik L, const M& model, const D& data, const M_Matrix<E>& param, double beta)const
  {
    return get_mcmc_step(L,model,data,L.get_mcmc_Dpost(model,data,param),beta);
  }

  template<class E>
  mcmc_step<E,myPropDistribution> get_mcmc_step(const D_Lik L, const M& model, const D& data, const M_Matrix<E>& param, double beta,double landa)const
  {
    return get_mcmc_step(L,model,data,L.get_mcmc_Dpost(model,data,param),beta,landa);
  }

  template<class E>
  mcmc_step<E,myPropDistribution> get_mcmc_step(const D_Lik L, const M& model, const D& data, mcmc_Dpost<E> p,double beta)const
  {
    mcmc_step<E,myPropDistribution> out(p,beta);
    return update_mcmc_step(L,model,data,out,beta);
  }

  template<class E>
  mcmc_step<E,myPropDistribution> get_mcmc_step(const D_Lik L, const M& model, const D& data, mcmc_Dpost<E> p,double beta,double landa)const
  {
    mcmc_step<E,myPropDistribution> out(p,beta);
    return update_mcmc_step(L,model,data,out,beta,landa);
  }


  template<class E>
  mcmc_step<E,myPropDistribution> get_mcmc_step(const D_Lik L, const M& model, const D& data, const M_Matrix<E>& param, mcmc_post<E> p,double beta)const
  {
    return get_mcmc_step(L,model,data,L.get_mcmc_Dpost(model,data,param,p),beta);
  }

  template<class E>
  mcmc_step<E,myPropDistribution> get_mcmc_step(const D_Lik L, const M& model, const D& data, const M_Matrix<E>& param, mcmc_post<E> p,double beta,double landa)const
  {
    return get_mcmc_step(L,model,data,L.get_mcmc_Dpost(model,data,param,p),beta,landa);
  }



  template<class E>
  mcmc_step<E,myPropDistribution>& update_mcmc_step(const D_Lik L, const M& model, const D& data, mcmc_step<E,myPropDistribution>& out,double beta, AP t_,double slogL_max,std::size_t ndts_max)const
  {
    out.beta=beta;


    double landa=0;
    std::size_t iloop=0;
    if (out.isValid)
      {
        LM_logL<E> LM_logL=update_landa( out,landa,beta);
        M_Matrix<double> next;
        mcmc_post<E> cand;
        if (!LM_logL.Hinv.empty())
          {
            next=out.param+LM_logL.d;
            cand=L.get_mcmc_Post(model,data,next,slogL_max,ndts_max);
          }
        while
            ((LM_logL.Hinv.empty()
              ||!t_(LM_logL.exp_delta_next_logL
                    ,cand.mlogbPL(beta)-out.logbPL()
                    ,out.param.size())
              )&&iloop<maxLoop_)

          {
            landa=landa0_*std::pow(v_,iloop);
            LM_logL=update_landa( out,landa, beta);
            next=out.param+LM_logL.d;
            cand=L.get_mcmc_Post(model,data,next,slogL_max,ndts_max);
            ++iloop;
          }
        if (LM_logL.Hinv.empty())
          out.isValid=false;
        else
          out.proposed=myPropDistribution(MultivariateGaussian<E>
                                          (next,LM_logL.Hinv,LM_logL.H),landa);
      }
    return
        out;
  }


  template<class E>
  mcmc_step<E,myPropDistribution>& update_mcmc_step(const D_Lik /*L*/, const M& /*model*/, const D& /* data*/, mcmc_step<E,myPropDistribution>& out,double beta, double landa)const
  {
    out.beta=beta;
    if (out.isValid)
      {
        LM_logL<E> LM_logL=update_landa( out,landa,beta);
        if (!LM_logL.Hinv.empty())
          {
            M_Matrix<E> next=out.param+LM_logL.d;
            out.proposed=myPropDistribution(MultivariateGaussian<E>
                                            (next,LM_logL.Hinv,LM_logL.H),landa);
          }
        else
          out.isValid=false;
      }
    return
        out;
  }




};



template<
    class D
    ,class M
    , class D_Lik
    , class myPropDistribution>
class LevenbergMarquardt_step<D,M,D_Lik,myPropDistribution,Landa>
{
public:


  //  double optimal_acceptance_rate_=0.574; // Handbook of Markov Chain Montecarlo p.100
  // optimal acceptance rate for Metropolis-Adjusted Langevin Algorithm

  //  double optimal_acceptance_rate_=0.234; // Handbook of Markov Chain Montecarlo p.96
  // optimal acceptance rate for Metropolis Algotithm


public:






  template <class E>
  struct LM_logL:public D_logL<E>
  {
    double logL;
    M_Matrix<E> Hinv;
    M_Matrix<E> d;
    double exp_delta_next_logL;
  };
  template <class E>
  static LM_logL<E> update_landa(const mcmc_Dpost<E>& postL,const Landa& landa,double beta)
  {
    LM_logL<E> out;
    out.G=postL.D_prior.G+postL.D_lik.G*beta;
    out.H=postL.D_prior.H+postL.D_lik.H*beta;

    /*uto Hp=M_Matrix<E>(M_Matrix<E>::FULL,postL.D_prior.H);
        auto Ht=M_Matrix<E>(M_Matrix<E>::FULL,postL.D_lik.H);


        auto test=out.H-(Hp+Ht*beta);
        */
    out.logL=postL.mlogbPL(beta);
    for (std::size_t i=0; i<out.H.nrows(); ++i)
      //    out.H(i,i)=out.H(i,i)+postL.D_lik.H(i,i)*beta*landa;
      // this alternative does not work
      out.H(i,i)=out.H(i,i)*(1+landa.getValue());
    out.Hinv=inv(out.H).first;
    if (!out.Hinv.empty())
      {
        out.d=-(out.G*out.Hinv);
        out.exp_delta_next_logL=-0.5*fullSum(multTransp(out.d,out.G));
      }
    return out;
  }

  template<class E>
  static mcmc_step<E,myPropDistribution> get_mcmc_step(const D_Lik L, const M& model, const D& data, const M_Matrix<E>& param,const Landa& landa, double beta, std::size_t iscout, double slogL_max, std::size_t ndts_max)
  {
    return get_mcmc_step(L,model,data,L.get_mcmc_Dpost(model,data,param,slogL_max,ndts_max),landa,beta, iscout);
  }


  template<class E>
  static mcmc_step<E,myPropDistribution> get_mcmc_step(const D_Lik L, const M& model, const D& data, mcmc_Dpost<E> p,const Landa& landa,double beta, std::size_t iscout)
  {
    mcmc_step<E,myPropDistribution> out(p,beta,iscout);
    return update_mcmc_step(L,model,data,out,landa,beta);
  }

  template<class E>
  static mcmc_step<E,myPropDistribution> get_mcmc_step(const D_Lik L, const M& model, const D& data, const M_Matrix<E>& param, mcmc_post<E> p,const  Landa& landa,double beta, std::size_t iscout)
  {
    return get_mcmc_step(L,model,data,L.get_mcmc_Dpost(model,data,param,p),landa ,beta, iscout);
  }





  template<class E>
  static  mcmc_step<E,myPropDistribution> update_mcmc_step(const D_Lik , const M& , const D& , mcmc_step<E,myPropDistribution> out,const Landa& landa,double beta)
  {
    out.beta=beta;
    if (out.isValid)
      {
        LM_logL<E> LM_logL=update_landa( out,landa,beta);
        if (!LM_logL.Hinv.empty())
          {
            M_Matrix<E> next=out.param+LM_logL.d;
            out.proposed=myPropDistribution
                (MultivariateGaussian<E>
                 (next,LM_logL.Hinv,LM_logL.H),
                 landa,
                 LM_logL.exp_delta_next_logL);
            if (!out.proposed.isValid())
              {
                out.isValid=false;
                std::cerr<<"\ninvalid cholesky\n";
              }
          }
        else
          out.isValid=false;
      }
    return
        out;
  }
};



template<class mcmc>
class Tempered_Evidence_Evaluation
{
public:

  static double Evidence_pair(const std::pair<Beta,std::vector<mcmc>>& sample)
  {
    return Evidence(sample.first,sample.second);
  }
  static double Evidence(const Beta& mybeta,const std::vector<mcmc>& sample)
  {
    auto n=mybeta.getValue().size();
    double beta0=mybeta.getValue()[0];
    double sum=0;
    double sumdb=0;
    double logLik0=sample[0].logLik();
    for (std::size_t i=1; i<sample.size(); ++i)
      {
        double beta=mybeta.getValue()[i];
        double db=beta0-beta;
        double logLik=sample[i].logLik();
        sum+=db*(logLik0+logLik)/2;
        sumdb+=db;
        if ((i==n-1)&&(beta>0))
          {
            double db0=beta;
            sum+=(db0*(1+db0/db/2))*logLik-sqr(db0)/db/2*logLik0;
            sumdb+=db0;
          }
        beta0=beta;
        logLik0=logLik;
      }
    return sum;

  }


  static double Evidence(std::size_t i0,const Beta& mybeta,const std::vector<mcmc>& sample)
  {
    //  el verdadero logLik es logLik* beta[i0]
    // el verdadero beta[i] es beta[i]/beta[i0]
    double betaRef=mybeta.getValue()[i0];
    auto n=mybeta.getValue().size();
    double beta0=mybeta.getValue()[i0]/betaRef;
    double sum=0;
    double sumdb=0;
    double logLik0=sample[i0].logLik()*betaRef;
    if (i0+1==sample.size())return logLik0;
    else
      {
        for (std::size_t i=i0+1; i<sample.size(); ++i)
          {
            double beta=mybeta.getValue()[i]/betaRef;
            double db=beta0-beta;
            double logLik=sample[i].logLik()*betaRef;
            sum+=db*(logLik0+logLik)/2;
            sumdb+=db;
            if ((i==n-1)&&(beta>0))
              {
                double db0=beta;
                sum+=(db0*(1+db0/db/2))*logLik-sqr(db0)/db/2*logLik0;
                sumdb+=db0;
              }
            beta0=beta;
            logLik0=logLik;
          }
        return sum;
      }
  }

  static std::vector<double> Evidence_of_beta(const Beta& mybeta,const std::vector<mcmc>& sample)
  {
    std::vector<double> out(mybeta.getValue().size());
    for (std::size_t i=0; i<out.size(); ++i)
      {
        out[i]=Evidence(i,mybeta,sample);
      }
    return out;
  }

};





template
<   class Adaptive_parameterized,
    class D
    , class M
    ,  class D_Lik//=Poisson_DLikelihood
    , class my_PropD//=LM_MultivariateGaussian
    , class AP//=Landa
    ,template<
      class
      ,  class
      , class
      , class
      , class
      >
    class propDistStep=LevenbergMarquardt_step
    >
class Metropolis_Hastings_mcmc
{
public:
  typedef my_PropD pDist;

    struct test
    {
      test()=default;

      template<class E>
      static std::ostream& put(std::ostream& os,const mcmc_step<E,pDist>& sLik
                               ,const mcmc_step<E,pDist>& cLik)
      {
        if (cLik.isValid)
          {
            double logPcurrent=sLik.mlogbPL();
            double logPcandidate=cLik.mlogbPL();


            double logQforward=sLik.proposed.logP(cLik.param);
            double logQbackward=cLik.proposed.logP(sLik.param);
            double logForward=logPcandidate-logQforward;
            double logBackward=logPcurrent-logQbackward;

            os<<logForward<<" "<<logBackward<<" ";
          }
        return os;
      }

      template<class E>
      static bool accept(const mcmc_step<E,pDist>& sLik
                         ,const mcmc_step<E,pDist>& cLik
                         ,std::mt19937_64 &mt,
                         double& dHd,
                         double& logPcurrent,
                         double& logPcandidate,
                         double& logChi2forward,
                         double& logChi2backward,
                         double& logDetCurrent,
                         double& logDetCandidate)
      {
        logPcurrent=sLik.logbPL(mt);
        logPcandidate=cLik.logbPL(mt);

        logChi2forward=sLik.proposed.chi2(cLik.param);
        logChi2backward=cLik.proposed.chi2(sLik.param);
        logDetCandidate=cLik.proposed.logDetCov();
        logDetCurrent=sLik.proposed.logDetCov();
        auto logQforward=sLik.proposed.logP(cLik.param);
        auto logQbackward=cLik.proposed.logP(sLik.param);
        auto d=cLik.param-sLik.param;
        auto H=sLik.D_lik.H*sLik.beta+sLik.D_prior.H;
        dHd=0.5*fullSum(quadraticForm_B_A_BT(H,d));

        if (!cLik.isValid)
          {
            return false;
          }
        else
          {

            if (!std::isfinite(logQforward)|| !std::isfinite(logQbackward)
                ||!std::isfinite(logPcandidate))
              {
                return false;
              }
            else
              {
                double logA=logPcandidate-logQforward-(logPcurrent-logQbackward);
                double A=std::min(1.0,exp(logA));
                std::uniform_real_distribution<double> u(0,1);

                double r=u(mt);
                bool accept_=r<A;
                if (accept_)
                  {
                    return true;
                  }
                else
                  {
                    return false;
                  }
              }
          }
      }



      template<class E>
      static bool accept_swap(const mcmc_step<E,pDist>& sLik
                              ,const mcmc_step<E,pDist>& cLik
                              ,std::mt19937_64 &mt, double& s_dHd, double& c_dHd)
      {
        if (!cLik.isValid)
          {
            return false;
          }
        else
          {
            double logPcurrent=sLik.mlogbPL();
            double logPcandidate=cLik.mlogbPL();

            double logPSwapcurrent=sLik.logbPLb(cLik.beta,mt);
            double logPSwapcandidate=cLik.logbPLb(sLik.beta,mt);


            double logA=logPSwapcandidate+logPSwapcurrent-(logPcurrent+logPcandidate);
            logA=(cLik.beta-sLik.beta)*(sLik.logL(mt)-cLik.logL(mt));

            double A=std::min(1.0,exp(logA));
            std::uniform_real_distribution<double> u(0,1);
            auto d=cLik.param-sLik.param;
            auto s_H=sLik.D_lik.H*sLik.beta+sLik.D_prior.H+cLik.D_lik.H*sLik.beta+cLik.D_prior.H;
            auto c_H=sLik.D_lik.H*cLik.beta+sLik.D_prior.H+cLik.D_lik.H*cLik.beta+cLik.D_prior.H;
            s_dHd=0.25*fullSum(quadraticForm_B_A_BT(s_H,d));
            c_dHd=0.25*fullSum(quadraticForm_B_A_BT(c_H,d));

            double r=u(mt);
            bool accept_=r<A;
            if (accept_)
              {
                return true;
              }
            else
              {
                return false;
              }
          }
      }



    };


  template<class E>
  static void tempered_step
  (std::size_t nskip
   ,const propDistStep<D,M,D_Lik,my_PropD,AP>& LM_Lik
   ,const D_Lik& lik
   ,const M& model
   ,const D& data
   ,std::vector<mcmc_step<E,pDist>>& sDist
   ,std::vector<Adaptive_parameterized>& pars
   ,Master_Adaptive_Beta_New& aBeta
   ,std::vector<double>& dHd
   ,std::vector<double>& logPcandidate
   ,std::vector<double>& logPcurrent
   ,std::vector<double>& logChiforward
   ,std::vector<double>& logChibackward
   ,std::vector<double>& logDetCurrent
   ,std::vector<double>& logDetCandidate
   , double p_Tjump
   ,std::vector<std::mt19937_64>& mt
   ,std::size_t isamples
   , bool does_stdout
   ,double slogL_max,
   std::size_t ndts_max
   , std::ostream& os
   , const std::chrono::steady_clock::time_point& startTime
   , double& timeOpt)
  {
    std::vector<bool> out(sDist.size());
    aBeta.push_step();
    std::size_t isubSamples=0;

    double p_Tjump_cycle_desired=0.5;


    double pp_Tjump=p_Tjump;
    std::size_t nsubSamples=1;
    if (p_Tjump<p_Tjump_cycle_desired)
      {
        nsubSamples=p_Tjump_cycle_desired/p_Tjump;
        pp_Tjump=p_Tjump*nsubSamples;
      }
    while(isubSamples<nskip)
      {
        nsubSamples=std::min(nsubSamples,nskip-isubSamples);

        std::vector<std::stringstream> ss(sDist.size()+1);

#pragma omp parallel for
        for(std::size_t i=0; i<sDist.size(); ++i)
          {
            std::size_t isubSample_i=isubSamples;

            for(std::size_t ii=0; ii<nsubSamples; ++ii)
              {
                AP landa;
                mcmc_step<E,pDist> cDist;
                landa=pars[i].sample(mt[i]);
                sDist[i]=LM_Lik.update_mcmc_step
                    (lik,model,data,sDist[i],landa,aBeta.getBeta().getValue()[i]);
                M_Matrix<E> c=sDist[i].proposed.sample(mt[i]);
                //	sDist[i].proposed.autoTest(mt[i],500);
                cDist=LM_Lik.get_mcmc_step
                    (lik,model,data,c,landa,aBeta.getBeta().getValue()[i],sDist[i].iscout,slogL_max,ndts_max);
                //		std::cerr<<"logL "<<cDist.logLikelihood;
                //		std::cerr<<"("<<cDist.slogLik()<<","<<cDist.dts_size()<<") ";
                //		if (!cDist.isValid)
                //		    {
                //			std::cerr<<"rejection :"<<landa<<"\n";
                //			pars[i].push_rejection(landa);
                //		    }


                if (test::accept(sDist[i],cDist,mt[i],dHd[i],
                                 logPcandidate[i],logPcurrent[i],logChiforward[i],
                                 logChibackward[i],logDetCurrent[i],logDetCandidate[i]))
                  {
                    // aBeta.new_acceptance(i,sDist[i],cDist);
                    pars[i].push_acceptance(landa);
                    sDist[i]=cDist;
                    out[i] =true;
                  }
                else
                  {
                    //  aBeta.new_rjection(i,sDist[i]);
                    out[i]=false;
                    pars[i].push_rejection(landa);
                  }
                if (i==0)
                  {
                    auto tnow=std::chrono::steady_clock::now();
                    auto d=tnow-startTime;
                    double t0=timeOpt;
                    timeOpt=1.0e-6*std::chrono::duration_cast<std::chrono::microseconds>(d).count()/60.0;
                    auto timeIter_=60*(timeOpt-t0);
                    double evidence=Tempered_Evidence_Evaluation<mcmc_step<E,pDist>>::Evidence(aBeta.getBeta(),sDist);
                    ss[i]<<"isample::"<<isamples<<"\t"<<"isubSample::"<<isubSample_i<<"\t"<<timeOpt<<"\t"<<timeIter_<<"\t"<<"Evidence\t"<<evidence<<"\n";
                    isubSample_i++;

                  }
                put(i,ss[i+1],sDist,out,dHd,logPcandidate,logPcurrent,logChiforward, logChibackward,logDetCurrent,logDetCandidate);

              }



          }


        for (std::size_t i=0; i<ss.size(); ++i)
          os<<ss[i].str()<<"\n";

        if (does_stdout)
          for (std::size_t i=0; i<ss.size(); ++i)
            std::cout<<ss[i].str()<<"\n";

#pragma omp parallel for
        for(std::size_t i=1; i<sDist.size(); i+=2)
          {
            std::uniform_real_distribution<> u;
            double r=u(mt[i]);
            double s_dHd, c_dHd;
            if (r<pp_Tjump)
              {
                if (test::accept_swap(sDist[i-1],
                                      sDist[i],
                                      mt[i], s_dHd, c_dHd))
                  {
                    std::swap(sDist[i-1],sDist[i]);
                    std::swap(sDist[i-1].beta,sDist[i].beta);
                    //   os<<"swap::\t"<<sDist[i-1].beta<<"\t"<<sDist[i].beta<<"\n";
                    //   if (does_stdout)
                    //     std::cout<<"swap::\t"<<sDist[i-1].beta<<"\t"<<sDist[i].beta<<"\n";
                    //   aBeta.push_acceptance(sDist[i-1].beta,sDist[i].beta);
                  }
                // else
                //   aBeta.push_rejection(sDist[i-1].beta,sDist[i].beta);
                ;
              }
          }

#pragma omp parallel for

        for(std::size_t i=2; i<sDist.size(); i+=2)
              {
            std::uniform_real_distribution<> u;
            double r=u(mt[i]);
            double s_dHd, c_dHd;
            if (r<pp_Tjump)
              {
                if (test::accept_swap(sDist[i-1],
                                      sDist[i],
                                      mt[i], s_dHd, c_dHd))
                  {
                    std::swap(sDist[i-1],sDist[i]);
                    std::swap(sDist[i-1].beta,sDist[i].beta);
                    //   os<<"swap::\t"<<sDist[i-1].beta<<"\t"<<sDist[i].beta<<"\n";
                    //   if (does_stdout)
                    //     std::cout<<"swap::\t"<<sDist[i-1].beta<<"\t"<<sDist[i].beta<<"\n";
                    //   aBeta.push_acceptance(sDist[i-1].beta,sDist[i].beta);
                  }
                // else
                //   aBeta.push_rejection(sDist[i-1].beta,sDist[i].beta);
                ;
              }
          }
        isubSamples+=nsubSamples;
      }

  }



  template<class E>
  static
  void save_state(std::ostream& os,
                  const std::vector<mcmc_step<E,pDist>>& sDist
                  ,const std::vector<Adaptive_parameterized>& pars
                  ,const Master_Adaptive_Beta_New& aBeta
                  ,std::size_t isamples
                  , std::size_t i_sim,
                  const std::vector<std::mt19937_64>& mts )
  {
    os<<"sDist\n";
    for (std::size_t i=0; i<sDist.size(); ++i)
      {
        os<<sDist[i].iscout<<"\t"<<sDist[i].beta<<"\n";
        os<<sDist[i].param<<"\n";
      }
    os<<"pars\n";
    os<<pars<<"\n";
    os<<"aBeta\n";
    os<<aBeta<<"\n";
    os<<"isamples\n";
    os<< isamples<<"\n";
    os<<"i_sim\n";
    os<< i_sim<<"\n";
    os<<"mts\n";
    os<< mts<<"\n";
  };

  template<class E>
  static
  bool load_state(std::istream& is,
                  std::vector<mcmc_step<E,pDist>>& sDist
                  ,std::vector<Adaptive_parameterized>& pars
                  ,Master_Adaptive_Beta_New& aBeta
                  ,std::size_t& isamples
                  ,std::size_t& i_sim
                  ,std::vector<std::mt19937_64>& mts)
  {
    std::string line;
    std::getline(is,line);
    for (std::size_t i=0; i<sDist.size(); ++i)
      {
        if (! (is>>sDist[i].iscout>>sDist[i].beta))
          return false;
        std::getline(is,line);

        if (!(is>>sDist[i].param))
          return false;
        std::getline(is,line);
      }
    std::getline(is,line);
    if (!(is>>pars))
      return false;
    std::getline(is,line);
    std::getline(is,line);
    if (!(is>>aBeta))
      return false;
    std::getline(is,line);
    std::getline(is,line);
    if (!(is>> isamples))
      return false;
    std::getline(is,line);
    std::getline(is,line);
    if (!(is>> i_sim))
      return false;
    std::getline(is,line);
    std::getline(is,line);
    if (!(is>> mts))
      return false;
    std::getline(is,line);
    return true;
  };



};



template<
    class Ad,
    class D
    ,   class M
    , class D_Lik//=Poisson_DLikelihood
    , class my_PropD//=LM_MultivariateGaussian
    , class AP//=trust_region
    ,template<
      class
      ,  class
      ,  class
      , class
      , class
      >
    class propDist//=LevenbergMarquardt_step
    ,template<
      class
      ,class
      ,  class
      , class
      ,class
      ,class
      ,template<
        class
        , class
        ,  class
        , class
        , class
        >
      class
      >
    class MH=Metropolis_Hastings_mcmc>

class Template_Tempering_mcmc
{
public:
  template<class E>
  using mystep= mcmc_step<E,my_PropD>;

  template<class E>
  static void
  run
  (const MH<Ad,D,M,D_Lik,my_PropD,AP,propDist>& mcmc
   ,const propDist<D,M,D_Lik,my_PropD,AP>& LMLik
   ,const D_Lik& lik
   ,const M& model
   ,const D& data
   , const Ad& landa_Dist0
   ,const Master_Adaptive_Beta_New& beta0
   , double maxTime
   , std::size_t nsamples
   ,std::size_t nskip
   ,std::size_t nAdapt
   ,double pTjump
   ,double slogL_max
   ,std::size_t ndts_max
   ,std::mt19937_64::result_type seed
   ,std::string EviName
   ,std::string EviNameLog0
   ,const std::chrono::steady_clock::time_point& startTime
   ,double& timeOpt
   ,std::size_t maxSimFileSize
   ,bool does_stdout
   , const std::string& state_file)
  {
    std::mt19937_64 mt;
    mt.seed(seed);
    std::size_t i_sim=0;

    std::string f_par_name0=EviName+"_par.txt";
    std::string f_logL_name0=EviName+"_logL.txt";
    std::string f_fit_name0=EviName+"_fit.txt";
    std::string f_sim_name0=EviName+"_sim.txt";

    std::string f_state_name=EviName+"_state.txt";

    std::string f_par_name=f_par_name0    +"."+leadingZeroZero(i_sim);
    std::string f_logL_name=f_logL_name0    +"."+leadingZeroZero(i_sim);
    std::string f_fit_name=f_fit_name0    +"."+leadingZeroZero(i_sim);
    std::string f_sim_name=f_sim_name0    +"."+leadingZeroZero(i_sim);
    std::string f_log_name=EviNameLog0    +"."+leadingZeroZero(i_sim);

    std::ofstream os;
    std::ofstream f_par;
    std::ofstream f_logL;
    std::ofstream f_sim;
    std::ofstream f_fit;






    Master_Adaptive_Beta_New beta(beta0);
    auto n=beta.size();
    std::vector<mystep<E>> sDists(n);

    std::vector<Ad> pars(n,landa_Dist0);
    std::vector<double> dHd(beta.size());
    std::vector<double> logPcandidate(beta.size());
    std::vector<double> logPcurrent(beta.size());
    std::vector<double> logChiforward(beta.size());
    std::vector<double> logChibackward(beta.size());
    std::vector<double> logDetCurrent(beta.size());
    std::vector<double> logDetCandidate(beta.size());

    std::uniform_int_distribution<typename std::mt19937_64::result_type> useed;
    std::vector<std::mt19937_64> mts(n);
    std::size_t o=0;

    bool isContinuation =!state_file.empty();
    if (isContinuation)
      {
        std::ifstream f_state;
        f_state.open(state_file.c_str(), std::ifstream::in );
        mcmc.load_state(f_state,sDists,pars,beta,o,i_sim,mts);

        f_state.close();
        f_par_name=f_par_name0    +"."+leadingZeroZero(i_sim);
        f_logL_name=f_logL_name0    +"."+leadingZeroZero(i_sim);
        f_fit_name=f_fit_name0    +"."+leadingZeroZero(i_sim);
        f_sim_name=f_sim_name0    +"."+leadingZeroZero(i_sim);
        f_log_name=EviNameLog0    +"."+leadingZeroZero(i_sim);



        os.open(f_log_name.c_str(), std::ofstream::out | std::ofstream::trunc);
        f_par.open(f_par_name.c_str(), std::ofstream::out | std::ofstream::trunc);
        f_logL.open(f_logL_name.c_str(), std::ofstream::out | std::ofstream::trunc);
        f_sim.open(f_sim_name.c_str(), std::ofstream::out | std::ofstream::trunc);
        f_fit.open(f_fit_name.c_str(), std::ofstream::out | std::ofstream::trunc);
        f_par<<std::setprecision(std::numeric_limits<double>::digits10 + 1);

        if (o>nsamples/10)
          {
            for (std::size_t i=0; i<n;++i)
              {
                std::cerr<<"\n"<<beta.getBeta().getValue()[i]<<"\t";
                pars[i].actualize();
              }
          }
        else
          {
            for (std::size_t i=0; i<n;++i)
              {
                std::cerr<<"\n"<<beta.getBeta().getValue()[i]<<"\t";
                pars[i].actualize(nAdapt*nskip);
              }
          }
        for (std::size_t i=0; i<beta.size(); ++i)
          {
            double beta0=beta.getBeta().getValue()[i];
            AP landa;
            mcmc_post<E> postL;
            postL=lik.get_mcmc_Post(model,data,sDists[i].param,slogL_max,ndts_max);
            // landa=pars[i].sample(mts[i]);
            sDists[i]=LMLik.get_mcmc_step(lik,model,data,sDists[i].param,postL,landa,beta0,sDists[i].iscout);
          }

      }

    else{



        os.open(f_log_name.c_str(), std::ofstream::out | std::ofstream::app);
        f_par.open(f_par_name.c_str(), std::ofstream::out | std::ofstream::app);
        f_logL.open(f_logL_name.c_str(), std::ofstream::out | std::ofstream::app);
        f_sim.open(f_sim_name.c_str(), std::ofstream::out | std::ofstream::app);
        f_fit.open(f_fit_name.c_str(), std::ofstream::out | std::ofstream::app);
        f_par<<std::setprecision(std::numeric_limits<double>::digits10 + 1);

        std::stringstream ss;

        ss<<"model\t";
        ss<<"seed\t";
        ss<<"time\t";
        ss<<"nsteps\t";
        ss<<"nsample\t";
        ss<<"Evidence\t";

        auto s=mystep<E>{};
        s.writelogLHeaderDataFrame(ss);


        f_logL<<ss.str()<<"\tbetaEvidence"<<std::endl;
        f_par<<ss.str()<<"\t";
        f_sim<<ss.str()<<"\t";
        f_fit<<ss.str()<<"\t";


        s.writeParamHeaderDataFrame(f_par,model);

        //auto param=model.getPrior();
        //param.writeHeaderDataFrame(f_par);
        f_par<<std::endl;


        for (std::size_t i=0; i<n; ++i)
          mts[i].seed(useed(mt));

        for (std::size_t i=0; i<beta.size(); ++i)
          {
            M_Matrix<E> pinit;
            double beta0=beta.getBeta().getValue()[i];
            std::size_t ntrialsj=0;
            pars[i].actualize();

            bool isvalid=false;
            AP landa;
            while(!isvalid)
              {
                std::size_t ntrialsi=0;
                mcmc_post<E> postL;
                while(!postL.isValid)
                  {
                    pinit=lik.sample(model,data,mts[i]);
                    postL=lik.get_mcmc_Post(model,data,pinit,slogL_max,ndts_max);
                    ++ntrialsi;
                  }
                landa=pars[i].sample(mts[i]);
                sDists[i]=LMLik.get_mcmc_step(lik,model,data,pinit,postL,landa,beta0,i);
                isvalid=sDists[i].isValid;
                ++ntrialsj;
              }
            // LMLik.update_mcmc_step(lik,model,data,sDists[i],landa,beta0);
          }
        sDists[0].writeSimulationHeaderDataFrame(f_sim,data,model);
        f_sim<<std::endl;
        sDists[0].writeFitHeaderDataFrame(f_fit,data,model);
        f_fit<<std::endl;

      }



    while (o<nsamples&&timeOpt<maxTime*60)
      {
        mcmc.tempered_step
            (nskip,LMLik,lik,model,data,sDists,pars,beta,dHd,
             logPcandidate,logPcurrent,logChiforward,logChibackward,logDetCurrent,logDetCandidate,pTjump,mts,o,does_stdout,slogL_max,ndts_max,os,startTime,timeOpt);


        o++;


        auto tnow=std::chrono::steady_clock::now();
        auto d=tnow-startTime;
        double t0=1.0e-6*std::chrono::duration_cast<std::chrono::microseconds>(d).count()/60.0;
        double evidence=Tempered_Evidence_Evaluation<mystep<E>>::Evidence(beta.getBeta(),sDists);
        auto evidences=Tempered_Evidence_Evaluation<mystep<E>>::
                                                              Evidence_of_beta(beta.getBeta(),sDists);

        std::stringstream ss1;
        ss1<<model.id()<<"\t";
        ss1<<seed<<"\t";
        ss1<<t0<<"\t";
        ss1<<o*sDists.size()<<"\t";
        ss1<<o<<"\t";
        ss1<<evidence<<"\t";

//        std::vector<std::stringstream> s_logL(sDists.size());
//        std::vector<std::stringstream> s_par(sDists.size());
//        std::vector<std::stringstream> s_fit(sDists.size());
//        std::vector<std::stringstream> s_sim(sDists.size());

//#pragma omp parallel for
        for (std::size_t i=0; i<sDists.size(); ++i)
          {

            (sDists[i].writelogLRowDataFrame(f_logL,ss1.str()))<<"\t"<<evidences[i]<<"\n";
            (sDists[i].writeParamRowDataFrame(f_par,model,ss1.str()))<<"\n";
            (sDists[i].writeYfitRowDataFrame(f_fit,ss1.str()))<<"\n";
          }


        if (o%2==0)
          {
//#pragma omp parallel for
            for (std::size_t i=0; i<sDists.size(); ++i)
              {
                (sDists[i].writeSimulationRowDataFrame
                    (f_sim,data,model,ss1.str()))<<"\n";
              }
          }


//        for (std::size_t i=0; i<sDists.size(); ++i)
//          {
//            f_logL<<s_logL[i].str();
//            f_par<<s_par[i].str();
//            f_fit<<s_fit[i].str();
//            f_sim<<s_sim[i].str();

//          }

        f_logL.flush();
         f_par.flush();
         f_fit.flush();

          f_sim.flush();
          os.flush();
        std::size_t f_sim_pos=
            f_sim.tellp()+f_logL.tellp()+f_fit.tellp()+f_par.tellp()+os.tellp();
        if (f_sim_pos>maxSimFileSize)
          {
            f_sim.close();
            f_fit.close();
            f_logL.close();
            f_par.close();
            os.close();

            std::cerr<<f_sim_name<<"is completed !!\n";
            std::cerr<<f_fit_name<<"is completed !!\n";
            std::cerr<<f_logL_name<<"is completed !!\n";
            std::cerr<<f_log_name<<"is completed !!\n";
            std::cerr<<f_par_name<<"is completed !!\n";

            ++i_sim;

            std::ofstream f_state;
            std::string f_tmp=f_state_name+".tmp";
            f_state.open(f_tmp.c_str(), std::ofstream::out | std::ofstream::trunc);
            f_state<<std::setprecision(std::numeric_limits<double>::digits10 + 1);
            mcmc.save_state(f_state,sDists,pars,beta,o,i_sim,mts);
            f_state.close();
            std::remove(f_state_name.c_str());


            rename_done(f_sim_name);
            rename_done(f_logL_name);
            rename_done(f_log_name);
            rename_done(f_par_name);
            rename_done(f_fit_name);
            std::rename(f_tmp.c_str(),f_state_name.c_str());

            f_sim_name=f_sim_name0    +"."+leadingZeroZero(i_sim);
            f_sim.open(f_sim_name.c_str(), std::ofstream::out | std::ofstream::app);

            f_log_name=EviNameLog0+"."+leadingZeroZero(i_sim);
            os.open(f_log_name.c_str(), std::ofstream::out | std::ofstream::app);


            f_par_name=f_par_name0    +"."+leadingZeroZero(i_sim);
            f_par.open(f_par_name.c_str(), std::ofstream::out | std::ofstream::app);

            f_fit_name=f_fit_name0    +"."+leadingZeroZero(i_sim);;
            f_fit.open(f_fit_name.c_str(), std::ofstream::out | std::ofstream::app);

            f_logL_name=f_logL_name0    +"."+leadingZeroZero(i_sim);;
            f_logL.open(f_logL_name.c_str(), std::ofstream::out | std::ofstream::app);

            std::cerr<<f_sim_name<<"is opened !!\n";
            std::cerr<<f_par_name<<"is opened !!\n";
            std::cerr<<f_logL_name<<"is opened !!\n";
            std::cerr<<f_fit_name<<"is opened !!\n";
            std::cerr<<f_log_name<<"is opened !!\n";
          }

        if (o>nsamples/10)
          {
#pragma omp parallel for
            for (std::size_t i=0; i<n;++i)
              {
                //  std::cerr<<"\n"<<beta.getBeta().getValue()[i]<<"\t";
                pars[i].actualize();
              }
          }
        else
          {
#pragma omp parallel for
            for (std::size_t i=0; i<n;++i)
              {
                //  std::cerr<<"\n"<<beta.getBeta().getValue()[i]<<"\t";
                pars[i].actualize(nAdapt*nskip);
              }
          }


      }




  }







};








} // namespace opt


#endif // MYOPTIMIZATION_H
