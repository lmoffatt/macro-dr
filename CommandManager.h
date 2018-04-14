#ifndef COMMANDMANAGER
#define COMMANDMANAGER
#include <string>
#include <map>
#include "myCommandManagement.h"
#include"myInputSerializer.h"

#include <experimental/type_traits>
/*! @file CommandManager.h   Management of the commands
 *
 *
 * */









inline std::string& removeComments(std::string& line)
{
  auto pos=line.find("//");
  if (pos!=line.npos)
    line.erase(pos);
  return line;
}


inline std::string& replaceLabel(std::string& line,
                                 const std::string& label,
                                 const std::string& replacement)
{
  std::size_t i=0;
  while (i!=line.npos)
    {
      i=line.find(label,i);
      if (i!=line.npos)
        line.replace(i,label.size(),replacement);
    }
  return line;
}


inline std::string& replaceLabel(std::string& line,
                                 const std::vector<std::string>& label,
                                 const std::vector<std::string>& replacement)
{
  for (std::size_t j=0; j<label.size(); ++j)
    {
      std::size_t i=0;
      while (i!=line.npos)
        {
          i=line.find(label[j],i);
          if (i!=line.npos)
            line.replace(i,label[j].size(),replacement[j]);
        }
    }
  return line;
}





struct has_ClassName_tag{};

struct has_not_ClassName_tag{};


template<class, typename T=void>
struct has_ClassName_traits
{
  typedef has_not_ClassName_tag tag;
};

template<class T>
struct has_ClassName_traits<T,decltype (T::ClassName())>
{
  typedef has_ClassName_tag tag;

};



template<typename>
std::string ClassName_imp(has_not_ClassName_tag)
{
  return "unknown";
}



template<typename T>
auto ClassName_imp(has_ClassName_tag)->decltype (T::ClassName())
{
  return T::ClassName();
}





inline
bool read_from_stream(std::istream&,std::ostream&,... )
{
  return true;
}




template<typename T>
auto read_from_stream(std::istream& is,std::ostream& /*logstream*/,T& x)
->decltype (bool(is>>x))
{

  if (is>>x)
    return true;
  else
    return false;
}

template<typename T>
auto read_from_stream(std::istream& is,std::ostream& /*logstream*/,T*& x)
->decltype (bool(is>>*x))
{

  return bool(is>>*x);
}


template<typename T>
auto read_from_stream(std::istream& is,std::ostream& logstream,T*& x)
->decltype (x->read(std::string(),is,logstream))
{
  std::string s;
  return x->read(s,is,logstream);
}

template<typename T>
auto read_from_stream(std::istream& is,std::ostream& logstream,T& x)
->decltype (x.read(std::string(),is,logstream))
{
  std::string s;
  return x.read(s,is,logstream);
}


inline
bool write_to_stream(std::ostream&,std::ostream&,... )
{
  return true;
}




template<typename T>
auto write_to_stream(std::ostream& os,std::ostream& /*logstream*/,T& x)
->decltype (bool(os<<x))
{

  if (os<<x) return true;
  else return false;
}

template<typename T>
auto write_to_stream(std::ostream& os,std::ostream& /*logstream*/,T*& x)
->decltype (bool(os<<*x))
{

  return bool(os<<*x);
}


template<typename T>
auto write_to_stream(std::ostream& os,std::ostream& logstream,T*& x)
->decltype (x->write(std::string(),os,logstream))
{
  std::string s;
  return x->write(s,os,logstream);
}

template<typename T>
auto write_to_stream(std::ostream& os,std::ostream& logstream,T& x)
->decltype (x.write(std::string(),os,logstream))
{
  std::string s;
  return x.write(s,os,logstream);
}







template<class T>
class Cls
{
public:
  static std::string name()
  {
    return ClassName_imp<T>(typename has_ClassName_traits<T>::tag());
  }
  static
  bool read(std::istream& is,T& x,std::ostream& logstream)
  {
    return read_from_stream(is,logstream,x);
  }

  static
  bool write (std::ostream& os, const T& x, std::ostream& logstream)
  {
    return write_to_stream(os,logstream,x);

  }


};



template<class T>
class Cls<T*>
{
public:
  static std::string name()
  {
    return ClassName_imp<T>(typename has_ClassName_traits<T>::tag());
  }
  static
  bool read(std::istream& is,T*& x,std::ostream& logstream)
  {
    return read_from_stream(is,logstream,x);
  }
  static
  bool write (std::ostream& os, const T* x, std::ostream& logstream)
  {
    return write_to_stream(os,logstream,x);

  }

};





template<class Command_Manager>
class myScript
{
public:
  myScript(Command_Manager* cm):cm_(cm){}

  int run(char* filename, std::ostream &logs)
  {
    std::ifstream f(filename);
    while (f)
      {
        std::string line;
        safeGetline(f,line);
        removeComments(line);
        std::cerr<<line;
        cm_->execute(line,logs);
      }
    return 0;

  }


  int runDefine(const std::string& filename,
                const std::vector<std::__cxx11::string> &label,
                const std::vector<std::__cxx11::string> &valueInplace,
                std::ostream &logs)
  {
    std::ifstream f(filename);
    while (f)
      {
        std::string line;
        safeGetline(f,line);
        removeComments(line);

        replaceLabel(line,label,valueInplace);
        cm_->execute(line,logs);
      }
    return 0;

  }


private:
  Command_Manager* cm_;


};









template<typename T>
class WriteIt
{
public:
  static
  void apply(const T& x,
             const std::string& idname,
             std::ostream* f,
             std::ostream* logs)
  {
    *f<<idname<<"\n"<<Cls<T>::name()<<"\n";
    Cls<T>::write(*f,x,*logs);
    *f<<"\n";
  }

};


template<typename,class=void>
struct doesWriteDataFrame: std::false_type{};

template< class ... > using void_t = void;


template<typename T>
struct doesWriteDataFrame<T, void_t<decltype(std::declval<T>().writeDataFrame(std::declval<std::ostream&>())) >> : std::true_type{};



template <bool, typename T>
struct writeDataFrame
{
  static void write(const T&,std::ostream*, std::ostream* log)
  {
    *log<<Cls<T>::name()<<" writeDataFrame is not implemented\n";

  }

};

template <typename T>
struct writeDataFrame<true,T>
{
  static void write(const T& x,std::ostream* f, std::ostream*/* log*/)
  {
    x.writeDataFrame(*f);
    *f<<"\n";
  }

};

template <typename T>
struct writeDataFrame<true,T*>
{
  static void write(const T* x,std::ostream* f, std::ostream* /*log*/)
  {
    x->writeDataFrame(*f);
    *f<<"\n";

  }

};





template <typename T>
class DataFrameIt
{
public:
  static void apply(const T& x,std::ostream* f,std::ostream* log)
  {
    writeDataFrame<doesWriteDataFrame<std::remove_pointer_t<T>>::value,T>::write(x,f,log);
  }

};








template<class Cm>
struct read
{
  void operator()(Cm* cm,const std::string& filename, std::ostream* logs  )const
  {

  }

};

template<class Cm>
struct write{
  bool operator()(Cm* cm,const std::string& idname, std::ostream* logs,  std::string fname, bool append) const


  {
    if (fname.empty())
      fname=idname+"_out.txt";
    std::ofstream f;
    if (append)
      f.open(fname, std::ofstream::app | std::ofstream::out);
    else
      f.open(fname,  std::ofstream::out);

    if (!f)
      {
        *logs<<"could not open file "+fname+" for writing";
        return false;
      }
    else
      {
        f<<idname<<"\n";
        cm->apply(Co<WriteIt>(),idname,idname,&f,logs);
        if (f)
          {
            *logs<<idname<<" written in file "<<fname;
            f.close();
            return true;
          }
        else
          {
            *logs<<idname<<" something wrong when written in file "<<fname;
            f.close();
            return false;
          }

      }
  }

};

template<class Cm>
struct dataFrame
{
  bool operator()(Cm* cm,const std::string& idname, std::ostream* logs,  std::string fname) const
  {
    if (fname.empty())
      fname=idname+"_data_frame.txt";
    std::ofstream f;
    f.open(fname,  std::ofstream::out);

    if (!f)
      {
        *logs<<"could not open file "+fname+" for writing";
        return false;
      }
    else
      {
        cm->apply(Co<DataFrameIt>(),idname,&f,logs);
        if (f)
          {
            *logs<<idname<<" written in file "<<fname;
            f.close();
            return true;
          }
        else
          {
            *logs<<idname<<" something wrong when written in file "<<fname;
            f.close();
            return false;
          }

      }
  }
};




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






#endif // COMMANDMANAGER

