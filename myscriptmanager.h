#ifndef MYSCRIPTMANAGER_H
#define MYSCRIPTMANAGER_H

#include "mySerializer.h"
#include <fstream>
#include <iostream>
#include <vector>

inline std::string& removeComments(std::string& line)
{
  auto pos=line.find("//");
  if (pos!=line.npos)
    line.erase(pos);
  return line;
}


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
        std::cerr<<line<<std::endl;
        cm_->execute(line,logs);
      }
    return 0;

  }


  int runDefine(const std::string& filename,
                const std::vector<std::string> &/*label*/,
                const std::vector<std::string> &/*valueInplace*/,
                std::ostream &logs)
  {
    std::ifstream f(filename);
    while (f)
      {
        std::string line;
        safeGetline(f,line);
        removeComments(line);

   //     replaceLabel(line,label,valueInplace);
        cm_->execute(line,logs);
      }
    return 0;

  }


private:
  Command_Manager* cm_;


};



#endif // MYSCRIPTMANAGER_H
