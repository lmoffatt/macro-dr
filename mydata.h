#ifndef MYDATA_H
#define MYDATA_H

#include "mytypetraits.h"

template<class ...> class DataArray;


template<class variable,class...coordinates>
class DataArray<variable,Cs<coordinates...>>
{


private:


};



#endif // MYDATA_H
