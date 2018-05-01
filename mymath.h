#ifndef MYMATH_H
#define MYMATH_H

#include <cmath>


std::size_t base2_floor(std::size_t x)
{
    unsigned n=0;
    while (x>>=1) ++n;
    return 1<<n;
}




#endif // MYMATH_H
