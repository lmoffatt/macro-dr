#ifndef MYMATH_H
#define MYMATH_H

#include <cmath>


std::size_t base2_floor(std::size_t x)
{
    unsigned n=0;
    while (x>>=1) ++n;
    return 1<<n;
}

inline constexpr double PI = 3.14159265358979323846;

inline constexpr double mynan = std::numeric_limits<double>::quiet_NaN();


#endif // MYMATH_H
