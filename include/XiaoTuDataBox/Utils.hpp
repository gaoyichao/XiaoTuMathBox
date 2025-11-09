#ifndef XTMB_UTILS_H
#define XTMB_UTILS_H

#include <stdint.h>
#include <cassert>
#include <iostream>

#define XTLog(ost) (ost) << __FILE__ << ":" << __LINE__ << ": "

inline bool XOR(bool x, bool y)
{
    return x ^ y;
}

#endif
