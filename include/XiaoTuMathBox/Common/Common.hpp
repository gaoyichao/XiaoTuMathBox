#ifndef XTMB_COMMON_H
#define XTMB_COMMON_H

#define SMALL_VALUE 1e-12

namespace xiaotu {

    template <typename T>
    int Sign(T val) {
        return (T(0) < val) - (val < T(0));
    }
    
}


#endif
