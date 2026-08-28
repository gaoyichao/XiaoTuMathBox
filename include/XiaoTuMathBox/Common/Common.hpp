#ifndef XTMB_COMMON_H
#define XTMB_COMMON_H

#define SMALL_VALUE 1e-12

namespace xiaotu {

    template <typename T>
    int Sign(T val) {
        return (T(0) < val) - (val < T(0));
    }

    /**
     * @brief 判定 c 是否在 a,b 之间
     */
    template <typename T>
    bool InRange(T c, T a, T b)
    {
        if (a < b)
            return (a < c) && (c < b);
        else
            return (b < c) && (c < a);
    }
    
    //! @brief 各种类型萃取器的声明, 需要自行提供特化类
    template<typename T>
    struct Traits;

}


#endif
