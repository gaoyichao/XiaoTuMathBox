#ifndef XTMB_COMMON_H
#define XTMB_COMMON_H

#define SMALL_VALUE 1e-12

namespace xiaotu {

    /**
     * @brief 获取数值的符号
     */
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

    /**
     * @brief 计算一个数组中绝对值最大的元素索引
     */
    template <typename VectorLike>
    size_t IdxOfMaxAbs(VectorLike const & list)
    {
        size_t re = 0;
        auto max_abs = std::abs(list[0]);
        for (size_t i = 1; i < list.size(); i++) {
            auto abs = std::abs(list[i]);
            if (abs > max_abs) {
                re = i;
                max_abs = abs;
            }
        }
        return re;
    }
    
    //! @brief 各种类型萃取器的声明, 需要自行提供特化类
    template<typename T>
    struct Traits;

}


#endif
