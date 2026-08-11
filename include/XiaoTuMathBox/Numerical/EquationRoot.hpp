#ifndef XTMB_NUMERICAL_EQUATION_ROOT
#define XTMB_NUMERICAL_EQUATION_ROOT

#include <functional>
#include <XiaoTuMathBox/Common/Common.hpp>

namespace xiaotu {

    //! @brief 二分法求解方程 f(x) = 0 在闭区间 [a, b] 上的根
    //!
    //! @param [in] f 目标函数
    //! @param [in] a 闭区间起点
    //! @param [in] b 闭区间终点
    //! @param [in] max_iter 最大迭代次数
    //! @param [in] tol 终止迭代时的区间长度
    template <typename DataType>
    DataType Bisection(std::function<DataType(DataType)> f, DataType a, DataType b,
                       int max_iter = 100, DataType tol = SMALL_VALUE)
    {
        DataType A = f(a);
        DataType B = f(b);

        if (0 == A)
            return a;
        if (0 == B)
            return b;

        assert(a < b);
        assert(Sign(A) * Sign(B) < 0);

        DataType re = b;
        for (int i = 0; i < max_iter; ++i) {
            DataType half_length = 0.5 * std::abs(b - a);
            re = a + half_length;
            B = f(re);
            if (0 == B || half_length < tol)
                return re;

            if (Sign(A) * Sign(B) > 0) {
                A = B;
                a = re;
            } else {
                b = re;
            }
        }

        return re;
    }


    //! @brief 牛顿法求解方程 f(x) = 0
    //!
    //! @param [in] f 目标函数
    //! @param [in] df f(x) 的导函数
    //! @param [in] x0 迭代初值
    //! @param [in] max_iter 最大迭代次数
    //! @param [in] tol 终止迭代时的区间长度
    template <typename DataType>
    DataType NewtonRaphson(std::function<DataType(DataType)> f,
                           std::function<DataType(DataType)> df, 
                           DataType x0,
                           int max_iter = 100, DataType tol = SMALL_VALUE)
    {
        DataType re = x0;
        for (int i = 0; i < max_iter; ++i) {
            DataType y = f(re);
            if (0 == y)
                return re;

            DataType dydx = df(re);
            assert(0 != dydx); // 是否需要 std::abs(dfdx) < SMALL_VALUE ?

            re = x0 - y / dydx;
            if (std::abs(re - x0) < tol)
                return re;

            x0 = re;
        }

        return re;
    }


    //! @brief 割线法求解方程 f(x) = 0
    //!
    //! @param [in] f 目标函数
    //! @param [in] x0 迭代初值
    //! @param [in] x1 迭代初值
    //! @param [in] max_iter 最大迭代次数
    //! @param [in] tol 终止迭代时的区间长度
    template <typename DataType>
    DataType SecantRoot(std::function<DataType(DataType)> f, DataType x0, DataType x1,
                        int max_iter = 100, DataType tol = SMALL_VALUE)
    {
        DataType y0 = f(x0);
        DataType y1 = f(x1);

        for (int i = 2; i < max_iter; ++i) {
            DataType dx = x1 - x0;
            DataType dy = y1 - y0;
            assert(0 != dy);

            DataType x = x1 - y1 * dx / dy;
            if (std::abs(x - x1) < tol)
                return x;

            x0 = x1;
            y0 = y1;
            x1 = x;
            y1 = f(x);
        }

        return x1;
    }



}

#endif
