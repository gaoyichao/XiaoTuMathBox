#ifndef XTMB_NUMERICAL_EQUATION_ROOT
#define XTMB_NUMERICAL_EQUATION_ROOT

#include <functional>
#include <XiaoTuMathBox/Common/Common.hpp>

namespace xiaotu {

    template <typename T>
    int Sign(T val) {
        return (T(0) < val) - (val < T(0));
    }

    //! @brief 二分法求解方程 f(x) = 0 在闭区间 [a, b] 上的根
    //!
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


}

#endif
