#ifndef XTMB_NUMERICAL_EQUATION_ROOT
#define XTMB_NUMERICAL_EQUATION_ROOT

#include <functional>
#include <complex>
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

        if (0 == y0)
            return x0;
        if (0 == y1)
            return x1;

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

    //! @brief 试位法求解方程 f(x) = 0
    //!
    //! @param [in] f 目标函数
    //! @param [in] x0 迭代初值
    //! @param [in] x1 迭代初值
    //! @param [in] max_iter 最大迭代次数
    //! @param [in] tol 终止迭代时的区间长度
    template <typename DataType>
    DataType FalsePosition(std::function<DataType(DataType)> f, DataType x0, DataType x1,
                           int max_iter = 100, DataType tol = SMALL_VALUE)
    {
        DataType y0 = f(x0);
        DataType y1 = f(x1);

        if (0 == y0)
            return x0;
        if (0 == y1)
            return x1;

        assert(Sign(y0) * Sign(y1) < 0);

        for (int i = 2; i < max_iter; ++i) {
            DataType dx = x1 - x0;
            DataType dy = y1 - y0;
            assert(0 != dy);

            DataType x = x1 - y1 * dx / dy;
            if (std::abs(x - x1) < tol)
                return x;

            DataType y = f(x);
            if(Sign(y) * Sign(y1) < 0) {
                x0 = x1;
                y0 = y1;
            }
            x1 = x;
            y1 = y;
        }

        return x1;
    }

    //! @brief Dekker 求解方程 f(x) = 0
    //!
    //! @param [in] f 目标函数
    //! @param [in] x0 迭代初值
    //! @param [in] x1 迭代初值
    //! @param [in] max_iter 最大迭代次数
    //! @param [in] tol 终止迭代时的区间长度
    template <typename DataType>
    DataType DekkerRoot(std::function<DataType(DataType)> f, DataType x0, DataType x1,
                           int max_iter = 100, DataType tol = SMALL_VALUE)
    {
        DataType y0 = f(x0);
        DataType y1 = f(x1);

        if (0 == y0)
            return x0;
        if (0 == y1)
            return x1;

        assert(Sign(y0) * Sign(y1) < 0);

        for (int i = 2; i < max_iter; ++i) {
            // 保证 |y1| <= |y0|, 即 x1 是当前更佳近似
            if (std::abs(y0) < std::abs(y1)) {
                std::swap(x0, x1);
                std::swap(y0, y1);
            }

            DataType dx = x1 - x0;
            DataType dy = y1 - y0;

            DataType m = 0.5 * (x0 + x1);
            DataType s = (0 == dy) ? m : (x1 - y1 * dx / dy);
            DataType x = InRange<DataType>(s, m, x1) ? s : m;

            if (std::abs(x - x1) < tol)
                return x;

            DataType y = f(x);
            if (0 == y)
                return x;

            if(Sign(y) * Sign(y1) < 0) {
                x0 = x;
                y0 = y;
            } else {
                x1 = x;
                y1 = y;
            }
        }
        return x1;
    }


    //! @brief Brent 求解方程 f(x) = 0
    //!
    //! https://en.wikipedia.org/wiki/Brent's_method
    //!
    //! @param [in] f 目标函数
    //! @param [in] a 迭代初值
    //! @param [in] b 迭代初值
    //! @param [in] max_iter 最大迭代次数
    //! @param [in] tol 终止迭代时的区间长度
    template <typename DataType>
    DataType BrentRoot(std::function<DataType(DataType)> f, DataType a, DataType b,
                           int max_iter = 100, DataType tol = SMALL_VALUE)
    {
        DataType fa = f(a);
        DataType fb = f(b);

        if (0 == fa)
            return a;
        if (0 == fb)
            return fb;
        assert(Sign(fa) * Sign(fb) < 0);

        // 保证 |y1| <= |y0|, 即 x1 是当前更佳近似
        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }

        // b_{k-1}
        DataType c = a;
        DataType fc = fa;
        // b_{k-2}
        DataType d = a;
        DataType fd = fa;

        // 标记上次迭代是否采用二分法
        bool mflag = true;
        // 下一轮候选点
        DataType s = a;

        for (int i = 2; i < max_iter; ++i) {
            DataType m = 0.5 * (a + b);
            if (fa != fc && fb != fc) {
                s = a * fb * fc / ((fa - fb) * (fa - fc))
                  + b * fa * fc / ((fb - fa) * (fb - fc))
                  + c * fa * fb / ((fc - fa) * (fc - fb));
            } else if (fb != fa) {
                s = b - fb * (b - a) / (fb - fa);
            } else {
                s = m;
            }

            if (!InRange<DataType>(s, 0.75*a + 0.25 * b, b) ||
                ( mflag && (std::abs(s - b) >= 0.5 * std::abs(b - c))) ||
                (!mflag && (std::abs(s - b) >= 0.5 * std::abs(c-d))) ||
                ( mflag && (std::abs(b - c) < tol)) ||
                (!mflag && (std::abs(c - d) < tol))) {
                s = m;
                mflag = true;
            } else {
                mflag = false;
            }
            

            if (std::abs(s - b) < tol)
                return s;
            DataType fs = f(s);
            if (0 == fs)
                return s;

            d = c;  // 第一轮时 d 并未赋值, mflag 为 true, d 还未参与计算。
            fd = fc;
            c = b;
            fc = fb;

            if(Sign(fa) * Sign(fs) < 0) {
                b = s;
                fb = fs;
            } else {
                a = s;
                fa = fs;
            }

            if (std::abs(fa) < std::abs(fb)) {
                std::swap(a, b);
                std::swap(fa, fb);
            }
        }

        return b;
    }

    //! @brief 一元二次多项式方程的根 
    //!
    //! @param [in] a,b,c 方程的系数 \(ax^2 + bx + c = 0\)
    //! @param [out] x0 复数形式的根
    //! @param [out] x1 复数形式的根
    //! @return 是否为两个实根
    template <typename DataType>
    bool QuadraticRoot(DataType a, DataType b, DataType c,
                       std::complex<DataType> & x0, std::complex<DataType> & x1)
    {
        DataType v = b * b - 4 * a * c;

        if (v >= 0) {
            v = std::sqrt(v);
            // 两个实根
            if (b < 0) {
                x0 = (-b + v) / a * 0.5;
                x1 = -2 * c / (b - v);
            } else {
                x0 = (-b - v) / a * 0.5;
                x1 = -2 * c / (b + v);
            }
            return true;
        } else {
            DataType real = -0.5 * b / a;
            DataType img = std::sqrt(-v) / a * 0.5;

            x0.real(real);
            x0.imag(img);

            x1.real(real);
            x1.imag(-img);
            return false;
        }

    }

}

#endif
