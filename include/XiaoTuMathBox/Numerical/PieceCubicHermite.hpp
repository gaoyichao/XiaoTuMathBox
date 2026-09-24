/************************************************************************
 * 
 * 分段三次 Hermite 插值
 * 
 * https://gaoyichao.com/Xiaotu/?book=数值计算&title=分段三次多项式插值
 * 
 ***********************************************************************/
#ifndef XTMB_NUMERICAL_NAIVE_PIECE_CUBIC_HERMITE_H
#define XTMB_NUMERICAL_NAIVE_PIECE_CUBIC_HERMITE_H


namespace xiaotu {

    /**
     * @brief 估计采样点的导数
     * 
     * @param [in] x 采样点 x 列表
     * @param [in] y 对应 x 列表的采样值
     * @return 对应 x 列表的一阶导数值
     */
    template <typename Scalar>
    std::vector<Scalar> FritschCarlsonDerivatives(
            std::vector<Scalar> const & x, 
            std::vector<Scalar> const & y)
    {
        size_t n = x.size();
        std::vector<Scalar> d(n, 0);

        std::vector<Scalar> h(n - 1), s(n - 1);
        for (size_t i = 0; i < n - 1; ++i) {
            h[i] = x[i + 1] - x[i];
            s[i] = (y[i + 1] - y[i]) / h[i];
        }

        for (size_t i = 1; i < n - 1; ++i) {
            if (s[i - 1] * s[i] <= 0) {
                d[i] = 0;
            } else {
                Scalar w1 = 2.0 * h[i] + h[i - 1];
                Scalar w2 = h[i] + 2.0 * h[i - 1];
                d[i] = (w1 + w2) / (w1 / s[i - 1] + w2 / s[i]);
            }
        }

        if (n >= 2) {
            d[0] = ((2.0 * h[0] + h[1]) * s[0] - h[0] * s[1]) / (h[0] + h[1]);
            if (d[0] * s[0] <= 0) {
                d[0] = 0;
            } else if (s[0] * s[1] > 0 && std::abs(d[0]) > 3.0 * std::abs(s[0])) {
                d[0] = 3.0 * s[0];
            }

            size_t last = n - 1;
            d[last] = ((2.0 * h[last - 1] + h[last - 2]) * s[last - 1] - h[last - 1] * s[last - 2]) / (h[last - 1] + h[last - 2]);
            if (d[last] * s[last - 1] <= 0) {
                d[last] = 0;
            } else if (s[last - 1] * s[last - 2] > 0 && std::abs(d[last]) > 3.0 * std::abs(s[last - 1])) {
                d[last] = 3.0 * s[last - 1];
            }
        }
        return d;
    }


    template <typename Scalar>
    class PieceCubicHermite
    {
        public:
             /**
             * @brief 构造函数
             * 
             * @param [in] x 采样点 x 列表
             * @param [in] y 对应 x 列表的采样值
             * @param [in] dy 对应 x 列表的一阶导数值
             */
            PieceCubicHermite(
                    std::vector<Scalar> const & x,
                    std::vector<Scalar> const & y,
                    std::vector<Scalar> const & dy)
                : mXs(x), mYs(y), mDYs(dy)
            {
                assert(x.size() >= 2);
                assert(x.size() == y.size());
                assert(x.size() == dy.size());
        
                for (size_t i = 1; i < mXs.size(); ++i)
                    assert(mXs[i] > mXs[i-1]);
            }

            /**
             * @brief 计算给定 x 处的朴素分段三次 Hermite 插值
             */
            Scalar Evaluate(Scalar x) const
            {
                if (x <= mXs.front())
                    return mYs.front();
                if (x >= mXs.back())
                    return mYs.back();

                // 定位 x 所述区间 [x_i, x_{i+1}]
                int i = 1;
                for (; i < mXs.size(); ++i) {
                    if (x == mXs[i])
                        return mYs[i];
                    if (x < mXs[i])
                        break;
                }
                i--;

                Scalar x0 = mXs[i];
                Scalar x1 = mXs[i + 1];
                Scalar y0 = mYs[i];
                Scalar y1 = mYs[i + 1];
                Scalar d0 = mDYs[i];
                Scalar d1 = mDYs[i + 1];

                // 归一化映射：将局部区间映射到标准区间 t \in [0, 1]
                Scalar h = x1 - x0;
                Scalar t = (x - x0) / h;

                Scalar t2 = t * t;
                Scalar t3 = t2 * t;

                Scalar h00 =  2.0 * t3 - 3.0 * t2 + 1.0; // 控制左端点值
                Scalar h10 =        t3 - 2.0 * t2 + t;   // 控制左端点斜率
                Scalar h01 = -2.0 * t3 + 3.0 * t2;       // 控制右端点值
                Scalar h11 =        t3 -       t2;       // 控制右端点斜率

                return h00 * y0 + h10 * h * d0 + h01 * y1 + h11 * h * d1;
            }


            /**
             * @brief 计算给定 x 处的朴素分段三次 Hermite 插值
             */
            Scalar operator()(Scalar const & x) const { return Evaluate(x); }


        private:
            //! 采样点 x 坐标
            std::vector<Scalar> mXs;
            //! 采样值 y
            std::vector<Scalar> mYs;
            //! 采样点一阶导数值
            std::vector<Scalar> mDYs;
    };

}


#endif
