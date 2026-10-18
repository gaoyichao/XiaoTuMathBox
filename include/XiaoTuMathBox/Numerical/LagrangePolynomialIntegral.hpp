/************************************************************************
 * 
 * 基于拉格朗日多项式的数值积分
 * 
 * https://gaoyichao.com/Xiaotu/?book=数值计算&title=Lagrange多项式数值积分
 * 
 ***********************************************************************/
#ifndef XTMB_NUMERICAL_LAGRANGE_POLYNOMIAL_INTEGRAL_H
#define XTMB_NUMERICAL_LAGRANGE_POLYNOMIAL_INTEGRAL_H


namespace xiaotu {

    template <typename Scalar>
    class LagrangePolynomialIntegral
    {
        public:
            /**
             * @brief 构造函数
             * 
             * @param [in] a 积分区间 [a,b]
             * @param [in] b
             * @param [in] n 采样点数量
             * @param [in] closed 是否闭式采样, 区间端点 a,b 是否为采样点
             */
            LagrangePolynomialIntegral(Scalar a, Scalar b, size_t n, bool closed = true)
            {
                mXs.resize(n);

                if (closed) {
                    double step = (b - a) / (n-1);
                    for (int i = 0; i < n; i++)
                        mXs[i] = a + i * step;
                } else {
                    double step = (b - a) / (n-3);
                    for (int i = 0; i < n; i++)
                        mXs[i] = a + step + i * step;
                }
                InitTable(a, b);
            }

            /**
             * @brief 构造函数
             * 
             * @param [in] a 积分区间 [a,b]
             * @param [in] b
             * @param [in] x 采样点
             */
            LagrangePolynomialIntegral(Scalar a, Scalar b, std::vector<Scalar> const & x)
                : mXs(x)
            {
                assert(x.size() > 1);
                InitTable(a, b);
            }

            /**
             * @brief 数值积分
             * 
             * @param [in] y 各采样点函数值
             */
            Scalar Integral(std::vector<Scalar> const & y)
            {
                assert(y.size() == mCs.size());
                Scalar re = 0;
                for (size_t i = 0; i < y.size(); i++) {
                    re += mCs[i] * y[i];
                }
                return re;
            }

            /**
             * @brief 数值积分
             * 
             * @param [in] f 被积函数
             */
            Scalar Integral(std::function<Scalar(Scalar)> f)
            {
                Scalar re = 0;
                for (size_t i = 0; i < mCs.size(); i++) {
                    re += mCs[i] * f(mXs[i]);
                }
                return re;
            }

        private:

            /**
             * @brief 初始化积分系数表
             */
            void InitTable(Scalar a, Scalar b)
            {
                mCs.resize(mXs.size());
                for (size_t i = 0; i < mXs.size(); i++) {
                    auto p = LagrangeBasis(mXs, i).Integral();
                    mCs[i] = p(b) - p(a);
                }
            }

        private:
            // 采样点
            std::vector<Scalar> mXs;
            //! 积分系数
            std::vector<Scalar> mCs;
    };
}

#endif