/************************************************************************
 * 
 * Gaussian-Legendre 积分
 * 
 * https://gaoyichao.com/Xiaotu/?book=数值计算&title=高斯积分
 * 
 ***********************************************************************/
#ifndef XTMB_NUMERICAL_GAUSSIAN_INTEGRATION
#define XTMB_NUMERICAL_GAUSSIAN_INTEGRATION

#include <vector>
#include <functional>
#include <XiaoTuMathBox/Common/Common.hpp>


namespace xiaotu {

    template <typename Scalar>
    class GaussianIntegration
    {
        public:
            /**
             * @brief 构造函数
             * 
             * @param [in] n > 1 勒让德多项式阶数
             */
            GaussianIntegration(int n)
                : _n(n)
            {
                assert(_n > 1);

                _roots.resize(_n);
                _weights.resize(_n);
                int m = (_n + 1) / 2;

                auto P = LegendrePolynomial<Scalar>(_n);
                for (int i = 0; i < m; ++i) {
                    // 契比雪夫初值估计
                    Scalar z = std::cos(static_cast<Scalar>(M_PI) * (i + 0.75) / (_n + 0.5));
                    Scalar r = P.NewtonRaphson(z);

                    Scalar pxi, dpxi;
                    P.Horner(r, pxi, dpxi);

                    _roots[i] = -r;
                    _roots[_n - 1 - i] = r;

                    // 权重公式
                    Scalar w = 2.0 / ((1.0 - r * r) * dpxi * dpxi);
                    _weights[i] = w;
                    _weights[_n - 1 - i] = w;
                }
            }

            /**
             * @brief 对函数 f 在区间 [a,b] 上积分
             * 
             * @param [in] f 被积函数
             * @param [in] a 积分区间
             * @param [in] b
             */
            Scalar Evaluate(std::function<Scalar(Scalar)> && f, Scalar a, Scalar b) const
            {
                Scalar integral = 0;
                
                Scalar slope = (b - a) / 2.0;
                Scalar intercept = (b + a) / 2.0;
                for (int i = 0; i < _n; ++i) {
                    // 将标准区间 [-1, 1] 的节点映射到实际区间 [a, b]
                    Scalar x_mapped = slope * _roots[i] + intercept;
                    integral += _weights[i] * f(x_mapped);
                }

                return integral * slope;
            }

            /**
             * @brief 对函数 f 在区间 [a,b] 上积分
             * 
             * @param [in] f 被积函数
             * @param [in] a 积分区间
             * @param [in] b
             */
            Scalar operator()(std::function<Scalar(Scalar)> && f, Scalar a, Scalar b) const
            {
                return Evaluate(std::move(f), a, b);
            }

            friend std::ostream & operator << (std::ostream & s, GaussianIntegration const & g)
            {
                s << "--- " << g._n << " 点高斯-勒让德积分表 ---" << std::endl;
                for (int i = 0; i < g._n; ++i) {
                    s << "节点 x[" << i << "] = " << g._roots[i] 
                      << "\t 权重 w[" << i << "] = " << g._weights[i] << std::endl;
                }
                s << "-------------------------------------" << std::endl;
                return s;
            }

        private:
            //! 勒让德多项式阶数
            int _n;
            //! 勒让德多项式的零点
            std::vector<Scalar> _roots;
            //! 积分权重
            std::vector<Scalar> _weights;
    };

}


#endif

