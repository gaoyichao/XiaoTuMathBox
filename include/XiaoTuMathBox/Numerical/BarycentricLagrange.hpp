/************************************************************************
 * 
 * 重心拉格朗日多项式(Barycentric Lagrange Interpolation)
 * 
 * http://www.lmpt.univ-tours.fr/~nicolis/Licence_NEW/08-09/barycentric.pdf
 * https://gaoyichao.com/Xiaotu/?book=数值计算&title=重心拉格朗日插值
 * 
 ***********************************************************************/


#ifndef XTMB_NUMERICAL_BARYCENTRICLAGRANGE_H
#define XTMB_NUMERICAL_BARYCENTRICLAGRANGE_H


#include <vector>
#include <functional>
#include <XiaoTuMathBox/Common/Common.hpp>

namespace xiaotu {

    template <typename Scalar>
    class BarycentricLagrange
    {
        public:
            /**
             * @brief 构造函数
             * 
             * @param [in] x 采样点 x 列表
             * @param [in] y 对应 x 列表的采样值
             */
            BarycentricLagrange(std::vector<Scalar> const & x, std::vector<Scalar> const & y) 
                : mXs(x), mYs(y)
            {
                assert(mXs.size() == mYs.size());
                InitWeights();
            }

            /**
             * @brief 计算 lagrange 插值
             */
            Scalar Evaluate(Scalar x) const
            {
                Scalar numerator = 0;
                Scalar denominator = 0;

                for (int i = 0; i < mXs.size(); ++i) {
                    if (std::abs(x - mXs[i]) < SMALL_VALUE)
                        return mYs[i];

                    Scalar term = mWs[i] / (x - mXs[i]);
                    numerator += term * mYs[i];
                    denominator += term;
                }
                return numerator / denominator;
            }

            /**
             * @brief 计算 lagrange 插值
             */
            Scalar operator()(Scalar const & x) const { return Evaluate(x); }

            /**
             * @brief 添加新采样点
             * 
             * @param [in] x 新采样点 x 坐标
             * @param [in] y 新采样点 y 坐标
             */
            BarycentricLagrange & AddPoint(Scalar x, Scalar y)
            {
                size_t n = mXs.size();
                // 更新现有的权重
                for (int i = 0; i < n; ++i) {
                    assert(std::abs(mXs[i] - x) > SMALL_VALUE);
                    mWs[i] /= (mXs[i] - x);
                }

                // 计算新采样点权重
                Scalar w = 1;
                for (int i = 0; i < n; ++i)
                    w /= (x - mXs[i]);
                
                mXs.push_back(x);
                mYs.push_back(y);
                mWs.push_back(w);

                NormalizeWeights();
                return *this;
            }

        private:

            /**
             * @brief 初始化权重
             */
            void InitWeights()
            {
                size_t n = mXs.size();
                mWs.resize(n);
                
                // w_j = 1 / \prod_{k \neq i} (x_j - x_k)
                for (size_t i = 0; i < n; ++i) {
                    Scalar w = 1;
                    for (int k = 0; k < n; ++k) {
                        if (k == i)
                            continue;
                        w *= (mXs[i] - mXs[k]);
                    }
                    mWs[i] = 1 / w;
                }
                NormalizeWeights();
            }


            void NormalizeWeights()
            {
                size_t n = mXs.size();
                // 为了数值稳定性, 规避绝对值恨大的权重
                size_t idx = IdxOfMaxAbs(mWs);
                Scalar max = mWs[idx];
                for (size_t i = 0; i < n; ++i) {
                    mWs[i] /= max;
                }
            }

        private:
            std::vector<Scalar> mXs;
            std::vector<Scalar> mYs;
            std::vector<Scalar> mWs;
    };


}


#endif
