/************************************************************************
 * 
 * 三次样条插值
 * 
 * https://gaoyichao.com/Xiaotu/?book=数值计算&title=三次样条插值
 * 
 ***********************************************************************/
#ifndef XTMB_NUMERICAL_CUBIC_SPLINE_H
#define XTMB_NUMERICAL_CUBIC_SPLINE_H

#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

namespace xiaotu {


    template <typename Scalar>
    class CubicSpline
    {
        public:
             /**
             * @brief 自然三次样条构造函数(Nature Cubic Spline)
             * 
             * @param [in] x 采样点 x 列表，严格升序排列
             * @param [in] y 对应 x 列表的采样值
             */
            CubicSpline(
                    std::vector<Scalar> const & x,
                    std::vector<Scalar> const & y)
                : mXs(x), mAs(y)
            {
                //! 三次样条插值至少需要 3 个点
                assert(mXs.size() >= 3);
                assert(mXs.size() == mAs.size());
                //!  严格升序排列
                for (size_t i = 1; i < mXs.size(); ++i)
                    assert(mXs[i] > mXs[i-1]);

                
                size_t n_points = mXs.size();
                size_t n = n_points - 1;
                mBs.resize(n+1);
                mCs.resize(n+1);
                mDs.resize(n+1);

                std::vector<Scalar> h(n);
                for (size_t i = 0; i < n; ++i)
                    h[i] = mXs[i+1] - mXs[i];

                InitNatureSpline(h);
            }

            /**
             * @brief 固定三次样条构造函数(Clamped Cubic Spline)
             * 
             * @param [in] x 采样点 x 列表，严格升序排列
             * @param [in] y 对应 x 列表的采样值
             * @param [in] df_x0 起点的一阶导数 f'(x_0)
             * @param [in] df_xn 重点的一阶导数 f'(x_n)
             */
            CubicSpline(
                    std::vector<Scalar> const & x,
                    std::vector<Scalar> const & y,
                    Scalar df_x0, Scalar df_xn)
                : mXs(x), mAs(y)
            {
                //! 三次样条插值至少需要 3 个点
                assert(mXs.size() >= 3);
                assert(mXs.size() == mAs.size());
                //!  严格升序排列
                for (size_t i = 1; i < mXs.size(); ++i)
                    assert(mXs[i] > mXs[i-1]);

                
                size_t n_points = mXs.size();
                size_t n = n_points - 1;
                mBs.resize(n+1);
                mCs.resize(n+1);
                mDs.resize(n+1);

                std::vector<Scalar> h(n);
                for (size_t i = 0; i < n; ++i)
                    h[i] = mXs[i+1] - mXs[i];

                InitClampedSpline(h, df_x0, df_xn);
            }

            /**
             * @brief 计算给定 x 处的三次样条插值
             */
            Scalar Evaluate(Scalar x) const
            {
                if (x <= mXs.front())
                    return mAs.front();
                if (x >= mXs.back())
                    return mAs.back();

                // 定位 x 所属区间 [x_i, x_{i+1}]
                int i = 1;
                for (; i < mXs.size(); ++i) {
                    if (x == mXs[i])
                        return mAs[i];
                    if (x < mXs[i])
                        break;
                }
                i--;

                Scalar dx = x - mXs[i];
                Scalar dx2 = dx * dx;
                Scalar dx3 = dx2 * dx;

                return mAs[i] + mBs[i] * dx + mCs[i] * dx * dx + mDs[i] * dx * dx * dx;
            }

            /**
             * @brief 计算给定 x 处的三次样条插值
             */
            Scalar operator()(Scalar const & x) const { return Evaluate(x); }

        private:

            /**
             * @brief 按照自然三次样条构造
             * 
             * @param [in] h 各个子区间长度
             */
            void InitNatureSpline(std::vector<Scalar> const & h)
            {
                size_t n = h.size();
                //! 构造矩阵 H
                DMatrix<Scalar> H(n+1, n+1);
                H(0,0) = 1;
                for (size_t i = 1; i < n; ++i) {
                    H(i, i-1) = h[i-1];
                    H(i, i) = 2 * (h[i-1] + h[i]);
                    H(i, i+1) = h[i];
                }
                H(n,n) = 1;
                //! 构造矩阵 beta
                DMatrix<Scalar> beta(n+1, 1);
                beta(0) = 0;
                for (size_t i = 1; i < n; ++i) {
                    beta(i) = 3 *((mAs[i+1]-mAs[i])/h[i] - (mAs[i]-mAs[i-1])/h[i-1]);
                }
                beta(n) = 0;

                //! LU 分解求 c_i
                LU lu(H.View());
                auto c = DMatrixView<Scalar>(mCs.data(), n+1, 1);
                lu.Solve(beta, c);

                //! 计算 mBs, mDs
                for (size_t i = 0; i < n; ++i) {
                    mDs[i] = (mCs[i+1] - mCs[i]) / (3 * h[i]);
                    mBs[i] = (mAs[i+1] - mAs[i]) / h[i] - (mCs[i+1] + 2 * mCs[i]) * h[i] / 3;
                }
            }

            /**
             * @brief 按照固定三次样条构造
             * 
             * @param [in] h 各个子区间长度
             * @param [in] df_x0 起点的一阶导数 f'(x_0)
             * @param [in] df_xn 重点的一阶导数 f'(x_n)
             */
            void InitClampedSpline(std::vector<Scalar> const & h,
                    Scalar df_x0, Scalar df_xn)
            {
                size_t n = h.size();
                //! 构造矩阵 H
                DMatrix<Scalar> H(n+1, n+1);
                H(0,0) = 2 * h[0]; H(0,1) = h[0];
                for (size_t i = 1; i < n; ++i) {
                    H(i, i-1) = h[i-1];
                    H(i, i) = 2 * (h[i-1] + h[i]);
                    H(i, i+1) = h[i];
                }
                H(n,n-1) = h[n-1]; H(n,n) = 2 * h[n-1];
                //! 构造矩阵 beta
                DMatrix<Scalar> beta(n+1, 1);
                beta(0) = 3 *((mAs[1] - mAs[0]) / h[0] - df_x0);
                for (size_t i = 1; i < n; ++i) {
                    beta(i) = 3 *((mAs[i+1]-mAs[i])/h[i] - (mAs[i]-mAs[i-1])/h[i-1]);
                }
                beta(n) = 3 *(df_xn - (mAs[n] - mAs[n-1]) / h[n-1]);

                //! LU 分解求 c_i
                LU lu(H.View());
                auto c = DMatrixView<Scalar>(mCs.data(), n+1, 1);
                lu.Solve(beta, c);

                //! 计算 mBs, mDs
                for (size_t i = 0; i < n; ++i) {
                    mDs[i] = (mCs[i+1] - mCs[i]) / (3 * h[i]);
                    mBs[i] = (mAs[i+1] - mAs[i]) / h[i] - (mCs[i+1] + 2 * mCs[i]) * h[i] / 3;
                }
            }

        private:
            std::vector<Scalar> mXs; // 节点的 x 坐标
            std::vector<Scalar> mAs; // 常数项系数 (y_i)
            std::vector<Scalar> mBs; // 一次项系数
            std::vector<Scalar> mCs; // 二次项系数
            std::vector<Scalar> mDs; // 三次项系数
    };

}


#endif
