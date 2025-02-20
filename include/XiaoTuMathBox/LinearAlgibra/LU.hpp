#ifndef XTMB_LA_LU_H
#define XTMB_LA_LU_H

#include <cassert>
#include <vector>

namespace xiaotu::math {

    //! @brief LU 分解, 给定一个 n x n 的矩阵，根据行主元进行重新排序，
    //! 同时对重排的矩阵进行 LU 分解。结果记录在成员变量 mLU 中。
    //!
    //! | \alpha_00                               | | \beta_00 \beta_01 \beta_02 \beta_03 |   | a_00 a_01 a_02 a_03 |
    //! | \alpha_10 \alpha_11                     | |          \beta_11 \beta_12 \beta_13 | = | a_10 a_11 a_12 a_13 |
    //! | \alpha_20 \alpha_21 \alpha_22           | |                   \beta_22 \beta_23 |   | a_20 a_21 a_22 a_23 |
    //! | \alpha_30 \alpha_31 \alpha_32 \alpha_33 | |                            \beta_33 |   | a_30 a_31 a_32 a_33 |
    //!
    //! 其中 \alpha 阵的对角线被约束为 1，所以 mLU 可以把两个矩阵合到一起:
    //!
    //! |  \beta_00  \beta_01  \beta_02 \beta_03 |
    //! | \alpha_10  \beta_11  \beta_12 \beta_13 |
    //! | \alpha_20 \alpha_21  \beta_22 \beta_23 |
    //! | \alpha_30 \alpha_31 \alpha_32 \beta_33 |
    template <typename MatrixViewA>
    class LU {
        public:
            typedef typename MatrixViewA::Scalar Scalar;

            LU(MatrixViewA const & a)
                : mBuffer(a.NumDatas()),
                  mLU(mBuffer.data()),
                  mSwapTimes(0),
                  mSwapIndex(a.Rows())
            {
                assert(a.Rows() == a.Cols());
                mLU.Assign(a);
                Decompose();
            }

            //! @brief 求解方程组 Ax=b, b 和 x 可以是同一个对象
            //! 
            //! @param [in] b 方程右侧的列向量
            //! @param [out] x 对应 b 中每一列的解
            template <typename MatrixViewB>
            void Solve(MatrixViewB const & b, MatrixViewB & x)
            {
                assert(b.Rows() == mLU.Rows());
            
                const int N = mLU.Rows();
                const int M = b.Cols();
                x.Assign(b);

                // 遍历每一列
                for (int cidx = 0; cidx < M; ++cidx) {
                    for (int ridx = 0; ridx < N; ++ridx) {
                        int r = mSwapIndex[ridx];
                        Scalar y = x(r, cidx);
                        x(r, cidx) = x(ridx, cidx);

                        for (int k = 0; k < ridx; k++)
                            y -= mLU(ridx, k) * x(k, cidx);

                        x(ridx, cidx) = y;
                    }

                    for (int ridx = N-1; ridx >= 0; --ridx) {
                        Scalar s = x(ridx, cidx);
                        for (int k = ridx+1; k < N; ++k)
                            s -= mLU(ridx, k) * x(k, cidx);
                        x(ridx, cidx) = s / mLU(ridx, ridx);
                    }
                }
            }

            void Inverse(MatrixViewA inv)
            {
                inv.Identity();
                Solve(inv, inv);
            }


        private:
            void Decompose()
            {
                const int N = mLU.Rows();
                // 记录各行主元的倒数
                std::vector<Scalar> vlot_inv(N);
                for (int ridx = 0; ridx < N; ridx++) {
                    auto max_abs = 0;
                    for (int cidx = 0; cidx < N; cidx++) {
                        auto abs = std::abs(mLU(ridx, cidx));
                        if (abs > max_abs)
                            max_abs = abs;
                    }
                    if (max_abs < SMALL_VALUE)
                        throw std::runtime_error("奇异矩阵");
                    vlot_inv[ridx] = 1.0 / max_abs;
                }

                // 遍历构造 LU
                mSwapTimes = 0;
                for (int count = 0; count < N; count++) {
                    // 查找当前列的最大值作为主元
                    auto max_abs = 0;
                    int rmax = count;
                    for (int ridx = count; ridx < N; ridx++) {
                        auto scaled = vlot_inv[ridx] * std::abs(mLU(ridx, count));
                        if (scaled > max_abs) {
                            max_abs = scaled;
                            rmax = ridx;
                        }
                    }
                    // 若当前行不是主元则交换
                    if (rmax != count) {
                        mLU.RowSwap(rmax, count);
                        mSwapTimes++;
                        std::swap(vlot_inv[rmax], vlot_inv[count]);
                    }
                    mSwapIndex[count] = rmax;
                    if (std::abs(mLU(count, count)) < SMALL_VALUE)
                        throw std::runtime_error("奇异矩阵");
                    // 更新剩余子阵
                    for (int ridx = count + 1; ridx < N; ridx++) {
                        auto alpha = mLU(ridx, count) / mLU(count, count);
                        mLU(ridx, count) = alpha;
                        for (int cidx = count + 1; cidx < N; cidx++) {
                            auto beta = mLU(ridx, cidx) - alpha * mLU(count, cidx);
                            mLU(ridx, cidx) = beta;
                        }
                    }
                }
            }

        public:
            //! @brief 获取 LU 分解矩阵
            MatrixViewA & operator() () { return mLU; }
            MatrixViewA const & operator() () const { return mLU; }

        private:
            std::vector<Scalar> mBuffer;
            MatrixViewA mLU;
            //! 记录下交换行索引
            int mSwapTimes;
            std::vector<int> mSwapIndex;
    };

}


#endif
