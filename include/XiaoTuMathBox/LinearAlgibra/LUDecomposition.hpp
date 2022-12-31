#ifndef XTMB_LA_LU_DECOMPOSITION_H
#define XTMB_LA_LU_DECOMPOSITION_H

#include <cassert>
#include <vector>

#include <XiaoTuMathBox/LinearAlgibra/MatrixView.hpp>

namespace xiaotu {
namespace math {

    //! @brief LU 分解, 给定一个 n x n 的矩阵，根据行主元进行重新排序，
    //! 同时对重排的矩阵进行 LU 分解。结果记录在成员变量 mLU 中。
    //!
    //! | a_00 a_01 a_02 a_03 |   | \alpha_00                               | | \beta_00 \beta_01 \beta_02 \beta_03 |
    //! | a_10 a_11 a_12 a_13 | = | \alpha_10 \alpha_11                     | |          \beta_11 \beta_12 \beta_13 |
    //! | a_20 a_21 a_22 a_23 |   | \alpha_20 \alpha_21 \alpha_22           | |                   \beta_22 \beta_23 |
    //! | a_30 a_31 a_32 a_33 |   | \alpha_30 \alpha_31 \alpha_32 \alpha_33 | |                            \beta_33 |
    //!
    //! 其中 \alpha 阵的对角线被约束为 1，所以 mLU 可以把两个矩阵合到一起:
    //!
    //! |  \beta_00  \beta_01  \beta_02 \beta_03 |
    //! | \alpha_10  \beta_11  \beta_12 \beta_13 |
    //! | \alpha_20 \alpha_21  \beta_22 \beta_23 |
    //! | \alpha_30 \alpha_31 \alpha_32 \beta_33 |
    //!
    //! 
    template <typename T>
    class LU_Decompose {
        public:
            typedef T  DataType;
            typedef T* DataTypePtr;

        public:
            LU_Decompose(MatrixView<T> & A, MatrixView<T> & lu)
                : mLU(lu), mSwapIndex(A.Rows())
            {
                assert(A.Rows() == A.Cols());
                int n = mLU.Rows();

                // 记录各行主元的数组的倒数
                std::vector<T> vlot_inv(n);
                for (int ridx = 0; ridx < n; ridx++) {
                    T max_abs = 0;
                    for (int cidx = 0; cidx < n; cidx++) {
                        T abs = std::abs(mLU(ridx, cidx));
                        if (abs > max_abs)
                            max_abs = abs;
                    }
                    if (0 == max_abs)
                        throw "奇异矩阵";
                    vlot_inv[ridx] = 1.0 / max_abs;
                }

                // 遍历构造LU
                mSwapTimes = 0;
                for (int count = 0; count < n; count++) {
                    // 查找剩余子阵中主元的最大者
                    T max_abs = 0;
                    int rmax = count;
                    for (int ridx = count; ridx < n; ridx++) {
                        T scaled = vlot_inv[ridx] * std::abs(mLU(ridx, count));
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
                    if (0 == mLU(count, count))
                        throw "奇异矩阵";
                    // 更新剩余子阵
                    for (int ridx = count + 1; ridx < n; ridx++) {
                        T alpha = mLU(ridx, count) / mLU(count, count);
                        mLU(ridx, count) = alpha;
                        for (int cidx = count + 1; cidx < n; cidx++) {
                            T beta = mLU(ridx, cidx) - alpha * mLU(count, cidx);
                            mLU(count, cidx) = beta;
                        }
                    }
                }
            }

        public:
            MatrixView<T> const & LU() { return mLU; }

        private:
            MatrixView<T> mLU;
            //! 记录下交换行索引
            std::vector<int> mSwapIndex;
            int mSwapTimes;
    };

}
}

#endif
