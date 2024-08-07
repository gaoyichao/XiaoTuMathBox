#ifndef XTMB_LA_EQUATION_ELIMINATION_H
#define XTMB_LA_EQUATION_ELIMINATION_H

#include <stdint.h>
#include <cassert>
#include <vector>
#include <iostream>

#include <XiaoTuMathBox/LinearAlgibra/MatrixView.hpp>

namespace xiaotu {
namespace math {

    //! @brief Gauss-Jordan 消元发求解线性方程组 A x = b
    //!
    //! @param [inout] A 线性方程组的系数矩阵, 输出 A 的逆
    //! @param [inout] b 线性方程组右侧的向量们, 输出解们 x
    template <typename MatViewA, typename MatViewB>
    void GaussJordanEliminate(MatViewA & A, MatViewB * b)
    {
        typedef typename MatViewA::Scalar Scalar;

        assert(A.Rows() == A.Cols());
        assert(nullptr == b || A.Rows() == b->Rows());

        int n = A.Rows();
        int m = (nullptr == b) ? 0 : b->Cols();
        // 用于记录发现 pivoting 的行列坐标
        std::vector<int> idxc(n), idxr(n);
        // 用于记录已经消元的行列
        std::vector<int> pivflag(n, 0);

        for (int t = 0; t < n; t++) {
            int pivRow = -1;
            int pivCol = -1;
            Scalar max = 0;
            // 找出未消元的子阵中最大值作为 pivoting
            for (int r = 0; r < n; r++) {
                if (1 == pivflag[r])
                    continue;
                for (int c = 0; c < n; c++) {
                    if (1 == pivflag[c])
                        continue;
                    Scalar abs = std::abs(A(r, c));
                    if (abs > max) {
                        max = abs;
                        pivRow = r;
                        pivCol = c;
                    }
                }
            }
            pivflag[pivCol]++;
            idxc[t] = pivCol;
            idxr[t] = pivRow;

            if (0.0 == A(pivCol, pivRow))
                throw "奇异矩阵";

            // 交换 pivoting 到对角线上
            if (pivRow != pivCol) {
                A.RowSwap(pivRow, pivCol);
                if (nullptr != b)
                    b->RowSwap(pivRow, pivCol);
            }
            Scalar pinv = 1.0 / max;
            A(pivCol, pivCol) = 1;

            // 归一化 pivoting 所在行
            for (int c = 0; c < n; c++)
                A(pivCol, c) *= pinv;
            if (nullptr != b) {
                for (int c = 0; c < m; c++)
                    b->At(pivCol, c) *= pinv;
            }

            // 对其余各行消元
            for (int r = 0; r < n; r++) {
                if (r == pivCol)
                    continue;
                Scalar dum = A(r, pivCol);
                A(r, pivCol) = 0.0;

                for (int c = 0; c < n; c++)
                    A(r, c) -= A(pivCol, c) * dum;
                if (nullptr != b) {
                    for (int c = 0; c < m; c++)
                        (*b)(r, c) -= (*b)(pivCol, c) * dum;
                }
            }
        }

        // 恢复逆矩阵
        for (int idx = n - 1; idx >= 0; idx--) {
            if (idxc[idx] == idxr[idx])
                continue;
            A.ColSwap(idxc[idx], idxr[idx]);
        }
    }

    //! @brief Gauss-Jordan 消元
    //!
    //! @param [inout] A 方阵, 输出 A 的逆
    template <typename MatViewA>
    inline void GaussJordanEliminate(MatViewA & A)
    {
        GaussJordanEliminate<MatViewA, MatViewA>(A, nullptr);
    }

}
}

#endif
