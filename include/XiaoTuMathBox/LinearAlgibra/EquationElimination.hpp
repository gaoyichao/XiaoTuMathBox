#ifndef XTMB_LA_EQUATION_ELIMINATION_H
#define XTMB_LA_EQUATION_ELIMINATION_H

#include <stdint.h>
#include <cassert>
#include <vector>
#include <iostream>

#include <XiaoTuMathBox/LinearAlgibra/Matrix.hpp>

namespace xiaotu {
namespace math {

    //! @brief Gauss-Jordan 消元发求解线性方程组 A x = b
    //!
    //! @param [inout] A 线性方程组的系数矩阵, 输出 A 的逆
    //! @param [inout] b 线性方程组右侧的向量们, 输出解们 x
    template <typename T>
    void GaussJordanEliminate(Matrix<T> & A, Matrix<T> * b = nullptr)
    {
        assert(A.Rows() == A.Cols());
        assert(nullptr == b || A.Rows() == b->Rows());

        int n = A.Rows();
        int m = (nullptr == b) ? 0 : b->Cols();
        // 记录发现 pivoting 的行列坐标
        std::vector<int> idxc(n), idxr(n);
        // 记录已经消元的行列
        std::vector<int> pivflag(n, 0);

        for (int t = 0; t < n; t++) {
            int pivRow = -1;
            int pivCol = -1;
            T max = 0;
            // 找出未消元的子阵中最大值作为 pivoting
            for (int ridx = 0; ridx < n; ridx++) {
                if (1 == pivflag[ridx])
                    continue;
                for (int cidx = 0; cidx < n; cidx++) {
                    if (1 == pivflag[cidx])
                        continue;
                    T abs = std::abs(A[ridx][cidx]);
                    if (abs > max) {
                        max = abs;
                        pivRow = ridx;
                        pivCol = cidx;
                    }
                }
            }
            pivflag[pivCol]++;
            idxc[t] = pivCol;
            idxr[t] = pivRow;

            if (0.0 == A[pivCol][pivRow])
                throw "奇异矩阵";

            // 交换 pivoting 到对角线上
            if (pivRow != pivCol) {
                A.RowSwap(pivRow, pivCol);
                if (nullptr != b)
                    b->RowSwap(pivRow, pivCol);
            }
            T pinv = 1.0 / A[pivCol][pivCol];
            A[pivCol][pivCol] = 1;

            // 归一化 pivoting 所在行
            for (int cidx = 0; cidx < n; cidx++)
                A[pivCol][cidx] *= pinv;
            for (int cidx = 0; cidx < m; cidx++)
                (*b)[pivCol][cidx] *= pinv;

            // 对其余各行消元
            for (int ridx = 0; ridx < n; ridx++) {
                if (ridx == pivCol)
                    continue;
                T dum = A[ridx][pivCol];
                A[ridx][pivCol] = 0.0;

                for (int cidx = 0; cidx < n; cidx++)
                    A[ridx][cidx] -= A[pivCol][cidx] * dum;
                for (int cidx = 0; cidx < m; cidx++)
                    (*b)[ridx][cidx] -= (*b)[pivCol][cidx] * dum;
            }
        }

        // 恢复逆矩阵
        for (int idx = n - 1; idx >= 0; idx--) {
            if (idxc[idx] == idxr[idx])
                continue;
            A.ColSwap(idxc[idx], idxr[idx]);
        }
    }

}
}

#endif
