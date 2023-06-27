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
    template <typename DerivedA, typename DerivedB = MatrixView<double, 4, 4>>
    void GaussJordanEliminate(MatrixViewBase<DerivedA> & A, MatrixViewBase<DerivedB> * b = nullptr)
    {
        typedef typename traits<DerivedA>::Scalar Scalar;

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
            Scalar max = 0;
            // 找出未消元的子阵中最大值作为 pivoting
            for (int ridx = 0; ridx < n; ridx++) {
                if (1 == pivflag[ridx])
                    continue;
                for (int cidx = 0; cidx < n; cidx++) {
                    if (1 == pivflag[cidx])
                        continue;
                    Scalar abs = std::abs(A.template At(ridx, cidx));
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

            if (0.0 == A.template At(pivCol, pivRow))
                throw "奇异矩阵";

            // 交换 pivoting 到对角线上
            if (pivRow != pivCol) {
                A.RowSwap(pivRow, pivCol);
                if (nullptr != b)
                    b->RowSwap(pivRow, pivCol);
            }
            Scalar pinv = 1.0 / A.template At(pivCol, pivCol);
            A.template At(pivCol, pivCol) = 1;

            // 归一化 pivoting 所在行
            for (int cidx = 0; cidx < n; cidx++)
                A.template At(pivCol, cidx) *= pinv;
            for (int cidx = 0; cidx < m; cidx++)
                (*b).template At(pivCol, cidx) *= pinv;

            // 对其余各行消元
            for (int ridx = 0; ridx < n; ridx++) {
                if (ridx == pivCol)
                    continue;
                Scalar dum = A.template At(ridx, pivCol);
                A.template At(ridx, pivCol) = 0.0;

                for (int cidx = 0; cidx < n; cidx++)
                    A.template At(ridx, cidx) -= A.template At(pivCol, cidx) * dum;
                for (int cidx = 0; cidx < m; cidx++)
                    (*b).template At(ridx, cidx) -= (*b).template At(pivCol, cidx) * dum;
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
