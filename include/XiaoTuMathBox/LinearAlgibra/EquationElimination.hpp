#ifndef XTMB_LA_EQUATION_ELIMINATION_H
#define XTMB_LA_EQUATION_ELIMINATION_H

#include <stdint.h>
#include <cassert>
#include <vector>
#include <iostream>


namespace xiaotu::math {

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

            if (std::abs(A(pivCol, pivRow)) < SMALL_VALUE)
                throw std::runtime_error("奇异矩阵");

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

    //! @brief 高斯行消元, 得到阶梯矩阵
    //!
    //! @param [inout] A 消元矩阵
    //! @return 极大线性无关组
    template <typename MatViewA>
    std::vector<int> GaussRowEliminate(MatViewA & A)
    {
        typedef typename MatViewA::Scalar Scalar;
        int cols = A.Cols();
        int rows = A.Rows();
        std::vector<int> max_indep_set;

        for (int cidx = 0; cidx < cols; cidx++) {
            if (max_indep_set.size() >= rows)
                break;

            // 按行查找主元
            int pivot_row = max_indep_set.size();
            for (int ridx = max_indep_set.size(); ridx < rows; ridx++) {
                if (std::abs(A(ridx, cidx)) > std::abs(A(pivot_row, cidx)))
                    pivot_row = ridx;
            }
            Scalar pabs = std::abs(A(pivot_row, cidx));
            if (pabs <= SMALL_VALUE)
                continue;

            // 交换行
            if (pivot_row != max_indep_set.size())
                A.RowSwap(max_indep_set.size(), pivot_row);
            // 消元
            Scalar pinv = 1.0 / A(max_indep_set.size(), cidx);
            for (int i = cidx; i < cols; i++) {
                A(max_indep_set.size(), i) *= pinv;
            }
            for (int ridx = 0; ridx < rows; ridx++) {
                if (max_indep_set.size() == ridx)
                    continue;
                Scalar factor = A(ridx, cidx);
                for (int i = cidx; i < cols; i++)
                    A(ridx, i) -= A(max_indep_set.size(), i) * factor;
            }
            max_indep_set.push_back(cidx);
        }

        return max_indep_set;
    }    

    //! @brief 根据阶梯矩阵写出齐次线性方程组的解空间
    //!
    //! @param [in] 最大线性无关组列索引
    //! @param [in] 已经过高斯消元后的阶梯矩阵
    //! @return 解空间
    template <typename MatViewA>
    DMatrix<typename MatViewA::Scalar>
    SolveSpace(std::vector<int> const & indep_set, MatViewA const & A)
    {
        int s_rows = A.Cols();
        int s_cols = A.Cols() - indep_set.size();
        DMatrix<typename MatViewA::Scalar> re(s_rows, s_cols);

        int s_c = 0;
        int last_indep_idx = 0;
        for (int c = 0; c < s_rows; c++) {
            if (c == indep_set[last_indep_idx]) {
                last_indep_idx++;
                continue;
            }

            re(c, s_c) = 1;
            for (int r = 0; r < last_indep_idx; r++) {
                re(indep_set[r], s_c) = -A(r, c);
            }
            s_c++;
        }
        return re;
    }

    //! @brief 查找矩阵 A 中最大线性无关列向量组
    //!
    //! @param [in] A 待查找矩阵
    //! @return 最大线性无关向量组
    template <typename MatViewA>
    std::vector<int> FindMaximalIndepSet(MatViewA const & A)
    {
        DMatrix<typename MatViewA::Scalar> tmp(A.Rows(), A.Cols());
        tmp = A;
        return GaussRowEliminate(tmp);
    }

    //! @brief 查找矩阵 A 中最大线性无关列向量组
    //!
    //! @param [in] A 待查找矩阵
    //! @return 最大线性无关向量组
    template <typename MatViewA>
    int Rank(MatViewA const & A)
    {
        DMatrix<typename MatViewA::Scalar> tmp(A.Rows(), A.Cols());
        tmp = A;
        return GaussRowEliminate(tmp).size();
    }


}

#endif
