#ifndef XTMB_LA_MATRIX_EIGEN_VALUE_H
#define XTMB_LA_MATRIX_EIGEN_VALUE_H

#include <cassert>
#include <algorithm>
#include <vector>
#include <cmath>
#include <iostream>

namespace xiaotu::math {


    //! @brief 幂法求解特征值
    //!
    //! @param [in] 目标矩阵 nxn
    //! @param [in|out] 输入幂法迭代的初始向量，输出特征向量的近似向量
    //! @param [in] eps 当 re 的迭代变化量小于 eps 时终止迭代
    //! @param [in] max_iter 最大迭代次数
    //! @return 绝对值最大的那个特征值
    template <typename MatrixA, typename VectorV>
    typename MatrixA::Scalar
    PowerIterate(MatrixA const & A, VectorV & v, typename MatrixA::Scalar eps = SMALL_VALUE, int max_iter = 100)
    {
        assert(A.Rows() == A.Cols());
        assert(v.Rows() == A.Rows());

        typename MatrixA::Scalar re = 0;
        v = v / v.Norm();
        for (int i = 0; i < max_iter; i++) {
            auto w = A * v;
            auto r = v.Dot(w);
            v = w / w.Norm();

            auto delta = std::abs(r - re);
            re = r;
            if (delta < eps)
                break;
        }

        return re;
    }


    //! @brief 逆幂法求解特征值
    //!
    //! @param [in] 目标矩阵 nxn
    //! @param [in|out] 输入幂法迭代的初始向量，输出特征向量的近似向量
    //! @param [in] eps 当 re 的迭代变化量小于 eps 时终止迭代
    //! @param [in] max_iter 最大迭代次数
    //! @return 绝对值最小的那个特征值
    template <typename MatrixA, typename VectorV>
    typename MatrixA::Scalar
    InversePowerIterate(MatrixA const & A, VectorV & v, typename MatrixA::Scalar eps = SMALL_VALUE, int max_iter = 100)
    {
        using Scalar = typename MatrixA::Scalar;
        assert(A.Rows() == A.Cols());
        assert(v.Rows() == A.Rows());

        Scalar re = 0;
        LU lu(A);

        for (int i = 0; i < max_iter; i++) {
            lu.Solve(v, v);
            v = v / v.Norm();
            auto r = v.Dot(A * v);

            auto delta = std::abs(r - re);
            re = r;
            if (delta < eps)
                break;
        }

        return re;
    }




}

#endif
