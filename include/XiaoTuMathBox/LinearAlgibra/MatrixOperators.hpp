#ifndef XTMB_LA_MATRIX_OPERATORS_H
#define XTMB_LA_MATRIX_OPERATORS_H

#include <cassert>
#include <vector>
#include <iostream>

namespace xiaotu::math {

    //! @brief 矩阵 A 的转置
    //!
    //! @return 矩阵尺寸是否合法
    template <typename MatIn, typename MatOut>
    bool Transpose(MatIn const & A, MatOut & T)
    {
        if (A.Rows() != T.Cols() || A.Cols() != T.Rows())
            return false;

        int m = A.Rows();
        int n = A.Cols();
        for (int ridx = 0; ridx < m; ++ridx)
            for (int cidx = 0; cidx < n; ++cidx)
                T(cidx, ridx) = A(ridx, cidx);

        return true;
    }

    //! @brief 矩阵的乘法 Re = AB
    //!
    //! 适用于 MatrixView, Matrix
    //!
    //! @param [in] A 矩阵 A
    //! @param [in] B 矩阵 B
    //! @param [out] R R = AB
    //!
    //! @return 矩阵尺寸是否合法
    template <typename MatrixA, typename MatrixB, typename MatrixRe>
    bool Multiply(MatrixA const & A, MatrixB const & B, MatrixRe & R)
    {
        if (A.Cols() != B.Rows() || R.Rows() != A.Rows() || R.Cols() != B.Cols())
            return false;

        int l = A.Cols();
        int m = R.Rows();
        int n = R.Cols();
        for (int ridx = 0; ridx < m; ++ridx)
            for (int cidx = 0; cidx < n; ++cidx) {
                R(ridx, cidx) = 0;
                for (int k = 0; k < l; ++k)
                    R(ridx, cidx) += A(ridx, k) * B(k, cidx);
            }

        return true;
    }

    //! @brief 矩阵的乘法 Re = AB
    template <typename MatrixA, typename MatrixB>
    Matrix<typename MatrixA::Scalar, MatrixA::NumRows, MatrixB::NumCols>
    operator * (MatrixA const & A, MatrixB const B)
    {
        Matrix<typename MatrixA::Scalar, MatrixA::NumRows, MatrixB::NumCols> re;
        bool success = Multiply(A, B, re);
        assert(success);
        return re;
    }
}

#endif
