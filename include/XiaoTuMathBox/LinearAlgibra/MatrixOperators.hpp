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

    //! @brief 矩阵的数乘 Re = aA
    //!
    //! 适用于 MatrixView, Matrix
    //!
    //! @param [in] a 数值 a
    //! @param [in] A 矩阵 A
    //! @param [out] R R = aA
    //!
    //! @return 矩阵尺寸是否合法
    template <typename Scalar, typename MatrixA, typename MatrixRe>
    bool ScalarMultiply(Scalar const & a, MatrixA const & A, MatrixRe & R)
    {
        if (A.Cols() != R.Cols() || A.Rows() != R.Rows())
            return false;

        int m = R.Rows();
        int n = R.Cols();
        for (int ridx = 0; ridx < m; ++ridx)
            for (int cidx = 0; cidx < n; ++cidx)
                R(ridx, cidx) = a * A(ridx, cidx);

        return true;
    }
    
    //! @brief 矩阵的加法 Re = A + B
    //!
    //! 适用于 MatrixView, Matrix
    //!
    //! @param [in] A 矩阵 A
    //! @param [in] B 矩阵 B
    //! @param [out] R R = A + B
    //!
    //! @return 矩阵尺寸是否合法
    template <typename MatrixA, typename MatrixB, typename MatrixRe>
    bool Add(MatrixA const & A, MatrixB const & B, MatrixRe & R)
    {
        if (A.Cols() != B.Cols() || A.Rows() != B.Rows() ||
            R.Cols() != B.Cols() || R.Rows() != B.Rows())
            return false;

        int m = R.Rows();
        int n = R.Cols();
        for (int ridx = 0; ridx < m; ++ridx)
            for (int cidx = 0; cidx < n; ++cidx)
                R(ridx, cidx) = A(ridx, cidx) + B(ridx, cidx);

        return true;
    }

    //! @brief 矩阵的减法 Re = A - B
    //!
    //! 适用于 MatrixView, Matrix
    //!
    //! @param [in] A 矩阵 A
    //! @param [in] B 矩阵 B
    //! @param [out] R R = A - B
    //!
    //! @return 矩阵尺寸是否合法
    template <typename MatrixA, typename MatrixB, typename MatrixRe>
    bool Sub(MatrixA const & A, MatrixB const & B, MatrixRe & R)
    {
        if (A.Cols() != B.Cols() || A.Rows() != B.Rows() ||
            R.Cols() != B.Cols() || R.Rows() != B.Rows())
            return false;

        int m = R.Rows();
        int n = R.Cols();
        for (int ridx = 0; ridx < m; ++ridx)
            for (int cidx = 0; cidx < n; ++cidx)
                R(ridx, cidx) = A(ridx, cidx) - B(ridx, cidx);

        return true;
    }

    //! @brief 矩阵的乘法 Re = AB
    template <typename MatrixA, typename MatrixB>
    DMatrix<typename MatrixA::Scalar>
    operator * (MatrixA const & A, MatrixB const & B)
    {
        DMatrix<typename MatrixA::Scalar> re(A.Rows(), B.Cols());
        bool success = Multiply(A, B, re);
        assert(success);
        return re;
    }


    //! @brief 矩阵的乘法 Re = aA
    template <typename Matrix>
    Matrix operator * (typename Matrix::Scalar const & a, Matrix const & A)
    {
        Matrix re;
        bool success = ScalarMultiply(a, A, re);
        assert(success);
        return re;
    }

    //! @brief 矩阵的乘法 A *= a
    template <typename Matrix>
    bool operator *= (Matrix & A, typename Matrix::Scalar const & a)
    {
        return ScalarMultiply(a, A, A);
    }

    //! @brief 矩阵的乘法 Re = aA
    template <typename Matrix>
    Matrix operator * (Matrix const & A, typename Matrix::Scalar const & a)
    {
        Matrix re;
        bool success = ScalarMultiply(a, A, re);
        assert(success);
        return re;
    }

    //! @brief 矩阵的加法 Re = A + B
    template <typename MatrixA, typename MatrixB>
    Matrix<typename MatrixA::Scalar, MatrixA::NumRows, MatrixB::NumCols>
    operator + (MatrixA const & A, MatrixB const & B)
    {
        Matrix<typename MatrixA::Scalar, MatrixA::NumRows, MatrixB::NumCols> re;
        bool success = Add(A, B, re);
        assert(success);
        return re;
    }

    //! @brief 矩阵的加法 A += B
    template <typename MatrixA, typename MatrixB>
    MatrixA operator += (MatrixA & A, MatrixB const & B)
    {
        bool success = Add(A, B, A);
        assert(success);
        return A;
    }

    //! @brief 矩阵的减法 Re = A - B
    template <typename MatrixA, typename MatrixB>
    Matrix<typename MatrixA::Scalar, MatrixA::NumRows, MatrixB::NumCols>
    operator - (MatrixA const & A, MatrixB const & B)
    {
        Matrix<typename MatrixA::Scalar, MatrixA::NumRows, MatrixB::NumCols> re;
        bool success = Sub(A, B, re);
        assert(success);
        return re;
    }

    //! @brief 矩阵的减法 A -= B
    template <typename MatrixA, typename MatrixB>
    MatrixA operator -= (MatrixA & A, MatrixB const & B)
    {
        bool success = Sub(A, B, A);
        assert(success);
        return A;
    }

}

#endif
