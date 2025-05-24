#ifndef XTMB_LA_MATRIX_OPERATORS_H
#define XTMB_LA_MATRIX_OPERATORS_H

#include <cassert>
#include <algorithm>
#include <vector>
#include <cmath>
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
    template <typename MatrixA, typename MatrixB, bool AIsMatrix = MatrixA::IsMatrix, bool BIsMatrix = MatrixB::IsMatrix>
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
    DMatrix<typename Matrix::Scalar>
    operator * (typename Matrix::Scalar const & a, Matrix const & A)
    {
        DMatrix<typename Matrix::Scalar> re(A.Rows(), A.Cols());
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
    DMatrix<typename Matrix::Scalar>
    operator * (Matrix const & A, typename Matrix::Scalar const & a)
    {
        DMatrix<typename Matrix::Scalar> re(A.Rows(), A.Cols());
        bool success = ScalarMultiply(a, A, re);
        assert(success);
        return re;
    }

    template <typename Matrix>
    DMatrix<typename Matrix::Scalar>
    operator / (Matrix const & A, typename Matrix::Scalar const & a)
    {
        auto a_inv = 1 / a;
        return A * a_inv;
    }

    //! @brief 矩阵的加法 Re = A + B
    template <typename MatrixA, typename MatrixB, bool AIsMatrix = MatrixA::IsMatrix, bool BIsMatrix = MatrixB::IsMatrix>
    DMatrix<typename MatrixA::Scalar>
    operator + (MatrixA const & A, MatrixB const & B)
    {
        DMatrix<typename MatrixA::Scalar> re(A.Rows(), B.Cols());
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
    template <typename MatrixA, typename MatrixB, bool AIsMatrix = MatrixA::IsMatrix, bool BIsMatrix = MatrixB::IsMatrix>
    DMatrix<typename MatrixA::Scalar>
    operator - (MatrixA const & A, MatrixB const & B)
    {
        DMatrix<typename MatrixA::Scalar> re(A.Rows(), B.Cols());
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


    /////////////////////////////////////////////////////////////////////////
    //
    // 一些子阵的操作
    //
    /////////////////////////////////////////////////////////////////////////

    //! @brief 获取 [begin, end] 中除 list 之外的元素
    //!
    //! @param [in] list 将要排除的索引列表
    //! @param [in] begin 全集的起始
    //! @param [in] end 全集的结束
    template <typename T>
    std::vector<T> Complement(std::vector<T> const & list, T begin, T end)
    {
        assert(begin < end);

        std::vector<T> b_list = list;
        std::sort(b_list.begin(), b_list.end());

        size_t idx = 0;
        std::vector<T> re;
        for (T c = begin; c < end; c++) {
            if (c == b_list[idx]) {
                idx++;
                continue;
            }
            re.push_back(c);
        }
        return re;
    }

    //! @brief 获取矩阵 A 中的指定列
    //!
    //! @param [in] row_list 指定列索引列表
    //! @param [in] A 目标矩阵 A
    //! @param [out] R 已经完成内存分配的结果缓存
    template <typename MatrixA, typename MatrixRe>
    void SelectCol(std::vector<int> const & col_list, MatrixA const & A, MatrixRe & R)
    {
        assert(A.Rows() == R.Rows());
        assert(col_list.size() == R.Cols());

        for (size_t idx = 0; idx < col_list.size(); idx++) {
            int c = col_list[idx];
            for (int r = 0; r < A.Rows(); r++)
                R(r, idx) = A(r, c);
        }
    }

    //! @brief 获取矩阵 A 中的指定列
    //!
    //! @param [in] row_list 指定列索引列表
    //! @param [in] A 目标矩阵 A
    //! @return 指定列的矩阵
    template <typename MatrixA>
    DMatrix<typename MatrixA::Scalar>
    SelectCol(std::vector<int> const & col_list, MatrixA const & A)
    {
        DMatrix<typename MatrixA::Scalar> re(A.Rows(), col_list.size());
        SelectCol(col_list, A, re);
        return re;
    }


    //! @brief 获取矩阵 A 中的除指定列外的子阵
    //!
    //! @param [in] list 排除列索引列表
    //! @param [in] A 目标矩阵 A
    //! @return  A 中的除指定列外的子阵
    template <typename MatrixA>
    DMatrix<typename MatrixA::Scalar>
    SelectNotCol(std::vector<int> const & list, MatrixA const & A)
    {
        std::vector<int> col_list = Complement(list, 0, A.Cols());
        return SelectCol(col_list, A);
    }

    /////////////////////////////////////////////////////////////////////////
    //
    // 一些特别的计算
    //
    /////////////////////////////////////////////////////////////////////////

    //! @brief Gram-Schmidt 标准正交化
    //!
    //! @param [in] col_vectors 按照列向量的形式排列的基向量
    //! @param [out] ortho 标准正交基
    template <typename MatViewIn, typename MatViewOut>
    void GramSchmidt(MatViewIn const & col_vectors, MatViewOut & ortho)
    {
        int num = col_vectors.Cols();
        int dim = col_vectors.Rows();
        typedef typename Traits<MatViewIn>::Scalar Scalar;
        std::vector<Scalar> norms(num);

        for (int k = 0; k < num; k++) {
            auto vec_k = col_vectors.SubMatrix(0, k, dim, 1);
            auto ort_k = ortho.SubMatrix(0, k, dim, 1);
            ort_k = vec_k;

            for (int i = 0; i < k; i++) {
                auto ort_i = ortho.SubMatrix(0, i, dim, 1);
                auto proj = vec_k.Dot(ort_i) / norms[i];
                ort_k = ort_k - proj * ort_i;
            }
            norms[k] = ort_k.Dot(ort_k);
        }

        for (int k = 0; k < num; k++) {
            auto ort_k = ortho.SubMatrix(0, k, dim, 1);
            ort_k = ort_k / std::sqrt(norms[k]);
        }
    }

    //! @brief 指定需要保留的元素索引，计算 Householder 向量
    //!
    //! @param [in] x 目标向量
    //! @param [out] v Householder 向量
    //! @param [in] n 需要保留的元素索引，缺省为 0
    template <typename VectorIn, typename VectorOut>
    void HouseholderVector(VectorIn const & x, VectorOut & v, int n = 0)
    {
        assert(n < x.NumDatas());
        typedef typename Traits<VectorIn>::Scalar Scalar;
        Scalar x_norm = x.Norm();
        Scalar alpha = -std::copysign(x_norm, x(n));

        v = x;
        v(n) -= alpha;
    }

    //! @brief 指定需要保留的元素索引，计算 Householder 向量
    //!
    //! @param [in] x 目标向量
    //! @param [in] n 需要保留的元素索引，缺省为 0
    //! @return Householder 向量
    template <typename VectorIn>
    DMatrix<typename VectorIn::Scalar>
    HouseholderVector(VectorIn const & x, int n = 0)
    {
        DMatrix<typename VectorIn::Scalar> re(x.NumDatas(), 1);
        HouseholderVector(x, re, n);
        return re;
    }

    //! @brief 构建关于 v 的 Householder 矩阵
    //!
    //! Householder 矩阵的几何意义是，将向量 x 关于一个垂直于向量 v 的超平面的镜面反射 Hx
    //! 在 QR 分解中，常用来消除 x 中初指定元素外的其它元素。
    //!
    //! @param [in] v 参考向量
    //! @param [out] H 通过 v 构造的 Householder 矩阵
    template <typename VectorV, typename MatrixH>
    void HouseholderMatrix(VectorV const & v, MatrixH & H)
    {
        assert(H.Rows() == H.Cols());
        assert(v.Rows() == H.Rows());

        typedef typename Traits<VectorV>::Scalar Scalar;

        Scalar v_norm = v.SquaredNorm();
        auto vvT = v * v.Transpose();

        H = DMatrix<Scalar>::Eye(H.Rows(), H.Cols()) - 2 / v_norm * vvT;
    }


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


}

#endif
