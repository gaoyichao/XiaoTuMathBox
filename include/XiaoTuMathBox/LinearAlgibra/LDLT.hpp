#ifndef XTMB_LA_LDLT_H
#define XTMB_LA_LDLT_H

#include <cassert>
#include <vector>
#include <cmath>
#include <exception>


namespace xiaotu::math {

    //! @brief 给定一个 n x n 的对称正定矩阵 A 进行 Cholesky 分解。A = L * L^T
    //! 结果记录在成员变量 mLD 中。
    //!
    //! | 1                | | d0         | | 1    l_10 l_20 l_30 |   | a_00 a_10 a_20 a_30 |
    //! | l_10 1           | |    d1      | |      1    l_21 l_31 | = | a_10 a_11 a_21 a_31 |
    //! | l_20 l_21 1      | |       d2   | |           1    l_32 |   | a_20 a_21 a_22 a_32 |
    //! | l_30 l_31 l_32 1 | |          d3| |                1    |   | a_30 a_31 a_32 a_33 |

    template <typename MatrixViewA>
    class LDLT {
        public:
            typedef typename MatrixViewA::Scalar Scalar;
            constexpr static int N = MatrixViewA::NumRows;

            LDLT(MatrixViewA const & a)
            {
                assert(a.Rows() == a.Cols());
                mL.Assign(a);
                mD.resize(N);
                Decompose();
            }

            //! @brief 求解方程组 Ax=b, b 和 x 可以是同一个对象
            //! 
            //! @param [in] b 方程右侧的列向量
            //! @param [out] x 对应 b 中每一列的解
            template <typename MatrixViewB>
            void Solve(MatrixViewB const & b, MatrixViewB & x)
            {
                assert(b.Rows() == mL.Rows());
                const int M = b.Cols();
                x.Assign(b);

                // 遍历每一列
                for (int cidx = 0; cidx < M; ++cidx) {
                    for (int ridx = 0; ridx < N; ++ridx) {
                        Scalar y = x(ridx, cidx);
                        for (int k = 0; k < ridx; k++)
                            y -= mL(ridx, k) * x(k, cidx);
                        x(ridx, cidx) = y;
                    }

                    for (int ridx = N-1; ridx >= 0; --ridx) {
                        Scalar s = x(ridx, cidx) / mD[ridx];
                        for (int k = ridx+1; k < N; ++k)
                            s -= mL(k, ridx) * x(k, cidx);
                        x(ridx, cidx) = s;
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
                for (int ridx = 0; ridx < N; ++ridx) {
                    for (int cidx = ridx; cidx < N; ++cidx) {
                        Scalar sum = mL(ridx, cidx);
                        for (int k = ridx - 1; k >= 0; --k)
                            sum -= mL(ridx, k) * mD[k] * mL(cidx, k);
                        if (ridx == cidx) {
                            if (sum <= 0.0)
                                throw std::runtime_error("非正定矩阵");
                            mD[ridx] = sum;
                            mL(ridx, ridx) = 1;
                        } else {
                            mL(cidx, ridx) = sum / mD[ridx];
                        }
                    }
                }

                for (int ridx = 0; ridx < N; ++ridx)
                    for (int cidx = 0; cidx < ridx; ++cidx)
                        mL(cidx, ridx) = 0;
            }

        public:
            //! @brief 获取 L 矩阵
            Matrix<Scalar, N, N> L() { return mL; }

            //! @brief 获取 D 矩阵
            Matrix<Scalar, N, N> D()
            {
                auto re = Matrix<Scalar, N, N>::Eye();
                for (int ridx = 0; ridx < N; ++ridx)
                    re(ridx, ridx) = mD[ridx];
                return re;
            }

            //! @brief 获取 L 的转置矩阵
            Matrix<Scalar, N, N> LT()
            {
                auto  re = Matrix<Scalar, N, N>::Eye();

                for (int ridx = 0; ridx < N; ++ridx)
                    for (int cidx = 0; cidx < ridx; ++cidx)
                        re(cidx, ridx) = mL(ridx, cidx);
                return re;
            }

        private:
            Matrix<Scalar, N, N> mL;
            std::vector<Scalar> mD;
    };

}


#endif
