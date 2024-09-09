#ifndef XTMB_LA_CHOLESKY_H
#define XTMB_LA_CHOLESKY_H

#include <cassert>
#include <vector>
#include <cmath>
#include <exception>

namespace xiaotu::math {

    //! @brief 给定一个 n x n 的对称正定矩阵 A 进行 Cholesky 分解。A = L * L^T
    //! 结果记录在成员变量 ml 中。
    //!
    //! | l_00                | | l_00 l_10 l_20 l_30 |   | a_00 a_10 a_20 a_30 |
    //! | l_10 l_11           | |      l_11 l_21 l_31 | = | a_10 a_11 a_21 a_31 |
    //! | l_20 l_21 l_22      | |           l_22 l_32 |   | a_20 a_21 a_22 a_32 |
    //! | l_30 l_31 l_32 l_33 | |                l_33 |   | a_30 a_31 a_32 a_33 |

    template <typename MatrixViewA>
    class Cholesky {
        public:
            typedef typename MatrixViewA::Scalar Scalar;

            Cholesky(MatrixViewA const & a)
                : mBuffer(a.NumDatas()),
                  ml(mBuffer.data())
            {
                assert(a.Rows() == a.Cols());
                ml.Assign(a);
                Decompose();
            }

            //! @brief 求解方程组 Ax=b, b 和 x 可以是同一个对象
            //! 
            //! @param [in] b 方程右侧的列向量
            //! @param [out] x 对应 b 中每一列的解
            template <typename MatrixViewB>
            void Solve(MatrixViewB const & b, MatrixViewB & x)
            {
                assert(b.Rows() == ml.Rows());
            
                const int N = ml.Rows();
                const int M = b.Cols();
                x.Assign(b);

                // 遍历每一列
                for (int cidx = 0; cidx < M; ++cidx) {
                    for (int ridx = 0; ridx < N; ++ridx) {
                        Scalar y = x(ridx, cidx);
                        for (int k = 0; k < ridx; k++)
                            y -= ml(ridx, k) * x(k, cidx);
                        x(ridx, cidx) = y / ml(ridx, ridx);
                    }

                    for (int ridx = N-1; ridx >= 0; --ridx) {
                        Scalar s = x(ridx, cidx);
                        for (int k = ridx+1; k < N; ++k)
                            s -= ml(k, ridx) * x(k, cidx);
                        x(ridx, cidx) = s / ml(ridx, ridx);
                    }
                }
            }

            void Inverse(MatrixViewA  & inv)
            {
                inv.Identity();
                Solve(inv, inv);
            }

        private:
            void Decompose()
            {
                const int N = ml.Rows();

                for (int ridx = 0; ridx < N; ++ridx) {
                    for (int cidx = ridx; cidx < N; ++cidx) {
                        Scalar sum = ml(ridx, cidx);
                        for (int k = ridx - 1; k >= 0; --k)
                            sum -= ml(ridx, k) * ml(cidx, k);
                        if (ridx == cidx) {
                            if (sum <= 0.0)
                                throw std::runtime_error("非正定矩阵");
                            ml(ridx, ridx) = std::sqrt(sum);
                        } else {
                            ml(cidx, ridx) = sum / ml(ridx, ridx);
                        }
                    }
                }

                for (int ridx = 0; ridx < N; ++ridx)
                    for (int cidx = 0; cidx < ridx; ++cidx)
                        ml(cidx, ridx) = 0;
            }


        public:
            //! @brief 获取 L 分解矩阵
            MatrixViewA & operator() () { return ml; }
            MatrixViewA const & operator() () const { return ml; }

        private:
            std::vector<Scalar> mBuffer;
            MatrixViewA ml;


    };
}

#endif
