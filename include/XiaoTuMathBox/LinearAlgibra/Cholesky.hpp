#ifndef XTMB_LA_CHOLESKY_H
#define XTMB_LA_CHOLESKY_H

#include <cassert>
#include <vector>

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

        private:
            void Decompose()
            {

            }

        private:
            std::vector<Scalar> mBuffer;
            MatrixViewA ml;


    };
}

#endif
