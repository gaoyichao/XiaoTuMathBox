#ifndef XTMB_LA_QR_GRAM_SCHMIDT_H
#define XTMB_LA_QR_GRAM_SCHMIDT_H

#include <cassert>
#include <vector>
#include <cmath>

namespace xiaotu::math {

    //! @brief 基于 Gram-Schmidt 的 QR 分解
    //!
    //! 将输入的矩阵，分解为一个正交矩阵 Q 和一个上三角矩阵 R
    //!
    //! 由于数值上的不稳定，工程上基本不用该算法
    //! 本实现只是出于教学目的编写的，只支持对方阵的 QR 分解
    template <typename MatViewIn>
    class QR_GramSchmidt {
        public:
            typedef typename MatViewIn::Scalar Scalar;

            QR_GramSchmidt(MatViewIn const & a)
                : mQ(a.Rows(), a.Cols()), mR(a.Rows(), a.Cols())
            {
                assert(a.Rows() == a.Cols());

                mQ.Zeroing();
                mR.Zeroing();

                Decompose(a);
            }

        private:

            void Decompose(MatViewIn const & a)
            {
                int num = a.Cols();
                auto vec_0 = a.ColView(0);
                auto ort_0 = mQ.ColView(0);

                mR(0, 0) = std::sqrt(vec_0.Dot(vec_0));
                ort_0 = vec_0 / mR(0, 0);

                for (int k = 1; k < num; k++) {
                    auto vec_k = a.ColView(k);
                    auto ort_k = mQ.ColView(k);
                    ort_k = vec_k;

                    for (int i = 0; i < k; i++) {
                        auto ort_i = mQ.ColView(i);
                        mR(i, k) = vec_k.Dot(ort_i);
                        ort_k = ort_k - mR(i, k) * ort_i;
                    }
                    mR(k, k) = std::sqrt(ort_k.Dot(ort_k));
                    ort_k = ort_k / mR(k, k);
                }
            }

        public:
            DMatrix<Scalar> const & Q() { return mQ; }
            DMatrix<Scalar> const & R() { return mR; }

        private:
            DMatrix<Scalar> mQ;
            DMatrix<Scalar> mR;

    };
}

#endif


