#ifndef XTMB_LA_QR_HOUSEHOLDER_H
#define XTMB_LA_QR_HOUSEHOLDER_H

namespace xiaotu::math {

    //! @brief 基于 Householder 的 QR 分解
    //!
    //! 将输入的矩阵，分解为一个正交矩阵 Q 和一个上三角矩阵 R
    template <typename MatViewIn>
    class QR_Householder {
        public:
            typedef typename MatViewIn::Scalar Scalar;

            QR_Householder(MatViewIn const & a)
                : mQ(a.Rows(), a.Rows()), mR(a.Rows(), a.Cols())
            {
                mQ.Identity();
                mR = a;

                Decompose(a);
            }

        private:

            void Decompose(MatViewIn const & a)
            {
                int m = a.Rows();
                int n = a.Cols();
                
                for (int k = 0; k < (m - 1); k++) {
                    auto A_k = mR.SubMatrix(k, k, m - k, n - k);
                    auto a_1 = A_k.ColView(0);

                    auto H = DMatrix<Scalar>::Eye(m, m);
                    auto H_k = H.SubMatrix(k, k, m - k, m - k);

                    auto v = HouseholderVector(a_1);
                    HouseholderMatrix(v, H_k);

                    A_k = H_k * A_k;
                    mQ = mQ * H;
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

