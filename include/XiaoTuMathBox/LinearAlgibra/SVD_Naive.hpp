#ifndef XTMB_LA_SVD_Naive_H
#define XTMB_LA_SVD_Naive_H

#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

namespace xiaotu::math {

    /**
     * @brief Golub-Kahan 算法,进行 SVD 分解
     * 
     * \Sigma = U^T A V
     */
    template <typename MatViewIn>
    class SVD_Naive {
        public:
            typedef typename MatViewIn::Scalar Scalar;
            typedef DMatrix<Scalar> Mat;

            /**
             * @brief 默认构造函数
             */
            SVD_Naive(MatViewIn const & a, int max_iter, Scalar tolerance)
            {
                mSigma.Resize(a.Rows(), a.Cols()) = a;
                mUT.Resize(a.Rows(), a.Rows()).Identity();
                mV.Resize(a.Cols(), a.Cols()).Identity();

                mSigma.Bidiagonal(&mUT, &mV);

                int m = mSigma.Rows();
                int n = mSigma.Cols();
                int p = m < n ? m : n;

                DMatrix<Scalar> B(p, p);
                B = mSigma.SubMatrix(0, 0, p, p);
                auto BTB = B.Transpose() * B;

                EigenImplicitQR<DMatrix<Scalar>> qr;
                qr.Iterate(BTB, max_iter, tolerance);

                auto BQ = B * qr.Q().Transpose();
                QR_Householder pqr(BQ);
                mSigma.SubMatrix(0, 0, p, p) = pqr.R();

                auto ut = mUT.SubMatrix(0, 0, p, m);
                auto v = mV.SubMatrix(0, 0, n, p);
                ut = pqr.Q().Transpose() * ut;
                v = (v * qr.Q().Transpose());
            }

            DMatrix<Scalar> const Sigma() { return mSigma; }
            DMatrix<Scalar> const UT() { return mUT; }
            DMatrix<Scalar> const V() { return mV; }

        private:
            DMatrix<Scalar> mUT;
            DMatrix<Scalar> mSigma;
            DMatrix<Scalar> mV;
    };

}

#endif


