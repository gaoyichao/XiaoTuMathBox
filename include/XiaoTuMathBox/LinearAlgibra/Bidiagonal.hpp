#ifndef XTMB_LA_BI_DIAGONAL_H
#define XTMB_LA_BI_DIAGONAL_H

#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

namespace xiaotu::math {

    /**
     * @brief 对输入矩阵进行二对角化, B = U^T A V, B 是二对角矩阵
     */
    template <typename MatViewIn>
    class Bidiagonal {
        public:
            typedef typename MatViewIn::Scalar Scalar;

        public:

            /**
             * @brief 默认构造函数
             */
            Bidiagonal()
            {}

            /**
             * @brief 构造函数
             * 
             * @param [in] A 待二教化的矩阵
             */
            Bidiagonal(MatViewIn const & A)
                : mB(A.Rows(), A.Cols()),
                  mUT(A.Rows(), A.Rows()),
                  mV(A.Cols(), A.Cols())
            {
                mB = A;
                mUT.Identity();
                mV.Identity();

                if (A.Rows() >= A.Cols())
                    UpperBidiagonal();
                else
                    LowerBidiagonal();
            }
        private:

            void UpperBidiagonal()
            {
                int m = mB.Rows();
                int n = mB.Cols();

                for (int k = 0; k < n ; k++) {
                    auto A_k = mB.SubMatrix(k, k, m - k, n - k);
                    auto u = HouseholderVector(A_k.Col(0));
                    auto H = DMatrix<Scalar>::Eye(m-k, m-k);
                    HouseholderMatrix(u, H);
                    A_k = H * A_k;

                    auto UT_k = mUT.SubMatrix(k, 0, m - k, m);
                    UT_k = H * UT_k;

                    if (k < (n - 1)) {
                        auto v = HouseholderVector(A_k.SubMatrix(0, 1, 1, n-k-1));
                        auto vH = DMatrix<Scalar>::Eye(n-k, n-k);
                        auto H_k = vH.SubMatrix(1, 1, n-k-1, n- k-1);
                        HouseholderMatrix_Row(v, H_k);
                        A_k = A_k * vH;

                        auto V_k = mV.SubMatrix(0, k, n, n - k);
                        V_k = V_k * vH;
                    }
                }
            }

            void LowerBidiagonal()
            {
                int m = mB.Rows();
                int n = mB.Cols();

                for (int k = 0; k < m ; k++) {
                    auto A_k = mB.SubMatrix(k, k, m - k, n - k);
                    auto v = HouseholderVector(A_k.Row(0));
                    auto H = DMatrix<Scalar>::Eye(n-k, n-k);
                    HouseholderMatrix_Row(v, H);
                    A_k = A_k * H;

                    auto V_k = mV.SubMatrix(0, k, n, n-k);
                    V_k = V_k * H;

                    if (k < (m - 1)) {
                        auto u = HouseholderVector(A_k.SubMatrix(1, 0, m-k-1, 1));
                        auto uH = DMatrix<Scalar>::Eye(m-k, m-k);
                        auto H_k = uH.SubMatrix(1, 1, m-k-1, m-k-1);
                        HouseholderMatrix(u, H_k);

                        A_k = uH * A_k;

                        auto UT_k = mUT.SubMatrix(k, 0, m-k, m);
                        UT_k = uH * UT_k;
                    }
                }
            }

        public:
            DMatrix<Scalar> const & B() { return mB; }
            DMatrix<Scalar> const & UT() { return mUT; }
            DMatrix<Scalar> const & V() { return mV; }

        private:
            //! @brief B = U^T A V. 输出的二对角矩阵
            DMatrix<Scalar> mB;
            //! @brief B = U^T A V
            DMatrix<Scalar> mUT;
            //! @brief B = U^T A V
            DMatrix<Scalar> mV;
    };

}

#endif
