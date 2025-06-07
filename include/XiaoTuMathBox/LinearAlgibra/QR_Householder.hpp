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
                    auto a_1 = A_k.Col(0);

                    auto H = DMatrix<Scalar>::Eye(m, m);
                    auto H_k = H.SubMatrix(k, k, m - k, m - k);

                    auto v = HouseholderVector(a_1);
                    HouseholderMatrix(v, H_k);

                    A_k = H_k * A_k;

                    // 为了保证 QR 分解唯一，添加约束 R 的对角元素都是正数
                    // https://gaoyichao.com/Xiaotu/?book=algebra&title=QR算法的收敛原理
                    if (A_k(0, 0) < 0) {
                        auto a_k_1 = A_k.Row(0);
                        a_k_1 *= -1.0;
                        auto h_k_1 = H_k.Col(0);
                        h_k_1 *= -1.0;
                    }

                    mQ = mQ * H;
                }
                
                // 为了保证 QR 分解唯一，添加约束 R 的对角元素都是正数
                // https://gaoyichao.com/Xiaotu/?book=algebra&title=QR算法的收敛原理
                Scalar r_mm = mR(m-1, m-1);
                if (r_mm < 0) {
                    auto a_k_1 = mR.Row(m-1);
                    a_k_1 *= -1.0;
                    auto h_k_1 = mQ.Col(m-1);
                    h_k_1 *= -1.0;
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

