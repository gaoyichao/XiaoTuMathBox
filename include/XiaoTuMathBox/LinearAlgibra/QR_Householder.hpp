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

            /**
             * @brief 默认构造函数
             */
            QR_Householder()
            {}

            /**
             * @brief 构造函数
             * 
             * @param [in] a 待分解矩阵
             */
            QR_Householder(MatViewIn const & a)
                : mQ(a.Rows(), a.Rows()), mR(a.Rows(), a.Cols())
            {
                mQ.Identity();
                mR = a;
                __Decompose__();
            }

            /**
             * @brief 执行 QR 分解
             * 
             * 针对 QR 迭代的优化, 减少迭代过程中重新申请内存的次数
             * 
             * @param [in] a 待分解矩阵
             */
            void Decompose(MatViewIn const & a)
            {
                if (mR.Rows() != a.Rows() || mR.Cols() != a.Cols()) {
                    mQ.Resize(a.Rows(), a.Rows());
                    mR.Resize(a.Rows(), a.Cols());
                }

                mQ.Identity();
                mR = a;
                __Decompose__();
            }

        private:

            /**
             * @brief 具体的分解实现
             */
            void __Decompose__()
            {
                int m = mR.Rows();
                int n = mR.Cols();
                
                for (int k = 0; k < (m - 1); k++) {
                    auto A_k = mR.SubMatrix(k, k, m - k, n - k);
                    auto a_1 = A_k.Col(0);

                    auto H = DMatrix<Scalar>::Eye(m, m);
                    auto H_k = H.SubMatrix(k, k, m - k, m - k);

                    auto v = HouseholderVector(a_1);
                    HouseholderMatrix(v, H_k);

                    A_k = H_k * A_k;
                    mQ = mQ * H;
                }
                
                FixDiag();
            }
            /**
             * @brief 修正对角元素，保证它们为正数
             *
             * 为了保证 QR 分解唯一，添加约束 R 的对角元素都是正数
             * https://gaoyichao.com/Xiaotu/?book=algebra&title=QR算法的收敛原理
             */
            void FixDiag()
            {
                int rows = mR.Rows();
                int cols = mR.Cols();
                int n = cols < rows ? cols : rows;

                for (int i = 0; i < n; i++) {
                    if (mR(i, i) < 0) {
                        auto r = mR.Row(i);
                        r *= -1.0;
                        auto q = mQ.Col(i);
                        q *= -1.0;
                    }
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

