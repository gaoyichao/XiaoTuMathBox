#ifndef XTMB_LA_EIGEN_NAIVE_QR_H
#define XTMB_LA_EIGEN_NAIVE_QR_H

/**
 * 各种关于 QR 算法的实现，用于求解矩阵的特征值
 */


namespace xiaotu {

    //! @brief 最基础的 QR 算法, 出于教学目的编写的, 效率较低
    //!
    //! 将输入的矩阵 A，转换成 A = U T U^T
    //! 其中矩阵 U 为一系列正交矩阵 Q 的右乘
    //! T 为收敛后的矩阵，它的对角线就是特征值
    template <typename MatViewIn>
    class EigenNaiveQR {
        public:
            typedef typename MatViewIn::Scalar Scalar;

            /**
             * @brief 构造函数
             * 
             * @param [in] a 目标矩阵
             * @param [in] max_iter 最大迭代次数
             */
            EigenNaiveQR(MatViewIn const & a, int max_iter = 100)
                : mT(a.Rows(), a.Cols()), mU(a.Rows(), a.Cols())
            {
                assert(a.Rows() == a.Cols());
                mT = a;
                mU.Identity();

                Iterate(max_iter);
            }

        private:

            /**
             * @brief QR 迭代
             * 
             * @param [in] max_iter 最大迭代次数
             */
            void Iterate(int max_iter)
            {
                QR_Householder<DMatrix<Scalar>> qr;
                for (int i = 0; i < max_iter; i++) {
                    qr.Decompose(mT);
                    mT = qr.R() * qr.Q();
                    mU = mU * qr.Q();
                }
            }


        public:
            DMatrix<Scalar> const & T() { return mT; }
            DMatrix<Scalar> const & U() { return mU; }

            /**
             * @brief 从矩阵 mT 中获取对角元素，即特征值
             */
            std::vector<Scalar> EigenValues()
            {
                int n = mT.Rows();
                std::vector<Scalar> re(n);

                for (int i = 0; i < n; i++)
                    re[i] = mT(i, i);

                return re;
            }

        private:
            DMatrix<Scalar> mT;
            DMatrix<Scalar> mU;
    };

}

#endif

