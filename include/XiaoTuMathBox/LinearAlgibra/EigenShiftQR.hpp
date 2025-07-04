#ifndef XTMB_LA_EIGEN_SHIFT_QR_H
#define XTMB_LA_EIGEN_SHIFT_QR_H

namespace xiaotu::math {

    /**
     * @brief 带位移的 QR 算法, 出于教学目的编写的, 效率较低
     * 
     * 将输入的矩阵 A，转换成 A = U T U^T
     * 其中矩阵 U 为一系列正交矩阵 Q 的右乘
     * T 为收敛后的矩阵，它的对角线就是特征值
     */
    template <typename MatViewIn>
    class EigenShiftQR {
        public:
            typedef typename MatViewIn::Scalar Scalar;

            /**
             * @brief 默认构造函数
             */
            EigenShiftQR()
            {
            }

            /**
             * @brief 构造函数
             * 
             * @param [in] a 目标矩阵
             * @param [in] max_iter 最大迭代次数
             * @param [in] tolerance 结束迭代的容忍度
             */
            EigenShiftQR(MatViewIn const & a, int max_iter = 100, Scalar tolerance = SMALL_VALUE)
            {
                Iterate(a, max_iter, tolerance);
            }

            /**
             * @brief 带位移的 QR 迭代
             * 
             * @param [in] a 目标矩阵
             * @param [in] max_iter 最大迭代次数
             * @param [in] tolerance 结束迭代的容忍度
             * @return 迭代次数
             */
            int Iterate(MatViewIn const & a, int max_iter, Scalar tolerance)
            {
                assert(a.Rows() == a.Cols());

                int n = a.Rows();
                mT.Resize(n, n);

                DMatrix<Scalar> tmp0 = a;
                DMatrix<Scalar> tmp1 = DMatrix<Scalar>::Zero(n, n);
                DMatrix<Scalar> * pA0 = &tmp0;
                DMatrix<Scalar> * pA1 = &tmp1;
                QR_Householder<DMatrix<Scalar>> qr;
                Scalar offset = (*pA0)(n - 1, n - 1);

                int i = 0;
                for (; i < max_iter; i++) {
                    SubDiagScalar(*pA0, offset);
                    qr.Decompose(*pA0);
                    (*pA1) = qr.R() * qr.Q();

                    if (IsConverged(*pA0, *pA1, tolerance)) {
                        AddDiagScalar(*pA1, offset);
                        mT = *pA1;
                        return i;
                    }

                    AddDiagScalar(*pA1, offset);
                    std::swap(pA0, pA1);
                    offset = (*pA0)(n - 1, n - 1);
                }

                mT = *pA0;
                return i;
            }

            /**
             * @brief 最基础的 QR 迭代实现, 用于对比带位移的 QR 算法
             * 
             * @param [in] a 目标矩阵
             * @param [in] max_iter 最大迭代次数
             * @param [in] tolerance 结束迭代的容忍度
             * @return 迭代次数
             */
            int NaiveIterate(MatViewIn const & a, int max_iter, Scalar tolerance)
            {
                assert(a.Rows() == a.Cols());
                mT.Resize(a.Rows(), a.Cols());

                DMatrix<Scalar> tmp0 = a;
                DMatrix<Scalar> tmp1 = DMatrix<Scalar>::Zero(a.Rows(), a.Cols());
                DMatrix<Scalar> * pA0 = &tmp0;
                DMatrix<Scalar> * pA1 = &tmp1;
                QR_Householder<DMatrix<Scalar>> qr;

                int i = 0;
                for (; i < max_iter; i++) {
                    qr.Decompose(*pA0);
                    (*pA1) = qr.R() * qr.Q();
                    std::swap(pA0, pA1);

                    if (IsConverged(*pA0, *pA1, tolerance))
                        break;
                }

                mT = *pA0;
                return i;
            }

        private:

            /**
             * @brief 判定算法是否收敛
             */
            bool IsConverged(DMatrix<Scalar> const & A0, DMatrix<Scalar> const & A1, Scalar tolerance)
            {
                int num = A0.NumDatas();
                for (int i = 0; i < num; ++i) {
                    Scalar diff = std::abs(A0(i) - A1(i));
                    if (diff >= tolerance)
                        return false;
                }
                return true;
            }

        public:
            DMatrix<Scalar> const & T() { return mT; }

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
    };

}

#endif

