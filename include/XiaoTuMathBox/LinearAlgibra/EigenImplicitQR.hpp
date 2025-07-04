#ifndef XTMB_LA_EIGEN_IMPLICIT_QR_H
#define XTMB_LA_EIGEN_IMPLICIT_QR_H

namespace xiaotu::math {


    /**
     * @brief 双位移-对角块-隐 QR 算法
     */
    template <typename MatViewIn>
    class EigenImplicitQR {
        public:
            typedef typename MatViewIn::Scalar Scalar;

            /**
             * @brief 默认构造函数
             */
            EigenImplicitQR()
            {
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
                mT = *pA0;
                return i;
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
