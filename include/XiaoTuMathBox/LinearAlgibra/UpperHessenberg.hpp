#ifndef XTMB_LA_UPPER_HESSENBERG_H
#define XTMB_LA_UPPER_HESSENBERG_H

namespace xiaotu {


    //! @brief 将输入矩阵，相似变换成上 Hessenberg 矩阵
    //!
    //! 出于教学的目的而保留
    //! https://gaoyichao.com/Xiaotu/?book=algebra&title=加快QR算法的收敛过程
    //!
    //! | a_11 a_12 a_13 a_14 |    | a_11 a_12 a_13 a_14 |
    //! | a_21 a_22 a_23 a_24 |    | a_21 a_22 a_23 a_24 |
    //! | a_31 a_32 a_33 a_34 | => | 0    a_32 a_33 a_34 |
    //! | a_41 a_42 a_43 a_44 |    | 0    0    a_43 a_44 |
    //!
    //! A = Q^T * H * Q
    template <typename MatViewIn>
    class UpperHessenberg {
        public:
            typedef typename MatViewIn::Scalar Scalar;

        public:

            /**
             * @brief 默认构造函数
             */
            UpperHessenberg()
            {}


            /**
             * @brief 构造函数
             * 
             * @param [in] a 待分解矩阵
             */
            UpperHessenberg(MatViewIn const & a)
                : mQ(a.Rows(), a.Rows()), mH(a.Rows(), a.Cols())
            {
                mQ.Identity();
                mH = a;
                mH.UpperHessenbergByHouseholder(&mQ);
            }

            /**
             * @brief 通过 Householder 完成分解
             * 
             * @param [in] a 待分解矩阵
             */
            void ByHouseholder(MatViewIn const & a)
            {
                if (mH.Rows() != a.Rows() || mH.Cols() != a.Cols()) {
                    mQ.Resize(a.Rows(), a.Rows());
                    mH.Resize(a.Rows(), a.Cols());
                }
                mQ.Identity();
                mH = a;
                mH.UpperHessenbergByHouseholder(&mQ);
            }

            //! @todo 通过 Givens 完成分解

        public:
            DMatrix<Scalar> const & Q() { return mQ; }
            DMatrix<Scalar> const & H() { return mH; }

        private:
            DMatrix<Scalar> mQ;
            DMatrix<Scalar> mH;
    };
}

#endif
