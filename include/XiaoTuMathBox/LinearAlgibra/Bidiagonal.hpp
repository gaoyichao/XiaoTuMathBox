#ifndef XTMB_LA_BI_DIAGONAL_H
#define XTMB_LA_BI_DIAGONAL_H

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
                    mB.UpperBidiagonal(&mUT, &mV);
                else
                    mB.LowerBidiagonal(&mUT, &mV);
            }

        public:
            //! @brief B = U^T A V. 输出的二对角矩阵
            DMatrix<Scalar> const & B() { return mB; }
            //! @brief B = U^T A V
            DMatrix<Scalar> const & UT() { return mUT; }
            //! @brief B = U^T A V
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
