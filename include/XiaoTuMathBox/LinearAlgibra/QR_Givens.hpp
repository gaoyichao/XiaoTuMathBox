#ifndef XTMB_LA_QR_GIVENS_H
#define XTMB_LA_QR_GIVENS_H

namespace xiaotu::math {

    //! @brief 基于 Givens 的 QR 分解
    //!
    //! 将输入的矩阵，分解为一个正交矩阵 Q 和一个上三角矩阵 R
    template <typename MatViewIn>
    class QR_Givens {
        public:
            typedef typename MatViewIn::Scalar Scalar;

            /**
             * @brief 默认构造函数
             */
            QR_Givens()
            {}

            /**
             * @brief 构造函数
             * 
             * @param [in] a 待分解矩阵
             */
            QR_Givens(MatViewIn const & a)
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
                int rows = mR.Rows();
                int cols = mR.Cols();
                
                for (int cidx = 0; cidx < cols; cidx++) {
                    for (int ridx = cidx+1; ridx < rows; ridx++) {
                        if (std::abs(mR(ridx, cidx)) < SMALL_VALUE)
                            continue;
                        Givens<Scalar> G(cidx, ridx, mR(cidx, cidx), mR(ridx, cidx));
                        G.LeftApplyOn(mR);
                        G.TRightApplyOn(mQ);
                    }
                }
                
                // 为了保证 QR 分解唯一，添加约束 R 的对角元素都是正数
                // https://gaoyichao.com/Xiaotu/?book=algebra&title=QR算法的收敛原理
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

