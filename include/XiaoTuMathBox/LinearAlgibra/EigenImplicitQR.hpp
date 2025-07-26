#ifndef XTMB_LA_EIGEN_IMPLICIT_QR_H
#define XTMB_LA_EIGEN_IMPLICIT_QR_H

namespace xiaotu::math {


    /**
     * @brief 隐式 QR 迭代, 带偏移，对角分块
     *
     * QAQ^T = Sigma, A = Q^T Sigma Q
     */
    template <typename MatViewIn>
    class EigenImplicitQR {
        public:
            typedef typename MatViewIn::Scalar Scalar;
            typedef DMatrix<Scalar> Mat;

            /**
             * @brief 默认构造函数
             */
            EigenImplicitQR()
            {
            }
            /**
             * @brief 带位移、对角分块、隐式 QR 迭代
             * 
             * @param [in] a 目标矩阵
             * @param [in] max_iter 最大迭代次数
             * @param [in] tolerance 结束迭代的容忍度
             * @param [in] first_off 若非 nullptr，将指定第一次迭代时的偏移量
             * @param [in] keep_q 是否需要保留正交矩阵 Q, QAQ^T = Sigma, A = Q^T Sigma Q
             * @return 迭代次数
             */
            int Iterate(MatViewIn const & a, int max_iter, Scalar tolerance,
                        Scalar * first_off = nullptr,
                        bool keep_q = true)
            {
                assert(a.Rows() == a.Cols());
                int n = a.Rows();

                mSigma.Resize(n, n) = a;
                mQ.Resize(n, n).Identity();
                UpperHessenbergByHouseholder(mSigma, mQ);

                std::vector<MatrixSubView<Mat>> parts0;
                std::vector<MatrixSubView<Mat>> parts1;
                std::vector<MatrixSubView<Mat>> * pParts0 = &parts0;
                std::vector<MatrixSubView<Mat>> * pParts1 = &parts1;
                pParts0->push_back(mSigma.SubMatrix(0, 0, n, n));


                std::vector<MatrixSubView<Mat>> partQ0;
                std::vector<MatrixSubView<Mat>> partQ1;
                std::vector<MatrixSubView<Mat>> * pPartQ0 = &partQ0;
                std::vector<MatrixSubView<Mat>> * pPartQ1 = &partQ1;
                pPartQ0->push_back(mQ.SubMatrix(0, 0, n, n));

                int i = 0;
                for (; i < max_iter; i++) {
                    if (pParts0->empty())
                        break;
                    for (int idx = 0; idx < pParts0->size(); idx++) {
                        auto & a0 = (*pParts0)[idx];
                        int n = a0.Rows();
                        Scalar offset = (0 == i && nullptr != first_off)
                                      ? *first_off
                                      : a0(n - 1, n - 1);

                        SubDiagScalar(a0, offset);
                        MatrixSubView<Mat> * pQ0 = keep_q
                                                 ? &(*pPartQ0)[idx]
                                                 : nullptr;
                        ImplicitQR(a0, pQ0);

                        AddDiagScalar(a0, offset);
                        Partition(a0, *pParts1, pQ0, pPartQ1);
                    }

                    std::swap(pParts0, pParts1);
                    std::swap(pPartQ0, pPartQ1);
                    pParts1->clear();
                    pPartQ1->clear();
                }

                return i;
            }


        private:

            /**
             * @brief 直接写出 2x2 矩阵的特征值
             *
             * @return true - 复特征值, false 实特征值
             */
            bool Solve2x2(MatrixSubView<Mat> const & a1)
            {
                assert(a1.Rows() == 2);
                Scalar tr = a1(0, 0) + a1(1, 1);
                Scalar det = a1(0, 0) * a1(1,1) - a1(0,1) * a1(1,0);
                Scalar k = tr * tr - 4 * det;
                if (k < 0) {
                    std::cout << "存在复特征值" << std::endl;
                    Scalar sk = std::sqrt(-k);
                    std::cout << 0.5 * tr << " + " << 0.5 * sk << "i" << std::endl;
                    std::cout << 0.5 * tr << " - " << 0.5 * sk << "i" << std::endl;
                    return true;
                } else {
                    Scalar sk = std::sqrt(k);
                    eigens_.push_back(0.5 * (tr + sk));
                    eigens_.push_back(0.5 * (tr - sk));
                    return false;
                }
            }

            /**
             * @brief 按照下次对角线是否为 0 对 a0 进行分割
             * 
             * @param [in] a0 分割对象 
             * @param [out] parts 输出分割列表
             * @return 分割的对角块数量
             */
            int Partition(MatrixSubView<Mat> & a0,
                          std::vector<MatrixSubView<Mat>> & parts,
                          MatrixSubView<Mat> * pQ0,
                          std::vector<MatrixSubView<Mat>> * partQ)
            {
                int n = a0.Rows();
                int start = 0;

                for (int i = 0; i < (n-1); i++) {
                    if (std::abs(a0(i+1, i)) < SMALL_VALUE) {
                        int m = i + 1 - start;
                        bool split = m > 2;
                        if (2 == m)
                            split = !Solve2x2(a0.SubMatrix(start, start, m, m));

                        if (split) {
                            parts.push_back(a0.SubMatrix(start, start, m, m));
                            if (nullptr != pQ0 && nullptr != partQ)
                                partQ->push_back(pQ0->SubMatrix(start, 0, m, pQ0->Cols()));
                        }
                        start = i+1;
                    }
                }
                {
                    int m = n - start;
                    bool split = m > 2;
                    if (2 == m)
                        split = !Solve2x2(a0.SubMatrix(start, start, m, m));

                    if (split) {
                        parts.push_back(a0.SubMatrix(start, start, m, m));
                        if (nullptr != pQ0 && nullptr != partQ)
                            partQ->push_back(pQ0->SubMatrix(start, 0, m, pQ0->Cols()));
                    }
                }
                return parts.size();
            }


            /**
             * @brief 对子阵 H 进行隐式 QR 迭代
             */
            void ImplicitQR(MatrixSubView<Mat> & H, MatrixSubView<Mat> * pQ)
            {
                Givens<Scalar> G(0, 1, H(0, 0), H(1, 0));
                G.LeftApplyOn(H);
                G.TRightApplyOn(H);
                if (nullptr != pQ)
                    G.LeftApplyOn(*pQ);

                int n = H.Cols() - 2;
                for (int i = 0; i < n; i++) {
                    Scalar bugle = H(i + 2, i);
                    if (std::abs(bugle) < SMALL_VALUE)
                        break;

                    Givens<Scalar> G(i+1, i+2, H(i+1, i), H(i+2, i));
                    G.LeftApplyOn(H);
                    G.TRightApplyOn(H);
                    if (nullptr != pQ)
                        G.LeftApplyOn(*pQ);
                }
            }

        public:
            DMatrix<Scalar> const & Sigma() { return mSigma; }
            DMatrix<Scalar> const & Q() { return mQ; }

            /**
             * @brief 从矩阵 mT 中获取对角元素，即特征值
             */
            std::vector<Scalar> EigenValues() const
            {
                int n = mSigma.Rows();
                std::vector<Scalar> re(n);

                for (int i = 0; i < n; i++)
                    re[i] = mSigma(i, i);

                return re;
            }

        private:
            std::vector<Scalar> eigens_;
            DMatrix<Scalar> mSigma;
            DMatrix<Scalar> mQ;
    };

}


#endif
