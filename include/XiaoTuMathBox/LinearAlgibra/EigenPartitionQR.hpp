#ifndef XTMB_LA_EIGEN_IMPLICIT_QR_H
#define XTMB_LA_EIGEN_IMPLICIT_QR_H

namespace xiaotu::math {


    /**
     * @brief 对角线分块 QR 算法, 出于教学目的编写的, 效率较低
     */
    template <typename MatViewIn>
    class EigenPartitionQR {
        public:
            typedef typename MatViewIn::Scalar Scalar;
            typedef DMatrix<Scalar> Mat;

            /**
             * @brief 默认构造函数
             */
            EigenPartitionQR()
            {
            }
            /**
             * @brief 带位移的 QR 迭代
             * 
             * @param [in] a 目标矩阵
             * @param [in] max_iter 最大迭代次数
             * @param [in] tolerance 结束迭代的容忍度
             * @param [in] first_off 若非 nullptr，将指定第一次迭代时的偏移量
             * @return 迭代次数
             */
            int Iterate(MatViewIn const & a, int max_iter, Scalar tolerance, Scalar * first_off = nullptr)
            {
                assert(a.Rows() == a.Cols());
                int n = a.Rows();

                UpperHessenberg h(a);
                DMatrix<Scalar> tmp0 = h.H();
                DMatrix<Scalar> tmp1 = DMatrix<Scalar>::Zero(n, n);
                // std::cout << "Hessenberg:" << h.H() << std::endl;

                std::vector<MatrixSubView<Mat>> parts0;
                std::vector<MatrixSubView<Mat>> parts1;
                std::vector<MatrixSubView<Mat>> parts2;
                std::vector<MatrixSubView<Mat>> parts3;
                std::vector<MatrixSubView<Mat>> * pParts0 = &parts0;
                std::vector<MatrixSubView<Mat>> * pParts1 = &parts1;
                std::vector<MatrixSubView<Mat>> * pParts2 = &parts2;
                std::vector<MatrixSubView<Mat>> * pParts3 = &parts3;
                pParts0->push_back(tmp0.SubMatrix(0, 0, n, n));
                pParts1->push_back(tmp1.SubMatrix(0, 0, n, n));

                int i = 0;
                QR_Householder<DMatrix<Scalar>> qr;
                for (; i < max_iter; i++) {
                    if (pParts0->empty())
                        break;
                    // std::cout << "----------" << std::endl;
                    for (int idx = 0; idx < pParts0->size(); idx++) {
                        auto & a0 = (*pParts0)[idx];
                        auto & a1 = (*pParts1)[idx];
                        int n = a0.Rows();
                        Scalar offset = (0 == i && nullptr != first_off)
                                      ? *first_off
                                      : a0(n - 1, n - 1);

                        // std::cout << "-- " << idx << ":" << n << a0;

                        SubDiagScalar(a0, offset);
                        qr.Decompose(a0);

                        a1 = qr.R() * qr.Q();

                        AddDiagScalar(a1, offset);
                        int np = Partition(a0, a1, *pParts2, *pParts3);
                        // std::cout << "np:" << np << std::endl;
                    }

                    std::swap(pParts0, pParts3);
                    std::swap(pParts1, pParts2);
                    pParts2->clear();
                    pParts3->clear();
                }

                return i;
            }


            /**
             * @brief 从矩阵 mT 中获取对角元素，即特征值
             */
            std::vector<Scalar> const & EigenValues() const
            {
                return eigens_;
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
             * @brief 参考 a1, 按照下次对角线是否为 0 对 a0/a1 进行分割
             * 
             * @param [in] 分割对象 a0
             * @param [in] 分割对象 a1
             * @param [out] 输出分割列表 part3 对应 a0
             * @param [out] 输出分割列表 part4 对应 a1
             * @return 分割的对角块数量
             */
            int Partition(MatrixSubView<Mat> & a0,
                          MatrixSubView<Mat> & a1,
                          std::vector<MatrixSubView<Mat>> & parts2,
                          std::vector<MatrixSubView<Mat>> & parts3)
            {
                int n = a1.Rows();
                int start = 0;

                for (int i = 0; i < (n-1); i++) {
                    if (std::abs(a1(i+1, i)) < SMALL_VALUE) {
                        int m = i + 1 - start;
                        if (1 == m) {
                            eigens_.push_back(a1(start, start));
                        } else if (2 == m) {
                            Solve2x2(a1.SubMatrix(start, start, m, m));
                        } else {
                            parts2.push_back(a0.SubMatrix(start, start, m, m));
                            parts3.push_back(a1.SubMatrix(start, start, m, m));
                        }
                        start = i+1;
                    }
                }
                int m = n - start;
                if (1 == m) {
                    eigens_.push_back(a1(start, start));
                } else if (2 == m) {
                    Solve2x2(a1.SubMatrix(start, start, m, m));
                } else {
                    parts2.push_back(a0.SubMatrix(start, start, m, m));
                    parts3.push_back(a1.SubMatrix(start, start, m, m));
                }
                return parts3.size();
            }


        private:
            std::vector<Scalar> eigens_;
    };

}


#endif
