#ifndef XTMB_LA_SVD_GKR_H
#define XTMB_LA_SVD_GKR_H

namespace xiaotu {

    /**
     * @brief Golub-Kahan-Reinsch 算法,进行 SVD 分解
     * 
     * \Sigma = U^T A V
     */
    template <typename MatViewIn>
    class SVD_GKR {
        public:
            typedef typename MatViewIn::Scalar Scalar;
            typedef DMatrix<Scalar> Mat;
            typedef MatrixSubView<Mat> SubView;

            /**
             * @brief 默认构造函数
             */
            SVD_GKR(MatViewIn const & a, bool keepU, bool keepV)
                : mKeepU(keepU), mKeepV(keepV)
            {
                int m = a.Rows();
                int n = a.Cols();
                int p = (m < n) ? m : n;

                mSigma.Resize(m, n) = a;
                mSigmaPingPang.Init(std::move(mSigma.SubMatrix(0, 0, p, p)));
                if (keepU) {
                    mUT.Resize(m, m).Identity();
                    mUTPingPang.Init(std::move(mUT.SubMatrix(0, 0, p, m)));
                }
                if (keepV) {
                    mV.Resize(n, n).Identity();
                    mVPingPang.Init(std::move(mV.SubMatrix(0, 0, n, p)));
                }
                mSigma.Bidiagonal(keepU ? &mUT : nullptr, keepV ? &mV  : nullptr);

                {
                    auto & sigma = mSigmaPingPang.Ping()[0];
                    SubView * pu = mKeepU ? &(mUTPingPang.Ping()[0]) : nullptr;
                    SubView * pv = mKeepV ? &(mVPingPang.Ping()[0]) : nullptr;

                    if (mSigma.Rows() >= mSigma.Cols()) {
                        UpperScanDiagonal(sigma, pu, pv);
                        Partition(sigma, pu, pv, true);
                    }
                    else {
                        LowerScanDiagonal(sigma, pu, pv);
                        Partition(sigma, pu, pv, false);
                    }
                    PingPang();
                }

            }

            int Iterate(int max_iter, Scalar tolerance)
            {
                mAbsEps = tolerance;
                int num = 0;

                if (mSigma.Rows() >= mSigma.Cols())
                    num = UpperIterate(max_iter, tolerance);
                else
                    num = LowerIterate(max_iter, tolerance);

                Sort();
                return num;
            }

        private:
            /**
             * @brief 对奇异值从大到小降序排列
             */
            void Sort()
            {
                int m = mSigma.Rows();
                int n = mSigma.Cols();
                int p = (m < n) ? m : n;

                std::vector<int> idx(p);
                for (int i = 0; i < p; i++)
                    idx[i] = i;

                std::sort(idx.begin(), idx.end(), [this](int i, int j){
                    return std::abs(mSigma(i,i)) > std::abs(mSigma(j,j));
                });


                auto sigma = mSigma.SubMatrix(0, 0, p, p);
                auto ut = mUT.SubMatrix(0, 0, p, m);
                auto v = mV.SubMatrix(0, 0, n, p);

                Permutation perm(idx);
                perm.LeftApplyOn(sigma);
                if (mKeepU)
                    perm.LeftApplyOn(ut);
                perm.RightApplyOn(sigma);
                if (mKeepV)
                    perm.RightApplyOn(v);
            }

            /**
             * @brief 上二对角阵的迭代
             */
            int UpperIterate(int max_iter, Scalar tolerance)
            {
                int i = 0;
                for (; i < max_iter; i++) {
                    auto & sigma_list = mSigmaPingPang.Ping();
                    if (sigma_list.empty())
                        break;

                    for (int idx = 0; idx < sigma_list.size(); idx++) {
                        auto & a0 = sigma_list[idx];
                        SubView * pu = mKeepU ? &(mUTPingPang.Ping()[idx]) : nullptr;
                        SubView * pv = mKeepV ? &(mVPingPang.Ping()[idx]) : nullptr;

                        int p = a0.Rows();
                        if (UpperIsAnnZero(a0, pv)) {
                            p = p - 1;
                            a0.Shrink(0, 0, p, p);
                            if (nullptr != pu)
                                pu->Shrink(0, 0, p, mSigma.Cols());
                            if (nullptr != pv)
                                pv->Shrink(0, 0, mSigma.Rows(), p);
                        }
                        if (p < 2)
                            continue;

                        UpperImplicitIter(a0, pu, pv);

                        UpperScanDiagonal(a0, pu, pv);
                        Partition(a0, pu, pv, true);
                    }

                    PingPang();
                }


                return i;
            }

            /**
             * @brief 扫描对角线，消除对角线为 0 的那一行
             * 
             * @param [in|out] sigma 对角矩阵块
             * @param [out] u 若非 nullptr 则记录变换过程，左乘 Givens 变换
             */
            void UpperScanDiagonal(MatrixSubView<Mat> & sigma,
                              MatrixSubView<Mat> * u,
                              MatrixSubView<Mat> * v)
            {
                int p = sigma.Rows();
                for (int i = 0; i < p; i++) {
                    if (std::abs(sigma(i, i)) > mAbsEps)
                        continue;

                    for (int j = i+1; j < p; j++) {
                        if (std::abs(sigma(i, j)) < mAbsEps)
                            break;

                        Givens<Scalar> G(j, i, sigma(j, j), sigma(i, j));
                        G.LeftApplyOn(sigma);
                        if (nullptr != u)
                            G.LeftApplyOn(*u);
                    }
                }
            }

            /**
             * @brief 检查右下角的元素，若为 0 则消除最后一列
             * 
             * @param [in|out] sigma 对角矩阵块
             * @param [out] v 若非 nullptr 则记录列消解过程, 右乘 Givens 变换
             * @return 是否执行消解操作
             */
            bool UpperIsAnnZero(MatrixSubView<Mat> & sigma, MatrixSubView<Mat> * v)
            {
                int p = sigma.Rows();
                if (std::abs(sigma(p-1, p-1)) > mAbsEps)
                    return false;

                p = p - 1;
                for (int i = p-1; i >= 0; i--) {
                    Givens<Scalar> G(p, i, sigma(i, i), sigma(i, p));
                    G.RightApplyOn(sigma);
                    if (nullptr != v)
                        G.RightApplyOn(*v);
                }
                return true;
            }

            /**
             * @brief 对上二对角矩阵 a 计算 Rayleigh 偏移
             */
            Scalar UpperRayleighShift(MatrixSubView<Mat> const & a)
            {
                int p = a.Rows() - 1;
                return a(p, p) * a(p, p) + a(p-1, p) * a(p-1, p);
            }
            /**
             * @brief 对上二对角矩阵 a 计算 Wilkinson 偏移
             *
             * 相比于 RayleighShift 收敛性更好一些
             */
            Scalar UpperWilkinsonShift(MatrixSubView<Mat> const & a)
            {
                int p = a.Rows() - 1;
                Scalar d_n1 = a(p-1, p-1);
                Scalar d_n = a(p, p);
                Scalar f_n2 = a(p-2, p-1);
                Scalar f_n1 = a(p-1, p);

                Scalar t11 = d_n1 * d_n1 + f_n2 * f_n2;
                Scalar t12 = d_n1 * f_n1;
                Scalar t21 = f_n1 * d_n1;
                Scalar t22 = d_n * d_n + f_n1 * f_n1;

                Scalar tr = t11 + t22;
                Scalar det = t11 * t22 - t12 * t21;
                Scalar sk = sqrt(tr * tr - 4 * det);

                Scalar e1 = 0.5 * (tr + sk);
                Scalar e2 = 0.5 * (tr - sk);

                Scalar diff1 = std::abs(e1 - d_n);
                Scalar diff2 = std::abs(e2 - d_n);
                return diff1 < diff2 ? e1 : e2;
            }

            void UpperImplicitIter(MatrixSubView<Mat> & B,
                              MatrixSubView<Mat> * pU,
                              MatrixSubView<Mat> * pV)
            {
                int n = B.Rows();
                //Scalar mu = UpperRayleighShift(B);
                Scalar mu = UpperWilkinsonShift(B);
                Scalar y = B(0, 0) * B(0, 0) - mu;
                Scalar z = B(0, 0) * B(0, 1);

                Givens<Scalar> G1(1, 0, y, z);
                G1.RightApplyOn(B);
                if (nullptr != pV)
                    G1.RightApplyOn(*pV);


                for (int k = 0; k < n-1; k++) {
                    y = B(k,   k);
                    z = B(k+1, k);
                    Givens<Scalar> G(k, k+1, y, z);
                    G.LeftApplyOn(B);
                    if (nullptr != pU)
                        G.LeftApplyOn(*pU);

                    if (k < n-2) {
                        y = B(k, k+1);
                        z = B(k, k+2);

                        Givens<Scalar> G1(k+2, k+1, y, z);
                        G1.RightApplyOn(B);
                        if (nullptr != pV)
                            G1.RightApplyOn(*pV);
                    }
                }
            }


        private:
            /**
             * @brief 下二对角阵的迭代
             */
            int LowerIterate(int max_iter, Scalar tolerance)
            {
                int i = 0;
                for (; i < max_iter; i++) {
                    auto & sigma_list = mSigmaPingPang.Ping();
                    if (sigma_list.empty())
                        break;

                    for (int idx = 0; idx < sigma_list.size(); idx++) {
                        auto & a0 = sigma_list[idx];
                        SubView * pu = mKeepU ? &(mUTPingPang.Ping()[idx]) : nullptr;
                        SubView * pv = mKeepV ? &(mVPingPang.Ping()[idx]) : nullptr;

                        int p = a0.Rows();
                        if (LowerIsAnnZero(a0, pu)) {
                            p = p - 1;
                            a0.Shrink(0, 0, p, p);
                            if (nullptr != pu)
                                pu->Shrink(0, 0, p, mSigma.Cols());
                            if (nullptr != pv)
                                pv->Shrink(0, 0, mSigma.Rows(), p);
                        }
                        if (p < 2)
                            continue;

                        LowerImplicitIter(a0, pu, pv);
                        LowerScanDiagonal(a0, pu, pv);
                        Partition(a0, pu, pv, false);
                    }

                    PingPang();
                }


                return i;
            }

            /**
             * @brief 检查右下角的元素，若为 0 则消除最后一列
             * 
             * @param [in|out] sigma 对角矩阵块
             * @param [out] u 若非 nullptr 则记录列消解过程, 左乘 Givens 变换
             * @return 是否执行消解操作
             */
            bool LowerIsAnnZero(MatrixSubView<Mat> & sigma, MatrixSubView<Mat> * u)
            {
                int p = sigma.Rows();
                if (std::abs(sigma(p-1, p-1)) > mAbsEps)
                    return false;

                p = p - 1;
                for (int i = p-1; i >= 0; i--) {
                    Givens<Scalar> G(i, p, sigma(i, i), sigma(p, i));
                    G.LeftApplyOn(sigma);
                    if (nullptr != u)
                        G.LeftApplyOn(*u);
                }
                return true;
            }

            /**
             * @brief 扫描对角线，消除对角线为 0 的那一列
             * 
             * @param [in|out] sigma 对角矩阵块
             * @param [out] u 若非 nullptr 则记录变换过程，右乘 Givens 变换
             */
            void LowerScanDiagonal(MatrixSubView<Mat> & sigma,
                              MatrixSubView<Mat> * u,
                              MatrixSubView<Mat> * v)
            {
                int p = sigma.Rows();
                for (int i = 0; i < p; i++) {
                    if (std::abs(sigma(i, i)) > mAbsEps)
                        continue;

                    for (int j = i+1; j < p; j++) {
                        if (std::abs(sigma(j, i)) < mAbsEps)
                            break;

                        Givens<Scalar> G(i, j, sigma(j, j), sigma(j, i));
                        G.RightApplyOn(sigma);
                        if (nullptr != v)
                            G.RightApplyOn(*v);
                    }
                }
            }

            /**
             * @brief 对下二对角矩阵 a 计算 Wilkinson 偏移
             *
             * 相比于 RayleighShift 收敛性更好一些
             */
            Scalar LowerWilkinsonShift(MatrixSubView<Mat> const & a)
            {
                int p = a.Rows() - 1;
                Scalar d_n1 = a(p-1, p-1);
                Scalar d_n = a(p, p);
                Scalar f_n2 = a(p-1, p-2);
                Scalar f_n1 = a(p, p-1);

                Scalar t11 = d_n1 * d_n1 + f_n2 * f_n2;
                Scalar t12 = d_n1 * f_n1;
                Scalar t21 = f_n1 * d_n1;
                Scalar t22 = d_n * d_n + f_n1 * f_n1;

                Scalar tr = t11 + t22;
                Scalar det = t11 * t22 - t12 * t21;
                Scalar sk = sqrt(tr * tr - 4 * det);

                Scalar e1 = 0.5 * (tr + sk);
                Scalar e2 = 0.5 * (tr - sk);

                Scalar diff1 = std::abs(e1 - d_n);
                Scalar diff2 = std::abs(e2 - d_n);
                return diff1 < diff2 ? e1 : e2;
            }

            void LowerImplicitIter(MatrixSubView<Mat> & B,
                              MatrixSubView<Mat> * pU,
                              MatrixSubView<Mat> * pV)
            {
                int n = B.Rows();
                Scalar mu = LowerWilkinsonShift(B);
                Scalar y = B(0, 0) * B(0, 0) - mu;
                Scalar z = B(0, 0) * B(1, 0);

                Givens<Scalar> G1(0, 1, y, z);
                G1.LeftApplyOn(B);
                if (nullptr != pU)
                    G1.LeftApplyOn(*pU);


                for (int k = 0; k < n-1; k++) {
                    y = B(k, k);
                    z = B(k, k+1);
                    Givens<Scalar> G(k+1, k, y, z);
                    G.RightApplyOn(B);
                    if (nullptr != pV)
                        G.RightApplyOn(*pV);

                    if (k < n-2) {
                        y = B(k+1, k);
                        z = B(k+2, k);

                        Givens<Scalar> G1(k+1, k+2, y, z);
                        G1.LeftApplyOn(B);
                        if (nullptr != pU)
                            G1.LeftApplyOn(*pU);
                    }
                
                }
            }

        public:
            DMatrix<Scalar> const Sigma() { return mSigma; }
            DMatrix<Scalar> const UT() { return mUT; }
            DMatrix<Scalar> const V() { return mV; }

        private:
            DMatrix<Scalar> mUT;
            DMatrix<Scalar> mSigma;
            DMatrix<Scalar> mV;

        public:
            //! @brief 绝对精度
            Scalar mAbsEps{SMALL_VALUE};
            //! @brief 是否需要保留 UT
            bool mKeepU{false};
            //! @brief 是否需要保留 V
            bool mKeepV{false};

        private:

            /**
             * @brief 对二对角矩阵按照次对角线元素为 0 对 a0 进行分割
             * 
             * @param [in] a0 分割对象 
             * @param [out] parts 输出分割列表
             * @return 分割的对角块数量
             */
            int Partition(MatrixSubView<Mat> & sigma,
                          MatrixSubView<Mat> * pU0,
                          MatrixSubView<Mat> * pV0,
                          bool upper)
            {
                int p = sigma.Rows();
                int m = mSigma.Rows();
                int n = mSigma.Cols();
                int start = 0;

                auto & partS = mSigmaPingPang.Pang();
                auto & partU = mUTPingPang.Pang();
                auto & partV = mVPingPang.Pang();
                for (int i = 0; i < (p-1); i++) {
                    bool need_part = (upper && std::abs(sigma(i, i+1)) < mAbsEps) ||
                                     (!upper && std::abs(sigma(i+1, i)) < mAbsEps);
                    if (need_part) {
                        int size = i + 1 - start;
                        if (size > 1) {
                            partS.push_back(sigma.SubMatrix(start, start, size, size));
                            if (nullptr != pU0)
                                partU.push_back(pU0->SubMatrix(start, 0, size, m));
                            if (nullptr != pV0)
                                partV.push_back(pV0->SubMatrix(0, start, n, size));
                        }
                        start = i+1;
                    }
                }
                {
                    int size = p - start;
                    if (size > 1) {
                        partS.push_back(sigma.SubMatrix(start, start, size, size));
                        if (nullptr != pU0)
                            partU.push_back(pU0->SubMatrix(start, 0, size, m));
                        if (nullptr != pV0)
                            partV.push_back(pV0->SubMatrix(0, start, n, size));
                    }
                }
                return partS.size();
            }

            /**
             * @brief 交换乒乓队列
             */
            void PingPang()
            {
                mSigmaPingPang.Swap();
                if (mKeepU)
                    mUTPingPang.Swap();
                if (mKeepV)
                    mVPingPang.Swap();
            }

            //! @brief mUT 的乒乓队列
            PingPangView<Mat> mUTPingPang;
            //! @brief mSigma 的乒乓队列
            PingPangView<Mat> mSigmaPingPang;
            //! @brief mV 的乒乓队列
            PingPangView<Mat> mVPingPang;
    };

}

#endif


