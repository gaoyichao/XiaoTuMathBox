#ifndef XTMB_LA_SVD_GKR_H
#define XTMB_LA_SVD_GKR_H

#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

namespace xiaotu::math {

    /**
     * @brief Golub-Kahan-Reinsch 算法,进行 SVD 分解
     * 
     * \Sigma = U^T A V
     * 
     * 只处理瘦高形矩阵(m >= n)
     */
    template <typename MatViewIn>
    class SVD_UpperGKR {
        public:
            typedef typename MatViewIn::Scalar Scalar;
            typedef DMatrix<Scalar> Mat;
            typedef MatrixSubView<Mat> SubView;

            /**
             * @brief 默认构造函数
             */
            SVD_UpperGKR(MatViewIn const & a, bool keepU, bool keepV)
                : mKeepU(keepU), mKeepV(keepV)
            {
                assert(a.Rows() >= a.Cols());

                int m = a.Rows();
                int n = a.Cols();

                mSigma.Resize(m, n) = a;
                mUT.Resize(m, m).Identity();
                mV.Resize(n, n).Identity();

                mSigmaPingPang.Init(std::move(mSigma.SubMatrix(0, 0, n, n)));
                mUTPingPang.Init(std::move(mUT.SubMatrix(0, 0, n, m)));
                mVPingPang.Init(std::move(mV.SubMatrix(0, 0, n, n)));

                mSigma.UpperBidiagonal(&mUT, &mV);
            }

            int Iterate(int max_iter, Scalar tolerance)
            {
                mAbsEps = tolerance;
                int m = mSigma.Rows();
                int n = mSigma.Cols();
                int p = m < n ? m : n;

                {
                    auto & sigma = mSigmaPingPang.Ping()[0];
                    SubView * pu = mKeepU ? &(mUTPingPang.Ping()[0]) : nullptr;
                    SubView * pv = mKeepV ? &(mVPingPang.Ping()[0]) : nullptr;

                    ScanDiagonal(sigma, pu, pv);
                    UpperPartition(sigma, pu, pv);
                    PingPang();
                }

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
                        if (IsAnnZero(a0, pv)) {
                            p = p - 1;
                            a0.Shrink(0, 0, p, p);
                            if (nullptr != pu)
                                pu->Shrink(0, 0, p, mSigma.Cols());
                            if (nullptr != pv)
                                pv->Shrink(0, 0, mSigma.Rows(), p);
                        }
                        if (p < 2)
                            continue;

                        ImplicitIter(a0, pu, pv);

                        ScanDiagonal(a0, pu, pv);
                        UpperPartition(a0, pu, pv);
                        PingPang();
                    }
                }


                return i;
            }

        private:

            /**
             * @brief 扫描对角线，消除对角线为 0 的那一行
             * 
             * @param [in|out] sigma 对角矩阵块
             * @param [out] u 若非 nullptr 则记录变换过程，左乘 Givens 变换
             */
            void ScanDiagonal(MatrixSubView<Mat> & sigma,
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
            bool IsAnnZero(MatrixSubView<Mat> & sigma, MatrixSubView<Mat> * v)
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
             * @brief 对上二对角矩阵按照次对角线元素为 0 对 a0 进行分割
             * 
             * @param [in] a0 分割对象 
             * @param [out] parts 输出分割列表
             * @return 分割的对角块数量
             */
            int UpperPartition(MatrixSubView<Mat> & sigma,
                               MatrixSubView<Mat> * pU0,
                               MatrixSubView<Mat> * pV0)
            {
                int p = sigma.Rows();
                int m = mSigma.Rows();
                int n = mSigma.Cols();
                int start = 0;

                auto & partS = mSigmaPingPang.Pang();
                auto & partU = mUTPingPang.Pang();
                auto & partV = mVPingPang.Pang();
                for (int i = 0; i < (p-1); i++) {
                    if (std::abs(sigma(i, i+1)) < mAbsEps) {
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
             * @brief 交换乒乓对接
             */
            void PingPang()
            {
                mSigmaPingPang.Swap();
                if (mKeepU)
                    mUTPingPang.Swap();
                if (mKeepV)
                    mVPingPang.Swap();
            }

            /**
             * @brief 对上二对角矩阵 a 计算 Rayleigh 偏移
             */
            Scalar RayleighShift(MatrixSubView<Mat> const & a)
            {
                int p = a.Rows() - 1;
                return a(p, p) * a(p, p) + a(p-1, p) * a(p-1, p);
            }

            Scalar WilkinsonShift(MatrixSubView<Mat> const & a)
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

            void ImplicitIter(MatrixSubView<Mat> & B,
                              MatrixSubView<Mat> * pU,
                              MatrixSubView<Mat> * pV)
            {
                int n = B.Rows();
                //Scalar mu = RayleighShift(B);
                Scalar mu = WilkinsonShift(B);
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
            //! @brief mUT 的乒乓队列
            PingPangView<Mat> mUTPingPang;
            //! @brief mSigma 的乒乓队列
            PingPangView<Mat> mSigmaPingPang;
            //! @brief mV 的乒乓队列
            PingPangView<Mat> mVPingPang;
    };

}

#endif


