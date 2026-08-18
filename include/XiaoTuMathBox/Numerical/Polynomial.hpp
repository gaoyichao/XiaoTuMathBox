#ifndef XTMB_NUMERICAL_POLYNOMIAL
#define XTMB_NUMERICAL_POLYNOMIAL

#include <vector>
#include <functional>
#include <XiaoTuMathBox/Common/Common.hpp>


namespace xiaotu {


    template <typename DataType>
    class Polynomial {

        public:

            Polynomial() : mCoeffs(1, 0) {}
            Polynomial(Polynomial const & p) : mCoeffs(p.mCoeffs) {}

            /**
             * @brief 构造函数
             * 
             * @param [in] coefs 多项式系数, 一个 n 次多项式需要有 n+1 个系数
             * @param [in] descending coefs 是否降序排列，默认否
             * @param [in] tol 去除高阶零系数的阈值
             */
            Polynomial(std::vector<DataType> const & coefs, bool descending = false, DataType const & tol = SMALL_VALUE)
            {
                if (coefs.empty())
                    mCoeffs.push_back(0);

                if (!descending) {
                    mCoeffs.assign(coefs.begin(), coefs.end());
                } else {
                    mCoeffs.assign(coefs.rbegin(), coefs.rend());
                }

                Normalize(tol);
            }

            /**
             * @brief 去除高阶零系数
             */
            Polynomial & Normalize(DataType const & tol = SMALL_VALUE)
            {
                while (mCoeffs.size() > 1 && std::abs(mCoeffs.back()) < tol) {
                    mCoeffs.pop_back();
                }
                return *this;
            }

        public:

            //! @brief 多项式的次数
            int Degree() const { return mCoeffs.size() - 1; }

            //! @brief 获取某次项系数
            DataType const & operator[](size_t n) const
            {
                assert(n < mCoeffs.size());
                return mCoeffs[n];
            }

            //! @brief 计算多项式的值， 嵌套乘法
            DataType Evaluate(DataType const & x) const
            {
                DataType re = 0;
                for (int i = Degree(); i >= 0; --i) {
                    re = re * x + mCoeffs[i];
                }
                return re;
            }

            /**
             * @brief 通过嵌套乘法计算 P(x_0), P'(x_0)
             * 
             * @param [in] x0 参考点
             * @param [out] Px0 多项式值 \(P(x_0)\)
             * @param [out] DPx0 一阶导数值 \(P'(x_0)\)
             */
            void Horner(DataType const & x0, DataType & Px0, DataType & DPx0)
            {
                size_t n = Degree();
                DataType bk = mCoeffs[n];
                DPx0 = mCoeffs[n];

                for (size_t k = (n-1); k >= 1; --k) {
                    // b_k =   a_k  + b_{k+1} x_0
                    bk = mCoeffs[k] + bk * x0;
                    // Q(x) = b_nx^{n-1} + b_{n-1}x^{n-2} + \cdots + b_2 x + b_1
                    // Q(x) = \left(\cdots\left(b_n x + b_{n-1}\right)x + \cdots + b_2\right)x + b_1
                    DPx0 = DPx0 * x0 + bk;
                }
                Px0 = mCoeffs[0] + bk * x0;
            }

            /**
             * @brief 一元二次多项式方程的根 \(ax^2 + bx + c = 0\)
             *
             * @param [out] x0 复数形式的根
             * @param [out] x1 复数形式的根
             * @return 是否为两个实根
             */
            bool QuadraticRoot(std::complex<DataType> & x0, std::complex<DataType> & x1)
            {
                return xiaotu::QuadraticRoot<DataType>(mCoeffs[2], mCoeffs[1], mCoeffs[0], x0, x1);
            }


        public:

            DataType operator()(DataType const & x) const { return Evaluate(x); }

            friend std::ostream & operator << (std::ostream & s, Polynomial const & m)
            {
                if (0 == m.Degree()) {
                    s << m.mCoeffs[0];
                    return s;
                }

                bool first = true;
                for (int i = m.Degree(); i >= 0; --i) {
                    auto const & ai = m.mCoeffs[i];
                    if (std::abs(ai) < SMALL_VALUE)
                        continue;

                    if (first) {
                        if (ai < 0)
                            s << "-";
                    } else {
                        s << (ai < 0 ? " - " : " + ");
                    }

                    auto abs_ai = std::abs(ai);
                    if (0 == i) {
                        s << abs_ai;
                    } else if (1 == i) {
                        if (std::abs(abs_ai - 1) > SMALL_VALUE)
                            s << abs_ai;
                        s << "x";
                    } else {
                        if (std::abs(abs_ai - 1) > SMALL_VALUE)
                            s << abs_ai;
                        s << "x^" << i;
                    }
                    first = false;
                }
                
                if (first)
                    s << "0";
                return s;
            }

        private:
            //! @brief 升序排列的多项式系数
            std::vector<DataType> mCoeffs;
    };


}




#endif
