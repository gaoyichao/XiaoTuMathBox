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
             * @brief 拷贝赋值
             */
            Polynomial & operator = (Polynomial const & p)
            {
                mCoeffs.assign(p.mCoeffs.begin(), p.mCoeffs.end());
                return *this;
            }

            static Polynomial Zero()
            {
                return Polynomial();
            }

            static Polynomial One()
            {
                return Polynomial({1.0});
            }

        public:

            /**
             * @brief 去除最高阶零系数
             */
            Polynomial & Normalize(DataType const & tol = SMALL_VALUE)
            {
                while (mCoeffs.size() > 1 && std::abs(mCoeffs.back()) < tol) {
                    mCoeffs.pop_back();
                }
                return *this;
            }

            /**
             * @brief 多项式的次数
             */
            int Degree() const { return mCoeffs.size() - 1; }

            /**
             * @brief 判定是否为 0 元
             */
            bool IsZero(DataType const & tol = SMALL_VALUE) const
            {
                return (0 == Degree() && std::abs(mCoeffs[0]) < tol);
            }

            /**
             * @brief 计算多项式的值， 嵌套乘法
             */
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
             * @brief 综合除法 P(x) = (x - x0) Q(x) + b0
             * 
             * @param [in] x0 参考点
             * @param [out] Q 商式 Q(x)
             * @param [out] b0 余数
             */
            void SyntheticDivide(DataType const & x0, Polynomial & Q, DataType & b0) const
            {
                size_t n = this->Degree();
                
                if (n == 0) {
                    Q = Polynomial::Zero();
                    b0 = mCoeffs[0];
                    return;
                }

                auto & q_coeffs = Q.mCoeffs;
                q_coeffs.resize(n);
                std::fill(q_coeffs.begin(), q_coeffs.end(), 0);

                DataType bk = mCoeffs[n];
                q_coeffs[n - 1] = bk;

                for (size_t k = (n-1); k >= 1; --k) {
                    bk = mCoeffs[k] + bk * x0;
                    q_coeffs[k - 1] = bk;
                }

                b0 = mCoeffs[0] + bk * x0;
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

            /**
             * @brief 多项式带余除法 A(x) = B(x) Q(x) + R(x)
             * 
             * @param [in] divisor 除式 B(x)
             * @param [out] quotient 商式 Q(x)
             * @param [out] remainder 余式 R(x)
             * @return 余式的次数
             */
            size_t Divide(Polynomial const & divisor, Polynomial & quotient, Polynomial& remainder) const
            {
                assert(!divisor.IsZero());
                remainder = *this;

                // 如果被除式次数小于除式次数，商为0，余数就是被除式本身
                if (this->Degree() < divisor.Degree()) {
                    quotient = Polynomial::Zero();
                    return remainder.Degree();
                }

                auto & q_coeffs = quotient.mCoeffs;
                auto & r_coeffs = remainder.mCoeffs;
                auto const & b_coeffs = divisor.mCoeffs;

                // 商的最高可能次数 = 被除式次数 - 除式次数
                size_t q_deg = this->Degree() - divisor.Degree();
                q_coeffs.resize(q_deg + 1);
                std::fill(q_coeffs.begin(), q_coeffs.end(), 0);

                size_t b_deg = divisor.Degree();
                double b_lead = b_coeffs[b_deg];

                while ((r_coeffs.size() - 1) >= b_deg) {
                    size_t r_deg = r_coeffs.size() - 1;
                    double r_lead = r_coeffs[r_deg];

                    if (std::abs(r_lead) < SMALL_VALUE) {
                        r_coeffs.pop_back();
                        continue;
                    }

                    size_t pow_diff = r_deg - b_deg;
                    DataType coeff_quotient = r_lead / b_lead;
                    q_coeffs[pow_diff] = coeff_quotient;
                    for (size_t i = 0; i <= b_deg; ++i) {
                        r_coeffs[i + pow_diff] -= coeff_quotient * b_coeffs[i];
                    }

                    r_coeffs.pop_back();
                }

                if (r_coeffs.empty())
                    r_coeffs.push_back(0);

                quotient.Normalize();
                remainder.Normalize();

                return remainder.Degree();
            }

        public:

            DataType GetCoeff(size_t i) const
            {
                if (i < mCoeffs.size())
                    return mCoeffs[i];
                return 0;
            }

            /**
             * @brief 获取某次项系数
             */
            DataType const & operator[](size_t i) const
            {
                assert(i < mCoeffs.size());
                return mCoeffs[i];
            }

            /**
             * @brief 获取某次项系数
             */
            DataType & operator[](size_t i)
            {
                assert(i < mCoeffs.size());
                return mCoeffs[i];
            }

            /**
             * @brief 计算多项式的值
             */
            DataType operator()(DataType const & x) const { return Evaluate(x); }

            ////////////////////////////////////////////////////////
            //
            //  ==, !=
            //
            ////////////////////////////////////////////////////////

            bool operator == (Polynomial const & other) const
            {
                size_t deg = std::max(this->Degree(), other.Degree());
                
                for (size_t i = 0; i <= deg; ++i) {
                    auto a = this->GetCoeff(i);
                    auto b = other.GetCoeff(i);
                    if (std::abs(a - b) >= SMALL_VALUE)
                        return false;
                }

                return true;
            }

            bool operator != (Polynomial const & other) const
            {
                return !(*this == other);
            }

            friend bool operator == (DataType const & a, Polynomial const & b)
            {
                if (b.Degree() > 0)
                    return false;
                if (std::abs(a - b.GetCoeff(0)) >= SMALL_VALUE)
                    return false;
                return true;
            }

            friend bool operator != (DataType const & a, Polynomial const & b)
            {
                return !(a == b);
            }

            friend bool operator == (Polynomial const & a, DataType const & b)
            {
                return b == a;
            }

            friend bool operator != (Polynomial const & a, DataType const & b)
            {
                return !(b == a);
            }

            ////////////////////////////////////////////////////////
            //
            //  c = a + b
            //
            ////////////////////////////////////////////////////////

            friend Polynomial operator + (Polynomial const & a, Polynomial const & b)
            {
                Polynomial re;
                size_t degree = std::max(a.Degree(), b.Degree());

                re.mCoeffs.resize(degree + 1);
                for (size_t i = 0; i <= degree; ++i) {
                    re.mCoeffs[i] = a.GetCoeff(i) + b.GetCoeff(i);
                }
                
                re.Normalize();
                return re;
            }

            friend Polynomial operator + (DataType const & a, Polynomial const & b)
            {
                Polynomial re = b;
                re.mCoeffs[0] += a;
                return re;
            }

            friend Polynomial operator + (Polynomial const & a, DataType const & b)
            {
                return b + a;
            }

            ////////////////////////////////////////////////////////
            //
            //  c = a * b
            //
            ////////////////////////////////////////////////////////

            friend Polynomial operator * (Polynomial const & a, Polynomial const & b)
            {
                if (a.IsZero() || b.IsZero())
                    return Polynomial::Zero();

                size_t degree = a.Degree() + b.Degree();
                Polynomial re;

                re.mCoeffs.resize(degree + 1, 0);
                for (size_t i = 0; i <= a.Degree(); ++i) {
                    for (size_t j = 0; j <= b.Degree(); ++j) {
                        re.mCoeffs[i + j] += a.mCoeffs[i] * b.mCoeffs[j];
                    }
                }

                re.Normalize();
                return re;
            }

            friend Polynomial operator * (DataType const & a, Polynomial const & b)
            {
                if (std::abs(a) < SMALL_VALUE || b.IsZero())
                    return Polynomial::Zero();

                Polynomial re = b;
                for (size_t i = 0; i < re.mCoeffs.size(); ++i)
                    re.mCoeffs[i] *= a;

                re.Normalize();
                return re;
            }

            friend Polynomial operator * (Polynomial const & a, DataType const & b)
            {
                return b * a;
            }

            ////////////////////////////////////////////////////////
            //
            //  c = a - b
            //
            ////////////////////////////////////////////////////////

            friend Polynomial operator - (Polynomial const & a, Polynomial const & b)
            {
                Polynomial re;
                size_t degree = std::max(a.Degree(), b.Degree());

                re.mCoeffs.resize(degree + 1);
                for (size_t i = 0; i <= degree; ++i) {
                    re.mCoeffs[i] = a.GetCoeff(i) - b.GetCoeff(i);
                }
                
                re.Normalize();
                return re;
            }

            friend Polynomial operator - (DataType const & a, Polynomial const & b)
            {
                Polynomial re = -1 * b;
                re.mCoeffs[0] += a;
                return re;
            }

            friend Polynomial operator - (Polynomial const & a, DataType const & b)
            {
                Polynomial re = b;
                re.mCoeffs[0] += a;
                return re;
            }



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
