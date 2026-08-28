#ifndef XTMB_NUMERICAL_POLYNOMIAL_BASE_H
#define XTMB_NUMERICAL_POLYNOMIAL_BASE_H

#include <vector>
#include <functional>
#include <XiaoTuMathBox/Common/Common.hpp>

namespace xiaotu {

    template <typename Derived>
    class PolynomialBase
    {
        public:
            typedef Traits<Derived> Trait;
            typedef typename Traits<Derived>::Scalar Scalar;

        public:
            /**
             * @brief 获取某次项系数
             */
            inline Scalar const & operator[](size_t i) const { return At(i); }

            /**
             * @brief 获取某次项系数
             */
            inline Scalar & operator[](size_t i) { return At(i); }

            /**
             * @brief 获取指定次数项的系数，超过最高次时返回 0
             */
            Scalar Coeff(size_t i) const
            {
                if (i <= Degree())
                    return At(i);

                return Trait::zero;
            }

            /**
             * @brief 去除最高阶零系数
             */
            Derived & Normalize()
            {
                while (Degree() > 0 && std::abs(At(Degree())) < Trait::small_value) {
                    PopBack();
                }
                return derived();
            }

            /**
             * @brief 判定是否为 0 元
             */
            bool IsZero() const
            {
                return (0 == Degree() && Trait::Abs(At(0)) < Trait::small_value);
            }

            /**
             * @brief 计算多项式的值， 嵌套乘法
             */
            Scalar Evaluate(Scalar const & x) const
            {
                Scalar re = 0;
                for (int i = Degree(); i >= 0; --i) {
                    re = re * x + At(i);
                }
                return re;
            }

            /**
             * @brief 计算多项式的值
             */
            Scalar operator()(Scalar const & x) const { return Evaluate(x); }

            /**
             * @brief 通过嵌套乘法计算 P(x_0), P'(x_0)
             * 
             * @param [in] x0 参考点
             * @param [out] Px0 多项式值 \(P(x_0)\)
             * @param [out] DPx0 一阶导数值 \(P'(x_0)\)
             */
            void Horner(Scalar const & x0, Scalar & Px0, Scalar & DPx0)
            {
                size_t n = Degree();
                Scalar bk = At(n);
                DPx0 = At(n);

                for (size_t k = (n-1); k >= 1; --k) {
                    // b_k =   a_k  + b_{k+1} x_0
                    bk = At(k) + bk * x0;
                    // Q(x) = b_nx^{n-1} + b_{n-1}x^{n-2} + \cdots + b_2 x + b_1
                    // Q(x) = \left(\cdots\left(b_n x + b_{n-1}\right)x + \cdots + b_2\right)x + b_1
                    DPx0 = DPx0 * x0 + bk;
                }
                Px0 = At(0) + bk * x0;
            }

            /**
             * @brief 综合除法 P(x) = (x - x0) Q(x) + b0
             * 
             * @param [in] x0 参考点
             * @param [out] Q 商式 Q(x)
             * @param [out] b0 余数
             */
            void SyntheticDivide(Scalar const & x0, Derived & Q, Scalar & b0) const
            {
                size_t n = this->Degree();
                
                if (n == 0) {
                    Q = Derived::Zero();
                    b0 = At(0);
                    return;
                }

                Q.ReAlloc(n - 1);

                Scalar bk = At(n);
                Q[n - 1] = bk;

                for (size_t k = (n-1); k >= 1; --k) {
                    bk = At(k) + bk * x0;
                    Q[k - 1] = bk;
                }

                b0 = At(0) + bk * x0;
            }

            /**
             * @brief 多项式带余除法 A(x) = B(x) Q(x) + R(x)
             * 
             * @param [in] divisor 除式 B(x)
             * @param [out] quotient 商式 Q(x)
             * @param [out] remainder 余式 R(x)
             * @return 余式的次数
             */
            size_t Divide(Derived const & divisor, Derived & quotient, Derived & remainder) const
            {
                assert(!divisor.IsZero());
                remainder = derived();

                // 如果被除式次数小于除式次数，商为0，余数就是被除式本身
                if (this->Degree() < divisor.Degree()) {
                    quotient = Derived::Zero();
                    return remainder.Degree();
                }

                // 商的最高可能次数 = 被除式次数 - 除式次数
                size_t q_deg = derived().Degree() - divisor.Degree();
                quotient.ReAlloc(q_deg);

                size_t b_deg = divisor.Degree();
                Scalar b_lead = divisor[b_deg];

                while (remainder.Degree() >= b_deg) {
                    size_t r_deg = remainder.Degree();
                    Scalar r_lead = remainder[r_deg];

                    if (Trait::Abs(r_lead) < Trait::small_value) {
                        remainder.PopBack();
                        continue;
                    }

                    size_t pow_diff = r_deg - b_deg;
                    Scalar coeff_quotient = r_lead / b_lead;
                    quotient[pow_diff] = coeff_quotient;
                    for (size_t i = 0; i <= b_deg; ++i) {
                        remainder[i + pow_diff] -= coeff_quotient * divisor[i];
                    }

                    remainder.PopBack();
                }

                quotient.Normalize();
                remainder.Normalize();

                return remainder.Degree();
            }


            Scalar NewtonRaphson(Scalar x0, int max_iter = 100)
            {
                Scalar re = x0;
                for (int i = 0; i < max_iter; ++i) {
                    Scalar y, dydx;
                    this->Horner(re, y, dydx);

                    if (Trait::Abs(y) < Trait::small_value)
                        return re;
                    assert(Trait::zero != dydx);

                    re = x0 - y / dydx;
                    if (std::abs(re - x0) < Trait::small_value)
                        return re;
                    x0 = re;
                }

                return re;
            }

        public:
            /**
             * @brief 获取多项式次数
             */
            inline int Degree() const { return derived().Degree(); }

            /**
             * @brief 重置内存并填充 0
             */ 
            inline Derived & ReAlloc(size_t degree) { return derived().ReAlloc(degree); }

            /**
             * @brief 删除最高次项
             */
            inline Derived & PopBack() { return derived().PopBack(); }

            /**
             * @brief 添加最高次项
             */
            inline Derived & PushBack(Scalar an) { return derived().PushBack(an); }

            /**
             * @brief 获取指定次数项的系数, 直接访问系数内存, 派生类需要保证对象可用
             */
            inline Scalar const & At(size_t i) const { return derived().At(i); }

            /**
             * @brief 获取指定次数项的系数, 直接访问系数内存, 派生类需要保证对象可用
             */
            inline Scalar & At(size_t i) { return derived().At(i); }

        private:
            Derived & derived() { return *static_cast<Derived*>(this); }
            Derived const & derived() const { return *static_cast<const Derived*>(this); }
    };

}


#endif
