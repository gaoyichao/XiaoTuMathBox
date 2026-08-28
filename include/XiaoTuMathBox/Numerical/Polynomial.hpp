#ifndef XTMB_NUMERICAL_POLYNOMIAL
#define XTMB_NUMERICAL_POLYNOMIAL

#include <vector>
#include <functional>
#include <XiaoTuMathBox/Common/Common.hpp>


namespace xiaotu {


    template <typename _Scalar>
    struct Traits<Polynomial<_Scalar>> {
        typedef _Scalar Scalar;
        static Scalar Abs(Scalar v) { return std::abs(v); }
        constexpr static Scalar zero = 0;
        constexpr static Scalar small_value = SMALL_VALUE;
    };

    template <typename Scalar>
    class Polynomial : public PolynomialBase<Polynomial<Scalar>>
    {
        public:
            typedef PolynomialBase<Polynomial> Base;
            typedef Traits<Polynomial> Trait;

            using Base::Coeff;
            using Base::IsZero;
            using Base::Normalize;
            using Base::NewtonRaphson;

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
            Polynomial(std::vector<Scalar> const & coefs, bool descending = false, Scalar const & tol = SMALL_VALUE)
            {
                if (coefs.empty())
                    mCoeffs.push_back(0);

                if (!descending) {
                    mCoeffs.assign(coefs.begin(), coefs.end());
                } else {
                    mCoeffs.assign(coefs.rbegin(), coefs.rend());
                }

                Normalize();
            }

            /**
             * @brief 拷贝赋值
             */
            Polynomial & operator = (Polynomial const & p)
            {
                mCoeffs.assign(p.mCoeffs.begin(), p.mCoeffs.end());
                return *this;
            }

            /**
             * @brief 构造零元
             */
            static Polynomial Zero()
            {
                return Polynomial();
            }

            /**
             * @brief 构造单位元
             */
            static Polynomial One()
            {
                return Polynomial({1.0});
            }

        public:

            /**
             * @brief 重置内存
             * 
             * @param [in] degree 支持的最高次数
             */
            Polynomial & ReAlloc(size_t degree)
            {
                mCoeffs.resize(degree + 1);
                std::fill(mCoeffs.begin(), mCoeffs.end(), Trait::zero);
                return *this;
            }

            /**
             * @brief 删除最高次项
             */
            Polynomial & PopBack()
            {
                if (1 == mCoeffs.size())
                    mCoeffs[0] = Trait::zero;
                else
                    mCoeffs.pop_back();
                return *this;
            }

            /**
             * @brief 添加最高次项
             */
            Polynomial & PushBack(Scalar an)
            {
                mCoeffs.push_back(an);
                return *this;
            }

            /**
             * @brief 多项式的次数
             */
            int Degree() const { return mCoeffs.size() - 1; }

            /**
             * @brief 获取指定次数项的系数
             */
            Scalar const & At(size_t i) const { return mCoeffs[i]; }

            /**
             * @brief 获取指定次数项的系数
             */
            Scalar & At(size_t i) { return mCoeffs[i]; }


            friend std::ostream & operator << (std::ostream & s, Polynomial const & m)
            {
                if (0 == m.Degree()) {
                    s << m[0];
                    return s;
                }

                bool first = true;
                for (int i = m.Degree(); i >= 0; --i) {
                    auto const & ai = m[i];
                    if (Trait::Abs(ai) < Trait::small_value)
                        continue;

                    if (first) {
                        if (ai < 0)
                            s << "-";
                    } else {
                        s << (ai < 0 ? " - " : " + ");
                    }

                    auto abs_ai = Trait::Abs(ai);
                    if (0 == i) {
                        s << abs_ai;
                    } else if (1 == i) {
                        if (Trait::Abs(abs_ai - 1)  < Trait::small_value)
                            s << abs_ai;
                        s << "x";
                    } else {
                        if (Trait::Abs(abs_ai - 1)  < Trait::small_value)
                            s << abs_ai;
                        s << "x^" << i;
                    }
                    first = false;
                }
                
                if (first)
                    s << "0";
                return s;
            }

            /**
             * @brief 一元二次多项式方程的根 \(ax^2 + bx + c = 0\)
             *
             * @param [out] x0 复数形式的根
             * @param [out] x1 复数形式的根
             * @return 是否为两个实根
             */
            bool QuadraticRoot(std::complex<Scalar> & x0, std::complex<Scalar> & x1)
            {
                return xiaotu::QuadraticRoot<Scalar>(At(2), At(1), At(0), x0, x1);
            }

        private:
            //! @brief 升序排列的多项式系数
            std::vector<Scalar> mCoeffs;
    };


    template <typename _Scalar>
    struct Traits<Polynomial<std::complex<_Scalar>>> {
        typedef std::complex<_Scalar> Complex;
        typedef Complex Scalar;

        constexpr static Scalar zero = Complex(0, 0);

        static _Scalar Abs(Scalar v) { return std::abs(v); }
        constexpr static _Scalar small_value = SMALL_VALUE;
    };


    template <typename Scalar>
    class Polynomial<std::complex<Scalar>> : public PolynomialBase<Polynomial<std::complex<Scalar>>>
    {
        public:
            typedef std::complex<Scalar> Complex;
            typedef PolynomialBase<Polynomial> Base;
            typedef Traits<Polynomial> Trait;

            using Base::Coeff;
            using Base::IsZero;
            using Base::Normalize;
            using Base::NewtonRaphson;

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
            Polynomial(std::vector<Complex> const & coefs, bool descending = false, Scalar const & tol = SMALL_VALUE)
            {
                if (coefs.empty())
                    mCoeffs.push_back(Trait::zero);

                if (!descending) {
                    mCoeffs.assign(coefs.begin(), coefs.end());
                } else {
                    mCoeffs.assign(coefs.rbegin(), coefs.rend());
                }

                Normalize();
            }

            /**
             * @brief 拷贝赋值
             */
            Polynomial & operator = (Polynomial const & p)
            {
                mCoeffs.assign(p.mCoeffs.begin(), p.mCoeffs.end());
                return *this;
            }

            /**
             * @brief 构造零元
             */
            static Polynomial Zero()
            {
                return Polynomial();
            }

            /**
             * @brief 构造单位元
             */
            static Polynomial One()
            {
                return Polynomial({Complex(1, 0)});
            }

        public:

            /**
             * @brief 重置内存
             * 
             * @param [in] degree 支持的最高次数
             */
            Polynomial & ReAlloc(size_t degree)
            {
                mCoeffs.resize(degree + 1);
                std::fill(mCoeffs.begin(), mCoeffs.end(), Trait::zero);
                return *this;
            }

            /**
             * @brief 删除最高次项
             */
            Polynomial & PopBack()
            {
                if (1 == mCoeffs.size())
                    mCoeffs[0] = Trait::zero;
                else
                    mCoeffs.pop_back();
                return *this;
            }

            /**
             * @brief 添加最高次项
             */
            Polynomial & PushBack(Scalar an)
            {
                mCoeffs.push_back(an);
                return *this;
            }


            /**
             * @brief 多项式的次数
             */
            int Degree() const { return mCoeffs.size() - 1; }

            /**
             * @brief 获取指定次数项的系数
             */
            Complex const & At(size_t i) const { return mCoeffs[i]; }

            /**
             * @brief 获取指定次数项的系数
             */
            Complex & At(size_t i) { return mCoeffs[i]; }


            friend std::ostream & operator << (std::ostream & s, Polynomial const & m)
            {
                if (0 == m.Degree()) {
                    s << m[0];
                    return s;
                }

                bool first = true;
                for (int i = m.Degree(); i >= 0; --i) {
                    auto const & ai = m[i];
                    if (Trait::Abs(ai) < Trait::small_value)
                        continue;

                    if (!first)
                        s << " + ";

                    if (0 == i) {
                        s << ai;
                    } else if (1 == i) {
                        s << ai << "x";
                    } else {
                        s << ai << "x^" << i;
                    }
                    first = false;
                }
                
                if (first)
                    s << Trait::zero;
                return s;
            }

            /**
             * @brief 一元二次多项式方程的根 \(ax^2 + bx + c = 0\)
             *
             * @param [out] x0 复数形式的根
             * @param [out] x1 复数形式的根
             * @return 是否为两个实根
             */
            bool QuadraticRoot(Complex & x0, Complex & x1)
            {
                Complex v = std::sqrt(At(1) * At(1) - 4.0 * At(2) * At(0));
                Complex tmp0 = At(1) + v;
                Complex tmp1 = At(1) - v;
                
                if (Trait::Abs(tmp0) < Trait::Abs(tmp1)) {
                    x0 = -tmp1 / At(0) * 0.5;
                    x1 = -2.0 * At(0) / tmp1;
                } else {
                    x0 = -2.0 * At(0) / tmp0;
                    x1 = -tmp0 / At(0) * 0.5;
                }

                return std::abs(x0.imag()) < Trait::small_value &&
                       std::abs(x1.imag()) < Trait::small_value;
            }

            /**
             * @brief 计算多项式方程的所有根
             * 
             * @param [out] roots 所有根
             */
            void AllRoots(std::vector<Complex> & roots)
            {
                using namespace std::placeholders;
                Complex x0 = -1;
                Complex x1 = 0;
                Complex x2 = 1;
                Complex b0;

                Polynomial P, Q;
                Polynomial *p_ptr = &P;
                Polynomial *q_ptr = &Q;
                *p_ptr = *this;

                while (p_ptr->Degree() > 1) {
                
                    Complex r = MullerRoot<Scalar>(
                        std::bind(&Polynomial::Evaluate, p_ptr, _1),
                        x0, x1, x2);

                    r = NewtonRaphson(r);

                    roots.push_back(r);

                    p_ptr->SyntheticDivide(r, *q_ptr, b0);
                    std::swap(p_ptr, q_ptr);
                }

                if (1 == p_ptr->Degree()) {
                    Complex r = -p_ptr->At(0) / p_ptr->At(1);
                    r = NewtonRaphson(r);
                    roots.push_back(r);
                }
            }


        private:
            //! @brief 升序排列的多项式系数
            std::vector<Complex> mCoeffs;
    };


}




#endif
