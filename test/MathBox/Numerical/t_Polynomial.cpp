#include <iostream>

#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Numerical/Numerical.hpp>

#include <cmath>
#include <gtest/gtest.h>



TEST(Polynomial, Creative)
{
    {
        xiaotu::Polynomial<double> P;
        EXPECT_EQ(P.Degree(), 0);
        XTLog(std::cout) << P << std::endl;
    }

    {
        xiaotu::Polynomial<double> P({1.0, 2.0, 3.0});
        EXPECT_EQ(P.Degree(), 2);
        EXPECT_DOUBLE_EQ(P(1.0), 1.0 + 2.0 + 3.0);
        EXPECT_DOUBLE_EQ(P(2.0), 1.0 + 2.0 * 2.0 + 3.0 * 2.0 * 2.0);
        XTLog(std::cout) << P << std::endl;
    }

    {
        xiaotu::Polynomial<double> P({1.0, 2.0, 3.0}, true);
        EXPECT_EQ(P.Degree(), 2);
        EXPECT_DOUBLE_EQ(P(1.0), 1.0 + 2.0 + 3.0);
        EXPECT_DOUBLE_EQ(P(2.0), 1.0 * 2.0 * 2.0 + 2.0 * 2.0 + 3.0 );
        XTLog(std::cout) << P << std::endl;
    }

    {
        xiaotu::Polynomial<double> P({0.0, 1.0, 2.0, 3.0});
        EXPECT_EQ(P.Degree(), 3);
        EXPECT_DOUBLE_EQ(P(1.0), 1.0 + 2.0 + 3.0);
        EXPECT_DOUBLE_EQ(P(2.0), 1.0 * 2.0 + 2.0 * 2.0 * 2.0 + 3.0 * 2.0 * 2.0 * 2.0);
        XTLog(std::cout) << P << std::endl;
    }

    {
        xiaotu::Polynomial<double> P({0.0, 1.0, 2.0, 3.0}, true);
        EXPECT_EQ(P.Degree(), 2);
        EXPECT_DOUBLE_EQ(P(1.0), 1.0 + 2.0 + 3.0);
        EXPECT_DOUBLE_EQ(P(2.0), 1.0 * 2.0 * 2.0 + 2.0 * 2.0 + 3.0);
        XTLog(std::cout) << P << std::endl;
    }

    {
        auto zero = xiaotu::Polynomial<double>::Zero();
        XTLog(std::cout) << "零元: " << zero << std::endl;

        auto one = xiaotu::Polynomial<double>::One();
        XTLog(std::cout) << "单位元: " << one << std::endl;
    }
}


TEST(Polynomial, Copy)
{
    xiaotu::Polynomial<double> P1({0.0, 1.0, 2.0, 3.0});
    xiaotu::Polynomial<double> P2({0.0, 1.0, 2.0});
    
    XTLog(std::cout) << "P1: " << P1 << std::endl;
    XTLog(std::cout) << "P2: " << P2 << std::endl;
    
    P1 = P2;
    
    XTLog(std::cout) << "P1: " << P1 << std::endl;
    XTLog(std::cout) << "P2: " << P2 << std::endl;
    
    P2[0] = 3.0;
    
    XTLog(std::cout) << "P1: " << P1 << std::endl;
    XTLog(std::cout) << "P2: " << P2 << std::endl;
}


TEST(Polynomial, Horner)
{
    xiaotu::Polynomial<double> P({0.0, 1.0, 2.0, 3.0});
    EXPECT_EQ(P.Degree(), 3);
    EXPECT_DOUBLE_EQ(P(1.0), 1.0 + 2.0 + 3.0);
    EXPECT_DOUBLE_EQ(P(2.0), 1.0 * 2.0 + 2.0 * 2.0 * 2.0 + 3.0 * 2.0 * 2.0 * 2.0);
    XTLog(std::cout) << "P(x)  = " << P << std::endl;

    xiaotu::Polynomial<double> DP({1.0, 2.0 * 2.0, 3.0 * 3.0});
    EXPECT_EQ(DP.Degree(), 2);
    XTLog(std::cout) << "P'(x) = " << DP << std::endl;

    {
        double p, dp;
        P.Horner(1.0, p, dp);
        EXPECT_DOUBLE_EQ(P(1.0), p);
        EXPECT_DOUBLE_EQ(DP(1.0), dp);
        XTLog(std::cout) << "P(1.0) = " << p << ", P'(1.0) = " << dp << std::endl;
    }

    {
        double p, dp;
        P.Horner(2.0, p, dp);
        EXPECT_DOUBLE_EQ(P(2.0), p);
        EXPECT_DOUBLE_EQ(DP(2.0), dp);
        XTLog(std::cout) << "P(2.0) = " << p << ", P'(2.0) = " << dp << std::endl;
    }
}


TEST(Polynomial, QuadraticRoot)
{
    std::complex<double> x0, x1;
    {
        xiaotu::Polynomial<double> p({1.0, -10000.0001, 1}, true);
        bool re = p.QuadraticRoot(x0, x1);
        EXPECT_TRUE(re);
        EXPECT_TRUE(std::abs(x0 - 10000.0) < SMALL_VALUE || std::abs(x0 - 0.0001) < SMALL_VALUE);
        EXPECT_TRUE(std::abs(x1 - 10000.0) < SMALL_VALUE || std::abs(x1 - 0.0001) < SMALL_VALUE);
        EXPECT_TRUE(x0 != x1);
    }

    {
        xiaotu::Polynomial<double> p({1.0, 0.0, 1}, true);
        bool re = p.QuadraticRoot(x0, x1);
        EXPECT_FALSE(re);
        EXPECT_TRUE(x0.real() == 0);
        EXPECT_TRUE(x1.real() == 0);

        EXPECT_TRUE(std::abs(x0.imag() - 1.0) < SMALL_VALUE || std::abs(x0.imag() + 1.0) < SMALL_VALUE);
        EXPECT_TRUE(std::abs(x1.imag() - 1.0) < SMALL_VALUE || std::abs(x1.imag() + 1.0) < SMALL_VALUE);
    }

    {
        xiaotu::Polynomial<double> p({2.0, 2.0, 1}, true);
        bool re = p.QuadraticRoot(x0, x1);
        EXPECT_FALSE(re);

        XTLog(std::cout) << x0 << x1 << std::endl;
    }
}


TEST(Polynomial, Plus)
{
    {
        xiaotu::Polynomial<double> P0({1.0, 2.0, 3.0});
        xiaotu::Polynomial<double> P1({1.0, 2.0, 3.0});
        auto P2 = P0 + P1;
        XTLog(std::cout) << P2 << std::endl;

        P2 = P2 + 5;
        XTLog(std::cout) << P2 << std::endl;

        P2 = -5 + P2;
        XTLog(std::cout) << P2 << std::endl;

        P2 = P2 - P1;
        XTLog(std::cout) << P2 << std::endl;
    }

    {
        xiaotu::Polynomial<double> P0({1.0, 2.0, 3.0});
        auto P2 = 2 - P0;
        XTLog(std::cout) << P2 << std::endl;
    }
}


TEST(Polynomial, Mult)
{
    {
        auto P0 = xiaotu::Polynomial<double>::Zero();
        auto P1 = xiaotu::Polynomial<double>::Zero();
        auto P2 = P0 * P1;
        XTLog(std::cout) << P2 << std::endl;
    }

    {
        xiaotu::Polynomial<double> P0({1.0, 2.0, 3.0});
        xiaotu::Polynomial<double> P1({1.0, 2.0, 3.0});
        auto P2 = P0 * P1;
        XTLog(std::cout) << P2 << std::endl;
    }

    {
        xiaotu::Polynomial<double> P0({1.0, 2.0, 3.0});
        auto P2 = 2 * P0;
        XTLog(std::cout) << P2 << std::endl;

        P2 = P2 * 0.5;
        XTLog(std::cout) << P2 << std::endl;
    }

}


TEST(Polynomial, Divide)
{
    {
        xiaotu::Polynomial<double> a({1.0, 2.0, 3.0});
        xiaotu::Polynomial<double> b({1.0, 2.0, 3.0});

        auto q = xiaotu::Polynomial<double>::Zero();
        auto r = xiaotu::Polynomial<double>::Zero();
        auto r_deg = a.Divide(b, q, r);
        
        XTLog(std::cout) << "r_deg:" << r_deg << std::endl;
        XTLog(std::cout) << "q:" << q << std::endl;
        XTLog(std::cout) << "r:" << r << std::endl;
    }

    XTLog(std::cout) << "-----------------------" << std::endl;

    {
        xiaotu::Polynomial<double> a({1.0, 2.0, 3.0});
        xiaotu::Polynomial<double> b({1.0, 2.0});

        auto q = xiaotu::Polynomial<double>::Zero();
        auto r = xiaotu::Polynomial<double>::Zero();
        auto r_deg = a.Divide(b, q, r);
        
        XTLog(std::cout) << "r_deg:" << r_deg << std::endl;
        XTLog(std::cout) << "q:" << q << std::endl;
        XTLog(std::cout) << "r:" << r << std::endl;

        auto aa = q * b + r;
        XTLog(std::cout) << "a:" << a << std::endl;
        XTLog(std::cout) << "aa:" << aa << std::endl;
    }

    XTLog(std::cout) << "-----------------------" << std::endl;

    {
        xiaotu::Polynomial<double> a({1.0, 2.0, 3.0});

        auto q = xiaotu::Polynomial<double>::Zero();
        double b0;
        a.SyntheticDivide(2.0, q, b0);
        
        XTLog(std::cout) << "a:" << a << std::endl;
        XTLog(std::cout) << "q:" << q << std::endl;
        XTLog(std::cout) << "r:" << b0 << std::endl;

        auto aa = q * xiaotu::Polynomial<double>({-2.0, 1}) + b0;
        XTLog(std::cout) << "aa:" << aa << std::endl;
    }

}


TEST(Polynomial, Equal)
{
    {
        xiaotu::Polynomial<double> a({1.0, 2.0, 3.0});
        xiaotu::Polynomial<double> b({1.0, 2.0, 3.0});

        EXPECT_TRUE(a == b);
        EXPECT_EQ(a, b);
    }

    {
        xiaotu::Polynomial<double> a({1.0, 2.0, 3.0});
        xiaotu::Polynomial<double> b({4.0, 2.0, 3.0});

        EXPECT_TRUE(a != b);
        EXPECT_NE(a, b);
    }

    {
        auto zero = xiaotu::Polynomial<double>::Zero();
        EXPECT_TRUE(0 == zero);
        EXPECT_EQ(0, zero);
        EXPECT_TRUE(zero == 0);
        EXPECT_EQ(zero, 0);
    }

    {
        auto one = xiaotu::Polynomial<double>::One();
        EXPECT_TRUE(1 == one);
        EXPECT_EQ(1, one);
        EXPECT_TRUE(one == 1);
        EXPECT_EQ(one, 1);
    }

    {
        xiaotu::Polynomial<double> a({1.0, 2.0, 3.0});
        EXPECT_FALSE(1 == a);
        EXPECT_NE(1, a);
        EXPECT_FALSE(a == 1);
        EXPECT_NE(a, 1);
    }
}


/////////////////////////////////////////////////////////////////

TEST(ComplexPolynomial, Creative)
{
    {
        xiaotu::Polynomial<std::complex<double>> P;
        EXPECT_EQ(P.Degree(), 0);
        XTLog(std::cout) << P << std::endl;
    }


    {
        xiaotu::Polynomial<std::complex<double>> P({1.0, 2.0, 3.0});
        EXPECT_EQ(P.Degree(), 2);
        EXPECT_DOUBLE_EQ(P(1.0).real(), 1.0 + 2.0 + 3.0);
        EXPECT_DOUBLE_EQ(P(2.0).real(), 1.0 + 2.0 * 2.0 + 3.0 * 2.0 * 2.0);
        XTLog(std::cout) << P << std::endl;
    }

    {
        auto zero = xiaotu::Polynomial<std::complex<double>>::Zero();
        XTLog(std::cout) << "零元: " << zero << std::endl;

        auto one = xiaotu::Polynomial<std::complex<double>>::One();
        XTLog(std::cout) << "单位元: " << one << std::endl;
    }
}

TEST(ComplexPolynomial, Copy)
{
    xiaotu::Polynomial<std::complex<double>> P1({0.0, 1.0, 2.0, 3.0});
    xiaotu::Polynomial<std::complex<double>> P2({0.0, 1.0, 2.0});
    
    XTLog(std::cout) << "P1: " << P1 << std::endl;
    XTLog(std::cout) << "P2: " << P2 << std::endl;
    
    P1 = P2;
    
    XTLog(std::cout) << "P1: " << P1 << std::endl;
    XTLog(std::cout) << "P2: " << P2 << std::endl;
    
    P2[0] = std::complex<double>(3.0, 2.9);
    
    XTLog(std::cout) << "P1: " << P1 << std::endl;
    XTLog(std::cout) << "P2: " << P2 << std::endl;
}

TEST(ComplexPolynomial, Horner)
{
    using Complex = std::complex<double>;
    xiaotu::Polynomial<Complex> P({0.0, 1.0, 2.0, 3.0});
    EXPECT_EQ(P.Degree(), 3);
    EXPECT_DOUBLE_EQ(P(1.0).real(), 1.0 + 2.0 + 3.0);
    EXPECT_DOUBLE_EQ(P(2.0).real(), 1.0 * 2.0 + 2.0 * 2.0 * 2.0 + 3.0 * 2.0 * 2.0 * 2.0);
    XTLog(std::cout) << "P(x)  = " << P << std::endl;

    xiaotu::Polynomial<Complex> DP({1.0, 2.0 * 2.0, 3.0 * 3.0});
    EXPECT_EQ(DP.Degree(), 2);
    XTLog(std::cout) << "P'(x) = " << DP << std::endl;

    {
        Complex p, dp;
        Complex x(1.0, 2.0);
        P.Horner(x, p, dp);
        EXPECT_EQ(P(x), p);
        EXPECT_EQ(DP(x), dp);
        XTLog(std::cout) << "P(" << x << ") = " << p << ", P' = " << dp << std::endl;
    }

    {
        Complex p, dp;
        Complex x(3.0, 2.0);

        P.Horner(x, p, dp);
        EXPECT_EQ(P(x), p);
        EXPECT_EQ(DP(x), dp);
        XTLog(std::cout) << "P(" << x << ") = " << p << ", P' = " << dp << std::endl;
    }

}

TEST(ComplexPolynomial, QuadraticRoot)
{
    using Complex = std::complex<double>;
    Complex x0, x1;
    {
        xiaotu::Polynomial<Complex> p({1.0, -10000.0001, 1}, true);
        bool re = p.QuadraticRoot(x0, x1);

        XTLog(std::cout) << "x0:" << x0 << std::endl;
        XTLog(std::cout) << "x1:" << x1 << std::endl;
        EXPECT_TRUE(re);
        EXPECT_TRUE(std::abs(x0 - Complex(10000.0, 0)) < 1e-6 || std::abs(x0 - 0.0001) < 1e-6);
        EXPECT_TRUE(std::abs(x1 - 10000.0) < 1e-6 || std::abs(x1 - 0.0001) < 1e-6);
        EXPECT_TRUE(x0 != x1);
    }

    {
        xiaotu::Polynomial<Complex> p({1.0, 0.0, 1.0}, true);
        bool re = p.QuadraticRoot(x0, x1);

        XTLog(std::cout) << "x0:" << x0 << std::endl;
        XTLog(std::cout) << "x1:" << x1 << std::endl;
        EXPECT_FALSE(re);
    }
}

TEST(ComplexPolynomial, Plus)
{
    using Complex = std::complex<double>;
    {
        xiaotu::Polynomial<Complex> P0({1.0, 2.0, 3.0});
        xiaotu::Polynomial<Complex> P1({1.0, 2.0, 3.0});
        auto P2 = P0 + P1;
        XTLog(std::cout) << P2 << std::endl;

        P2 = P2 + 5;
        XTLog(std::cout) << P2 << std::endl;

        P2 = -5 + P2;
        XTLog(std::cout) << P2 << std::endl;

        P2 = P2 - P1;
        XTLog(std::cout) << P2 << std::endl;
    }

    {
        xiaotu::Polynomial<Complex> P0({1.0, 2.0, 3.0});
        auto P2 = 2 - P0;
        XTLog(std::cout) << P2 << std::endl;
    }
}

TEST(ComplexPolynomial, Mult)
{
    using Complex = std::complex<double>;
    {
        auto P0 = xiaotu::Polynomial<Complex>::Zero();
        auto P1 = xiaotu::Polynomial<Complex>::Zero();
        auto P2 = P0 * P1;
        XTLog(std::cout) << P2 << std::endl;
    }

    {
        xiaotu::Polynomial<Complex> P0({1.0, 2.0, 3.0});
        xiaotu::Polynomial<Complex> P1({1.0, 2.0, 3.0});
        auto P2 = P0 * P1;
        XTLog(std::cout) << P2 << std::endl;
    }

    {
        xiaotu::Polynomial<Complex> P0({1.0, 2.0, 3.0});
        auto P2 = 2 * P0;
        XTLog(std::cout) << P2 << std::endl;

        P2 = P2 * 0.5;
        XTLog(std::cout) << P2 << std::endl;
    }

}


TEST(ComplexPolynomial, Divide)
{
    using Complex = std::complex<double>;
    {
        xiaotu::Polynomial<Complex> a({1.0, 2.0, 3.0});
        xiaotu::Polynomial<Complex> b({1.0, 2.0, 3.0});

        auto q = xiaotu::Polynomial<Complex>::Zero();
        auto r = xiaotu::Polynomial<Complex>::Zero();
        auto r_deg = a.Divide(b, q, r);
        
        XTLog(std::cout) << "r_deg:" << r_deg << std::endl;
        XTLog(std::cout) << "q:" << q << std::endl;
        XTLog(std::cout) << "r:" << r << std::endl;
    }

    XTLog(std::cout) << "-----------------------" << std::endl;

    {
        xiaotu::Polynomial<Complex> a({1.0, 2.0, 3.0});
        xiaotu::Polynomial<Complex> b({1.0, 2.0});

        auto q = xiaotu::Polynomial<Complex>::Zero();
        auto r = xiaotu::Polynomial<Complex>::Zero();
        auto r_deg = a.Divide(b, q, r);
        
        XTLog(std::cout) << "r_deg:" << r_deg << std::endl;
        XTLog(std::cout) << "q:" << q << std::endl;
        XTLog(std::cout) << "r:" << r << std::endl;

        auto aa = q * b + r;
        XTLog(std::cout) << "a:" << a << std::endl;
        XTLog(std::cout) << "aa:" << aa << std::endl;
    }

    XTLog(std::cout) << "-----------------------" << std::endl;

    {
        xiaotu::Polynomial<Complex> a({1.0, 2.0, 3.0});

        auto q = xiaotu::Polynomial<Complex>::Zero();
        Complex b0;
        a.SyntheticDivide(2.0, q, b0);
        
        XTLog(std::cout) << "a:" << a << std::endl;
        XTLog(std::cout) << "q:" << q << std::endl;
        XTLog(std::cout) << "r:" << b0 << std::endl;

        auto aa = q * xiaotu::Polynomial<Complex>({-2.0, 1}) + b0;
        XTLog(std::cout) << "aa:" << aa << std::endl;
    }
}


TEST(ComplexPolynomial, Equal)
{
    using Complex = std::complex<double>;
    {
        xiaotu::Polynomial<Complex> a({1.0, 2.0, 3.0});
        xiaotu::Polynomial<Complex> b({1.0, 2.0, 3.0});

        EXPECT_TRUE(a == b);
        EXPECT_EQ(a, b);
    }

    {
        xiaotu::Polynomial<Complex> a({1.0, 2.0, 3.0});
        xiaotu::Polynomial<Complex> b({4.0, 2.0, 3.0});

        EXPECT_TRUE(a != b);
        EXPECT_NE(a, b);
    }

    {
        auto zero = xiaotu::Polynomial<Complex>::Zero();
        EXPECT_TRUE(0 == zero);
        EXPECT_EQ(0, zero);
        EXPECT_TRUE(zero == 0);
        EXPECT_EQ(zero, 0);
    }

    {
        auto one = xiaotu::Polynomial<Complex>::One();
        EXPECT_TRUE(1 == one);
        EXPECT_EQ(1, one);
        EXPECT_TRUE(one == 1);
        EXPECT_EQ(one, 1);
    }

    {
        xiaotu::Polynomial<Complex> a({1.0, 2.0, 3.0});
        EXPECT_FALSE(1 == a);
        EXPECT_NE(1, a);
        EXPECT_FALSE(a == 1);
        EXPECT_NE(a, 1);
    }
}


TEST(ComplexPolynomial, AllRoots)
{
    using Complex = std::complex<double>;
    {
        xiaotu::Polynomial<Complex> a({-1.0, 0.0, 0.0, 1.0});
        std::vector<Complex> re;

        a.AllRoots(re);
        EXPECT_EQ(a.Degree(), re.size());

        for (size_t i = 0; i < re.size(); i++) {
            EXPECT_TRUE(std::abs(a(re[i])) < SMALL_VALUE);
            XTLog(std::cout) << "[" << i << "] " << re[i] << std::endl;
        }
    }

    {
        xiaotu::Polynomial<Complex> a({-8, 4, 6, -5, 1});
        std::vector<Complex> re;

        a.AllRoots(re);
        EXPECT_EQ(a.Degree(), re.size());

        for (size_t i = 0; i < re.size(); i++) {
            EXPECT_TRUE(std::abs(a(re[i])) < SMALL_VALUE);
            XTLog(std::cout) << "[" << i << "] " << re[i] << std::endl;
        }
    }

    Complex ha(1.0, 1.0);
    std::cout << std::abs(ha) << std::endl;
}
