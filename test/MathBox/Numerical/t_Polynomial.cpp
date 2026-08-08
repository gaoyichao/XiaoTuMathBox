#include <iostream>

#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Numerical/Polynomial.hpp>

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
        XTLog(std::cout) << "P(1.0) = " << p << "P'(1.0) = " << dp << std::endl;
    }

    {
        double p, dp;
        P.Horner(2.0, p, dp);
        EXPECT_DOUBLE_EQ(P(2.0), p);
        EXPECT_DOUBLE_EQ(DP(2.0), dp);
        XTLog(std::cout) << "P(2.0) = " << p << "P'(2.0) = " << dp << std::endl;
    }

}




