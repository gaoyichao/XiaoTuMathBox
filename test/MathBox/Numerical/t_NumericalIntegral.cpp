#include <iostream>

#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Numerical/Numerical.hpp>

#include <gtest/gtest.h>


TEST(GaussianIntegration, LegendrePolynomial)
{
    {
        auto P0 = xiaotu::LegendrePolynomial<double>(0);
        XTLog(std::cout) << P0 << std::endl;

        auto P1 = xiaotu::LegendrePolynomial<double>(1);
        XTLog(std::cout) << P1 << std::endl;

        auto P2 = xiaotu::LegendrePolynomial<double>(2);
        XTLog(std::cout) << P2 << std::endl;

        auto P3 = xiaotu::LegendrePolynomial<double>(3);
        XTLog(std::cout) << P3 << std::endl;

        auto P4 = xiaotu::LegendrePolynomial<double>(4);
        XTLog(std::cout) << P4 << std::endl;
    }
}

TEST(GaussianIntegration, init)
{
    xiaotu::GaussianIntegration<double> gaussian(3);

    XTLog(std::cout) << gaussian << std::endl;
}

