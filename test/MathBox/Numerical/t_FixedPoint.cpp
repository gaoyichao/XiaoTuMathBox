#include <iostream>

#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Numerical/FixedPoint.hpp>

#include <cmath>
#include <gtest/gtest.h>



TEST(FixedPoint, Naive)
{
    {
        double root = xiaotu::NaiveFixedPoint<double>(
        [](double x) {
            return std::sqrt(0.5 * (x*x*x - x + 2));
        }, 1.01);
        EXPECT_TRUE(std::abs(root - 1.0) < SMALL_VALUE);
    }

}

TEST(FixedPoint, AitkenSteffensen)
{
    {
        double root = xiaotu::AitkenSteffensen<double>(
        [](double x) {
            return std::sqrt(0.5 * (x*x*x - x + 2));
        }, 1.01);
        EXPECT_TRUE(std::abs(root - 1.0) < SMALL_VALUE);
    }

}



