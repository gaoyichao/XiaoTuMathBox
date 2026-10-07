#include <iostream>

#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Numerical/Numerical.hpp>

#include <gtest/gtest.h>


TEST(NumericalDifference, ForwardBackward)
{
    auto fn = [](double x) {
        return std::log(x);
    };


    {
        auto f1 = xiaotu::ForwardDifference<double>(fn, 1.0, 0.1);
        auto f2 = xiaotu::BackwardDifference<double>(fn, 1.0, -0.1);
        EXPECT_DOUBLE_EQ(f1, f2);
    }

    {
        auto f1 = xiaotu::ForwardDifference<double>(fn, 1.0, -0.1);
        auto f2 = xiaotu::BackwardDifference<double>(fn, 1.0, 0.1);
        EXPECT_DOUBLE_EQ(f1, f2);
    }
}



