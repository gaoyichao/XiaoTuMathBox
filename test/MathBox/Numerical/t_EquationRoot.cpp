#include <iostream>

#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Numerical/EquationRoot.hpp>

#include <gtest/gtest.h>

double douniwan(double x)
{
    return 0.5 * x;
}


TEST(EquationRoot, Bisection)
{
    double root = xiaotu::Bisection<double>(douniwan, -1.0, 1.0);
    XTLog(std::cout) << root << std::endl;

}

