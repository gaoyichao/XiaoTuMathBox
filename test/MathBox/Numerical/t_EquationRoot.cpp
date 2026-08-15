#include <iostream>

#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Numerical/EquationRoot.hpp>

#include <cmath>
#include <gtest/gtest.h>


TEST(EquationRoot, Bisection)
{
    {
        double root = xiaotu::Bisection<double>(
        [](double x) {
            return 0.5 * x;
        }, -1.0, 1.0);
        EXPECT_TRUE(std::abs(root - 0) < SMALL_VALUE);
    }

    {
        double root = xiaotu::Bisection<double>(
        [](double x) {
            return x * x - 4;
        }, 0.0, 3.0, 100, 1e-15);
        EXPECT_DOUBLE_EQ(root, 2);
    }

}

TEST(EquationRoot, NewtonRaphson)
{
    {
        double root = xiaotu::NewtonRaphson<double>(
        [](double x) {
            return 0.5 * x;
        }, 
        [](double x) {
            return 0.5;
        }, 
        1.0);
        EXPECT_TRUE(std::abs(root - 0) < SMALL_VALUE);
    }

    {
        double root = xiaotu::NewtonRaphson<double>(
        [](double x) {
            return x * x - 4;
        }, 
        [](double x) {
            return 2 * x;
        }, 
        1.0, 100, 1e-15);
        EXPECT_DOUBLE_EQ(root, 2);
    }


    {
        double root = xiaotu::NewtonRaphson<double>(
        [](double x) {
            return std::cos(x) - x;
        }, 
        [](double x) {
            return -std::sin(x) - 1;
        }, 
        0.0, 100, 1e-15);
        EXPECT_TRUE(std::abs(root - 0.7390851332) < 1e-9);
    }


}

TEST(EquationRoot, SecantRoot)
{
    {
        double root = xiaotu::SecantRoot<double>(
        [](double x) {
            return 0.5 * x;
        }, -1.0, 1.0);
        EXPECT_TRUE(std::abs(root - 0) < SMALL_VALUE);
    }

    {
        double root = xiaotu::SecantRoot<double>(
        [](double x) {
            return x * x - 4;
        }, 0.0, 3.0, 100, 1e-15);
        EXPECT_DOUBLE_EQ(root, 2);
    }

    {
        double root = xiaotu::SecantRoot<double>(
        [](double x) {
            return x * x - 10000.0001 * x + 1;
        }, -1.0, 1.0, 100, 1e-15);
        EXPECT_DOUBLE_EQ(root, 0.0001);
    }

    {
        double root = xiaotu::SecantRoot<double>(
        [](double x) {
            return x * x - 10000.0001 * x + 1;
        }, 10001.0, 9999.0, 100, 1e-15);
        EXPECT_DOUBLE_EQ(root, 10000);
    }

}


TEST(EquationRoot, FalsePosition)
{
    {
        double root = xiaotu::FalsePosition<double>(
        [](double x) {
            return 0.5 * x;
        }, -1.0, 1.0);
        EXPECT_TRUE(std::abs(root - 0) < SMALL_VALUE);
    }

    {
        double root = xiaotu::FalsePosition<double>(
        [](double x) {
            return x * x - 4;
        }, 0.0, 3.0, 100, 1e-15);
        EXPECT_DOUBLE_EQ(root, 2);
    }


    {
        double root = xiaotu::FalsePosition<double>(
        [](double x) {
            return x * x - 10000.0001 * x + 1;
        }, -1.0, 1.0, 100, 1e-15);
        EXPECT_DOUBLE_EQ(root, 0.0001);
    }

    {
        double root = xiaotu::FalsePosition<double>(
        [](double x) {
            return x * x - 10000.0001 * x + 1;
        }, 10001.0, 9999.0, 100, 1e-15);
        EXPECT_DOUBLE_EQ(root, 10000);
    }

}


TEST(EquationRoot, DekkerRoot)
{
    {
        double root = xiaotu::DekkerRoot<double>(
        [](double x) {
            return 0.5 * x;
        }, -1.0, 1.0);
        EXPECT_TRUE(std::abs(root - 0) < SMALL_VALUE);
    }

    {
        double root = xiaotu::DekkerRoot<double>(
        [](double x) {
            return x * x - 4;
        }, 0.0, 3.0, 100, 1e-15);
        EXPECT_DOUBLE_EQ(root, 2);
    }


    {
        double root = xiaotu::DekkerRoot<double>(
        [](double x) {
            return x * x - 10000.0001 * x + 1;
        }, -1.0, 1.0, 100, 1e-15);
        EXPECT_DOUBLE_EQ(root, 0.0001);
    }

    {
        double root = xiaotu::DekkerRoot<double>(
        [](double x) {
            return x * x - 10000.0001 * x + 1;
        }, 10001.0, 9999.0, 100, 1e-15);
        EXPECT_DOUBLE_EQ(root, 10000);
    }

}


