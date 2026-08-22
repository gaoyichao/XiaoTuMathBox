#include <iostream>

#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Numerical/Numerical.hpp>

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
        EXPECT_TRUE(std::abs(root - 2) < SMALL_VALUE);
    }

    {
        auto func = [](double x) {
            return std::pow(2, -x) - x;
        };
        double root = xiaotu::Bisection<double>(func, 0.0, 1.0);
        EXPECT_TRUE(std::abs(func(root)) < SMALL_VALUE);
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
        auto func = [](double x) {
            return std::cos(x) - x;
        };
        auto dfunc = [](double x) {
            return -std::sin(x) - 1;
        };

        double root = xiaotu::NewtonRaphson<double>(func, dfunc, 0.0, 100, 1e-15);
        EXPECT_DOUBLE_EQ(func(root), 0);
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

TEST(EquationRoot, BrentRoot)
{
    {
        double root = xiaotu::BrentRoot<double>(
        [](double x) {
            return 0.5 * x;
        }, -1.0, 1.0);
        EXPECT_TRUE(std::abs(root - 0) < SMALL_VALUE);
    }

    {
        double root = xiaotu::BrentRoot<double>(
        [](double x) {
            return x * x - 4;
        }, 0.0, 3.0, 100, 1e-15);
        EXPECT_DOUBLE_EQ(root, 2);
    }


    {
        double root = xiaotu::BrentRoot<double>(
        [](double x) {
            return x * x - 10000.0001 * x + 1;
        }, -1.0, 1.0, 100, 1e-15);
        EXPECT_DOUBLE_EQ(root, 0.0001);
    }

    {
        double root = xiaotu::BrentRoot<double>(
        [](double x) {
            return x * x - 10000.0001 * x + 1;
        }, 10001.0, 9999.0, 100, 1e-15);
        EXPECT_DOUBLE_EQ(root, 10000);
    }

    {
        auto func = [](double x) {
            return std::pow(2, -x) - x;
        };
        double root = xiaotu::BrentRoot<double>(func, 0.0, 1.0);
        EXPECT_DOUBLE_EQ(func(root), 0);
    }

}

TEST(EquationRoot, QuadraticRoot)
{
    std::complex<double> x0, x1;

    {
        bool re = xiaotu::QuadraticRoot<double>(1, -10000.0001, 1, x0, x1);
        EXPECT_TRUE(re);
        EXPECT_TRUE(std::abs(x0 - 10000.0) < SMALL_VALUE || std::abs(x0 - 0.0001) < SMALL_VALUE);
        EXPECT_TRUE(std::abs(x1 - 10000.0) < SMALL_VALUE || std::abs(x1 - 0.0001) < SMALL_VALUE);
        EXPECT_TRUE(x0 != x1);
    }

    {
        bool re = xiaotu::QuadraticRoot<double>(1, 0, 1, x0, x1);
        EXPECT_FALSE(re);
        EXPECT_TRUE(x0.real() == 0);
        EXPECT_TRUE(x1.real() == 0);

        EXPECT_TRUE(std::abs(x0.imag() - 1.0) < SMALL_VALUE || std::abs(x0.imag() + 1.0) < SMALL_VALUE);
        EXPECT_TRUE(std::abs(x1.imag() - 1.0) < SMALL_VALUE || std::abs(x1.imag() + 1.0) < SMALL_VALUE);
    }
}

TEST(EquationRoot, MullerRoot)
{
    {
        double root = xiaotu::MullerRoot<double>(
        [](double x) {
            return 0.5 * x;
        }, -1.0, 1.0, 2.0);
        EXPECT_TRUE(std::abs(root - 0) < SMALL_VALUE);
    }

    {
        double root = xiaotu::MullerRoot<double>(
        [](double x) {
            return x * x - 4;
        }, 0.0, 3.0, 1.0, 100, 1e-15);
        EXPECT_DOUBLE_EQ(root, 2);
    }


    {
        double root = xiaotu::MullerRoot<double>(
        [](double x) {
            return x * x - 10000.0001 * x + 1;
        }, -1.0, 1.0, 2.0, 100, 1e-15);
        EXPECT_DOUBLE_EQ(root, 0.0001);
    }

    {
        double root = xiaotu::MullerRoot<double>(
        [](double x) {
            return x * x - 10000.0001 * x + 1;
        }, 10001.0, 9999.0, 9990.0, 100, 1e-15);
        EXPECT_DOUBLE_EQ(root, 10000);
    }


    {
        auto func = [](double x) {
            return std::pow(2, -x) - x;
        };
        double root = xiaotu::MullerRoot<double>(func, 0.0, 1.0, 2.0);
        EXPECT_DOUBLE_EQ(func(root), 0);
    }
}

TEST(EquationRoot, ComplexMullerRoot)
{
    {
        using Complex = std::complex<double>;
        Complex root = xiaotu::MullerRoot<double>(
            [](Complex x) {
                return 0.5 * x;
            },
            Complex(-1.0, 0), Complex(1.0, 0), Complex(2.0, 0)
        );
        EXPECT_TRUE(std::abs(root - 0.0) < SMALL_VALUE);
    }

    {
        using Complex = std::complex<double>;
        auto func = [](Complex x) {
            return x * x + 1.0;
        };
        Complex root = xiaotu::MullerRoot<double>(func,
            Complex(-1.0, 0), Complex(1.0, 0), Complex(2.0, 0)
        );
        EXPECT_TRUE(std::abs(func(root)) < SMALL_VALUE);
    }

}


