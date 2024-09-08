#include <iostream>

#include <XiaoTuMathBox/Utils.hpp>
#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

#include <gtest/gtest.h>

#include <memory>
#include <vector>
#include <cmath>

using namespace xiaotu::math;


TEST(LinearAlgibra, GaussJordanEliminate)
{
    Matrix<double, 3, 3> A = {
        9, -3, 1,
        1,  1, 1,
        4,  2, 1
    };
    Matrix<double, 3, 1> b = { 20, 0, 10 };

    GaussJordanEliminate(A, &b);
    XTLog(std::cout) << "A = " << A << std::endl;
    XTLog(std::cout) << "b = " << b << std::endl;

    EXPECT_TRUE(std::abs(b(0) -   3)  < 1e-9);
    EXPECT_TRUE(std::abs(b(1) -   1)  < 1e-9);
    EXPECT_TRUE(std::abs(b(2) - (-4)) < 1e-9);
}

TEST(LinearAlgibra, LU)
{
    Matrix<double, 3, 3> A = {
        9, -3, 1,
        1,  1, 1,
        4,  2, 1
    };

    LU lu(A.View());
    XTLog(std::cout) << "lu = " << lu() << std::endl;

    Matrix<double, 3, 1> b = { 20, 0, 10 };
    Matrix<double, 3, 1> x;

    lu.Solve(b, x);
    XTLog(std::cout) << "x = " << x << std::endl;
    XTLog(std::cout) << "b = " << b << std::endl;

    Matrix<double, 3, 3> Ainv;
    lu.Inverse(Ainv.View());
    XTLog(std::cout) << "A^{-1} = " << Ainv << std::endl;

    Matrix<double, 3, 1> c = Ainv * b;
    XTLog(std::cout) << "c = " << c << std::endl;
}

TEST(LinearAlgibra, douniwan)
{
    Matrix<double, 3, 1> b = { 20, 0, 10 };
    Matrix<double, 1, 3> bT = { 20, 0, 10 };

    auto bbt = b * bT;
    XTLog(std::cout) << "bbt = " << bbt << std::endl;

    auto btb = bT * b;
    XTLog(std::cout) << "btb = " << btb << std::endl;

    Matrix<double, 3, 3> A;
    A.Identity();

    XTLog(std::cout) << "btA = " << bT * A << std::endl;
    XTLog(std::cout) << "btAb = " << bT * A * b << std::endl;

}

