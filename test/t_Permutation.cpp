#include <iostream>

#include <XiaoTuMathBox/Utils.hpp>
#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

#include <gtest/gtest.h>

#include <memory>
#include <vector>
#include <cmath>

using namespace xiaotu::math;


TEST(Permutation, RowMatrix)
{
    Matrix<double, 3, 2> A = {
        0, 3,
        1, 4,
        2, 5
    };

    Permutation p({2, 1, 0});
    auto pA = p.LeftMatrix<double>() * A;

    EXPECT_DOUBLE_EQ(pA(0, 0), A(2, 0));
    EXPECT_DOUBLE_EQ(pA(1, 0), A(1, 0));
    EXPECT_DOUBLE_EQ(pA(2, 0), A(0, 0));

    EXPECT_DOUBLE_EQ(pA(0, 1), A(2, 1));
    EXPECT_DOUBLE_EQ(pA(1, 1), A(1, 1));
    EXPECT_DOUBLE_EQ(pA(2, 1), A(0, 1));

    {
        auto ppA = p * A;
        EXPECT_DOUBLE_EQ(pA(0, 0), ppA(0, 0));
        EXPECT_DOUBLE_EQ(pA(1, 0), ppA(1, 0));
        EXPECT_DOUBLE_EQ(pA(2, 0), ppA(2, 0));

        EXPECT_DOUBLE_EQ(pA(0, 1), ppA(0, 1));
        EXPECT_DOUBLE_EQ(pA(1, 1), ppA(1, 1));
        EXPECT_DOUBLE_EQ(pA(2, 1), ppA(2, 1));
    }

    p.LeftApplyOn(A);

    EXPECT_DOUBLE_EQ(pA(0, 0), A(0, 0));
    EXPECT_DOUBLE_EQ(pA(1, 0), A(1, 0));
    EXPECT_DOUBLE_EQ(pA(2, 0), A(2, 0));

    EXPECT_DOUBLE_EQ(pA(0, 1), A(0, 1));
    EXPECT_DOUBLE_EQ(pA(1, 1), A(1, 1));
    EXPECT_DOUBLE_EQ(pA(2, 1), A(2, 1));

    XTLog(std::cout) << A << std::endl;
    XTLog(std::cout) << A.Rows() << std::endl;
}

TEST(Permutation, ColMatrix)
{
    Matrix<double, 2, 3> A = {
        0, 1, 2,
        3, 4, 5
    };

    Permutation p({2, 1, 0});
    auto Ap = A * p.RightMatrix<double>();

    EXPECT_DOUBLE_EQ(Ap(0, 0), A(0, 2));
    EXPECT_DOUBLE_EQ(Ap(0, 1), A(0, 1));
    EXPECT_DOUBLE_EQ(Ap(0, 2), A(0, 0));

    EXPECT_DOUBLE_EQ(Ap(1, 0), A(1, 2));
    EXPECT_DOUBLE_EQ(Ap(1, 1), A(1, 1));
    EXPECT_DOUBLE_EQ(Ap(1, 2), A(1, 0));

    {
        auto App = A * p;
        EXPECT_DOUBLE_EQ(Ap(0, 0), App(0, 0));
        EXPECT_DOUBLE_EQ(Ap(0, 1), App(0, 1));
        EXPECT_DOUBLE_EQ(Ap(0, 2), App(0, 2));

        EXPECT_DOUBLE_EQ(Ap(1, 0), App(1, 0));
        EXPECT_DOUBLE_EQ(Ap(1, 1), App(1, 1));
        EXPECT_DOUBLE_EQ(Ap(1, 2), App(1, 2));
    }

    p.RightApplyOn(A);

    EXPECT_DOUBLE_EQ(Ap(0, 0), A(0, 0));
    EXPECT_DOUBLE_EQ(Ap(0, 1), A(0, 1));
    EXPECT_DOUBLE_EQ(Ap(0, 2), A(0, 2));

    EXPECT_DOUBLE_EQ(Ap(1, 0), A(1, 0));
    EXPECT_DOUBLE_EQ(Ap(1, 1), A(1, 1));
    EXPECT_DOUBLE_EQ(Ap(1, 2), A(1, 2));

    XTLog(std::cout) << A << std::endl;
    XTLog(std::cout) << A.Rows() << std::endl;
}

