#include <iostream>

#include <XiaoTuMathBox/Utils.hpp>
#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>
#include <XiaoTuMathBox/LinearAlgibra/Bidiagonal.hpp>

#include <gtest/gtest.h>

#include <memory>
#include <vector>
#include <cmath>

using namespace xiaotu::math;

TEST(SVD, UpperBidiagonal)
{
    Matrix<double, 7, 4> A = {
        2, 1, 1, 1,
        1, 2, 1, 1,
        1, 1, 2, 1,
        1, 1, 1, 2,
        1, 1, 1, 2,
        1, 2, 1, 2,
        2, 2, 1, 2,
    };

    Bidiagonal bidiag(A);
    EXPECT_TRUE(bidiag.B().IsUpperBiDiagonal());

    XTLog(std::cout) << "B = " << bidiag.B() << std::endl;

    auto _A_ = bidiag.UT().Transpose() *
                  bidiag.B() *
                  bidiag.V().Transpose();
    EXPECT_TRUE(_A_ == A);

    XTLog(std::cout) << "_A_ = " << _A_ << std::endl;
}

TEST(SVD, LowerBidiagonal)
{
    Matrix<double, 4, 7> A = {
        2, 1, 1, 1, 1, 1, 1, 
        1, 2, 1, 1, 1, 1, 1,
        1, 1, 2, 1, 1, 1, 1,
        1, 1, 1, 2, 1, 1, 1,
    };

    Bidiagonal bidiag(A);
    EXPECT_TRUE(bidiag.B().IsLowerBiDiagonal());

    XTLog(std::cout) << "B = " << bidiag.B() << std::endl;

    auto _A_ = bidiag.UT().Transpose() *
                  bidiag.B() *
                  bidiag.V().Transpose();
    EXPECT_TRUE(_A_ == A);

    XTLog(std::cout) << "_A_ = " << _A_ << std::endl;
}

