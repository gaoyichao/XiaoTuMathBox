#include <iostream>

#include <XiaoTuMathBox/Utils.hpp>
#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

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

    {
        Bidiagonal bidiag(A);
        EXPECT_TRUE(bidiag.B().IsUpperBiDiagonal());

        XTLog(std::cout) << "B = " << bidiag.B() << std::endl;

        auto _A_ = bidiag.UT().Transpose() *
                      bidiag.B() *
                      bidiag.V().Transpose();
        EXPECT_TRUE(_A_ == A);

        XTLog(std::cout) << "_A_ = " << _A_ << std::endl;
    }

    XTLog(std::cout) << "__A__ = " << A.UpperBidiagonal() << std::endl;
}

TEST(SVD, LowerBidiagonal)
{
    Matrix<double, 4, 7> A = {
        2, 1, 1, 1, 1, 1, 1, 
        1, 2, 1, 1, 1, 1, 1,
        1, 1, 2, 1, 1, 1, 1,
        1, 1, 1, 2, 1, 1, 1,
    };

    {
        Bidiagonal bidiag(A);
        EXPECT_TRUE(bidiag.B().IsLowerBiDiagonal());

        XTLog(std::cout) << "B = " << bidiag.B() << std::endl;

        auto _A_ = bidiag.UT().Transpose() *
                      bidiag.B() *
                      bidiag.V().Transpose();
        EXPECT_TRUE(_A_ == A);

        XTLog(std::cout) << "_A_ = " << _A_ << std::endl;
    }

    XTLog(std::cout) << "__A__ = " << A.LowerBidiagonal() << std::endl;
}



TEST(SVD, SVD_Naive)
{
    Matrix<double, 7, 5> A = {
        2, 1, 1, 1, 3,
        1, 2, 1, 1, 3,
        1, 1, 2, 1, 3,
        1, 1, 1, 2, 3,
        1, 1, 1, 2, 3,
        1, 2, 1, 2, 3,
        2, 2, 1, 2, 3,
    };

    {
        SVD_Naive svd(A, 1000, SMALL_VALUE);
        auto _sigma_ = svd.Sigma();
        XTLog(std::cout) << "_Sigma_ = " << _sigma_.Truncate() << std::endl;
        
        auto ha = svd.UT() * A * svd.V();
        XTLog(std::cout) << "_ha_ = " << ha.Truncate() << std::endl;

        auto _A_ = svd.UT().Transpose() * svd.Sigma() * svd.V().Transpose();
        XTLog(std::cout) << "_A_ = " << _A_.Truncate() << std::endl;
    }
}

TEST(SVD, SVD_GKR)
{
    Matrix<double, 7, 5> A = {
        2, 1, 1, 1, 3,
        1, 2, 1, 1, 3,
        1, 1, 2, 1, 3,
        1, 1, 1, 2, 3,
        1, 1, 1, 2, 3,
        1, 2, 1, 2, 3,
        2, 2, 1, 2, 3,
    };

    {
        SVD_GKR svd(A, true, true);
        int n = svd.Iterate(1000, SMALL_VALUE);
        XTLog(std::cout) << "迭代次数:" << n << std::endl;

        auto _sigma_ = svd.Sigma();
        XTLog(std::cout) << "_Sigma_ = " << _sigma_.Truncate() << std::endl;
        
        auto ha = svd.UT() * A * svd.V();
        XTLog(std::cout) << "_ha_ = " << ha.Truncate() << std::endl;

        auto _A_ = svd.UT().Transpose() * svd.Sigma() * svd.V().Transpose();
        XTLog(std::cout) << "_A_ = " << _A_.Truncate() << std::endl;
    }
}


TEST(SVD, SVD_LowerGKR)
{
    Matrix<double, 5, 7> A = {
        2, 1, 1, 1, 3, 9, 1,
        1, 2, 1, 1, 3, 9, 1,
        1, 1, 2, 1, 3, 9, 1,
        1, 1, 1, 2, 3, 9, 1,
        1, 1, 1, 2, 3, 9, 1,
    };
    {
        auto _at_ = A.Transpose();
        XTLog(std::cout) << "_at_:" << _at_ << std::endl;

        SVD_GKR svd(_at_, true, true);
        int n = svd.Iterate(1000, SMALL_VALUE);
        XTLog(std::cout) << "迭代次数:" << n << std::endl;

        auto _sigma_ = svd.Sigma();
        XTLog(std::cout) << "_Sigma_ = " << _sigma_.Truncate() << std::endl;
        
    }

    {
        SVD_GKR svd(A, true, true);
        int n = svd.Iterate(1000, SMALL_VALUE);
        XTLog(std::cout) << "迭代次数:" << n << std::endl;

        auto _sigma_ = svd.Sigma();
        XTLog(std::cout) << "_Sigma_ = " << _sigma_.Truncate() << std::endl;
        
        auto ha = svd.UT() * A * svd.V();
        XTLog(std::cout) << "_ha_ = " << ha.Truncate() << std::endl;

        auto _A_ = svd.UT().Transpose() * svd.Sigma() * svd.V().Transpose();
        XTLog(std::cout) << "_A_ = " << _A_.Truncate() << std::endl;
    }
}



