#include <iostream>

#include <XiaoTuMathBox/Utils.hpp>
#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

#include <gtest/gtest.h>

#include <memory>
#include <vector>
#include <cmath>

using namespace xiaotu::math;

TEST(QR, GramSchmidt)
{
    Matrix<double, 3, 3> A = {
        9, -3, 1,
        1,  1, 1,
        4,  2, 1
    };

    QR_GramSchmidt qr(A);
    XTLog(std::cout) << "Q = " << qr.Q() << std::endl;
    XTLog(std::cout) << "R = " << qr.R() << std::endl;

    auto ortho_0 = qr.Q().SubMatrix(0, 0, 3, 1);
    auto ortho_1 = qr.Q().SubMatrix(0, 1, 3, 1);
    auto ortho_2 = qr.Q().SubMatrix(0, 2, 3, 1);

    EXPECT_TRUE(std::abs(ortho_0.Dot(ortho_1)) < 1e-9);
    EXPECT_TRUE(std::abs(ortho_0.Dot(ortho_2)) < 1e-9);
    EXPECT_TRUE(std::abs(ortho_1.Dot(ortho_2)) < 1e-9);

    auto eye = qr.Q().Transpose() * qr.Q();
    XTLog(std::cout) << "eye = " << eye << std::endl;

    auto a = qr.Q() * qr.R();
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            EXPECT_TRUE(std::abs(a(i, j) - A(i, j)) < 1e-9);
}

TEST(QR, HouseholderMatrix)
{
    Vector<double, 3> x = { 1, 2, 3 };
    Matrix<double, 3, 3> H;
    double x_norm = x.Norm();

    {
        auto v = HouseholderVector(x);
        HouseholderMatrix(v, H);
        auto y = H * x;

        EXPECT_DOUBLE_EQ(x_norm, std::abs(y(0)));
        EXPECT_TRUE(std::abs(y(1)) < 1e-12);
        EXPECT_TRUE(std::abs(y(2)) < 1e-12);
    }

    {
        auto v = HouseholderVector(x, 1);
        HouseholderMatrix(v, H);
        auto y = H * x;

        EXPECT_TRUE(std::abs(y(0)) < 1e-12);
        EXPECT_DOUBLE_EQ(x_norm, std::abs(y(1)));
        EXPECT_TRUE(std::abs(y(2)) < 1e-12);
    }

    {
        auto v = HouseholderVector(x, 2);
        HouseholderMatrix(v, H);
        auto y = H * x;

        EXPECT_TRUE(std::abs(y(0)) < 1e-12);
        EXPECT_TRUE(std::abs(y(1)) < 1e-12);
        EXPECT_DOUBLE_EQ(x_norm, std::abs(y(2)));
    }
}

TEST(QR, Householder)
{
    Matrix<double, 3, 4> A = {
        9, -3, 1, 4,
        1,  1, 1, 4,
        4,  2, 1, 4
    };

    QR_Householder qr(A);
    XTLog(std::cout) << "Q = " << qr.Q() << std::endl;
    XTLog(std::cout) << "R = " << qr.R() << std::endl;

    EXPECT_TRUE(qr.R()(0, 0) > 0);
    EXPECT_TRUE(qr.R()(1, 1) > 0);
    EXPECT_TRUE(qr.R()(2, 2) > 0);

    auto ortho_0 = qr.Q().SubMatrix(0, 0, 3, 1);
    auto ortho_1 = qr.Q().SubMatrix(0, 1, 3, 1);
    auto ortho_2 = qr.Q().SubMatrix(0, 2, 3, 1);

    EXPECT_TRUE(std::abs(ortho_0.Dot(ortho_1)) < 1e-9);
    EXPECT_TRUE(std::abs(ortho_0.Dot(ortho_2)) < 1e-9);
    EXPECT_TRUE(std::abs(ortho_1.Dot(ortho_2)) < 1e-9);

    auto eye = qr.Q().Transpose() * qr.Q();
    XTLog(std::cout) << "eye = " << eye << std::endl;

    auto a = qr.Q() * qr.R();
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++)
            EXPECT_TRUE(std::abs(a(i, j) - A(i, j)) < 1e-9);
    XTLog(std::cout) << "a = " << a << std::endl;
}

TEST(QR, UpperHessenberg)
{
    Matrix<double, 4, 4> A = {
        9, -3, 1, 4,
        9,  1, 1, 4,
        4,  2, 1, 4,
        4,  2, 9, 7
    };

    UpperHessenberg h(A);
    XTLog(std::cout) << "Q = " << h.Q() << std::endl;
    XTLog(std::cout) << "H = " << h.H() << std::endl;

    auto ha = h.Q().Transpose() * h.H() * h.Q();
    XTLog(std::cout) << "ha = " << ha << std::endl;
}

TEST(QR, Givens)
{
    Vector<double, 3> x = { 1, 2, 3 };
    {
        auto G = Givens<double>(0, 1, x);
        auto v = G * x;

        EXPECT_TRUE(std::abs(v(1)) < SMALL_VALUE);
        EXPECT_DOUBLE_EQ(3.0, v(2));
        XTLog(std::cout) << x << std::endl;
    }

    {
        EXPECT_DOUBLE_EQ(1.0, x(0));
        EXPECT_DOUBLE_EQ(2.0, x(1));
        EXPECT_DOUBLE_EQ(3.0, x(2));

        auto G = Givens<double>(0, 1, x);
        G.LeftApplyOn(x);

        EXPECT_TRUE(std::abs(x(1)) < SMALL_VALUE);
        EXPECT_DOUBLE_EQ(3.0, x(2));
        XTLog(std::cout) << x << std::endl;
    }

    {
        x = { 1, 2, 3 };
        RowVector<double, 3> y = { 1, 2, 3 };
        auto G = Givens<double>(0, 1, x);
        XTLog(std::cout) << "c:" << G.c() << std::endl;
        XTLog(std::cout) << "s:" << G.s() << std::endl;

        G.TRightApplyOn(y);

        EXPECT_TRUE(std::abs(y(1)) < SMALL_VALUE);
        EXPECT_DOUBLE_EQ(3.0, y(2));
        XTLog(std::cout) << y << std::endl;
    }
}

TEST(QR, GivensMatrix)
{
    Vector<double, 3> x = { 1, 2, 3 };

    {
        auto G = GivensMatrix(0, 1, x);
        auto y = G * x;

        EXPECT_TRUE(std::abs(y(1)) < SMALL_VALUE);
        EXPECT_DOUBLE_EQ(3.0, y(2));
        XTLog(std::cout) << y << std::endl;
    }

    {
        auto G = GivensMatrix(2, 1, x);
        auto y = G * x;

        EXPECT_TRUE(std::abs(y(1)) < SMALL_VALUE);
        EXPECT_DOUBLE_EQ(1.0, y(0));
        XTLog(std::cout) << y << std::endl;
    }


    {
        RowVector<double, 3> y = { 1, 2, 3 };
        auto G = Givens<double>(0, 1, y);
        XTLog(std::cout) << "y * G:" << y * G << std::endl;

        G.TRightApplyOn(y);

        EXPECT_TRUE(std::abs(y(1)) < SMALL_VALUE);
        EXPECT_DOUBLE_EQ(3.0, y(2));
        XTLog(std::cout) << y << std::endl;
    }
}






