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


TEST(QR, EigenNaiveQR)
{
    Matrix<double, 3, 3> A = {
        2, 1, 1,
        1, 2, 1,
        1, 1, 2
    };

    EigenNaiveQR qr(A, 20);
    auto eigen_values = qr.EigenValues();
    std::sort(eigen_values.begin(), eigen_values.end(), [](double a, double b){ return a > b; });
    EXPECT_TRUE(std::abs(eigen_values[0] - 4.0) < 1e-9);
    EXPECT_TRUE(std::abs(eigen_values[1] - 1.0) < 1e-9);
    EXPECT_TRUE(std::abs(eigen_values[2] - 1.0) < 1e-9);

    XTLog(std::cout) << "T = " << qr.T() << std::endl;
}

TEST(QR, EigenShiftQR)
{
    Matrix<double, 5, 5> A = {
        3, 1, 1, 1, 1,
        1, 2, 1, 1, 1,
        1, 1, 2, 1, 1,
        1, 1, 1, 2, 1,
        1, 1, 1, 1, 3,
    };

    {
        EigenShiftQR<Matrix<double, 5, 5>> qr;
        int n = qr.NaiveIterate(A, 1000, SMALL_VALUE);
        XTLog(std::cout) << "朴素迭代次数: " << n << std::endl;

        auto eigen_values = qr.EigenValues();
        std::sort(eigen_values.begin(), eigen_values.end(), [](double a, double b){ return a > b; });
        XTLog(std::cout) << "T = " << qr.T() << std::endl;

        MatrixView<double, 5, 1> lambdas(eigen_values.data());
        XTLog(std::cout) << "lambdas = " << lambdas << std::endl;
    }

    {
        UpperHessenberg h(A);

        EigenShiftQR<Matrix<double, 5, 5>> qr;
        int n = qr.NaiveIterate(h.H(), 1000, SMALL_VALUE);
        XTLog(std::cout) << "Upper Hessenberg 迭代次数: " << n << std::endl;

        auto eigen_values = qr.EigenValues();
        std::sort(eigen_values.begin(), eigen_values.end(), [](double a, double b){ return a > b; });
        XTLog(std::cout) << "T = " << qr.T() << std::endl;

        MatrixView<double, 5, 1> lambdas(eigen_values.data());
        XTLog(std::cout) << "lambdas = " << lambdas << std::endl;
    }

    {
        UpperHessenberg h(A);

        EigenShiftQR<Matrix<double, 5, 5>> qr;
        int n = qr.Iterate(h.H(), 1000, SMALL_VALUE);
        XTLog(std::cout) << "Upper Hessenberg 偏移迭代次数: " << n << std::endl;

        auto eigen_values = qr.EigenValues();
        std::sort(eigen_values.begin(), eigen_values.end(), [](double a, double b){ return a > b; });
        XTLog(std::cout) << "T = " << qr.T() << std::endl;

        MatrixView<double, 5, 1> lambdas(eigen_values.data());
        XTLog(std::cout) << "lambdas = " << lambdas << std::endl;
    }
}






