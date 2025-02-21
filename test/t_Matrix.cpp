#include <iostream>

#include <XiaoTuMathBox/Utils.hpp>
#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

#include <gtest/gtest.h>

#include <memory>
#include <vector>
#include <cmath>

#include <Eigen/Eigen>

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

TEST(LinearAlgibra, Cholesky)
{
    std::vector<double> _A_(9);
    MatrixView<double, 3, 3> A(_A_.data());
    MatrixView<double, 3, 3, EAlignType::eRowMajor> AT(_A_.data());
    A = {
        9, -3, 1,
        1,  1, 1,
        4,  2, 1
    };

    Matrix<double, 3, 3> D = {
        3, 0, 0,
        0, 2, 0,
        0, 0, 1
    };
    Matrix<double, 3, 3> ATDA = AT * D * A;
    XTLog(std::cout) << "ATDA = " << ATDA << std::endl;

    Cholesky cholesky(ATDA.View());
    XTLog(std::cout) << "cholesky = " << cholesky() << std::endl;

    MatrixView<double, 3, 3, EAlignType::eRowMajor> lt(cholesky().StorBegin());
    XTLog(std::cout) << "lt = " << lt << std::endl;
    XTLog(std::cout) << "ha = " << cholesky() * lt << std::endl;

    Matrix<double, 3, 3> ATDA_inv;
    cholesky.Inverse(ATDA_inv.View());
    XTLog(std::cout) << "ATDA_inv = " << ATDA_inv << std::endl;
    XTLog(std::cout) << "ha = " << ATDA_inv * ATDA << std::endl;
}

TEST(LinearAlgibra, LDLT)
{
    std::vector<double> _A_(9);
    MatrixView<double, 3, 3> A(_A_.data());
    MatrixView<double, 3, 3, EAlignType::eRowMajor> AT(_A_.data());
    A = {
        9, -3, 1,
        1,  1, 1,
        4,  2, 1
    };

    Matrix<double, 3, 3> D = {
        3, 0, 0,
        0, 2, 0,
        0, 0, 1
    };
    Matrix<double, 3, 3> ATDA = AT * D * A;
    XTLog(std::cout) << "ATDA = " << ATDA << std::endl;

    LDLT ldlt(ATDA.View());
    auto l = ldlt.L();
    auto d = ldlt.D();
    auto lt = ldlt.LT();
    XTLog(std::cout) << "ldlt = " << l * d * lt << std::endl;

    Matrix<double, 3, 3> ATDA_inv;
    ldlt.Inverse(ATDA_inv.View());
    XTLog(std::cout) << "ATDA_inv = " << ATDA_inv << std::endl;
    XTLog(std::cout) << "ha = " << ATDA_inv * ATDA << std::endl;
}

TEST(LinearAlgibra, Operations)
{
    Vector<double, 3> a = { 20, 0.1, 10 };
    Vector<double, 3> b = { 20, 0.1, 10 };

    auto c = a + b;
    RowVector<double, 3> c_T = c.Transpose();
    XTLog(std::cout) << "c  = " << c << std::endl;
    XTLog(std::cout) << "cT = " << c_T << std::endl;
    XTLog(std::cout) << "cT *c = " << c_T * c << std::endl;
    EXPECT_DOUBLE_EQ(2000.04, c.Dot(c));

    auto d = c - b;
    EXPECT_TRUE(std::abs(d(0) - a(0)) < 1e-9);
    EXPECT_TRUE(std::abs(d(1) - a(1)) < 1e-9);
    EXPECT_TRUE(std::abs(d(2) - a(2)) < 1e-9);

    auto e = 2.0 * b;
    EXPECT_TRUE(std::abs(e(0) - c(0)) < 1e-9);
    EXPECT_TRUE(std::abs(e(1) - c(1)) < 1e-9);
    EXPECT_TRUE(std::abs(e(2) - c(2)) < 1e-9);

    e -= b;
    EXPECT_TRUE(std::abs(e(0) - a(0)) < 1e-9);
    EXPECT_TRUE(std::abs(e(1) - a(1)) < 1e-9);
    EXPECT_TRUE(std::abs(e(2) - a(2)) < 1e-9);

    e += b;
    EXPECT_TRUE(std::abs(e(0) - c(0)) < 1e-9);
    EXPECT_TRUE(std::abs(e(1) - c(1)) < 1e-9);
    EXPECT_TRUE(std::abs(e(2) - c(2)) < 1e-9);

    e *= 0.5;
    EXPECT_TRUE(std::abs(e(0) - a(0)) < 1e-9);
    EXPECT_TRUE(std::abs(e(1) - a(1)) < 1e-9);
    EXPECT_TRUE(std::abs(e(2) - a(2)) < 1e-9);

    Matrix<double, 3, 3> A = {
        1,  2, 3,
        1,  1, 1,
        4,  2, 1
    };
    auto B = Matrix<double, 3, 3>::Eye() * 2;
    auto C = A + A;
    auto D = A * B;
    for (int idx = 0; idx < A.NumDatas(); idx++)
        EXPECT_TRUE(std::abs(C(idx) - D(idx)) < 1e-9);

    A = A * B;
    for (int idx = 0; idx < A.NumDatas(); idx++)
        EXPECT_TRUE(std::abs(C(idx) - A(idx)) < 1e-9);

}

TEST(LinearAlgibra, VMatrix)
{
    {
        VMatrix<double, 3, 3> A = {
            9, -3, 1,
            1,  1, 1,
            4,  2, 1
        };
        EXPECT_EQ(sizeof(A), sizeof(std::vector<double>));
    }

    {
        VMatrix<double, 3, 3, eRowMajor> A = {
            9, -3, 1,
            1,  1, 1,
            4,  2, 1
        };
        EXPECT_EQ(sizeof(A), sizeof(std::vector<double>));
    }

    {
        AMatrix<double, 3, 3> A = {
            9, -3, 1,
            1,  1, 1,
            4,  2, 1
        };
        EXPECT_EQ(72, sizeof(double) * A.NumDatas());
    }

    {
        AMatrix<double, 3, 3> A = {
             1, 2, 3,
             0, 1, 0,
             0, 0, 1
        };
        EXPECT_DOUBLE_EQ(4, A.Norm());
        EXPECT_DOUBLE_EQ(16, A.SquaredNorm());
        EXPECT_DOUBLE_EQ(1.0, A.Normalize().Norm());

        XTLog(std::cout) << A << std::endl;
    }
}

