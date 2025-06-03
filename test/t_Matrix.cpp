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

TEST(LinearAlgibra, GramSchmidt)
{
    Matrix<double, 3, 3> A = {
        9, -3, 1,
        1,  1, 1,
        4,  2, 1
    };

    Matrix<double, 3, 3> ortho;
    GramSchmidt(A, ortho);
    XTLog(std::cout) << "A = " << A << std::endl;
    XTLog(std::cout) << "ortho = " << ortho << std::endl;

    auto ortho_0 = ortho.SubMatrix(0, 0, 3, 1);
    auto ortho_1 = ortho.SubMatrix(0, 1, 3, 1);
    auto ortho_2 = ortho.SubMatrix(0, 2, 3, 1);

    EXPECT_TRUE(std::abs(ortho_0.Dot(ortho_1)) < 1e-9);
    EXPECT_TRUE(std::abs(ortho_0.Dot(ortho_2)) < 1e-9);
    EXPECT_TRUE(std::abs(ortho_1.Dot(ortho_2)) < 1e-9);

    auto eye = ortho.Transpose() * ortho;
    XTLog(std::cout) << "eye = " << eye << std::endl;
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

    auto c = Ainv * b;
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
    auto ATDA = AT * D * A;
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
    auto ATDA = AT * D * A;
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
    EXPECT_TRUE(std::abs(d(0) - a(0)) < SMALL_VALUE);
    EXPECT_TRUE(std::abs(d(1) - a(1)) < SMALL_VALUE);
    EXPECT_TRUE(std::abs(d(2) - a(2)) < SMALL_VALUE);

    auto e = 2.0 * b;
    EXPECT_TRUE(std::abs(e(0) - c(0)) < SMALL_VALUE);
    EXPECT_TRUE(std::abs(e(1) - c(1)) < SMALL_VALUE);
    EXPECT_TRUE(std::abs(e(2) - c(2)) < SMALL_VALUE);

    e -= b;
    EXPECT_TRUE(std::abs(e(0) - a(0)) < SMALL_VALUE);
    EXPECT_TRUE(std::abs(e(1) - a(1)) < SMALL_VALUE);
    EXPECT_TRUE(std::abs(e(2) - a(2)) < SMALL_VALUE);

    e += b;
    EXPECT_TRUE(std::abs(e(0) - c(0)) < SMALL_VALUE);
    EXPECT_TRUE(std::abs(e(1) - c(1)) < SMALL_VALUE);
    EXPECT_TRUE(std::abs(e(2) - c(2)) < SMALL_VALUE);

    e *= 0.5;
    EXPECT_TRUE(std::abs(e(0) - a(0)) < SMALL_VALUE);
    EXPECT_TRUE(std::abs(e(1) - a(1)) < SMALL_VALUE);
    EXPECT_TRUE(std::abs(e(2) - a(2)) < SMALL_VALUE);

    EXPECT_TRUE(Saxpy(-2, a, e));
    EXPECT_TRUE(std::abs(e(0) + a(0)) < SMALL_VALUE);
    EXPECT_TRUE(std::abs(e(1) + a(1)) < SMALL_VALUE);
    EXPECT_TRUE(std::abs(e(2) + a(2)) < SMALL_VALUE);


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

TEST(LinearAlgibra, Gaxpy)
{
    Matrix<double, 3, 2> A = {
        1.0, 2.0,
        3.0, 4.0,
        5.0, 6.0
    };
    Vector<double, 2> x = { 7.0, 8.0 };
    Vector<double, 3> y = { 9.0, 10.0, 11.0 };

    EXPECT_TRUE(Gaxpy(A, x, y));
    EXPECT_DOUBLE_EQ(32, y(0));
    EXPECT_DOUBLE_EQ(63, y(1));
    EXPECT_DOUBLE_EQ(94, y(2));
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
        EXPECT_DOUBLE_EQ(4, A.PNorm(2));
        EXPECT_DOUBLE_EQ(3, A.InftyNorm());
        EXPECT_DOUBLE_EQ(16, A.SquaredNorm());
        EXPECT_DOUBLE_EQ(1.0, A.Normalize().Norm());

        XTLog(std::cout) << A << std::endl;
    }
}


TEST(LinearAlgibra, PowerIterate)
{
    Matrix<double, 3, 3> A = {
        2, 1, 1,
        1, 2, 1,
        1, 1, 2
    };
    Vector<double, 3> v = {
        0.333333, 0.3333333, 0.3333333
    };

    double lambda = PowerIterate(A, v, 1e-20);

    XTLog(std::cout) << lambda << std::endl;
    XTLog(std::cout) << v << std::endl;
}


TEST(LinearAlgibra, InversePowerIterate)
{
    Matrix<double, 3, 3> A = {
        2, 1, 1,
        1, 2, 1,
        1, 1, 2
    };
    Vector<double, 3> v = {
        0.333333, 0.3333333, 0.3333333
    };

    double lambda = InversePowerIterate(A, v, 1e-20);

    XTLog(std::cout) << lambda << std::endl;
    XTLog(std::cout) << v << std::endl;
}


TEST(LinearAlgibra, OffInvPowerIterate)
{
    Matrix<double, 3, 3> A = {
        2, 1, 1,
        1, 2, 1,
        1, 1, 2
    };
    Vector<double, 3> v = {
        0.333333, 0.3333333, 0.3333333
    };

    double lambda = OffInvPowerIterate(A, v, 4);

    XTLog(std::cout) << lambda << std::endl;
    XTLog(std::cout) << v << std::endl;
    XTLog(std::cout) << v.Rows() << std::endl;
    XTLog(std::cout) << v.Cols() << std::endl;
}
