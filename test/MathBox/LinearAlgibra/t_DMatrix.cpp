#include <iostream>

#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

#include <gtest/gtest.h>

#include <memory>
#include <vector>
#include <cmath>

#include <Eigen/Eigen>

using namespace xiaotu;


TEST(LinearAlgibra, DMatrix)
{
    {
        DMatrix<double> A(3, 3);
        A = {
            9, -3, 1,
            1,  1, 1,
            4,  2, 1
        };

        VMatrix<double, 3, 3> v = A;
        XTLog(std::cout) << v << std::endl;

        AMatrix<double, 3, 3> a = A;
        XTLog(std::cout) << a << std::endl;
    }

    {
        VMatrix<double, 3, 3> A = {
            9, -3, 1,
            1,  1, 1,
            4,  2, 1
        };

        DMatrix<double> d(3, 3);
        d = A;
        XTLog(std::cout) << d << std::endl;
    }

    {
        DMatrix<double> A;
        XTLog(std::cout) << A << std::endl;

        A.Resize(3, 3);
        XTLog(std::cout) << A << std::endl;
    }
}

TEST(LinearAlgibra, GaussJordanEliminate)
{
    DMatrix<double> A(3, 3);
    A = {
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
    DMatrix<double> A(3, 3);
    A = {
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

    DMatrix<double> Ainv(3, 3);
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

    DMatrix<double> D(3, 3);
    D = {
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

    DMatrix<double> D(3, 3);
    D = {
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

    DMatrix<double> ATDA_inv(3, 3);
    ldlt.Inverse(ATDA_inv.View());
    XTLog(std::cout) << "ATDA_inv = " << ATDA_inv << std::endl;
    XTLog(std::cout) << "ha = " << ATDA_inv * ATDA << std::endl;
}

TEST(LinearAlgibra, Operations)
{
    DMatrix<double> A(3, 3);
    A = {
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

TEST(LinearAlgibra, GaussRowEliminate)
{
    {
        DMatrix<double> A(3, 3);
        A = {
            1,  2, 0,
            2,  4, 0,
            0,  0, 1
        };
        auto max_indep_set = GaussRowEliminate(A);
        EXPECT_EQ(2, max_indep_set.size());
        XTLog(std::cout) << A << std::endl;
    }
    {
        DMatrix<double> A(3, 3);
        A = {
            1,  2, 0,
            2,  4, 0,
            0,  0, 1
        };
        auto max_indep_set = FindMaximalIndepSet(A);
        EXPECT_EQ(2, max_indep_set.size());        
        XTLog(std::cout) << A << std::endl;
        XTLog(std::cout) << max_indep_set[0] << std::endl;
        XTLog(std::cout) << max_indep_set[1] << std::endl;
    }
    {
        DMatrix<double> A(3, 3);
        A = {
            1,  2, 0,
            2,  4, 0,
            0,  0, 1
        };
        EXPECT_EQ(2, Rank(A));
        XTLog(std::cout) << A << std::endl;
    }
    {
        DMatrix<double> A(4, 4);
        A = {
            1, 2, 3, 4,
            1, 2, 0, -5,
            2, 4, -3, -19,
            3, 6, -3, -24
        };
        auto max_indep_set = GaussRowEliminate(A);
        EXPECT_EQ(2, max_indep_set.size());
        XTLog(std::cout) << A << std::endl;

        XTLog(std::cout) << SolveSpace(max_indep_set, A) << std::endl;
    }

}

TEST(LinearAlgibra, QR_Householder)
{
    Matrix<double, 3, 4> A = {
        9, -3, 1, 4,
        1,  1, 1, 4,
        4,  2, 1, 4
    };

    QR_Householder<Matrix<double, 3, 4>> qr;
    qr.Decompose(A);
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

