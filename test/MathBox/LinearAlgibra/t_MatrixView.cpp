#include <iostream>

#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

#include <gtest/gtest.h>

#include <memory>
#include <vector>
#include <cmath>

using namespace xiaotu;

TEST(MatrixView, Assign)
{
    int len = 6;
    std::vector<double> buffer(len);

    {
        MatrixView<double, 2, 3> m(buffer.data());
        m = { 1.0, 3.0, 5.0,
            2.0, 4.0, 6.0 };
        for (int i = 0; i < len; i++)
            EXPECT_DOUBLE_EQ(1.0 + i, buffer[i]);

        m = { { 10.0, 20.0 },
            { 30.0 } };
        EXPECT_DOUBLE_EQ(10, buffer[0]);
        EXPECT_DOUBLE_EQ(30, buffer[1]);
        EXPECT_DOUBLE_EQ(20, buffer[2]);
    }

    {
        MatrixView<double, 2, 3, eRowMajor> m(buffer.data());
        m = { 1.0, 2.0, 3.0,
              4.0, 5.0, 6.0 };
        for (int i = 0; i < len; i++)
            EXPECT_DOUBLE_EQ(1.0 + i, buffer[i]);

        m = { { 10.0, 20.0 },
              { 30.0 } };
        EXPECT_DOUBLE_EQ(10, buffer[0]);
        EXPECT_DOUBLE_EQ(20, buffer[1]);
        EXPECT_DOUBLE_EQ(3, buffer[2]);
    }

    {
        DMatrixView<double, eRowMajor> h(buffer.data(), 2, 3);
        MatrixView<double, 2, 3, eRowMajor> m(h);

        XTLog(std::cout) << "m = " << m << std::endl;
    }

    {
        const DMatrixView<double, eRowMajor> h(buffer.data(), 2, 3);
        MatrixView<const double, 2, 3, eRowMajor> m(h);

        XTLog(std::cout) << "m = " << m << std::endl;
    }

}

TEST(MatrixView, Swap)
{
    int len = 16;
    std::vector<double> buffer(len);

    MatrixView<double, 4, 4> m(buffer.data());
    for (int ridx = 0; ridx < 4; ridx++) {
        for (int cidx = 0; cidx < 4; cidx++) {
            m(ridx, cidx) = ridx + 0.1 * cidx;
        }
    }

    m.RowSwap(0, 2);
    for (int cidx = 0; cidx < 4; cidx++) {
        EXPECT_DOUBLE_EQ(0 + 0.1 * cidx, m(2, cidx));
        EXPECT_DOUBLE_EQ(2 + 0.1 * cidx, m(0, cidx));
    }

    m.RowSwap(0, 2);
    m.ColSwap(0, 2);
    for (int ridx = 0; ridx < 4; ridx++) {
        EXPECT_DOUBLE_EQ(ridx + 0.1 * 0, m(ridx, 2));
        EXPECT_DOUBLE_EQ(ridx + 0.1 * 2, m(ridx, 0));
    }
}

TEST(MatrixView, GaussJordanEliminate)
{
    std::vector<double> _A_(9);
    MatrixView<double, 3, 3> A(_A_.data());
    A = {
        9, -3, 1,
        1,  1, 1,
        4,  2, 1
    };

    std::vector<double> _b_(3);
    MatrixView<double, 3, 1> b(_b_.data());
    b = { 20, 0, 10 };

    GaussJordanEliminate(A, &b);
    XTLog(std::cout) << "A = " << A << std::endl;
    XTLog(std::cout) << "b = " << b << std::endl;

    EXPECT_TRUE(std::abs(b(0) -   3)  < 1e-9);
    EXPECT_TRUE(std::abs(b(1) -   1)  < 1e-9);
    EXPECT_TRUE(std::abs(b(2) - (-4)) < 1e-9);
}

TEST(MatrixView, LU)
{
    std::vector<double> _A_(9);
    MatrixView<double, 3, 3> A(_A_.data());
    A = {
        9, -3, 1,
        1,  1, 1,
        4,  2, 1
    };

    LU lu(A);
    XTLog(std::cout) << "lu = " << lu() << std::endl;

    std::vector<double> _b_(3);
    std::vector<double> _x_(3);
    MatrixView<double, 3, 1> b(_b_.data());
    MatrixView<double, 3, 1> x(_x_.data());
    b = { 20, 0, 10 };

    lu.Solve(b, x);
    XTLog(std::cout) << "x = " << x << std::endl;
    XTLog(std::cout) << "b = " << b << std::endl;

    std::vector<double> _Ainv_(9);
    MatrixView<double, 3, 3> Ainv(_Ainv_.data());
    lu.Inverse(Ainv);
    XTLog(std::cout) << "A^{-1} = " << Ainv << std::endl;


    std::vector<double> _c_(3);
    MatrixView<double, 3, 1> c(_c_.data());
    bool re = Multiply(Ainv, b, c);
    EXPECT_TRUE(re);
    XTLog(std::cout) << "c = " << c << std::endl;
}

TEST(MatrixView, Cholesky)
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

TEST(MatrixView, SubView)
{
    int len = 6;
    std::vector<double> buffer(len);
    MatrixView<double, 2, 3> m(buffer.data());

    {
        m = { 1.0, 3.0, 5.0,
            2.0, 4.0, 6.0 };
        for (int i = 0; i < len; i++)
            EXPECT_DOUBLE_EQ(1.0 + i, buffer[i]);

        m = { { 10.0, 20.0 },
            { 30.0 } };
        EXPECT_DOUBLE_EQ(10, buffer[0]);
        EXPECT_DOUBLE_EQ(30, buffer[1]);
        EXPECT_DOUBLE_EQ(20, buffer[2]);

        auto sub = m.SubMatrix(0, 0, 1, 3);
        EXPECT_DOUBLE_EQ(10, sub(0));
        EXPECT_DOUBLE_EQ(20, sub(1));
        EXPECT_DOUBLE_EQ(5, sub(2));
    }

    {
        MatrixView<double, 2, 3> const & mm = m;
        auto sub = mm.SubMatrix(0, 0, 1, 3);
        XTLog(std::cout) << sub << std::endl;
        XTLog(std::cout) << sub.Dot(sub) << std::endl;

        Matrix<double, 2, 3> A = {
            9, -3, 1,
            1,  1, 1
        };

        auto sub_A = A.SubMatrix(1, 0, 1, 3);
        sub_A = sub;
        XTLog(std::cout) << A << std::endl;

        auto sub_mm = mm.SubMatrix(1, 0, 1, 3);
        sub_A << sub_mm;
        XTLog(std::cout) << A << std::endl;
    }
}

TEST(MatrixView, Col_Row)
{
    std::vector<double> buffer(64);
    MatrixView<double, 8, 8> m(buffer.data());

    m << 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
         21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0,
         31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0,
         41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0,
         51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0,
         61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0,
         71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0,
         81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0;

    {
        MatrixColView cols(m, {1, 3, 5, 7, 1});
        MatrixRowView rows(cols, {1, 3, 5, 7});
        XTLog(std::cout) << "cols-rows:" << std::endl;
        XTLog(std::cout) << rows << std::endl;
    }

    {
        MatrixConstColView const_cols(m, {1, 3, 5, 7});
        MatrixConstRowView const_rows(const_cols, {1, 3, 5, 7});
        XTLog(std::cout) << "cols-rows:" << std::endl;
        XTLog(std::cout) << const_rows << std::endl;
    }

}


