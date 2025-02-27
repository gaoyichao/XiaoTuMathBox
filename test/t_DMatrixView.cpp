#include <iostream>

#include <XiaoTuMathBox/Utils.hpp>
#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

#include <gtest/gtest.h>

#include <memory>
#include <vector>
#include <cmath>

using namespace xiaotu::math;

TEST(DMatrixView, Assign)
{
    int len = 6;
    std::vector<double> buffer(len);

    {
        DMatrixView<double> m(buffer.data(), 2, 3);
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
        DMatrixView<double, eRowMajor> m(buffer.data(), 2, 3);
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
        MatrixView<double, 2, 3, eRowMajor> m(buffer.data());
        DMatrixView<double, eRowMajor> h = m;

        XTLog(std::cout) << "h = " << h << std::endl;
    }

    {
        const MatrixView<double, 2, 3, eRowMajor> m(buffer.data());
        DMatrixView<const double, eRowMajor> h = m;

        XTLog(std::cout) << "h = " << h << std::endl;
    }
}

TEST(DMatrixView, Swap)
{
    int len = 16;
    std::vector<double> buffer(len);

    DMatrixView<double> m(buffer.data(), 4, 4);
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

TEST(DMatrixView, GaussJordanEliminate)
{
    std::vector<double> _A_(9);
    DMatrixView<double> A(_A_.data(), 3, 3);
    A = {
        9, -3, 1,
        1,  1, 1,
        4,  2, 1
    };

    std::vector<double> _b_(3);
    DMatrixView<double> b(_b_.data(), 3, 1);
    b = { 20, 0, 10 };

    GaussJordanEliminate(A, &b);
    XTLog(std::cout) << "A = " << A << std::endl;
    XTLog(std::cout) << "b = " << b << std::endl;

    EXPECT_TRUE(std::abs(b(0) -   3)  < 1e-9);
    EXPECT_TRUE(std::abs(b(1) -   1)  < 1e-9);
    EXPECT_TRUE(std::abs(b(2) - (-4)) < 1e-9);
}

TEST(DMatrixView, LU)
{
    std::vector<double> _A_(9);
    DMatrixView<double> A(_A_.data(), 3, 3);
    A = {
        9, -3, 1,
        1,  1, 1,
        4,  2, 1
    };

    LU lu(A);
    XTLog(std::cout) << "lu = " << lu() << std::endl;

    std::vector<double> _b_(3);
    std::vector<double> _x_(3);
    DMatrixView<double> b(_b_.data(), 3, 1);
    DMatrixView<double> x(_x_.data(), 3, 1);
    b = { 20, 0, 10 };

    lu.Solve(b, x);
    XTLog(std::cout) << "x = " << x << std::endl;
    XTLog(std::cout) << "b = " << b << std::endl;

    std::vector<double> _Ainv_(9);
    DMatrixView<double> Ainv(_Ainv_.data(), 3, 3);
    lu.Inverse(Ainv);
    XTLog(std::cout) << "A^{-1} = " << Ainv << std::endl;


    std::vector<double> _c_(3);
    DMatrixView<double> c(_c_.data(), 3, 1);
    bool re = Multiply(Ainv, b, c);
    EXPECT_TRUE(re);
    XTLog(std::cout) << "c = " << c << std::endl;
}









