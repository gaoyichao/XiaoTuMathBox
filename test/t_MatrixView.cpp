#include <iostream>

#include <XiaoTuMathBox/Utils.hpp>
#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

#include <gtest/gtest.h>

#include <memory>
#include <vector>
#include <cmath>

using namespace xiaotu::math;

TEST(MatrixView, Assign)
{
    int len = 6;
    std::vector<double> buffer(len);

    MatrixView<double, 2, 3> m(buffer.data());
    m = { 1.0, 2.0, 3.0,
           4.0, 5.0, 6.0 };
    XTLog(std::cout) << "m = " << m << std::endl;

    m = { { 10.0, 20.0 },
           { 30.0 } };
    XTLog(std::cout) << "m = " << m << std::endl;

    std::cout << buffer[0] << std::endl;
    std::cout << buffer[1] << std::endl;
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
        EXPECT_EQ(0 + 0.1 * cidx, m(2, cidx));
        EXPECT_EQ(2 + 0.1 * cidx, m(0, cidx));
    }

    m.RowSwap(0, 2);
    m.ColSwap(0, 2);
    for (int ridx = 0; ridx < 4; ridx++) {
        EXPECT_EQ(ridx + 0.1 * 0, m(ridx, 2));
        EXPECT_EQ(ridx + 0.1 * 2, m(ridx, 0));
    }
}

TEST(LinearAlgibra, GaussJordanEliminate)
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

