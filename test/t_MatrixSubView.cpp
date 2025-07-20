#include <iostream>

#include <XiaoTuMathBox/Utils.hpp>
#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

#include <gtest/gtest.h>

#include <memory>
#include <vector>
#include <cmath>

using namespace xiaotu::math;

TEST(MatrixSubView, Assign)
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

TEST(MatrixSubView, SubSubView)
{
    auto all = DMatrix<double>::Zero(10, 10);
    for (int i = 0; i < 10; i++)
        all(i, i) = 1.0 * i;

    {
        const auto sub_5_10 = all.SubMatrix(5, 5, 5, 5);
        XTLog(std::cout) << "------------------------" << std::endl;
        XTLog(std::cout) << "sub_5_10(0, 0):" << sub_5_10(0, 0) << std::endl;

        XTLog(std::cout) << "----" << std::endl;
        auto sub_0_3 = sub_5_10.SubMatrix(0, 0, 3, 3);
        XTLog(std::cout) << "sub_0_3(0, 0):" << sub_0_3(0, 0) << std::endl;

        XTLog(std::cout) << "----" << std::endl;
        auto sub_0_1 = sub_0_3.SubMatrix(0, 0, 1, 1);
        XTLog(std::cout) << "sub_0_1(0, 0):" << sub_0_1(0, 0) << std::endl;

    }
}


TEST(MatrixSubView, SubSub)
{
    Matrix<double, 5, 5> A = {
        0, 1, 2, 3, 4,
        5, 6, 7, 8, 9,
        10, 11, 12, 13, 14,
        15, 16, 17, 18, 19,
        20, 21, 22, 23, 24,
    };

    {
        auto sub = A.SubMatrix(0, 0, 1, 3);
        EXPECT_DOUBLE_EQ(0, sub(0));
        EXPECT_DOUBLE_EQ(1, sub(1));
        EXPECT_DOUBLE_EQ(2, sub(2));
        XTLog(std::cout) << sub << std::endl;
    }

    {
        auto sub = A.SubMatrix(1, 1, 3, 3);
        EXPECT_DOUBLE_EQ(6, sub(0, 0));
        EXPECT_DOUBLE_EQ(7, sub(0, 1));
        EXPECT_DOUBLE_EQ(8, sub(0, 2));

        EXPECT_DOUBLE_EQ(11, sub(1, 0));
        EXPECT_DOUBLE_EQ(12, sub(1, 1));
        EXPECT_DOUBLE_EQ(13, sub(1, 2));

        EXPECT_DOUBLE_EQ(16, sub(2, 0));
        EXPECT_DOUBLE_EQ(17, sub(2, 1));
        EXPECT_DOUBLE_EQ(18, sub(2, 2));

        auto ssub = sub.SubMatrix(1, 1, 2, 2);
        EXPECT_DOUBLE_EQ(12, ssub(0, 0));
        EXPECT_DOUBLE_EQ(13, ssub(0, 1));
        EXPECT_DOUBLE_EQ(17, ssub(1, 0));
        EXPECT_DOUBLE_EQ(18, ssub(1, 1));

    }
}

TEST(MatrixSubView, ColSwap)
{
    Matrix<double, 5, 5> A = {
        0, 1, 2, 3, 4,
        5, 6, 7, 8, 9,
        10, 11, 12, 13, 14,
        15, 16, 17, 18, 19,
        20, 21, 22, 23, 24,
    };

    auto sub = A.SubMatrix(0, 0, 1, 3);
    EXPECT_DOUBLE_EQ(0, sub(0));
    EXPECT_DOUBLE_EQ(1, sub(1));
    EXPECT_DOUBLE_EQ(2, sub(2));
    XTLog(std::cout) << sub << std::endl;

    A.ColSwap(0, 2);
    EXPECT_DOUBLE_EQ(2, sub(0));
    EXPECT_DOUBLE_EQ(1, sub(1));
    EXPECT_DOUBLE_EQ(0, sub(2));
    XTLog(std::cout) << sub << std::endl;
}




