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



