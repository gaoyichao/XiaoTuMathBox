#include <iostream>

#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

#include <gtest/gtest.h>

#include <memory>
#include <vector>
#include <cmath>

using namespace xiaotu::math;


TEST(MatrixComma, MatrixView)
{
    int len = 6;
    std::vector<double> buffer(len);

    MatrixView<double, 2, 3> m(buffer.data());

    m << 1.9, 2.1, 3.3,
         4.5, 5.2, 6.2;

    XTLog(std::cout) << m << std::endl;
}

TEST(MatrixComma, Matrix)
{
    Matrix<double, 3, 3> A = {
        9, -3, 1,
        1,  1, 1,
        4,  2, 1
    };
    Matrix<double, 3, 1> b = { 20, 0, 10 };

    Matrix<double, 6, 4> Ab;
    Ab << A, b,
          b, A;

    XTLog(std::cout) << Ab << std::endl;
}


