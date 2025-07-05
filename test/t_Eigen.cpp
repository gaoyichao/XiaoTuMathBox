#include <iostream>

#include <XiaoTuMathBox/Utils.hpp>
#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

#include <gtest/gtest.h>

#include <memory>
#include <vector>
#include <cmath>

using namespace xiaotu::math;

TEST(Eigen, ParitionQR)
{
    Matrix<double, 5, 5> A = {
        3, 1, 1, 1, 1,
        1, 2, 1, 1, 1,
        1, 1, 2, 1, 1,
        1, 1, 1, 2, 1,
        1, 1, 1, 1, 3,
    };

    {
        EigenPartitionQR<Matrix<double, 5, 5>> qr;
        int n = qr.Iterate(A, 25, SMALL_VALUE);
        XTLog(std::cout) << "迭代次数: " << n << std::endl;

        auto eigen_values = qr.EigenValues();
        std::sort(eigen_values.begin(), eigen_values.end(), [](double a, double b){ return a > b; });
        MatrixView<double, 5, 1> lambdas(eigen_values.data());
        XTLog(std::cout) << "lambdas = " << lambdas << std::endl;
    }
}


TEST(Eigen, ParitionQR_Complex_Eigen)
{
    Matrix<double, 5, 5> A = {
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1,
       -1, 0, 0, 0, 0,
    };

    {
        double first_off = -1.0;
        EigenPartitionQR<Matrix<double, 5, 5>> qr;
        int n = qr.Iterate(A, 1000, SMALL_VALUE, &first_off);
        XTLog(std::cout) << "迭代次数: " << n << std::endl;

        auto eigen_values = qr.EigenValues();
        std::sort(eigen_values.begin(), eigen_values.end(), [](double a, double b){ return a > b; });
        MatrixView<double, 5, 1> lambdas(eigen_values.data());
        XTLog(std::cout) << "lambdas = " << lambdas << std::endl;
    }
}


