#include <iostream>

#include <XiaoTuMathBox/Utils.hpp>
#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

#include <gtest/gtest.h>

#include <memory>
#include <vector>
#include <cmath>

using namespace xiaotu::math;


TEST(Eigen, EigenNaiveQR)
{
    Matrix<double, 3, 3> A = {
        2, 1, 1,
        1, 2, 1,
        1, 1, 2
    };

    EigenNaiveQR qr(A, 20);
    auto eigen_values = qr.EigenValues();
    std::sort(eigen_values.begin(), eigen_values.end(), [](double a, double b){ return a > b; });
    EXPECT_TRUE(std::abs(eigen_values[0] - 4.0) < 1e-9);
    EXPECT_TRUE(std::abs(eigen_values[1] - 1.0) < 1e-9);
    EXPECT_TRUE(std::abs(eigen_values[2] - 1.0) < 1e-9);

    XTLog(std::cout) << "T = " << qr.T() << std::endl;
}

TEST(Eigen, EigenShiftQR)
{
    Matrix<double, 5, 5> A = {
        3, 1, 1, 1, 1,
        1, 2, 1, 1, 1,
        1, 1, 2, 1, 1,
        1, 1, 1, 2, 1,
        1, 1, 1, 1, 3,
    };

    {
        EigenShiftQR<Matrix<double, 5, 5>> qr;
        int n = qr.NaiveIterate(A, 1000, SMALL_VALUE);
        XTLog(std::cout) << "朴素迭代次数: " << n << std::endl;

        auto eigen_values = qr.EigenValues();
        std::sort(eigen_values.begin(), eigen_values.end(), [](double a, double b){ return a > b; });
        XTLog(std::cout) << "T = " << qr.T() << std::endl;

        MatrixView<double, 5, 1> lambdas(eigen_values.data());
        XTLog(std::cout) << "lambdas = " << lambdas << std::endl;
    }

    {
        UpperHessenberg h(A);

        EigenShiftQR<Matrix<double, 5, 5>> qr;
        int n = qr.NaiveIterate(h.H(), 1000, SMALL_VALUE);
        XTLog(std::cout) << "Upper Hessenberg 迭代次数: " << n << std::endl;

        auto eigen_values = qr.EigenValues();
        std::sort(eigen_values.begin(), eigen_values.end(), [](double a, double b){ return a > b; });
        XTLog(std::cout) << "T = " << qr.T() << std::endl;

        MatrixView<double, 5, 1> lambdas(eigen_values.data());
        XTLog(std::cout) << "lambdas = " << lambdas << std::endl;
    }

    {
        UpperHessenberg h(A);

        EigenShiftQR<Matrix<double, 5, 5>> qr;
        int n = qr.Iterate(h.H(), 1000, SMALL_VALUE);
        XTLog(std::cout) << "Upper Hessenberg 偏移迭代次数: " << n << std::endl;

        auto eigen_values = qr.EigenValues();
        std::sort(eigen_values.begin(), eigen_values.end(), [](double a, double b){ return a > b; });
        XTLog(std::cout) << "T = " << qr.T() << std::endl;

        MatrixView<double, 5, 1> lambdas(eigen_values.data());
        XTLog(std::cout) << "lambdas = " << lambdas << std::endl;
    }
}


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




