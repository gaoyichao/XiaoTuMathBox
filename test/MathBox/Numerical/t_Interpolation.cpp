#include <iostream>

#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Numerical/Numerical.hpp>

#include <cmath>
#include <gtest/gtest.h>



TEST(Interpolation, LagrangePoly)
{
    {
        std::vector<double> x_nodes = {1.0, 2.0, 4.0};
        std::vector<double> y_nodes = {1.0, 4.0, 16.0};

        for (size_t i = 0; i < x_nodes.size(); i++) {
            double y = xiaotu::LagrangeInterpolation(x_nodes, y_nodes, x_nodes[i]);
            EXPECT_DOUBLE_EQ(y_nodes[i], y);
        }
    }

    {
        std::vector<double> x_nodes = {1.0, 2.0, 4.0};
        std::vector<double> y_nodes = {1.0, 4.0, 10.0};
        auto poly = xiaotu::LagrangePolynomial<xiaotu::Polynomial<double>>(x_nodes, y_nodes);
        XTLog(std::cout) << poly << std::endl;

        for (size_t i = 0; i < x_nodes.size(); i++) {
            double y = poly(x_nodes[i]);
            auto diff = std::abs(y_nodes[i] - y);
            XTLog(std::cout) << "[" << i << "]: " << diff << std::endl;
            EXPECT_TRUE(diff < SMALL_VALUE);
        }
    }

}

TEST(Interpolation, BarycentricLagrange)
{

    {
        std::vector<double> x_nodes = {1.0, 2.0, 3.0,  4.0};
        std::vector<double> y_nodes = {1.0, 4.0, 9.0, 16.0};
        auto bl = xiaotu::BarycentricLagrange(x_nodes, y_nodes);

        for (size_t i = 0; i < x_nodes.size(); i++) {
            double y = bl(x_nodes[i]);
            EXPECT_DOUBLE_EQ(y_nodes[i], y);
        }

        for (double x = 0.0; x < 5.0; x += 0.1) {
            double diff = x * x - bl(x);
            EXPECT_TRUE(std::abs(diff) < SMALL_VALUE);
        }

        EXPECT_DOUBLE_EQ(1.5 * 1.5, bl(1.5));
        EXPECT_DOUBLE_EQ(2.5 * 2.5, bl(2.5));
        EXPECT_DOUBLE_EQ(3.5 * 3.5, bl(3.5));
        EXPECT_DOUBLE_EQ(4.5 * 4.5, bl(4.5));
    }
}

TEST(Interpolation, NewtonDividedDifference)
{

    {
        std::vector<double> x_nodes = {1.0, 2.0, 3.0,  4.0};
        std::vector<double> y_nodes = {1.0, 4.0, 9.0, 16.0};
        auto ndd = xiaotu::NewtonDividedDifference(x_nodes, y_nodes);
        XTLog(std::cout) << ndd << std::endl;

        for (size_t i = 0; i < x_nodes.size(); i++) {
            double y = ndd(x_nodes[i]);
            EXPECT_DOUBLE_EQ(y_nodes[i], y);
        }

        for (double x = 0.0; x < 5.0; x += 0.1) {
            double diff = x * x - ndd(x);
            EXPECT_TRUE(std::abs(diff) < SMALL_VALUE);
        }

        EXPECT_DOUBLE_EQ(1.5 * 1.5, ndd(1.5));
        EXPECT_DOUBLE_EQ(2.5 * 2.5, ndd(2.5));
        EXPECT_DOUBLE_EQ(3.5 * 3.5, ndd(3.5));
        EXPECT_DOUBLE_EQ(4.5 * 4.5, ndd(4.5));
    }

    {
        std::vector<double> x_nodes = {1.0, 2.0, 3.0};
        std::vector<double> y_nodes = {1.0, 4.0, 9.0}; 
        
        auto ndd = xiaotu::NewtonDividedDifference(x_nodes, y_nodes);
        XTLog(std::cout) << ndd << std::endl;
        XTLog(std::cout) << "x=1.5: " << ndd(1.5) << " Expected: " << (1.5*1.5) << std::endl;
        
        // 动态新增采样点
        ndd.AddPoint(4.0, 16.0);
        XTLog(std::cout) << "AddPoint(4.0, 16.0):" << ndd << std::endl;
        XTLog(std::cout) << "x=2.5: " << ndd(2.5) << " Expected: " << (2.5*2.5) << std::endl;
        
        ndd.AddPoint(5.0, 25.0);
        XTLog(std::cout) << "AddPoint(5.0, 25.0):" << ndd << std::endl;
        XTLog(std::cout) << "x=3.5: " << ndd(3.5) << " Expected: " << (3.5*3.5) << std::endl;
    }
}

