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
        std::vector<double> x_nodes = {1.0, 2.0, 4.0};
        std::vector<double> y_nodes = {1.0, 4.0, 10.0};
        auto poly = xiaotu::BarycentricLagrange(x_nodes, y_nodes);

        for (size_t i = 0; i < x_nodes.size(); i++) {
            double y = poly(x_nodes[i] + 2 * SMALL_VALUE);
            // auto diff = std::abs(y_nodes[i] - y);
            XTLog(std::cout) << "[" << i << "]: " << y << std::endl;
            // EXPECT_TRUE(diff < SMALL_VALUE);
            // EXPECT_DOUBLE_EQ(y_nodes[i], y);
        }
    }
}



