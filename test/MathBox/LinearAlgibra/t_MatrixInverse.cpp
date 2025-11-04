#include <iostream>

#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

#include <gtest/gtest.h>

#include <memory>
#include <vector>
#include <cmath>

using namespace xiaotu::math;

TEST(LinearAlgibra, LU)
{
    double rad = M_PI * 0.15;
    double c = std::cos(rad);
    double s = std::sin(rad);

    Matrix<double, 3, 3> A = {
        c, -s, 0,
        s,  c, 0,
        0,  0, 1
    };

    auto A_inv = A.InverseMat();
    auto eye = A * A_inv;
    EXPECT_TRUE(eye.IsIdentity());
    
    XTLog(std::cout) << "A = " << A << std::endl;
    XTLog(std::cout) << "Inverse = " << A_inv << std::endl;
    XTLog(std::cout) << "eye = " << eye << std::endl;
}



