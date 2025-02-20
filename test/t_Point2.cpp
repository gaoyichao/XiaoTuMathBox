#include <iostream>

#include <XiaoTuMathBox/Utils.hpp>
#include <XiaoTuMathBox/Geometry/Point2.hpp>

#include <gtest/gtest.h>

#include <memory>
#include <vector>
#include <cmath>

using namespace xiaotu::math;



TEST(Geometry, Point2)
{
    {
        Point2<double> p = { 9, -3};
        XTLog(std::cout) << "p = " << p << std::endl;
        XTLog(std::cout) << "sizeof(p): " << sizeof(p) << std::endl;
    }

}

