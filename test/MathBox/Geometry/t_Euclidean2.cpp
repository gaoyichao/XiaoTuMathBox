#include <iostream>
#include <XiaoTuMathBox/Geometry/Geometry.hpp>

#include <gtest/gtest.h>

TEST(Euclidean2, PointAndVector)
{
    using namespace xiaotu::math;

    Point2<double> p0(3.14159, 1.41421);
    EXPECT_DOUBLE_EQ(3.14159, p0.x());
    EXPECT_DOUBLE_EQ(1.41421, p0.y());

    Point2<double> p1(0, 1.41421);
    Vector2<double> v(p1 - p0);
    EXPECT_DOUBLE_EQ(-3.14159, v.x());
    EXPECT_DOUBLE_EQ(0.0, v.y());
 
    v = p1 - p0;
    v.Normalize();
    EXPECT_DOUBLE_EQ(-1, v.x());
}


TEST(Euclidean2, RadianAngles)
{
    using namespace xiaotu::math;
    Vector2<double> v0(0, 1);
    Vector2<double> v1(0, 1);

    EXPECT_DOUBLE_EQ(1.0, Cos(v0, v1));
    EXPECT_DOUBLE_EQ(1.0, Cos(v1, v0));
    EXPECT_DOUBLE_EQ(0.0, Sin(v0, v1));

    v1 = { 0.0, -1.0 };
    EXPECT_DOUBLE_EQ(-1.0, Cos(v0, v1));
    EXPECT_DOUBLE_EQ(-1.0, Cos(v1, v0));

    v1 = { 1.0, 0.0 };
    EXPECT_DOUBLE_EQ(0.0, Cos(v0, v1));
    EXPECT_DOUBLE_EQ(-1.0, Sin(v0, v1));
    EXPECT_DOUBLE_EQ(1.0, Sin(v1, v0));

    EXPECT_DOUBLE_EQ(90.0, 180 * Radian(v1, v0) / M_PI);
    EXPECT_DOUBLE_EQ(-90.0, 180 * Radian(v0, v1) / M_PI);

    Point2<double> p0(0, 0);
    Point2<double> p1(0, 1);
    Point2<double> p2(1, 1);

    EXPECT_DOUBLE_EQ(90.0, 180 * Radian(p2, p1, p0) / M_PI);
    EXPECT_DOUBLE_EQ(-90.0, 180 * Radian(p0, p1, p2) / M_PI);
}




