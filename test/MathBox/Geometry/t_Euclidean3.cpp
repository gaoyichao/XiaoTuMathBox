#include <iostream>
#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Geometry/Geometry.hpp>

#include <gtest/gtest.h>

TEST(Euclidean3, PointAndVector)
{
    using namespace xiaotu::math;

    Point3<double> p0(3.14159, 1.41421, 2.71828);
    EXPECT_DOUBLE_EQ(3.14159, p0.x());
    EXPECT_DOUBLE_EQ(1.41421, p0.y());
    EXPECT_DOUBLE_EQ(2.71828, p0.z());

    Point3<double> p1(0, 1.41421, 0);
    Vector3<double> v = p1 - p0;
    EXPECT_DOUBLE_EQ(-3.14159, v.x());
    EXPECT_DOUBLE_EQ(0.0, v.y());
    EXPECT_DOUBLE_EQ(-2.71828, v.z());
 
    v = p1 - p0;
    v.Normalize();
    EXPECT_DOUBLE_EQ(1, v.Norm());
}

TEST(Euclidean3, Line)
{
    using namespace xiaotu::math;

    Point3<double> p0(1, 1, 0);
    Point3<double> p1(1, 2, 0);
    Line3<double> line0 = Join(p0, p1);
    XTLog(std::cout) << "line0: " << line0 << std::endl;

    {
        Line3<double> line1 = Join(p1, p0);
        EXPECT_TRUE(line0 == line1);
        EXPECT_TRUE(OnLine(p0, line1));
        EXPECT_TRUE(OnLine(p1, line1));
        EXPECT_TRUE(!line1.IsInfinity());
        XTLog(std::cout) << line1 << std::endl;
    }

    {
        Point3<double> p3(1, 9, 0);
        Point3<double> p4(1, 2, 0);

        Line3<double> line1 = Join(p3, p4);
        EXPECT_TRUE(line0 == line1);
        EXPECT_TRUE(OnLine(p3, line0));
        EXPECT_TRUE(OnLine(p4, line0));
        EXPECT_TRUE(!line1.IsInfinity());
        XTLog(std::cout) << line1 << std::endl;
    }


    {
        Point3<double> p3(1, 9, 1);
        Point3<double> p4(1, 2, 0);

        Line3<double> line1 = Join(p3, p4);
        EXPECT_TRUE(line0 != line1);
        EXPECT_TRUE(!OnLine(p3, line0));
        EXPECT_TRUE(OnLine(p3, line1));
        EXPECT_TRUE(!line1.IsInfinity());
        XTLog(std::cout) << line1 << std::endl;
    }


}

TEST(Euclidean3, Plane)
{
    using namespace xiaotu::math;
    
    {
        Plane3<double> plane1({0, 0, 1}, -1);
        XTLog(std::cout) << plane1 << std::endl;
        EXPECT_TRUE(OnPlane({0, 0, 1}, plane1));
        EXPECT_TRUE(plane1.IsValid());

        Plane3<double> plane2({0, 0, 1}, {0, 1, 1});
        XTLog(std::cout) << plane2 << std::endl;
        EXPECT_TRUE(plane1 == plane2);
        EXPECT_TRUE(plane2.IsValid());
    }
}

TEST(Euclidean3, PlaneMeetPlane)
{
    using namespace xiaotu::math;

    Plane3<double> plane1({0, 0, 1}, {0, 1, 1});
    Plane3<double> plane2({0, 0, 1}, {0, 1, 0});
    
    Line3<double> line1 = Meet(plane1, plane2);
    EXPECT_TRUE(line1.IsInfinity());
    XTLog(std::cout) << line1 << std::endl;
}

 

