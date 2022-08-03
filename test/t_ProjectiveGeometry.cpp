#include <iostream>
#include <XiaoTuMathBox/HomoPoint2.hpp>
#include <XiaoTuMathBox/HomoLine2.hpp>
#include <XiaoTuMathBox/HomoUtils2.hpp>

#include <gtest/gtest.h>

TEST(Homogeneous, HomoPoint2)
{
    using namespace xiaotu::math;
    EXPECT_EQ(sizeof(Eigen::Vector3d), sizeof(HomoPoint2<double>));

    HomoPoint2<double> p1(3.14159, 1.41421, 2.71828);
    EXPECT_EQ(3.14159, p1.x());
    EXPECT_EQ(1.41421, p1.y());
    EXPECT_EQ(2.71828, p1.k());
    EXPECT_EQ(&p1[0], &p1.x());
    EXPECT_EQ(&p1[1], &p1.y());
    EXPECT_EQ(&p1[2], &p1.k());

    p1.SetValue(0.1, 1.0, 2.0);
    EXPECT_EQ(0.1, p1[0]);
    EXPECT_EQ(1.0, p1[1]);
    EXPECT_EQ(2.0, p1[2]);

    HomoPoint2<double> p2 = p1.Normalize();
    EXPECT_EQ(1.0, p2.k());
    EXPECT_EQ(1.0, p1.k());

    p1 << 0.1, 1.0, 2.0;
    p2 = p1.Normalization();
    EXPECT_EQ(1.0, p2.k());
    EXPECT_EQ(2.0, p1.k());

    EXPECT_FALSE(p1.IsInfinity());
    EXPECT_FALSE(p2.IsInfinity());

    p1 << 0.1, 1.0, 0.0;
    EXPECT_TRUE(p1.IsInfinity());
}

TEST(Homogeneous, HomoLine2)
{
    using namespace xiaotu::math;
    EXPECT_EQ(sizeof(Eigen::Vector3d), sizeof(HomoLine2<double>));

    HomoLine2<double> l1(0.1, 1.0, 0.0);
    EXPECT_EQ(0.1, l1.a());
    EXPECT_EQ(1.0, l1.b());
    EXPECT_EQ(0.0, l1.c());
    EXPECT_EQ(&l1[0], &l1.a());
    EXPECT_EQ(&l1[1], &l1.b());
    EXPECT_EQ(&l1[2], &l1.c());

    l1.SetValue(0.1, 1.0, 2.0);
    EXPECT_EQ(0.1, l1[0]);
    EXPECT_EQ(1.0, l1[1]);
    EXPECT_EQ(2.0, l1[2]);

    HomoLine2<double> l2 = l1.Normalize();
    EXPECT_EQ(1.0, l2.c());
    EXPECT_EQ(1.0, l1.c());

    l1 << 0.1, 1.0, 2.0;
    l2 = l1.Normalization();
    EXPECT_EQ(1.0, l2.c());
    EXPECT_EQ(2.0, l1.c());

    EXPECT_FALSE(l1.IsInfinity());
    EXPECT_FALSE(l2.IsInfinity());

    l1 << 0.0, 0.0, 1.0;
    EXPECT_TRUE(l1.IsInfinity());
}

TEST(Homogeneous, HomoUtils2)
{
    using namespace xiaotu::math;

    HomoPoint2<double> p0, p1;
    p0 << 1.0, 0.0, 1.0;
    p1 << 2.0, 0.0, 1.0;
    HomoLine2<double> l = Collinear(p0, p1);

    HomoPoint2<double> p2(3.0, 0.0, 1.0);
    EXPECT_TRUE(OnLine(p2, l));
    p2 << 1.0, 0.1, 1.0;
    EXPECT_FALSE(OnLine(p2, l));

    HomoLine2<double> l2(1.0, 0.0, 0.0);
    HomoPoint2<double> p = Intersection(l, l2);
    p2 << 0.0, 0.0, 1.0;
    EXPECT_EQ(p, p2);


}


