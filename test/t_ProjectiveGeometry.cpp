#include <iostream>

#include <XiaoTuMathBox/ProjectiveGeometry.h>

#include <gtest/gtest.h>

TEST(Homogeneous, Point2D)
{
    using namespace xiaotu::math;

    HomoPoint2D bienao(0.1, 1.0, 2.0);
    EXPECT_EQ(&bienao[0], &bienao.x);
    EXPECT_EQ(&bienao[1], &bienao.y);
    EXPECT_EQ(&bienao[2], &bienao.k);

    HomoPoint2D douniwan = bienao.Normalize();
    EXPECT_EQ(1.0, bienao.k);
    EXPECT_EQ(0.05, bienao.x);
    EXPECT_EQ(0.5, bienao.y);
    EXPECT_EQ(1.0, douniwan.k);
    EXPECT_EQ(0.05, douniwan.x);
    EXPECT_EQ(0.5, douniwan.y);

    bienao.SetValue(0.1, 1.0, 2.0);
    EXPECT_TRUE(bienao == douniwan);

    bienao.SetValue(1.0, 0.0, 1.0);
    douniwan.SetValue(0.0, 1.0, 1.0);
    HomoLine2D l1 = bienao.JoinLine(douniwan);
    HomoLine2D l2 = douniwan.JoinLine(bienao);
    EXPECT_EQ(l1, l2);

    HomoPoint2D p0 = HomoPoint2D::Infinity();
    EXPECT_TRUE(p0.IsInfinity());
}

TEST(Homogeneous, Line2D)
{
    using namespace xiaotu::math;

    HomoLine2D l(0.1, 1.0, 0.0);
    EXPECT_EQ(&l[0], &l.a);
    EXPECT_EQ(&l[1], &l.b);
    EXPECT_EQ(&l[2], &l.c);

    HomoPoint2D p(0.0, 0.0, 1.0);
    EXPECT_TRUE(IsPointOnLine(p, l));

    HomoLine2D l0(1.1, 1.0, 1.0);
    p = l.Intersection(l0);
    EXPECT_TRUE(IsPointOnLine(p, l0));
    EXPECT_TRUE(IsPointOnLine(p, l));
    p.Normalize();
    EXPECT_EQ(1.0, p.k);

    HomoPoint2D p0 = l0.Intersection(l);
    EXPECT_TRUE(IsPointOnLine(p0, l0));
    EXPECT_TRUE(IsPointOnLine(p0, l));
    p0.Normalize();
    EXPECT_EQ(1.0, p0.k);

    EXPECT_TRUE(p == p0);
    EXPECT_EQ(p, p0);
}

TEST(Homogeneous, IdealPoints)
{
    using namespace xiaotu::math;

    HomoLine2D l0(0.1, 1.0, 0.0);
    HomoLine2D l1(0.1, 1.0, 1.0);
    EXPECT_TRUE(HomoPoint2D(l1.Intersection(l0)).IsInfinity());

    HomoLine2D l2 = HomoLine2D::Infinity();
    EXPECT_TRUE(l2.IsInfinity());
    EXPECT_TRUE(HomoPoint2D(l2.Intersection(l0)).IsInfinity());
    EXPECT_TRUE(HomoPoint2D(l2.Intersection(l1)).IsInfinity());

    HomoPoint2D p0 = l0.Intersection(l2);
    EXPECT_EQ(l0.a, -p0.y);
    EXPECT_EQ(l0.b, p0.x);
    EXPECT_EQ(0.0, p0.k);
}


