#include <iostream>

#include <XiaoTuMathBox/Geometry/Geometry.hpp>
#include <XiaoTuMathBox/Geometry/Polygon.hpp>

#include <gtest/gtest.h>

using namespace xiaotu;

TEST(Polygon, TriaArea)
{
    Point2<double> a{0, 0};
    Point2<double> b{1, 0};
    Point2<double> c{0.5, 1};

    {
        auto area = Area2(a, b, c);
        EXPECT_DOUBLE_EQ(1.0, area);

        Polygon<double> P{a, b, c};
        EXPECT_DOUBLE_EQ(1.0, P.Area2());
        EXPECT_DOUBLE_EQ(0.5, P.Area());
    }

    {
        auto area = Area2(c, b, a);
        EXPECT_DOUBLE_EQ(-1.0, area);

        Polygon<double> P{c, b, a};
        EXPECT_DOUBLE_EQ(-1.0, P.Area2());
        EXPECT_DOUBLE_EQ(-0.5, P.Area());
    }

    {
        Polygon<double> P{b, a};
        EXPECT_DOUBLE_EQ(0.0, P.Area2());
    }
}

TEST(Polygon, QuadArea)
{
    Point2<double> a{0, 0};
    Point2<double> b{1, 0};
    Point2<double> c{0.5, 1};
    Point2<double> d{0, 1};

    {
        Polygon<double> P{a, b, c, d};
        EXPECT_DOUBLE_EQ(1.5, P.Area2());
        EXPECT_DOUBLE_EQ(0.75, P.Area());
    }

    {
        Polygon<double> P{d, c, b, a};
        EXPECT_DOUBLE_EQ(-1.5, P.Area2());
        EXPECT_DOUBLE_EQ(-0.75, P.Area());
    }
}

TEST(Polygon, SegmentIntersection)
{
    {
        Point2<double> a{0, 0};
        Point2<double> b{1, 1};
        Point2<double> c{0.5, 0.5};
        Point2<double> d{0, 1};
        EXPECT_TRUE(Intersect(a, b, c, d));
    }

    {
        Point2<double> a{0, 0};
        Point2<double> b{1, 1};
        Point2<double> c{0.5, 0.6};
        Point2<double> d{0, 1};
        EXPECT_FALSE(Intersect(a, b, c, d));
    }

    {
        Point2<double> a{0, 0};
        Point2<double> b{1, 1};
        Point2<double> c{0.5, 0.4};
        Point2<double> d{0, 1};
        EXPECT_TRUE(Intersect(a, b, c, d));
    }
}




