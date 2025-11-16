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


TEST(Polygon, Diagnalie)
{
    Point2<double> a{0, 0};
    Point2<double> b{1, 0};
    Point2<double> c{0.5, 1};
    Point2<double> d{0, 1};

    {
        Polygon<double> P{a, b, c, d};
        EXPECT_FALSE(P.Diagnal(0, 1)); // a, b
        EXPECT_TRUE(P.Diagnal(0, 2));  // a, c
        EXPECT_FALSE(P.Diagnal(0, 3));  // a, d
    }

}

TEST(Polygon, Triangulate_0)
{
    Polygon<double> P{{0, 0}, {1, 0}, {0.5, 1}, {0, 1}};

    std::vector<size_t> tran_list;
    P.Triangulate(tran_list);
    
    EXPECT_EQ(6, tran_list.size());
    DMatrixView<size_t, eRowMajor> h(tran_list.data(), tran_list.size() / 3, 3);
    XTLog(std::cout) << h << std::endl;
}

TEST(Polygon, Triangulate_1)
{
    Polygon<double> P{{0, 0}, {1, 0}, {1, -10}, {10, 10}};

    std::vector<size_t> tran_list;
    P.Triangulate(tran_list);
    
    EXPECT_EQ(6, tran_list.size());
    DMatrixView<size_t, eRowMajor> h(tran_list.data(), tran_list.size() / 3, 3);
    XTLog(std::cout) << h << std::endl;
}

TEST(Polygon, Triangulate_2)
{

    Polygon<double> P{
        {0,  0},  {10, 7},  {12, 3},  {20, 8},  {13, 17},
        {10, 12}, {12, 14}, {14, 9},  {8,  10}, {6,  14},
        {10, 15}, {7,  17}, {0,  16}, {1,  13}, {3,  15},
        {5,  8},  {-2, 9},  {5,  5}
    };

    std::vector<size_t> tran_list;
    P.Triangulate(tran_list);
    
    EXPECT_EQ(0, tran_list.size() % 3);
    DMatrixView<size_t, eRowMajor> h(tran_list.data(), tran_list.size() / 3, 3);
    XTLog(std::cout) << "三角形数量: " << tran_list.size() / 3 << std::endl;
    XTLog(std::cout) << h << std::endl;
}

