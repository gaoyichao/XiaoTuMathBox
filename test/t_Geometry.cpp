#include <iostream>
#include <XiaoTuMathBox/Equation.h>
#include <XiaoTuMathBox/Point2D.h>
#include <XiaoTuMathBox/Line2D.h>
#include <XiaoTuMathBox/Circle.h>

#include <XiaoTuMathBox/CircularValue.h>

#include <gtest/gtest.h>

TEST(Circle, intersections)
{
    xiaotu::math::QuadraticEquation<double> equa(1, 0, -1);
    std::cout << equa.n << std::endl;
    std::cout << equa.x1 << std::endl;
    std::cout << equa.x2 << std::endl;

    xiaotu::math::Point2D p1(2.0, 1.0);
    xiaotu::math::Point2D p2(2.0, 2.0);
    xiaotu::math::Vector2D vec;
    xiaotu::math::UnionVector2D(p1, p2, vec);

    xiaotu::math::Line2D l2(p1, vec);
    std::cout << l2 << std::endl;

    xiaotu::math::Point2D p3(1, 1.5);
    std::vector<xiaotu::math::Point2D> intersections;
    int n = l2.IntersectCircle(p3, 1.5, intersections);

    std::cout << n << std::endl;
    std::cout << "p1:" << intersections[0] << std::endl;
    std::cout << "p2:" << intersections[1] << std::endl;

    intersections.clear();
    n = xiaotu::math::CirclesIntersection(p1, 1, p2, 1, intersections);
    std::cout << n << std::endl;
    std::cout << "p1:" << intersections[0] << std::endl;
    std::cout << "p2:" << intersections[1] << std::endl;
}

TEST(vector, angle)
{
    using namespace xiaotu::math;
    Vector2D v1, v2;
    v1 << 1, 0;
    v2 << 0, 1;
    EXPECT_LT(std::abs(Angle(v1, v2) - 90 * DEGREE_TO_RAD), 1e-6);
    EXPECT_LT(std::abs(Angle(v2, v1) + 90 * DEGREE_TO_RAD), 1e-6);

    v2 << 1, 1;
    EXPECT_LT(std::abs(Angle(v1, v2) - 45 * DEGREE_TO_RAD), 1e-6);
    EXPECT_LT(std::abs(Angle(v2, v1) + 45 * DEGREE_TO_RAD), 1e-6);

    v2 << 2, 2;
    EXPECT_LT(std::abs(Angle(v1, v2) - 45 * DEGREE_TO_RAD), 1e-6);
    EXPECT_LT(std::abs(Angle(v2, v1) + 45 * DEGREE_TO_RAD), 1e-6);

    v2 << -2, 2;
    EXPECT_LT(std::abs(Angle(v1, v2) - 135 * DEGREE_TO_RAD), 1e-6);
    EXPECT_LT(std::abs(Angle(v2, v1) + 135 * DEGREE_TO_RAD), 1e-6);

    v2 << -2, -2;
    EXPECT_LT(std::abs(Angle(v1, v2) + 135 * DEGREE_TO_RAD), 1e-6);
    EXPECT_LT(std::abs(Angle(v2, v1) - 135 * DEGREE_TO_RAD), 1e-6);

    v2 << 1, 0;
    EXPECT_LT(std::abs(Angle(v1, v2) - 0 * DEGREE_TO_RAD), 1e-6);
    EXPECT_LT(std::abs(Angle(v2, v1) + 0 * DEGREE_TO_RAD), 1e-6);

    v2 << -1, 0;
    std::cout << Angle(v1, v2) << std::endl;
    std::cout << Angle(v2, v1) << std::endl;
    EXPECT_LT(std::abs(Angle(v1, v2) - 180 * DEGREE_TO_RAD), 1e-6);
    EXPECT_LT(std::abs(Angle(v2, v1) + 180 * DEGREE_TO_RAD), 1e-6);
}

TEST(point, xxxside)
{
    using namespace xiaotu::math;
    Point2D sp, ep, p;
    sp << 0, 0;
    ep << 1, 1;

    p << 0, 1;
    EXPECT_TRUE(OnLeftSide(sp, ep, p));
    EXPECT_TRUE(OnRightSide(ep, sp, p));

    p << 2, 2;
    std::cout << __CrossProduct__(sp, ep, p) << std::endl;
    EXPECT_FALSE(OnLeftSide(sp, ep, p));
    EXPECT_FALSE(OnRightSide(sp, ep, p));
}

TEST(CircularValue, bienao)
{
    using namespace xiaotu::math;

    std::cout << SignedRadRange::lower_bound << std::endl;
    std::cout << SignedRadRange::upper_bound << std::endl;

    CircularValue<SignedRadRange> cirval;
    std::cout << cirval.GetLowerBound() << std::endl;
    std::cout << cirval.GetUpperBound() << std::endl;
    std::cout << cirval.GetZero() << std::endl;
    std::cout << cirval.GetRange() << std::endl;
    std::cout << cirval.GetHalfRange() << std::endl;

    std::cout << CircularValue<SignedRadRange>::ShortestDist(-3.14, 3.1) << std::endl;
}

