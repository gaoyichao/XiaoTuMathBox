#include <iostream>
#include <XiaoTuMathBox/Equation.h>
#include <XiaoTuMathBox/Point2D.h>
#include <XiaoTuMathBox/Line2D.h>
#include <XiaoTuMathBox/Circle.h>

int main(int argc, char * argv[])
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



    return 0;
}

