#include <XiaoTuMathBox/Circle.h>


namespace xiaotu {
namespace math {
    int CirclesIntersection(Point2D const & c1, double r1,
                            Point2D const & c2, double r2,
                            std::vector<Point2D> & intersections)
    {
        intersections.clear();
        double d = Distance(c1, c2);

        // 两圆相隔太远, 不相交
        if (d > (r1 + r2))
            return 0;

        // 两圆圆心基本重叠，若半径一样则有无穷多解
        if (d < 1e-9)
            return 0;

        // 大圆套小圆, 不相交
        if (d < std::abs(r1 - r2))
            return 0;

        double inv_d = 1.0 / d;
        // 两圆相切
        if (d == (r1 + r2)) {
            intersections.push_back(c1 + r1 * inv_d * (c2 - c1));
            return 1;
        }

        // 两个交点分别为p1, p2, 弦p1p2与c1c2的交点为p3
        // h = | p1 - p2 | / 2
        // a = | c1 - p3 |
        // b = | c2 - p3 |
        // d = a + b
        // a^2 + h^2 = r1^2
        // b^2 + h^2 = r2^2
        double r1r1 = r1 * r1;
        double r2r2 = r2 * r2;
        double a = (r1r1 - r2r2 + d * d) * 0.5 * inv_d;
        double h = sqrt(r1r1 - a * a);

        Point2D diff_c = c2 - c1;
        Point2D p3 = c1 + a * inv_d * diff_c;
        intersections.push_back(Point2D(p3.x() + h * inv_d * diff_c.y(), p3.y() - h * inv_d * diff_c.x()));
        intersections.push_back(Point2D(p3.x() - h * inv_d * diff_c.y(), p3.y() + h * inv_d * diff_c.x()));

        return 2;
    }
}
}

