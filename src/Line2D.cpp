#include <XiaoTuMathBox/Line2D.h>
#include <XiaoTuMathBox/Equation.h>


namespace xiaotu {
namespace math {

    /*
     * IntersectCircle - 与圆的交点
     * 
     * 把直线的参数方程直接带入圆方程，得到关于u的一元二次方程，解出u即可
     *   x = x_0 + u dx
     *   y = y_0 + u dy 
     *   (x - m)^2 + (y - n)^2 = r^2
     */
    int Line2D::IntersectCircle(Point2D const & cp, double r, std::vector<Point2D> & intersections)
    {
        intersections.clear();

        Point2D dir = Dir();

        Point2D p1 = RefPoint();
        Point2D p2 = p1 + dir;
        Point2D const & p3 = cp;

        // au^2 + bu + c = 0
        double a = dir.squaredNorm();
        double b = 2 * (dir.x() * (p1.x() - p3.x()) + dir.y() * (p1.y() - p3.y()));
        double c = p3.squaredNorm() + p1.squaredNorm() - 2 * (p3.x() * p1.x() + p3.y() * p1.y()) - r * r;

        QuadraticEquation<double> equa(a, b, c);
        if (0 == equa.n)
            return 0;

        if (1 == equa.n) {
            intersections.push_back(p1 + equa.x1 * dir);
        } else {
            intersections.push_back(p1 + equa.x1 * dir);
            intersections.push_back(p1 + equa.x2 * dir);
        }

        return intersections.size();
    }

}
}

