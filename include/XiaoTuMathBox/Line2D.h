#ifndef XTMB_LINE2D_H
#define XTMB_LINE2D_H

#include <cmath>
#include <iostream>
#include <vector>

#include <XiaoTuMathBox/Point2D.h>

namespace xiaotu {
namespace math {

    /*
     * ax + by + c = 0
     */
    class Line2D {
        public:
            Line2D() : a(0), b(0), c(0) {}

            Line2D(double _a, double _b, double _c) : a(_a), b(_b), c(_c) {}

            Line2D(Line2D const & l) : a(l.a), b(l.b), c(l.c) {}

            /*
             * Line2D - 根据一个点和方向给出直线方程
             * 
             * 需注意，Point2D和Vector2D是同一种数据类型
             */
            Line2D(Point2D const & p, Vector2D const & dir)
            {
                a = -dir.y();
                b = dir.x();
                c = -a * p.x() - b * p.y();
            }

            Line2D & operator = (Line2D const & l)
            {
                a = l.a;
                b = l.b;
                c = l.c;
                return *this;
            }

            /*
             * IntersectCircle - 与圆的交点
             * 
             * @cp: 圆心坐标
             * @r: 圆半径
             * @intersections: 交点坐标
             * return: 交点数量
             */
            int IntersectCircle(Point2D const & cp, double r, std::vector<Point2D> & intersections);

            double y(double _x) const { return (-c - a * _x) / b; }
            double x(double _y) const { return (-c - b * _y) / a; }

            Point2D RefPoint() const
            {
                if (b == 0)
                    return Point2D(-c / a, 0);

                return Point2D(0, -c / b);
            }

            Vector2D Dir() const { return Point2D(b, -a); }
            
            friend std::ostream & operator << (std::ostream & stream, Line2D const & line)
            {
                stream << line.a << "x + " << line.b << "y + " << line.c << " = 0";
                return stream;
            }

        public:
            double a;
            double b;
            double c;
    };

}
}

#endif

