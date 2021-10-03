#ifndef XTMB_CIRCLE_H
#define XTMB_CIRCLE_H

#include <cmath>
#include <iostream>
#include <vector>

#include <XiaoTuMathBox/Point2D.h>

namespace xiaotu {
namespace math {


    /*
     * CirclesIntersection - 两圆的交点
     */
    int CirclesIntersection(Point2D const & c1, double r1,
                            Point2D const & c2, double r2,
                            std::vector<Point2D> & intersections);

}
}

#endif

