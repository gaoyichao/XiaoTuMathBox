#ifndef XTMB_POINT2D_H
#define XTMB_POINT2D_H

#include <cmath>
#include <Eigen/Eigen>

namespace xiaotu {
namespace math {

    typedef Eigen::Vector2d Point2D;
    typedef Eigen::Vector2d Vector2D;

    inline double Distance(Point2D const & p1, Point2D const & p2)
    {
        Vector2D diff = p2 - p1;
        return diff.norm();
    }

    void UnionVector2D(Point2D const & from, Point2D const & to, Vector2D & vec);
}
}

#endif

