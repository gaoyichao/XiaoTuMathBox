#ifndef XTMB_PROJECTIVE_GEOMETRY_H
#define XTMB_PROJECTIVE_GEOMETRY_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Eigen>

#include <XiaoTuMathBox/HomoPoint2D.h>
#include <XiaoTuMathBox/HomoLine2D.h>

namespace xiaotu {
namespace math {

    inline bool IsPointOnLine(HomoPoint2D const & p, HomoLine2D const & l, double tolerance = 1e-9)
    {
        return (std::fabs(p.transpose() * l) < tolerance);
    }

    bool operator == (HomoPoint2D const & p1, HomoPoint2D const & p2)
    {
        if (std::fabs(p1.k - p2.k) < 1e-9) {
            return (std::fabs(p1.x - p2.x) < 1e-9) &&
                   (std::fabs(p1.y - p2.y) < 1e-9);
        } else {
            double k1_inv = 1.0 / p1.k;
            double k2_inv = 1.0 / p2.k;
            return (std::fabs(p1.x * k1_inv - p2.x * k2_inv) < 1e-9) &&
                   (std::fabs(p1.y * k1_inv - p2.y * k2_inv) < 1e-9);
        }
    }

    bool operator == (HomoLine2D const & l1, HomoLine2D const & l2)
    {
        if (std::fabs(l1.c - l2.c) < 1e-9) {
            return (std::fabs(l1.a - l2.a) < 1e-9) &&
                   (std::fabs(l1.b - l2.b) < 1e-9);
        } else {
            double c1_inv = 1.0 / l1.c;
            double c2_inv = 1.0 / l2.c;
            return (std::fabs(l1.a * c1_inv - l2.a * c2_inv) < 1e-9) &&
                   (std::fabs(l1.b * c1_inv - l2.b * c2_inv) < 1e-9);
        }

    }
}
}


#endif
