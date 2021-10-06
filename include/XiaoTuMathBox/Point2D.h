#ifndef XTMB_POINT2D_H
#define XTMB_POINT2D_H

#include <cmath>
#include <Eigen/Eigen>

namespace xiaotu {
namespace math {

    constexpr double CFG_SMALL_VALUE = 1.0e-6;
    constexpr double RAD_TO_DEGREE = 180.0 / 3.1415926;
    constexpr double DEGREE_TO_RAD = 3.1415926 / 180.0;

    typedef Eigen::Vector2d Point2D;
    typedef Eigen::Vector2d Vector2D;

    inline double Distance(Point2D const & p1, Point2D const & p2)
    {
        Vector2D diff = p2 - p1;
        return diff.norm();
    }

    /*
     * Angle - 计算两个向量之间的夹角
     */
    inline double Angle(Vector2D const & v1, Vector2D const & v2)
    {
        return std::atan2(v1.x() * v2.y() - v1.y() * v2.x(),
                          v1.x() * v2.x() + v1.y() * v2.y());
    }

    /*
     * __CrossProduct - 计算 sp->p 叉积 sp->ep, 主要用于判定点 p 在线段 sp->ep 的左侧还是右侧
     */
    inline double __CrossProduct__(Point2D const & sp, Point2D const & ep, Point2D const & p)
    {
        double dx0 = sp.x() -  p.x();
        double dy0 = sp.y() -  p.y();
        double dx1 = ep.x() - sp.x();
        double dy1 = ep.y() - sp.y();
 
        return dx0 * dy1 - dx1 * dy0;
    }
    /*
     * OnRightSide - 判定点 p 是否在线段 sp->ep 的右侧
     */
    inline bool OnRightSide(Point2D const & sp, Point2D const & ep, Point2D const & p)
    {
        return __CrossProduct__(sp, ep, p) < 0;
    }
    /*
     * OnLeftSide - 判定点 p 是否在线段 sp->ep 的左侧
     */
    inline bool OnLeftSide(Point2D const & sp, Point2D const & ep, Point2D const & p)
    {
        return __CrossProduct__(sp, ep, p) > 0;
    }


    void UnionVector2D(Point2D const & from, Point2D const & to, Vector2D & vec);
}
}

#endif

