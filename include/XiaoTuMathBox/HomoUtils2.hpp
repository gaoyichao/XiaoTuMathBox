#ifndef XTMB_HOMOUTILS2_H
#define XTMB_HOMOUTILS2_H

#include <XiaoTuMathBox/HomoPoint2.hpp>
#include <XiaoTuMathBox/HomoLine2.hpp>

namespace xiaotu {
namespace math {


    //! @brief 两点共线
    template <typename DataType>
    inline HomoLine2<DataType> Collinear(HomoPoint2<DataType> const & p1, HomoPoint2<DataType> const & p2)
    {
        return p1.cross(p2);
    }

    //! @brief 两线交点
    template <typename DataType>
    inline HomoPoint2<DataType> Intersection(HomoLine2<DataType> const & l1, HomoLine2<DataType> const & l2)
    {
        return l1.cross(l2);
    }

    //! @brief 判定点在线上
    template <typename DataType>
    inline bool OnLine(HomoPoint2<DataType> const & p, HomoLine2<DataType> const & l,  double tolerance = 1e-9)
    {
        return (std::fabs(p.transpose() * l) < tolerance);
    }

    //! @brief 相同点判定
    template <typename DataType>
    inline bool Equal(HomoPoint2<DataType> const & p1, HomoPoint2<DataType> const & p2, double tolerance = 1e-9)
    {
        if (std::fabs(p1.k() - p2.k()) < tolerance) {
            return (std::fabs(p1.x() - p2.x()) < tolerance) &&
                   (std::fabs(p1.y() - p2.y()) < tolerance);
        } else {
            double k1_inv = 1.0 / p1.k();
            double k2_inv = 1.0 / p2.k();
            return (std::fabs(p1.x() * k1_inv - p2.x() * k2_inv) < tolerance) &&
                   (std::fabs(p1.y() * k1_inv - p2.y() * k2_inv) < tolerance);
        }

    }

    //! @brief 相同点判定
    template <typename DataType>
    inline bool operator == (HomoPoint2<DataType> const & p1, HomoPoint2<DataType> const & p2)
    {
        return Equal(p1, p2);
    }

}
}

#endif
