/********************************************************************************************************
 * 
 * 3维欧式空间下的一些几何工具
 *
 **************************************************************************** GAO YiChao 2025.0826 *****/
#ifndef XTMB_GEO_GEOMETRY_H
#error "请勿直接引用 Euclidean2.hpp, 请使用 #include <XiaoTuMathBox/Geometry/Geometry.hpp>"
#endif

#include <XiaoTuMathBox/Geometry/Geometry.hpp>

namespace xiaotu::math {

    //! @brief 相同直线判定
    template <typename DataType>
    bool Equal(Line3<DataType> const & l1, Line3<DataType> const & l2, DataType tolerance = SMALL_VALUE)
    {
        auto n1 = l1.Normalization();
        auto n2 = l2.Normalization();

        if ((n1.Direction() - n2.Direction()).IsZero(tolerance) && (n1.Moment() - n2.Moment()).IsZero(tolerance))
            return true;
        if ((n1.Direction() + n2.Direction()).IsZero(tolerance) && (n1.Moment() + n2.Moment()).IsZero(tolerance))
            return true;
        return false;
    }

    //! @brief 相同直线判定
    template <typename DataType>
    bool operator == (Line3<DataType> const & l1, Line3<DataType> const & l2)
    {
        return Equal(l1, l2);
    }

    //! @brief 不同直线判定
    template <typename DataType>
    bool operator != (Line3<DataType> const & l1, Line3<DataType> const & l2)
    {
        return !Equal(l1, l2);
    }

    //! @brief 两点共线 
    template <typename DataType>
    inline Line3<DataType> Join(Point3<DataType> const & p0, Point3<DataType> const & p1)
    {
        Vector3<DataType> dir = p1 - p0;
        Vector3<DataType> mem = p0.Cross(dir);
        return Line3<DataType>(dir, mem);
    }


    //! @brief 判定点在线上
    //!
    //! @param [in] p 目标点
    //! @param [in] l 目标直线
    //! @param [in] tolerance 判定容忍度
    //! @return true 在，false 不在
    template <typename DataType>
    bool OnLine(Point3<DataType> const & p, Line3<DataType> const & l, DataType tolerance = SMALL_VALUE)
    {
        auto m = p.Cross(l.Direction());
        return Equal(m, l.Moment(), tolerance);
    }



}

