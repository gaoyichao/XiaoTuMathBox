/*********************************************************************************************
 * 
 * 2维欧式空间下的一些几何工具
 *
 ***************************************************************** GAO YiChao 2025.0819 *****/
#ifndef XTMB_GEO_GEOMETRY_H
#error "请勿直接引用 Euclidean2.hpp, 请使用 #include <XiaoTuMathBox/Geometry/Geometry.hpp>"
#endif

/////////////////////////////////////////////////////////////////////////
//
// 在《Computational Geometry in C》的第一章中，给出了这些函数
// 原文是按照整型实现的。个人认为，大多数计算场景都是浮点数。
// 所以这里虽然用模板类来实现，但并未认真考虑过整型的问题。
// http://gaoyichao.com/Xiaotu/?book=计算几何_in_C&title=index
//
/////////////////////////////////////////////////////////////////////////
namespace xiaotu {

    /**
     * @brief 给定三个点，计算三角形 a->b->c 的面积的两倍
     */
    template <typename DataType>
    DataType Area2(Point2<DataType> const & a,
                   Point2<DataType> const & b,
                   Point2<DataType> const & c)
    {
        return (b.x() - a.x()) * (c.y() - a.y()) - (c.x() - a.x()) * (b.y() - a.y());
    }

    /**
     * @brief 给定三个点，计算三角形 a->b->c 的面积
     */
    template <typename DataType>
    DataType Area(Point2<DataType> const & a,
                  Point2<DataType> const & b,
                  Point2<DataType> const & c)
    {
        return 0.5 * Area2(a, b, c);
    }

    /**
     * @brief 判定点 c 是否在线段 ab 的左侧
     */
    template <typename DataType>
    bool Left(Point2<DataType> const & a,
              Point2<DataType> const & b,
              Point2<DataType> const & c,
              DataType tolerance = SMALL_VALUE)
    {
        return Area2(a,b,c) > tolerance;
    }

    /**
     * @brief 判定点 c 是否在线段 ab 的左侧，或者三点共线
     */
    template <typename DataType>
    bool LeftOn(Point2<DataType> const & a,
                Point2<DataType> const & b,
                Point2<DataType> const & c,
                DataType tolerance = SMALL_VALUE)
    {
        return Area2(a,b,c) >= -tolerance;
    }

    /**
     * @brief 判定三点是否共线
     */
    template <typename DataType>
    bool Collinear(Point2<DataType> const & a,
                Point2<DataType> const & b,
                Point2<DataType> const & c,
                DataType tolerance = SMALL_VALUE)
    {
        return std::abs(Area2(a,b,c)) <= tolerance;
    }


    /**
     * @brief 判定线段 ab 与 cd 是否正常相交，即，交点在两个线段的内部
     */
    template <typename DataType>
    bool IntersectProp(Point2<DataType> const & a,
                       Point2<DataType> const & b,
                       Point2<DataType> const & c,
                       Point2<DataType> const & d,
                       DataType tolerance = SMALL_VALUE)
    {
        DataType abc = Area2(a,b,c);
        DataType abd = Area2(a,b,d);
        DataType cda = Area2(c,d,a);
        DataType cdb = Area2(c,d,b);

        //! 任意三个端点共线
        if (std::abs(abc) <= tolerance || std::abs(abd) <= tolerance ||
            std::abs(cda) <= tolerance || std::abs(cdb) <= tolerance)
            return false;

        return XOR(abc > tolerance, abd > tolerance) &&
               XOR(cda > tolerance, cdb > tolerance);
    }

    /**
     * @brief 判定点 c 是否在 ab 之间
     */
    template <typename DataType>
    bool Between(Point2<DataType> const & a,
                 Point2<DataType> const & b,
                 Point2<DataType> const & c,
                 DataType tolerance = SMALL_VALUE)
    {
        if (!Collinear(a,b,c, tolerance))
            return false;

        DataType dx = std::abs(a.x() - b.x());
        DataType dy = std::abs(a.y() - b.y());
        if (dx > dy) {
            return ((a.x() <= c.x()) && (c.x() <= b.x())) || ((a.x() >= c.x()) && (c.x() >= b.x()));
        } else {
            return ((a.y() <= c.y()) && (c.y() <= b.y())) || ((a.y() >= c.y()) && (c.y() >= b.y()));
        }
    }

    /**
     * @brief 判定线段 ab 与 cd 是否相交
     * 
     * 1. 交点在两个线段的内部
     * 2. 一条线段的一个端点位于另一条线段的两个端点之间
     */
    template <typename DataType>
    bool Intersect(Point2<DataType> const & a,
                   Point2<DataType> const & b,
                   Point2<DataType> const & c,
                   Point2<DataType> const & d,
                   DataType tolerance = SMALL_VALUE)
    {
        if (IntersectProp(a,b,c,d,tolerance))
            return true;
        
        if (Between(a,b,c) || Between(a,b,d) || Between(c,d,a) || Between(c,d,b))
            return true;

        return false;
    }
}


namespace xiaotu {

    //! @brief 计算两个点之间的距离
    template <typename DataType>
    inline DataType Distance(Point2<DataType> const & p0, Point2<DataType> const & p1)
    {
        Vector2<DataType> v0(p1 - p0);
        return v0.Norm();
    }

    //! @brief 计算两个向量夹角的余弦值
    template <typename DataType>
    inline DataType Cos(Vector2<DataType> const & v0,  Vector2<DataType> const & v1)
    {
        DataType dot = v0.Dot(v1);
        return dot / v0.Norm() / v1.Norm();
    }

    //! @brief 计算两个向量夹角的正弦值
    template <typename DataType>
    inline DataType Sin(Vector2<DataType> const & v0,  Vector2<DataType> const & v1)
    {
        DataType cross = v0.Cross(v1);
        return cross / v0.Norm() / v1.Norm();
    }

    //! @brief 计算两个向量夹角的弧度值
    template <typename DataType>
    inline DataType Radian(Vector2<DataType> const & v0,  Vector2<DataType> const & v1)
    {
        DataType dot = v0.Dot(v1);
        DataType cross = v0.Cross(v1);
        return std::atan2(cross, dot);
    }

    //! @brief 计算两个向量夹角的余弦值, p0 -> p1, p1 -> p2
    template <typename DataType>
    inline DataType Cos(Point2<DataType> const & p0, Point2<DataType> const & p1, Point2<DataType> const & p2)
    {
        Vector2<DataType> v0(p1 - p0);
        Vector2<DataType> v1(p2 - p1);
        return Cos(v0, v1);
    }

    //! @brief 计算两个向量夹角的正弦值, p0 -> p1, p1 -> p2
    template <typename DataType>
    inline DataType Sin(Point2<DataType> const & p0, Point2<DataType> const & p1, Point2<DataType> const & p2)
    {
        Vector2<DataType> v0(p1 - p0);
        Vector2<DataType> v1(p2 - p1);
        return Sin(v0, v1);
    }


    //! @brief 计算两个向量夹角的弧度值, p0 -> p1, p1 -> p2
    template <typename DataType>
    inline DataType Radian(Point2<DataType> const & p0, Point2<DataType> const & p1, Point2<DataType> const & p2)
    {
        Vector2<DataType> v0(p1 - p0);
        Vector2<DataType> v1(p2 - p1);
        return Radian(v0, v1);
    }

}






