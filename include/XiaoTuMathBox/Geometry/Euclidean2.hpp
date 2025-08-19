/********************************************************************************************************
 * 
 * 2维欧式空间下的一些几何工具
 *
 **************************************************************************** GAO YiChao 2025.0819 *****/
#ifndef XTMB_GEO_GEOMETRY_H
#error "请勿直接引用 Euclidean2.hpp, 请使用 #include <XiaoTuMathBox/Geometry/Geometry.hpp>"
#endif


namespace xiaotu::math {

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
