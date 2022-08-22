#ifndef XTMB_EUCLIDEAN2_H
#define XTMB_EUCLIDEAN2_H

#include <XiaoTuMathBox/Point2.hpp>
#include <XiaoTuMathBox/Vector2.hpp>


namespace xiaotu {
namespace math {

    //! @brief 计算两个点之间的距离
    template <typename DataType>
    inline DataType Distance(Point2<DataType> const & p0, Point2<DataType> const & p1)
    {
        Vector2<DataType> v0(p1 - p0);
        return v0.Norm();
    }
    //! @brief 计算两个点之间的距离
    template <typename DataType, typename PointType>
    inline DataType Distance(PointType const & p0, PointType const & p1)
    {
        Vector2<DataType> v0(p1.x - p0.x, p1.y - p0.y);
        return v0.Norm();
    }

    //! @brief 计算两个向量夹角的余弦值
    template <typename DataType>
    inline DataType Cos(Vector2<DataType> const & v0,  Vector2<DataType> const & v1)
    {
        DataType dot = v0.transpose() * v1;
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
        DataType dot = v0.transpose() * v1;
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

    //! @brief 计算两个向量夹角的余弦值, p0 -> p1, p1 -> p2
    template <typename DataType, typename PointType>
    inline DataType Cos(PointType const & p0, PointType const & p1, PointType const & p2)
    {
        Vector2<DataType> v0(p1.x - p0.x, p1.y - p0.y);
        Vector2<DataType> v1(p2.x - p1.x, p2.y - p1.y);
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

    //! @brief 计算两个向量夹角的正弦值, p0 -> p1, p1 -> p2
    template <typename DataType, typename PointType>
    inline DataType Sin(PointType const & p0, PointType const & p1, PointType const & p2)
    {
        Vector2<DataType> v0(p1.x - p0.x, p1.y - p0.y);
        Vector2<DataType> v1(p2.x - p1.x, p2.y - p1.y);
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

    //! @brief 计算两个向量夹角的弧度值, p0 -> p1, p1 -> p2
    template <typename DataType, typename PointType>
    inline DataType Radian(PointType const & p0, PointType const & p1, PointType const & p2)
    {
        Vector2<DataType> v0(p1.x - p0.x, p1.y - p0.y);
        Vector2<DataType> v1(p2.x - p1.x, p2.y - p1.y);
        return Radian(v0, v1);
    }


}
}


#endif
