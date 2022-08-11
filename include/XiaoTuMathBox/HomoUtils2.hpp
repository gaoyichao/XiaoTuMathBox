#ifndef XTMB_HOMOUTILS2_H
#define XTMB_HOMOUTILS2_H

#include <XiaoTuMathBox/HomoPoint2.hpp>
#include <XiaoTuMathBox/HomoLine2.hpp>
#include <XiaoTuMathBox/HomoConic2.hpp>
#include <XiaoTuMathBox/Projective2.hpp>

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

    //! @brief 判定点在直线上
    //!
    //! @param [in] p 目标点
    //! @param [in] l 目标直线
    //! @param [in] tolerance 判定容忍度
    //! @return true 在，false 不在
    template <typename DataType>
    inline bool OnLine(HomoPoint2<DataType> const & p, HomoLine2<DataType> const & l,  DataType tolerance = 1e-9)
    {
        return (std::fabs(p.transpose() * l) < tolerance);
    }

    //! @brief 判定点在圆锥曲线上
    //!
    //! @param [in] p 目标点
    //! @param [in] c 目标圆锥曲线
    //! @param [in] tolerance 判定容忍度
    //! @return true 在，false 不在
    template <typename DataType>
    inline bool OnLine(HomoPoint2<DataType> const & p, HomoConic2<DataType> const & c, DataType tolerance = 1e-9)
    {
        return (std::fabs(p.transpose() * c * p) < tolerance);
    }

    //! @brief 相同点判定
    template <typename DataType>
    inline bool Equal(HomoPoint2<DataType> const & p1, HomoPoint2<DataType> const & p2, DataType tolerance = 1e-9)
    {
        Eigen::Matrix<DataType, 3, 1> n1 = p1.Normalization();
        Eigen::Matrix<DataType, 3, 1> n2 = p2.Normalization();
        return (std::fabs(n1[0] - n2[0]) < tolerance) &&
               (std::fabs(n1[1] - n2[1]) < tolerance) &&
               (std::fabs(n1[2] - n2[2]) < tolerance);
    }

    //! @brief 相同点判定
    template <typename DataType>
    inline bool operator == (HomoPoint2<DataType> const & p1, HomoPoint2<DataType> const & p2)
    {
        return Equal(p1, p2);
    }

    //! @brief 相同直线判定
    template <typename DataType>
    inline bool Equal(HomoLine2<DataType> const & l1, HomoLine2<DataType> const & l2, DataType tolerance = 1e-9)
    {
        Eigen::Matrix<DataType, 3, 1> n1 = l1.Normalization();
        Eigen::Matrix<DataType, 3, 1> n2 = l2.Normalization();
        return (std::fabs(n1[0] - n2[0]) < tolerance) &&
               (std::fabs(n1[1] - n2[1]) < tolerance) &&
               (std::fabs(n1[2] - n2[2]) < tolerance);
    }

    //! @brief 相同直线判定
    template <typename DataType>
    inline bool operator == (HomoLine2<DataType> const & l1, HomoLine2<DataType> const & l2)
    {
        return Equal(l1, l2);
    }

    //! @brief 相同圆锥曲线判定
    template <typename DataType>
    inline bool Equal(HomoConic2<DataType> const & c1, HomoConic2<DataType> const & c2, DataType tolerance = 1e-9)
    {
        HomoConic2<DataType> n1 = c1.Normalization();
        HomoConic2<DataType> n2 = c2.Normalization();
        return (std::fabs(n1.a() - n2.a()) < tolerance) &&
               (std::fabs(n1.b() - n2.b()) < tolerance) &&
               (std::fabs(n1.c() - n2.c()) < tolerance) &&
               (std::fabs(n1.d() - n2.d()) < tolerance) &&
               (std::fabs(n1.e() - n2.e()) < tolerance) &&
               (std::fabs(n1.f() - n2.f()) < tolerance);
    }

    //! @brief 相同圆锥曲线判定
    template <typename DataType>
    inline bool operator == (HomoConic2<DataType> const & c1, HomoConic2<DataType> const & c2)
    {
        return Equal(c1, c2);
    }


}
}

#endif
