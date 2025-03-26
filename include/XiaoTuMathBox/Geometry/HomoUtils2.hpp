#ifndef XTMB_HOMOUTILS2_H
#define XTMB_HOMOUTILS2_H

#include <XiaoTuMathBox/Geometry/HomoPoint2.hpp>
#include <XiaoTuMathBox/Geometry/HomoLine2.hpp>
#include <XiaoTuMathBox/Geometry/HomoConic2.hpp>
// #include <XiaoTuMathBox/Projective2.hpp>

namespace xiaotu {
namespace math {


    //! @brief 两点共线
    //!
    //! @param [in] p1 目标点1
    //! @param [in] p2 目标点2
    //! @return 直线对象
    template <typename DataType>
    inline HomoLine2<DataType> Collinear(HomoPoint2<DataType> const & p1, HomoPoint2<DataType> const & p2)
    {
        return p1.Cross(p2);
    }

    //! @brief 两线交点
    template <typename DataType>
    inline HomoPoint2<DataType> Intersection(HomoLine2<DataType> const & l1, HomoLine2<DataType> const & l2)
    {
        return l1.Cross(l2);
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
        return (std::fabs(p.Dot(l)) < tolerance);
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
        return (std::fabs((p.Dot(c * p))) < tolerance);
    }

    //! @brief 计算点到直线的距离
    //!
    //! @param [in] p 目标点
    //! @param [in] l 目标直线
    //! @return 点到直线的距离
    template <typename DataType>
    inline DataType Distance(HomoPoint2<DataType> const & p, HomoLine2<DataType> const & l)
    {
        AMatrix<DataType, 3, 1> nl = l.Normalization();
        AMatrix<DataType, 3, 1> np = p.Normalization();
        DataType dis = np.Dot(nl);
        return dis;
    }

    //! @brief 相同点判定
    template <typename DataType>
    inline bool Equal(HomoPoint2<DataType> const & p1, HomoPoint2<DataType> const & p2, DataType tolerance = 1e-9)
    {
        auto n1 = p1.Normalization();
        auto n2 = p2.Normalization();
        return (std::fabs(n1(0) - n2(0)) < tolerance) &&
               (std::fabs(n1(1) - n2(1)) < tolerance) &&
               (std::fabs(n1(2) - n2(2)) < tolerance);
    }

    //! @brief 相同点判定
    template <typename DataType>
    inline bool operator == (HomoPoint2<DataType> const & p1, HomoPoint2<DataType> const & p2)
    {
        return Equal(p1, p2);
    }

    //! @brief 不同点判定
    template <typename DataType>
    inline bool operator != (HomoPoint2<DataType> const & p1, HomoPoint2<DataType> const & p2)
    {
        return !Equal(p1, p2);
    }

    //! @brief 相同直线判定
    template <typename DataType>
    inline bool Equal(HomoLine2<DataType> const & l1, HomoLine2<DataType> const & l2, DataType tolerance = 1e-9)
    {
        auto n1 = l1.Normalization();
        auto n2 = l2.Normalization();
        return (std::fabs(n1(0) - n2(0)) < tolerance) &&
               (std::fabs(n1(1) - n2(1)) < tolerance) &&
               (std::fabs(n1(2) - n2(2)) < tolerance);
    }

    //! @brief 相同直线判定
    template <typename DataType>
    inline bool operator == (HomoLine2<DataType> const & l1, HomoLine2<DataType> const & l2)
    {
        return Equal(l1, l2);
    }

    //! @brief 不同直线判定
    template <typename DataType>
    inline bool operator != (HomoLine2<DataType> const & l1, HomoLine2<DataType> const & l2)
    {
        return !Equal(l1, l2);
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

    //! @brief 不同圆锥曲线判定
    template <typename DataType>
    inline bool operator != (HomoConic2<DataType> const & c1, HomoConic2<DataType> const & c2)
    {
        return !Equal(c1, c2);
    }


}
}

#endif
