#ifndef XTMB_HOMOUTILS3_H
#define XTMB_HOMOUTILS3_H

#include <XiaoTuMathBox/HomoPoint3.hpp>
#include <XiaoTuMathBox/HomoPlane3.hpp>
#include <XiaoTuMathBox/HomoLine3.hpp>
#include <XiaoTuMathBox/Projective3.hpp>

namespace xiaotu {
namespace math {


    //! @brief 点线共面
    template <typename DataType>
    inline HomoPlane3<DataType> Coplanar(HomoPoint3<DataType> const & p, HomoLine3<DataType> const & l)
    {
        return HomoPlane3<DataType>(l.DualForm() * p);
    }

    //! @brief 两点共线 
    template <typename DataType>
    inline HomoLine3<DataType> Collinear(HomoPoint3<DataType> const & p1, HomoPoint3<DataType> const & p2)
    {
        return HomoLine3<DataType>(p1 * p2.transpose() - p2 * p1.transpose());
    }

    //! @brief 两面交线
    template <typename DataType>
    inline HomoLine3<DataType> Intersection(HomoPlane3<DataType> const & p1, HomoPlane3<DataType> const & p2)
    {
        HomoLine3<DataType> dual(p1 * p2.transpose() - p2 * p1.transpose());
        return dual.DualForm();
    }

    //! @brief 线面交点
    template <typename DataType>
    inline HomoPoint3<DataType> Intersection(HomoLine3<DataType> const & l, HomoPlane3<DataType> const & pi)
    {
        return HomoPoint3<DataType>(l * pi);
    }


    //! @brief 判定点在线上
    //!
    //! @param [in] po 目标点
    //! @param [in] l 目标直线
    //! @param [in] tolerance 判定容忍度
    //! @return true 在，false 不在
    template <typename DataType>
    inline bool OnLine(HomoPoint3<DataType> const & po, HomoLine3<DataType> const & l, DataType tolerance = 1e-9)
    {
        Eigen::Matrix<DataType, 4, 1> re = l.DualForm() * po;
        return (std::abs(re[0]) < tolerance) &&
               (std::abs(re[1]) < tolerance) &&
               (std::abs(re[2]) < tolerance) &&
               (std::abs(re[3]) < tolerance);
    }

    //! @brief 判定线在平面上
    //!
    //! @param [in] l  目标直线
    //! @param [in] pi 目标平面
    //! @param [in] tolerance 判定容忍度
    //! @return true 在，false 不在
    template <typename DataType>
    inline bool OnPlane(HomoLine3<DataType> const & l, HomoPlane3<DataType> const & pi,  DataType tolerance = 1e-9)
    {
        return !(Intersection(l, pi).IsValid(tolerance));
    }


}
}

#endif
