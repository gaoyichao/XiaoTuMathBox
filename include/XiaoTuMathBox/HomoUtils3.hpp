#ifndef XTMB_HOMOUTILS3_H
#define XTMB_HOMOUTILS3_H

#include <XiaoTuMathBox/HomoPoint3.hpp>
#include <XiaoTuMathBox/HomoPlane3.hpp>
#include <XiaoTuMathBox/HomoLine3.hpp>
#include <XiaoTuMathBox/Projective3.hpp>

namespace xiaotu {
namespace math {

    //! @brief 三点共面
    template <typename DataType>
    inline HomoPlane3<DataType> Coplanar(HomoPoint3<DataType> const & p1, HomoPoint3<DataType> const & p2, HomoPoint3<DataType> const & p3)
    {
        DataType d234 = p1[1] * p2[2] * p3[3] + p2[1] * p3[2] * p1[3] + p3[1] * p1[2] * p2[3]
                      - p3[1] * p2[2] * p1[3] - p2[1] * p1[2] * p3[3] + p1[1] * p3[2] * p2[3];
        DataType d134 = p1[0] * p2[2] * p3[3] + p2[0] * p3[2] * p1[3] + p3[0] * p1[2] * p2[3]
                      - p3[0] * p2[2] * p1[3] - p2[0] * p1[2] * p3[3] + p1[0] * p3[2] * p2[3];
        DataType d124 = p1[0] * p2[1] * p3[3] + p2[0] * p3[1] * p1[3] + p3[0] * p1[1] * p2[3]
                      - p3[0] * p2[1] * p1[3] - p2[0] * p1[1] * p3[3] + p1[0] * p3[1] * p2[3];
        DataType d123 = p1[0] * p2[1] * p3[2] + p2[0] * p3[1] * p1[2] + p3[0] * p1[1] * p2[2]
                      - p3[0] * p2[1] * p1[2] - p2[0] * p1[1] * p3[2] + p1[0] * p3[1] * p2[2];

        return HomoPlane3<DataType>(d234, -d134, d124, -d123);
    }

    //! @brief 三点共面, SVD解零空间
    template <typename DataType>
    inline HomoPlane3<DataType> CoplanarSVD(HomoPoint3<DataType> const & p1, HomoPoint3<DataType> const & p2, HomoPoint3<DataType> const & p3)
    {
        Eigen::Matrix<DataType, 3, 4> m;
        m << p1.transpose(),
             p2.transpose(),
             p3.transpose();
        Eigen::JacobiSVD<Eigen::Matrix<DataType, 3, 4>> svd(m, Eigen::ComputeFullV);
        return HomoPlane3<DataType>(svd.matrixV().col(3));
    }

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

    //! @brief 三面交点 
    template <typename DataType>
    inline HomoPoint3<DataType> Intersection(HomoPlane3<DataType> const & p1, HomoPlane3<DataType> const & p2, HomoPlane3<DataType> const & p3)
    {
        DataType d234 = p1[1] * p2[2] * p3[3] + p2[1] * p3[2] * p1[3] + p3[1] * p1[2] * p2[3]
                      - p3[1] * p2[2] * p1[3] - p2[1] * p1[2] * p3[3] + p1[1] * p3[2] * p2[3];
        DataType d134 = p1[0] * p2[2] * p3[3] + p2[0] * p3[2] * p1[3] + p3[0] * p1[2] * p2[3]
                      - p3[0] * p2[2] * p1[3] - p2[0] * p1[2] * p3[3] + p1[0] * p3[2] * p2[3];
        DataType d124 = p1[0] * p2[1] * p3[3] + p2[0] * p3[1] * p1[3] + p3[0] * p1[1] * p2[3]
                      - p3[0] * p2[1] * p1[3] - p2[0] * p1[1] * p3[3] + p1[0] * p3[1] * p2[3];
        DataType d123 = p1[0] * p2[1] * p3[2] + p2[0] * p3[1] * p1[2] + p3[0] * p1[1] * p2[2]
                      - p3[0] * p2[1] * p1[2] - p2[0] * p1[1] * p3[2] + p1[0] * p3[1] * p2[2];

        return HomoPoint3<DataType>(d234, -d134, d124, -d123);
    }

    //! @brief 三面交点, SVD解零空间
    template <typename DataType>
    inline HomoPoint3<DataType> IntersectionSVD(HomoPlane3<DataType> const & p1, HomoPlane3<DataType> const & p2, HomoPlane3<DataType> const & p3)
    {
        Eigen::Matrix<DataType, 3, 4> m;
        m << p1.transpose(),
             p2.transpose(),
             p3.transpose();
        Eigen::JacobiSVD<Eigen::Matrix<DataType, 3, 4>> svd(m, Eigen::ComputeFullV);
        return HomoPoint3<DataType>(svd.matrixV().col(3));
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


    //! @brief 判定点在平面上
    //!
    //! @param [in] po 目标点
    //! @param [in] pi 目标平面
    //! @param [in] tolerance 判定容忍度
    //! @return true 在，false 不在
    template <typename DataType>
    inline bool OnPlane(HomoPoint3<DataType> const & po, HomoPlane3<DataType> const & pi,  DataType tolerance = 1e-9)
    {
        return (std::abs(po.transpose() * pi) < tolerance);
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

    //! @brief 相同点判定
    template <typename DataType>
    inline bool Equal(HomoPoint3<DataType> const & p1, HomoPoint3<DataType> const & p2, DataType tolerance = 1e-9)
    {
        Eigen::Matrix<DataType, 4, 1> n1 = p1.Normalization();
        Eigen::Matrix<DataType, 4, 1> n2 = p2.Normalization();
        return (std::abs(n1[0] - n2[0]) < tolerance) &&
               (std::abs(n1[1] - n2[1]) < tolerance) &&
               (std::abs(n1[2] - n2[2]) < tolerance) &&
               (std::abs(n1[3] - n2[3]) < tolerance);
    }

    //! @brief 相同点判定
    template <typename DataType>
    inline bool operator == (HomoPoint3<DataType> const & p1, HomoPoint3<DataType> const & p2)
    {
        return Equal(p1, p2);
    }

    //! @brief 相同平面判定
    template <typename DataType>
    inline bool Equal(HomoPlane3<DataType> const & p1, HomoPlane3<DataType> const & p2, DataType tolerance = 1e-9)
    {
        Eigen::Matrix<DataType, 4, 1> n1 = p1.Normalization();
        Eigen::Matrix<DataType, 4, 1> n2 = p2.Normalization();
        return (std::abs(n1[0] - n2[0]) < tolerance) &&
               (std::abs(n1[1] - n2[1]) < tolerance) &&
               (std::abs(n1[2] - n2[2]) < tolerance) &&
               (std::abs(n1[3] - n2[3]) < tolerance);
    }

    //! @brief 相同平面判定
    template <typename DataType>
    inline bool operator == (HomoPlane3<DataType> const & p1, HomoPlane3<DataType> const & p2)
    {
        return Equal(p1, p2);
    }




}
}

#endif
