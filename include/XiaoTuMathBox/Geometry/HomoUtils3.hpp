#ifndef XTMB_GEO_GEOMETRY_H
#error "请勿直接引用 HomoUtils3.hpp, 请使用 #include <XiaoTuMathBox/Geometry/Geometry.hpp>"
#endif

/////////////////////////////////////////////////////////////////////////
//
// 相同点、平面的判定
//
/////////////////////////////////////////////////////////////////////////

namespace xiaotu::math {

    //! @brief 相同点判定
    template <typename DataType>
    inline bool Equal(HomoPoint3<DataType> const & p1, HomoPoint3<DataType> const & p2, DataType tolerance = SMALL_VALUE)
    {
        auto n1 = p1.Normalization();
        auto n2 = p2.Normalization();
        return (std::abs(n1(0) - n2(0)) < tolerance) &&
               (std::abs(n1(1) - n2(1)) < tolerance) &&
               (std::abs(n1(2) - n2(2)) < tolerance) &&
               (std::abs(n1(3) - n2(3)) < tolerance);
    }

    //! @brief 相同点判定
    template <typename DataType>
    inline bool operator == (HomoPoint3<DataType> const & p1, HomoPoint3<DataType> const & p2)
    {
        return Equal(p1, p2);
    }

    //! @brief 不同点判定
    template <typename DataType>
    inline bool operator != (HomoPoint3<DataType> const & p1, HomoPoint3<DataType> const & p2)
    {
        return !Equal(p1, p2);
    }

    //! @brief 相同平面判定
    template <typename DataType>
    inline bool Equal(HomoPlane3<DataType> const & p1, HomoPlane3<DataType> const & p2, DataType tolerance = SMALL_VALUE)
    {
        auto n1 = p1.Normalization();
        auto n2 = p2.Normalization();
        return (std::abs(n1(0) - n2(0)) < tolerance) &&
               (std::abs(n1(1) - n2(1)) < tolerance) &&
               (std::abs(n1(2) - n2(2)) < tolerance) &&
               (std::abs(n1(3) - n2(3)) < tolerance);
    }

    //! @brief 相同平面判定
    template <typename DataType>
    inline bool operator == (HomoPlane3<DataType> const & p1, HomoPlane3<DataType> const & p2)
    {
        return Equal(p1, p2);
    }

    //! @brief 不同平面判定
    template <typename DataType>
    inline bool operator != (HomoPlane3<DataType> const & p1, HomoPlane3<DataType> const & p2)
    {
        return !Equal(p1, p2);
    }

}

/////////////////////////////////////////////////////////////////////////
//
// 点与平面的关系
//
/////////////////////////////////////////////////////////////////////////

namespace xiaotu::math {

    //! @brief 三点共面, SVD解零空间, 不能保证法线的方向
    template <typename DataType>
    inline HomoPlane3<DataType> CoplanarSVD(HomoPoint3<DataType> const & p1, HomoPoint3<DataType> const & p2, HomoPoint3<DataType> const & p3)
    {
        AMatrix<DataType, 3, 4> m;
        m << p1.Transpose(),
             p2.Transpose(),
             p3.Transpose();
        
        SVD_GKR svd(m, false, true);
        svd.Iterate(1000, SMALL_VALUE);
        return svd.V().Col(3);
    }

    //! @brief 三点共面
    template <typename DataType>
    inline HomoPlane3<DataType> Coplanar(HomoPoint3<DataType> const & p1, HomoPoint3<DataType> const & p2, HomoPoint3<DataType> const & p3)
    {
        DataType d234 = p1(1) * p2(2) * p3(3) + p2(1) * p3(2) * p1(3) + p3(1) * p1(2) * p2(3)
                      - p3(1) * p2(2) * p1(3) - p2(1) * p1(2) * p3(3) + p1(1) * p3(2) * p2(3);
        DataType d134 = p1(0) * p2(2) * p3(3) + p2(0) * p3(2) * p1(3) + p3(0) * p1(2) * p2(3)
                      - p3(0) * p2(2) * p1(3) - p2(0) * p1(2) * p3(3) + p1(0) * p3(2) * p2(3);
        DataType d124 = p1(0) * p2(1) * p3(3) + p2(0) * p3(1) * p1(3) + p3(0) * p1(1) * p2(3)
                      - p3(0) * p2(1) * p1(3) - p2(0) * p1(1) * p3(3) + p1(0) * p3(1) * p2(3);
        DataType d123 = p1(0) * p2(1) * p3(2) + p2(0) * p3(1) * p1(2) + p3(0) * p1(1) * p2(2)
                      - p3(0) * p2(1) * p1(2) - p2(0) * p1(1) * p3(2) + p1(0) * p3(1) * p2(2);

        return HomoPlane3<DataType>(d234, -d134, d124, -d123);
    }

    //! @brief 判定点在平面上
    //!
    //! @param [in] po 目标点
    //! @param [in] pi 目标平面
    //! @param [in] tolerance 判定容忍度
    //! @return true 在，false 不在
    template <typename DataType>
    inline bool OnPlane(HomoPoint3<DataType> const & po, HomoPlane3<DataType> const & pi,  DataType tolerance = SMALL_VALUE)
    {
        return (std::abs(po.Dot(pi)) < tolerance);
    }

    //! @brief 三面交点, SVD解零空间
    template <typename DataType>
    inline HomoPoint3<DataType> IntersectionSVD(HomoPlane3<DataType> const & p1, HomoPlane3<DataType> const & p2, HomoPlane3<DataType> const & p3)
    {
        AMatrix<DataType, 3, 4> m;
        m << p1.Transpose(),
             p2.Transpose(),
             p3.Transpose();
 
        SVD_GKR svd(m, false, true);
        svd.Iterate(1000, SMALL_VALUE);
        return svd.V().Col(3);
    }

    //! @brief 三面交点 
    template <typename DataType>
    inline HomoPoint3<DataType> Intersection(HomoPlane3<DataType> const & p1, HomoPlane3<DataType> const & p2, HomoPlane3<DataType> const & p3)
    {
        DataType d234 = p1(1) * p2(2) * p3(3) + p2(1) * p3(2) * p1(3) + p3(1) * p1(2) * p2(3)
                      - p3(1) * p2(2) * p1(3) - p2(1) * p1(2) * p3(3) + p1(1) * p3(2) * p2(3);
        DataType d134 = p1(0) * p2(2) * p3(3) + p2(0) * p3(2) * p1(3) + p3(0) * p1(2) * p2(3)
                      - p3(0) * p2(2) * p1(3) - p2(0) * p1(2) * p3(3) + p1(0) * p3(2) * p2(3);
        DataType d124 = p1(0) * p2(1) * p3(3) + p2(0) * p3(1) * p1(3) + p3(0) * p1(1) * p2(3)
                      - p3(0) * p2(1) * p1(3) - p2(0) * p1(1) * p3(3) + p1(0) * p3(1) * p2(3);
        DataType d123 = p1(0) * p2(1) * p3(2) + p2(0) * p3(1) * p1(2) + p3(0) * p1(1) * p2(2)
                      - p3(0) * p2(1) * p1(2) - p2(0) * p1(1) * p3(2) + p1(0) * p3(1) * p2(2);

        return HomoPoint3<DataType>(d234, -d134, d124, -d123);
    }





}


