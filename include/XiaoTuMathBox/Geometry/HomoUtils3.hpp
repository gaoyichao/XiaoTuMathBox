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
    bool Equal(HomoPoint3<DataType> const & p1, HomoPoint3<DataType> const & p2, DataType tolerance = SMALL_VALUE)
    {
        auto n1 = p1.Normalization();
        auto n2 = p2.Normalization();
        return (n1 - n2).IsZero(tolerance);
    }

    //! @brief 相同点判定
    template <typename DataType>
    bool operator == (HomoPoint3<DataType> const & p1, HomoPoint3<DataType> const & p2)
    {
        return Equal(p1, p2);
    }

    //! @brief 不同点判定
    template <typename DataType>
    bool operator != (HomoPoint3<DataType> const & p1, HomoPoint3<DataType> const & p2)
    {
        return !Equal(p1, p2);
    }

    //! @brief 相同平面判定
    template <typename DataType>
    bool Equal(HomoPlane3<DataType> const & p1, HomoPlane3<DataType> const & p2, DataType tolerance = SMALL_VALUE)
    {
        auto n1 = p1.Normalization();
        auto n2 = p2.Normalization();
        return (n1 - n2).IsZero(tolerance) || (n1 + n2).IsZero(tolerance);
    }

    //! @brief 相同平面判定
    template <typename DataType>
    bool operator == (HomoPlane3<DataType> const & p1, HomoPlane3<DataType> const & p2)
    {
        return Equal(p1, p2);
    }

    //! @brief 不同平面判定
    template <typename DataType>
    bool operator != (HomoPlane3<DataType> const & p1, HomoPlane3<DataType> const & p2)
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
    HomoPlane3<DataType> JoinSVD(HomoPoint3<DataType> const & p1, HomoPoint3<DataType> const & p2, HomoPoint3<DataType> const & p3)
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
    HomoPlane3<DataType> Join(HomoPoint3<DataType> const & p1, HomoPoint3<DataType> const & p2, HomoPoint3<DataType> const & p3)
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
    bool OnPlane(HomoPoint3<DataType> const & po, HomoPlane3<DataType> const & pi,  DataType tolerance = SMALL_VALUE)
    {
        return (std::abs(po.Dot(pi)) < tolerance);
    }

    //! @brief 三面交点, SVD解零空间
    template <typename DataType>
    HomoPoint3<DataType> MeetSVD(HomoPlane3<DataType> const & p1, HomoPlane3<DataType> const & p2, HomoPlane3<DataType> const & p3)
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
    HomoPoint3<DataType> Meet(HomoPlane3<DataType> const & p1, HomoPlane3<DataType> const & p2, HomoPlane3<DataType> const & p3)
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


/////////////////////////////////////////////////////////////////////////
//
// 普吕克直线与点、平面的关系
//
/////////////////////////////////////////////////////////////////////////

namespace xiaotu::math {

    //! @brief 两点共线 
    template <typename DataType>
    HomoLine3<DataType> Join(HomoPoint3<DataType> const & p1, HomoPoint3<DataType> const & p2)
    {
        return HomoLine3<DataType>(p1 * p2.Transpose() - p2 * p1.Transpose());
    }

    //! @brief 两面交线
    template <typename DataType>
    HomoLine3<DataType> Meet(HomoPlane3<DataType> const & p1, HomoPlane3<DataType> const & p2)
    {
        HomoLine3<DataType> dual(p1 * p2.Transpose() - p2 * p1.Transpose());
        return dual.DualForm();
    }

    //! @brief 线面交点
    //!
    //! @param [in] l 直线
    //! @param [in] pi 平面
    //! @return 交点
    template <typename DataType>
    HomoPoint3<DataType> Meet(HomoLine3<DataType> const & l, HomoPlane3<DataType> const & pi)
    {
        return HomoPoint3<DataType>(l * pi);
    }

    //! @brief 点线共面
    //!
    //! @param [in] p 点
    //! @param [in] l 直线
    //! @return 共面
    template <typename DataType>
    HomoPlane3<DataType> Join(HomoPoint3<DataType> const & p, HomoLine3<DataType> const & l)
    {
        return HomoPlane3<DataType>(l.DualForm() * p);
    }

    //! @brief 判定点在线上
    //!
    //! @param [in] po 目标点
    //! @param [in] l 目标直线
    //! @param [in] tolerance 判定容忍度
    //! @return true 在，false 不在
    template <typename DataType>
    bool OnLine(HomoPoint3<DataType> const & po, HomoLine3<DataType> const & l, DataType tolerance = SMALL_VALUE)
    {
        auto re = l.DualForm() * po;
        return re.IsZero(tolerance);
    }


    //! @brief 判定两条直线是否共面
    template <typename DataType>
    bool AreCoplanar(HomoLine3<DataType> const & l1, HomoLine3<DataType> const & l2, DataType tolerance = SMALL_VALUE)
    {
        double tmp = l1.l12() * l2.l34() + l2.l12() * l1.l34() + l1.l13() * l2.l42()
                   + l2.l13() * l1.l42() + l1.l14() * l2.l23() + l2.l14() * l1.l23();
        return (std::abs(tmp) < tolerance);
    }

    //! @brief 判定线在平面上
    //!
    //! @param [in] l  目标直线
    //! @param [in] pi 目标平面
    //! @param [in] tolerance 判定容忍度
    //! @return true 在，false 不在
    template <typename DataType>
    bool OnPlane(HomoLine3<DataType> const & l, HomoPlane3<DataType> const & pi,  DataType tolerance = SMALL_VALUE)
    {
        return Meet(l, pi).IsZero(tolerance);
    }



}



