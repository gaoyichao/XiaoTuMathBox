#ifndef XTMB_GEO_GEOMETRY_H
#error "请勿直接引用 HomoUtils2.hpp, 请使用 #include <XiaoTuMathBox/Geometry/Geometry.hpp>"
#endif

/////////////////////////////////////////////////////////////////////////
//
// 相同点、直线、圆锥曲线的判定
//
/////////////////////////////////////////////////////////////////////////
namespace xiaotu {

    //! @brief 相同点判定
    template <typename DataType>
    inline bool Equal(HomoPoint2<DataType> const & p1, HomoPoint2<DataType> const & p2, DataType tolerance = SMALL_VALUE)
    {
        auto n1 = p1.Normalization();
        auto n2 = p2.Normalization();
        return (n1 - n2).IsZero(tolerance);
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
    inline bool Equal(HomoLine2<DataType> const & l1, HomoLine2<DataType> const & l2, DataType tolerance = SMALL_VALUE)
    {
        auto n1 = l1.Normalization();
        auto n2 = l2.Normalization();
        return (n1 - n2).IsZero(tolerance) || (n1 + n2).IsZero(tolerance);
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
    inline bool Equal(HomoConic2<DataType> const & c1, HomoConic2<DataType> const & c2, DataType tolerance = SMALL_VALUE)
    {
        HomoConic2<DataType> n1 = c1.Normalization();
        HomoConic2<DataType> n2 = c2.Normalization();
        return (n1 - n2).IsZero(tolerance);
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

/////////////////////////////////////////////////////////////////////////
//
// 点与直线的关系
//
/////////////////////////////////////////////////////////////////////////
namespace xiaotu {

    //! @brief 两点共线
    //!
    //! @param [in] p1 目标点1
    //! @param [in] p2 目标点2
    //! @return 直线对象
    template <typename DataType>
    inline HomoLine2<DataType> Join(HomoPoint2<DataType> const & p1, HomoPoint2<DataType> const & p2)
    {
        return p1.Cross(p2);
    }

    //! @brief 两线交点
    template <typename DataType>
    inline HomoPoint2<DataType> Meet(HomoLine2<DataType> const & l1, HomoLine2<DataType> const & l2)
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
    inline bool OnLine(HomoPoint2<DataType> const & p, HomoLine2<DataType> const & l,  DataType tolerance = SMALL_VALUE)
    {
        return (std::fabs(p.Dot(l)) < tolerance);
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
}

/////////////////////////////////////////////////////////////////////////
//
// 点与圆锥曲线的关系
//
/////////////////////////////////////////////////////////////////////////
namespace xiaotu {

    //! @brief 判定点在圆锥曲线上
    //!
    //! @param [in] p 目标点
    //! @param [in] c 目标圆锥曲线
    //! @param [in] tolerance 判定容忍度
    //! @return true 在，false 不在
    template <typename DataType>
    inline bool OnConic(HomoPoint2<DataType> const & p, HomoConic2<DataType> const & c, DataType tolerance = SMALL_VALUE)
    {
        return (std::fabs((p.Dot(c * p))) < tolerance);
    }

    //! @brief 通过 5 个点求解线性方程组获得圆锥曲线
    template <typename DataType>
    HomoConic2<DataType> Join(HomoPoint2<DataType> const & p0,
                                 HomoPoint2<DataType> const & p1,
                                 HomoPoint2<DataType> const & p2,
                                 HomoPoint2<DataType> const & p3,
                                 HomoPoint2<DataType> const & p4)
    {
        AMatrix<DataType, 5, 6> m;
        m << HomoConic2<DataType>::PointEquation(p0),
             HomoConic2<DataType>::PointEquation(p1),
             HomoConic2<DataType>::PointEquation(p2),
             HomoConic2<DataType>::PointEquation(p3),
             HomoConic2<DataType>::PointEquation(p4);
        
        SVD_GKR svd(m, false, true);
        svd.Iterate(1000, SMALL_VALUE);
        auto re = svd.V().Col(5);
        
        // auto max_indep_set = GaussRowEliminate(m);
        // auto re = SolveSpace(max_indep_set, m);
        
        return HomoConic2(re(0), re(1), re(2), re(3), re(4), re(5));
    }

    //! @brief 通过 5 个点求解线性方程组获得圆锥曲线
    //! 
    //! @param [in] ps 5 个点的起始地址
    template <typename DataType>
    HomoConic2<DataType> Join(HomoPoint2<DataType> const * ps)
    {
        return Join(ps[0], ps[1], ps[2], ps[3], ps[4]);
    }


}


