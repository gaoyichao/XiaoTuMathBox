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




}

