#ifndef XTMB_GEO_DECLARATIONS_H
#define XTMB_GEO_DECLARATIONS_H

#include <cmath>
#include <iostream>
#include <vector>

#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

/////////////////////////////////////////////////////////////////////////
//
// 二维欧式空间
//
/////////////////////////////////////////////////////////////////////////
namespace xiaotu::math {

    /**
     * @brief 二维欧式空间下的向量
     */
    template <typename DataType>
    class Vector2;

    /**
     * @brief 二维欧式空间下的点
     */
    template <typename DataType>
    using Point2 = Vector2<DataType>;
}

/////////////////////////////////////////////////////////////////////////
//
// 二维射影空间
//
/////////////////////////////////////////////////////////////////////////
namespace xiaotu::math {
    /**
     * @brief 二维射影空间下的点
     */
    template <typename DataType>
    class HomoPoint2;

    /**
     * @brief 二维射影空间下的直线
     */
    template <typename DataType>
    class HomoLine2;

    /**
     * @brief 二维射影空间下的圆锥曲线, 3x3对称矩阵
     */
    template <typename DataType>
    class HomoConic2;
}


/////////////////////////////////////////////////////////////////////////
//
// 三维射影空间
//
/////////////////////////////////////////////////////////////////////////
namespace xiaotu::math {
    /**
     * @brief 三维射影空间下的点
     */
    template <typename DataType>
    class HomoPoint3;
}



#endif


