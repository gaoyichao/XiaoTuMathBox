/********************************************************************************************************
 * 
 * 在射影空间\(\mathbb{P}^2\)中，点和线对偶(dual) 
 * 
 * https://gaoyichao.com/Xiaotu/?book=几何&title=2D射影空间中的点直线和圆锥曲线
 *
 **************************************************************************** GAO YiChao 2022.0803 *****/
#ifndef XTMB_HOMOUTILS2_H
#error "请勿直接引用 HomoPoint2.hpp, 请使用 #include <XiaoTuMathBox/HomoUtils2.hpp>"
#endif

#ifndef XTMB_HOMOPOINT2_H
#define XTMB_HOMOPOINT2_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Eigen>

namespace xiaotu {
namespace math {

    template <typename DataType>
    class HomoLine2;

    template <typename DataType>
    class HomoPoint2 : public Eigen::Matrix<DataType, 3, 1>
    {
        typedef Eigen::Matrix<DataType, 3, 1> EigenVector;
        
        public:
            HomoPoint2()
            {
                SetValue(0.0, 0.0, 1.0);
            }

            HomoPoint2(DataType _x, DataType _y, DataType _k = 1.0)
            {
                SetValue(_x, _y, _k);
            }

            HomoPoint2(EigenVector const & v)
            {
                SetValue(v[0], v[1], v[2]);
            }

            HomoPoint2(HomoPoint2 const & p)
            {
                SetValue(p.x(), p.y(), p.k());
            }

            HomoPoint2 & operator = (HomoPoint2 const & p)
            {
                SetValue(p.x(), p.y(), p.k());
                return *this;
            }

            HomoPoint2 & operator = (EigenVector const & v)
            {
                SetValue(v[0], v[1], v[2]);
                return *this;
            }

            inline void SetValue(DataType _x, DataType _y, DataType _k)
            {
                x() = _x;
                y() = _y;
                k() = _k;
            }

            HomoPoint2 & Normalize()
            {
                if (0 == k())
                    return *this;

                DataType k_inv = 1.0 / k();
                x() = x() * k_inv;
                y() = y() * k_inv;
                k() = 1.0;
                return *this;
            }

            EigenVector Normalization() const
            {
                if (0 == k())
                    return *this;

                DataType k_inv = 1.0 / k();
                return *this * k_inv;
            }

            inline bool IsInfinity()
            {
                return 0 == k();
            }

            static HomoPoint2 Infinity()
            {
                return HomoPoint2(1, 0, 0);
            }

        public:
            DataType & x() { return this->data()[0]; }
            DataType & y() { return this->data()[1]; }
            DataType & k() { return this->data()[2]; }
            DataType const & x() const { return this->data()[0]; }
            DataType const & y() const { return this->data()[1]; }
            DataType const & k() const { return this->data()[2]; }
    };

}
}

#endif
