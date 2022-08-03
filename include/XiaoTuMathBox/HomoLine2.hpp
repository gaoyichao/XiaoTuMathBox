/********************************************************************************************************
 * 
 * 在射影空间\(\mathbb{P}^2\)中，点和线对偶(dual) 
 * 
 * https://gaoyichao.com/Xiaotu/?book=几何&title=2D几何空间和变换
 *
 **************************************************************************** GAO YiChao 2022.0803 *****/
#ifndef XTMB_HOMOLINE2_H
#define XTMB_HOMOLINE2_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Eigen>

namespace xiaotu {
namespace math {

    template <typename DataType>
    class HomoPoint2;

    template <typename DataType>
    class HomoLine2 : public Eigen::Matrix<DataType, 3, 1>
    {
        typedef Eigen::Matrix<DataType, 3, 1> EigenVector;
        public:
            HomoLine2()
            {
                SetValue(0.0, 0.0, 1.0);
            }

            HomoLine2(DataType _a, DataType _b, DataType _c)
            {
                SetValue(_a, _b, _c);
            }

            HomoLine2(Eigen::Vector3d const & v)
            {
                SetValue(v[0], v[1], v[2]);
            }

            HomoLine2(HomoLine2 const & p)
            {
                SetValue(p.a(), p.b(), p.c());
            }

            HomoLine2 & operator = (HomoLine2 const & p)
            {
                SetValue(p.a(), p.b(), p.c());
                return *this;
            }

            HomoLine2 & operator = (EigenVector const & v)
            {
                SetValue(v[0], v[1], v[2]);
                return *this;
            }

            inline void SetValue(DataType _a, DataType _b, DataType _c)
            {
                a() = _a;
                b() = _b;
                c() = _c;
            }

            HomoLine2 & Normalize()
            {
                if (0 == c())
                    return *this;

                DataType c_inv = 1.0 / c();
                a() = a() * c_inv;
                b() = b() * c_inv;
                c() = 1.0;
                return *this;
            }

            EigenVector Normalization() const
            {
                if (0 == c())
                    return *this;

                DataType c_inv = 1.0 / c();
                return *this * c_inv;
            }

            static HomoLine2 Infinity()
            {
                return HomoLine2(0, 0, 1);
            }

            inline bool IsInfinity()
            {
                return c() != 0 && a() == 0 && b() == 0;
            }

        public:
            DataType & a() { return this->data()[0]; }
            DataType & b() { return this->data()[1]; }
            DataType & c() { return this->data()[2]; }
            DataType const & a() const { return this->data()[0]; }
            DataType const & b() const { return this->data()[1]; }
            DataType const & c() const { return this->data()[2]; }
    };

}
}


#endif
