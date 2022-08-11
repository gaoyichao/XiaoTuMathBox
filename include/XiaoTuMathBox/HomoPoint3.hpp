/********************************************************************************************************
 * 
 * 在射影空间\(\mathbb{P}^3\)中，点和面对偶(dual) 
 * 
 * https://gaoyichao.com/Xiaotu/?book=几何&title=index
 *
 **************************************************************************** GAO YiChao 2022.0810 *****/
#ifndef XTMB_HOMOUTILS3_H
#error "请勿直接引用 HomoPoint3.hpp, 请使用 #include <XiaoTuMathBox/HomoUtils3.hpp>"
#endif

#ifndef XTMB_HOMOPOINT3_H
#define XTMB_HOMOPOINT3_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Eigen>

namespace xiaotu {
namespace math {

    template <typename DataType>
    class HomoPoint3 : public Eigen::Matrix<DataType, 4, 1>
    {
        typedef Eigen::Matrix<DataType, 4, 1> EigenVector;
        public:
            HomoPoint3()
            {
                SetValue(0.0, 0.0, 0.0, 1.0);
            }

            HomoPoint3(DataType _x, DataType _y, DataType _z, DataType _k = 1.0)
            {
                SetValue(_x, _y, _z, _k);
            }

            HomoPoint3(EigenVector const & v)
            {
                SetValue(v[0], v[1], v[2], v[3]);
            }

            HomoPoint3(HomoPoint3 const & p)
            {
                SetValue(p.x(), p.y(), p.z(), p.k());
            }

            HomoPoint3 & operator = (HomoPoint3 const & p)
            {
                SetValue(p.x(), p.y(), p.z(), p.k());
                return *this;
            }

            HomoPoint3 & operator = (EigenVector const & v)
            {
                SetValue(v[0], v[1], v[2], v[3]);
                return *this;
            }

            inline void SetValue(DataType _x, DataType _y, DataType _z, DataType _k)
            {
                x() = _x;
                y() = _y;
                z() = _z;
                k() = _k;
            }

            HomoPoint3 & Normalize()
            {
                if (0 == k())
                    return *this;

                DataType k_inv = 1.0 / k();
                x() = x() * k_inv;
                y() = y() * k_inv;
                z() = z() * k_inv;
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

            inline bool IsValid(DataType tolerance = 1e-9) const
            {
                return (std::abs(x()) >= tolerance) ||
                       (std::abs(y()) >= tolerance) ||
                       (std::abs(z()) >= tolerance) ||
                       (std::abs(k()) >= tolerance);
            }

            inline bool IsInfinity()
            {
                return 0 == k();
            }

            static HomoPoint3 Infinity()
            {
                return HomoPoint3(1, 0, 0, 0);
            }


        public:
            DataType & x() { return this->data()[0]; }
            DataType & y() { return this->data()[1]; }
            DataType & z() { return this->data()[2]; }
            DataType & k() { return this->data()[3]; }
            DataType const & x() const { return this->data()[0]; }
            DataType const & y() const { return this->data()[1]; }
            DataType const & z() const { return this->data()[2]; }
            DataType const & k() const { return this->data()[3]; }
    };
}
}

#endif
