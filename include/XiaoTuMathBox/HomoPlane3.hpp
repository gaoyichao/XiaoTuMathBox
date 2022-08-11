/********************************************************************************************************
 * 
 * 在射影空间\(\mathbb{P}^3\)中，点和面对偶(dual) 
 * 
 * https://gaoyichao.com/Xiaotu/?book=几何&title=index
 *
 **************************************************************************** GAO YiChao 2022.0810 *****/
#ifndef XTMB_HOMOUTILS3_H
#error "请勿直接引用 HomoPlane3.hpp, 请使用 #include <XiaoTuMathBox/HomoUtils3.hpp>"
#endif


#ifndef XTMB_HOMOPLANE3_H
#define XTMB_HOMOPLANE3_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Eigen>

namespace xiaotu {
namespace math {

    template <typename DataType>
    class HomoPlane3 : public Eigen::Matrix<DataType, 4, 1>
    {
        typedef Eigen::Matrix<DataType, 4, 1> EigenVector;
        public:
            HomoPlane3()
            {
                SetValue(0.0, 0.0, 0.0, 1.0);
            }

            HomoPlane3(DataType _a, DataType _b, DataType _c, DataType _d = 1.0)
            {
                SetValue(_a, _b, _c, _d);
            }

            HomoPlane3(EigenVector const & v)
            {
                SetValue(v[0], v[1], v[2], v[3]);
            }

            HomoPlane3(HomoPlane3 const & p)
            {
                SetValue(p.a(), p.b(), p.c(), p.d());
            }

            HomoPlane3 & operator = (HomoPlane3 const & p)
            {
                SetValue(p.a(), p.b(), p.c(), p.d());
                return *this;
            }

            HomoPlane3 & operator = (EigenVector const & v)
            {
                SetValue(v[0], v[1], v[2], v[3]);
                return *this;
            }

            inline void SetValue(DataType _a, DataType _b, DataType _c, DataType _d)
            {
                a() = _a;
                b() = _b;
                c() = _c;
                d() = _d;
            }

            HomoPlane3 & Normalize()
            {
                DataType k = 0.0; 
                k = (std::abs(a()) > k) ? std::abs(a()) : k;
                k = (std::abs(b()) > k) ? std::abs(b()) : k;
                k = (std::abs(c()) > k) ? std::abs(c()) : k;
                k = (std::abs(d()) > k) ? std::abs(d()) : k;
                if (0 == k)
                    return *this;

                DataType k_inv = 1.0 / k;
                a() = a() * k_inv;
                b() = b() * k_inv;
                c() = c() * k_inv;
                d() = d() * k_inv;
                return *this;
            }

            EigenVector Normalization() const
            {
                DataType k = 0.0; 
                k = (std::abs(a()) > k) ? std::abs(a()) : k;
                k = (std::abs(b()) > k) ? std::abs(b()) : k;
                k = (std::abs(c()) > k) ? std::abs(c()) : k;
                k = (std::abs(d()) > k) ? std::abs(d()) : k;
                if (0 == k)
                    return *this;

                DataType k_inv = 1.0 / k;
                return *this * k_inv;
            }

            inline bool IsInfinity()
            {
                return d() != 0 && a() == 0 && b() == 0 && c() == 0;
            }

            static HomoPlane3 Infinity()
            {
                return HomoPlane3(0, 0, 0, 1);
            }

        public:
            DataType & a() { return this->data()[0]; }
            DataType & b() { return this->data()[1]; }
            DataType & c() { return this->data()[2]; }
            DataType & d() { return this->data()[3]; }
            DataType const & a() const { return this->data()[0]; }
            DataType const & b() const { return this->data()[1]; }
            DataType const & c() const { return this->data()[2]; }
            DataType const & d() const { return this->data()[3]; }
    };

}
}

#endif
