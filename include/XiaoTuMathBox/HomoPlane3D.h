/********************************************************************************************************
 * 
 * 在射影空间\(\mathbb{P}^2\)中，点和线对偶(dual) 
 * 在射影空间\(\mathbb{P}^3\)中，点和面对偶(dual) 
 *
 **************************************************************************** GAO YiChao 2021.1029 *****/
#ifndef XTMB_HOMOPLANE3D_H
#define XTMB_HOMOPLANE3D_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Eigen>

#include <XiaoTuMathBox/HomoPoint3D.h>

namespace xiaotu {
namespace math {
    class HomoPoint3D;
    class HomoPlane3D;

    Eigen::Vector4d Intersection(HomoPlane3D const & pi1, HomoPlane3D const & pi2, HomoPlane3D const & pi3);

    /*
     * ax + by + cz + d = 0
     */
    class HomoPlane3D : public Eigen::Vector4d
    {
        public:
            HomoPlane3D()
                : a(Eigen::Vector4d::data()[0]),
                  b(Eigen::Vector4d::data()[1]),
                  c(Eigen::Vector4d::data()[2]),
                  d(Eigen::Vector4d::data()[3])
            {
                SetValue(0.0, 0.0, 0.0, 1.0);
            }

            HomoPlane3D(double _a, double _b, double _c, double _d)
                : a(Eigen::Vector4d::data()[0]),
                  b(Eigen::Vector4d::data()[1]),
                  c(Eigen::Vector4d::data()[2]),
                  d(Eigen::Vector4d::data()[3])
            {
                SetValue(_a, _b, _c, _d);
            }

            HomoPlane3D(Eigen::Vector4d const & v)
                : a(Eigen::Vector4d::data()[0]),
                  b(Eigen::Vector4d::data()[1]),
                  c(Eigen::Vector4d::data()[2]),
                  d(Eigen::Vector4d::data()[3])
            {
                SetValue(v[0], v[1], v[2], v[3]);
            }

            HomoPlane3D(HomoPlane3D const & p)
                : a(Eigen::Vector4d::data()[0]),
                  b(Eigen::Vector4d::data()[1]),
                  c(Eigen::Vector4d::data()[2]),
                  d(Eigen::Vector4d::data()[3])
            {
                SetValue(p.a, p.b, p.c, p.d);
            }

            /*
             * 构造函数 - 三个不共线的点确定一个平面
             */
            HomoPlane3D(HomoPoint3D const & po1, HomoPoint3D const & po2, HomoPoint3D const & po3);

            HomoPlane3D & operator = (HomoPlane3D const & p)
            {
                SetValue(p.a, p.b, p.c, p.d);
                return *this;
            }

            HomoPlane3D & operator = (Eigen::Vector4d const & v)
            {
                SetValue(v[0], v[1], v[2], v[3]);
                return *this;
            }

            friend bool operator == (HomoPlane3D const & p1, HomoPlane3D const & p2);

            inline void SetValue(double _a, double _b, double _c, double _d)
            {
                a = _a;
                b = _b;
                c = _c;
                d = _d;
            }

            HomoPlane3D & Normalize()
            {
                if (0 == d)
                    return *this;

                double d_inv = 1.0 / d;
                a = a * d_inv;
                b = b * d_inv;
                c = c * d_inv;
                d = 1.0;

                return *this;
            }

            Eigen::Vector4d Normalization() const
            {
                if (0 == d)
                    return *this;

                double d_inv = 1.0 / d;
                return *this * d_inv;
            }

            static HomoPlane3D Infinity()
            {
                return HomoPlane3D(0, 0, 0, 1);
            }

            inline bool IsInfinity()
            {
                return d == 1 && a == 0 && b == 0 && c == 0;
            }

            friend Eigen::Vector4d Intersection(HomoPlane3D const & pi1, HomoPlane3D const & pi2, HomoPlane3D const & pi3);

        public:
            double & a;
            double & b;
            double & c;
            double & d;
    };
 
}
}

#endif
