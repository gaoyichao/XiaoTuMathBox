/********************************************************************************************************
 * 
 * 在射影空间\(\mathbb{P}^2\)中，点和线对偶(dual) 
 * 在射影空间\(\mathbb{P}^3\)中，点和面对偶(dual) 
 *
 **************************************************************************** GAO YiChao 2021.1029 *****/
#ifndef XTMB_HOMOPOINT3D_H
#define XTMB_HOMOPOINT3D_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Eigen>

namespace xiaotu {
namespace math {
    class HomoPlane3D;

    class HomoPoint3D : public Eigen::Vector4d
    {
        public:
            HomoPoint3D()
                : x(Eigen::Vector4d::data()[0]),
                  y(Eigen::Vector4d::data()[1]),
                  z(Eigen::Vector4d::data()[2]),
                  k(Eigen::Vector4d::data()[3])
            {
                SetValue(0.0, 0.0, 0.0, 1.0);
            }

            HomoPoint3D(double _x, double _y, double _z, double _k = 1.0)
                : x(Eigen::Vector4d::data()[0]),
                  y(Eigen::Vector4d::data()[1]),
                  z(Eigen::Vector4d::data()[2]),
                  k(Eigen::Vector4d::data()[3])
            {
                SetValue(_x, _y, _z, _k);
            }

            HomoPoint3D(Eigen::Vector4d const & v)
                : x(Eigen::Vector4d::data()[0]),
                  y(Eigen::Vector4d::data()[1]),
                  z(Eigen::Vector4d::data()[2]),
                  k(Eigen::Vector4d::data()[3])
            {
                SetValue(v[0], v[1], v[2], v[3]);
            }

            HomoPoint3D(HomoPoint3D const & p)
                : x(Eigen::Vector4d::data()[0]),
                  y(Eigen::Vector4d::data()[1]),
                  z(Eigen::Vector4d::data()[2]),
                  k(Eigen::Vector4d::data()[3])
            {
                SetValue(p.x, p.y, p.z, p.k);
            }

            /*
             * 构造函数 - 三个平面的交点
             */
            HomoPoint3D(HomoPlane3D const & pi1, HomoPlane3D const & pi2, HomoPlane3D const & pi3);

            HomoPoint3D & operator = (HomoPoint3D const & p)
            {
                SetValue(p.x, p.y, p.z, p.k);
                return *this;
            }

            HomoPoint3D & operator = (Eigen::Vector4d const & v)
            {
                SetValue(v[0], v[1], v[2], v[3]);
                return *this;
            }

            friend bool operator == (HomoPoint3D const & p1, HomoPoint3D const & p2);

            inline void SetValue(double _x, double _y, double _z, double _k)
            {
                x = _x;
                y = _y;
                z = _z;
                k = _k;
            }

            HomoPoint3D & Normalize()
            {
                if (0 == k)
                    return *this;

                double k_inv = 1.0 / k;
                x = x * k_inv;
                y = y * k_inv;
                z = z * k_inv;
                k = 1.0;
                return *this;
            }

            Eigen::Vector4d Normalization() const
            {
                if (0 == k)
                    return *this;

                double k_inv = 1.0 / k;
                return *this * k_inv;
            }

            static const HomoPoint3D Infinity;
            static const HomoPoint3D XInfinity;
            static const HomoPoint3D YInfinity;
            static const HomoPoint3D ZInfinity;

            inline bool IsInfinity()
            {
                return 0 == k;
            }

            //inline Eigen::Vector3d JoinLine(HomoPoint3D const p)
            //{
            //    return this->cross(p);
            //}
        
        public:
            double & x;
            double & y;
            double & z;
            double & k;
    };
}
}

#endif

