/********************************************************************************************************
 * 
 * 在射影空间\(\mathbb{P}^2\)中，点和线对偶(dual) 
 * 在射影空间\(\mathbb{P}^3\)中，点和面对偶(dual) 
 *
 **************************************************************************** GAO YiChao 2021.1029 *****/
#ifndef XTMB_HOMOPOINT2D_H
#define XTMB_HOMOPOINT2D_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Eigen>

namespace xiaotu {
namespace math {

    class HomoLine2D;

    class HomoPoint2D : public Eigen::Vector3d
    {
        public:
            HomoPoint2D()
                : x(Eigen::Vector3d::data()[0]),
                  y(Eigen::Vector3d::data()[1]),
                  k(Eigen::Vector3d::data()[2])
            {
                SetValue(0.0, 0.0, 1.0);
            }

            HomoPoint2D(double _x, double _y, double _k = 1.0)
                : x(Eigen::Vector3d::data()[0]),
                  y(Eigen::Vector3d::data()[1]),
                  k(Eigen::Vector3d::data()[2])
            {
                SetValue(_x, _y, _k);
            }

            HomoPoint2D(Eigen::Vector3d const & v)
                : x(Eigen::Vector3d::data()[0]),
                  y(Eigen::Vector3d::data()[1]),
                  k(Eigen::Vector3d::data()[2])
            {
                SetValue(v[0], v[1], v[2]);
            }

            HomoPoint2D(HomoPoint2D const & p)
                : x(Eigen::Vector3d::data()[0]),
                  y(Eigen::Vector3d::data()[1]),
                  k(Eigen::Vector3d::data()[2])
            {
                SetValue(p.x, p.y, p.k);
            }

            HomoPoint2D(HomoLine2D const & l1, HomoLine2D const & l2);

            HomoPoint2D & operator = (HomoPoint2D const & p)
            {
                SetValue(p.x, p.y, p.k);
                return *this;
            }

            HomoPoint2D & operator = (Eigen::Vector3d const & v)
            {
                SetValue(v[0], v[1], v[2]);
                return *this;
            }

            friend bool operator == (HomoPoint2D const & p1, HomoPoint2D const & p2);

            inline void SetValue(double _x, double _y, double _k)
            {
                x = _x;
                y = _y;
                k = _k;
            }

            HomoPoint2D & Normalize()
            {
                if (0 == k)
                    return *this;

                double k_inv = 1.0 / k;
                x = x * k_inv;
                y = y * k_inv;
                k = 1.0;
                return *this;
            }

            Eigen::Vector3d Normalization() const
            {
                if (0 == k)
                    return *this;

                double k_inv = 1.0 / k;
                return *this * k_inv;
            }

            static HomoPoint2D Infinity()
            {
                return HomoPoint2D(0, 0, 0);
            }

            inline bool IsInfinity()
            {
                return 0 == k;
            }

            inline Eigen::Vector3d JoinLine(HomoPoint2D const p) const
            {
                return this->cross(p);
            }

        public:
            double & x;
            double & y;
            double & k;
    };
}
}

#endif
