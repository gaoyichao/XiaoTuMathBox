#ifndef XTMB_HOMOPOINT2D_H
#define XTMB_HOMOPOINT2D_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Eigen>

namespace xiaotu {
namespace math {

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

            HomoPoint2D(double _x, double _y)
                : x(Eigen::Vector3d::data()[0]),
                  y(Eigen::Vector3d::data()[1]),
                  k(Eigen::Vector3d::data()[2])
            {
                SetValue(_x, _y, 1.0);
            }

            HomoPoint2D(double _x, double _y, double _k)
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

                x = x / k;
                y = y / k;
                k = 1.0;
                return *this;
            }

            static HomoPoint2D Infinity()
            {
                return HomoPoint2D(1, 1, 0);
            }

            inline Eigen::Vector3d JoinLine(HomoPoint2D const p)
            {
                return this->cross(p);
            }

            inline bool IsInfinity()
            {
                return 0 == k;
            }
        
        public:
            double & x;
            double & y;
            double & k;
    };
}
}

#endif
