#ifndef XTMB_HOMOLINE2D_H
#define XTMB_HOMOLINE2D_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Eigen>

#include <XiaoTuMathBox/HomoPoint2D.h>

namespace xiaotu {
namespace math {

    class HomoLine2D : public Eigen::Vector3d
    {
        public:
            HomoLine2D()
                : a(Eigen::Vector3d::data()[0]),
                  b(Eigen::Vector3d::data()[1]),
                  c(Eigen::Vector3d::data()[2])
            {
                SetValue(0.0, 0.0, 1.0);
            }

            HomoLine2D(double _a, double _b, double _c)
                : a(Eigen::Vector3d::data()[0]),
                  b(Eigen::Vector3d::data()[1]),
                  c(Eigen::Vector3d::data()[2])
            {
                SetValue(_a, _b, _c);
            }

            HomoLine2D(Eigen::Vector3d const & v)
                : a(Eigen::Vector3d::data()[0]),
                  b(Eigen::Vector3d::data()[1]),
                  c(Eigen::Vector3d::data()[2])
            {
                SetValue(v[0], v[1], v[2]);
            }

            HomoLine2D(HomoLine2D const & p)
                : a(Eigen::Vector3d::data()[0]),
                  b(Eigen::Vector3d::data()[1]),
                  c(Eigen::Vector3d::data()[2])
            {
                SetValue(p.a, p.b, p.c);
            }

            HomoLine2D & operator = (HomoLine2D const & p)
            {
                SetValue(p.a, p.b, p.c);
                return *this;
            }

            HomoLine2D & operator = (Eigen::Vector3d const & v)
            {
                SetValue(v[0], v[1], v[2]);
                return *this;
            }

            friend bool operator == (HomoLine2D const & l1, HomoLine2D const & l2);

            inline void SetValue(double _a, double _b, double _c)
            {
                a = _a;
                b = _b;
                c = _c;
            }

            HomoLine2D & Normalize()
            {
                if (0 == c)
                    return *this;

                a = a / c;
                b = b / c;
                c = 1.0;

                return *this;
            }

            static HomoLine2D Infinity()
            {
                return HomoLine2D(0, 0, 1);
            }

            inline Eigen::Vector3d Intersection(HomoLine2D const & l)
            {
                return this->cross(l);
            }

            inline bool IsInfinity()
            {
                return c == 1 && a == 0 && b == 0;
            }
        
        public:
            double & a;
            double & b;
            double & c;
    };
}
}

#endif
