#ifndef XTMB_HOMOCONIC2D_H
#define XTMB_HOMOCONIC2D_H

#include <cmath>
#include <iostream>
#include <vector>

#include <XiaoTuMathBox/EigenExt.h>

namespace xiaotu {
namespace math {

    /*
     * ax^2 + bxy + cy^2 + dx + ey + f = 0
     */
    class HomoConic2D : public Eigen::Matrix3d
    {
        public:
            HomoConic2D() {}

            HomoConic2D(double _a, double _b, double _c,
                        double _d, double _e, double _f)
            {
                SetValue(_a, _b, _c,
                         _d, _e, _f);
            }

            HomoConic2D(double const * v)
            {
                SetValue(v);
            }

            inline void SetValue(double const * v)
            {
                SetValue(v[0], v[1], v[2],
                         v[3], v[4], v[5]);
            }

            inline void SetValue(double _a, double _b, double _c,
                                 double _d, double _e, double _f)
            {
                SetA(_a); SetB(_b); SetC(_c);
                SetD(_d); SetE(_e); SetF(_f);
            }

            inline void From5Points(HomoPoint2D const * ps)
            {
                Eigen::Matrix5x6d m;
                m << HomoConic2D::Expand(ps[0]).transpose(),
                     HomoConic2D::Expand(ps[1]).transpose(),
                     HomoConic2D::Expand(ps[2]).transpose(),
                     HomoConic2D::Expand(ps[3]).transpose(),
                     HomoConic2D::Expand(ps[4]).transpose();

                Eigen::JacobiSVD<Eigen::Matrix5x6d> svd(m, Eigen::ComputeFullV);
                SetValue(svd.matrixV().col(5).data());
            }

            inline static Eigen::Vector6d Expand(HomoPoint2D const p)
            {
                Eigen::Vector6d re;
                re << p.x * p.x,
                      p.x * p.y,
                      p.y * p.y,
                      p.x * p.k,
                      p.y * p.k,
                      p.k * p.k;
                return re;
            }

            /*
             * Tangent -  过点 @p 的切线
             */
            inline Eigen::Vector3d Tangent(HomoPoint2D const & p) const
            {
                return (*this) * p;
            }

            inline static Eigen::Matrix3d Transform(HomoConic2D const & c, Eigen::Matrix3d const & H)
            {
                Eigen::Matrix3d H_inv = H.inverse();
                return H_inv.transpose() * c * H_inv;
            }

            inline HomoConic2D & Normalize()
            {
                if (0 == f())
                    return *this;

                double f_inv = 1.0 / f();
                SetA(a() * f_inv); SetB(b() * f_inv); SetC(c() * f_inv);
                SetD(d() * f_inv); SetE(e() * f_inv); SetF(1.0);

                return *this;
            }

            inline HomoConic2D & Transform(Eigen::Matrix3d const & H)
            {
                *(Eigen::Matrix3d*)this = Transform(*this, H);
                return *this;
            }

        public:
            inline double a() const { return this->data()[0]; }
            inline void SetA(double _a) { this->data()[0] = _a; }

            inline double b() const { return this->data()[1] + this->data()[3]; }
            inline void SetB(double _b)
            {
                double half_b = 0.5 * _b;
                this->data()[1] = half_b;
                this->data()[3] = half_b;
            }

            inline double c() const { return this->data()[4]; }
            inline void SetC(double _c) { this->data()[4] = _c; }

            inline double d() const { return this->data()[2] + this->data()[6]; }
            inline void SetD(double _d)
            {
                double half_d = 0.5 * _d;
                this->data()[2] = half_d;
                this->data()[6] = half_d;
            }

            inline double e() const { return this->data()[5] + this->data()[7]; }
            inline void SetE(double _e)
            {
                double half_e = 0.5 * _e;
                this->data()[5] = half_e;
                this->data()[7] = half_e;
            }

            inline double f() const { return this->data()[8]; }
            inline void SetF(double _f) { this->data()[8] = _f; }
    };

}
}

#endif
