/********************************************************************************************************
 * 
 * 在射影空间\(\mathbb{P}^2\)中的圆锥曲线
 * 
 * https://gaoyichao.com/Xiaotu/?book=几何&title=2D射影空间中的点直线和圆锥曲线
 *
 **************************************************************************** GAO YiChao 2022.0803 *****/
#ifndef XTMB_HOMOCONIC2_H
#define XTMB_HOMOCONIC2_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Eigen>

#include <XiaoTuMathBox/HomoPoint2.hpp>
#include <XiaoTuMathBox/HomoLine2.hpp>

namespace xiaotu {
namespace math {

    //! 圆锥曲线, 3X3对称矩阵
    template <typename DataType>
    class HomoConic2 : public Eigen::Matrix<DataType, 3, 3>
    {
        typedef Eigen::Matrix<DataType, 3, 3> Matrix3;
        public:
            HomoConic2()
            {
                SetValue(0, 0, 0,
                         0, 0, 1);
            }

            HomoConic2(DataType _a, DataType _b, DataType _c,
                       DataType _d, DataType _e, DataType _f)
            {
                SetValue(_a, _b, _c,
                         _d, _e, _f);
            }

            HomoConic2(Matrix3 const & m)
            {
                *this << m;
            }

            HomoConic2(HomoConic2 const & p)
            {
                *this << p;
            }

            HomoConic2 & operator = (HomoConic2 const & p)
            {
                *this << p;
                return *this;
            }

            HomoConic2 & operator = (Matrix3 const & m)
            {
                *this << m;
                return *this;
            }

            inline void SetValue(DataType _a, DataType _b, DataType _c,
                                 DataType _d, DataType _e, DataType _f)
            {
                DataType b_2 = 0.5 * _b;
                DataType d_2 = 0.5 * _d;
                DataType e_2 = 0.5 * _e;
                *this <<  _a, b_2, d_2,
                         b_2,  _c, e_2,
                         d_2, e_2,  _f;
            }

            inline void SetValue(DataType const * v)
            {
                SetValue(v[0], v[1], v[2], v[3], v[4], v[5]);
            }


            inline void From5Points(HomoPoint2<DataType> const * ps)
            {
                Eigen::Matrix<DataType, 5, 6> m;
                m << HomoConic2::PointEquation(ps[0]).transpose(),
                     HomoConic2::PointEquation(ps[1]).transpose(),
                     HomoConic2::PointEquation(ps[2]).transpose(),
                     HomoConic2::PointEquation(ps[3]).transpose(),
                     HomoConic2::PointEquation(ps[4]).transpose();

                Eigen::JacobiSVD<Eigen::Matrix<DataType, 5, 6>> svd(m, Eigen::ComputeFullV);
                SetValue(svd.matrixV().col(5).data());
            }

            HomoConic2 & Normalize()
            {
                DataType k = 0.0; 
                k = (std::abs(a()) > k) ? std::abs(a()) : k;
                k = (std::abs(b()) > k) ? std::abs(b()) : k;
                k = (std::abs(c()) > k) ? std::abs(c()) : k;
                k = (std::abs(d()) > k) ? std::abs(d()) : k;
                k = (std::abs(e()) > k) ? std::abs(e()) : k;
                k = (std::abs(f()) > k) ? std::abs(f()) : k;

                if (0 == k)
                    return *this;

                DataType k_inv = 1.0 / k;
                *this = *this * k_inv;
                return *this;
            }

            Matrix3 Normalization() const
            {
                DataType k = 0.0; 
                k = (std::abs(a()) > k) ? std::abs(a()) : k;
                k = (std::abs(b()) > k) ? std::abs(b()) : k;
                k = (std::abs(c()) > k) ? std::abs(c()) : k;
                k = (std::abs(d()) > k) ? std::abs(d()) : k;
                k = (std::abs(e()) > k) ? std::abs(e()) : k;
                k = (std::abs(f()) > k) ? std::abs(f()) : k;

                if (0 == k)
                    return *this;

                DataType k_inv = 1.0 / k;
                return *this * k_inv;
            }

            //! @brief 经过 p 的切线, 需由调用者保证点 p 在曲线上
            //!
            //! @param [in] p 圆锥曲线上一点
            inline Eigen::Matrix<DataType, 3, 1> Tangent(HomoPoint2<DataType> const & p) const
            {
                return (*this) * p;
            }
 
            //! @brief 按照圆锥曲线方程展开点，主要用于根据5个点构建圆锥曲线
            //!
            //! @param [in] p 目标点
            inline static Eigen::Matrix<DataType, 6, 1> PointEquation(HomoPoint2<DataType> const & p)
            {
                assert(0 != p.k());
                DataType kk_inv = p.k() * p.k();
                kk_inv = 1.0 / kk_inv;

                Eigen::Matrix<DataType, 6, 1> re;
                re << p.x() * p.x(), p.x() * p.y(), p.y() * p.y(),
                      p.x() * p.k(), p.y() * p.k(), p.k() * p.k();

                return re * kk_inv;
            }

        public:
            inline DataType a() const { return this->data()[0]; }
            inline DataType b() const { return this->data()[1] + this->data()[3]; }
            inline DataType c() const { return this->data()[4]; }
            inline DataType d() const { return this->data()[2] + this->data()[6]; }
            inline DataType e() const { return this->data()[5] + this->data()[7]; }
            inline DataType f() const { return this->data()[8]; }
    };
}
}

#endif
