#ifndef XTMB_HOMOUTILS3_H
#error "请勿直接引用 Projective3.hpp, 请使用 #include <XiaoTuMathBox/HomoUtils3.hpp>"
#endif

#ifndef XTMB_PERSPECTIVE3_H
#define XTMB_PERSPECTIVE3_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Eigen>

namespace xiaotu {
namespace math {

    //! 射影变换, 4X4矩阵 
    template <typename DataType>
    class Projective3 : public Eigen::Matrix<DataType, 4, 4>
    {
        typedef Eigen::Matrix<DataType, 4, 4> Matrix4;
        public:
            Projective3()
            {
                *this << 1, 0, 0, 0,
                         0, 1, 0, 0,
                         0, 0, 1, 0,
                         0, 0, 0, 1;
            }

            Projective3(Matrix4 const & m)
            {
                *this << m;
            }

            Projective3(Projective3 const & p)
            {
                *this << p;
            }

            Projective3 & operator = (Projective3 const & p)
            {
                *this << p;
                return *this;
            }

            Projective3 & operator = (Matrix4 const & m)
            {
                *this << m;
                return *this;
            }

            //! @brief 对点做射影变换
            inline Eigen::Matrix<DataType, 4, 1> ApplyOn(HomoPoint3<DataType> const & p)
            {
                return *this * p;
            }

            //! @brief 对平面做射影变换
            inline Eigen::Matrix<DataType, 4, 1> ApplyOn(HomoPlane3<DataType> const & p)
            {
                Eigen::Matrix<DataType, 4, 4> H_inv = this->inverse();
                return H_inv.transpose() * p;
            }

            //! @brief 对直线做射影变换
            inline Eigen::Matrix<DataType, 4, 4> ApplyOn(HomoLine3<DataType> const & l)
            {
                Eigen::Matrix<DataType, 4, 4> H_inv = this->inverse();
                return (*this) * l * this->transpose();
            }



    };
}
}

#endif
