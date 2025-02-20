#ifndef XTMB_Vector2_H
#define XTMB_Vector2_H


#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Eigen>

namespace xiaotu {
namespace math {


    template <typename DataType>
    class Vector2 : public Eigen::Matrix<DataType, 2, 1>
    {
        typedef Eigen::Matrix<DataType, 2, 1> EigenVector;
        public:
            Vector2()
            {
                SetValue(0.0, 0.0);
            }

            Vector2(DataType _x, DataType _y)
            {
                SetValue(_x, _y);
            }

            Vector2(EigenVector const & v)
            {
                SetValue(v[0], v[1]);
            }

            Vector2(Vector2 const & p)
            {
                SetValue(p.x(), p.y());
            }

            Vector2 & operator = (Vector2 const & p)
            {
                SetValue(p.x(), p.y());
                return *this;
            }

            Vector2 & operator = (EigenVector const & v)
            {
                SetValue(v[0], v[1]);
                return *this;
            }

            inline void SetValue(DataType _x, DataType _y)
            {
                x() = _x;
                y() = _y;
            }

            inline DataType Norm() const
            {
                return std::sqrt(x() * x() + y() * y());
            }

            inline DataType Cross(Vector2 const & v) const
            {
                return x() * v.y() - v.x() * y();
            }

            inline EigenVector Normalization() const
            {
                DataType factor = 1.0 / Norm();
                return *this * factor;
            }

            inline Vector2 & Normalize()
            {
                *this << Normalization();
                return *this;
            }

        public:
            inline DataType & x() { return this->data()[0]; }
            inline DataType & y() { return this->data()[1]; }
            inline DataType const & x() const { return this->data()[0]; }
            inline DataType const & y() const { return this->data()[1]; }
    };
}
}

#endif
