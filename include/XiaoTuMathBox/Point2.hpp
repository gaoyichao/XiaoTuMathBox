#ifndef XTMB_POINT2_H
#define XTMB_POINT2_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Eigen>

namespace xiaotu {
namespace math {


    template <typename DataType>
    class Point2 : public Eigen::Matrix<DataType, 2, 1>
    {
        typedef Eigen::Matrix<DataType, 2, 1> EigenVector;
        public:
            Point2()
            {
                SetValue(0.0, 0.0);
            }

            Point2(DataType _x, DataType _y)
            {
                SetValue(_x, _y);
            }

            Point2(EigenVector const & v)
            {
                SetValue(v[0], v[1]);
            }

            Point2(Point2 const & p)
            {
                SetValue(p.x(), p.y());
            }

            Point2 & operator = (Point2 const & p)
            {
                SetValue(p.x(), p.y());
                return *this;
            }

            Point2 & operator = (EigenVector const & v)
            {
                SetValue(v[0], v[1]);
                return *this;
            }

            inline void SetValue(DataType _x, DataType _y)
            {
                x() = _x;
                y() = _y;
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
