#ifndef XTMB_GEO_POINT2_H
#define XTMB_GEO_POINT2_H

#include <cmath>
#include <iostream>
#include <vector>

#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

namespace xiaotu::math {


    template <typename DataType>
    class Point2 : public AMatrix<DataType, 2, 1>
    {
        typedef AMatrix<DataType, 2, 1> Base;
        public:
            Point2()
            {
                SetValue(0.0, 0.0);
            }

            Point2(DataType _x, DataType _y)
            {
                SetValue(_x, _y);
            }

            Point2(Base const & v)
            {
                SetValue(v(0), v(1));
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

            Point2 & operator = (Base const & v)
            {
                SetValue(v(0), v(1));
                return *this;
            }

            inline void SetValue(DataType _x, DataType _y)
            {
                x() = _x;
                y() = _y;
            }

        public:
            inline DataType & x() { return this->At(0); }
            inline DataType & y() { return this->At(1); }
            inline DataType const & x() const { return this->At(0); }
            inline DataType const & y() const { return this->At(1); }
    };
}

#endif
