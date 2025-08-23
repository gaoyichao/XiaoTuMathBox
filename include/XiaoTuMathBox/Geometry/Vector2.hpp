#ifndef XTMB_GEO_GEOMETRY_H
#error "请勿直接引用 VectorPoint2.hpp, 请使用 #include <XiaoTuMathBox/Geometry/Geometry.hpp>"
#endif


#include <cmath>
#include <iostream>
#include <vector>

namespace xiaotu::math {


    /**
     * @brief 二维欧式空间下的向量
     */
    template <typename DataType>
    class Vector2 : public AMatrix<DataType, 2, 1>
    {
        typedef AMatrix<DataType, 2, 1> Base;
        using Base::View;
        using MatrixBase = typename AMatrix<DataType, 2, 1>::Base;

        public:
            Vector2()
            {
                SetValue(0.0, 0.0);
            }

            Vector2(DataType _x, DataType _y)
            {
                SetValue(_x, _y);
            }

            template <typename M, bool MIsMatrix=M::IsMatrix>
            Vector2(M const & v)
            {
                SetValue(v(0), v(1));
            }

            Vector2(std::initializer_list<DataType> && li)
            {
                auto view = View();
                view = std::move(li);
            }

            template <typename M, bool MIsMatrix=M::IsMatrix>
            Vector2 & operator = (M const & v)
            {
                SetValue(v(0), v(1));
                return *this;
            }

            Vector2 & operator = (std::initializer_list<DataType> && li)
            {
                auto view = View();
                view = std::move(li);
                return *this;
            }

            inline void SetValue(DataType _x, DataType _y)
            {
                x() = _x;
                y() = _y;
            }

            inline DataType Cross(Vector2 const & v) const
            {
                return x() * v.y() - v.x() * y();
            }


        public:
            inline DataType & x() { return this->At(0); }
            inline DataType & y() { return this->At(1); }
            inline DataType const & x() const { return this->At(0); }
            inline DataType const & y() const { return this->At(1); }
    };
}

