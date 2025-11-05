#ifndef XTMB_GEO_GEOMETRY_H
#error "请勿直接引用 Vector3.hpp, 请使用 #include <XiaoTuMathBox/Geometry/Geometry.hpp>"
#endif

#include <cmath>
#include <iostream>
#include <vector>

namespace xiaotu {

    /**
     * @brief 三维欧式空间下的向量
     */
    template <typename DataType>
    class Vector3 : public AMatrix<DataType, 3, 1>
    {
        typedef AMatrix<DataType, 3, 1> Base;
        using Base::View;

        public:
            Vector3()
            {
                SetValue(0, 0, 0);
            }

            Vector3(DataType _x, DataType _y, DataType _z)
            {
                SetValue(_x, _y, _z);
            }

            template <typename M, bool MIsMatrix=M::IsMatrix>
            Vector3(M const & v)
            {
                SetValue(v(0), v(1), v(2));
            }

            Vector3(std::initializer_list<DataType> && li)
            {
                auto view = View();
                view = std::move(li);
            }

            template <typename M, bool MIsMatrix=M::IsMatrix>
            Vector3 & operator = (M const & v)
            {
                SetValue(v(0), v(1), v(2));
                return *this;
            }

            Vector3 & operator = (std::initializer_list<DataType> && li)
            {
                auto view = View();
                view = std::move(li);
                return *this;
            }

            inline void SetValue(DataType _x, DataType _y, DataType _z)
            {
                x() = _x;
                y() = _y;
                z() = _z;
            }

            Vector3 Cross(Vector3 const & v) const
            {
                Vector3 re;
                re << y() * v.z() - z() * v.y(),
                      z() * v.x() - x() * v.z(),
                      x() * v.y() - y() * v.x();
                return re;
            }

        public:
            inline DataType & x() { return this->At(0); }
            inline DataType & y() { return this->At(1); }
            inline DataType & z() { return this->At(2); }
            inline DataType const & x() const { return this->At(0); }
            inline DataType const & y() const { return this->At(1); }
            inline DataType const & z() const { return this->At(2); }
    };
}

