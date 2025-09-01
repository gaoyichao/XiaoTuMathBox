#ifndef XTMB_GEO_GEOMETRY_H
#error "请勿直接引用 Plane3.hpp, 请使用 #include <XiaoTuMathBox/Geometry/Geometry.hpp>"
#endif

#include <cmath>
#include <iostream>
#include <vector>

namespace xiaotu::math {
    /**
     * @brief 三维欧式空间下的直线 ax + by + cz + d = 0
     */
    template <typename DataType>
    class Plane3
    {
        public:
            /**
             * @brief 构造函数, 一般方程
             * 
             * @param [in] norm 平面的法向量
             * @param [in] d 平面方程的常数项
             */
            Plane3(Vector3<DataType> const & norm, DataType d)
            {
                norm_ << norm;
                d_ = d;
                Normalize();
            }

            /**
             * @brief 构造函数, 点法式
             * 
             * @param [in] norm 平面的法向量
             * @param [in] point 平面上一点
             */
            Plane3(Vector3<DataType> const & norm, Point3<DataType> const & point)
            {
                norm_ << norm;
                d_ = -norm.Dot(point);
                Normalize();
            }

            /**
             * @brief 在欧式空间中，平面的法向量为 0 不具有几何意义
             */
            bool IsValid()
            {
                return !norm_.IsZero();
            }

            /**
             * @brief 归一化，不修改对象本身
             */
            Plane3 Normalization() const
            {
                return Plane3(norm_, d_);
            }
            /**
             * @brief 归一化，修改对象本身
             */
            Plane3 & Normalize()
            {
                if (norm_.IsZero())
                    return *this;
                DataType s = 1.0 / norm_.Norm();
                norm_ = s * norm_;
                d_ = s * d_;
                return *this;
            }

        public:
            friend std::ostream & operator << (std::ostream & os, Plane3 const & p)
            {
                os << "(" << p.a() << ", " << p.b() << ", " << p.c()
                   << " | " << p.d() << ")";
                return os;
            }

            DataType & a() { return norm_.At(0); }
            DataType & b() { return norm_.At(1); }
            DataType & c() { return norm_.At(2); }
            DataType & d() { return d_; }
            DataType const & a() const { return norm_.At(0); }
            DataType const & b() const { return norm_.At(1); }
            DataType const & c() const { return norm_.At(2); }
            DataType const & d() const { return d_; }

            /**
             * @brief 获取平面的法向量
             */
            Vector3<DataType> & Normal() { return norm_; }

            /**
             * @brief 获取平面的法向量
             */
            Vector3<DataType> const & Normal() const { return norm_; }
 
        private:
            //! @brief 平面的法向量(a, b, c)
            Vector3<DataType> norm_;
            //! @brief 平面方程的常数项 ax + by + cz + d = 0
            DataType d_;
    };

}

