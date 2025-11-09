#ifndef XTMB_GEO_GEOMETRY_H
#error "请勿直接引用 Line3.hpp, 请使用 #include <XiaoTuMathBox/Geometry/Geometry.hpp>"
#endif

#include <cmath>
#include <iostream>
#include <vector>

namespace xiaotu {


    /**
     * @brief 三维欧式空间下的直线
     */
    template <typename DataType>
    class Line3
    {
        public:
 
            /**
             * @brief 直线的普吕克坐标
             * 
             * @param [in] mom 直线到原点的矩
             * @param [in] dir 直线的方向向量
             */
            Line3(Vector3<DataType> const & mom, Vector3<DataType> const & dir)
            {
                mom_ = mom;
                dir_ = dir;
                Normalize();
            }

            inline bool IsInfinity()
            {
                return dir_.IsZero() && (!mom_.IsZero());
            }

            /**
             * @brief 对于一条直线，它的 6 个参数不能全为 0
             */
            bool IsValid()
            {
                return (!dir_.IsZero()) || (!mom_.IsZero());
            }

            /**
             * @brief 归一化，不修改对象本身
             */
            Line3 Normalization() const
            {
                return Line3(mom_, dir_);
            }
            /**
             * @brief 归一化，修改对象本身
             */
            Line3 & Normalize()
            {
                if (dir_.IsZero()) {
                    if (!mom_.IsZero())
                        mom_.Normalize();
                } else {
                    DataType s = 1.0 / dir_.Norm();
                    mom_ = s * mom_;
                    dir_ = s * dir_;
                }
                return *this;
            }

            friend std::ostream & operator << (std::ostream & os, Line3 const & line)
            {
                os << "(m:" << line.mom_(0) << ", " << line.mom_(1) << ", " << line.mom_(2)
                   << " |l: " << line.dir_(0) << ", " << line.dir_(1) << ", " << line.dir_(2)
                   << ")";
                return os;
            }

            Vector3<DataType> & Moment() { return mom_; }
            Vector3<DataType> const & Moment() const { return mom_; }

            Vector3<DataType> & Direction() { return dir_; }
            Vector3<DataType> const & Direction() const { return dir_; }

        private:
            //! @brief 矩，\boldsymbol{p} \times \boldsymbol{l}
            Vector3<DataType> mom_;
            //! @brief 直线的方向 \boldsymbol{l}
            Vector3<DataType> dir_;
    };

}

