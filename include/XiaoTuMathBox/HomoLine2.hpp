/********************************************************************************************************
 * 
 * 在射影空间\(\mathbb{P}^2\)中，点和线对偶(dual) 
 * 
 * https://gaoyichao.com/Xiaotu/?book=几何&title=2D射影空间中的点直线和圆锥曲线
 *
 **************************************************************************** GAO YiChao 2022.0803 *****/
#ifndef XTMB_HOMOUTILS2_H
#error "请勿直接引用 HomoLine2.hpp, 请使用 #include <XiaoTuMathBox/HomoUtils2.hpp>"
#endif

#ifndef XTMB_HOMOLINE2_H
#define XTMB_HOMOLINE2_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Eigen>

namespace xiaotu {
namespace math {

    template <typename DataType>
    class HomoPoint2;

    template <typename DataType>
    class HomoLine2 : public Eigen::Matrix<DataType, 3, 1>
    {
        typedef Eigen::Matrix<DataType, 3, 1> EigenVector;
        public:
            HomoLine2()
            {
                SetValue(0.0, 0.0, 1.0);
            }

            HomoLine2(DataType _a, DataType _b, DataType _c)
            {
                SetValue(_a, _b, _c);
            }

            HomoLine2(EigenVector const & v)
            {
                SetValue(v[0], v[1], v[2]);
            }

            HomoLine2(HomoLine2 const & p)
            {
                SetValue(p.a(), p.b(), p.c());
            }

            HomoLine2 & operator = (HomoLine2 const & p)
            {
                SetValue(p.a(), p.b(), p.c());
                return *this;
            }

            HomoLine2 & operator = (EigenVector const & v)
            {
                SetValue(v[0], v[1], v[2]);
                return *this;
            }

            inline void SetValue(DataType _a, DataType _b, DataType _c)
            {
                a() = _a;
                b() = _b;
                c() = _c;
            }

            //! @brief 归一化，不修改对象本身
            //!
            //! 直线方程: \(ax + by + c = 0\)
            //! 点到直线的距离公式:
            //! $$
            //!     dis = \frac{|ax + by + c|}{\sqrt{a^2 + b^2}}
            //! $$
            //!
            //! @param [in] tolerance 分母为0的判定容忍度
            EigenVector Normalization(DataType tolerance = 1e-9) const
            {
                DataType k = std::sqrt(a() * a() + b() * b());
                EigenVector re;
                if (k < tolerance) {
                    re << 0, 0, 1;
                } else {
                    DataType k_inv = 1.0 / k;
                    re = *this * k_inv;
                }
                return re;
            }
            //! @brief 归一化，修改对象本身
            //!
            //! @param [in] tolerance 分母为0的判定容忍度
            HomoLine2 & Normalize(DataType tolerance = 1e-9)
            {
                (*this) << Normalization(tolerance);
                return *this;
            }

            static HomoLine2 Infinity()
            {
                return HomoLine2(0, 0, 1);
            }

            inline bool IsInfinity()
            {
                return c() != 0 && a() == 0 && b() == 0;
            }

        public:
            DataType & a() { return this->data()[0]; }
            DataType & b() { return this->data()[1]; }
            DataType & c() { return this->data()[2]; }
            DataType const & a() const { return this->data()[0]; }
            DataType const & b() const { return this->data()[1]; }
            DataType const & c() const { return this->data()[2]; }
    };

}
}


#endif
