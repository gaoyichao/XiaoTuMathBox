/********************************************************************************************************
 * 
 * 在射影空间\(\mathbb{P}^2\)中，点和线对偶(dual) 
 * 
 * https://gaoyichao.com/Xiaotu/?book=几何&title=2D射影空间中的点直线和圆锥曲线
 *
 **************************************************************************** GAO YiChao 2022.0803 *****/
#ifndef XTMB_GEO_GEOMETRY_H
#error "请勿直接引用 HomoLine2.hpp, 请使用 #include <XiaoTuMathBox/Geometry/Geometry.hpp>"
#endif

#include <cmath>
#include <iostream>
#include <vector>

namespace xiaotu {


    /**
     * @brief 二维摄影空间下的直线
     */
    template <typename DataType>
    class HomoLine2 : public AMatrix<DataType, 3, 1>
    {
        typedef AMatrix<DataType, 3, 1> AVector;
        using AVector::View;

        public:
            HomoLine2()
            {
                SetValue(0.0, 0.0, 1.0);
            }

            HomoLine2(DataType _a, DataType _b, DataType _c)
            {
                SetValue(_a, _b, _c);
            }

            template <typename Matrix, bool IsMatrix = Matrix::IsMatrix>
            HomoLine2(Matrix const & v)
            {
                SetValue(v(0), v(1), v(2));
            }

            HomoLine2(HomoLine2 const & p)
            {
                SetValue(p.a(), p.b(), p.c());
            }

            HomoLine2(std::initializer_list<DataType> && li)
            {
                auto view = View();
                view = std::move(li);
            }

            HomoLine2 & operator = (HomoLine2 const & p)
            {
                SetValue(p.a(), p.b(), p.c());
                return *this;
            }

            template <typename Matrix, bool IsMatrix = Matrix::IsMatrix>
            HomoLine2 & operator = (Matrix const & v)
            {
                SetValue(v(0), v(1), v(2));
                return *this;
            }

            HomoLine2 & operator = (std::initializer_list<DataType> && li)
            {
                auto view = View();
                view = std::move(li);
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
            //! @param [in] tolerance 分母为0的判定容忍度
            AVector Normalization(DataType tolerance = SMALL_VALUE) const
            {
                DataType norm = this->Norm();
                if (norm < tolerance) {
                    return {0, 0, 1};
                } else {
                    DataType norm_inv = 1.0 / norm;
                    return (*this) * norm_inv;
                }
            }

            //! @brief 归一化，修改对象本身
            //!
            //! @param [in] tolerance 分母为0的判定容忍度
            HomoLine2 & Normalize(DataType tolerance = SMALL_VALUE)
            {
                *this = Normalization(tolerance);
                return *this;
            }

            inline bool IsInfinity()
            {
                return c() != 0 && a() == 0 && b() == 0;
            }

            static HomoLine2 Infinity()
            {
                return HomoLine2(0, 0, 1);
            }

        public:
            DataType & a() { return this->At(0); }
            DataType & b() { return this->At(1); }
            DataType & c() { return this->At(2); }
            DataType const & a() const { return this->At(0); }
            DataType const & b() const { return this->At(1); }
            DataType const & c() const { return this->At(2); }
    };

}

