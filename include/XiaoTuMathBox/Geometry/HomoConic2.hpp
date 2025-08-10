/********************************************************************************************************
 * 
 * 在射影空间\(\mathbb{P}^2\)中的圆锥曲线
 * 
 * https://gaoyichao.com/Xiaotu/?book=几何&title=2D射影空间中的点直线和圆锥曲线
 *
 **************************************************************************** GAO YiChao 2022.0803 *****/
#ifndef XTMB_GEO_HOMOUTILS2_H
#error "请勿直接引用 HomoConic2.hpp, 请使用 #include <XiaoTuMathBox/Geometry/HomoUtils2.hpp>"
#endif

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>


namespace xiaotu::math {

    /**
     * @brief 圆锥曲线, 3x3对称矩阵
     */
    template <typename DataType>
    class HomoConic2 : public AMatrix<DataType, 3, 3>
    {
        typedef AMatrix<DataType, 3, 3> AMatrix3;

        public:
            HomoConic2()
            {
                SetValue(0, 0, 0, 0, 0, 1);
            }

            HomoConic2(DataType _a, DataType _b, DataType _c,
                       DataType _d, DataType _e, DataType _f)
            {
                SetValue(_a, _b, _c, _d, _e, _f);
            }

            HomoConic2(HomoConic2 const & p)
            {
                this->Assign(p);
            }

            template <typename Mat, bool IsMatrix = Mat::IsMatrix>
            HomoConic2(Mat const & mv)
            {
                this->Assign(mv);
            }

            HomoConic2 & operator = (HomoConic2 const & p)
            {
                this->Assign(p);
                return *this;
            }

            template <typename Mat, bool IsMatrix = Mat::IsMatrix>
            HomoConic2 & operator = (Mat const & mv)
            {
                this->Assign(mv);
                return *this;
            }

            inline void SetValue(DataType _a, DataType _b, DataType _c,
                                 DataType _d, DataType _e, DataType _f)
            {
                DataType b_2 = 0.5 * _b;
                DataType d_2 = 0.5 * _d;
                DataType e_2 = 0.5 * _e;
                this->At(0, 0) =  _a; this->At(0, 1) = b_2; this->At(0, 2) = d_2;
                this->At(1, 0) = b_2; this->At(1, 1) =  _c; this->At(1, 2) = e_2;
                this->At(2, 0) = d_2; this->At(2, 1) = e_2; this->At(2, 2) =  _f;
            }

            /**
             * @brief 按照圆锥曲线方程展开点，主要用于根据5个点构建圆锥曲线
             * 
             *  @param [in] p 目标点
             */
            inline static AMatrix<DataType, 1, 6> PointEquation(HomoPoint2<DataType> const & p)
            {
                assert(0 != p.k());

                DataType kk_inv = p.k() * p.k();
                kk_inv = 1.0 / kk_inv;

                AMatrix<DataType, 1, 6> re;
                re << p.x() * p.x(), p.x() * p.y(), p.y() * p.y(),
                      p.x() * p.k(), p.y() * p.k(), p.k() * p.k();

                return re * kk_inv;
            }

            /**
             * @brief 经过 p 的切线, 需由调用者保证点 p 在曲线上
             * 
             *  @param [in] p 圆锥曲线上一点
             */
            inline AMatrix<DataType, 3, 1> Tangent(HomoPoint2<DataType> const & p) const
            {
                return (*this) * p;
            }

            /**
             * @brief 正则化
             */
            AMatrix<DataType, 3, 3> Normalization() const
            {
                DataType k = 0.0; 
                k = (std::abs(a()) > k) ? std::abs(a()) : k;
                k = (std::abs(b()) > k) ? std::abs(b()) : k;
                k = (std::abs(c()) > k) ? std::abs(c()) : k;
                k = (std::abs(d()) > k) ? std::abs(d()) : k;
                k = (std::abs(e()) > k) ? std::abs(e()) : k;
                k = (std::abs(f()) > k) ? std::abs(f()) : k;

                if (0 == k)
                    return *this;

                DataType k_inv = 1.0 / k;
                return *this * k_inv;
            }

            /**
             * @brief 正则化，修改对象本身
             */
            HomoConic2 & Normalize()
            {
                *this = Normalization();
                return *this;
            }

            /**
             * @brief 判定是否满足圆锥曲线方程的形式
             */
            bool Valid() const
            {
                return this->IsSymmetric();
            }

        public:
            inline DataType a() const { return this->At(0, 0); }
            inline DataType b() const { return this->At(0, 1) + this->At(1, 0); }
            inline DataType c() const { return this->At(1, 1); }
            inline DataType d() const { return this->At(0, 2) + this->At(2, 0); }
            inline DataType e() const { return this->At(1, 2) + this->At(2, 1); }
            inline DataType f() const { return this->At(2, 2); }
    };

}

