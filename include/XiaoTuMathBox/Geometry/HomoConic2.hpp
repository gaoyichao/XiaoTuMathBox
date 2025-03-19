/********************************************************************************************************
 * 
 * 在射影空间\(\mathbb{P}^2\)中的圆锥曲线
 * 
 * https://gaoyichao.com/Xiaotu/?book=几何&title=2D射影空间中的点直线和圆锥曲线
 *
 **************************************************************************** GAO YiChao 2022.0803 *****/
// #ifndef XTMB_HOMOUTILS2_H
// #error "请勿直接引用 HomoConic2.hpp, 请使用 #include <XiaoTuMathBox/HomoUtils2.hpp>"
// #endif

#ifndef XTMB_HOMOCONIC2_H
#define XTMB_HOMOCONIC2_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>


namespace xiaotu::math {

    //! 圆锥曲线, 3x3对称矩阵
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


        public:
            inline DataType a() const { return this->At(0, 0); }
            inline DataType b() const { return this->At(0, 1) + this->At(1, 0); }
            inline DataType c() const { return this->At(1, 1); }
            inline DataType d() const { return this->At(0, 2) + this->At(2, 0); }
            inline DataType e() const { return this->At(1, 2) + this->At(2, 1); }
            inline DataType f() const { return this->At(2, 2); }
    };

}

#endif

