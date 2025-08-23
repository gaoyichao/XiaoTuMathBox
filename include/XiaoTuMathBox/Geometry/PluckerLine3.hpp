/********************************************************************************************************
 * 
 * Plucker 矩阵形式的线, 两点共线, 两面交线
 * 
 * 在 3-space 中有4个自由度
 * 
 * 在 MVG 一书中提到表示直线的三种方法:
 * 1. 两点或者两直线的生成子空间 W^T = \lambda A + \miu B
 * 
 *              | A^T |
 *    W_{2*4} = | B^T |
 *    
 *    W 的右零空间则是直线或者点的pencil
 *    
 * 2. Plücker 矩阵
 *   
 *    L = A * B^T - B * A^T
 *    
 *    该阵虽然是一个4*4的阵，但它的秩只有2，有4个自由度
 *    该阵与点或者面AB的选择无关
 *  
 * 3. Plücker line coordinates
 *  
 *    L = { l12, l13, l14, l23, l42, l34 }
 *
 * 
 * https://gaoyichao.com/Xiaotu/?book=几何&title=index
 *
 **************************************************************************** GAO YiChao 2022.0810 *****/
#ifndef XTMB_GEO_GEOMETRY_H
#error "请勿直接引用 PluckerLine3.hpp, 请使用 #include <XiaoTuMathBox/Geometry/Geometry.hpp>"
#endif


namespace xiaotu::math {
 
    /**
     * @brief 三维射影空间下的普吕克直线
     */
    template <typename DataType>
    class PluckerLine3 : public AMatrix<DataType, 4, 4>
    {
        typedef AMatrix<DataType, 4, 4> Mat4;
        typedef AMatrix<DataType, 4, 1> Vec4;

        public:
            PluckerLine3()
            {
                this->Zeroing();
            }

            PluckerLine3(DataType _l12, DataType _l13, DataType _l14,
                         DataType _l23, DataType _l42, DataType _l34)
            {
                SetValue(_l12, _l13, _l14, _l23, _l42, _l34);
            }

            template <typename Matrix, bool IsMatrix = Matrix::IsMatrix>
            PluckerLine3(Matrix const & m)
            {
                *this << m;
            }

            PluckerLine3(PluckerLine3 const & l)
            {
                *this << l;
            }

            PluckerLine3 & operator = (PluckerLine3 const & l)
            {
                *this << l;
                return *this;
            }

            template <typename Matrix, bool IsMatrix = Matrix::IsMatrix>
            PluckerLine3 & operator = (Matrix const & m)
            {
                *this << m;
                return *this;
            }

            inline PluckerLine3 & SetValue(DataType _l12, DataType _l13, DataType _l14,
                                           DataType _l23, DataType _l42, DataType _l34)
            {
                this->At(0, 0) = (DataType)0;
                this->At(1, 1) = (DataType)0;
                this->At(2, 2) = (DataType)0;
                this->At(3, 3) = (DataType)0;

                SetL12(_l12); SetL13(_l13); SetL14(_l14);
                SetL23(_l23); SetL42(_l42); SetL34(_l34);
                return *this;
            }

            inline bool IsValid(DataType tolerance = 1e-9) const
            {
                DataType tmp = l12() * l34() + l13() * l42() + l14() * l23();
                return (std::abs(tmp) < tolerance);
            }


            inline PluckerLine3 DualForm() const
            {
                PluckerLine3 dual;
                dual <<    0,  l34(),  l42(),  l23(),
                        -l34(),    0,  l14(), -l13(),
                        -l42(), -l14(),    0,  l12(),
                        -l23(),  l13(), -l12(),    0;
                return dual;
            }

            inline Vec4 Direction() const
            {
                return { (*this)(0, 3), (*this)(1, 3), (*this)(2, 3), (*this)(3, 3) };
            }

        public:
            PluckerLine3 & SetL12(DataType const & d)
            {
                this->At(0, 1) = d;
                this->At(1, 0) = -d;
                return *this;
            }

            PluckerLine3 & SetL13(DataType const & d)
            {
                this->At(0, 2) = d;
                this->At(2, 0) = -d;
                return *this;
            }

            PluckerLine3 & SetL14(DataType const & d)
            {
                this->At(0, 3) = d;
                this->At(3, 0) = -d;
                return *this;
            }

            PluckerLine3 & SetL23(DataType const & d)
            {
                this->At(1, 2) = d;
                this->At(2, 1) = -d;
                return *this;
            }

            PluckerLine3 & SetL42(DataType const & d)
            {
                this->At(3, 1) = d;
                this->At(1, 3) = -d;
                return *this;
            }

            PluckerLine3 & SetL34(DataType const & d)
            {
                this->At(2, 3) = d;
                this->At(3, 2) = -d;
                return *this;
            }

        public:
            DataType l12() const { return this->At(0, 1); }
            DataType l13() const { return this->At(0, 2); }
            DataType l14() const { return this->At(0, 3); }
            DataType l23() const { return this->At(1, 2); }
            DataType l42() const { return this->At(3, 1); }
            DataType l34() const { return this->At(2, 3); }
    };
}
