/********************************************************************************************************
 * 
 * Plucker 矩阵形式的线, 两点共线
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
#ifndef XTMB_HOMOUTILS3_H
#error "请勿直接引用 HomoLine3.hpp, 请使用 #include <XiaoTuMathBox/HomoUtils3.hpp>"
#endif

#ifndef XTMB_HOMOLINE3_H
#define XTMB_HOMOLINE3_H

namespace xiaotu {
namespace math {
 
    template <typename DataType>
    class HomoLine3 : public Eigen::Matrix<DataType, 4, 4>
    {
        typedef Eigen::Matrix<DataType, 4, 4> Matrix4;
        typedef Eigen::Matrix<DataType, 4, 1> Vector4;
        public:
            HomoLine3()
            {
                this->setZero();
            }

            HomoLine3(DataType _l12, DataType _l13, DataType _l14,
                      DataType _l23, DataType _l42, DataType _l34)
            {
                SetValue(_l12, _l13, _l14, _l23, _l42, _l34);
            }

            HomoLine3(Matrix4 const & m)
            {
                *this << m;
            }

            HomoLine3(HomoLine3 const & l)
            {
                *this << l;
            }

            HomoLine3 & operator = (HomoLine3 const & l)
            {
                *this << l;
                return *this;
            }

            HomoLine3 & operator = (Matrix4 const & m)
            {
                *this << m;
                return *this;
            }

            inline HomoLine3 & SetValue(DataType _l12, DataType _l13, DataType _l14,
                                 DataType _l23, DataType _l42, DataType _l34)
            {
                SetL12(_l12); SetL13(_l13); SetL14(_l14);
                SetL23(_l23); SetL42(_l42); SetL34(_l34);
                return *this;
            }

            inline bool IsValid(DataType tolerance = 1e-9) const
            {
                DataType tmp = l12() * l34() + l13() * l42() + l14() * l23();
                return (std::abs(tmp) < tolerance);
            }

            inline bool AreCoplanar(HomoLine3 const & l, DataType tolerance = 1e-9) const
            {
                double tmp = l12() * l.l34() + l.l12() * l34() + l13() * l.l42()
                           + l.l13() * l42() + l14() * l.l23() + l.l14() * l23();

                return (std::abs(tmp) < tolerance);
            }


            inline HomoLine3 DualForm() const
            {
                HomoLine3 dual;
                dual <<    0,  l34(),  l42(),  l23(),
                        -l34(),    0,  l14(), -l13(),
                        -l42(), -l14(),    0,  l12(),
                        -l23(),  l13(), -l12(),    0;
                return dual;
            }

            inline Vector4 Direction() const
            {
                Vector4 re((*this)(0, 3), (*this)(1, 3), (*this)(2, 3), (*this)(3, 3));
                return re;
            }

        public:
            HomoLine3 & SetL12(DataType const & d)
            {
                this->data()[1 * 4 + 0] = d;
                this->data()[0 * 4 + 1] = -d;
                return *this;
            }

            HomoLine3 & SetL13(DataType const & d)
            {
                this->data()[2 * 4 + 0] = d;
                this->data()[0 * 4 + 2] = -d;
                return *this;
            }

            HomoLine3 & SetL14(DataType const & d)
            {
                this->data()[3 * 4 + 0] = d;
                this->data()[0 * 4 + 3] = -d;
                return *this;
            }

            HomoLine3 & SetL23(DataType const & d)
            {
                this->data()[2 * 4 + 1] = d;
                this->data()[1 * 4 + 2] = -d;
                return *this;
            }

            HomoLine3 & SetL42(DataType const & d)
            {
                this->data()[1 * 4 + 3] = d;
                this->data()[3 * 4 + 1] = -d;
                return *this;
            }

            HomoLine3 & SetL34(DataType const & d)
            {
                this->data()[3 * 4 + 2] = d;
                this->data()[2 * 4 + 3] = -d;
                return *this;
            }

        public:
            DataType l12() const { return this->data()[1 * 4 + 0]; }
            DataType l13() const { return this->data()[2 * 4 + 0]; }
            DataType l14() const { return this->data()[3 * 4 + 0]; }
            DataType l23() const { return this->data()[2 * 4 + 1]; }
            DataType l42() const { return this->data()[1 * 4 + 3]; }
            DataType l34() const { return this->data()[3 * 4 + 2]; }
    };
}
}

#endif
