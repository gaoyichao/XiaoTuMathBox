/********************************************************************************************************
 * 
 * 几何意义上：
 * 射影映射 <i>projective</i> 是一个从 \(\mathbb{P}^2\) 到其自身空间的一种可逆映射，不妨用符号\(h\)表示，
 * 满足条件，若三点 \(\bold{x}_1, \bold{x}_2, \bold{x}_3\) 共线，
 * 则\(h(\bold{x}_1), h(\bold{x}_2), h(\bold{x}_3)\)也共线
 * 
 * 代数意义上：
 * 映射 \(h\) 是空间\(\mathbb{P}^2\)上的一个射影映射 <==> 存在一个非奇异的3x3的矩阵 H, 使得 \(h(x) = Hx\)
 * 
 * 射影映射也称保线映射(collineation), 构成一个 group
 * 1. 射影映射的逆也是射影映射
 * 2. 射影映射的组合还是射影映射
 * 
 *            H_e          H_s                            H_a                        H_p
 *          isometries --> similarity transformations --> affine transformations --> projective transformations
 *          等距           相似                           仿射                       射影
 * 不变量:  欧氏距离       形状                           平行关系                   共线关系
 *          3Dof           4Dof                           6Dof                       8Dof
 *          刚体                                                                     model vanishing points
 *                                                                                         灭点
 *                                                                                         
 * A projective transformation can be decomposed into a chain of transformations 
 * H = H_s H_a H_p
 * 
 **************************************************************************** GAO YiChao 2022.0805 *****/
#ifndef XTMB_HOMOUTILS2_H
#error "请勿直接引用 Projective2.hpp, 请使用 #include <XiaoTuMathBox/HomoUtils2.hpp>"
#endif

#ifndef XTMB_PERSPECTIVE2_H
#define XTMB_PERSPECTIVE2_H

namespace xiaotu {
namespace math {

    //! 射影变换, 3X3矩阵 
    template <typename DataType>
    class Projective2 : public Eigen::Matrix<DataType, 3, 3>
    {
        typedef Eigen::Matrix<DataType, 3, 3> Matrix3;
        public:
            Projective2()
            {
                *this << 1, 0, 0,
                         0, 1, 0,
                         0, 0, 1;
            }

            Projective2(Matrix3 const & m)
            {
                *this << m;
            }

            Projective2(Projective2 const & p)
            {
                *this << p;
            }

            Projective2 & operator = (Projective2 const & p)
            {
                *this << p;
                return *this;
            }

            Projective2 & operator = (Matrix3 const & m)
            {
                *this << m;
                return *this;
            }

            //! @brief 对点做射影变换
            inline Eigen::Matrix<DataType, 3, 1> ApplyOn(HomoPoint2<DataType> const & p)
            {
                return *this * p;
            }

            //! @brief 对线做射影变换
            inline Eigen::Matrix<DataType, 3, 1> ApplyOn(HomoLine2<DataType> const & l)
            {
                Eigen::Matrix<DataType, 3, 3> H_inv = this->inverse();
                return H_inv.transpose() * l;
            }

            //! @brief 对圆锥曲线做射影变换
            inline Eigen::Matrix<DataType, 3, 3> ApplyOn(HomoConic2<DataType> const & c)
            {
                Eigen::Matrix<DataType, 3, 3> H_inv = this->inverse();
                return H_inv.transpose() * c * H_inv;
            }
    };
}
}

#endif
