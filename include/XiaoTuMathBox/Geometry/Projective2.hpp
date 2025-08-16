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
#ifndef XTMB_GEO_GEOMETRY_H
#error "请勿直接引用 Projective2.hpp, 请使用 #include <XiaoTuMathBox/Geometry/Geometry.hpp>"
#endif

namespace xiaotu::math {

 
    /**
     * @brief 射影变换, 3X3矩阵
     * 
     * 出于教学目的编写，只是为了说明射影变换如何用到点、直线、圆锥曲线上的。
     * 工程中通常直接进行矩阵运算，基本不需要专门封装一个数据结构出来。
     */
    template <typename DataType>
    class Projective2 : public AMatrix<DataType, 3, 3>
    {
        public:
            typedef AMatrix<DataType, 3, 3> Matrix3;

        public:
            /**
             * @brief 默认构造
             */
            Projective2()
            {
                *this << 1, 0, 0,
                         0, 1, 0,
                         0, 0, 1;
            }

            /**
             * @brief 拷贝构造，深度拷贝
             */
            Projective2(Projective2 const & p)
            {
                *this << p;
            }

            /**
             * @brief 从 3x3 的矩阵中构造
             */
            template <typename Derived>
            Projective2(MatrixBase<Derived> const & m)
            {
                *this << m;
            }

            /**
             * @brief 拷贝复制
             */
            Projective2 & operator = (Projective2 const & p)
            {
                *this << p;
                return *this;
            }

            /**
             * @brief 对 3x3 的矩阵执行拷贝操作
             */
            template <typename Derived>
            Projective2 & operator = (MatrixBase<Derived> const & m)
            {
                *this << m;
                return *this;
            }

            //! @brief 对点做射影变换
            HomoPoint2<DataType> ApplyOn(HomoPoint2<DataType> & p)
            {
                return *this * p;
            }

            //! @brief 对线做射影变换
            HomoLine2<DataType> ApplyOn(HomoLine2<DataType> const & l)
            {
                return this->InverseMat().Transpose() * l;
            }

            //! @brief 对圆锥曲线做射影变换
            inline HomoConic2<DataType> ApplyOn(HomoConic2<DataType> const & c)
            {
                auto H_inv = this->InverseMat();
                return H_inv.Transpose() * c * H_inv;
            }
    };
}

