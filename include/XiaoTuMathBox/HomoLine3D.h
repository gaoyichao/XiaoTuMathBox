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
 *  2. Plücker 矩阵
 *   
 *    L = A * B^T - B * A^T
 *    
 *    该阵虽然是一个4*4的阵，但它的秩只有2，有4个自由度
 *    该阵与点或者面AB的选择无关
 *  
 *  3. Plücker line coordinates
 *  
 *    L = { l12, l13, l14, l23, l42, l34 }
 *
 **************************************************************************** GAO YiChao 2021.1029 *****/
#ifndef XTMB_HOMOLINE3D_H
#define XTMB_HOMOLINE3D_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Eigen>

#include <XiaoTuMathBox/HomoPoint3D.h>
#include <XiaoTuMathBox/HomoPlane3D.h>

namespace xiaotu {
namespace math {
    class HomoPoint3D;
    class HomoPlane3D;
 
    class HomoLine3D : public Eigen::Matrix4d
    {
        public:
            HomoLine3D()
                : l12(Eigen::Matrix4d::data()[1 * 4 + 0]),
                  l13(Eigen::Matrix4d::data()[2 * 4 + 0]),
                  l14(Eigen::Matrix4d::data()[3 * 4 + 0]),
                  l23(Eigen::Matrix4d::data()[2 * 4 + 1]),
                  l42(Eigen::Matrix4d::data()[1 * 4 + 3]),
                  l34(Eigen::Matrix4d::data()[3 * 4 + 2])
            {
                this->setZero();
            }

            HomoLine3D(HomoPoint3D const & p1, HomoPoint3D const & p2);

            inline void DualForm(Eigen::Matrix4d & dual) const
            {
                dual <<    0,  l34,  l42,  l23,
                        -l34,    0,  l14, -l13,
                        -l42, -l14,    0,  l12,
                        -l23,  l13, -l12,    0;
            }

            inline static Eigen::Matrix4d Transform(HomoLine3D const & L, Eigen::Matrix4d const & H)
            {
                return H * L * H.transpose();
            }

            inline HomoLine3D & Transform(Eigen::Matrix4d const & H)
            {
                *(Eigen::Matrix4d*)this = Transform(*this, H);
                return *this;
            }

            inline bool IsValid() const
            {
                double tmp = l12 * l34 + l13 * l42 + l14 * l23;
                return (std::abs(tmp) < 1e-9);
            }

            inline bool Coplanar(HomoLine3D const & l) const
            {
                double tmp = l12 * l.l34 + l.l12 * l34 + l13 * l.l42
                           + l.l13 * l42 + l14 * l.l23 + l.l14 * l23;

                return (std::abs(tmp) < 1e-9);
            }

            inline Eigen::Vector4d Intersection(HomoPlane3D const & pi) const
            {
                return (*this * pi);
            }

            inline Eigen::Vector4d JoinPlane(HomoPoint3D const p) const
            {
                Eigen::Matrix4d dual;
                DualForm(dual);
                return dual * p;
            }


        public:
            double const & l12;
            double const & l13;
            double const & l14;
            double const & l23;
            double const & l42;
            double const & l34;
    };

}
}

#endif

