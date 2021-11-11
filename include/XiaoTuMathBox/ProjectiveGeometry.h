/********************************************************************************************************
 * 
 * 几何意义上：
 * 射影映射 <i>projective</i> 是一个从 \(\mathbb{P}^2\) 到其自身空间的一种可逆映射，不妨用符号\(h\)表示，
 * 即，若三点 \(\bold{x}_1, \bold{x}_2, \bold{x}_3\) 共线，
 * 则，\(h(\bold{x}_1), h(\bold{x}_2), h(\bold{x}_3)\)也共线
 * 
 * 代数意义上：
 * 映射 \(h\) 是空间\(\mathbb{P}^2\)上的一个射影映射 <==> 存在一个非奇异的3x3的矩阵 H, 使得 \(h(x) = Hx\)
 * 
 * 射影映射也称保线映射(collineation), 构成一个 group
 * 1. 射影映射的逆也是射影映射
 * 2. 射影映射的组合还是射影映射
 * 
 *                         H_s                            H_a                        H_p
 *          isometries --> similarity transformations --> affine transformations --> projective transformations
 *          等距           相似                           仿射                       射影
 * 不变量:  欧氏距离       形状                           平行关系
 *          3Dof           4Dof                           6Dof                       8Dof
 *          刚体                                                                     model vanishing points
 *                                                                                         灭点
 *                                                                                         
 * A projective transformation can be decomposed into a chain of transformations 
 * H = H_s H_a H_p
 * 
 **************************************************************************** GAO YiChao 2021.1029 *****/
#ifndef XTMB_PROJECTIVE_GEOMETRY_H
#define XTMB_PROJECTIVE_GEOMETRY_H

#include <cmath>
#include <iostream>
#include <vector>

#include <XiaoTuMathBox/EigenExt.h>

#include <XiaoTuMathBox/HomoPoint2D.h>
#include <XiaoTuMathBox/HomoLine2D.h>
#include <XiaoTuMathBox/HomoConic2D.h>

#include <XiaoTuMathBox/HomoPoint3D.h>
#include <XiaoTuMathBox/HomoPlane3D.h>
#include <XiaoTuMathBox/HomoLine3D.h>

namespace xiaotu {
namespace math {

    inline bool PointOn(HomoPoint2D const & p, HomoLine2D const & l, double tolerance = 1e-9)
    {
        return (std::fabs(p.transpose() * l) < tolerance);
    }

    inline bool PointOn(HomoPoint2D const & p, HomoConic2D const & c, double tolerance = 1e-9)
    {
        return (std::fabs(p.transpose() * c * p) < tolerance);
    }

    inline bool PointOn(HomoPoint3D const & po, HomoPlane3D const & pi, double tolerance = 1e-9)
    {
        return (std::fabs(po.transpose() * pi) < tolerance);
    }


    inline bool operator == (HomoPoint2D const & p1, HomoPoint2D const & p2)
    {
        if (std::fabs(p1.k - p2.k) < 1e-9) {
            return (std::fabs(p1.x - p2.x) < 1e-9) &&
                   (std::fabs(p1.y - p2.y) < 1e-9);
        } else {
            double k1_inv = 1.0 / p1.k;
            double k2_inv = 1.0 / p2.k;
            return (std::fabs(p1.x * k1_inv - p2.x * k2_inv) < 1e-9) &&
                   (std::fabs(p1.y * k1_inv - p2.y * k2_inv) < 1e-9);
        }
    }

    inline bool operator == (HomoLine2D const & l1, HomoLine2D const & l2)
    {
        if (std::fabs(l1.c - l2.c) < 1e-9) {
            return (std::fabs(l1.a - l2.a) < 1e-9) &&
                   (std::fabs(l1.b - l2.b) < 1e-9);
        } else {
            double c1_inv = 1.0 / l1.c;
            double c2_inv = 1.0 / l2.c;
            return (std::fabs(l1.a * c1_inv - l2.a * c2_inv) < 1e-9) &&
                   (std::fabs(l1.b * c1_inv - l2.b * c2_inv) < 1e-9);
        }
    }

    inline bool operator == (HomoPoint3D const & p1, HomoPoint3D const & p2)
    {
        if (std::fabs(p1.k - p2.k) < 1e-9) {
            return (std::fabs(p1.x - p2.x) < 1e-9) &&
                   (std::fabs(p1.y - p2.y) < 1e-9) &&
                   (std::fabs(p1.z - p2.z) < 1e-9);
        } else {
            double k1_inv = 1.0 / p1.k;
            double k2_inv = 1.0 / p2.k;
            return (std::fabs(p1.x * k1_inv - p2.x * k2_inv) < 1e-9) &&
                   (std::fabs(p1.y * k1_inv - p2.y * k2_inv) < 1e-9) &&
                   (std::fabs(p1.z * k1_inv - p2.z * k2_inv) < 1e-9);
        }
    }

    inline bool operator == (HomoPlane3D const & p1, HomoPlane3D const & p2)
    {
        if (std::fabs(p1.d - p2.d) < 1e-9) {
            return (std::fabs(p1.a - p2.a) < 1e-9) &&
                   (std::fabs(p1.b - p2.b) < 1e-9) &&
                   (std::fabs(p1.c - p2.c) < 1e-9);
        } else {
            double d1_inv = 1.0 / p1.d;
            double d2_inv = 1.0 / p2.d;
            return (std::fabs(p1.a * d1_inv - p2.a * d2_inv) < 1e-9) &&
                   (std::fabs(p1.b * d1_inv - p2.b * d2_inv) < 1e-9) &&
                   (std::fabs(p1.c * d1_inv - p2.c * d2_inv) < 1e-9);
        }
    }

}
}


#endif
