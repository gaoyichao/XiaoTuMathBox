/********************************************************************************************************
 * 
 * Plucker 矩阵形式的线, 两点共线
 * 
 * 在 3-space 中有4个自由度
 *
 **************************************************************************** GAO YiChao 2021.1029 *****/
#ifndef XTMB_HOMOLINE3D_H
#define XTMB_HOMOLINE3D_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Eigen>

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

            //inline void SetValue(double _l12, double _l13, double _l14,
            //                     double _l23, double _l42, double _l34)
            //{
            //    l12 = _l12;
            //    l13 = _l13;
            //    l14 = _l14;
            //    l23 = _l23;
            //    l42 = _l42;
            //    l34 = _l34;
            //}

        public:
            double & l12;
            double & l13;
            double & l14;
            double & l23;
            double & l42;
            double & l34;
    };

}
}

#endif

