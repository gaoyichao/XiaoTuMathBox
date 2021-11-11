#include <XiaoTuMathBox/ProjectiveGeometry.h>
#include <iostream>

namespace xiaotu {
namespace math {

    const HomoPoint3D HomoPoint3D::Infinity = HomoPoint3D(0, 0, 0, 0);
    const HomoPoint3D HomoPoint3D::XInfinity = HomoPoint3D(1, 0, 0, 0);
    const HomoPoint3D HomoPoint3D::YInfinity = HomoPoint3D(0, 1, 0, 0);
    const HomoPoint3D HomoPoint3D::ZInfinity = HomoPoint3D(0, 0, 1, 0);


    /*
     * HomoPoint3D - 构造函数, 由三个不同的平面交一个点
     */
    HomoPoint3D::HomoPoint3D(HomoPlane3D const & pi1, HomoPlane3D const & pi2, HomoPlane3D const & pi3)
        : x(Eigen::Vector4d::data()[0]),
          y(Eigen::Vector4d::data()[1]),
          z(Eigen::Vector4d::data()[2]),
          k(Eigen::Vector4d::data()[3])
    {
        *this << Intersection(pi1, pi2, pi3);
    }


}
}




