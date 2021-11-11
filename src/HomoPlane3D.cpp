#include <XiaoTuMathBox/ProjectiveGeometry.h>

namespace xiaotu {
namespace math {

    /*
     * HomoPlane3D - 构造函数, 由三个不同线的点构造一个平面
     */
    HomoPlane3D::HomoPlane3D(HomoPoint3D const & po1, HomoPoint3D const & po2, HomoPoint3D const & po3)
        : a(Eigen::Vector4d::data()[0]),
          b(Eigen::Vector4d::data()[1]),
          c(Eigen::Vector4d::data()[2]),
          d(Eigen::Vector4d::data()[3])
    {
        Eigen::Vector3d x1 = po1.Normalization().block(0, 0, 3, 1);
        Eigen::Vector3d x2 = po2.Normalization().block(0, 0, 3, 1);
        Eigen::Vector3d x3 = po3.Normalization().block(0, 0, 3, 1);

        Eigen::Vector3d x13, x23;
        x13 << (x1 - x3);
        x23 << (x2 - x3);

        *this << x13.cross(x23),
                 -x3.transpose() * x1.cross(x2);
        Normalize();
    }

    /*
     * Intersection - 仨平面的交点
     */
    Eigen::Vector4d Intersection(HomoPlane3D const & pi1, HomoPlane3D const & pi2, HomoPlane3D const & pi3)
    {
        double x = pi1.b * pi2.c * pi3.d + pi2.b * pi3.c * pi1.d + pi3.b * pi1.c * pi2.d
                 - pi3.b * pi2.c * pi1.d - pi2.b * pi1.c * pi3.d + pi1.b * pi3.c * pi2.d;
        double y = pi1.a * pi2.c * pi3.d + pi2.a * pi3.c * pi1.d + pi3.a * pi1.c * pi2.d
                 - pi3.a * pi2.c * pi1.d - pi2.a * pi1.c * pi3.d + pi1.a * pi3.c * pi2.d;
        double z = pi1.a * pi2.b * pi3.d + pi2.a * pi3.b * pi1.d + pi3.a * pi1.b * pi2.d
                 - pi3.a * pi2.b * pi1.d - pi2.a * pi1.b * pi3.d + pi1.a * pi3.b * pi2.d;
        double k = pi1.a * pi2.b * pi3.c + pi2.a * pi3.b * pi1.c + pi3.a * pi1.b * pi2.c
                 - pi3.a * pi2.b * pi1.c - pi2.a * pi1.b * pi3.c + pi1.a * pi3.b * pi2.c;
        y = -y;
        k = -k;

        return Eigen::Vector4d(x, y, z, k);
    }

}
}




