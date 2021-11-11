#include <XiaoTuMathBox/ProjectiveGeometry.h>

namespace xiaotu {
namespace math {

    HomoLine3D::HomoLine3D(HomoPoint3D const & p1, HomoPoint3D const & p2)
        :  l12(Eigen::Matrix4d::data()[1 * 4 + 0]),
          l13(Eigen::Matrix4d::data()[2 * 4 + 0]),
          l14(Eigen::Matrix4d::data()[3 * 4 + 0]),
          l23(Eigen::Matrix4d::data()[2 * 4 + 1]),
          l42(Eigen::Matrix4d::data()[1 * 4 + 3]),
          l34(Eigen::Matrix4d::data()[3 * 4 + 2])
    {
        *this << (p1 * p2.transpose() - p2 * p1.transpose());
    }



}
}
