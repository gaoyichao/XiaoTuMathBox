/********************************************************************************************************
 * 
 * 在射影空间\(\mathbb{P}^2\)中，点和线对偶(dual) 
 * 在射影空间\(\mathbb{P}^3\)中，点和面对偶(dual) 
 *
 **************************************************************************** GAO YiChao 2021.1029 *****/
#include <XiaoTuMathBox/ProjectiveGeometry.h>

namespace xiaotu {
namespace math {

    HomoPoint2D::HomoPoint2D(HomoLine2D const & l1, HomoLine2D const & l2)
        : x(Eigen::Vector3d::data()[0]),
          y(Eigen::Vector3d::data()[1]),
          k(Eigen::Vector3d::data()[2])
    {
        *this = l1.cross(l2);
    }


}
}


