/********************************************************************************************************
 * 
 * 在射影空间\(\mathbb{P}^2\)中，点和线对偶(dual) 
 * 在射影空间\(\mathbb{P}^3\)中，点和面对偶(dual) 
 *
 **************************************************************************** GAO YiChao 2021.1029 *****/
#include <XiaoTuMathBox/ProjectiveGeometry.h>

namespace xiaotu {
namespace math {

    HomoLine2D::HomoLine2D(HomoPoint2D const & p1, HomoPoint2D const & p2)
        : a(Eigen::Vector3d::data()[0]),
          b(Eigen::Vector3d::data()[1]),
          c(Eigen::Vector3d::data()[2])
    {
        *this = p1.JoinLine(p2);
    }


}
}


