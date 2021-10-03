#include <XiaoTuMathBox/Point2D.h>


namespace xiaotu {
namespace math {

    void UnionVector2D(Point2D const & from, Point2D const & to, Vector2D & vec)
    {
        vec = to - from;
        vec.normalize();
    }

}
}

