#ifndef XTMB_EQUATION_H
#define XTMB_EQUATION_H

#include <cmath>

namespace xiaotu {
namespace math {

    /*
     * QuadraticEquation - 二元一次方程组
     * 
     * ax^2 + bx + c = 0
     */
    template <typename T>
    class QuadraticEquation {
        public:
            QuadraticEquation(T const & _a, T const & _b, T const & _c)
                : a(_a), b(_b), c(_c), x1(0), x2(0), n(0)
            {
                T re = b * b - 4 * a * c;
                if (re < 0) {
                    n = 0;
                } else if (0 == re) {
                    x1 = x2 = -b / a * 0.5;
                    n = 1;
                } else {
                    x1 = (-b + std::sqrt(re)) / a * 0.5;
                    x2 = (-b - std::sqrt(re)) / a * 0.5;
                    n = 2;
                }
            }
    
        public:
            T a;
            T b;
            T c;
    
            T x1;
            T x2;
            int n;
    };

}
}

#endif

