#ifndef XTMB_NUMERICAL_BANACHFIXEDPOINT
#define XTMB_NUMERICAL_BANACHFIXEDPOINT


#include <functional>
#include <XiaoTuMathBox/Common/Common.hpp>

namespace xiaotu {

    template <typename DataType>
    DataType NaiveFixedPoint(std::function<DataType(DataType)> g, DataType p0,
                       int max_iter = 100, DataType tol = SMALL_VALUE)
    {
        for (int i = 0; i < max_iter; ++i) {
            DataType p = g(p0);
            if (std::abs(p - p0) < tol)
                return p;
            p0 = p;
        }
        return p0;
    }



    template <typename DataType>
    DataType AitkenSteffensen(std::function<DataType(DataType)> g, DataType p0,
                       int max_iter = 100, DataType tol = SMALL_VALUE)
    {
        for (int i = 0; i < max_iter; ++i) {
            DataType p1 = g(p0);
            DataType p2 = g(p1);
            if (std::abs(p2 - p1) < tol)
                return p2;

            DataType p = p0 - (p1 - p0)*(p1 - p0)/(p2 - 2*p1 + p0);
            if (std::abs(p - p0) < tol)
                return p;
            p0 = p;
        }
        return p0;
    }


    
}


#endif
