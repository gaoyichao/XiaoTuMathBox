#include <iostream>

#include <XiaoTuMathBox/Geometry/Geometry.hpp>

#include <gtest/gtest.h>

using namespace xiaotu;

std::ostream & operator << (std::ostream & os, Polygon<double>::Vertex const & v)
{
    os << v.ToString() << ":" << *v << std::endl;
    return os;
}


TEST(Polygon, Empty)
{

    Polygon<double> p;
    EXPECT_EQ(0, p.NumVertices());

    {
        auto & v = p.Add(1.0, 2.0);
        EXPECT_EQ(1, p.NumVertices());
        std::cout << v << std::endl;
    }

    {
        auto & v = p.Add(2.0, 2.0);
        EXPECT_EQ(2, p.NumVertices());
        std::cout << v << std::endl;
        std::cout << v.Next() << std::endl;
        std::cout << v.Prev() << std::endl;
    }


}

