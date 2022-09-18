#include <iostream>
#include <XiaoTuMathBox/LinearAlgibra/Matrix.hpp>

#include <gtest/gtest.h>

TEST(LinearAlgibra, Matrix)
{
    using namespace xiaotu::math;

    std::cout << "sizeof Matrix:" << sizeof(Matrix<double>) << std::endl;
    std::cout << "sizeof Matrix:" << sizeof(Matrix<int>) << std::endl;
    std::cout << "sizeof Matrix:" << sizeof(Matrix<float>) << std::endl;
}


