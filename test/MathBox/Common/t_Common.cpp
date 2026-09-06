#include <iostream>

#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Common/Common.hpp>
#include <vector>

#include <cmath>
#include <gtest/gtest.h>



TEST(Common, IdxOfMaxAbs)
{
    std::vector<double> haha = { 0.1, 1.0, 0.5, -9.0, -4.0, 4.0 };
    size_t idx = xiaotu::IdxOfMaxAbs(haha);
    XTLog(std::cout) << idx << ":" << haha[idx] << std::endl;
}