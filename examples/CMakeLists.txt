cmake_minimum_required(VERSION 3.5.0)
project(XiaoTuMathBoxExample)

add_executable(douniwan ./douniwan.cpp)
target_link_libraries(douniwan XiaoTuMathBox)

