cmake_minimum_required(VERSION 3.5.0)
project(XiaoTuMathBoxExample)

add_executable(RungePhenomenon ./RungePhenomenon.cpp)
target_link_libraries(RungePhenomenon XiaoTuMathBox)

