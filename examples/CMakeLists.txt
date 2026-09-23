cmake_minimum_required(VERSION 3.5.0)
project(XiaoTuMathBoxExample)

find_package(Matplot++ REQUIRED)

add_executable(RungePhenomenon ./RungePhenomenon.cpp)
target_link_libraries(RungePhenomenon XiaoTuMathBox matplot)

add_executable(NaivePieceCubicHermite ./NaivePieceCubicHermite.cpp)
target_link_libraries(NaivePieceCubicHermite XiaoTuMathBox matplot)

