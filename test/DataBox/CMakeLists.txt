cmake_minimum_required(VERSION 3.5.0)
project(TestXiaoTuMathBox)

set(TEST_DEPENDS gtest gtest_main pthread rt)

function(build_test_case NAME EXEC TSRC)
    add_executable(${EXEC} ${TSRC})
    target_link_libraries(${EXEC} ${TEST_DEPENDS} gflags)
    add_test(NAME ${NAME} COMMAND ${EXEC})
endfunction()

build_test_case(VectorList    t_VectorList    ./t_VectorList.cpp)

