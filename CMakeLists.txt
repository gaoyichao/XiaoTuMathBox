cmake_minimum_required(VERSION 3.5.0)
project(XiaoTuMathBox)
set(PACKAGE_VERSION "0.0.1")

string(REGEX REPLACE "([0-9]+).*$" "\\1" PACKAGE_MAJOR_VERSION ${PACKAGE_VERSION})
string(REGEX REPLACE "[0-9]+\\.([0-9]+).*$" "\\1" PACKAGE_MINOR_VERSION ${PACKAGE_VERSION})
string(REGEX REPLACE "[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1" PACKAGE_PATCH_VERSION ${PACKAGE_VERSION})

add_compile_options(-std=c++17)
set(CMAKE_BUILD_TYPE DEBUG)

message("${PROJECT_NAME}")
message("${PROJECT_SOURCE_DIR}")
message("${PROJECT_BINARY_DIR}")
message("${CMAKE_INSTALL_PREFIX}")

include_directories(
    /home/gyc/local/include/eigen3
    ${PROJECT_SOURCE_DIR}/include
)

############################################################################
#
# libXiaoTuMathBox.so
#
############################################################################

file(GLOB PROJECT_SRC_FILES
          ${PROJECT_SOURCE_DIR}/src/*cpp
)

add_library(${PROJECT_NAME} SHARED
    ${PROJECT_SRC_FILES}
)
target_link_libraries(${PROJECT_NAME} pthread)

############################################################################
#
# utils
#
############################################################################

# echo 服务器 - 原始 poll 接口实现
#add_executable(u_raw_poll_echo utils/u_raw_poll_echo.cpp)
#target_link_libraries(u_raw_poll_echo ${PROJECT_NAME})

############################################################################
#
# tests 
#
############################################################################

enable_testing()

set(TEST_DEPENDS gtest gtest_main pthread rt)

function(build_test_case NAME EXEC TSRC)
    message(STATUS "build_test_case:"  ${NAME})
    message(STATUS "EXEC = " ${EXEC})
    message(STATUS "TSRC = " ${TSRC})

    add_executable(${EXEC} ${TSRC})
    target_link_libraries(${EXEC} ${TEST_DEPENDS} ${PROJECT_NAME} gflags)
    add_test(NAME ${NAME} COMMAND ${EXEC})
endfunction()

# build_test_case(ProjectiveGeometry3 t_ProjectiveGeometry3 ./test/t_ProjectiveGeometry3.cpp)

build_test_case(MatrixView t_MatrixView ./test/t_MatrixView.cpp)
build_test_case(MatrixSubView t_MatrixSubView ./test/t_MatrixSubView.cpp)
build_test_case(Matrix t_Matrix ./test/t_Matrix.cpp)
build_test_case(DMatrixView t_DMatrixView ./test/t_DMatrixView.cpp)
build_test_case(DMatrix t_DMatrix ./test/t_DMatrix.cpp)
build_test_case(MatrixComma t_MatrixComma ./test/t_MatrixComma.cpp)

build_test_case(QR t_QR ./test/t_QR.cpp)
build_test_case(Eigen t_Eigen ./test/t_Eigen.cpp)

build_test_case(Euclidean2 t_Euclidean2 ./test/t_Euclidean2.cpp)
build_test_case(ProjectiveGeometry2 t_ProjectiveGeometry2 ./test/t_ProjectiveGeometry2.cpp)


############################################################################
#
# documentation 
#
############################################################################

include(cmake/build_doxygen.cmake)
build_doxygen()

############################################################################
#
# install 
#
############################################################################

include(cmake/install_package.cmake)
install_package()

INSTALL(TARGETS ${PROJECT_NAME}
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

file(GLOB PROJECT_INC_FILES
          ${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME}/*.h
          ${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME}/*.hpp
)
INSTALL(FILES ${PROJECT_INC_FILES}
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/${PROJECT_NAME}
)

function(install_include_subset SUBDIR)
    message(STATUS "install_include_subset:"  ${SUBDIR})
    file(GLOB SUBDIR_INC_FILES
        ${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME}/${SUBDIR}/*.h
        ${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME}/${SUBDIR}/*.hpp
    )

    INSTALL(FILES ${SUBDIR_INC_FILES}
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include/${PROJECT_NAME}/${SUBDIR}
    )
endfunction()


add_custom_target("uninstall" COMMENT "Uninstall installed files")
add_custom_command(
    TARGET "uninstall"
    POST_BUILD
    COMMENT "Uninstall files with install_manifest.txt"
    COMMAND xargs rm -vf < install_manifest.txt || echo Nothing in
            install_manifest.txt to be uninstalled!
)

#######################################################################
#
# deb 打包 
#
#######################################################################
include(./cmake/deb.cmake)
include(CPack)

