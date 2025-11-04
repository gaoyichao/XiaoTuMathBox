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
# tests 
#
############################################################################

enable_testing()
add_subdirectory(./test/MathBox)
add_subdirectory(./test/DataBox)

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

# INSTALL(TARGETS ${PROJECT_NAME}
#     LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
#     ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
# )

INSTALL(DIRECTORY "${PROJECT_SOURCE_DIR}/include"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/"
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.hpp"
)

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

