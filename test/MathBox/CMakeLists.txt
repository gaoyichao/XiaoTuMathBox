cmake_minimum_required(VERSION 3.5.0)
project(TestXiaoTuMathBox)

set(TEST_DEPENDS gtest gtest_main pthread rt)

function(build_test_case NAME EXEC TSRC)
    add_executable(${EXEC} ${TSRC})
    target_link_libraries(${EXEC} ${TEST_DEPENDS} gflags)
    add_test(NAME ${NAME} COMMAND ${EXEC})
endfunction()

build_test_case(MatrixView    t_MatrixView    ./LinearAlgibra/t_MatrixView.cpp)
build_test_case(MatrixSubView t_MatrixSubView ./LinearAlgibra/t_MatrixSubView.cpp)
build_test_case(Matrix        t_Matrix        ./LinearAlgibra/t_Matrix.cpp)
build_test_case(Permutation   t_Permutation   ./LinearAlgibra/t_Permutation.cpp)
build_test_case(DMatrixView   t_DMatrixView   ./LinearAlgibra/t_DMatrixView.cpp)
build_test_case(DMatrix       t_DMatrix       ./LinearAlgibra/t_DMatrix.cpp)
build_test_case(MatrixComma   t_MatrixComma   ./LinearAlgibra/t_MatrixComma.cpp)
build_test_case(MatrixInverse t_MatrixInverse ./LinearAlgibra/t_MatrixInverse.cpp)
build_test_case(QR            t_QR            ./LinearAlgibra/t_QR.cpp)
build_test_case(Eigen         t_Eigen         ./LinearAlgibra/t_Eigen.cpp)
build_test_case(SVD           t_SVD           ./LinearAlgibra/t_SVD.cpp)

build_test_case(Euclidean2    t_Euclidean2    ./Geometry/t_Euclidean2.cpp)
build_test_case(Euclidean3    t_Euclidean3    ./Geometry/t_Euclidean3.cpp)
build_test_case(Projective2   t_Projective2   ./Geometry/t_Projective2.cpp)
build_test_case(Projective3   t_Projective3   ./Geometry/t_Projective3.cpp)

build_test_case(Polygon   t_Polygon   ./Geometry/t_Polygon.cpp)


