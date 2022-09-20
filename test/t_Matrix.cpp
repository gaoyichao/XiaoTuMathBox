#include <iostream>
#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

#include <gtest/gtest.h>
#include <Eigen/Eigen>

TEST(LinearAlgibra, Matrix)
{
    using namespace xiaotu::math;

    Matrix<double> m;
    EXPECT_TRUE(m.Empty());
    EXPECT_EQ(0, m.Storage());
    EXPECT_EQ(Matrix<double>::EAlloc_No, m.GetAllocFlags());

    m.Resize(3, 3);
    EXPECT_EQ(9, m.Storage());
    EXPECT_EQ(Matrix<double>::EAlloc_All, m.GetAllocFlags());

    m.Resize(4, 4);
    EXPECT_EQ(16, m.Storage());
    EXPECT_EQ(Matrix<double>::EAlloc_All, m.GetAllocFlags());

    m = {{2.0, 3.0},
         {1.0, 4.1}};
    EXPECT_EQ(16, m.Storage());

    m.Resize(4, 4);
    EXPECT_EQ(16, m.Storage());

    for (int ridx = 0; ridx < 4; ridx++) {
        for (int cidx = 0; cidx < 4; cidx++) {
            m[ridx][cidx] = ridx + 0.1 * cidx;
        }
    }

    m.RowSwap(0, 2);
    for (int cidx = 0; cidx < 4; cidx++) {
        EXPECT_EQ(0 + 0.1 * cidx, m[2][cidx]);
        EXPECT_EQ(2 + 0.1 * cidx, m[0][cidx]);
    }

    m.RowSwap(0, 2);
    m.ColSwap(0, 2);
    for (int ridx = 0; ridx < 4; ridx++) {
        EXPECT_EQ(ridx + 0.1 * 0, m[ridx][2]);
        EXPECT_EQ(ridx + 0.1 * 2, m[ridx][0]);
    }
}


TEST(LinearAlgibra, GaussJordanEliminate)
{
    using namespace xiaotu::math;

    Matrix<double> A(4, 4, 0.0);
    for (int idx = 0; idx < 4; idx++) {
        A[idx][idx] = idx + 1;
    }
    GaussJordanEliminate(A);
    for (int idx = 0; idx < 4; idx++) {
        EXPECT_EQ(1.0 / (idx + 1), A[idx][idx]);
    }
    
    double s = std::sin(0.25 * M_PI);
    double c = std::cos(0.25 * M_PI);
    A.Resize(3, 3);
    A = { c, -s, 0,
          s,  c, 0,
          0,  0, 1 };

    // x 的解是 b 中向量顺时针旋转 45 度
    Matrix<double> b(3, 2);
    b = { { s, 1 },
          { s, 0 },
          { 1, 1 } };

    GaussJordanEliminate(A, &b);
    EXPECT_EQ(A[0][0], c);  EXPECT_EQ(A[0][1], s); EXPECT_EQ(A[0][2], 0);
    EXPECT_EQ(A[1][0], -s); EXPECT_EQ(A[1][1], c); EXPECT_EQ(A[1][2], 0);
    EXPECT_EQ(A[2][0], 0);  EXPECT_EQ(A[2][1], 0); EXPECT_EQ(A[2][2], 1);

    EXPECT_TRUE(std::abs(b[0][0] - 1) < 1e-9);
    EXPECT_TRUE(std::abs(b[1][0] - 0) < 1e-9);
    EXPECT_TRUE(std::abs(b[2][0] - 1) < 1e-9);

    EXPECT_TRUE(std::abs(b[0][1] - s) < 1e-9);
    EXPECT_TRUE(std::abs(b[1][1] + s) < 1e-9);
    EXPECT_TRUE(std::abs(b[2][1] - 1) < 1e-9);
}

TEST(LinearAlgibra, LU)
{
    using namespace xiaotu::math;

    double s = std::sin(0.25 * M_PI);
    double c = std::cos(0.25 * M_PI);
    Matrix<double> A = { { c, -s, 0 },
                         { s,  c, 0 },
                         { 0,  0, 1 } };

    LU_Decompose<double> lu(A);
    XTLog(std::cout) << "A:" << A;
    XTLog(std::cout) << "lu:" << lu.LU();

}
