#include <iostream>

#include <XiaoTuMathBox/Utils.hpp>
#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

#include <memory>
#include <gtest/gtest.h>
#include <Eigen/Eigen>

#include <vector>


TEST(LinearAlgibra, MatrixView)
{
    using namespace xiaotu::math;

    int len = 16;
    std::vector<double> buffer(len);
    for (int i = 0; i < len; i++)
        buffer[i] = 1.0 + i;

    MatrixView<double, 4, 4> m(buffer.data());

    for (int ridx = 0; ridx < 4; ridx++) {
        for (int cidx = 0; cidx < 4; cidx++) {
            m(ridx, cidx) = ridx + 0.1 * cidx;
        }
    }

    m.RowSwap(0, 2);
    for (int cidx = 0; cidx < 4; cidx++) {
        EXPECT_EQ(0 + 0.1 * cidx, m(2, cidx));
        EXPECT_EQ(2 + 0.1 * cidx, m(0, cidx));
    }

    m.RowSwap(0, 2);
    m.ColSwap(0, 2);
    for (int ridx = 0; ridx < 4; ridx++) {
        EXPECT_EQ(ridx + 0.1 * 0, m(ridx, 2));
        EXPECT_EQ(ridx + 0.1 * 2, m(ridx, 0));
    }

    MatrixView<double, 2, 3> m2(m.Ptr());
    m2 = { 1.0, 2.0, 3.0, 4.0, 5.0 , 6.0};
    XTLog(std::cout) << "m2 = " << m2 << std::endl;

    m2 = { { 10.0, 20.0 },
          { 30.0 } };
    XTLog(std::cout) << "m2 = " << m2 << std::endl;
}

TEST(LinearAlgibra, MatrixView_Dynamic)
{
    using namespace xiaotu::math;

    int len = 16;
    std::vector<double> buffer(len);
    for (int i = 0; i < len; i++)
        buffer[i] = 1.0 + i;

    MatrixView<double, Dynamic, Dynamic> m(buffer.data(), 4, 4);

    for (int ridx = 0; ridx < 4; ridx++) {
        for (int cidx = 0; cidx < 4; cidx++) {
            m(ridx, cidx) = ridx + 0.1 * cidx;
        }
    }

    m.RowSwap(0, 2);
    for (int cidx = 0; cidx < 4; cidx++) {
        EXPECT_EQ(0 + 0.1 * cidx, m(2, cidx));
        EXPECT_EQ(2 + 0.1 * cidx, m(0, cidx));
    }

    m.RowSwap(0, 2);
    m.ColSwap(0, 2);
    for (int ridx = 0; ridx < 4; ridx++) {
        EXPECT_EQ(ridx + 0.1 * 0, m(ridx, 2));
        EXPECT_EQ(ridx + 0.1 * 2, m(ridx, 0));
    }

    MatrixView<double, Dynamic, Dynamic> m2(m);
    m2.Reshape(2, 3);
    m2 = { 1.0, 2.0, 3.0, 4.0, 5.0 , 6.0};
    XTLog(std::cout) << "m2 = " << m2 << std::endl;

    m2 = { { 10.0, 20.0 },
          { 30.0 } };
    XTLog(std::cout) << "m2 = " << m2 << std::endl;
}



TEST(LinearAlgibra, GaussJordanEliminate)
{
    using namespace xiaotu::math;

    int len = 16;
    std::vector<double> _A_(len);
    for (int i = 0; i < len; i++)
        _A_[i] = 0.0;

    MatrixView<double, 4, 4> A(_A_.data());

    for (int idx = 0; idx < 4; idx++) {
        A(idx, idx) = idx + 1;
    }

    GaussJordanEliminate(A);
    for (int idx = 0; idx < 4; idx++) {
        EXPECT_EQ(1.0 / (idx + 1), A(idx, idx));
    }
    
    double s = std::sin(0.25 * M_PI);
    double c = std::cos(0.25 * M_PI);
    MatrixView<double, 3, 3> A1(_A_.data());
    A1 = { c, -s, 0,
           s,  c, 0,
           0,  0, 1 };

    // x 的解是 b 中向量顺时针旋转 45 度
    std::vector<double> _b_(6);
    MatrixView<double, 3, 2> b(_b_.data());
    b = { s, 1,
          s, 0,
          1, 1  };

    GaussJordanEliminate(A1, &b);
    EXPECT_EQ(A1(0, 0),  c); EXPECT_EQ(A1(0, 1), s); EXPECT_EQ(A1(0, 2), 0);
    EXPECT_EQ(A1(1, 0), -s); EXPECT_EQ(A1(1, 1), c); EXPECT_EQ(A1(1, 2), 0);
    EXPECT_EQ(A1(2, 0),  0); EXPECT_EQ(A1(2, 1), 0); EXPECT_EQ(A1(2, 2), 1);

    EXPECT_TRUE(std::abs(b(0, 0) - 1) < 1e-9);
    EXPECT_TRUE(std::abs(b(1, 0) - 0) < 1e-9);
    EXPECT_TRUE(std::abs(b(2, 0) - 1) < 1e-9);

    EXPECT_TRUE(std::abs(b(0, 1) - s) < 1e-9);
    EXPECT_TRUE(std::abs(b(1, 1) + s) < 1e-9);
    EXPECT_TRUE(std::abs(b(2, 1) - 1) < 1e-9);
}

//TEST(LinearAlgibra, LU)
//{
//    using namespace xiaotu::math;
//
//    std::vector<double> _A_(9);
//    std::vector<double> _LU_(9);
//    MatrixView<double> A(_A_.data(), _A_.size());
//    MatrixView<double> LU(_LU_.data(), _LU_.size());
//    A.Reshape(3, 3);
//    LU.Reshape(3, 3);
//
//    double s = std::sin(0.25 * M_PI);
//    double c = std::cos(0.25 * M_PI);
//    A = { {c, -s, 0 },
//          {s,  c, 0 },
//          {0,  0, 1 }};
//    LU = { c, -s, 0,
//          s,  c, 0,
//          0,  0, 1 };
//
//
//    LU_Decompose<double> lu(A, LU);
//    XTLog(std::cout) << "A:" << A;
//    XTLog(std::cout) << "lu:" << lu.LU();
//
//}

