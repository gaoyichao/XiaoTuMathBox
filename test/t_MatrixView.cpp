#include <iostream>
#include <XiaoTuMathBox/LinearAlgibra/MatrixViewBase.hpp>
#include <XiaoTuMathBox/LinearAlgibra/MatrixView.hpp>
#include <XiaoTuMathBox/LinearAlgibra/LinearAlgibra.hpp>

#include <gtest/gtest.h>
#include <Eigen/Eigen>

#include <vector>

TEST(LinearAlgibra, MatrixViewBase)
{
    using namespace xiaotu::math;

    int len = 4;
    std::vector<float> buffer(len);
    for (int i = 0; i < len; i++)
        buffer[i] = 1.0 + i;

    MatrixViewBase mv((uint8_t*)buffer.data(), sizeof(float) * buffer.size());

    mv.Reshape(1, 4, sizeof(float));
    EXPECT_EQ(1.0f, *(mv.Ptr<float>(0, 0)));
    EXPECT_EQ(2.0f, *(mv.Ptr<float>(0, 1)));
    EXPECT_EQ(3.0f, *(mv.Ptr<float>(0, 2)));
    EXPECT_EQ(4.0f, *(mv.Ptr<float>(0, 3)));

    mv.Reshape(4, 1, sizeof(float));
    EXPECT_EQ(1.0f, mv.At<float>(0, 0));
    EXPECT_EQ(2.0f, mv.At<float>(1, 0));
    EXPECT_EQ(3.0f, mv.At<float>(2, 0));
    EXPECT_EQ(4.0f, mv.At<float>(3, 0));

    mv.Reshape(2, 2, sizeof(float));
    mv.At<float>(0, 0) = 1.0f;
    mv.At<float>(0, 1) = 1.1f;
    mv.At<float>(1, 0) = 1.2f;
    mv.At<float>(1, 1) = 1.3f;

    mv.Reshape(2, 2, sizeof(float), MatrixViewBase::eRowMajor);
    EXPECT_EQ(1.0f, mv.At<float>(0, 0));
    EXPECT_EQ(1.2f, mv.At<float>(0, 1));
    EXPECT_EQ(1.1f, mv.At<float>(1, 0));
    EXPECT_EQ(1.3f, mv.At<float>(1, 1));

    mv.Full(3.14f);
    EXPECT_EQ(3.14f, mv.At<float>(0, 0));
    EXPECT_EQ(3.14f, mv.At<float>(0, 1));
    EXPECT_EQ(3.14f, mv.At<float>(1, 0));
    EXPECT_EQ(3.14f, mv.At<float>(1, 1));
}


TEST(LinearAlgibra, MatrixView)
{
    using namespace xiaotu::math;

    int len = 16;
    std::vector<double> buffer(len);
    for (int i = 0; i < len; i++)
        buffer[i] = 1.0 + i;

    MatrixView<double> m(buffer.data(), buffer.size());
    m.Reshape(4, 4);

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
}

TEST(LinearAlgibra, GaussJordanEliminate)
{
    using namespace xiaotu::math;

    int len = 16;
    std::vector<double> _A_(len);
    for (int i = 0; i < len; i++)
        _A_[i] = 0.0;

    MatrixView<double> A(_A_.data(), _A_.size());
    A.Reshape(4, 4);

    for (int idx = 0; idx < 4; idx++) {
        A(idx, idx) = idx + 1;
    }

    GaussJordanEliminate(A);
    for (int idx = 0; idx < 4; idx++) {
        EXPECT_EQ(1.0 / (idx + 1), A(idx, idx));
    }
    
    double s = std::sin(0.25 * M_PI);
    double c = std::cos(0.25 * M_PI);
    A.Reshape(3, 3);
    A = { c, -s, 0,
          s,  c, 0,
          0,  0, 1 };

    // x 的解是 b 中向量顺时针旋转 45 度
    std::vector<double> _b_(6);
    MatrixView<double> b(_b_.data(), _b_.size());
    b.Reshape(3, 2);
    b = { s, 1,
          s, 0,
          1, 1  };

    GaussJordanEliminate(A, &b);
    EXPECT_EQ(A(0, 0),  c); EXPECT_EQ(A(0, 1), s); EXPECT_EQ(A(0, 2), 0);
    EXPECT_EQ(A(1, 0), -s); EXPECT_EQ(A(1, 1), c); EXPECT_EQ(A(1, 2), 0);
    EXPECT_EQ(A(2, 0),  0); EXPECT_EQ(A(2, 1), 0); EXPECT_EQ(A(2, 2), 1);

    EXPECT_TRUE(std::abs(b(0, 0) - 1) < 1e-9);
    EXPECT_TRUE(std::abs(b(1, 0) - 0) < 1e-9);
    EXPECT_TRUE(std::abs(b(2, 0) - 1) < 1e-9);

    EXPECT_TRUE(std::abs(b(0, 1) - s) < 1e-9);
    EXPECT_TRUE(std::abs(b(1, 1) + s) < 1e-9);
    EXPECT_TRUE(std::abs(b(2, 1) - 1) < 1e-9);
}

TEST(LinearAlgibra, LU)
{
    using namespace xiaotu::math;

    std::vector<double> _A_(9);
    std::vector<double> _LU_(9);
    MatrixView<double> A(_A_.data(), _A_.size());
    MatrixView<double> LU(_LU_.data(), _LU_.size());
    A.Reshape(3, 3);
    LU.Reshape(3, 3);

    double s = std::sin(0.25 * M_PI);
    double c = std::cos(0.25 * M_PI);
    A = { c, -s, 0,
          s,  c, 0,
          0,  0, 1 };
    LU = { c, -s, 0,
          s,  c, 0,
          0,  0, 1 };


    LU_Decompose<double> lu(A, LU);
    XTLog(std::cout) << "A:" << A;
    XTLog(std::cout) << "lu:" << lu.LU();

}

