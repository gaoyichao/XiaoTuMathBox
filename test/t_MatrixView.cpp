#include <iostream>

#include <XiaoTuMathBox/Utils.hpp>
#include <XiaoTuMathBox/LinearAlgibra/MatrixViewBase.hpp>
#include <XiaoTuMathBox/LinearAlgibra/MatrixView.hpp>
#include <XiaoTuMathBox/LinearAlgibra/Matrix.hpp>
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
    EXPECT_EQ(sizeof(float) * buffer.size(), mv.Bytes());

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

    m.Reshape(2, 3);
    m = { 1.0, 2.0, 3.0, 4.0, 5.0 , 6.0};
    XTLog(std::cout) << "m = " << m << std::endl;

    m = { { 10.0, 20.0 },
          { 30.0 } };
    XTLog(std::cout) << "m = " << m << std::endl;
}

TEST(LinearAlgibra, Matrix)
{
    using namespace xiaotu::math;

    Matrix<double> *pm = new Matrix<double>();
    EXPECT_TRUE(pm->Empty());
    EXPECT_EQ(0, pm->Storage());
    EXPECT_EQ(Matrix<double>::EAlloc_No, pm->GetAllocFlags());

    pm->Resize(4, 4);
    EXPECT_EQ(16, pm->Storage());
    EXPECT_EQ(4, pm->Rows());
    EXPECT_EQ(4, pm->Cols());
    EXPECT_EQ(16, pm->NumDatas());
    EXPECT_TRUE(Matrix<double>::EAlloc_DataBuffer & pm->GetAllocFlags());
    EXPECT_TRUE(Matrix<double>::EAlloc_View & pm->GetAllocFlags());

    pm->Resize(3, 3);
    EXPECT_EQ(16, pm->Storage());
    EXPECT_EQ(3, pm->Rows());
    EXPECT_EQ(3, pm->Cols());
    EXPECT_EQ(9, pm->NumDatas());
    EXPECT_TRUE(Matrix<double>::EAlloc_DataBuffer & pm->GetAllocFlags());
    EXPECT_TRUE(Matrix<double>::EAlloc_View & pm->GetAllocFlags());

    pm->Full(0.618);
    XTLog(std::cout) << "pm->" << *pm << std::endl;

    pm->Resize(3, 4);
    XTLog(std::cout) << "pm->" << *pm << std::endl;

    *pm = { { 2.0, 3.0, 4.0, 6.0 },
            { 3.0 } };
    XTLog(std::cout) << "pm->" << *pm << std::endl;

    Matrix<double> m = { {0.1, 0.2, 0.3, 0.4 }, { 1.1 } };
    XTLog(std::cout) << "bienao->" << m << std::endl;

    double vlist[4] = { 1, 2, 3, 4};
    Matrix<double> m1(2, 2, vlist);
    XTLog(std::cout) << "m1->" << m1 << std::endl;

    Matrix<double> m2(m1);
    m2 = { { 1.2, 2.1 }, {2.1, 2.2 }};
    XTLog(std::cout) << "m1->" << m1 << std::endl;
    XTLog(std::cout) << "m2->" << m2 << std::endl;

    m2.Resize(4, 1);
    XTLog(std::cout) << "m1->" << m1 << std::endl;
    XTLog(std::cout) << "m2->" << m2 << std::endl;

    bool nothrow = true;
    try {
        m2.Resize(4, 2);
        nothrow = true;
    } catch (NotDeepCopyException e) {
        nothrow = false;
    }
    EXPECT_FALSE(nothrow);

    Matrix<double> const & _m3 = m1;
    Matrix<double> m3(_m3);
    nothrow = true;
    try {
        m3.Resize(4, 2);
        nothrow = true;
    } catch (NotDeepCopyException e) {
        nothrow = false;
    }
    EXPECT_TRUE(nothrow);

    m3 = { { 3.2, 3.1 }, {3.1, 3.2 }};
    XTLog(std::cout) << "m1->" << m1 << std::endl;
    XTLog(std::cout) << "m3->" << m3 << std::endl;

    m3 = *pm;
    XTLog(std::cout) << "m3->" << m3 << std::endl;

    m3 = { { 1.0, 2.0, 3.0, 6.0 },
           { 2.0, 3.0, 4.0, 7.0 },
           { 3.0, 4.0, 5.0, 8.0 } };
    XTLog(std::cout) << "pm->" << *pm << std::endl;

    m3 = *(Matrix<double> const *)pm;
    m3 = { { 1.1, 2.0, 3.0, 6.0 },
           { 2.1, 3.0, 4.0, 7.0 },
           { 3.1, 4.0, 5.0, 8.0 } };
    XTLog(std::cout) << "pm->" << *pm << std::endl;


    delete pm;
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

