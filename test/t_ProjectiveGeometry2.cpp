#include <iostream>
#include <XiaoTuMathBox/Utils.hpp>
#include <XiaoTuMathBox/Geometry/HomoUtils2.hpp>

#include <gtest/gtest.h>

TEST(Homogeneous, HomoPoint2)
{
    using namespace xiaotu::math;
    EXPECT_EQ(sizeof(double) * 3, sizeof(HomoPoint2<double>));

    HomoPoint2<double> p1(3.14159, 1.41421, 2.71828);
    EXPECT_DOUBLE_EQ(3.14159, p1.x());
    EXPECT_DOUBLE_EQ(1.41421, p1.y());
    EXPECT_DOUBLE_EQ(2.71828, p1.k());
    EXPECT_EQ(&p1(0), &p1.x());
    EXPECT_EQ(&p1(1), &p1.y());
    EXPECT_EQ(&p1(2), &p1.k());

    {
        HomoPoint2<double> p2 = p1;
        EXPECT_TRUE(p2 == p1);
        EXPECT_EQ(p2, p1);

        p2 = p2 * 0.1;
        XTLog(std::cout) << p1 << std::endl;
        XTLog(std::cout) << p2 << std::endl;
        EXPECT_TRUE(p2 == p1);
        EXPECT_EQ(p2, p1);
    }

    {
        HomoPoint2<double> p2 = p1.Normalization();
        EXPECT_DOUBLE_EQ(1.0, p2.k());
        EXPECT_DOUBLE_EQ(2.71828, p1.k());

        EXPECT_TRUE(p2 == p1);
        EXPECT_EQ(p2, p1);
    }

    {

        p1 = { 0.1, 1.0, 2.0 };
        HomoPoint2<double> p2 = (1.0 / p1.k()) * p1;
        EXPECT_DOUBLE_EQ(1.0, p2.k());
        EXPECT_DOUBLE_EQ(2.0, p1.k());

        EXPECT_TRUE(p2 == p1);
        EXPECT_EQ(p2, p1);
    }

    {
        EXPECT_FALSE(p1.IsInfinity());
        p1 = { 0.1, 1.0, 0.0 };
        EXPECT_TRUE(p1.IsInfinity());
    }
}


TEST(Homogeneous, HomoLine2)
{
    using namespace xiaotu::math;
    EXPECT_EQ(sizeof(double) * 3, sizeof(HomoLine2<double>));

    HomoLine2<double> l1(3.14159, 1.41421, 2.71828);
    EXPECT_DOUBLE_EQ(3.14159, l1.a());
    EXPECT_DOUBLE_EQ(1.41421, l1.b());
    EXPECT_DOUBLE_EQ(2.71828, l1.c());
    EXPECT_EQ(&l1(0), &l1.a());
    EXPECT_EQ(&l1(1), &l1.b());
    EXPECT_EQ(&l1(2), &l1.c());

    {
        HomoLine2<double> l2 = l1;
        EXPECT_TRUE(l2 == l1);
        EXPECT_EQ(l2, l1);

        l2 = l2 * 0.1;
        XTLog(std::cout) << l1 << std::endl;
        XTLog(std::cout) << l2 << std::endl;
        EXPECT_TRUE(l2 == l1);
        EXPECT_EQ(l2, l1);
    }


    {
        HomoLine2<double> l2 = l1.Normalization();
        EXPECT_DOUBLE_EQ(1.0, l2.Norm());
        EXPECT_DOUBLE_EQ(2.71828, l1.c());

        EXPECT_TRUE(l2 == l1);
        EXPECT_EQ(l2, l1);
    }

    {
        l1.SetValue(0.1, 1.0, 2.0);
        EXPECT_EQ(0.1, l1(0));
        EXPECT_EQ(1.0, l1(1));
        EXPECT_EQ(2.0, l1(2));

        HomoLine2<double> l2 = l1.Normalize();
        EXPECT_TRUE(l2 == l1);
        EXPECT_EQ(l2, l1);
    }

    {
        EXPECT_FALSE(l1.IsInfinity());
        l1 = { 0.0, 0.0, 1.0 };
        EXPECT_TRUE(l1.IsInfinity());
    }
}

TEST(Homogeneous, HomoConic2)
{
    using namespace xiaotu::math;
    EXPECT_EQ(sizeof(double) * 3 * 3, sizeof(HomoConic2<double>));

    HomoConic2<double> conic0(1, 2, 3, 4, 5, 6);
    EXPECT_DOUBLE_EQ(1.0, conic0.a());
    EXPECT_DOUBLE_EQ(2.0, conic0.b());
    EXPECT_DOUBLE_EQ(3.0, conic0.c());
    EXPECT_DOUBLE_EQ(4.0, conic0.d());
    EXPECT_DOUBLE_EQ(5.0, conic0.e());
    EXPECT_DOUBLE_EQ(6.0, conic0.f());

    {
        HomoConic2<double> conic1(conic0);
        EXPECT_DOUBLE_EQ(1.0, conic1.a());
        EXPECT_DOUBLE_EQ(2.0, conic1.b());
        EXPECT_DOUBLE_EQ(3.0, conic1.c());
        EXPECT_DOUBLE_EQ(4.0, conic1.d());
        EXPECT_DOUBLE_EQ(5.0, conic1.e());
        EXPECT_DOUBLE_EQ(6.0, conic1.f());
        EXPECT_TRUE(conic0 == conic1);
        EXPECT_EQ(conic0, conic1);

        conic1.SetValue(1.0, 2.0, 3.0, 4.0, 5.0, 7.0);
        EXPECT_DOUBLE_EQ(7.0, conic1.f());
        EXPECT_TRUE(conic0 != conic1);
        EXPECT_NE(conic0, conic1);
    }

    {
        HomoConic2<double> conic1 = 0.1 * conic0;
        EXPECT_DOUBLE_EQ(0.1, conic1.a());
        EXPECT_DOUBLE_EQ(0.2, conic1.b());
        EXPECT_DOUBLE_EQ(0.3, conic1.c());
        EXPECT_DOUBLE_EQ(0.4, conic1.d());
        EXPECT_DOUBLE_EQ(0.5, conic1.e());
        EXPECT_DOUBLE_EQ(0.6, conic1.f());
        EXPECT_TRUE(conic0 == conic1);
        EXPECT_EQ(conic0, conic1);
    }

    {
        Matrix<double, 3, 3> A = {
            9, -3, 1,
            1,  1, 1,
            4,  2, 1
        };
        HomoConic2<double> conic1(A);
        XTLog(std::cout) << conic1 << std::endl;
        EXPECT_FALSE(conic1.Valid());

        conic1 = conic0;
        EXPECT_DOUBLE_EQ(1.0, conic1.a());
        EXPECT_DOUBLE_EQ(2.0, conic1.b());
        EXPECT_DOUBLE_EQ(3.0, conic1.c());
        EXPECT_DOUBLE_EQ(4.0, conic1.d());
        EXPECT_DOUBLE_EQ(5.0, conic1.e());
        EXPECT_DOUBLE_EQ(6.0, conic1.f());
        EXPECT_TRUE(conic1.Valid());

        conic1 = A;
        XTLog(std::cout) << conic1 << std::endl;
        EXPECT_TRUE(conic0.Valid());
        EXPECT_FALSE(conic1.Valid());
    }

    {
        HomoPoint2<double> ps[5];
        ps[0].SetValue(0, -1, 1);
        ps[1].SetValue(1, 0, 1);
        ps[2].SetValue(-1, 0, 1);
        ps[3].SetValue(2, 3, 1);
        ps[4].SetValue(-2, 3, 1);
        HomoConic2<double> conic1 = CoConic(ps);

        EXPECT_TRUE(OnConic(ps[0], conic1, SMALL_VALUE));
        EXPECT_TRUE(OnConic(ps[1], conic1, SMALL_VALUE));
        EXPECT_TRUE(OnConic(ps[2], conic1, SMALL_VALUE));
        EXPECT_TRUE(OnConic(ps[3], conic1, SMALL_VALUE));
        EXPECT_TRUE(OnConic(ps[4], conic1, SMALL_VALUE));

        conic1.Normalize();
        EXPECT_TRUE(OnConic(ps[0], conic1, SMALL_VALUE));
        EXPECT_TRUE(OnConic(ps[1], conic1, SMALL_VALUE));
        EXPECT_TRUE(OnConic(ps[2], conic1, SMALL_VALUE));
        EXPECT_TRUE(OnConic(ps[3], conic1, SMALL_VALUE));
        EXPECT_TRUE(OnConic(ps[4], conic1, SMALL_VALUE));

        HomoLine2<double> l = conic1.Tangent(ps[0]);
        EXPECT_EQ(l, HomoLine2<double>(0, 1, 1));

        EXPECT_TRUE(std::abs(l.a()) < SMALL_VALUE);
        EXPECT_DOUBLE_EQ(0.5, l.b());
        EXPECT_DOUBLE_EQ(0.5, l.c());
    }
}

TEST(Homogeneous, HomoUtils2)
{
    using namespace xiaotu::math;

    HomoPoint2<float> p0, p1;
    p0 << 1.0f, 0.0f, 1.0f;
    p1 << 2.0f, 0.0f, 1.0f;
    HomoLine2<float> l = Collinear(p0, p1);

    HomoPoint2<float> p2(3.0, 0.0, 1.0);
    EXPECT_TRUE(OnLine(p2, l));
    p2 << 1.0f, 0.1f, 1.0f;
    EXPECT_FALSE(OnLine(p2, l));

    HomoLine2<float> l2(1.0, 0.0, 0.0);
    HomoPoint2<float> p = Intersection(l, l2);
    p2 << 0.0f, 0.0f, 1.0f;
    EXPECT_EQ(p, p2);

    l << -1.0f, 0.0f, 1.0f;
    l2 << 0.0f, -1.0f, 1.0f;
    p = Intersection(l, l2);
    p2 << 1.0f, 1.0f, 1.0f;
    EXPECT_EQ(p, p2);

    p << 0.0f, 0.f, 1.0f;
    EXPECT_TRUE((Distance(p, l) - 1) < 1e-6);

    l *= 2;
    EXPECT_TRUE((Distance(p, l) - 1) < 1e-6);

    l << 1.0f, 1.0f, 1.0f;
    EXPECT_TRUE((Distance(p, l) - std::sqrt(2)) < 1e-6);
}

TEST(Homogeneous, Infinity)
{
    using namespace xiaotu::math;
    float a = 1.0;
    float b = 2.0;
    float c = 3.0;
    float _c = 5.0;

    //! 两条平行线的交点在无穷远处
    HomoLine2<float> l1(a, b, c);
    HomoLine2<float> l2(a, b, _c);
    HomoPoint2<float> p = Intersection(l1, l2);
    EXPECT_TRUE(p.IsInfinity());

    float k = _c - c;
    p = p / k;
    EXPECT_EQ(b, p.x());
    EXPECT_EQ(-a, p.y());
    EXPECT_EQ(0, p.k());
 
    // 正常的直线与无穷远处的直线的交点是无穷远处的点
    l1 << 0.0f, 0.0f, 1.0f;
    l2 << a, b, c;
    p = Intersection(l1, l2);
    EXPECT_TRUE(p.IsInfinity());
}

// TEST(Homogeneous, Projective2)
// {
//     using namespace xiaotu::math;
//     EXPECT_EQ(sizeof(Eigen::Matrix3d), sizeof(Projective2<double>));

//     Projective2<double> H;
//     H = Eigen::Matrix3d::Identity();
//     EXPECT_TRUE(H.isIdentity());

//     double rad = M_PI * 0.5;
//     double c = std::cos(rad);
//     double s = std::sin(rad);

//     H << c, -s, 0,
//          s,  c, 0,
//          0,  0, 1;

//     HomoPoint2<double> p0(1, 0, 1);
//     HomoPoint2<double> p1(1, 1, 1);
//     HomoPoint2<double> p2 = H.ApplyOn(p0);
//     EXPECT_EQ(p2, HomoPoint2<double>(0, 1, 1));

//     HomoLine2<double> l = Collinear(p0, p1);
//     HomoLine2<double> _l = H.ApplyOn(l);
//     EXPECT_TRUE(OnConic(p2, _l));

//     HomoPoint2<double> ps[5];
//     ps[0].SetValue(0, -1, 1);
//     ps[1].SetValue(1, 0, 1);
//     ps[2].SetValue(-1, 0, 1);
//     ps[3].SetValue(2, 3, 1);
//     ps[4].SetValue(-2, 3, 1);

//     HomoConic2<double> conic1;
//     conic1.From5Points(ps);
//     HomoConic2<double> conic2 = H.ApplyOn(conic1);

//     ps[0] = H.ApplyOn(ps[0]);
//     ps[1] = H.ApplyOn(ps[1]);
//     ps[2] = H.ApplyOn(ps[2]);
//     ps[3] = H.ApplyOn(ps[3]);
//     ps[4] = H.ApplyOn(ps[4]);

//     EXPECT_TRUE(OnConic(ps[0], conic2, 1e-6));
//     EXPECT_TRUE(OnConic(ps[1], conic2, 1e-6));
//     EXPECT_TRUE(OnConic(ps[2], conic2, 1e-6));
//     EXPECT_TRUE(OnConic(ps[3], conic2, 1e-6));
//     EXPECT_TRUE(OnConic(ps[4], conic2, 1e-6));

//     conic1.From5Points(ps);
//     EXPECT_EQ(conic1, conic2);
// }

