#include <iostream>

#include <XiaoTuMathBox/HomoUtils3.hpp>

#include <gtest/gtest.h>

TEST(Homogeneous, HomoPoint3)
{
    using namespace xiaotu::math;
    EXPECT_EQ(sizeof(Eigen::Vector4d), sizeof(HomoPoint3<double>));

    HomoPoint3<double> p1(3.14159, 1.41421, 2.71828, 0.5);
    EXPECT_EQ(3.14159, p1.x());
    EXPECT_EQ(1.41421, p1.y());
    EXPECT_EQ(2.71828, p1.z());
    EXPECT_EQ(0.5, p1.k());
    EXPECT_EQ(&p1[0], &p1.x());
    EXPECT_EQ(&p1[1], &p1.y());
    EXPECT_EQ(&p1[2], &p1.z());
    EXPECT_EQ(&p1[3], &p1.k());

    p1.SetValue(0.1, 1.0, 2.0, 3.0);
    EXPECT_EQ(0.1, p1[0]);
    EXPECT_EQ(1.0, p1[1]);
    EXPECT_EQ(2.0, p1[2]);
    EXPECT_EQ(3.0, p1[3]);

    HomoPoint3<double> p2 = p1.Normalize();
    EXPECT_EQ(1.0, p2.k());
    EXPECT_EQ(1.0, p1.k());

    p1 << 0.1, 1.0, 2.0, 3.0;
    p2 = p1.Normalization();
    EXPECT_EQ(1.0, p2.k());
    EXPECT_EQ(3.0, p1.k());

    EXPECT_FALSE(p1.IsInfinity());
    EXPECT_FALSE(p2.IsInfinity());

    p1 << 0.1, 1.0, 2.0, 0.0;
    EXPECT_TRUE(p1.IsInfinity());
}


TEST(Homogeneous, HomoPlane3)
{
    using namespace xiaotu::math;
    EXPECT_EQ(sizeof(Eigen::Vector4d), sizeof(HomoPlane3<double>));

    HomoPlane3<double> l1(0.1, 1.0, 0.0, 0.0);
    EXPECT_EQ(0.1, l1.a());
    EXPECT_EQ(1.0, l1.b());
    EXPECT_EQ(0.0, l1.c());
    EXPECT_EQ(0.0, l1.d());
    EXPECT_EQ(&l1[0], &l1.a());
    EXPECT_EQ(&l1[1], &l1.b());
    EXPECT_EQ(&l1[2], &l1.c());
    EXPECT_EQ(&l1[3], &l1.d());

    l1.SetValue(0.1, 1.0, 2.0, 1.0);
    EXPECT_EQ(0.1, l1[0]);
    EXPECT_EQ(1.0, l1[1]);
    EXPECT_EQ(2.0, l1[2]);
    EXPECT_EQ(1.0, l1[3]);

    HomoPlane3<double> l2 = l1.Normalize();
    EXPECT_EQ(1.0, l2.c());
    EXPECT_EQ(1.0, l1.c());

    l1 << 0.1, 1.0, 2.0, 1.0;
    l2 = l1.Normalization();
    EXPECT_EQ(1.0, l2.c());
    EXPECT_EQ(2.0, l1.c());

    EXPECT_FALSE(l1.IsInfinity());
    EXPECT_FALSE(l2.IsInfinity());

    l1 << 0.0, 0.0, 0.0, 1.0;
    EXPECT_TRUE(l1.IsInfinity());
}

TEST(Homogeneous, HomoLine3)
{
    using namespace xiaotu::math;
    EXPECT_EQ(sizeof(Eigen::Matrix4d), sizeof(HomoLine3<double>));

    HomoLine3<double> line;
    line << 11, 12, 13, 14,
            21, 22, 23, 24,
            31, 32, 33, 34,
            41, 42, 43, 44;
    EXPECT_EQ(12, line.l12());
    EXPECT_EQ(13, line.l13());
    EXPECT_EQ(14, line.l14());
    EXPECT_EQ(23, line.l23());
    EXPECT_EQ(42, line.l42());
    EXPECT_EQ(34, line.l34());

    line.SetValue(12, 13, 14, 23, 42, 34);
    EXPECT_EQ(12, line.l12());
    EXPECT_EQ(13, line.l13());
    EXPECT_EQ(14, line.l14());
    EXPECT_EQ(23, line.l23());
    EXPECT_EQ(42, line.l42());
    EXPECT_EQ(34, line.l34());

    EXPECT_EQ(-12, line(1, 0));
    EXPECT_EQ(-13, line(2, 0));
    EXPECT_EQ(-14, line(3, 0));
    EXPECT_EQ(-23, line(2, 1));
    EXPECT_EQ(-42, line(1, 3));
    EXPECT_EQ(-34, line(3, 2));

    EXPECT_FALSE(line.IsValid());
}


TEST(Homogeneous, HomoUtils3)
{
    using namespace xiaotu::math;

    HomoPlane3<double> pi(0.1, 1.0, 1.0, 0.0);
    EXPECT_EQ(&pi[0], &pi.a());
    EXPECT_EQ(&pi[1], &pi.b());
    EXPECT_EQ(&pi[2], &pi.c());
    EXPECT_EQ(&pi[3], &pi.d());

    HomoPoint3<double> po(0.0, 0.0, 0.0, 1.0);
    EXPECT_TRUE(OnPlane(po, pi));

    HomoPoint3<double> po1(1.0,  2.0, 0.0, 1.0);
    HomoPoint3<double> po2(2.0,  2.0, 0.0, 1.0);
    HomoPoint3<double> po3(2.0, -2.0, 0.0, 1.0);
    HomoPlane3<double> pi1 = Coplanar(po1, po2, po3);

    EXPECT_TRUE(OnPlane(po, pi1));
    po.SetValue(-9.0, 11.0, 0.0, 1.0);
    EXPECT_TRUE(OnPlane(po, pi1));
    EXPECT_TRUE(OnPlane(po1, pi1));
    EXPECT_TRUE(OnPlane(po2, pi1));
    EXPECT_TRUE(OnPlane(po3, pi1));

    pi = CoplanarSVD(po1, po2, po3);
    EXPECT_EQ(pi, pi1);
    EXPECT_TRUE(OnPlane(po, pi));

    pi1.SetValue(0.0, 0.0, 1.0, 0.0);
    HomoPlane3<double> pi2(1.0, 0.0, 0.0, 0.0);
    HomoPlane3<double> pi3(0.0, 1.0, 0.0, 0.0);

    po1 = Intersection(pi1, pi2, pi3);
    EXPECT_EQ(po1, HomoPoint3<double>(0, 0, 0, 1));

    po2 = IntersectionSVD(pi1, pi2, pi3);
    EXPECT_EQ(po1, po2);

    HomoLine3<double> line;
    po1 << 0.0, 0.0, 0.0, 1.0;
    po2 << 1.0, 0.0, 0.0, 0.0;
    po3 << 3.0, 0.0, 0.0, 1.0;
    line = Collinear(po1, po2);
    EXPECT_EQ(-1, line.l14());
    EXPECT_TRUE(line.IsValid());
    EXPECT_TRUE(OnLine(po3, line));
 
    pi1 << 0.0, 1.0, 0.0, 0.0;
    pi2 << 0.0, 0.0, 1.0, 0.0;
    HomoLine3<double> line2 = Intersection(pi1, pi2);

    EXPECT_TRUE(line.AreCoplanar(line2));
    EXPECT_TRUE(OnLine(po1, line2));
    EXPECT_TRUE(OnLine(po2, line2));

    EXPECT_TRUE(OnPlane(line, pi1));
    EXPECT_TRUE(OnPlane(line, pi2));

    pi2 << 1.0, 0.0, 0.0, -1.0;
    po3 = Intersection(line, pi2);
    EXPECT_EQ(po3, HomoPoint3<double>(1, 0, 0, 1));

    po3 << 0.0, 1.0, 0.0, 1.0;
    pi1 = Coplanar(po3, line);
    EXPECT_EQ(pi1, HomoPlane3<double>(0, 0, 2, 0));
}

TEST(Homogeneous, Projective3)
{
    using namespace xiaotu::math;
    EXPECT_EQ(sizeof(Eigen::Matrix4d), sizeof(Projective3<double>));

    Projective3<double> H;
    H = Eigen::Matrix4d::Identity();
    EXPECT_TRUE(H.isIdentity());

    double rad = M_PI * 0.5;
    double c = std::cos(rad);
    double s = std::sin(rad);

    H << c, -s, 0, 0,
         s,  c, 0, 0,
         0,  0, 1, 0,
         0,  0, 0, 1;
    HomoPoint3<double> p0(1, 0, 0, 1);
    HomoPoint3<double> p1(1, 1, 0, 1);
    HomoPoint3<double> p2(0, 1, 0, 1);
    HomoPoint3<double> p3 = H.ApplyOn(p0);
    EXPECT_EQ(p3, HomoPoint3<double>(0, 1, 0, 1));

    HomoPlane3<double> pi = Coplanar(p0, p1, p2);
    HomoPlane3<double> _pi = H.ApplyOn(pi);
    EXPECT_TRUE(OnPlane(p3, _pi));

    HomoLine3<double> line;
    p0 << 0.0, 0.0, 0.0, 1.0;
    p1 << 1.0, 0.0, 0.0, 0.0;
    p2 << 3.0, 0.0, 0.0, 1.0;
    line = Collinear(p0, p1);
    EXPECT_TRUE(line.IsValid());
    EXPECT_TRUE(OnLine(p2, line));
 
    HomoLine3<double> _line = H.ApplyOn(line);
    p2 = H.ApplyOn(p2);
    EXPECT_TRUE(_line.IsValid());
    EXPECT_TRUE(OnLine(p2, _line));
}

TEST(Homogeneous, Properties)
{
    using namespace xiaotu::math;
    EXPECT_EQ(sizeof(Eigen::Matrix4d), sizeof(Projective3<double>));

    HomoPlane3<double> pi1(1.0, 0.0, 0.0, 0.0);
    HomoPlane3<double> pi2(3.0, 0.0, 0.0, 0.0);
    HomoLine3<double> line1 = Intersection(pi1, pi2);

    EXPECT_TRUE(OnPlane(line1, HomoPlane3<double>(0, 0, 0, 1)));


}

