#include <iostream>

#include <XiaoTuMathBox/ProjectiveGeometry.h>

#include <gtest/gtest.h>

TEST(Homogeneous, Point2D)
{
    using namespace xiaotu::math;

    HomoPoint2D bienao(0.1, 1.0, 2.0);
    EXPECT_EQ(&bienao[0], &bienao.x);
    EXPECT_EQ(&bienao[1], &bienao.y);
    EXPECT_EQ(&bienao[2], &bienao.k);

    HomoPoint2D douniwan = bienao.Normalization();
    EXPECT_EQ(1.0, douniwan.k);
    EXPECT_EQ(2.0, bienao.k);
    EXPECT_TRUE(bienao == douniwan);
    bienao.Normalize();
    EXPECT_EQ(1.0, bienao.k);

    bienao.SetValue(1.0, 0.0, 1.0);
    douniwan.SetValue(0.0, 1.0, 1.0);
    HomoLine2D l1(bienao, douniwan);
    HomoLine2D l2 = douniwan.JoinLine(bienao);
    EXPECT_EQ(l1, l2);

    HomoPoint2D p0 = HomoPoint2D::Infinity();
    EXPECT_TRUE(p0.IsInfinity());
}

TEST(Homogeneous, Line2D)
{
    using namespace xiaotu::math;

    HomoLine2D l(0.1, 1.0, 0.0);
    EXPECT_EQ(&l[0], &l.a);
    EXPECT_EQ(&l[1], &l.b);
    EXPECT_EQ(&l[2], &l.c);

    HomoPoint2D p(0.0, 0.0, 1.0);
    EXPECT_TRUE(PointOn(p, l));

    HomoLine2D l0(1.1, 1.0, 1.0);
    p = l.Intersection(l0);
    EXPECT_TRUE(PointOn(p, l0));
    EXPECT_TRUE(PointOn(p, l));
    p.Normalize();
    EXPECT_EQ(1.0, p.k);

    HomoPoint2D p0 = l0.Intersection(l);
    EXPECT_TRUE(PointOn(p0, l0));
    EXPECT_TRUE(PointOn(p0, l));
    p0.Normalize();
    EXPECT_EQ(1.0, p0.k);

    EXPECT_TRUE(p == p0);
    EXPECT_EQ(p, p0);
}

TEST(Homogeneous, IdealPoints)
{
    using namespace xiaotu::math;

    HomoLine2D l0(0.1, 1.0, 0.0);
    HomoLine2D l1(0.1, 1.0, 1.0);
    EXPECT_TRUE(HomoPoint2D(l1.Intersection(l0)).IsInfinity());

    HomoLine2D l2 = HomoLine2D::Infinity();
    EXPECT_TRUE(l2.IsInfinity());
    EXPECT_TRUE(HomoPoint2D(l2.Intersection(l0)).IsInfinity());
    EXPECT_TRUE(HomoPoint2D(l2.Intersection(l1)).IsInfinity());

    HomoPoint2D p0 = l0.Intersection(l2);
    EXPECT_EQ(l0.a, -p0.y);
    EXPECT_EQ(l0.b, p0.x);
    EXPECT_EQ(0.0, p0.k);
}

TEST(Homogeneous, Conic2D)
{
    using namespace xiaotu::math;

    HomoConic2D conic0(1, 2, 3, 4, 5, 6);
    EXPECT_EQ(1.0, conic0.a());
    EXPECT_EQ(2.0, conic0.b());
    EXPECT_EQ(3.0, conic0.c());
    EXPECT_EQ(4.0, conic0.d());
    EXPECT_EQ(5.0, conic0.e());
    EXPECT_EQ(6.0, conic0.f());

    std::vector<double> datas = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };
    HomoConic2D conic1(datas.data());
    EXPECT_EQ(0.1, conic1.a());
    EXPECT_EQ(0.2, conic1.b());
    EXPECT_EQ(0.3, conic1.c());
    EXPECT_EQ(0.4, conic1.d());
    EXPECT_EQ(0.5, conic1.e());
    EXPECT_EQ(0.6, conic1.f());

    conic0.SetValue(1.0, 0, 0, 0, -1, 0);
    HomoPoint2D p(0, 0, 1);
    EXPECT_TRUE(PointOn(p, conic0));

    HomoPoint2D ps[5];
    ps[0].SetValue(0, -1, 1);
    ps[1].SetValue(1, 0, 1);
    ps[2].SetValue(-1, 0, 1);
    ps[3].SetValue(2, 3, 1);
    ps[4].SetValue(-2, 3, 1);

    conic1.From5Points(ps);
    EXPECT_TRUE(PointOn(ps[0], conic1));
    EXPECT_TRUE(PointOn(ps[1], conic1));
    EXPECT_TRUE(PointOn(ps[2], conic1));
    EXPECT_TRUE(PointOn(ps[3], conic1));
    EXPECT_TRUE(PointOn(ps[4], conic1));

    HomoLine2D l = conic1.Tangent(ps[0]);
    EXPECT_EQ(l, HomoLine2D(0, 1, 1));

    Eigen::Matrix3d H;
    H << 1, 0, 1,
         0, 1, 3,
         0, 0, 1;

    p = H * ps[0];
    conic1.Transform(H);
    EXPECT_TRUE(PointOn(p, conic1));
}


TEST(Homogeneous, Point3D)
{
    using namespace xiaotu::math;

    HomoPoint3D p0(0.1, 1.0, 2.0, 2.0);
    EXPECT_EQ(&p0[0], &p0.x);
    EXPECT_EQ(&p0[1], &p0.y);
    EXPECT_EQ(&p0[2], &p0.z);
    EXPECT_EQ(&p0[3], &p0.k);

    HomoPoint3D & p1 = p0.Normalize();
    EXPECT_EQ(1.0, p0.k);
    EXPECT_EQ(0.05, p0.x);
    EXPECT_EQ(0.5, p0.y);
    EXPECT_EQ(1.0, p0.z);
    EXPECT_EQ(1.0, p1.k);
    EXPECT_EQ(&p1, &p0);

    HomoPoint3D p2(0.1, 1.0, 2.0, 2.0);
    EXPECT_TRUE(p2 == p1);

    //p1.SetValue(1.0, 0.0, 1.0);
    //p2.SetValue(0.0, 1.0, 1.0);
    //HomoLine2D l1 = p0.JoinLine(p1);
    //HomoLine2D l2 = p1.JoinLine(p0);
    //EXPECT_EQ(l1, l2);

    HomoPoint3D p3 = HomoPoint3D::Infinity;
    EXPECT_TRUE(p3.IsInfinity());


}


TEST(Homogeneous, Plane3D)
{
    using namespace xiaotu::math;

    HomoPlane3D pi(0.1, 1.0, 1.0, 0.0);
    EXPECT_EQ(&pi[0], &pi.a);
    EXPECT_EQ(&pi[1], &pi.b);
    EXPECT_EQ(&pi[2], &pi.c);
    EXPECT_EQ(&pi[3], &pi.d);

    HomoPoint3D po(0.0, 0.0, 0.0, 1.0);
    EXPECT_TRUE(PointOn(po, pi));

    HomoPoint3D po1(1.0,  2.0, 0.0, 1.0);
    HomoPoint3D po2(2.0,  2.0, 0.0, 1.0);
    HomoPoint3D po3(2.0, -2.0, 0.0, 1.0);
    HomoPlane3D pi1(po1, po2, po3);

    EXPECT_TRUE(PointOn(po, pi1));
    po.SetValue(-9.0, 11.0, 0.0, 1.0);
    EXPECT_TRUE(PointOn(po, pi1));
    EXPECT_TRUE(PointOn(po1, pi1));
    EXPECT_TRUE(PointOn(po2, pi1));
    EXPECT_TRUE(PointOn(po3, pi1));

    pi1.SetValue(0.0, 0.0, 1.0, 0.0);
    HomoPlane3D pi2(1.0, 0.0, 0.0, 0.0);
    HomoPlane3D pi3(0.0, 1.0, 0.0, 0.0);

    po1 = Intersection(pi1, pi2, pi3);
    EXPECT_EQ(po1, HomoPoint3D(0, 0, 0, 1));
}

TEST(Homogeneous, Line3D)
{
    using namespace xiaotu::math;

    HomoLine3D line;
    line << 11, 12, 13, 14,
            21, 22, 23, 24,
            31, 32, 33, 34,
            41, 42, 43, 44;
    EXPECT_EQ(12, line.l12);
    EXPECT_EQ(13, line.l13);
    EXPECT_EQ(14, line.l14);
    EXPECT_EQ(23, line.l23);
    EXPECT_EQ(42, line.l42);
    EXPECT_EQ(34, line.l34);

    HomoPoint3D po1(0.0, 0.0, 0.0, 1.0);
    HomoPoint3D po2(1.0, 0.0, 0.0, 0.0);
    HomoLine3D line1(po1, po2);
    EXPECT_EQ(-1, line1.l14);
    EXPECT_TRUE(line1.IsValid());
    

    po1 << 1.0, 0.0, 0.0, 1.0;
    po2 << 2.0, 0.0, 1.0, 1.0;
    HomoLine3D line2(po1, po2);
    EXPECT_TRUE(line2.IsValid());


    po1 << 0.0,  2.0, 1.0, 1.0;
    po2 << 0.0,  2.0, 2.0, 1.0;
    HomoPoint3D po3(0.0, -2.0, 2.0, 1.0);
    HomoPlane3D pi1(po1, po2, po3);
    po3 = line1.Intersection(pi1);
    EXPECT_EQ(po3, HomoPoint3D(0, 0, 0, 1));

    po1 << 3.1415926,  2.0, 0.0, 1.0;
    po2 << 2.0, 99.0, 0.0, 1.0;
    pi1 = line1.JoinPlane(po2);
    EXPECT_TRUE(PointOn(po1, pi1));
}





