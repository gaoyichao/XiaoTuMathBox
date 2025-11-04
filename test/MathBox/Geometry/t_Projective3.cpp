#include <iostream>

#include <XiaoTuDataBox/Utils.hpp>
#include <XiaoTuMathBox/Geometry/Geometry.hpp>

#include <gtest/gtest.h>

TEST(Projective3, HomoPoint3)
{
    using namespace xiaotu::math;
    EXPECT_EQ(sizeof(double) * 4, sizeof(HomoPoint3<double>));

    HomoPoint3<double> p1(3.14159, 1.41421, 2.71828, 0.5);
    EXPECT_EQ(3.14159, p1.x());
    EXPECT_EQ(1.41421, p1.y());
    EXPECT_EQ(2.71828, p1.z());
    EXPECT_EQ(0.5, p1.k());
    EXPECT_EQ(&p1(0), &p1.x());
    EXPECT_EQ(&p1(1), &p1.y());
    EXPECT_EQ(&p1(2), &p1.z());
    EXPECT_EQ(&p1(3), &p1.k());

    {
        HomoPoint3<double> p2 = p1;
        EXPECT_TRUE(p2 == p1);
        EXPECT_EQ(p2, p1);

        p2 = p2 * 0.1;
        XTLog(std::cout) << p1 << std::endl;
        XTLog(std::cout) << p2 << std::endl;
        EXPECT_TRUE(p2 == p1);
        EXPECT_EQ(p2, p1);
    }

    {
        HomoPoint3<double> p2 = p1.Normalization();
        EXPECT_DOUBLE_EQ(1.0, p2.k());
        EXPECT_DOUBLE_EQ(0.5, p1.k());

        EXPECT_TRUE(p2 == p1);
        EXPECT_EQ(p2, p1);
    }

    {

        p1 = { 0.1, 1.0, 2.0, 3.0 };
        HomoPoint3<double> p2 = (1.0 / p1.k()) * p1;
        EXPECT_DOUBLE_EQ(1.0, p2.k());
        EXPECT_DOUBLE_EQ(3.0, p1.k());

        EXPECT_TRUE(p2 == p1);
        EXPECT_EQ(p2, p1);
    }

    {
        EXPECT_FALSE(p1.IsInfinity());

        p1 = { 0.1, 1.0, 10.0, 0.0 };
        EXPECT_TRUE(p1.IsInfinity());
        XTLog(std::cout) << p1 << std::endl;

        p1 << 0.1, 1.0, 2.0, 0.0;
        EXPECT_TRUE(p1.IsInfinity());
        XTLog(std::cout) << p1 << std::endl;
    }
}

TEST(Projective3, HomoPlane3)
{
    using namespace xiaotu::math;
    EXPECT_EQ(sizeof(double) * 4, sizeof(HomoPlane3<double>));

    HomoPlane3<double> l1(3.14159, 1.41421, 2.71828, 0.0);
    EXPECT_DOUBLE_EQ(3.14159, l1.a());
    EXPECT_DOUBLE_EQ(1.41421, l1.b());
    EXPECT_DOUBLE_EQ(2.71828, l1.c());
    EXPECT_EQ(&l1(0), &l1.a());
    EXPECT_EQ(&l1(1), &l1.b());
    EXPECT_EQ(&l1(2), &l1.c());
    EXPECT_EQ(&l1(3), &l1.d());

    {
        HomoPlane3<double> l2 = l1;
        EXPECT_TRUE(l2 == l1);
        EXPECT_EQ(l2, l1);

        l2 = l2 * 0.1;
        XTLog(std::cout) << l1 << std::endl;
        XTLog(std::cout) << l2 << std::endl;
        EXPECT_TRUE(l2 == l1);
        EXPECT_EQ(l2, l1);
    }


    {
        HomoPlane3<double> l2 = l1.Normalization();
        EXPECT_DOUBLE_EQ(1.0, l2.Norm());
        EXPECT_DOUBLE_EQ(2.71828, l1.c());

        EXPECT_TRUE(l2 == l1);
        EXPECT_EQ(l2, l1);
    }
    
    {
        l1.SetValue(0.1, 1.0, 2.0, 1.0);
        EXPECT_EQ(0.1, l1(0));
        EXPECT_EQ(1.0, l1(1));
        EXPECT_EQ(2.0, l1(2));
        EXPECT_EQ(1.0, l1(3));

        HomoPlane3<double> l2 = l1.Normalize();
        EXPECT_TRUE(l2 == l1);
        EXPECT_EQ(l2, l1);
    }

    {
        EXPECT_FALSE(l1.IsInfinity());
        l1 << 0.0, 0.0, 0.0, 1.0;
        EXPECT_TRUE(l1.IsInfinity());
    }
}

TEST(Projective3, HomoPoint3Plane3)
{
    using namespace xiaotu::math;

    HomoPoint3<double> po1(1.0,  2.0, 0.0, 1.0);
    HomoPoint3<double> po2(2.0,  2.0, 0.0, 1.0);
    HomoPoint3<double> po3(2.0, -2.0, 0.0, 1.0);
    
    HomoPlane3<double> pi1 = JoinSVD(po1, po2, po3);
    HomoPlane3<double> pi2 = Join(po1, po2, po3);
    XTLog(std::cout) << "pi1" << pi1 << std::endl;
    EXPECT_TRUE(pi1 == pi2);
    EXPECT_TRUE(OnPlane(po1, pi1));
    EXPECT_TRUE(OnPlane(po2, pi1));
    EXPECT_TRUE(OnPlane(po3, pi1));
    EXPECT_TRUE(OnPlane(po1, pi2));
    EXPECT_TRUE(OnPlane(po2, pi2));
    EXPECT_TRUE(OnPlane(po3, pi2));

    HomoPoint3<double> po(0.0, 0.0, 0.0, 1.0);
    EXPECT_TRUE(OnPlane(po, pi1));
    EXPECT_TRUE(OnPlane(po, pi2));
    po.SetValue(-9.0, 11.0, 0.0, 1.0);
    EXPECT_TRUE(OnPlane(po, pi1));
    EXPECT_TRUE(OnPlane(po, pi2));

                       pi1 = { 0.0, 0.0, 1.0, 0.0 };
                       pi2 = { 1.0, 0.0, 0.0, 0.0 };
    HomoPlane3<double> pi3 = { 0.0, 1.0, 0.0, 0.0 };

    po1 = MeetSVD(pi1, pi2, pi3);
    EXPECT_EQ(po1, HomoPoint3<double>(0, 0, 0, 1));
    XTLog(std::cout) << "po1" << po1 << std::endl;

    po2 = MeetSVD(pi1, pi2, pi3);
    EXPECT_EQ(po1, po2);
}

TEST(Projective3, HomoLine3)
{
    using namespace xiaotu::math;
    EXPECT_EQ(sizeof(double) * 4 * 4, sizeof(HomoLine3<double>));

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


TEST(Projective3, HomoLine3Point3Plane3)
{
     using namespace xiaotu::math;

     HomoLine3<double> line;
     HomoPoint3<double> po1({0.0, 0.0, 0.0, 1.0});
     HomoPoint3<double> po2({1.0, 0.0, 0.0, 0.0});
     HomoPoint3<double> po3({3.0, 0.0, 0.0, 1.0});
     line = Join(po1, po2);
     EXPECT_EQ(-1, line.l14());
     EXPECT_TRUE(line.IsValid());
     EXPECT_TRUE(OnLine(po3, line));

 
     HomoPlane3<double> pi1({0.0, 1.0, 0.0, 0.0});
     HomoPlane3<double> pi2({0.0, 0.0, 1.0, 0.0});
     HomoLine3<double> line2 = Meet(pi1, pi2);

     EXPECT_TRUE(AreCoplanar(line, line2));
     EXPECT_TRUE(OnLine(po1, line2));
     EXPECT_TRUE(OnLine(po2, line2));
     EXPECT_TRUE(OnPlane(line, pi1));
     EXPECT_TRUE(OnPlane(line, pi2));

     pi2 << 1.0, 0.0, 0.0, -1.0;
     po3 = Meet(line, pi2);
     EXPECT_EQ(po3, HomoPoint3<double>(1, 0, 0, 1));

     po3 << 0.0, 1.0, 0.0, 1.0;
     pi1 = Join(po3, line);
     EXPECT_EQ(pi1, HomoPlane3<double>(0, 0, 2, 0));
}

// TEST(Projective3, Projective3)
// {
//     using namespace xiaotu::math;
//     EXPECT_EQ(sizeof(Eigen::Matrix4d), sizeof(Projective3<double>));

//     Projective3<double> H;
//     H = Eigen::Matrix4d::Identity();
//     EXPECT_TRUE(H.isIdentity());

//     double rad = M_PI * 0.5;
//     double c = std::cos(rad);
//     double s = std::sin(rad);

//     H << c, -s, 0, 0,
//          s,  c, 0, 0,
//          0,  0, 1, 0,
//          0,  0, 0, 1;
//     HomoPoint3<double> p0(1, 0, 0, 1);
//     HomoPoint3<double> p1(1, 1, 0, 1);
//     HomoPoint3<double> p2(0, 1, 0, 1);
//     HomoPoint3<double> p3 = H.ApplyOn(p0);
//     EXPECT_EQ(p3, HomoPoint3<double>(0, 1, 0, 1));

//     HomoPlane3<double> pi = Join(p0, p1, p2);
//     HomoPlane3<double> _pi = H.ApplyOn(pi);
//     EXPECT_TRUE(OnPlane(p3, _pi));

//     HomoLine3<double> line;
//     p0 << 0.0, 0.0, 0.0, 1.0;
//     p1 << 1.0, 0.0, 0.0, 0.0;
//     p2 << 3.0, 0.0, 0.0, 1.0;
//     line = Join(p0, p1);
//     EXPECT_TRUE(line.IsValid());
//     EXPECT_TRUE(OnLine(p2, line));
 
//     HomoLine3<double> _line = H.ApplyOn(line);
//     p2 = H.ApplyOn(p2);
//     EXPECT_TRUE(_line.IsValid());
//     EXPECT_TRUE(OnLine(p2, _line));
// }

// TEST(Projective3, Properties)
// {
//     using namespace xiaotu::math;
//     EXPECT_EQ(sizeof(Eigen::Matrix4d), sizeof(Projective3<double>));

//     HomoPlane3<double> pi1(1.0, 0.0, 0.0, 0.0);
//     HomoPlane3<double> pi2(3.0, 0.0, 0.0, 0.0);
//     HomoLine3<double> line1 = Meet(pi1, pi2);

//     EXPECT_TRUE(OnPlane(line1, HomoPlane3<double>(0, 0, 0, 1)));


// }

