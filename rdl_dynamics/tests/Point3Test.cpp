
#include "UnitTestUtils.hpp"
#include <gtest/gtest.h>

using namespace RobotDynamics;
using namespace RobotDynamics::Math;

class Point3Test : public ::testing::Test
{
protected:

    virtual void SetUp()
    {
        std::srand(time(NULL));
    }

    virtual void TearDown()
    {

    }

    int nTests = 1000;
};

TEST_F(Point3Test, testConstCorrectness)
{
    RobotDynamics::Math::Point3d p(1.1,2.2,3.3);

    p.x() = 10.2;

    EXPECT_EQ(10.2,p.x());

    double& Y = p.y();

    Y = 8.2;

    EXPECT_EQ(p.y(),8.2);

    EXPECT_EQ(p.z(),3.3);

    double* ptr = p.data();

    EXPECT_EQ(ptr[0],p.x());
    EXPECT_EQ(ptr[1],p.y());
    EXPECT_EQ(ptr[2],p.z());
}

TEST_F(Point3Test, testCross)
{
    Point3d p(1.,0.,0.);
    Vector3d v_p(1.,0.,0.);
    Vector3d v(0.,1.,0.);


    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(v_p.cross(v),p.cross(v),unit_test_utils::TEST_PREC));
}

TEST_F(Point3Test, testDistanceSquared)
{
    for (int i = 0; i < nTests; i++)
    {
        double x = unit_test_utils::getRandomNumber<double>();
        double y = unit_test_utils::getRandomNumber<double>();
        double z = unit_test_utils::getRandomNumber<double>();

        Point3d point(x, y, z);

        double x2 = unit_test_utils::getRandomNumber<double>();
        double y2 = unit_test_utils::getRandomNumber<double>();
        double z2 = unit_test_utils::getRandomNumber<double>();

        Point3d point2(x2, y2, z2);

        double distanceSquared = pow((x2 - x), 2) + pow((y2 - y), 2) + pow((z2 - z), 2);

        EXPECT_NEAR(distanceSquared,point.distanceSquared(point2),unit_test_utils::TEST_PREC);
    }

    for (int i = 0; i < nTests; i++)
    {
        double x = unit_test_utils::getRandomNumber<double>();
        double y = unit_test_utils::getRandomNumber<double>();
        double z = unit_test_utils::getRandomNumber<double>();

        Point3d point(x, y, z);

        double x2 = unit_test_utils::getRandomNumber<double>();
        double y2 = unit_test_utils::getRandomNumber<double>();
        double z2 = unit_test_utils::getRandomNumber<double>();

        Point3d point2(x2, y2, z2);

        double dx,dy,dz;
        dx=x-x2;
        dy=y-y2;
        dz=z-z2;

        double distanceSquared = dx*dx+dy*dy+dz*dz;

        EXPECT_NEAR(distanceSquared,point.distanceSquared(point2),unit_test_utils::TEST_PREC);
    }
}

TEST_F(Point3Test, testDistance)
{
    for (int i = 0; i < nTests; i++)
    {
        double x = unit_test_utils::getRandomNumber<double>();
        double y = unit_test_utils::getRandomNumber<double>();
        double z = unit_test_utils::getRandomNumber<double>();

        Point3d point(x, y, z);

        double x2 = unit_test_utils::getRandomNumber<double>();
        double y2 = unit_test_utils::getRandomNumber<double>();
        double z2 = unit_test_utils::getRandomNumber<double>();

        Point3d point2(x2, y2, z2);

        double distance = sqrt(pow((x2 - x), 2) + pow((y2 - y), 2) + pow((z2 - z), 2));

        EXPECT_NEAR(distance,point.distance(point2),unit_test_utils::TEST_PREC);
    }

    for (int i = 0; i < nTests; i++)
    {
        double x = unit_test_utils::getRandomNumber<double>();
        double y = unit_test_utils::getRandomNumber<double>();
        double z = unit_test_utils::getRandomNumber<double>();

        Point3d point(x, y, z);

        double x2 = unit_test_utils::getRandomNumber<double>();
        double y2 = unit_test_utils::getRandomNumber<double>();
        double z2 = unit_test_utils::getRandomNumber<double>();

        Point3d point2(x2, y2, z2);

        double dx,dy,dz;
        dx=x-x2;
        dy=y-y2;
        dz=z-z2;

        double distance = sqrt(dx*dx+dy*dy+dz*dz);

        EXPECT_NEAR(distance,point.distance(point2),unit_test_utils::TEST_PREC);
    }
}

TEST_F(Point3Test, testDistanceL1)
{
    for (int i = 0; i < nTests; i++)
    {
        double x = unit_test_utils::getRandomNumber<double>();
        double y = unit_test_utils::getRandomNumber<double>();
        double z = unit_test_utils::getRandomNumber<double>();

        Point3d point(x, y, z);

        double x2 = unit_test_utils::getRandomNumber<double>();
        double y2 = unit_test_utils::getRandomNumber<double>();
        double z2 = unit_test_utils::getRandomNumber<double>();

        Point3d point2(x2, y2, z2);

        double distanceL1 = fabs(x2 - x) + fabs(y2 - y) + fabs(z2 - z);

        EXPECT_NEAR(distanceL1,point.distanceL1(point2),unit_test_utils::TEST_PREC);
    }

    for (int i = 0; i < nTests; i++)
    {
        double x = unit_test_utils::getRandomNumber<double>();
        double y = unit_test_utils::getRandomNumber<double>();
        double z = unit_test_utils::getRandomNumber<double>();

        Point3d point(x, y, z);

        double x2 = unit_test_utils::getRandomNumber<double>();
        double y2 = unit_test_utils::getRandomNumber<double>();
        double z2 = unit_test_utils::getRandomNumber<double>();

        Point3d point2(x2, y2, z2);

        double distanceL1 = fabs(x2 - x) + fabs(y2 - y) + fabs(z2 - z);

        EXPECT_NEAR(distanceL1,point.distanceL1(point2),unit_test_utils::TEST_PREC);
    }
}

TEST_F(Point3Test, testDistanceLInf)
{
    for (int i = 0; i < nTests; i++)
    {
        double x = unit_test_utils::getRandomNumber<double>();
        double y = unit_test_utils::getRandomNumber<double>();
        double z = unit_test_utils::getRandomNumber<double>();

        Point3d point(x, y, z);

        double x2 = unit_test_utils::getRandomNumber<double>();
        double y2 = unit_test_utils::getRandomNumber<double>();
        double z2 = unit_test_utils::getRandomNumber<double>();

        Point3d point2(x2, y2, z2);

        double distanceLInf = std::max(fabs(x2 - x), fabs(y2 - y));
        distanceLInf = std::max(distanceLInf, fabs(z2 - z));

        EXPECT_NEAR(distanceLInf,point.distanceLinf(point2),unit_test_utils::TEST_PREC);
    }

    for (int i = 0; i < nTests; i++)
    {
        double x = unit_test_utils::getRandomNumber<double>();
        double y = unit_test_utils::getRandomNumber<double>();
        double z = unit_test_utils::getRandomNumber<double>();

        Point3d point(x, y, z);

        double x2 = unit_test_utils::getRandomNumber<double>();
        double y2 = unit_test_utils::getRandomNumber<double>();
        double z2 = unit_test_utils::getRandomNumber<double>();

        Point3d point2(x2, y2, z2);

        double distanceLInf = std::max(fabs(x2 - x), fabs(y2 - y));
        distanceLInf = std::max(distanceLInf, fabs(z2 - z));

        EXPECT_NEAR(distanceLInf,point.distanceLinf(point2),unit_test_utils::TEST_PREC);
    }
}

TEST_F(Point3Test, testAdd1)
{
    for (int i = 0; i < 1000; i++)
    {
        double array1[3] = {unit_test_utils::getRandomNumber<double>(), unit_test_utils::getRandomNumber<double>(), unit_test_utils::getRandomNumber<double>()};
        double array2[3] = {unit_test_utils::getRandomNumber<double>(), unit_test_utils::getRandomNumber<double>(), unit_test_utils::getRandomNumber<double>()};

        Point3d point1(array1);
        Point3d point2(array2);

        Point3d point3(array1[0] + array2[0], array1[1] + array2[1], array1[2] + array2[2]);

        point1 += point2;

        EXPECT_TRUE(point1.epsilonEquals(point3, unit_test_utils::TEST_PREC));
    }

    for (int i = 0; i < 1000; i++)
    {
        double array1[3] = {unit_test_utils::getRandomNumber<double>(), unit_test_utils::getRandomNumber<double>(), unit_test_utils::getRandomNumber<double>()};
        double array2[3] = {unit_test_utils::getRandomNumber<double>(), unit_test_utils::getRandomNumber<double>(), unit_test_utils::getRandomNumber<double>()};

        Point3d point1(array1);
        Point3d point2(array2);

        Point3d point3(array1[0] + array2[0], array1[1] + array2[1], array1[2] + array2[2]);

        point1 += point2;

        EXPECT_TRUE(point1.epsilonEquals(point3, unit_test_utils::TEST_PREC));
    }
}

TEST_F(Point3Test, testAdd2)
{
    Point3d point1,point2;
    Point3d point3,point4;

    for (int i = 0; i < 1000; i++)
    {
        point1.set(unit_test_utils::getRandomNumber<double>(), unit_test_utils::getRandomNumber<double>(), unit_test_utils::getRandomNumber<double>());

        std::vector<double> vector(3);
        vector[0] = unit_test_utils::getRandomNumber<double>();
        vector[1] = unit_test_utils::getRandomNumber<double>();
        vector[2] = unit_test_utils::getRandomNumber<double>();

        point2.set(vector);

        Point3d tmpPoint;
        tmpPoint.set(vector[0] + point1.x(), vector[1] + point1.y(), vector[2] + point1.z());

        point2 = point2 + point1;

        EXPECT_TRUE(point2.epsilonEquals(tmpPoint, unit_test_utils::TEST_PREC));
    }

    for (int i = 0; i < 1000; i++)
    {
        point3.set(unit_test_utils::getRandomNumber<double>(), unit_test_utils::getRandomNumber<double>(), unit_test_utils::getRandomNumber<double>());

        std::vector<double> vector(3);
        vector[0] = unit_test_utils::getRandomNumber<double>();
        vector[1] = unit_test_utils::getRandomNumber<double>();
        vector[2] = unit_test_utils::getRandomNumber<double>();

        point4.set(vector);

        Point3d tmpPoint;
        tmpPoint.set(vector[0] + point3.x(), vector[1] + point3.y(), vector[2] + point3.z());

        point4 = point4 + point3;

        EXPECT_TRUE(point4.epsilonEquals(tmpPoint, unit_test_utils::TEST_PREC));
    }
}

TEST_F(Point3Test, testSubtract1)
{
    for (int i = 0; i < 1000; i++)
    {
        Point3d point1 = unit_test_utils::getRandomPoint3();
        Point3d point2 = unit_test_utils::getRandomPoint3();
        Point3d point3(point1.x() - point2.x(), point1.y() - point2.y(), point1.z() - point2.z());
        point1 -= point2;

        EXPECT_TRUE(point3.epsilonEquals(point1, 1e-12));
    }

    for (int i = 0; i < 1000; i++)
    {
        Point3d point1 = unit_test_utils::getRandomPoint3();
        Point3d point2 = unit_test_utils::getRandomPoint3();
        Point3d point3(point1.x() - point2.x(), point1.y() - point2.y(), point1.z() - point2.z());
        point1 -= point2;

        EXPECT_TRUE(point3.epsilonEquals(point1, unit_test_utils::TEST_PREC));
    }
}

TEST_F(Point3Test, testScale1)
{
    for (int i = 0; i < 1000; i++)
    {
        Point3d point1 = unit_test_utils::getRandomPoint3();
        Point3d point2 = point1;
        double scale = rand() % 100 - 50;

        point1 *= scale;

        EXPECT_NEAR(point1.x(), point2.x() * scale,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point1.y(), point2.y() * scale,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point1.z(), point2.z() * scale,unit_test_utils::TEST_PREC);
    }
}

TEST_F(Point3Test, testAbsoluteValue1)
{
    for (int i = 0; i < 1000; i++)
    {
        Point3d point1 = unit_test_utils::getRandomPoint3();
        Point3d point2 = point1;
        point1.absoluteValue();

        EXPECT_NEAR(point1.x(),fabs(point2.x()),unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point1.y(),fabs(point2.y()),unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point1.z(),fabs(point2.z()),unit_test_utils::TEST_PREC);
    }

    for (int i = 0; i < 1000; i++)
    {
        Point3d point1 = unit_test_utils::getRandomPoint3();
        Point3d point2 = point1;
        point1.absoluteValue();

        EXPECT_NEAR(point1.x(),fabs(point2.x()),unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point1.y(),fabs(point2.y()),unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point1.z(),fabs(point2.z()),unit_test_utils::TEST_PREC);
    }
}

TEST_F(Point3Test, testClampMinMax1)
{
    {
        Point3d point1(100., 200., 300.);
        Point3d point2 = point1;
        point2.clampMinMax(-100., -50.);

        EXPECT_NEAR(point2.x(),-50.,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.y(),-50.,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.z(),-50.,unit_test_utils::TEST_PREC);
    }

    {
        Point3d point1(100., 200., 300.);
        Point3d point2 = point1;
        point2.clampMinMax(-100., -50.);

        EXPECT_NEAR(point2.x(),-50.,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.y(),-50.,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.z(),-50.,unit_test_utils::TEST_PREC);
    }
}

TEST_F(Point3Test, testClampMinMax2)
{
    {
        Point3d point1(100, 200, 300);
        Point3d point2 = point1;
        point2.clampMinMax(-100, 200);

        EXPECT_NEAR(point2.x(),100,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.y(),200,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.z(),200,unit_test_utils::TEST_PREC);
    }

    {
        Point3d point1(100, 200, 300);
        Point3d point2 = point1;
        point2.clampMinMax(-100, 200);

        EXPECT_NEAR(point2.x(),100,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.y(),200,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.z(),200,unit_test_utils::TEST_PREC);
    }
}

TEST_F(Point3Test, testClampMinMax3)
{
    {
        Point3d point1(100, 200, 300);
        Point3d point2 = point1;
        point2.clampMinMax(201, 220);

        EXPECT_NEAR(point2.x(),201,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.y(),201,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.z(),220,unit_test_utils::TEST_PREC);
    }

    {
        Point3d point1(100, 200, 300);
        Point3d point2 = point1;
        point2.clampMinMax(201, 220);

        EXPECT_NEAR(point2.x(),201,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.y(),201,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.z(),220,unit_test_utils::TEST_PREC);
    }
}

TEST_F(Point3Test, testClampMin)
{
    {
        Point3d point1(100, 200, 300);
        Point3d point2 = point1;
        point2.clampMin(201);

        EXPECT_NEAR(point2.x(),201,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.y(),201,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.z(),300,unit_test_utils::TEST_PREC);
    }

    {
        Point3d point1(100, 200, 300);
        Point3d point2 = point1;
        point2.clampMin(201);

        EXPECT_NEAR(point2.x(), 201,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.y(), 201,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.z(), 300,unit_test_utils::TEST_PREC);
    }

    {
        Point3d point1(100, 200, 199);
        Point3d point2 = point1;
        point2.clampMin(201);

        EXPECT_NEAR(point2.x(), 201,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.y(), 201,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.z(), 201,unit_test_utils::TEST_PREC);
    }
}

TEST_F(Point3Test, testClampMax)
{
    {
        Point3d point1(100, 200, 300);
        Point3d point2 = point1;
        point2.clampMax(201);

        EXPECT_NEAR(point2.x(),100,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.y(),200,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.z(),201,unit_test_utils::TEST_PREC);
    }

    {
        Point3d point1(100, 200, 300);
        Point3d point2 = point1;
        point2.clampMax(201);

        EXPECT_NEAR(point2.x(),100,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.y(),200,unit_test_utils::TEST_PREC);
        EXPECT_NEAR(point2.z(),201,unit_test_utils::TEST_PREC);
    }
}

TEST_F(Point3Test, testEqualsEquals)
{
    {
        double x = unit_test_utils::getRandomNumber<double>();
        double y = unit_test_utils::getRandomNumber<double>();
        double z = unit_test_utils::getRandomNumber<double>();
        Point3d point1(x, y, z);
        Point3d point2 = point1;

        EXPECT_TRUE(point2==point1);
        EXPECT_FALSE(point2!=point1);

        point2.x()=point2.x()+0.1;
        EXPECT_FALSE(point2==point1);
    }

    {
        float x = unit_test_utils::getRandomNumber<float>();
        float y = unit_test_utils::getRandomNumber<float>();
        float z = unit_test_utils::getRandomNumber<float>();
        Point3d point1(x, y, z);
        Point3d point2 = point1;

        EXPECT_TRUE(point2==point1);
        EXPECT_FALSE(point2!=point1);
    }
}

TEST_F(Point3Test, testTemplateTypeInference)
{
    double x = unit_test_utils::getRandomNumber<double>();
    double y = unit_test_utils::getRandomNumber<double>();
    double z = unit_test_utils::getRandomNumber<double>();
    std::vector<double> v(3);
    v[0] = x;
    v[1] = y;
    v[2] = z;
    Point3d point1(v);
    Point3d point2(x,y,z);

    EXPECT_NEAR(point1.x(),point2.x(),unit_test_utils::TEST_PREC);
    EXPECT_NEAR(point1.y(),point2.y(),unit_test_utils::TEST_PREC);
    EXPECT_NEAR(point1.z(),point2.z(),unit_test_utils::TEST_PREC);

    Point3d point3 = unit_test_utils::getRandomPoint3();

    point1.set(point3);

    EXPECT_NEAR(point3.x(),point1.x(),unit_test_utils::TEST_PREC);
    EXPECT_NEAR(point3.y(),point1.y(),unit_test_utils::TEST_PREC);
    EXPECT_NEAR(point3.z(),point1.z(),unit_test_utils::TEST_PREC);

    Point3d point4(point3);

    EXPECT_NEAR(point3.x(),point4.x(),unit_test_utils::TEST_PREC);
    EXPECT_NEAR(point3.y(),point4.y(),unit_test_utils::TEST_PREC);
    EXPECT_NEAR(point3.z(),point4.z(),unit_test_utils::TEST_PREC);
}

TEST_F(Point3Test, SpatialTransformApplyOnPoint3)
{
    RobotDynamics::Math::Matrix3d m1,m2;

    double th = M_PI_2;
    RobotDynamics::Math::Matrix3d Mx = RobotDynamics::Math::Xrotx(th).E;
    RobotDynamics::Math::Matrix3d My = RobotDynamics::Math::Xroty(th).E;

    RobotDynamics::Math::Vector3d r1 = RobotDynamics::Math::Xtrans(RobotDynamics::Math::Vector3d(6.,5.,3.)).r;
    RobotDynamics::Math::Vector3d r2 = RobotDynamics::Math::Xtrans(RobotDynamics::Math::Vector3d(3.,-6.,4.)).r;

    Point3d p(3.,-1.,2.);

    RobotDynamics::Math::SpatialTransform X_1(Mx,r1);
    RobotDynamics::Math::SpatialTransform X_2(My,r2);

    p.transform(X_1);

    EXPECT_TRUE(p.epsilonEquals(Point3d(-3.,-1.,6.),unit_test_utils::TEST_PREC));

    p.transform(X_2);

    Point3d pTest(3.,-1.,2.);

    RobotDynamics::Math::SpatialTransform X_3 = X_2*X_1;

    pTest.transform(X_3);

    EXPECT_TRUE(pTest.epsilonEquals(Point3d(-2.,5.,-6.),unit_test_utils::TEST_PREC));
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    return RUN_ALL_TESTS();
}