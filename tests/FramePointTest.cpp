#include <gtest/gtest.h>
#include "rdl_dynamics/FramePoint.hpp"
#include "UnitTestUtils.hpp"

using namespace RobotDynamics;

class FramePointTest : public ::testing::Test
{
protected:

    virtual void SetUp()
    {
        std::srand( time(NULL) );
    }

    virtual void TearDown()
    {
    }

    std::shared_ptr<ReferenceFrame> root1 = ReferenceFrame::createARootFrame("root1");

    int nTests = 1000;

private:

};

TEST_F(FramePointTest, testChangeFrame)
{
    RobotDynamics::Math::SpatialTransform transform1;

    transform1.E = RobotDynamics::Math::Xrotx(M_PI_2).E;
    transform1.r.set(6.,5.,3.);

    std::shared_ptr<ReferenceFrame> frameA(new unit_test_utils::RandomUnchangingFrame("A", root1.get(), transform1,1));

    transform1.E = RobotDynamics::Math::Xroty(M_PI_2).E;
    transform1.r.set(3.,-6.,4.);

    std::shared_ptr<ReferenceFrame> frameB(new unit_test_utils::RandomUnchangingFrame("B", frameA.get(), transform1,2));

    transform1.E = RobotDynamics::Math::Xrotz(M_PI_2).E;
    transform1.r = RobotDynamics::Math::Vector3d(0.,-5.,0.);
    std::shared_ptr<ReferenceFrame> frameC(new unit_test_utils::RandomUnchangingFrame("C", frameB.get(), transform1,3));

    frameA->update();
    frameB->update();
    frameC->update();

    double x = 3.;
    double y = -1.;
    double z = 2.;

    FramePoint framePoint(root1.get(), x, y, z);

    EXPECT_NEAR(framePoint.x(),x,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(framePoint.y(),y,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(framePoint.z(),z,unit_test_utils::TEST_PREC);

    framePoint.changeFrame(root1.get());

    EXPECT_NEAR(framePoint.x(),x,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(framePoint.y(),y,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(framePoint.z(),z,unit_test_utils::TEST_PREC);

    framePoint.changeFrame(frameA.get());

    EXPECT_NEAR(framePoint.x(),-3.,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(framePoint.y(),-1.,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(framePoint.z(),6.,unit_test_utils::TEST_PREC);

    framePoint.changeFrame(frameB.get());

    EXPECT_NEAR(framePoint.x(),-2.,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(framePoint.y(),5.,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(framePoint.z(),-6.,unit_test_utils::TEST_PREC);

    framePoint.changeFrame(frameC.get());

    EXPECT_NEAR(framePoint.x(),10.,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(framePoint.y(),2.,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(framePoint.z(),-6.,unit_test_utils::TEST_PREC);

    framePoint.changeFrame(frameB.get());

    EXPECT_NEAR(framePoint.x(),-2.,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(framePoint.y(),5.,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(framePoint.z(),-6.,unit_test_utils::TEST_PREC);

    framePoint.changeFrame(frameA.get());

    EXPECT_NEAR(framePoint.x(),-3.,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(framePoint.y(),-1.,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(framePoint.z(),6.,unit_test_utils::TEST_PREC);

    framePoint.changeFrame(root1.get());

    EXPECT_NEAR(framePoint.x(),x,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(framePoint.y(),y,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(framePoint.z(),z,unit_test_utils::TEST_PREC);
}

TEST_F(FramePointTest, testChangeFrameAndCopy)
{
    RobotDynamics::Math::SpatialTransform transform1;

    transform1.E = RobotDynamics::Math::Xrotx(M_PI_2).E;
    transform1.r.set(6.,5.,3.);

    std::shared_ptr<ReferenceFrame> frameA(new unit_test_utils::RandomUnchangingFrame("A", root1.get(), transform1,1));

    transform1.E = RobotDynamics::Math::Xroty(M_PI_2).E;
    transform1.r.set(3.,-6.,4.);

    std::shared_ptr<ReferenceFrame> frameB(new unit_test_utils::RandomUnchangingFrame("B", frameA.get(), transform1,2));

    transform1.E = RobotDynamics::Math::Xrotz(M_PI_2).E;
    transform1.r = RobotDynamics::Math::Vector3d(0.,-5.,0.);
    std::shared_ptr<ReferenceFrame> frameC(new unit_test_utils::RandomUnchangingFrame("C", frameB.get(), transform1,3));

    frameA->update();
    frameB->update();
    frameC->update();

    double x = 3.;
    double y = -1.;
    double z = 2.;

    FramePoint framePoint(root1.get(), x, y, z);

    EXPECT_NEAR(framePoint.x(),x,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(framePoint.y(),y,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(framePoint.z(),z,unit_test_utils::TEST_PREC);

    FramePoint p1 = framePoint.changeFrameAndCopy(root1);

    EXPECT_STREQ(framePoint.getReferenceFrame()->getName().c_str(),root1->getName().c_str());
    EXPECT_NEAR(p1.x(),x,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(p1.y(),y,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(p1.z(),z,unit_test_utils::TEST_PREC);

    FramePoint p2 = p1;
    p1.changeFrameAndCopy(frameA,p2);

    EXPECT_NEAR(p2.x(),-3.,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(p2.y(),-1.,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(p2.z(),6.,unit_test_utils::TEST_PREC);
}

TEST_F(FramePointTest, testDistance)
{
    FramePoint framePoint1(root1.get(), 1, 2, 3);
    FramePoint framePoint2(root1.get(), -1, -2, -3);

    RobotDynamics::Math::SpatialTransform transform1;

    transform1.E = RobotDynamics::Math::Xrotx(M_PI_2).E;
    transform1.r.set(5.,0.,0.);

    std::shared_ptr<ReferenceFrame> frameA(new unit_test_utils::RandomUnchangingFrame("A", root1.get(), transform1,1));

    FramePoint framePoint3(frameA.get(),1.,2.,3.);

    EXPECT_TRUE(framePoint1.distance(framePoint2)==sqrt(36.+16.+4.));

    try
    {
        framePoint3.distance(framePoint2);
    }
    catch(RobotDynamics::ReferenceFrameException &e)
    {
        EXPECT_STREQ(e.what(),"Reference frames do not match!");
    }
}

TEST_F(FramePointTest, testSetIncludingFrameMismatchingFrames)
{
    FramePoint framePoint2;

    try
    {
        framePoint2.setIncludingFrame(0.1,0.2,0.3,nullptr);
    }
    catch(RobotDynamics::ReferenceFrameException &e)
    {
        EXPECT_STREQ(e.what(),"Reference frame is nullptr!");
    }
}

TEST_F(FramePointTest, testDistanceMixingTypes)
{
    FramePoint framePoint1(root1.get(), 1, 2, 3);
    FramePoint framePoint2(root1.get(), -1, -2, -3);

    RobotDynamics::Math::SpatialTransform transform1;

    transform1.E = RobotDynamics::Math::Xrotx(M_PI_2).E;
    transform1.r.set(5.,0.,0.);

    std::shared_ptr<ReferenceFrame> frameA(new unit_test_utils::RandomUnchangingFrame("A", root1.get(), transform1.inverse(),1));

    FramePoint framePoint3(frameA.get(),1.,2.,3.);

    EXPECT_TRUE(framePoint1.distance(framePoint2)==sqrt(36.+16.+4.));

    try
    {
        framePoint3.distance(framePoint2);
    }
    catch(RobotDynamics::ReferenceFrameException &e)
    {
        EXPECT_STREQ(e.what(),"Reference frames do not match!");
    }
}

TEST_F(FramePointTest, testDistanceL1)
{
    std::shared_ptr<FramePoint> framePoint1(new RobotDynamics::Math::FramePoint(root1.get(), 1., 2., 3.));
    std::shared_ptr<FramePoint> framePoint2(new RobotDynamics::Math::FramePoint(root1.get(), -1., -2., -3.));

    RobotDynamics::Math::SpatialTransform transform1;

    transform1 = RobotDynamics::Math::Xrotx(M_PI_2);
    transform1.r.set(5.,0.,0.);

    std::shared_ptr<ReferenceFrame> frameA(new unit_test_utils::RandomUnchangingFrame("A", root1.get(), transform1.inverse(),1));

    std::shared_ptr<FramePoint> framePoint3(new RobotDynamics::Math::FramePoint(frameA.get(),1.,2.,3.));

    EXPECT_TRUE(framePoint1->distanceL1(*framePoint2)==12.);

    try
    {
        framePoint3->distanceL1(*framePoint2);
    }
    catch(RobotDynamics::ReferenceFrameException &e)
    {
        EXPECT_STREQ(e.what(),"Reference frames do not match!");
    }
}

TEST_F(FramePointTest, testDistanceL1MixingTypes)
{
    std::shared_ptr<FramePoint> framePoint1(new RobotDynamics::Math::FramePoint(root1.get(), 1., 2., 3.));
    std::shared_ptr<FramePoint> framePoint2(new RobotDynamics::Math::FramePoint(root1.get(), -1., -2., -3.));

    RobotDynamics::Math::SpatialTransform transform1;

    transform1 = RobotDynamics::Math::Xrotx(M_PI_2);
    transform1.r.set(5.,0.,0.);

    std::shared_ptr<ReferenceFrame> frameA(new unit_test_utils::RandomUnchangingFrame("A", root1.get(), transform1.inverse(),1));

    std::shared_ptr<FramePoint> framePoint3(new RobotDynamics::Math::FramePoint(frameA.get(),1.,2.,3.));

    EXPECT_TRUE(framePoint1->distanceL1(*framePoint2)==12.);

    try
    {
        framePoint3->distanceL1(*framePoint2);
    }
    catch(RobotDynamics::ReferenceFrameException &e)
    {
        EXPECT_STREQ(e.what(),"Reference frames do not match!");
    }
}

TEST_F(FramePointTest, testDistanceLinf)
{
    std::shared_ptr<FramePoint> framePoint1(new RobotDynamics::Math::FramePoint(root1.get(), 1., 2., 3.));
    std::shared_ptr<FramePoint> framePoint2(new RobotDynamics::Math::FramePoint(root1.get(), -1., -2., -3.));

    RobotDynamics::Math::SpatialTransform transform1;

    transform1 = RobotDynamics::Math::Xrotx(M_PI_2);
    transform1.r.set(5.,0.,0.);

    std::shared_ptr<ReferenceFrame> frameA(new unit_test_utils::RandomUnchangingFrame("A", root1.get(), transform1.inverse(),1));

    std::shared_ptr<FramePoint> framePoint3(new RobotDynamics::Math::FramePoint(frameA.get(),1.,2.,3.));

    EXPECT_TRUE(framePoint1->distanceLinf(*framePoint2)==6.);

    try
    {
        framePoint3->distanceL1(*framePoint2);
    }
    catch(RobotDynamics::ReferenceFrameException &e)
    {
        EXPECT_STREQ(e.what(),"Reference frames do not match!");
    }
}

TEST_F(FramePointTest, testDistanceLinfMixingTypes)
{
    std::shared_ptr<FramePoint> framePoint1(new RobotDynamics::Math::FramePoint(root1.get(), 1., 2., 3.));
    std::shared_ptr<FramePoint> framePoint2(new RobotDynamics::Math::FramePoint(root1.get(), -1., -2., -3.));

    RobotDynamics::Math::SpatialTransform transform1;

    transform1 = RobotDynamics::Math::Xrotx(M_PI_2);
    transform1.r.set(5.,0.,0.);

    std::shared_ptr<ReferenceFrame> frameA(new unit_test_utils::RandomUnchangingFrame("A", root1.get(), transform1.inverse(),1));

    std::shared_ptr<FramePoint> framePoint3(new RobotDynamics::Math::FramePoint(frameA.get(),1.,2.,3.));

    EXPECT_TRUE(framePoint1->distanceLinf(*framePoint2)==6.);

    try
    {
        framePoint3->distanceL1(*framePoint2);
    }
    catch(RobotDynamics::ReferenceFrameException &e)
    {
        EXPECT_STREQ(e.what(),"Reference frames do not match!");
    }
}

TEST_F(FramePointTest, testVectorConstructor)
{
    std::vector<double> v(3);
    v[0] = 1.;
    v[1] = 2.;
    v[2] = 3.;

    std::shared_ptr<FramePoint> framePoint(new RobotDynamics::Math::FramePoint(root1.get(),v));

    EXPECT_EQ(framePoint->x(),v[0]);
    EXPECT_EQ(framePoint->y(),v[1]);
    EXPECT_EQ(framePoint->z(),v[2]);
}

TEST_F(FramePointTest, testAdd)
{
    std::vector<double> v(3);
    v[0] = 1.;
    v[1] = 2.;
    v[2] = 3.;

    FramePoint framePoint(root1.get(),v);

    framePoint*=3.;

    EXPECT_EQ(framePoint.x(),v[0]*3.);
    EXPECT_EQ(framePoint.y(),v[1]*3.);
    EXPECT_EQ(framePoint.z(),v[2]*3.);

    FramePoint p2(root1.get(),RobotDynamics::Math::Vector3d(-1.,-2.,-3.));

    FramePoint p3 = framePoint - p2;

    EXPECT_EQ(p3.x(),4.);
    EXPECT_EQ(p3.y(),8.);
    EXPECT_EQ(p3.z(),12.);

    FramePoint p4 = p3;
    p4.add(p2);

    EXPECT_EQ(framePoint.x(),p4.x());
    EXPECT_EQ(framePoint.y(),p4.y());
    EXPECT_EQ(framePoint.z(),p4.z());
}

TEST_F(FramePointTest, testOperatorOverloads)
{
    std::vector<double> v(3);
    v[0] = 1.;
    v[1] = 2.;
    v[2] = 3.;

    FramePoint framePoint(root1.get(),v);

    framePoint*=3.;

    EXPECT_EQ(framePoint.x(),v[0]*3.);
    EXPECT_EQ(framePoint.y(),v[1]*3.);
    EXPECT_EQ(framePoint.z(),v[2]*3.);

    FramePoint p2(root1.get(),RobotDynamics::Math::Vector3d(-1.,-2.,-3.));

    FramePoint p3 = framePoint - p2;

    EXPECT_EQ(p3.x(),4.);
    EXPECT_EQ(p3.y(),8.);
    EXPECT_EQ(p3.z(),12.);

    FramePoint p4 = p3+p2;

    EXPECT_EQ(framePoint.x(),p4.x());
    EXPECT_EQ(framePoint.y(),p4.y());
    EXPECT_EQ(framePoint.z(),p4.z());

    Vector3d v3d(13.,14.,15.);
    p4.set(1.,2.,3.);
    p4+=v3d;

    EXPECT_EQ(p4.x(),14.);
    EXPECT_EQ(p4.y(),16.);
    EXPECT_EQ(p4.z(),18.);

    p4-=v3d;

    EXPECT_EQ(p4.x(),1.);
    EXPECT_EQ(p4.y(),2.);
    EXPECT_EQ(p4.z(),3.);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    return RUN_ALL_TESTS();
}

