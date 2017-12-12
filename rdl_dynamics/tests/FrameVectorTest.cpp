#include <gtest/gtest.h>
#include "rdl_dynamics/FrameVector.hpp"
#include "UnitTestUtils.hpp"

using namespace RobotDynamics;

class FrameVectorTest : public ::testing::Test
{
protected:

    virtual void SetUp()
    {
        std::srand( time(NULL) );
    }

    virtual void TearDown()
    {
    }

    std::shared_ptr<ReferenceFrame> root = ReferenceFrame::createARootFrame("root1");
    std::shared_ptr<unit_test_utils::RandomUnchangingFrame> frame1 = unit_test_utils::RandomUnchangingFrame::create("frame1", root.get(),1);
    std::shared_ptr<unit_test_utils::RandomUnchangingFrame> frame2 = unit_test_utils::RandomUnchangingFrame::create("frame2", frame1.get(),2);

    int nTests = 1000;

private:

};

TEST_F(FrameVectorTest, operators)
{
    FrameVector v1(root.get(), 1.,2.,3.);
    Vector3d v2(2.,3.,4.);
    FrameVector v3;

    v3 = v1+v2;

    EXPECT_STREQ(v3.getReferenceFrame()->getName().c_str(),"root1");
    EXPECT_EQ(v3.x(),3.);
    EXPECT_EQ(v3.y(),5.);
    EXPECT_EQ(v3.z(),7.);

    EXPECT_STREQ(v1.getReferenceFrame()->getName().c_str(),"root1");
    EXPECT_EQ(v1.x(),1.);
    EXPECT_EQ(v1.y(),2.);
    EXPECT_EQ(v1.z(),3.);

    v3 = v1-v2;
    EXPECT_STREQ(v3.getReferenceFrame()->getName().c_str(),"root1");
    EXPECT_EQ(v3.x(),-1.);
    EXPECT_EQ(v3.y(),-1.);
    EXPECT_EQ(v3.z(),-1.);
}

TEST_F(FrameVectorTest, testConstructors)
{
    FrameVector frameVector(root.get(),1.,2.,3.);

    Eigen::Vector3d vector = frameVector;
    EXPECT_TRUE(vector(0)==1.);
    EXPECT_TRUE(vector(1)==2.);
    EXPECT_TRUE(vector(2)==3.);

    Eigen::Vector3d vector2(3.,2.,1.);
    FrameVector frameVector2(root.get(),vector2);

    Eigen::Vector3d vectorCheck = frameVector2;
    EXPECT_TRUE(vectorCheck(0)==3.);
    EXPECT_TRUE(vectorCheck(1)==2.);
    EXPECT_TRUE(vectorCheck(2)==1.);
}

TEST_F(FrameVectorTest, testSetIncludingFrame)
{
    FrameVector frameVector(root.get(),1.,2.,3.);

    frameVector.setIncludingFrame(0.1,0.2,0.3,frame1.get());

    EXPECT_STREQ(frameVector.getReferenceFrame()->getName().c_str(),"frame1");
    EXPECT_EQ(frameVector.x(),0.1);
    EXPECT_EQ(frameVector.y(),0.2);
    EXPECT_EQ(frameVector.z(),0.3);

    frameVector.setIncludingFrame(Vector3d(0.3,0.2,0.1),frame2.get());

    EXPECT_STREQ(frameVector.getReferenceFrame()->getName().c_str(),"frame2");
    EXPECT_EQ(frameVector.x(),0.3);
    EXPECT_EQ(frameVector.y(),0.2);
    EXPECT_EQ(frameVector.z(),0.1);

    try
    {
        frameVector.setIncludingFrame(0.1,0.2,0.3,nullptr);
    }
    catch(RobotDynamics::ReferenceFrameException &e)
    {
        EXPECT_STREQ(e.what(),"Reference frame cannot be nullptr!");
    }

    try
    {
        frameVector.setIncludingFrame(Vector3d(-0.1,-0.2,-0.3),nullptr);
    }
    catch(RobotDynamics::ReferenceFrameException &e)
    {
        EXPECT_STREQ(e.what(),"Reference frame cannot be nullptr!");
    }
}

TEST_F(FrameVectorTest, testDot)
{
    FrameVector frameVector1(frame1.get(),-1.,2.,-3.);
    FrameVector frameVector2(frame2.get(),1.,2.,3.);
    FrameVector frameVector3(frame1.get(),4.,5.,-6.);

    try
    {
        frameVector1.dot(frameVector2);
    }
    catch(ReferenceFrameException &e)
    {
        EXPECT_STREQ(e.what(),"Reference frames do not match!");
    }

    double value = frameVector1.dot(frameVector3);

    EXPECT_TRUE(value=24);
}

TEST_F(FrameVectorTest, testCross)
{
    Eigen::Vector3d result;
    Eigen::Vector3d expectedResult;

    for(int i = 0; i<nTests; i++)
    {
        Eigen::Vector3d v1(unit_test_utils::getRandomNumber<double>(),unit_test_utils::getRandomNumber<double>(),unit_test_utils::getRandomNumber<double>());
        Eigen::Vector3d v2(unit_test_utils::getRandomNumber<double>(),unit_test_utils::getRandomNumber<double>(),unit_test_utils::getRandomNumber<double>());
        Eigen::Vector3d v3(unit_test_utils::getRandomNumber<double>(),unit_test_utils::getRandomNumber<double>(),unit_test_utils::getRandomNumber<double>());

        FrameVector frameVector1(frame1.get(),v1);
        FrameVector frameVector2(frame2.get(),v2);
        FrameVector frameVector3(frame1.get(),v3);

        try
        {
            frameVector1.cross(frameVector2);
        }
        catch(ReferenceFrameException &e)
        {
            EXPECT_STREQ(e.what(),"Reference frames do not match!");
        }

        Eigen::Vector3d result = frameVector1.cross(frameVector3);
        Eigen::Vector3d expectedResult = v1.cross(v3);
        EXPECT_TRUE(result(0)==expectedResult(0));
        EXPECT_TRUE(result(1)==expectedResult(1));
        EXPECT_TRUE(result(2)==expectedResult(2));
    }
}

TEST_F(FrameVectorTest, testAngleBetweenVectors)
{
    FrameVector frameVector1(frame1.get(),2.,3.,1.);
    FrameVector frameVector2(frame1.get(),4.,1.,2.);


    double angle = frameVector1.getAngleBetweenVectors(frameVector2);
    double expectedResult = acos(13/(sqrt(14)*sqrt(21)));
    EXPECT_TRUE(angle==expectedResult);

    frameVector1.set(0.,0.,0.);

    EXPECT_NEAR(frameVector1.getAngleBetweenVectors(frameVector2),M_PI,1e-12);
}

TEST_F(FrameVectorTest, testChangeFrame)
{
    RobotDynamics::Math::SpatialTransform transform1;

    transform1.E = RobotDynamics::Math::Xrotx(M_PI_2).E;
    transform1.r = RobotDynamics::Math::Xtrans(RobotDynamics::Math::Vector3d(5.,0.,0.)).r;

    std::shared_ptr<ReferenceFrame> frameA(new unit_test_utils::RandomUnchangingFrame("A", root.get(), transform1,3));

    transform1.E = RobotDynamics::Math::Xroty(M_PI_2).E;
    transform1.r = RobotDynamics::Math::Xtrans(RobotDynamics::Math::Vector3d(5.,0.,0.)).r;

    std::shared_ptr<ReferenceFrame> frameB(new unit_test_utils::RandomUnchangingFrame("B", frameA.get(), transform1,4));

    transform1.E = RobotDynamics::Math::Xrotz(M_PI_2).E;
    transform1.r = RobotDynamics::Math::Xtrans(RobotDynamics::Math::Vector3d(5.,0.,0.)).r;

    std::shared_ptr<ReferenceFrame> frameC(new unit_test_utils::RandomUnchangingFrame("C", frameB.get(), transform1,5));

    double x = 3.0;
    double y = 1.0;
    double z = -9.0;

    FrameVector frameVector(frameC.get(), x, y, z);
    frameVector.changeFrame(frameB.get());

    EXPECT_NEAR(frameVector.x(),-1,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(frameVector.y(),3,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(frameVector.z(),-9,unit_test_utils::TEST_PREC);

    frameVector.changeFrame(frameA.get());

    EXPECT_NEAR(frameVector.x(),-9,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(frameVector.y(),3,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(frameVector.z(),1,unit_test_utils::TEST_PREC);
}

TEST_F(FrameVectorTest, testChangeFrameAndCopy)
{
    RobotDynamics::Math::SpatialTransform transform1;

    transform1.E = RobotDynamics::Math::Xrotx(M_PI_2).E;
    transform1.r = RobotDynamics::Math::Xtrans(RobotDynamics::Math::Vector3d(5.,0.,0.)).r;

    std::shared_ptr<ReferenceFrame> frameA(new unit_test_utils::RandomUnchangingFrame("A", root.get(), transform1,3));

    transform1.E = RobotDynamics::Math::Xroty(M_PI_2).E;
    transform1.r = RobotDynamics::Math::Xtrans(RobotDynamics::Math::Vector3d(5.,0.,0.)).r;

    std::shared_ptr<ReferenceFrame> frameB(new unit_test_utils::RandomUnchangingFrame("B", frameA.get(), transform1,4));

    transform1.E = RobotDynamics::Math::Xrotz(M_PI_2).E;
    transform1.r = RobotDynamics::Math::Xtrans(RobotDynamics::Math::Vector3d(5.,0.,0.)).r;

    std::shared_ptr<ReferenceFrame> frameC(new unit_test_utils::RandomUnchangingFrame("C", frameB.get(), transform1,5));

    double x = 3.0;
    double y = 1.0;
    double z = -9.0;

    FrameVector v;
    FrameVector frameVector(frameC.get(), x, y, z);
    v = frameVector.changeFrameAndCopy(frameB);

    EXPECT_NEAR(v.x(),-1,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(v.y(),3,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(v.z(),-9,unit_test_utils::TEST_PREC);

    FrameVector v2;
    v.changeFrameAndCopy(frameA,v2);

    EXPECT_NEAR(v2.x(),-9,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(v2.y(),3,unit_test_utils::TEST_PREC);
    EXPECT_NEAR(v2.z(),1,unit_test_utils::TEST_PREC);
}

TEST_F(FrameVectorTest, testVectorLength)
{
    Eigen::Vector3d result;
    Eigen::Vector3d expectedResult;

    for(int i = 0; i<nTests; i++)
    {
        Eigen::Vector3d v1(unit_test_utils::getRandomNumber<double>(),unit_test_utils::getRandomNumber<double>(),unit_test_utils::getRandomNumber<double>());

        FrameVector frameVector1(frame1.get(),v1);

        double result = frameVector1.norm();
        double expectedResult = v1.norm();
        EXPECT_TRUE(result==expectedResult);
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    return RUN_ALL_TESTS();
}