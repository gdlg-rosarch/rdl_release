

#include <gtest/gtest.h>
#include "rdl_dynamics/FrameOrientation.hpp"
#include "UnitTestUtils.hpp"

using namespace RobotDynamics;
using namespace RobotDynamics::Math;

class FrameOrientationTest : public ::testing::Test
{
public:
    FrameOrientationTest(){}

    void SetUp()
    {}

    void TearDown()
    {}
};

TEST_F(FrameOrientationTest, setters)
{
    ReferenceFrame f1("f1",ReferenceFrame::getWorldFrame().get(),SpatialTransform(),true,0);
    FrameOrientation orientation(&f1,Quaternion());

    Quaternion q(1.,2.,3.,4.);

    orientation.setOrientation(q);
    EXPECT_EQ(orientation.getOrientation().x(),1.);
    EXPECT_EQ(orientation.getOrientation().y(),2.);
    EXPECT_EQ(orientation.getOrientation().z(),3.);
    EXPECT_EQ(orientation.getOrientation().w(),4.);

    SpatialTransform X = Xrotz(0.1);
    Quaternion q2(Quaternion::fromMatrix(X.E));

    orientation.setOrientation(X.E);

    EXPECT_EQ(orientation.getOrientation().x(),Quaternion::fromMatrix(X.E).x());
    EXPECT_EQ(orientation.getOrientation().y(),Quaternion::fromMatrix(X.E).y());
    EXPECT_EQ(orientation.getOrientation().z(),Quaternion::fromMatrix(X.E).z());
    EXPECT_EQ(orientation.getOrientation().w(),Quaternion::fromMatrix(X.E).w());
}

TEST_F(FrameOrientationTest, changeFrame)
{
    SpatialTransform X1 = Xrotz(M_PI_2);
    ReferenceFrame f1("f1",ReferenceFrame::getWorldFrame().get(),X1,true,0);
    SpatialTransform X2 = Xrotz(-M_PI_2);
    ReferenceFrame f2("f2", &f1,X2,true,1);
    SpatialTransform X3 = Xrotz(-M_PI_2);
    ReferenceFrame f3("f3", &f2,X3,true,1);

    FrameOrientation orientation(&f2,Quaternion());
    orientation.changeFrame(&f1);

    Quaternion expected(Quaternion::fromMatrix(X2.E));
    EXPECT_TRUE(unit_test_utils::checkVector4dEpsilonClose(orientation.getOrientation(),expected,unit_test_utils::TEST_PREC));
    EXPECT_STREQ(orientation.getReferenceFrame()->getName().c_str(),"f1");

    orientation.changeFrame(ReferenceFrame::getWorldFrame());

    expected = Quaternion(0.,0.,0.,1.);
    EXPECT_TRUE(unit_test_utils::checkVector4dEpsilonClose(orientation.getOrientation(),expected,unit_test_utils::TEST_PREC));
    EXPECT_STREQ(orientation.getReferenceFrame()->getName().c_str(),"World");

    orientation.changeFrame(&f3);
    expected = Quaternion::fromMatrix(Xrotz(M_PI_2).E);
    EXPECT_TRUE(unit_test_utils::checkVector4dEpsilonClose(orientation.getOrientation(),expected,unit_test_utils::TEST_PREC));
    EXPECT_STREQ(orientation.getReferenceFrame()->getName().c_str(),"f3");
}

TEST_F(FrameOrientationTest, changeFrameAndCopy)
{
    SpatialTransform X1 = Xrotz(M_PI_2);
    std::shared_ptr<ReferenceFrame> f1(new ReferenceFrame("f1",ReferenceFrame::getWorldFrame().get(),X1,true,0));
    SpatialTransform X2 = Xrotz(-M_PI_2);
    std::shared_ptr<ReferenceFrame> f2(new ReferenceFrame("f2", f1.get(),X2,true,1));
    SpatialTransform X3 = Xrotz(-M_PI_2);
    std::shared_ptr<ReferenceFrame> f3(new ReferenceFrame("f3", f2.get(),X3,true,1));

    FrameOrientation orientation(f2.get(),Quaternion());
    FrameOrientation orientation_copy = orientation.changeFrameAndCopy(f1);

    Quaternion expected(Quaternion::fromMatrix(X2.E));
    EXPECT_TRUE(unit_test_utils::checkVector4dEpsilonClose(orientation_copy.getOrientation(),expected,unit_test_utils::TEST_PREC));
    EXPECT_STREQ(orientation.getReferenceFrame()->getName().c_str(),"f2");
    EXPECT_STREQ(orientation_copy.getReferenceFrame()->getName().c_str(),"f1");

    orientation = orientation_copy;
    orientation.changeFrameAndCopy(ReferenceFrame::getWorldFrame(),orientation_copy);

    expected = Quaternion(0.,0.,0.,1.);
    EXPECT_TRUE(unit_test_utils::checkVector4dEpsilonClose(orientation_copy.getOrientation(),expected,unit_test_utils::TEST_PREC));
    EXPECT_STREQ(orientation.getReferenceFrame()->getName().c_str(),"f1");
    EXPECT_STREQ(orientation_copy.getReferenceFrame()->getName().c_str(),"World");

    orientation = orientation_copy;
    orientation.changeFrameAndCopy(f3,orientation_copy);
    expected = Quaternion::fromMatrix(Xrotz(M_PI_2).E);
    EXPECT_TRUE(unit_test_utils::checkVector4dEpsilonClose(orientation_copy.getOrientation(),expected,unit_test_utils::TEST_PREC));
    EXPECT_STREQ(orientation.getReferenceFrame()->getName().c_str(),"World");
    EXPECT_STREQ(orientation_copy.getReferenceFrame()->getName().c_str(),"f3");
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    return RUN_ALL_TESTS();
}