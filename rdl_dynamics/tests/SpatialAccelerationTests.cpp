//
// Created by jordan on 10/10/16.
//
#include "UnitTestUtils.hpp"
#include "Fixtures.h"

#include "rdl_dynamics/Kinematics.h"

using namespace RobotDynamics::Math;

class SpatialAccelerationTests : public testing::Test
{
public:
    SpatialAccelerationTests()
    {

    }

    void SetUp()
    {

    }

    void TearDown()
    {

    }
};

TEST_F(FixedBase3DoFPlanar,testAdd)
{
    SpatialAcceleration a1(model->bodyFrames[1].get(),ReferenceFrame::getWorldFrame().get(),ReferenceFrame::getWorldFrame().get(),0.1,0.2,0.3,-0.1,-0.2,-0.3);
    SpatialAcceleration a2(model->bodyFrames[2].get(),model->bodyFrames[1].get(),ReferenceFrame::getWorldFrame().get(),0.1,0.2,0.3,-0.1,-0.2,-0.3);

    SpatialVector av_1 = a1;
    SpatialVector av_2 = a2;

    SpatialAcceleration a3 = a1+a2;

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(av_1+av_2,a3,unit_test_utils::TEST_PREC));

    SpatialAcceleration a4 = a3-a1;

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(a4,a2,unit_test_utils::TEST_PREC));

    a4 = a3-a2;

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(a4,a1,unit_test_utils::TEST_PREC));

    try
    {
        SpatialAcceleration a5(model->bodyFrames[2].get(),model->bodyFrames[1].get(),ReferenceFrame::getWorldFrame().get(),0.1,0.2,0.3,-0.1,-0.2,-0.3);
        a5-=a1;
    }
    catch(ReferenceFrameException &e)
    {
        EXPECT_STREQ("Reference frame mismatch during subtraction of SpatialAccelerations!",e.what());
    }
}

TEST_F(FixedBase3DoFPlanar,testChangeFrame)
{

    Q[0] = -0.25; Q[1] = 0.12; Q[2] = 1.1;
    QDot[0] = 0.; QDot[1] = 0.; QDot[2] = 0.;
    QDDot[0] = 0.; QDDot[1] = 0.; QDDot[2] = 0.;

    updateKinematics(*model, Q, QDot, QDDot);

    ReferenceFrame* newFrame = model->bodyFrames[2].get();
    ReferenceFrame* currentExpressedInFrame = ReferenceFrame::getWorldFrame().get();

    SpatialAcceleration a2(model->bodyFrames[2].get(),model->bodyFrames[1].get(),currentExpressedInFrame,0.1,0.2,0.3,-0.1,-0.2,-0.3);
    SpatialVector a2_v = a2;

    SpatialTransform X = currentExpressedInFrame->getTransformToDesiredFrame(newFrame);

    SpatialVector v_expected = X.apply(a2_v);
    a2.changeFrame(newFrame);

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(a2,v_expected,unit_test_utils::TEST_PREC));
}

TEST_F(FixedBase3DoFPlanar,testchangeFrameWithRelativeMotion)
{
    Q[0] = -0.25; Q[1] = 0.12; Q[2] = 1.1;
    QDot[0] = 0.1; QDot[1] = 0.2; QDot[2] = 0.3;
    QDDot[0] = 0.5; QDDot[1] = 0.4; QDDot[2] = 0.3;

    updateKinematics(*model, Q, QDot, QDDot);

    ReferenceFrame* newFrame = model->bodyFrames[2].get();
    ReferenceFrame* currentExpressedInFrame = ReferenceFrame::getWorldFrame().get();
    ReferenceFrame* bodyFrame = model->bodyFrames[2].get();
    ReferenceFrame* baseFrame = model->bodyFrames[1].get();

    SpatialAcceleration a2(bodyFrame,baseFrame,currentExpressedInFrame,0.1,0.2,0.3,-0.1,-0.2,-0.3);
    SpatialVector a2_v = a2;

    SpatialMotion m2(bodyFrame,baseFrame,currentExpressedInFrame,0.5,-0.1,-0,.11,0.23,0.9);
    SpatialMotion m1(currentExpressedInFrame,newFrame,currentExpressedInFrame,0.1,0.2,0.3,0.4,0.5,0.6);

    SpatialMotion expected = m1%m2;

    expected.set(expected.toMotionVector() + a2_v);

    SpatialTransform X = currentExpressedInFrame->getTransformToDesiredFrame(newFrame);
    expected.transform(X);

    a2.changeFrameWithRelativeMotion(newFrame,m1,m2);
    EXPECT_TRUE(a2.getReferenceFrame()==newFrame);

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(a2,expected,unit_test_utils::TEST_PREC));
}

TEST_F(FixedBase3DoFPlanar,testchangeFrameWithRelativeMotionSameFrame)
{

    Q[0] = -0.25; Q[1] = 0.12; Q[2] = 1.1;
    QDot[0] = 0.1; QDot[1] = -0.2; QDot[2] = -0.3;
    QDDot[0] = 1.; QDDot[1] = 2.; QDDot[2] = 3.;

    updateKinematics(*model, Q, QDot, QDDot);

    ReferenceFrame* newFrame = model->bodyFrames[2].get();
    ReferenceFrame* currentExpressedInFrame = ReferenceFrame::getWorldFrame().get();

    SpatialAcceleration a2(model->bodyFrames[2].get(),model->bodyFrames[1].get(),currentExpressedInFrame,0.1,0.2,0.3,-0.1,-0.2,-0.3);

    SpatialMotion m2(currentExpressedInFrame,newFrame,currentExpressedInFrame,0.5,-0.1,-0,.11,0.23,0.9);
    SpatialMotion m1(model->bodyFrames[2].get(),model->bodyFrames[1].get(),currentExpressedInFrame,0.1,0.2,0.3,0.4,0.5,0.6);

    SpatialAcceleration expected = a2;
    a2.changeFrameWithRelativeMotion(a2.getReferenceFrame(),m1,m2);

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(a2,expected,unit_test_utils::TEST_PREC));
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}