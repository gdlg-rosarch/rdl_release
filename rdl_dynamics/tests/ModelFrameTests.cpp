//
// Created by jordan on 8/14/16.
//

#include "gtest/gtest.h"
#include "Fixtures.h"
#include "rdl_dynamics/FramePoint.hpp"
#include "rdl_dynamics/Kinematics.h"
#include "rdl_dynamics/Dynamics.h"
#include "UnitTestUtils.hpp"
#include <math.h>

using namespace RobotDynamics::Math;
using namespace RobotDynamics;

TEST_F(FixedBase3DoFPlanar, testPlanar3LinkPendulum)
{
    Q[0] = 0.; Q[1] = 0.; Q[2] = 0.;
    QDot[0] = 0.; QDot[1] = 0.; QDot[2] = 0.;
    QDDot[0] = 0.; QDDot[1] = 0.; QDDot[2] = 0.;

    updateKinematics(*model, Q, QDot, QDDot);

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[1]->getTransformFromParent().r,Vector3d(0.,0.,0.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[2]->getTransformFromParent().r,Vector3d(0.,0.,-1.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[3]->getTransformFromParent().r,Vector3d(0.,0.,-2.),unit_test_utils::TEST_PREC));

    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Matrix3dIdentity,model->bodyFrames[0]->getTransformToParent().E,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Matrix3dIdentity,model->bodyFrames[1]->getTransformToParent().E,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Matrix3dIdentity,model->bodyFrames[2]->getTransformToParent().E,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Matrix3dIdentity,model->bodyFrames[3]->getTransformToParent().E,unit_test_utils::TEST_PREC));

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[0]->getTransformToParent().r,Vector3d(0,0,0),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[1]->getTransformToParent().r,Vector3d(0,0,0),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[2]->getTransformToParent().r,Vector3d(0,0,1),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[3]->getTransformToParent().r,Vector3d(0,0,2),unit_test_utils::TEST_PREC));

    FramePointd p(model->bodyFrames[model->mBodyNameMap["body_c"]].get(),0.,0.,-3.);
    FramePointd p_com(model->bodyCenteredFrames[model->mBodyNameMap["body_c"]].get(),0.,0.,0.);

    Vector3d expected_c_com = 1/(2*fixed_body.mMass+body_c.mMass)*((fixed_body.mCenterOfMass+X_fixed_c.r)*fixed_body.mMass + (X_fixed_c.r+X_fixed_c2.r + fixed_body.mCenterOfMass)*fixed_body.mMass + body_c.mCenterOfMass*body_c.mMass);
    // Need to make sure that COM frames are properly updated when fixed joint bodies are involed in which bodies mass properties are combined.
    EXPECT_NEAR(model->I[3].h[0]/model->I[3].m,expected_c_com[0],unit_test_utils::TEST_PREC);
    EXPECT_NEAR(model->I[3].h[1]/model->I[3].m,expected_c_com[1],unit_test_utils::TEST_PREC);
    EXPECT_NEAR(model->I[3].h[2]/model->I[3].m,expected_c_com[2],unit_test_utils::TEST_PREC);

    p.changeFrame(ReferenceFrame::getWorldFrame().get());
    FramePointd p_com_point(model->bodyFrames[model->mBodyNameMap["body_c"]].get(),expected_c_com);
    p_com.changeFrame(ReferenceFrame::getWorldFrame().get());
    p_com_point.changeFrame(model->worldFrame);

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(p.vec(),Vector3d(0.,0.,-6.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(p_com.vec(),p_com_point.vec(),unit_test_utils::TEST_PREC));

    p.changeFrame(model->bodyFrames[model->mBodyNameMap["body_b"]].get());

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(p.vec(),Vector3d(0.,0.,-5.),unit_test_utils::TEST_PREC));

    p.changeFrame(model->bodyFrames[model->mBodyNameMap["body_c"]].get());

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(p.vec(),Vector3d(0.,0.,-3.),unit_test_utils::TEST_PREC));

    p_com.changeFrame(model->bodyCenteredFrames[model->mBodyNameMap["body_c"]].get());

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(p_com.vec(),Vector3d(0.,0.,0.),unit_test_utils::TEST_PREC));

    double j1Angle = M_PI_2;
    double j2Angle = -M_PI_2;
    double j3Angle = 0.;

    Q[0] = j1Angle; Q[1] = j2Angle; Q[2] = j3Angle;
    QDot[0] = 0.; QDot[1] = 0.; QDot[2] = 0.;
    QDDot[0] = 0.; QDDot[1] = 0.; QDDot[2] = 0.;

    updateKinematicsCustom(*model, &Q, nullptr, nullptr);

    p_com_point.setIncludingFrame(expected_c_com,model->bodyFrames[model->mBodyNameMap["body_c"]].get());

    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Matrix3dIdentity,model->bodyFrames[0]->getTransformToParent().E,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Xroty(-j1Angle).E,model->bodyFrames[1]->getTransformToParent().E,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Xroty(-j2Angle).E,model->bodyFrames[2]->getTransformToParent().E,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Xroty(-j3Angle).E,model->bodyFrames[3]->getTransformToParent().E,unit_test_utils::TEST_PREC));

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[0]->getTransformToParent().r,Vector3d(0,0,0),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[1]->getTransformToParent().r,Vector3d(0,0,0),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[2]->getTransformToParent().r,Vector3d(1,0,0),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[3]->getTransformToParent().r,Vector3d(0,0,2),unit_test_utils::TEST_PREC));

    p.changeFrame(ReferenceFrame::getWorldFrame().get());
    p_com.changeFrame(ReferenceFrame::getWorldFrame().get());
    p_com_point.changeFrame(model->worldFrame);

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(p.vec(),Vector3d(-1,0,-5),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(p_com.vec(),p_com_point.vec(),unit_test_utils::TEST_PREC));

    p.changeFrame(model->bodyFrames[model->mBodyNameMap["body_b"]].get());

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(p.vec(),Vector3d(0,0,-5),unit_test_utils::TEST_PREC));

    p.changeFrame(model->bodyFrames[model->mBodyNameMap["body_c"]].get());

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(p.vec(),Vector3d(0,0,-3),unit_test_utils::TEST_PREC));

    p.changeFrame(model->bodyFrames[model->mBodyNameMap["body_a"]].get());

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(p.vec(),Vector3d(5,0,-1),unit_test_utils::TEST_PREC));

    FramePointd p2(model->bodyFrames[model->mBodyNameMap["body_c"]].get(),1.,1.,1.);

    try
    {
        p.add(p2);
    }
    catch(RobotDynamics::ReferenceFrameException &e)
    {
        EXPECT_STREQ(e.what(),"Reference frames do not match!");
    }

    EXPECT_STREQ(model->getBodyFrame("body_a")->getName().c_str(),"body_a");
    EXPECT_STREQ(model->getBodyFrame("body_b")->getName().c_str(),"body_b");
    EXPECT_STREQ(model->getBodyFrame("body_c")->getName().c_str(),"body_c");
    EXPECT_FALSE(model->getBodyFrame("body_d"));
}

TEST_F(FixedBase3DoFPlanar, testPlanar3LinkPendulumFixedBodies)
{
    Q[0] = 0.; Q[1] = 0.; Q[2] = 0.;
    QDot[0] = 0.; QDot[1] = 0.; QDot[2] = 0.;
    QDDot[0] = 0.; QDDot[1] = 0.; QDDot[2] = 0.;

    updateKinematics(*model, Q, QDot, QDDot);

    FramePointd pFixedFrameOrigin(model->fixedBodyFrames[0].get(),0.,0.,0.);
    FramePointd pFixedFrameOrigin_2(model->fixedBodyFrames[1].get(),0.,0.,0.);

    pFixedFrameOrigin.changeFrame(ReferenceFrame::getWorldFrame().get());
    pFixedFrameOrigin_2.changeFrame(ReferenceFrame::getWorldFrame().get());

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(Vector3d(0.02,0.01,-3.25),pFixedFrameOrigin.vec(),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(Vector3d(-0.02,-0.01,-2.75),pFixedFrameOrigin_2.vec(),unit_test_utils::TEST_PREC));

    double j1Angle = M_PI_2;
    double j2Angle = -M_PI_2;
    double j3Angle = 0.;

    Q[0] = j1Angle; Q[1] = j2Angle; Q[2] = j3Angle;
    QDot[0] = 0.; QDot[1] = 0.; QDot[2] = 0.;
    QDDot[0] = 0.; QDDot[1] = 0.; QDDot[2] = 0.;

    updateKinematics(*model, Q, QDot, QDDot);

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(Vector3d(0.02,0.01,-3.25),pFixedFrameOrigin.vec(),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(Vector3d(-0.02,-0.01,-2.75),pFixedFrameOrigin_2.vec(),unit_test_utils::TEST_PREC));

    pFixedFrameOrigin.setIncludingFrame(0.,0.,0.,model->fixedBodyFrames[0].get());
    pFixedFrameOrigin_2.setIncludingFrame(0.,0.,0.,model->fixedBodyFrames[1].get());

    pFixedFrameOrigin.changeFrame(ReferenceFrame::getWorldFrame());
    pFixedFrameOrigin_2.changeFrame(ReferenceFrame::getWorldFrame());
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(Vector3d((0.02-1),0.01,-2.25),pFixedFrameOrigin.vec(),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(Vector3d((-0.02-1),-0.01,-1.75),pFixedFrameOrigin_2.vec(),unit_test_utils::TEST_PREC));
}

// Kinda a weird test, but I had a bug where inverseDynamics didn't update the fixed body frames, so it gets a test.
TEST_F(FixedBase3DoFPlanar, testPlanar3LinkPendulumFixedBodiesUsingInverseDynamics)
{
    Q[0] = 0.; Q[1] = 0.; Q[2] = 0.;
    QDot[0] = 0.; QDot[1] = 0.; QDot[2] = 0.;
    QDDot[0] = 0.; QDDot[1] = 0.; QDDot[2] = 0.;

    inverseDynamics(*model, Q, QDot, QDDot,Tau);

    FramePointd pFixedFrameOrigin(model->fixedBodyFrames[0].get(),0.,0.,0.);
    FramePointd pFixedFrameOrigin_2(model->fixedBodyFrames[1].get(),0.,0.,0.);

    pFixedFrameOrigin.changeFrame(ReferenceFrame::getWorldFrame().get());
    pFixedFrameOrigin_2.changeFrame(ReferenceFrame::getWorldFrame().get());

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(Vector3d(0.02,0.01,-3.25),pFixedFrameOrigin.vec(),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(Vector3d(-0.02,-0.01,-2.75),pFixedFrameOrigin_2.vec(),unit_test_utils::TEST_PREC));

    double j1Angle = M_PI_2;
    double j2Angle = -M_PI_2;
    double j3Angle = 0.;

    Q[0] = j1Angle; Q[1] = j2Angle; Q[2] = j3Angle;
    QDot[0] = 0.; QDot[1] = 0.; QDot[2] = 0.;
    QDDot[0] = 0.; QDDot[1] = 0.; QDDot[2] = 0.;

    inverseDynamics(*model, Q, QDot, QDDot,Tau);

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(Vector3d(0.02,0.01,-3.25),pFixedFrameOrigin.vec(),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(Vector3d(-0.02,-0.01,-2.75),pFixedFrameOrigin_2.vec(),unit_test_utils::TEST_PREC));

    pFixedFrameOrigin.setIncludingFrame(0.,0.,0.,model->fixedBodyFrames[0].get());
    pFixedFrameOrigin_2.setIncludingFrame(0.,0.,0.,model->fixedBodyFrames[1].get());

    pFixedFrameOrigin.changeFrame(ReferenceFrame::getWorldFrame());
    pFixedFrameOrigin_2.changeFrame(ReferenceFrame::getWorldFrame());
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(Vector3d((0.02-1),0.01,-2.25),pFixedFrameOrigin.vec(),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(Vector3d((-0.02-1),-0.01,-1.75),pFixedFrameOrigin_2.vec(),unit_test_utils::TEST_PREC));
}

TEST_F(FixedBaseTwoChain6DoF3D,testPointChangeFrameWith2ChainModel)
{
    Q[0] = 0.; Q[1] = 0.; Q[2] = 0.; Q[3] = 0.; Q[4] = 0.; Q[5] = 0.;
    QDot[0] = 0.; QDot[1] = 0.; QDot[2] = 0.; QDot[3] = 0.; QDot[4] = 0.; QDot[5] = 0.;
    QDDot[0] = 0.; QDDot[1] = 0.; QDDot[2] = 0.; QDDot[3] = 0.; QDDot[4] = 0.; QDDot[5] = 0.;

    updateKinematics(*model, Q, QDot, QDDot);

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[1]->getTransformFromParent().r,Vector3d(0.,0.,0.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[2]->getTransformFromParent().r,Vector3d(0.,0.,-1.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[3]->getTransformFromParent().r,Vector3d(0.,-2.,0.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[4]->getTransformFromParent().r,Vector3d(0.,0.,0.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[5]->getTransformFromParent().r,Vector3d(2.,0.,0.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[6]->getTransformFromParent().r,Vector3d(0.,0.,1.),unit_test_utils::TEST_PREC));

    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Matrix3dIdentity,model->bodyFrames[0]->getTransformToParent().E,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Xrotx(0.).E,model->bodyFrames[1]->getTransformToParent().E,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Xrotx(-M_PI_2).E,model->bodyFrames[2]->getTransformToParent().E,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Xrotx(-M_PI_2).E,model->bodyFrames[3]->getTransformToParent().E,unit_test_utils::TEST_PREC));

    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Xrotx(0.).E,model->bodyFrames[4]->getTransformToParent().E,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Xrotz(-M_PI_2).E,model->bodyFrames[5]->getTransformToParent().E,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Xrotx(-M_PI_2).E,model->bodyFrames[6]->getTransformToParent().E,unit_test_utils::TEST_PREC));

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[0]->getTransformToParent().r,Vector3d(0.,0.,0.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[1]->getTransformToParent().r,Vector3d(0.,0.,0.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[2]->getTransformToParent().r,Vector3d(0.,1.,0.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[3]->getTransformToParent().r,Vector3d(0.,0.,-2.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[4]->getTransformToParent().r,Vector3d(0.,0.,0.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[5]->getTransformToParent().r,Vector3d(0.,2.,0.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[6]->getTransformToParent().r,Vector3d(0.,-1.,0.),unit_test_utils::TEST_PREC));

    Math::FramePointd p(model->bodyFrames[6].get(),3.,0.,0.);
    Math::FramePointd p_com(model->bodyCenteredFrames[6].get(),0.,0.,0.);

    p.changeFrame(model->bodyFrames[6].get());
    p_com.changeFrame(model->bodyCenteredFrames[6].get());
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(p.vec(),Vector3d(3.,0.,0.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(p_com.vec(),Vector3d(0.,0.,0.),unit_test_utils::TEST_PREC));

    p.changeFrame(ReferenceFrame::getWorldFrame().get());
    p_com.changeFrame(ReferenceFrame::getWorldFrame().get());
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(p.vec(),Vector3d(2.,3.,1.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(p_com.vec(),Vector3d(2.5,0.,1.),unit_test_utils::TEST_PREC));

    p.changeFrame(model->bodyFrames[3].get());
    p_com.changeFrame(model->bodyFrames[3].get());
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(p.vec(),Vector3d(2.,-3.,-4.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(p_com.vec(),Vector3d(2.5,0.,-4.),unit_test_utils::TEST_PREC));

    Q[0] = M_PI_2;

    updateKinematics(*model, Q, QDot, QDDot);

    Math::FramePointd p2(model->bodyFrames[6].get(),3.,0.,0.);

    p2.changeFrame(model->bodyFrames[3].get());

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(p2.vec(),Vector3d(2.,-1.,0.),unit_test_utils::TEST_PREC));
}

TEST_F(FloatingBaseWith2SingleDofJoints, testPointChangeFrames)
{
    Q[0] = 0.; Q[1] = 0.; Q[2] = 0.; Q[3] = 0.; Q[4] = 0.; Q[5] = 0.; Q[6] = 0.; Q[7] = 0.;Q[8] = 1.;
    QDot[0] = 0.; QDot[1] = 0.; QDot[2] = 0.; QDot[3] = 0.; QDot[4] = 0.; QDot[5] = 0.;
    QDDot[0] = 0.; QDDot[1] = 0.; QDDot[2] = 0.; QDDot[3] = 0.; QDDot[4] = 0.; QDDot[5] = 0.;

    updateKinematics(*model, Q, QDot, QDDot);

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[0]->getTransformFromParent().r,Vector3d(0.,0.,0.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[1]->getTransformFromParent().r,Vector3d(0.,0.,0.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[2]->getTransformFromParent().r,Vector3d(0.,0.,0.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[3]->getTransformFromParent().r,Vector3d(-1.5,0.,0.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[4]->getTransformFromParent().r,Vector3d(1.5,0.,0.),unit_test_utils::TEST_PREC));

    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Matrix3dIdentity,model->bodyFrames[0]->getTransformFromParent().E,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Matrix3dIdentity,model->bodyFrames[1]->getTransformFromParent().E,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Xtrans(Vector3d(0.,0.,0.)).E,model->bodyFrames[2]->getTransformFromParent().E,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Matrix3dIdentity,model->bodyFrames[3]->getTransformFromParent().E,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Xrotz(M_PI_2).E,model->bodyFrames[4]->getTransformFromParent().E,unit_test_utils::TEST_PREC));

    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Matrix3dIdentity,model->bodyFrames[0]->getTransformToParent().E,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Matrix3dIdentity,model->bodyFrames[1]->getTransformToParent().E,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Matrix3dIdentity,model->bodyFrames[2]->getTransformToParent().E,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Matrix3dIdentity,model->bodyFrames[3]->getTransformToParent().E,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(RobotDynamics::Math::Xrotz(-M_PI_2).E,model->bodyFrames[4]->getTransformToParent().E,unit_test_utils::TEST_PREC));

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[0]->getTransformToParent().r,Vector3d(0.,0.,0.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[1]->getTransformToParent().r,Vector3d(0.,0.,0.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[2]->getTransformToParent().r,Vector3d(0.,0.,0.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[3]->getTransformToParent().r,Vector3d(1.5,0.,0.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(model->bodyFrames[4]->getTransformToParent().r,Vector3d(0.,1.5,0.),unit_test_utils::TEST_PREC));

    FramePointd p1_end_effector(model->bodyFrames[model->mBodyNameMap["body_1"]].get(),0.,0.,-1.);
    FramePointd p2_end_effector(model->bodyFrames[model->mBodyNameMap["body_2"]].get(),0.,0.,1.);
    FramePointd root_origin(model->bodyFrames[model->mBodyNameMap["floating_body"]].get(),0.,0.,0.);
    FramePointd world_origin(ReferenceFrame::getWorldFrame().get(),0.,0.,0.);

    p1_end_effector.changeFrame(ReferenceFrame::getWorldFrame().get());
    p2_end_effector.changeFrame(model->bodyFrames[model->mBodyNameMap["body_1"]].get());
    root_origin.changeFrame(model->bodyFrames[model->mBodyNameMap["body_2"]].get());

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(p1_end_effector.vec(),Vector3d(-1.5,0.,-1.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(p2_end_effector.vec(),Vector3d(3.,0.,1.),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(root_origin.vec(),Vector3d(0.,1.5,0.),unit_test_utils::TEST_PREC));

    Q[0] = 1.; Q[1] = 2.; Q[2] = 3.;

    updateKinematics(*model, Q, QDot, QDDot);

    world_origin.changeFrame(model->bodyFrames[model->mBodyNameMap["floating_body"]].get());

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(world_origin.vec(),Vector3d(-1.,-2.,-3.),unit_test_utils::TEST_PREC));

    double a = M_PI_2;
    Q[0] = 1.; Q[1] = 2.; Q[2] = 3.; Q[3] = 0.; Q[4] = 0.; Q[5] = sin(a/2); Q[6] = 0.; Q[7] = 0.;Q[8] = cos(a/2);

    world_origin.setIncludingFrame(0.,0.,0.,ReferenceFrame::getWorldFrame().get());
    updateKinematics(*model, Q, QDot, QDDot);

    root_origin.changeFrame(ReferenceFrame::getWorldFrame().get());

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(root_origin.vec(),Vector3d(1.,2.,3.),unit_test_utils::TEST_PREC));

    world_origin.changeFrame(model->bodyFrames[model->mBodyNameMap["floating_body"]].get());

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(world_origin.vec(),Vector3d(-2.,1.,-3.),unit_test_utils::TEST_PREC));
}

TEST_F(FixedAndMovableJoint, testComFramesWithFixedJoint)
{
    EXPECT_TRUE(model_fixed->bodyCenteredFrames[1]->getTransformToParent().r[0]==-1.);
    EXPECT_TRUE(model_fixed->bodyCenteredFrames[1]->getTransformToParent().r[1]==-0.5);
    EXPECT_TRUE(model_fixed->bodyCenteredFrames[1]->getTransformToParent().r[2]==0.);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}