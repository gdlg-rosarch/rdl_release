//
// Created by jordan on 9/27/16.
//

#include "rdl_dynamics/rdl_math.hpp"
#include "UnitTestUtils.hpp"
#include "Fixtures.h"
#include <gtest/gtest.h>

using namespace RobotDynamics::Math;

class SpatialForceTests : public ::testing::Test
{
public:
    SpatialForceTests(){};

    void SetUp()
    {

    }

    void TearDown()
    {
        
    }
};

TEST_F(FixedBase3DoFPlanar,changeFrameAndCopy)
{
    SpatialForce v1(model->getBodyFrame("body_b"), 1.,2.,3.,4.,5.,6.);
    SpatialForce v2;

    randomizeStates();
    updateKinematicsCustom(*model,&Q,nullptr,nullptr);

    v2=v1.changeFrameAndCopy(ReferenceFrame::getWorldFrame());

    EXPECT_FALSE(v2.mx()==v1.mx());
    EXPECT_FALSE(v2.my()==v1.my());
    EXPECT_FALSE(v2.mz()==v1.mz());
    EXPECT_FALSE(v2.fx()==v1.fx());
//    EXPECT_FALSE(v2.fy()==v1.fy());  // its a rot about y axis, so this one isn't supposed to change
    EXPECT_FALSE(v2.fz()==v1.fz());

    v1.changeFrame(ReferenceFrame::getWorldFrame());

    EXPECT_TRUE(v2.mx()==v1.mx());
    EXPECT_TRUE(v2.my()==v1.my());
    EXPECT_TRUE(v2.mz()==v1.mz());
    EXPECT_TRUE(v2.fx()==v1.fx());
    EXPECT_TRUE(v2.fy()==v1.fy());
    EXPECT_TRUE(v2.fz()==v1.fz());

    SpatialForce v3(model->getBodyFrame("body_b"), 1.,2.,3.,4.,5.,6.);
    SpatialForce v4;

    v3.changeFrameAndCopy(ReferenceFrame::getWorldFrame(),v4);

    EXPECT_FALSE(v3.mx()==v4.mx());
    EXPECT_FALSE(v3.my()==v4.my());
    EXPECT_FALSE(v3.mz()==v4.mz());
    EXPECT_FALSE(v3.fx()==v4.fx());
//    EXPECT_FALSE(v3.fy()==v4.fy());  // its a rot about y axis, so this one isn't supposed to change
    EXPECT_FALSE(v3.fz()==v4.fz());

    v3.changeFrame(ReferenceFrame::getWorldFrame());

    EXPECT_TRUE(v3.mx()==v4.mx());
    EXPECT_TRUE(v3.my()==v4.my());
    EXPECT_TRUE(v3.mz()==v4.mz());
    EXPECT_TRUE(v3.fx()==v4.fx());
    EXPECT_TRUE(v3.fy()==v4.fy());
    EXPECT_TRUE(v3.fz()==v4.fz());
}

TEST_F(FixedBase3DoFPlanar, test_framed_accessors)
{
    SpatialForce v1(model->getBodyFrame("body_a"), 1.,2.,3.,4.,5.,6.);

    FrameVector fvl = v1.getFramedLinearPart();
    FrameVector fva = v1.getFramedAngularPart();

    EXPECT_STREQ(fvl.getReferenceFrame()->getName().c_str(),v1.getReferenceFrame()->getName().c_str());

    EXPECT_EQ(fva.x(),v1.mx());
    EXPECT_EQ(fva.y(),v1.my());
    EXPECT_EQ(fva.z(),v1.mz());

    EXPECT_EQ(fvl.x(),v1.fx());
    EXPECT_EQ(fvl.y(),v1.fy());
    EXPECT_EQ(fvl.z(),v1.fz());

    v1.setIncludingFrame(model->getBodyFrame("body_b"),ForceVector(0.1,-0.1,0.2,-0.2,0.3,-0.3));

    EXPECT_STREQ(v1.getReferenceFrame()->getName().c_str(),"body_b");
    EXPECT_EQ(v1.mx(),0.1);
    EXPECT_EQ(v1.my(),-0.1);
    EXPECT_EQ(v1.mz(),0.2);

    EXPECT_EQ(v1.fx(),-0.2);
    EXPECT_EQ(v1.fy(),0.3);
    EXPECT_EQ(v1.fz(),-0.3);

    v1.setIncludingFrame(model->getBodyFrame("body_c"),3.,4.,5.,6.,7.,8.);

    EXPECT_STREQ(v1.getReferenceFrame()->getName().c_str(),"body_c");
    EXPECT_EQ(v1.mx(),3.);
    EXPECT_EQ(v1.my(),4.);
    EXPECT_EQ(v1.mz(),5.);

    EXPECT_EQ(v1.fx(),6.);
    EXPECT_EQ(v1.fy(),7.);
    EXPECT_EQ(v1.fz(),8.);
}

TEST_F(FixedBase3DoFPlanar, testAdd)
{
    SpatialForce v1(model->getBodyFrame("body_a"), 0., 1., 0., 0., 0., 0.);
    SpatialForce v2(model->getBodyFrame("body_b"), 0., 1., 0., 0., 0., 0.);

    try
    {
        //Frame mismatches, should throw
        SpatialForce v3 = v1 + v2;
    }
    catch(ReferenceFrameException &e)
    {
        EXPECT_STREQ("Reference frames do not match!",e.what());
    }

    SpatialForce v4(model->getBodyFrame("body_b"), 0., 1., 0., 0., 0., 0.);

    SpatialForce v5 = v2 + v4;

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(v5, SpatialVector(0., 2., 0., 0., 0., 0.), unit_test_utils::TEST_PREC));

    SpatialForce v6 = v5;

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(v5, SpatialVector(0., 2., 0., 0., 0., 0.), unit_test_utils::TEST_PREC));

    SpatialForce v7(model->getBodyFrame("body_b"), -0.2, 0.3, -0.4, 0.5, -0.6, 0.7);

    EXPECT_EQ(v7.mx(), -0.2);
    EXPECT_EQ(v7.my(), 0.3);
    EXPECT_EQ(v7.mz(), -0.4);
    EXPECT_EQ(v7.fx(), 0.5);
    EXPECT_EQ(v7.fy(), -0.6);
    EXPECT_EQ(v7.fz(), 0.7);

    SpatialForce v8 = v6-v7;

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(v8, SpatialVector(0.2, 1.7, 0.4, -0.5, 0.6, -0.7), unit_test_utils::TEST_PREC));
}

TEST_F(SpatialForceTests, testChangeFrame)
{
    Vector3d rot(1.1, 1.2, 1.3);
    Vector3d trans(1.1, 1.2, 1.3);

    SpatialTransform X_st;
    X_st.r = trans;

    SpatialMatrix X_66_matrix(SpatialMatrix::Zero(6, 6));
    X_66_matrix = Xrotz_mat(rot[2]) * Xroty_mat(rot[1]) * Xrotx_mat(rot[0]) * Xtrans_mat(trans);
    X_st.E = X_66_matrix.block<3, 3>(0, 0);

    std::shared_ptr<ReferenceFrame> frameA(new ReferenceFrame("frameA", ReferenceFrame::getWorldFrame().get(), X_st, 0, 0));

    ForceVector v(1.1, 2.1, 3.1, 4.1, 5.1, 6.1);

    SpatialForce m(frameA.get(), v);

    v.transform(X_st.inverse());

    m.changeFrame(ReferenceFrame::getWorldFrame().get());

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(m, v, unit_test_utils::TEST_PREC));
}

TEST_F(SpatialForceTests, testFrameCorrectness)
{
    SpatialForce F(ReferenceFrame::getWorldFrame().get(),0.,0.,0.,1.,2.,3.);
    SpatialTransform X_z = Xrotz(M_PI_2);

    std::shared_ptr<ReferenceFrame> R_2(new ReferenceFrame("Frame2",ReferenceFrame::getWorldFrame().get(),X_z,true,1));

    F.changeFrame(R_2.get());

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(F,SpatialVector(0.,0.,0.,2.,-1.,3.),unit_test_utils::TEST_PREC));

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(F.transform_copy(R_2->getTransformToRoot()),SpatialVector(0.,0.,0.,1.,2.,3.),unit_test_utils::TEST_PREC));
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}