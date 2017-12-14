#include <gtest/gtest.h>

#include "UnitTestUtils.hpp"
#include "Fixtures.h"
#include "Human36Fixture.h"
#include "rdl_dynamics/rdl_utils.h"

using namespace std;
using namespace RobotDynamics;
using namespace RobotDynamics::Math;

struct RdlUtilsTests : public testing::Test
{

};

TEST_F(RdlUtilsTests, testGetDofName)
{
    SpatialVector rx(1.,0.,0.,0.,0.,0.);
    SpatialVector ry(0.,1.,0.,0.,0.,0.);
    SpatialVector rz(0.,0.,1.,0.,0.,0.);
    SpatialVector tx(0.,0.,0.,1.,0.,0.);
    SpatialVector ty(0.,0.,0.,0.,1.,0.);
    SpatialVector tz(0.,0.,0.,0.,0.,1.);

    EXPECT_STREQ("RX",Utils::getDofName(rx).c_str());
    EXPECT_STREQ("RY",Utils::getDofName(ry).c_str());
    EXPECT_STREQ("RZ",Utils::getDofName(rz).c_str());
    EXPECT_STREQ("TX",Utils::getDofName(tx).c_str());
    EXPECT_STREQ("TY",Utils::getDofName(ty).c_str());
    EXPECT_STREQ("TZ",Utils::getDofName(tz).c_str());

    SpatialVector c(1.,0.,0.,0.,0.,1.);

    ostringstream dof_stream(ostringstream::out);
    dof_stream << "custom (" << c.transpose() << ")";
    EXPECT_STREQ(dof_stream.str().c_str(),Utils::getDofName(c).c_str());
}

TEST_F(RdlUtilsTests, testGetBodyName)
{
    Model model;
    Body b1;
    b1.mIsVirtual=true;
    Joint j(JointType6DoF);
    j.mJointAxes[0]=SpatialVector(1.,0.,0.,0.,0.,0.);
    j.mJointAxes[1]=SpatialVector(0.,1.,0.,0.,0.,0.);
    j.mJointAxes[2]=SpatialVector(0.,0.,1.,0.,0.,0.);
    j.mJointAxes[3]=SpatialVector(0.,0.,0.,1.,0.,0.);
    j.mJointAxes[4]=SpatialVector(0.,0.,0.,0.,1.,0.);
    j.mJointAxes[5]=SpatialVector(0.,0.,0.,0.,0.,1.);
    unsigned int b_id=model.addBody (0,SpatialTransform(),j,b1,"b1");
    Body b2;
    unsigned int b_id2=model.addBody(3,SpatialTransform(),Joint(JointTypeRevoluteY),b2,"b2");

    EXPECT_STREQ("b2",Utils::getBodyName(model,b_id2).c_str());
    EXPECT_STREQ("",Utils::getBodyName(model,b_id).c_str());
}

TEST_F(FloatingBase12DoF, TestKineticEnergy)
{
    VectorNd q = VectorNd::Zero(model->q_size);
    VectorNd qdot = VectorNd::Zero(model->q_size);

    for(unsigned int i = 0; i < q.size(); i++)
    {
        q[i] = 0.1 * i;
        qdot[i] = 0.3 * i;
    }

    MatrixNd H = MatrixNd::Zero(model->q_size, model->q_size);
    compositeRigidBodyAlgorithm(*model, q, H, true);

    double kinetic_energy_ref = 0.5 * qdot.transpose() * H * qdot;
    double ke = Utils::calcKineticEnergy(*model, q, qdot);

    EXPECT_EQ(ke, kinetic_energy_ref);
}

TEST_F(RdlUtilsTests, TestPotentialEnergy)
{
    Model model;
    Matrix3d inertia = Matrix3d::Zero(3, 3);
    Body body(0.5, Vector3d(0., 0., 0.), inertia);
    Joint joint(SpatialVector(0., 0., 0., 1., 0., 0.), SpatialVector(0., 0., 0., 0., 1., 0.), SpatialVector(0., 0., 0., 0., 0., 1.));

    model.appendBody(Xtrans(Vector3d::Zero()), joint, body);

    VectorNd q = VectorNd::Zero(model.q_size);
    double potential_energy_zero = Utils::calcPotentialEnergy(model, q);
    EXPECT_EQ (0., potential_energy_zero);

    q[1] = 1.;
    double potential_energy_lifted = Utils::calcPotentialEnergy(model, q);
    EXPECT_EQ (4.905, potential_energy_lifted);
}

TEST_F(RdlUtilsTests, TestCOMSimple)
{
    Model model;
    Matrix3d inertia = Matrix3d::Zero(3, 3);
    Body body(123., Vector3d(0., 0., 0.), inertia);
    Joint joint(SpatialVector(0., 0., 0., 1., 0., 0.), SpatialVector(0., 0., 0., 0., 1., 0.), SpatialVector(0., 0., 0., 0., 0., 1.));

    model.appendBody(Xtrans(Vector3d::Zero()), joint, body);

    VectorNd q = VectorNd::Zero(model.q_size);
    VectorNd qdot = VectorNd::Zero(model.qdot_size);

    double mass;
    Vector3d com;
    Vector3d com_velocity;
    FramePointd p_com,pcom_2;
    FrameVector v_com;
    Utils::calcCenterOfMass(model, q, qdot, mass, com, &com_velocity);
    Utils::calcCenterOfMass(model, q, qdot, mass, p_com, &v_com);
    Utils::calcCenterOfMass(model, q, pcom_2);

    EXPECT_EQ (123., mass);
    EXPECT_TRUE(unit_test_utils::checkVector3dEq(Vector3d(0., 0., 0.), com));
    EXPECT_TRUE(unit_test_utils::checkVector3dEq(pcom_2.vec(),com));
    EXPECT_TRUE(unit_test_utils::checkVector3dEq(p_com.vec(),Vector3d(0.,0.,0.)));
    EXPECT_TRUE(unit_test_utils::checkVector3dEq(Vector3d(0., 0., 0.), com_velocity));

    q[1] = 1.;
    Utils::calcCenterOfMass(model, q, qdot, mass, com, &com_velocity);
    Utils::calcCenterOfMass(model,q,qdot,p_com);
    Utils::calcCenterOfMass(model, q, pcom_2);
    EXPECT_TRUE(unit_test_utils::checkVector3dEq(Vector3d(0., 1., 0.), p_com.vec()));
    EXPECT_TRUE(unit_test_utils::checkVector3dEq(pcom_2.vec(),com));
    EXPECT_TRUE(unit_test_utils::checkVector3dEq(Vector3d(0., 1., 0.), com));
    EXPECT_TRUE(unit_test_utils::checkVector3dEq(Vector3d(0., 0., 0.), com_velocity));

    qdot[1] = 1.;
    Utils::calcCenterOfMass(model, q, qdot, mass, com, &com_velocity);
    Utils::calcCenterOfMass(model,q,qdot,p_com);
    Utils::calcCenterOfMass(model, q, pcom_2);
    EXPECT_TRUE(unit_test_utils::checkVector3dEq(Vector3d(0., 1., 0.), com));
    EXPECT_TRUE(unit_test_utils::checkVector3dEq(pcom_2.vec(),com));
    EXPECT_TRUE(unit_test_utils::checkVector3dEq(Vector3d(0., 1., 0.), p_com.vec()));
    EXPECT_TRUE(unit_test_utils::checkVector3dEq(Vector3d(0., 1., 0.), com_velocity));
}

TEST_F(Human36, TestCOM)
{
    randomizeStates();

    double mass;
    Vector3d com;
    Vector3d com_velocity;
    Utils::calcCenterOfMass(*model, q, qdot, mass, com, &com_velocity);

    FramePointd p_com,pcom_2;
    FrameVector v_com;
    Utils::calcCenterOfMass(*model,q,qdot,p_com,&v_com);
    Utils::calcCenterOfMass(*model,q,pcom_2);

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(com,p_com.vec(),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(pcom_2.vec(),p_com.vec(),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(com_velocity,v_com,unit_test_utils::TEST_PREC));

    randomizeStates();

    Utils::calcCenterOfMass(*model_3dof, q, qdot, mass, com, &com_velocity);
    Utils::calcCenterOfMass(*model_3dof,q,qdot,mass,p_com,&v_com,nullptr);
    Utils::calcCenterOfMass(*model_3dof,q,pcom_2);
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(com,p_com.vec(),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(pcom_2.vec(),p_com.vec(),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(com_velocity,v_com,unit_test_utils::TEST_PREC));

    randomizeStates();

    Utils::calcCenterOfMass(*model_emulated, q, qdot, mass, com, &com_velocity);
    Utils::calcCenterOfMass(*model_emulated,q,qdot,p_com,&v_com);
    Utils::calcCenterOfMass(*model_emulated,q,pcom_2);
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(com,p_com.vec(),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(pcom_2.vec(),p_com.vec(),unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(com_velocity,v_com,unit_test_utils::TEST_PREC));
}

TEST_F(Human36, TestCOMCalcSubtreeMass)
{
    double subtreeMass,wholeRobotMass;
    Vector3d com,com_velocity;

    Utils::calcCenterOfMass(*model_emulated,q,qdot,wholeRobotMass,com,nullptr);
    subtreeMass = Utils::calcSubtreeMass(*model_emulated,0);

    EXPECT_EQ(subtreeMass,wholeRobotMass);

    subtreeMass = Utils::calcSubtreeMass(*model_emulated,body_id_emulated[BodyPelvis]+1);
    EXPECT_EQ(subtreeMass,10.3368+3.1609+1.001);
}

TEST_F(RdlUtilsTests, TestCOMJacobian)
{
    Model model;

    unsigned int id1 = model.addBody(0,Xtrans(Vector3d(0,1,0)),Joint(SpatialVector(1,0,0,0,0,0)),Body(2,Vector3d(0,1,0),Vector3d(0,0,0)));
    model.addBody(id1,Xtrans(Vector3d(0,1,0)),Joint(SpatialVector(1,0,0,0,0,0)),Body(3,Vector3d(0,1,0),Vector3d(0,0,0)));

    double wholeRobotMass;
    Vector3d com,com_velocity;

    MatrixNd J_com(3,model.qdot_size);
    J_com.setZero();

    VectorNd q(2);
    VectorNd qdot(2);

    q << 0.1,-0.4;
    qdot << 0.9,-0.125;

    Utils::calcCenterOfMassJacobian(model,q,J_com);

    Utils::calcCenterOfMass(model,q,qdot,wholeRobotMass,com,&com_velocity);
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(J_com*qdot,com_velocity,unit_test_utils::TEST_PREC));
}

TEST_F(Human36, TestCOMJacobianHuman36)
{
    randomizeStates();

    double wholeRobotMass;
    Vector3d com,com_velocity;

    MatrixNd J_com(3,model_emulated->qdot_size);
    J_com.setZero();

    Utils::calcCenterOfMassJacobian(*model,q,J_com);
    Utils::calcCenterOfMass(*model,q,qdot,wholeRobotMass,com,&com_velocity);

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(J_com*qdot,com_velocity,unit_test_utils::TEST_PREC));

    J_com.setZero();

    Utils::calcCenterOfMassJacobian(*model_3dof,q,J_com);

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(J_com*qdot,com_velocity,unit_test_utils::TEST_PREC));
}

TEST_F(RdlUtilsTests, TestAngularMomentumSimple)
{
    Model model;
    Matrix3d inertia = Matrix3d::Zero(3, 3);
    inertia(0, 0) = 1.1;
    inertia(1, 1) = 2.2;
    inertia(2, 2) = 3.3;

    Body body(0.5, Vector3d(1., 0., 0.), inertia);
    Joint joint(SpatialVector(1., 0., 0., 0., 0., 0.), SpatialVector(0., 1., 0., 0., 0., 0.), SpatialVector(0., 0., 1., 0., 0., 0.));

    model.appendBody(Xtrans(Vector3d(0., 0., 0.)), joint, body);

    VectorNd q = VectorNd::Zero(model.q_size);
    VectorNd qdot = VectorNd::Zero(model.qdot_size);

    double mass;
    Vector3d com;
    Vector3d angular_momentum;

    FramePointd p_com;
    FrameVector v_com, ang_momentum;

    qdot << 1., 0., 0.;
    Utils::calcCenterOfMass(model, q, qdot, mass, com, NULL, &angular_momentum);
    Utils::calcCenterOfMass(model, q, qdot, mass, p_com, nullptr, &ang_momentum);
    EXPECT_EQ (Vector3d(1.1, 0., 0.), angular_momentum);
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(Vector3d(1.1, 0., 0.), ang_momentum,unit_test_utils::TEST_PREC));

    qdot << 0., 1., 0.;
    Utils::calcCenterOfMass(model, q, qdot, mass, com, NULL, &angular_momentum);
    Utils::calcCenterOfMass(model, q, qdot, mass, p_com, nullptr, &ang_momentum);
    EXPECT_EQ (Vector3d(0., 2.2, 0.), angular_momentum);
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(Vector3d(0., 2.2, 0.), ang_momentum,unit_test_utils::TEST_PREC));

    qdot << 0., 0., 1.;
    Utils::calcCenterOfMass(model, q, qdot, mass, com, NULL, &angular_momentum);
    Utils::calcCenterOfMass(model, q, qdot, mass, p_com, nullptr, &ang_momentum);
    EXPECT_EQ (Vector3d(0., 0., 3.3), angular_momentum);
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(Vector3d(0., 0., 3.3), ang_momentum,unit_test_utils::TEST_PREC));
}

TEST_F (TwoArms12DoF, TestAngularMomentumSimple)
{
    double mass;
    Vector3d com;
    Vector3d angular_momentum;

    Utils::calcCenterOfMass(*model, q, qdot, mass, com, NULL, &angular_momentum);

    EXPECT_EQ (Vector3d(0., 0., 0.), angular_momentum);

    qdot[0] = 1.;
    qdot[1] = 2.;
    qdot[2] = 3.;

    Utils::calcCenterOfMass(*model, q, qdot, mass, com, NULL, &angular_momentum);

    // only a rough guess from test calculation
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(Vector3d(3.3, 2.54, 1.5), angular_momentum, 1.0e-1));

    qdot[3] = -qdot[0];
    qdot[4] = -qdot[1];
    qdot[5] = -qdot[2];


    Utils::calcCenterOfMass(*model, q, qdot, mass, com, NULL, &angular_momentum);

    EXPECT_TRUE (angular_momentum[0] == 0);
    EXPECT_TRUE (angular_momentum[1] < 0);
    EXPECT_TRUE (angular_momentum[2] == 0.);
}

TEST_F (TwoArms12DoF, calcCentroidalMomentumMatrix)
{
    double mass;
    Vector3d com,com_velocity;
    Vector3d angular_momentum;
    MatrixNd A(6,model->qdot_size);
    A.setZero();

    randomizeStates();

    Utils::calcCenterOfMass(*model, q, qdot, mass, com, &com_velocity, &angular_momentum);
    Utils::calcCentroidalMomentumMatrix(*model,q,A);

    SpatialVector m_exp;
    m_exp.setLinearPart(com_velocity*mass);
    m_exp.setAngularPart(angular_momentum);

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(SpatialVector(A*qdot),m_exp,unit_test_utils::TEST_PREC*10.));
}

TEST_F (RdlUtilsTests, calcCentroidalMomentumMatrix)
{
    Model model;

    Body b1(1.,Vector3d(0.,0.,-0.1),Vector3d(1.,1.,1.));
    Joint j1(JointTypeRevoluteX);

    unsigned int id = model.addBody(0,SpatialTransform(),j1,b1,"b1");

    model.addBody(id,Xtrans(Vector3d(0.,0.,-1.)),j1,b1,"b2");

    VectorNd Q(model.q_size);
    VectorNd QDot(model.qdot_size);
    Q.setZero();
    QDot.setZero();
    QDot[0] = 0.1;
    QDot[1] = 0.1;
    MatrixNd A(6,model.qdot_size);
    A.setZero();

    double mass;
    Vector3d com,com_velocity,ang_momentum;

    Utils::calcCenterOfMass(model,Q,QDot,mass,com,&com_velocity,&ang_momentum);

    SpatialVector m_exp;
    m_exp.setAngularPart(ang_momentum);
    m_exp.setLinearPart(com_velocity*mass);

    Utils::calcCentroidalMomentumMatrix(model,Q,A,true);

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(SpatialVector(A*QDot),m_exp,unit_test_utils::TEST_PREC*10.));
}

TEST_F (FixedBase3DoFPlanar, calcCentroidalMomentumMatrix)
{
    double mass;
    Vector3d com,com_velocity, ang_momentum;

    randomizeStates();

    MatrixNd A(6,model->qdot_size),G(3,model->qdot_size);
    A.setZero();

    Utils::calcCenterOfMass(*model, Q, QDot, mass, com,&com_velocity,&ang_momentum);
    Utils::calcCentroidalMomentumMatrix(*model,Q,A);

    SpatialVector m_exp;
    m_exp.setLinearPart(com_velocity*mass);
    m_exp.setAngularPart(ang_momentum);

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(SpatialVector(A*QDot),m_exp,unit_test_utils::TEST_PREC*10.));
}

TEST_F (Human36, calcCentroidalMomentumMatrix)
{
    double mass;
    Vector3d com,com_velocity, ang_momentum;

    randomizeStates();

    MatrixNd A(6,model->qdot_size),G(3,model->qdot_size);
    A.setZero();

    Utils::calcCenterOfMass(*model, q, qdot, mass, com,&com_velocity,&ang_momentum);
    Utils::calcCentroidalMomentumMatrix(*model,q,A);

    SpatialVector m_exp;
    m_exp.setLinearPart(com_velocity*mass);
    m_exp.setAngularPart(ang_momentum);

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(SpatialVector(A*qdot),m_exp,unit_test_utils::TEST_PREC*10.));

    A.setZero();

    Utils::calcCenterOfMass(*model_3dof, q, qdot, mass, com,&com_velocity,&ang_momentum);
    Utils::calcCentroidalMomentumMatrix(*model_3dof,q,A);

    m_exp.setLinearPart(com_velocity*mass);
    m_exp.setAngularPart(ang_momentum);

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(SpatialVector(A*qdot),m_exp,unit_test_utils::TEST_PREC*10.));

    A.setZero();

    Utils::calcCenterOfMass(*model_emulated, q, qdot, mass, com,&com_velocity,&ang_momentum);
    Utils::calcCentroidalMomentumMatrix(*model_emulated,q,A);

    m_exp.setLinearPart(com_velocity*mass);
    m_exp.setAngularPart(ang_momentum);

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(SpatialVector(A*qdot),m_exp,unit_test_utils::TEST_PREC*10.));
}

TEST_F (RdlUtilsTests, calcCentroidalMomentumMatrixDot)
{
    Model model;
    model.gravity.setZero();

    Body b1(1.,Vector3d(0.,0.,-0.1),Vector3d(1.,1.,1.));
    Joint jx(JointTypeRevoluteX);

    model.addBody(0,SpatialTransform(),jx,b1,"b1");
    model.appendBody(Xtrans(Vector3d(0.,0.,-1.)),jx,b1,"b2");

    VectorNd Q(model.q_size);
    VectorNd QDot(model.qdot_size);
    VectorNd QDDot(model.qdot_size);
    VectorNd Tau(model.qdot_size);
    Q.setZero();
    QDot.setZero();
    QDDot.setZero();
    QDot[0] = 0.1;
    QDot[1] = 0.1;
    Tau[0] = 4;
    Tau[1] = 5;
    Tau.setZero();
    MatrixNd A(6,model.qdot_size), ADot(6,model.qdot_size), ADot_num(6,model.qdot_size);
    A.setZero();
    ADot.setZero();

    double mass;
    Vector3d com,com_velocity,ang_momentum;

    Utils::calcCenterOfMass(model,Q,QDot,mass,com,&com_velocity,&ang_momentum);

    SpatialVector m_exp;
    m_exp.setAngularPart(ang_momentum);
    m_exp.setLinearPart(com_velocity*mass);

    Utils::calcCentroidalMomentumMatrix(model,Q,A,true);

    double h = 0.00000005;
    Math::VectorNd x_euler = Math::VectorNd::Zero(model.q_size+model.qdot_size);
    x_euler.setZero();
    Math::VectorNd x_rk4 = x_euler;

    unit_test_utils::integrateRk4(model,Q,QDot,x_rk4,Tau,h);
    Q = x_rk4.head(model.q_size);
    QDot = x_rk4.tail(model.qdot_size);

    MatrixNd A2(6,model.qdot_size);
    A2.setZero();
    Utils::calcCentroidalMomentumMatrix(model,Q,A2,true);
    Utils::calcCentroidalMomentumMatrixDot(model,Q,QDot,ADot,true);

    ADot_num = ((1.0/h)*(A2-A));

    EXPECT_TRUE(unit_test_utils::checkVectorNdEpsilonClose(ADot_num*QDot,ADot*QDot,1e-6));
}

TEST_F(FixedBase3DoFPlanar, calcCentroidalMomentumMatrixDot)
{
    randomizeStates();
    double h = 0.00000005;

    RobotDynamics::Math::MatrixNd A1(6,model->qdot_size),A2(6,model->qdot_size),ADot(6,model->qdot_size),ADot_num(6,model->qdot_size);
    A1.setZero();
    A2.setZero();
    ADot.setZero();
    ADot_num.setZero();
    RobotDynamics::Utils::calcCentroidalMomentumMatrix(*model,Q,A1);

    Math::VectorNd x_rk4 = Math::VectorNd::Zero(model->q_size+model->qdot_size);
    unit_test_utils::integrateRk4(*model,Q,QDot,x_rk4,Tau,h);

    // integrate
    Q = x_rk4.head(model->q_size);
    QDot = x_rk4.tail(model->qdot_size);

    RobotDynamics::Utils::calcCentroidalMomentumMatrixDot(*model,Q,QDot,ADot);

    RobotDynamics::Utils::calcCentroidalMomentumMatrix(*model,Q,A2);

    ADot_num = (1.0/h)*(A2-A1);

    EXPECT_TRUE(unit_test_utils::checkVectorNdEpsilonClose(ADot_num*QDot,ADot*QDot,1e-4));
}

TEST_F(Human36, calcCentroidalMomentumMatrixDot)
{
    randomizeStates();
    double h = 0.00000005;

    RobotDynamics::Math::MatrixNd A1(6,model_emulated->qdot_size),A2(6,model_emulated->qdot_size),ADot(6,model_emulated->qdot_size),ADot_num(6,model_emulated->qdot_size);
    A1.setZero();
    A2.setZero();
    ADot.setZero();
    ADot_num.setZero();
    RobotDynamics::Utils::calcCentroidalMomentumMatrix(*model_emulated,q,A1);

    Math::VectorNd x_rk4 = Math::VectorNd::Zero(model_emulated->q_size+model_emulated->qdot_size);
    unit_test_utils::integrateRk4(*model_emulated,q,qdot,x_rk4,tau,h);
    // integra
    q = x_rk4.head(model_emulated->q_size);
    qdot = x_rk4.tail(model_emulated->qdot_size);

    RobotDynamics::Utils::calcCentroidalMomentumMatrixDot(*model_emulated,q,qdot,ADot);

    RobotDynamics::Utils::calcCentroidalMomentumMatrix(*model_emulated,q,A2);

    ADot_num = (1.0/h)*(A2-A1);

    EXPECT_TRUE(unit_test_utils::checkVectorNdEpsilonClose(ADot_num*qdot,ADot*qdot,1e-4));
}

TEST_F(RdlUtilsTests, calcCentroidalMomentumMatrixDotSphericalJoint)
{
    Model model;
    Body b(1.,Vector3d(1.,1.,1.),Vector3d(1.,1.,1.));
    model.appendBody(SpatialTransform(),Joint(SpatialVector(0.,0.,0.,1.,0.,0.)),b);
    model.appendBody(SpatialTransform(),Joint(SpatialVector(0.,0.,0.,0.,1.,0.)),b);
    model.appendBody(SpatialTransform(),Joint(SpatialVector(0.,0.,0.,0.,0.,1.)),b);
    unsigned int id = model.appendBody(SpatialTransform(),JointTypeSpherical,b);
    RobotDynamics::Math::Quaternion quat(1.,2.,3.,4.);
    quat.normalize();
    double h = 0.00000005;
    VectorNd q(model.q_size),qdot(model.qdot_size),tau(model.qdot_size);
    q.setZero();
    qdot.setZero();
    tau.setZero();
    RobotDynamics::Math::MatrixNd A1(6,model.qdot_size),A2(6,model.qdot_size),ADot(6,model.qdot_size),ADot_num(6,model.qdot_size);
    model.SetQuaternion(id,quat,q);
    qdot[0] = 1.;
    qdot[1] = -0.3;
    qdot[2] = 0.115;
    qdot[3] = 0.36;
    qdot[4] = -1.1;
    qdot[5] = 2.2;
    A1.setZero();
    A2.setZero();
    ADot.setZero();
    ADot_num.setZero();
    RobotDynamics::Utils::calcCentroidalMomentumMatrix(model,q,A1);

    Math::VectorNd x_rk4 = Math::VectorNd::Zero(model.q_size+model.qdot_size);
    x_rk4.setZero();
    unit_test_utils::integrateEuler(model,q,qdot,x_rk4,tau,h);

    q = x_rk4.head(model.q_size);
    qdot = x_rk4.tail(model.qdot_size);

    RobotDynamics::Utils::calcCentroidalMomentumMatrixDot(model,q,qdot,ADot);

    RobotDynamics::Utils::calcCentroidalMomentumMatrix(model,q,A2);

    ADot_num = (1.0/h)*(A2-A1);

    EXPECT_TRUE(unit_test_utils::checkVectorNdEpsilonClose(ADot_num*qdot,ADot*qdot,1e-4));
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}