//
// Created by jordan on 4/20/17.
//

#include <gtest/gtest.h>
#include <Eigen/Dense>
#include "rdl_dynamics/Quaternion.h"
#include "UnitTestUtils.hpp"

struct QuaternionFixture : public testing::Test
{
    QuaternionFixture()
    {
        f = Eigen::IOFormat(Eigen::FullPrecision);
    }

    Eigen::IOFormat f;
};

TEST_F(QuaternionFixture, EigenConversion)
{
    RobotDynamics::Math::Quaternion q;
    Eigen::Quaterniond q_eig(1.,2.,3.,4.);

    EXPECT_EQ(q.x(),0.);
    EXPECT_EQ(q.y(),0.);
    EXPECT_EQ(q.z(),0.);
    EXPECT_EQ(q.w(),1.);

    q = q_eig;

    EXPECT_EQ(q.x(),2.);
    EXPECT_EQ(q.y(),3.);
    EXPECT_EQ(q.z(),4.);
    EXPECT_EQ(q.w(),1.);
}

TEST_F(QuaternionFixture, element_accessors)
{
    RobotDynamics::Math::Quaternion q1(0.2, 0.3, 0.4, 0.1);

    q1.x() = 1.1;
    q1.y() = 2.1;
    q1.z() = 3.1;
    q1.w() = 4.1;

    EXPECT_EQ(1.1,q1.x());
    EXPECT_EQ(2.1,q1.y());
    EXPECT_EQ(3.1,q1.z());
    EXPECT_EQ(4.1,q1.w());

    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(q1.getVectorPart(),Vector3d(1.1,2.1,3.1),unit_test_utils::TEST_PREC));
    EXPECT_EQ(q1.getScalarPart(),4.1);
}

TEST_F(QuaternionFixture, multiplication)
{
    RobotDynamics::Math::Quaternion q1(0.2, 0.3, 0.4, 0.1), q2(0.2, 0.1, 0.5, 0.3);
    q1.normalize();
    q2.normalize();

    RobotDynamics::Math::Quaternion q_exp(0.5554700788944518, 0.2338821384818745, 0.3800584750330459, -0.7016464154456233);

    EXPECT_TRUE(unit_test_utils::checkVector4dEpsilonClose(q1 * q2, q_exp, unit_test_utils::TEST_PREC));

    q1 *= q2;

    EXPECT_TRUE(unit_test_utils::checkVector4dEpsilonClose(q1, q_exp, unit_test_utils::TEST_PREC));
}

TEST_F(QuaternionFixture, slerp)
{
    RobotDynamics::Math::Quaternion q1(0.2, 0.3, 0.4, 0.1), q2(0.2, 0.1, 0.5, 0.3);
    q1.normalize();
    q2.normalize();

    RobotDynamics::Math::Quaternion q_exp(0.3633161731320073, 0.4782736431579663, 0.7599826545340369, 0.2483587031060483);

    EXPECT_TRUE(unit_test_utils::checkVector4dEpsilonClose(q_exp, q1.slerp(0.2, q2), unit_test_utils::TEST_PREC));

    q1 = q1 * -1.;
    q2 = q2 * -1.;

    EXPECT_TRUE(unit_test_utils::checkVector4dEpsilonClose(q_exp * -1., q1.slerp(0.2, q2), unit_test_utils::TEST_PREC));

    q1 = q1 * -1.;
    q2 = q2 * -1.;

    EXPECT_TRUE(unit_test_utils::checkVector4dEpsilonClose(q2, q1.slerp(1., q2), unit_test_utils::TEST_PREC));
}

TEST_F(QuaternionFixture, fromAxisAngle)
{
    RobotDynamics::Math::Vector3d v(0.5, 0.11, -0.3);
    v.normalize();
    Eigen::AngleAxisd aa(0.5, v);
    Eigen::Quaterniond q2(aa);

    RobotDynamics::Math::Quaternion q_exp(0.2084700340099281, 0.04586340748218418, -0.1250820204059569, 0.9689124217106447);
    RobotDynamics::Math::Quaternion q = RobotDynamics::Math::Quaternion::fromAxisAngle(v, 0.5);

    EXPECT_TRUE(unit_test_utils::checkVector4dEpsilonClose(q_exp, q, unit_test_utils::TEST_PREC));
}

TEST_F(QuaternionFixture, fromMatrix)
{
    RobotDynamics::Math::Matrix3d m = (RobotDynamics::Math::Xrotx(0.1) * RobotDynamics::Math::Xroty(0.1) * RobotDynamics::Math::Xrotz(-0.1)).E;

    Eigen::Quaterniond q2(m);

    RobotDynamics::Math::Quaternion q = RobotDynamics::Math::Quaternion::fromMatrix(m);
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(m, q.toMatrix(), unit_test_utils::TEST_PREC));

    m << 0., 1., 0., 1., 0., 0., 0., 0., -1.;

    q = RobotDynamics::Math::Quaternion::fromMatrix(m);
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(m, q.toMatrix(), unit_test_utils::TEST_PREC));

    m << -1., 0., 0., 0., 1., 0., 0., 0., -1.;

    q = RobotDynamics::Math::Quaternion::fromMatrix(m);
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(m, q.toMatrix(), unit_test_utils::TEST_PREC));

    m << 1., 0., 0., 0., -1., 0., 0., 0., -1.;

    q = RobotDynamics::Math::Quaternion::fromMatrix(m);
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(m, q.toMatrix(), unit_test_utils::TEST_PREC));

    m << -1., 0., 0., 0., -1., 0., 0., 0., 1.;

    q = RobotDynamics::Math::Quaternion::fromMatrix(m);
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(m, q.toMatrix(), unit_test_utils::TEST_PREC));

    m = (RobotDynamics::Math::Xrotx(0.1144) * RobotDynamics::Math::Xroty(-1.99331) * RobotDynamics::Math::Xrotz(-0.98011)).E;
    RobotDynamics::Math::Quaternion q_exp(-0.367164735122942, -0.7542393999032798, -0.2128590176200619, 0.5010030174893784);
    EXPECT_TRUE(unit_test_utils::checkVector4dEpsilonClose(q_exp,RobotDynamics::Math::Quaternion::fromMatrix(m),unit_test_utils::TEST_PREC));
}

TEST_F(QuaternionFixture, conjugate)
{
    RobotDynamics::Math::Quaternion q(0.1,0.2,-0.3,-0.4);
    EXPECT_TRUE(unit_test_utils::checkVector4dEpsilonClose(RobotDynamics::Math::Quaternion(-0.1,-0.2,0.3,-0.4),q.conjugate(),unit_test_utils::TEST_PREC));
}

TEST_F(QuaternionFixture, rotate)
{
    RobotDynamics::Math::Quaternion q1(0.1,0.2,-0.3,-0.4),q2(-0.1,-0.3,-0.2,0.);
    RobotDynamics::Math::Vector3d v(-0.1,-0.3,-0.2);
    q1.normalize();

    RobotDynamics::Math::Quaternion q_exp = q1.conjugate()*q2*q1;
    RobotDynamics::Math::Vector3d v_out = q1.rotate(v);

    EXPECT_NEAR(v_out.x(),q_exp.x(),unit_test_utils::TEST_PREC);
    EXPECT_NEAR(v_out.y(),q_exp.y(),unit_test_utils::TEST_PREC);
    EXPECT_NEAR(v_out.z(),q_exp.z(),unit_test_utils::TEST_PREC);
}

TEST_F(QuaternionFixture, euler)
{
    double th = 0.3;
    RobotDynamics::Math::Quaternion qx(std::sin(th/2.),0.,0.,std::cos(th/2.)),qy(0.,std::sin(th/2.),0.,std::cos(th/2.)),qz(0.,0.,std::sin(th/2.),std::cos(th/2.));

    EXPECT_TRUE(unit_test_utils::checkVector4dEpsilonClose(RobotDynamics::Math::Quaternion::fromXYZAngles(RobotDynamics::Math::Vector3d(0.3,0.3,0.3)),qx*qy*qz,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector4dEpsilonClose(RobotDynamics::Math::Quaternion::fromZYXAngles(RobotDynamics::Math::Vector3d(0.3,0.3,0.3)),qz*qy*qx,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkVector4dEpsilonClose(RobotDynamics::Math::Quaternion::fromYXZAngles(RobotDynamics::Math::Vector3d(0.3,0.3,0.3)),qy*qx*qz,unit_test_utils::TEST_PREC));
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}