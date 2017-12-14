
#include <gtest/gtest.h>
#include "rdl_dynamics/Momentum.hpp"
#include "UnitTestUtils.hpp"

using namespace RobotDynamics::Math;

class MomentumTests : public testing::Test
{
public:
    MomentumTests(){};

    void SetUp()
    {

    }

    void TearDown()
    {

    }
};

TEST_F(MomentumTests, testComparisonSpatialVectorCross)
{
    SpatialVector v(-1.,0.11,0.9,2.,0.56,-0.33);
    MotionVector m_v(v);
    RigidBodyInertia I(1.2,Vector3d(0.2,0.3,1.),1.,2.,3.,4.,5.,6.);
    RigidBodyInertia I_r(1.2,Vector3d(0.2,0.3,1.),1.,2.,3.,4.,5.,6.);
    Momentum m(I_r,m_v);
    MotionVector expected = m_v%m;
    SpatialVector v_sp_exp = crossf(v,m);

    SpatialVector v_out;
    v_out.setAngularPart(toTildeForm(v.getAngularPart())*m.getAngularPart()+toTildeForm(v.getLinearPart())*m.getLinearPart());//crossf(v,I * v);
    v_out.setLinearPart(toTildeForm(v.getAngularPart())*m.getLinearPart());

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(v_out,expected,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(v_sp_exp,expected,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(crossf(v)*m,expected,unit_test_utils::TEST_PREC));

    MotionVector m_2(0.1,-0.2,0.3,-0.4,-05,-0.6);

    expected = m%m_2;

    v_sp_exp = crossm(m,m_2);

    v_out.setAngularPart(toTildeForm(m.getAngularPart())*m_2.getAngularPart());
    v_out.setLinearPart(toTildeForm(m.getLinearPart())*m_2.getAngularPart() + toTildeForm(m.getAngularPart())*m_2.getLinearPart());

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(v_out,expected,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(v_sp_exp,expected,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(crossm(m)*m_2,expected,unit_test_utils::TEST_PREC));
}

TEST_F(MomentumTests, testConstructorsAndCompute)
{
    MotionVector v(1., 2., 3., 4., 5., 6.);
    RigidBodyInertia I(1., Vector3d(0.1, 1.1, 2.1), 0.1, 0.2, 0.3, 0.4, 0.5, 0.6);

    Momentum m = I*v;

    SpatialVector m_v = I.toMatrix() * v;

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(m,m_v,unit_test_utils::TEST_PREC));

    SpatialVector sv(-1.2,2.2,3.3,4.1,6.2,-0.8);

    Momentum mv(sv);

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(mv,sv,unit_test_utils::TEST_PREC));

    MotionVector mv2(sv);

    double kineticEnergy = m*mv2;

    EXPECT_EQ(kineticEnergy,m.dot(mv2)/2.);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}