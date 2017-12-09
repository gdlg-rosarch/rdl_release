#include <gtest/gtest.h>

#include "rdl_dynamics/rdl_mathutils.h"
#include "rdl_dynamics/Body.h"
#include "UnitTestUtils.hpp"

using namespace std;
using namespace RobotDynamics;
using namespace RobotDynamics::Math;

class RdlBodyTests : public testing::Test
{
public:
    RdlBodyTests()
    {
    };

    void SetUp()
    {

    }

    void TearDown()
    {

    }
};

/* Tests whether the spatial inertia matches the one specified by its
 * parameters
 */
TEST_F(RdlBodyTests, TestComputeSpatialInertiaFromAbsoluteRadiiGyration)
{
    Body body(1.1, Vector3d(1.5, 1.2, 1.3), Vector3d(1.4, 2., 3.));

    Matrix3d inertia_C(1.4, 0., 0., 0., 2., 0., 0., 0., 3.);

    SpatialMatrix reference_inertia(4.843, -1.98, -2.145, 0, -1.43, 1.32, -1.98, 6.334, -1.716, 1.43, 0, -1.65, -2.145, -1.716, 7.059, -1.32, 1.65, 0, 0, 1.43, -1.32, 1.1, 0, 0, -1.43, 0, 1.65, 0, 1.1, 0, 1.32, -1.65, 0, 0, 0, 1.1);

    RigidBodyInertia rbi = createFromMassComInertiaC(body.mMass,body.mCenterOfMass,body.mInertia);

    EXPECT_TRUE(unit_test_utils::checkSpatialMatrixEpsilonClose(reference_inertia, rbi.toMatrix(), unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(inertia_C, body.mInertia, unit_test_utils::TEST_PREC));
}

TEST_F(RdlBodyTests, TestJoinTwoBodiesWithZeroMassThrows)
{
    Body body1(0., Vector3d(1.5, 1.2, 1.3), Vector3d(1.4, 2., 3.));
    Body body2(0., Vector3d(1.5, 1.2, 1.3), Vector3d(1.4, 2., 3.));

    try
    {
        body1.join(RobotDynamics::Math::SpatialTransform(),body2);
        FAIL();
    }
    catch(RobotDynamics::RdlException &e)
    {
        EXPECT_STREQ(e.what(),"Error: cannot join bodies as both have zero mass!");
    }
}

TEST_F(RdlBodyTests, TestBodyConstructorMassComInertia)
{
    double mass = 1.1;
    Vector3d com(1.5, 1.2, 1.3);
    Matrix3d inertia_C(8.286, -3.96, -4.29, -3.96, 10.668, -3.432, -4.29, -3.432, 11.118);

    Body body(mass, com, inertia_C);

    SpatialMatrix reference_inertia(11.729, -5.94, -6.435, 0, -1.43, 1.32, -5.94, 15.002, -5.148, 1.43, 0, -1.65, -6.435, -5.148, 15.177, -1.32, 1.65, 0, 0, 1.43, -1.32, 1.1, 0, 0, -1.43, 0, 1.65, 0, 1.1, 0, 1.32, -1.65, 0, 0, 0, 1.1);

    RigidBodyInertia body_rbi = createFromMassComInertiaC(body.mMass, body.mCenterOfMass, body.mInertia);
    EXPECT_TRUE(unit_test_utils::checkSpatialMatrixEpsilonClose(reference_inertia, body_rbi.toMatrix(), unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkMatrix3dEpsilonClose(inertia_C, body.mInertia, unit_test_utils::TEST_PREC));
}

TEST_F(RdlBodyTests, TestBodyJoinNullbody)
{
    Body body(1.1, Vector3d(1.5, 1.2, 1.3), Vector3d(1.4, 2., 3.));
    Body nullbody(0., Vector3d(0., 0., 0.), Vector3d(0., 0., 0.));

    Body joined_body = body;
    joined_body.join(Xtrans(Vector3d(0., 0., 0.)), nullbody);

    RigidBodyInertia body_rbi(body.mMass, body.mCenterOfMass, body.mInertia);
    RigidBodyInertia joined_body_rbi(joined_body.mMass, joined_body.mCenterOfMass, joined_body.mInertia);

    ASSERT_EQ(1.1, body.mMass);
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(body.mCenterOfMass, joined_body.mCenterOfMass, unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkSpatialMatrixEpsilonClose(body_rbi.toMatrix(), joined_body_rbi.toMatrix(), unit_test_utils::TEST_PREC));
}

TEST_F(RdlBodyTests, TestBodyJoinTwoBodies)
{
    Body body_a(1.1, Vector3d(-1.1, 1.3, 0.), Vector3d(3.1, 3.2, 3.3));
    Body body_b(1.1, Vector3d(1.1, 1.3, 0.), Vector3d(3.1, 3.2, 3.3));

    Body body_joined(body_a);
    body_joined.join(Xtrans(Vector3d(0., 0., 0.)), body_b);

    RigidBodyInertia body_joined_rbi = createFromMassComInertiaC(body_joined.mMass, body_joined.mCenterOfMass, body_joined.mInertia);

    SpatialMatrix reference_inertia(9.918, 0, 0, 0, -0, 2.86, 0, 9.062, 0, 0, 0, -0, 0, 0, 12.98, -2.86, 0, 0, 0, 0, -2.86, 2.2, 0, 0, -0, 0, 0, 0, 2.2, 0, 2.86, -0, 0, 0, 0, 2.2);

    ASSERT_EQ(2.2, body_joined.mMass);
    EXPECT_TRUE(unit_test_utils::checkVector3dEq(Vector3d(0., 1.3, 0.), body_joined.mCenterOfMass));
    EXPECT_TRUE(unit_test_utils::checkSpatialMatrixEpsilonClose(reference_inertia, body_joined_rbi.toMatrix(), unit_test_utils::TEST_PREC));
}

TEST_F(RdlBodyTests, TestBodyJoinTwoBodiesDisplaced)
{
    Body body_a(1.1, Vector3d(-1.1, 1.3, 0.), Vector3d(3.1, 3.2, 3.3));
    Body body_b(1.1, Vector3d(0., 0., 0.), Vector3d(3.1, 3.2, 3.3));

    Body body_joined(body_a);
    body_joined.join(Xtrans(Vector3d(1.1, 1.3, 0.)), body_b);

    RigidBodyInertia body_joined_rbi = createFromMassComInertiaC(body_joined.mMass, body_joined.mCenterOfMass, body_joined.mInertia);

    SpatialMatrix reference_inertia(9.918, 0, 0, 0, -0, 2.86, 0, 9.062, 0, 0, 0, -0, 0, 0, 12.98, -2.86, 0, 0, 0, 0, -2.86, 2.2, 0, 0, -0, 0, 0, 0, 2.2, 0, 2.86, -0, 0, 0, 0, 2.2);

    ASSERT_EQ(2.2, body_joined.mMass);
    EXPECT_TRUE(unit_test_utils::checkVector3dEq(Vector3d(0., 1.3, 0.), body_joined.mCenterOfMass));
    EXPECT_TRUE(unit_test_utils::checkSpatialMatrixEpsilonClose(reference_inertia, body_joined_rbi.toMatrix(), unit_test_utils::TEST_PREC));


}

TEST_F(RdlBodyTests, TestBodyJoinTwoBodiesRotated)
{
    Body body_a(1.1, Vector3d(0., 0., 0.), Vector3d(3.1, 3.2, 3.3));
    Body body_b(1.1, Vector3d(0., 0., 0.), Vector3d(3.1, 3.3, 3.2));

    Body body_joined(body_a);
    body_joined.join(Xrotx(-M_PI * 0.5), body_b);

    RigidBodyInertia body_joined_rbi(body_joined.mMass, body_joined.mCenterOfMass, body_joined.mInertia);

    SpatialMatrix reference_inertia(6.2, 0., 0., 0., 0., 0., 0., 6.4, 0., 0., 0., 0., 0., 0., 6.6, 0., 0., 0., 0., 0., 0., 2.2, 0., 0., 0., 0., 0., 0., 2.2, 0., 0., 0., 0., 0., 0., 2.2);

    ASSERT_EQ(2.2, body_joined.mMass);
    EXPECT_TRUE(unit_test_utils::checkVector3dEq(Vector3d(0., 0., 0.), body_joined.mCenterOfMass));
    EXPECT_TRUE(unit_test_utils::checkSpatialMatrixEpsilonClose(reference_inertia, body_joined_rbi.toMatrix(), unit_test_utils::TEST_PREC));
}

TEST_F(RdlBodyTests, TestBodyJoinTwoBodiesRotatedAndTranslated)
{
    Body body_a(1.1, Vector3d(0., 0., 0.), Vector3d(3.1, 3.2, 3.3));
    Body body_b(1.1, Vector3d(-1., 1., 0.), Vector3d(3.2, 3.1, 3.3));

    Body body_joined(body_a);
    body_joined.join(Xrotz(M_PI * 0.5) * Xtrans(Vector3d(1., 1., 0.)), body_b);

    RigidBodyInertia body_joined_rbi(body_joined.mMass, body_joined.mCenterOfMass, body_joined.mInertia);

    SpatialMatrix reference_inertia(6.2, 0., 0., 0., 0., 0., 0., 6.4, 0., 0., 0., 0., 0., 0., 6.6, 0., 0., 0., 0., 0., 0., 2.2, 0., 0., 0., 0., 0., 0., 2.2, 0., 0., 0., 0., 0., 0., 2.2);

    ASSERT_EQ (2.2, body_joined.mMass);
    EXPECT_TRUE(unit_test_utils::checkVector3dEpsilonClose(Vector3d(0., 0., 0.), body_joined.mCenterOfMass, unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkSpatialMatrixEpsilonClose(reference_inertia, body_joined_rbi.toMatrix(), unit_test_utils::TEST_PREC));
}

TEST_F(RdlBodyTests, TestBodyConstructorRigidBodyInertiaMultiplyMotion)
{
    Body body(1.1, Vector3d(1.5, 1.2, 1.3), Vector3d(1.4, 2., 3.));

    RigidBodyInertia rbi = RigidBodyInertia(body.mMass, body.mCenterOfMass * body.mMass, body.mInertia);

    SpatialVector mv(1.1, 1.2, 1.3, 1.4, 1.5, 1.6);
    SpatialVector fv_matrix = rbi.toMatrix() * mv;
    SpatialVector fv_rbi = rbi * mv;

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(fv_matrix, fv_rbi, unit_test_utils::TEST_PREC));
}

TEST_F(RdlBodyTests, TestBodyConstructorRigidBodyInertia)
{
    Body body(1.1, Vector3d(1.5, 1.2, 1.3), Vector3d(1.4, 2., 3.));

    RigidBodyInertia rbi = RigidBodyInertia(body.mMass, body.mCenterOfMass * body.mMass, body.mInertia);
    SpatialMatrix spatial_inertia = rbi.toMatrix();

    EXPECT_TRUE(unit_test_utils::checkSpatialMatrixEpsilonClose(spatial_inertia, rbi.toMatrix(), unit_test_utils::TEST_PREC));
}

TEST_F(RdlBodyTests, TestBodyConstructorCopyRigidBodyInertia)
{
    Body body(1.1, Vector3d(1.5, 1.2, 1.3), Vector3d(1.4, 2., 3.));

    RigidBodyInertia rbi = RigidBodyInertia(body.mMass, body.mCenterOfMass * body.mMass, body.mInertia);

    RigidBodyInertia rbi_from_matrix;
    rbi_from_matrix.createFromMatrix(rbi.toMatrix());

    EXPECT_TRUE(unit_test_utils::checkSpatialMatrixEpsilonClose(rbi.toMatrix(), rbi_from_matrix.toMatrix(), unit_test_utils::TEST_PREC));
}


int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}