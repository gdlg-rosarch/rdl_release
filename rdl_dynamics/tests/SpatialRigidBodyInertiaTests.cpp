//
// Created by jordan on 9/28/16.
//

#include <gtest/gtest.h>
#include "rdl_dynamics/rdl_math.hpp"
#include "UnitTestUtils.hpp"
#include "Fixtures.h"

using namespace RobotDynamics::Math;

class SpatialRigidBodyInertiaTests : public testing::Test
{
public:
	SpatialRigidBodyInertiaTests()
	{

	}

	void SetUp()
	{

	}

	void TearDown()
	{

	}
};

TEST_F(FixedBase3DoFPlanar, testAdd)
{
    SpatialInertia Ibody = model->I[1];
    SpatialInertia I_2(model->bodyFrames[1].get(),3.2,Vector3d(-0.4,1.,0.1),1.2,0.1,0.2,0.3,0.4,0.5);

    SpatialMatrix I_expected = Ibody.toMatrix() + I_2.toMatrix();

    SpatialInertia I_added = Ibody + I_2;

    EXPECT_TRUE(unit_test_utils::checkSpatialMatrixEpsilonClose(I_added.toMatrix(),I_expected,unit_test_utils::TEST_PREC));
}

TEST_F(FixedBase3DoFPlanar, ChangeFrameTest)
{
    SpatialInertia Ibody = model->I[1];
    SpatialInertia Ibody_c = model->Ib_c[1];

    Ibody.changeFrame(model->bodyCenteredFrames[1].get());

    EXPECT_TRUE(unit_test_utils::checkSpatialMatrixEpsilonClose(Ibody.toMatrix(),Ibody_c.toMatrix(),unit_test_utils::TEST_PREC));

    Ibody = model->I[3];
    Ibody_c = model->Ib_c[3];

    Ibody_c.changeFrame(model->bodyFrames[3].get());

    EXPECT_TRUE(unit_test_utils::checkSpatialMatrixEpsilonClose(Ibody.toMatrix(),Ibody_c.toMatrix(),unit_test_utils::TEST_PREC));
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}