//
// Created by jordan on 9/25/16.
//

#include "UnitTestUtils.hpp"
#include <gtest/gtest.h>

using namespace RobotDynamics::Math;

class ForceVectorTests : public testing::Test
{
public:
    ForceVectorTests()
    {

    }

    void SetUp()
    {

    }

    void TearDown()
    {

    }
};

TEST_F(ForceVectorTests, testConstructorsAndAccessors)
{
    ForceVector f(1.1,2.2,3.3,4.4,5.5,6.6);

    EXPECT_EQ(f.mx(),1.1);
    EXPECT_EQ(f.my(),2.2);
    EXPECT_EQ(f.mz(),3.3);
    EXPECT_EQ(f.fx(),4.4);
    EXPECT_EQ(f.fy(),5.5);
    EXPECT_EQ(f.fz(),6.6);

    f.mx() = -1.1;
    f.my() = -2.2;
    f.mz() = -3.3;
    f.fx() = -4.4;
    f.fy() = -5.5;
    f.fz() = -6.6;

    EXPECT_EQ(f.mx(),-1.1);
    EXPECT_EQ(f.my(),-2.2);
    EXPECT_EQ(f.mz(),-3.3);
    EXPECT_EQ(f.fx(),-4.4);
    EXPECT_EQ(f.fy(),-5.5);
    EXPECT_EQ(f.fz(),-6.6);
}

TEST_F(ForceVectorTests, testTransform)
{
    SpatialVector v(1.1, 2.1, 3.1, 4.1, 5.1, 6.1);
    ForceVector f(v);
    ForceVector f2(v),f3(v);

    SpatialTransform X = Xrotx(0.1)*Xtrans(Vector3d(0.1,0.2,-0.3));
    f2 = X*f2;

    ForceVector f_tr = f.transform_copy(X);
    ForceVector f_exp;

    f_exp=X.toMatrixAdjoint()*ForceVector(v);
    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(f_exp,f_tr,unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(f_exp,f2,unit_test_utils::TEST_PREC));
}

TEST_F(ForceVectorTests, testTransformTranspose)
{
    Vector3d rot(1.1, 1.2, 1.3);
    Vector3d trans(1.1, 1.2, 1.3);

    SpatialTransform X_st;
    X_st.r = trans;

    SpatialMatrix X_66_matrix(SpatialMatrix::Zero(6, 6));
    X_66_matrix = Xrotz_mat(rot[2]) * Xroty_mat(rot[1]) * Xrotx_mat(rot[0]) * Xtrans_mat(trans);
    X_st.E = X_66_matrix.block<3, 3>(0, 0);

    SpatialVector v(1.1, 2.1, 3.1, 4.1, 5.1, 6.1);
    ForceVector f(v);
    ForceVector f2(v);
    SpatialVector v_66_res = X_66_matrix.transpose() * v;
    f.transformTranspose(X_st);
    f2 = X_st.toMatrixTranspose()*f2;

    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(v_66_res, f, unit_test_utils::TEST_PREC));
    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(v_66_res, f2, unit_test_utils::TEST_PREC));
}

TEST_F(ForceVectorTests, testOperatorOverloads)
{
    ForceVector f1(0.1,0.2,-0.3,-0.1,1.0,-3.);
    ForceVector f2(-0.1,-0.2,0.3,0.1,-1.0,3.);

    f1+=f2;
    EXPECT_TRUE(unit_test_utils::checkSpatialVectorsEpsilonClose(f1,SpatialVectorZero,unit_test_utils::TEST_PREC));
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}