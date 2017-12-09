//
// Created by jordan on 5/18/17.
//

#include <gtest/gtest.h>
#include <ros/package.h>

#include <rdl_dynamics/Dynamics.h>
#include <rdl_dynamics/Kinematics.h>
#include <rdl_dynamics/Model.h>
#include <rdl_urdfreader/urdfreader.h>

class UrdfReaderTests : public testing::Test
{
public:
    UrdfReaderTests()
    {
        srand(time(nullptr));
    }

    void SetUp()
    {

    }

    void TearDown()
    {

    }

    bool checkSpatialMatrixEpsilonClose(const RobotDynamics::Math::SpatialMatrix &t1,const RobotDynamics::Math::SpatialMatrix &t2,const double epsilon)
    {
        for(int i = 0;i<6;i++)
        {
            for(int j = 0; j<6; j++)
            {
                if(fabs(t1(i,j)-t2(i,j))>epsilon)
                {
                    return false;
                }
            }
        }

        return true;
    }

    bool checkSpatialVectorsEpsilonClose(const RobotDynamics::Math::SpatialVector &v1,const RobotDynamics::Math::SpatialVector &v2,const double epsilon)
    {
        for(int i = 0;i<6;i++)
        {
            if(fabs(v1(i)-v2(i))>epsilon)
            {
                return false;
            }
        }

        return true;
    }
};

TEST_F(UrdfReaderTests, testFixedArm)
{
    std::string file_path = ros::package::getPath("rdl_urdfreader") + "/tests/urdf/test_robot_arm.urdf";
    RobotDynamics::Model model;
    EXPECT_TRUE(RobotDynamics::Urdf::urdfReadFromFile(file_path.c_str(),&model,false,false));

    // First body in URDf is a fixed joint body, so it'll get mashed together with the root body
    EXPECT_EQ(model.mBodies[0].mMass,4.);
    EXPECT_TRUE(model.mBodies[0].mCenterOfMass.isApprox(RobotDynamics::Math::Vector3d(0., 0., 0.),1e-14));

    EXPECT_EQ(model.Ib_c[0].Ixx,0.0061063308908);
    EXPECT_EQ(model.Ib_c[0].Iyx,0.0);
    EXPECT_EQ(model.Ib_c[0].Izx,0.0);
    EXPECT_EQ(model.Ib_c[0].Iyy,0.0061063308908);
    EXPECT_EQ(model.Ib_c[0].Izy,0.0);
    EXPECT_EQ(model.Ib_c[0].Izz,0.01125);

    unsigned int id = model.GetBodyId("test_robot_shoulder_link");
    EXPECT_EQ(model.mBodies[id].mMass,7.778);
    EXPECT_TRUE(model.mBodies[id].mCenterOfMass.isApprox(RobotDynamics::Math::Vector3d(0., 0.01, 0.),1e-14));

    EXPECT_EQ(model.Ib_c[id].Ixx,0.0314743125769);
    EXPECT_EQ(model.Ib_c[id].Iyx,0.);
    EXPECT_EQ(model.Ib_c[id].Izx,0.);
    EXPECT_EQ(model.Ib_c[id].Iyy,0.0314743125769);
    EXPECT_EQ(model.Ib_c[id].Izy,0.);
    EXPECT_EQ(model.Ib_c[id].Izz,0.021875625);

    id = model.GetParentBodyId(model.GetBodyId("gripper_right_finger_link"));
    EXPECT_EQ(model.GetBodyId("gripper_right_knuckle_link"),id);

    EXPECT_EQ(model.mJoints[1].mJointType,RobotDynamics::JointTypePrismatic);
    EXPECT_EQ(model.GetBodyId("test_robot_shoulder_link"),model.lambda[model.GetBodyId("test_robot_upper_arm_link")]);
}

TEST_F(UrdfReaderTests, testNegative1DofJointAxes)
{
    std::string file_path = ros::package::getPath("rdl_urdfreader") + "/tests/urdf/test_robot_arm.urdf";
    std::string file_path_2 = ros::package::getPath("rdl_urdfreader") + "/tests/urdf/test_robot_arm_neg_jnt_axes.urdf";
    RobotDynamics::Model model,model_neg_axes;
    EXPECT_TRUE(RobotDynamics::Urdf::urdfReadFromFile(file_path.c_str(),&model,false,false));
    EXPECT_TRUE(RobotDynamics::Urdf::urdfReadFromFile(file_path_2.c_str(),&model_neg_axes,false,false));

    RobotDynamics::Math::VectorNd q = RobotDynamics::Math::VectorNd::Zero(model.q_size);
    RobotDynamics::Math::VectorNd q_neg = RobotDynamics::Math::VectorNd::Zero(model.q_size);
    RobotDynamics::Math::VectorNd qd = RobotDynamics::Math::VectorNd::Zero(model.qdot_size);
    RobotDynamics::Math::VectorNd qd_neg = RobotDynamics::Math::VectorNd::Zero(model.qdot_size);
    RobotDynamics::Math::VectorNd qdd = RobotDynamics::Math::VectorNd::Zero(model.qdot_size);
    RobotDynamics::Math::VectorNd qdd_neg = RobotDynamics::Math::VectorNd::Zero(model.qdot_size);
    RobotDynamics::Math::VectorNd tau = RobotDynamics::Math::VectorNd::Zero(model.qdot_size);
    RobotDynamics::Math::VectorNd tau_neg = RobotDynamics::Math::VectorNd::Zero(model.qdot_size);

    for (int i = 0; i < q.size(); i++)
    {
        q[i] = 0.4 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
        qd[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
        tau[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
        qdd[i] = 0.5 * M_PI * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    }

    q_neg = -q;
    qd_neg = -qd;
    qdd_neg = -qdd;
    tau_neg = -tau;

    RobotDynamics::updateKinematics(model,q,qd,qdd);
    RobotDynamics::updateKinematics(model_neg_axes,q_neg,qd_neg,qdd_neg);

    for(unsigned int i = 0; i<model.mBodies.size(); i++)
    {
        EXPECT_TRUE(checkSpatialVectorsEpsilonClose(model.v[i],model_neg_axes.v[i],1e-14));
        EXPECT_TRUE(checkSpatialVectorsEpsilonClose(model.a[i],model_neg_axes.a[i],1e-14));
    }

    RobotDynamics::forwardDynamics(model,q,qd,tau,qdd);
    RobotDynamics::forwardDynamics(model_neg_axes,q_neg,qd_neg,tau_neg,qdd_neg);

    EXPECT_TRUE(qdd.isApprox(-qdd_neg,1e-14));

    RobotDynamics::inverseDynamics(model,q,qd,qdd,tau);
    RobotDynamics::inverseDynamics(model_neg_axes,q_neg,qd_neg,qdd_neg,tau_neg);

    EXPECT_TRUE(tau.isApprox(-tau_neg,1e-14));
}

TEST_F(UrdfReaderTests, testJointBodyMap)
{
    std::string file_path = ros::package::getPath("rdl_urdfreader") + "/tests/urdf/test_robot_arm.urdf";
    std::map<std::string,std::string> jointBodyMap;

    ASSERT_TRUE(RobotDynamics::Urdf::parseJointBodyNameMapFromFile(file_path.c_str(),jointBodyMap));

    EXPECT_STREQ(jointBodyMap["test_robot_elbow_joint"].c_str(),"test_robot_forearm_link");
    EXPECT_TRUE(jointBodyMap.find("test_robot_ee_fixed_joint")==jointBodyMap.end()); // It's fixed, so shouldn't be here
    EXPECT_STREQ(jointBodyMap["test_robot_shoulder_pan_joint"].c_str(),"test_robot_shoulder_link");
}

TEST_F(UrdfReaderTests, testFloatingBaseRobot)
{
    RobotDynamics::Model model;
    std::string path_to_urdf = ros::package::getPath("rdl_urdfreader") + "/tests/urdf/floating_base_robot.urdf";
    EXPECT_TRUE(RobotDynamics::Urdf::urdfReadFromFile(path_to_urdf.c_str(),&model,true,false));

    EXPECT_NEAR(model.mBodies[2].mMass,4.,1e-14);
    EXPECT_TRUE(model.mBodies[2].mCenterOfMass.isApprox(RobotDynamics::Math::Vector3d(0.3, 0.2, 0.1),1e-14));

    EXPECT_NEAR(model.Ib_c[2].Ixx,0.1,1e-14);
    EXPECT_NEAR(model.Ib_c[2].Iyx,0.,1e-14);
    EXPECT_NEAR(model.Ib_c[2].Izx,0.,1e-14);
    EXPECT_NEAR(model.Ib_c[2].Iyy,0.2,1e-14);
    EXPECT_NEAR(model.Ib_c[2].Izy,0.,1e-14);
    EXPECT_NEAR(model.Ib_c[2].Izz,0.3,1e-14);

    EXPECT_EQ(model.mJoints[1].mJointType,RobotDynamics::JointTypeTranslationXYZ);
    EXPECT_EQ(model.mJoints[2].mJointType,RobotDynamics::JointTypeSpherical);
}

TEST_F(UrdfReaderTests, testRobotWithFloatingJoint)
{
    RobotDynamics::Model model;
    std::string path_to_urdf = ros::package::getPath("rdl_urdfreader") + "/tests/urdf/floating_joint_robot.urdf";
    EXPECT_TRUE(RobotDynamics::Urdf::urdfReadFromFile(path_to_urdf.c_str(),&model,true,false));

    EXPECT_EQ(model.mJoints[1].mJointType,RobotDynamics::JointTypeTranslationXYZ);
    EXPECT_EQ(model.mJoints[2].mJointType,RobotDynamics::JointTypeEulerXYZ);
}

TEST_F(UrdfReaderTests, testRobotSingleBodyFloatingBase)
{
    RobotDynamics::Model model;
    std::string path_to_urdf = ros::package::getPath("rdl_urdfreader") + "/tests/urdf/single_body_floating_base.urdf";
    ASSERT_TRUE(RobotDynamics::Urdf::urdfReadFromFile(path_to_urdf.c_str(),&model,true,false));
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    return RUN_ALL_TESTS();
}