/*
 * Original Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 *
 * RDL - Robot Dynamics Library
 * Modifications Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include "rdl_dynamics/rdl_utils.h"
#include "rdl_dynamics/Kinematics.h"

#include <iomanip>

namespace RobotDynamics
{

    namespace Utils
    {

        using namespace std;
        using namespace Math;

        string getDofName(const SpatialVector &joint_dof)
        {
            if (joint_dof == SpatialVector(1., 0., 0., 0., 0., 0.))
            {
                return "RX";
            }
            else if (joint_dof == SpatialVector(0., 1., 0., 0., 0., 0.))
            {
                return "RY";
            }
            else if (joint_dof == SpatialVector(0., 0., 1., 0., 0., 0.))
            {
                return "RZ";
            }
            else if (joint_dof == SpatialVector(0., 0., 0., 1., 0., 0.))
            {
                return "TX";
            }
            else if (joint_dof == SpatialVector(0., 0., 0., 0., 1., 0.))
            {
                return "TY";
            }
            else if (joint_dof == SpatialVector(0., 0., 0., 0., 0., 1.))
            {
                return "TZ";
            }

            ostringstream dof_stream(ostringstream::out);
            dof_stream << "custom (" << joint_dof.transpose() << ")";
            return dof_stream.str();
        }

        string getBodyName(const RobotDynamics::Model &model, unsigned int body_id)
        {
            if (model.mBodies[body_id].mIsVirtual)
            {
                // if there is not a unique child we do not know what to do...
                if (model.mu[body_id].size() != 1)
                {
                    return "";
                }

                return getBodyName(model, model.mu[body_id][0]);
            }

            return model.GetBodyName(body_id);
        }

        std::string getModelDOFOverview(const Model &model)
        {
            stringstream result("");

            unsigned int q_index = 0;
            for (unsigned int i = 1; i < model.mBodies.size(); i++)
            {
                if (model.mJoints[i].mDoFCount == 1)
                {
                    result << setfill(' ') << setw(3) << q_index << ": " << getBodyName(model, i) << "_"
                           << getDofName(model.S[i]) << endl;
                    q_index++;
                }
                else
                {
                    for (unsigned int j = 0; j < model.mJoints[i].mDoFCount; j++)
                    {
                        result << setfill(' ') << setw(3) << q_index << ": " << getBodyName(model, i) << "_"
                               << getDofName(model.mJoints[i].mJointAxes[j]) << endl;
                        q_index++;
                    }
                }
            }

            return result.str();
        }

        std::string printHierarchy(const RobotDynamics::Model &model, unsigned int body_index = 0, int indent = 0)
        {
            stringstream result("");

            for (int j = 0; j < indent; j++)
            {
                result << "  ";
            }

            result << getBodyName(model, body_index);

            if (body_index > 0)
            {
                result << " [ ";
            }

            while (model.mBodies[body_index].mIsVirtual)
            {
                if (model.mu[body_index].size() == 0)
                {
                    result << " end";
                    break;
                }
                else if (model.mu[body_index].size() > 1)
                {
                    cerr << endl << "Error: Cannot determine multi-dof joint as massless body with id " << body_index
                         << " (name: " << model.GetBodyName(body_index) << ") has more than one child:" << endl;
                    for (unsigned int ci = 0; ci < model.mu[body_index].size(); ci++)
                    {
                        cerr << "  id: " << model.mu[body_index][ci] << " name: "
                             << model.GetBodyName(model.mu[body_index][ci]) << endl;
                    }
                    abort();
                }

                result << getDofName(model.S[body_index]) << ", ";

                body_index = model.mu[body_index][0];
            }

            if (body_index > 0)
            {
                result << getDofName(model.S[body_index]) << " ]";
            }
            result << endl;

            unsigned int child_index = 0;
            for (child_index = 0; child_index < model.mu[body_index].size(); child_index++)
            {
                result << printHierarchy(model, model.mu[body_index][child_index], indent + 1);
            }

            // print fixed children
            for (unsigned int fbody_index = 0; fbody_index < model.mFixedBodies.size(); fbody_index++)
            {
                if (model.mFixedBodies[fbody_index].mMovableParent == body_index)
                {
                    for (int j = 0; j < indent + 1; j++)
                    {
                        result << "  ";
                    }

                    result << model.GetBodyName(model.fixed_body_discriminator + fbody_index) << " [fixed]" << endl;
                }
            }


            return result.str();
        }

        std::string getModelHierarchy(const Model &model)
        {
            stringstream result("");

            result << printHierarchy(model);

            return result.str();
        }

        std::string getNamedBodyOriginsOverview(Model &model)
        {
            stringstream result("");

            VectorNd Q(VectorNd::Zero(model.dof_count));
            updateKinematicsCustom(model, &Q, NULL, NULL);

            for (unsigned int body_id = 0; body_id < model.mBodies.size(); body_id++)
            {
                std::string body_name = model.GetBodyName(body_id);

                if (body_name.size() == 0)
                {
                    continue;
                }

                Vector3d position = model.bodyFrames[body_id]->getInverseTransformToRoot().r;

                result << body_name << ": " << position.transpose() << endl;
            }

            return result.str();
        }

        void calcCenterOfMass(Model               & model,
                              const Math::VectorNd& q,
                              Math::Vector3d      & com,
                              bool                  update_kinematics)
        {
            if (update_kinematics)
            {
                updateKinematicsCustom(model, &q, nullptr, nullptr);
            }

            for (size_t i = 1; i < model.mBodies.size(); i++)
            {
                model.Ic[i] = model.I[i];
            }

            RigidBodyInertia Itot;

            for (size_t i = model.mBodies.size() - 1; i > 0; i--)
            {
                unsigned int lambda = model.lambda[i];

                if (lambda != 0)
                {
                    model.Ic[lambda] = model.Ic[lambda] + model.Ic[i].transform_transpose_copy(model.bodyFrames[i]->getTransformFromParent());
                }
                else
                {
                    Itot = Itot + model.Ic[i].transform_transpose_copy(model.bodyFrames[i]->getTransformFromParent());
                }
            }

            com = Itot.h / Itot.m;
        }

        void calcCenterOfMass(Model               & model,
                              const Math::VectorNd& q,
                              FramePointd      & com,
                              bool                  update_kinematics)
        {
            Vector3d com_v;
            calcCenterOfMass(model,q,com_v,update_kinematics);
            com.setIncludingFrame(com_v,model.worldFrame.get());
        }

        void calcCenterOfMass(Model &model, const Math::VectorNd &q, const Math::VectorNd &qdot,
                                         double &mass, Math::FramePointd &com, Math::FrameVector *com_velocity,
                                         Math::FrameVector *angular_momentum, bool update_kinematics)
        {
            Vector3d comPosition, comVelocity, angularMomentum;

            calcCenterOfMass(model, q, qdot, mass, comPosition, &comVelocity, &angularMomentum, update_kinematics);

            ReferenceFrame *worldFrame = model.worldFrame.get();

            com.setIncludingFrame(comPosition, worldFrame);

            if (com_velocity)
            {
                com_velocity->setIncludingFrame(comVelocity, worldFrame);
            }

            if (angular_momentum)
            {
                angular_momentum->setIncludingFrame(angularMomentum, worldFrame);
            }
        }

        void calcCenterOfMass(Model &model, const Math::VectorNd &q, const Math::VectorNd &qdot,
                                         double &mass, Math::Vector3d &com, Math::Vector3d *com_velocity, Vector3d *angular_momentum,
                                         bool update_kinematics)
        {
            if (update_kinematics)
            {
                updateKinematicsCustom(model, &q, &qdot, NULL);
            }

            for (size_t i = 1; i < model.mBodies.size(); i++)
            {
                model.Ic[i] = model.I[i];
                model.hc[i] = model.Ic[i] * model.v[i];
            }

            RigidBodyInertia Itot;
            SpatialVector htot(SpatialVector::Zero(6));

            for (size_t i = model.mBodies.size() - 1; i > 0; i--)
            {
                unsigned int lambda = model.lambda[i];

                if (lambda != 0)
                {
                    model.Ic[lambda] = model.Ic[lambda] + model.Ic[i].transform_transpose_copy(model.bodyFrames[i]->getTransformFromParent());
                    model.hc[lambda] = model.hc[lambda] + model.bodyFrames[i]->getTransformFromParent().applyTranspose(model.hc[i]);
                }
                else
                {
                    Itot = Itot + model.Ic[i].transform_transpose_copy(model.bodyFrames[i]->getTransformFromParent());
                    htot = htot + model.bodyFrames[i]->getTransformFromParent().applyTranspose(model.hc[i]);
                }
            }

            mass = Itot.m;
            com = Itot.h / mass;

            if (com_velocity)
            {
                *com_velocity = Vector3d(htot[3] / mass, htot[4] / mass, htot[5] / mass);
            }

            if (angular_momentum)
            {
                htot = Xtrans(com).applyAdjoint(htot);
                angular_momentum->set(htot[0], htot[1], htot[2]);
            }
        }

        void calcCenterOfMass(Model &model, const Math::VectorNd &q, const Math::VectorNd &qdot,
                                         Math::FramePointd &com, Math::FrameVector *com_velocity, bool update_kinematics)
        {
            if (update_kinematics)
            {
                if (com_velocity)
                {
                    updateKinematicsCustom(model, &q, &qdot, nullptr);
                }
                else
                {
                    updateKinematicsCustom(model, &q, nullptr, nullptr);
                }
            }

            Vector3d p_com(0., 0., 0.), v_com(0., 0., 0.);
            double mass = 0;
            for (unsigned int i = 1; i < model.mBodies.size(); i++)
            {
                double bodyMass = model.mBodies[i].mMass;
                mass += bodyMass;
                p_com += bodyMass * model.bodyCenteredFrames[i]->getInverseTransformToRoot().r;

                if (com_velocity)
                {
                    v_com += bodyMass * (model.bodyFrames[i]->getTransformToRoot().E * (model.v[i].getLinearPart() - model.mBodies[i].mCenterOfMass.cross(model.v[i].getAngularPart())));
                }
            }

            com.setIncludingFrame(p_com / mass, model.worldFrame.get());

            if (com_velocity)
            {
                com_velocity->setIncludingFrame(v_com / mass, model.worldFrame.get());
            }
        }

        void
        calcCenterOfMassJacobian(Model &model, const Math::VectorNd &q, Math::MatrixNd &jCom, bool update_kinematics)
        {
            assert(jCom.cols() == model.qdot_size && jCom.rows() == 3);

            if (update_kinematics)
            {
                updateKinematicsCustom(model, &q, nullptr, nullptr);
            }

            Math::SpatialMatrix K;
            for (unsigned int i = 1; i < model.mBodies.size(); i++)
            {
                K.block<3,3>(3, 0) = -calcSubtreeCenterOfMassScaledByMass(model, i, q, false).toTildeForm();
                K.block<3,3>(3, 3) = Math::Matrix3dIdentity * calcSubtreeMass(model, i);
                if (model.mJoints[i].mJointType != JointTypeCustom)
                {
                    if (model.mJoints[i].mDoFCount == 1)
                    {
                        jCom.col(model.mJoints[i].q_index) = (K*model.S[i].transform_copy(model.bodyFrames[i]->getTransformToRoot())).block(3,0,3,1);
                    }
                    else if (model.mJoints[i].mDoFCount == 3)
                    {
                        jCom.block(0,model.mJoints[i].q_index,3,3) = (K*(model.bodyFrames[i]->getTransformToRoot().toMatrix()*model.multdof3_S[i])).block(3,0,3,3);
                    }
                }
                else
                {
                    jCom.block(0, model.mJoints[i].q_index, 3, model.mJoints[i].mDoFCount) = (K*(model.bodyFrames[i]->getTransformToRoot().toMatrix() * model.mCustomJoints[model.mJoints[i].custom_joint_index]->S)).block(3,0,3,model.mJoints[i].mDoFCount);
                }
            }

            jCom /= calcSubtreeMass(model, 0);
        }

        Math::FramePointd
        calcSubtreeCenterOfMassScaledByMass(Model &model, const unsigned int bodyId, const VectorNd &q, bool updateKinematics)
        {
            if (updateKinematics)
            {
                updateKinematicsCustom(model, &q, nullptr, nullptr);
            }

            std::vector<unsigned int> childBodyIds = model.mu[bodyId];

            Math::FramePointd comScaledByMass(model.worldFrame.get(), (model.bodyCenteredFrames[bodyId]->getInverseTransformToRoot().r * model.mBodies[bodyId].mMass));

            for (unsigned int j = 0; j < childBodyIds.size(); j++)
            {
                comScaledByMass += calcSubtreeCenterOfMassScaledByMass(model, childBodyIds[j], q, false);
            }

            return comScaledByMass;
        }

        double calcSubtreeMass(Model &model, const unsigned int bodyId)
        {
            std::vector<unsigned int> childBodyIds = model.mu[bodyId];

            double subtreeMass = model.mBodies[bodyId].mMass;

            for (unsigned int j = 0; j < childBodyIds.size(); j++)
            {
                subtreeMass += calcSubtreeMass(model, childBodyIds[j]);
            }

            return subtreeMass;
        }

        double calcPotentialEnergy(Model &model, const Math::VectorNd &q, bool update_kinematics)
        {
            double mass;
            Vector3d com;
            calcCenterOfMass(model, q, VectorNd::Zero(model.qdot_size), mass, com, NULL, NULL, update_kinematics);

            Vector3d g = -Vector3d(model.gravity[3], model.gravity[4], model.gravity[5]);

            return mass * com.dot(g);
        }

        double calcKineticEnergy(Model &model, const Math::VectorNd &q, const Math::VectorNd &qdot,
                                            bool update_kinematics)
        {
            if (update_kinematics)
            {
                updateKinematicsCustom(model, &q, &qdot, NULL);
            }

            double result = 0.;

            for (size_t i = 1; i < model.mBodies.size(); i++)
            {
                result += (0.5 * model.v[i].dot(model.I[i] * model.v[i]));
            }
            return result;
        }

        void calcCentroidalMomentumMatrix(Model &model, const Math::VectorNd &q, Math::MatrixNd &A,
                bool update_kinematics)
        {
            assert(A.cols() == model.qdot_size && A.rows() == 6);

            if(update_kinematics)
            {
                updateKinematicsCustom(model,&q,nullptr,nullptr);
            }

            Vector3d com;
            calcCenterOfMass(model,q,com,false);
            SpatialTransform X_com = Xtrans(com);

            // cppcheck-suppress variableScope
            unsigned int j;
            for(unsigned int i = 1; i < model.mBodies.size(); i++)
            {
                j = i;
                ReferenceFrame* bodyFrame = model.bodyFrames[i].get();
                while(j != 0)
                {
                    if(model.mJoints[j].mJointType != JointTypeCustom)
                    {
                        if(model.mJoints[j].mDoFCount == 1)
                        {
                            A.col(model.mJoints[j].q_index) += (model.I[i]*model.S[j].transform_copy(model.bodyFrames[j]->getTransformToDesiredFrame(bodyFrame))).transform_copy(bodyFrame->getTransformToRoot()).transform_copy(X_com);
                        }
                        else if(model.mJoints[j].mDoFCount == 3)
                        {
                            for(int k = 0;k<3;k++)
                            {
                                A.col(model.mJoints[j].q_index+k) += X_com.toMatrixAdjoint()*bodyFrame->getTransformToRoot().toMatrixAdjoint()*(model.I[i].toMatrix()*MotionVector(model.multdof3_S[j].col(k)).transform_copy(model.bodyFrames[j]->getTransformToDesiredFrame(bodyFrame)));
                            }
                        }
                    }
                    else if(model.mJoints[j].mJointType == JointTypeCustom)
                    {
                        unsigned int k = model.mJoints[j].custom_joint_index;

                        A.block(0, model.mJoints[j].q_index, 6, model.mCustomJoints[k]->mDoFCount) += X_com.toMatrixAdjoint()*bodyFrame->getTransformToRoot().toMatrixAdjoint()*(model.I[i].toMatrix()*model.bodyFrames[j]->getTransformToDesiredFrame(bodyFrame).toMatrix() * model.mCustomJoints[k]->S);
                    }

                    j = model.lambda[j];
                }
            }
        }

        void calcCentroidalMomentumMatrixDot(Model& model, const Math::VectorNd& q, const Math::VectorNd& qdot, Math::MatrixNd& Adot, bool update_kinematics)
        {
            assert(Adot.cols() == model.qdot_size && Adot.rows() == 6);

            if(update_kinematics)
            {
                updateKinematicsCustom(model,&q,&qdot,nullptr);
            }

            FramePointd com;
            FrameVector com_velocity;
            calcCenterOfMass(model,q,qdot,com,&com_velocity,false);
            SpatialTransform X_com = Xtrans(com.vec());
            MotionVector com_twist(0.,0.,0.,com_velocity.x(),com_velocity.y(),com_velocity.z());
            // cppcheck-suppress variableScope
            MotionVector v_tmp;

            // cppcheck-suppress variableScope
            unsigned int j;
            for(unsigned int i = 1; i < model.mBodies.size(); i++)
            {
                j = i;
                ReferenceFrame* bodyFrame = model.bodyFrames[i].get();

                v_tmp = -MotionVector(com_twist-model.v[i].transform_copy(bodyFrame->getTransformToRoot())).transform_copy(X_com);

                while(j != 0)
                {
                    if(model.mJoints[j].mJointType != JointTypeCustom)
                    {
                        if(model.mJoints[j].mDoFCount == 1)
                        {
                            Adot.col(model.mJoints[j].q_index) += v_tmp.crossf()*(model.I[i]*model.S[j].transform_copy(model.bodyFrames[j]->getTransformToDesiredFrame(bodyFrame))).transform_copy(bodyFrame->getTransformToRoot()).transform_copy(X_com) +
                                    (model.I[i]*MotionVector(model.S_o[j] + model.v[j]%model.S[j]).transform_copy(model.bodyFrames[j]->getTransformToDesiredFrame(bodyFrame))).transform_copy(bodyFrame->getTransformToRoot()).transform_copy(X_com);
                        }
                        else if(model.mJoints[j].mDoFCount == 3)
                        {
                            for(int k = 0;k<3;k++)
                            {
                                Adot.col(model.mJoints[j].q_index+k) += v_tmp.crossf()*X_com.toMatrixAdjoint()*bodyFrame->getTransformToRoot().toMatrixAdjoint()*(model.I[i].toMatrix()*MotionVector(model.multdof3_S[j].col(k)).transform_copy(model.bodyFrames[j]->getTransformToDesiredFrame(bodyFrame))) +
                                        X_com.toMatrixAdjoint()*bodyFrame->getTransformToRoot().toMatrixAdjoint()*
                                        (model.I[i].toMatrix()*MotionVector(model.multdof3_S_o[j].col(k) + model.v[j].crossm()*
                                        model.multdof3_S[j].col(k)).transform_copy(model.bodyFrames[j]->getTransformToDesiredFrame(bodyFrame)));
                            }
                        }
                    }
                    else if(model.mJoints[j].mJointType == JointTypeCustom)
                    {
                        unsigned int k = model.mJoints[j].custom_joint_index;

                        Adot.block(0, model.mJoints[j].q_index, 6, model.mCustomJoints[k]->mDoFCount) += v_tmp.crossf()*X_com.toMatrixAdjoint()*bodyFrame->getTransformToRoot().toMatrixAdjoint()*(model.I[i].toMatrix()*model.bodyFrames[j]->getTransformToDesiredFrame(bodyFrame).toMatrix() * model.mCustomJoints[k]->S) +
                                X_com.toMatrixAdjoint()*bodyFrame->getTransformToRoot().toMatrixAdjoint()*(model.I[i].toMatrix()*model.bodyFrames[j]->getTransformToDesiredFrame(bodyFrame).toMatrix() * (model.mCustomJoints[k]->S_o + model.v[j].crossm()*model.mCustomJoints[k]->S));
                    }

                    j = model.lambda[j];
                }
            }
        }
    }
}
