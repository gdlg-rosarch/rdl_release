/*
 * Original Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 *
 * RDL - Robot Dynamics Library
 * Modifications Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <limits>
#include <string.h>

#include "rdl_dynamics/Dynamics.h"

namespace RobotDynamics
{

    using namespace Math;

    RDL_DLLAPI void inverseDynamics(Model &model, const VectorNd &Q, const VectorNd &QDot, const VectorNd &QDDot,
            VectorNd &Tau, std::vector<ForceVector> *f_ext)
    {
        // Reset the velocity of the root body
        model.v[0].setZero();
        model.a[0].set(-model.gravity);

        for(unsigned int i = 1; i < model.mBodies.size(); i++)
        {
            unsigned int q_index = model.mJoints[i].q_index;
            unsigned int lambda = model.lambda[i];
            ReferenceFrame* bodyFrame = model.bodyFrames[i].get();

            jcalc(model, i, Q, QDot);

            model.v[i].set(model.v[lambda].transform_copy(bodyFrame->getTransformFromParent()) + model.v_J[i]);
            model.c[i] = model.c_J[i] + model.v[i]%model.v_J[i];

            if(model.mJoints[i].mJointType != JointTypeCustom)
            {
                if(model.mJoints[i].mDoFCount == 1)
                {
                    model.a[i].set(model.a[lambda].transform_copy(bodyFrame->getTransformFromParent()) + model.c[i] + model.S[i] * QDDot[q_index]);
                }
                else if(model.mJoints[i].mDoFCount == 3)
                {
                    model.a[i].set(bodyFrame->getTransformFromParent().apply(model.a[lambda]) + model.c[i] + model.multdof3_S[i] * Vector3d(QDDot[q_index], QDDot[q_index + 1], QDDot[q_index + 2]));
                }
            }
            else if(model.mJoints[i].mJointType == JointTypeCustom)
            {
                unsigned int k = model.mJoints[i].custom_joint_index;
                for(unsigned int z = 0; z < model.mCustomJoints[k]->mDoFCount; ++z)
                {
                    model.mCustomJoints[k]->tmp_ndof_vec[z] = QDDot[q_index + z];
                }

                model.a[i].set(bodyFrame->getTransformFromParent().apply(model.a[lambda]) + model.c[i] + model.mCustomJoints[k]->S * model.mCustomJoints[k]->tmp_ndof_vec);
            }

            if(!model.mBodies[i].mIsVirtual)
            {
                model.f_b[i] = model.I[i] * model.a[i] + model.v[i]%Momentum(model.I[i],model.v[i]);
            }
            else
            {
                model.f_b[i].setZero();
            }
        }

        for(unsigned int i = 0;i<model.fixedBodyFrames.size();i++)
        {
            model.fixedBodyFrames[i]->update();
        }

        if(f_ext != NULL)
        {
            for(unsigned int i = 1; i < model.mBodies.size(); i++)
            {
                model.f_b[i] -= (*f_ext)[i].transform_copy(model.bodyFrames[i]->getInverseTransformToRoot());
            }
        }

        for(unsigned int i = model.mBodies.size() - 1; i > 0; i--)
        {
            if(model.mJoints[i].mJointType != JointTypeCustom)
            {
                if(model.mJoints[i].mDoFCount == 1)
                {
                    Tau[model.mJoints[i].q_index] = model.S[i].dot(model.f_b[i]);
                }
                else if(model.mJoints[i].mDoFCount == 3)
                {
                    Tau.block<3, 1>(model.mJoints[i].q_index, 0) = model.multdof3_S[i].transpose() * model.f_b[i];
                }
            }
            else if(model.mJoints[i].mJointType == JointTypeCustom)
            {
                unsigned int k = model.mJoints[i].custom_joint_index;
                Tau.block(model.mJoints[i].q_index, 0, model.mCustomJoints[k]->mDoFCount, 1) = model.mCustomJoints[k]->S.transpose() * model.f_b[i];
            }

            if(model.lambda[i] != 0)
            {
                model.f_b[model.lambda[i]] = model.f_b[model.lambda[i]] + model.f_b[i].transformTranspose_copy(model.bodyFrames[i]->getTransformFromParent());
            }
        }
    }

    RDL_DLLAPI void nonlinearEffects(Model &model, const VectorNd &Q, const VectorNd &QDot, VectorNd &Tau)
    {
        model.v[0].setZero();
        model.a[0].set(-model.gravity);

        for(unsigned int i = 1; i < model.mBodies.size(); i++)
        {
            jcalc(model, i, Q, QDot);
        }

        for(unsigned int i = 0;i<model.fixedBodyFrames.size();i++)
        {
            model.fixedBodyFrames[i]->update();
        }

        for(unsigned int i = 1; i < model.mBodies.size(); i++)
        {
            if(model.lambda[i] == 0)
            {
                model.v[i].set(model.v_J[i]);
                model.a[i].set(-model.gravity.transform_copy(model.bodyFrames[i]->getTransformFromParent()));
            }
            else
            {
                model.v[i].set(model.v[model.lambda[i]].transform_copy(model.bodyFrames[i]->getTransformFromParent()) + model.v_J[i]);
                model.c[i] = model.c_J[i] + model.v[i]%model.v_J[i];
                model.a[i].set(model.a[model.lambda[i]].transform_copy(model.bodyFrames[i]->getTransformFromParent()) + model.c[i]);
            }

            if(!model.mBodies[i].mIsVirtual)
            {
                model.f_b[i].set(model.I[i] * model.a[i] + model.v[i]%(Momentum(model.I[i],model.v[i])));
            }
            else
            {
                model.f_b[i].setZero();
            }
        }

        for(unsigned int i = model.mBodies.size() - 1; i > 0; i--)
        {
            if(model.mJoints[i].mJointType != JointTypeCustom)
            {
                if(model.mJoints[i].mDoFCount == 1)
                {
                    Tau[model.mJoints[i].q_index] = model.S[i].dot(model.f_b[i]);
                }
                else if(model.mJoints[i].mDoFCount == 3)
                {
                    Tau.block<3, 1>(model.mJoints[i].q_index, 0) = model.multdof3_S[i].transpose() * model.f_b[i];
                }
            }
            else if(model.mJoints[i].mJointType == JointTypeCustom)
            {
                unsigned int k = model.mJoints[i].custom_joint_index;
                Tau.block(model.mJoints[i].q_index, 0, model.mCustomJoints[k]->mDoFCount, 1) = model.mCustomJoints[k]->S.transpose() * model.f_b[i];
            }

            if(model.lambda[i] != 0)
            {
                model.f_b[model.lambda[i]].set(model.f_b[model.lambda[i]] + model.f_b[i].transformTranspose_copy(model.bodyFrames[i]->getTransformFromParent()));
            }
        }
    }

    RDL_DLLAPI void compositeRigidBodyAlgorithm(Model &model, const VectorNd &Q, MatrixNd &H, bool update_kinematics)
    {
        assert (H.rows() == model.dof_count && H.cols() == model.dof_count);

        for(unsigned int i = 1; i < model.mBodies.size(); i++)
        {
            if(update_kinematics)
            {
                jcalc_X_lambda_S(model, i, Q);
            }
            model.Ic[i] = model.I[i];
        }

        for(unsigned int i = 0;i<model.fixedBodyFrames.size();i++)
        {
            model.fixedBodyFrames[i]->update();
        }

        for(unsigned int i = model.mBodies.size() - 1; i > 0; i--)
        {
            if(model.lambda[i] != 0)
            {
                model.Ic[model.lambda[i]] = model.Ic[model.lambda[i]] + model.Ic[i].transform_transpose_copy(model.bodyFrames[i]->getTransformFromParent());
            }

            unsigned int dof_index_i = model.mJoints[i].q_index;

            if(model.mJoints[i].mDoFCount == 1 && model.mJoints[i].mJointType != JointTypeCustom)
            {
                Momentum F(model.Ic[i],model.S[i]);
                H(dof_index_i, dof_index_i) = model.S[i].dot(F);

                unsigned int j = i;
                unsigned int dof_index_j = dof_index_i;

                while(model.lambda[j] != 0)
                {
                    F.transformTranspose(model.bodyFrames[j]->getTransformFromParent());
                    j = model.lambda[j];
                    dof_index_j = model.mJoints[j].q_index;

                    if(model.mJoints[j].mJointType != JointTypeCustom)
                    {
                        if(model.mJoints[j].mDoFCount == 1)
                        {
                            H(dof_index_i, dof_index_j) = F.dot(model.S[j]);
                            H(dof_index_j, dof_index_i) = H(dof_index_i, dof_index_j);
                        }
                        else if(model.mJoints[j].mDoFCount == 3)
                        {
                            Vector3d H_temp2 = (F.transpose() * model.multdof3_S[j]).transpose();

                            H.block<1, 3>(dof_index_i, dof_index_j) = H_temp2.transpose();
                            H.block<3, 1>(dof_index_j, dof_index_i) = H_temp2;
                        }
                    }
                    else if(model.mJoints[j].mJointType == JointTypeCustom)
                    {
                        unsigned int k = model.mJoints[j].custom_joint_index;
                        unsigned int dof = model.mCustomJoints[k]->mDoFCount;
                        VectorNd H_temp2 = (F.transpose() * model.mCustomJoints[k]->S).transpose();

                        H.block(dof_index_i, dof_index_j, 1, dof) = H_temp2.transpose();
                        H.block(dof_index_j, dof_index_i, dof, 1) = H_temp2;
                    }
                }
            }
            else if(model.mJoints[i].mDoFCount == 3 && model.mJoints[i].mJointType != JointTypeCustom)
            {
                Matrix63 F_63 = model.Ic[i].toMatrix() * model.multdof3_S[i];
                H.block<3, 3>(dof_index_i, dof_index_i) = model.multdof3_S[i].transpose() * F_63;

                unsigned int j = i;
                unsigned int dof_index_j = dof_index_i;

                while(model.lambda[j] != 0)
                {
                    F_63 = model.bodyFrames[j]->getTransformFromParent().toMatrixTranspose() * (F_63);
                    j = model.lambda[j];
                    dof_index_j = model.mJoints[j].q_index;

                    if(model.mJoints[j].mJointType != JointTypeCustom)
                    {
                        if(model.mJoints[j].mDoFCount == 1)
                        {
                            Vector3d H_temp2 = F_63.transpose() * (model.S[j]);

                            H.block<3, 1>(dof_index_i, dof_index_j) = H_temp2;
                            H.block<1, 3>(dof_index_j, dof_index_i) = H_temp2.transpose();
                        }
                        else if(model.mJoints[j].mDoFCount == 3)
                        {
                            Matrix3d H_temp2 = F_63.transpose() * (model.multdof3_S[j]);

                            H.block<3, 3>(dof_index_i, dof_index_j) = H_temp2;
                            H.block<3, 3>(dof_index_j, dof_index_i) = H_temp2.transpose();
                        }
                    }
                    else if(model.mJoints[j].mJointType == JointTypeCustom)
                    {
                        unsigned int k = model.mJoints[j].custom_joint_index;
                        unsigned int dof = model.mCustomJoints[k]->mDoFCount;

                        MatrixNd H_temp2 = F_63.transpose() * (model.mCustomJoints[k]->S);

                        H.block(dof_index_i, dof_index_j, 3, dof) = H_temp2;
                        H.block(dof_index_j, dof_index_i, dof, 3) = H_temp2.transpose();
                    }
                }
            }
            else if(model.mJoints[i].mJointType == JointTypeCustom)
            {
                unsigned int kI = model.mJoints[i].custom_joint_index;
                unsigned int dofI = model.mCustomJoints[kI]->mDoFCount;

                MatrixNd F_Nd = model.Ic[i].toMatrix() * model.mCustomJoints[kI]->S;

                H.block(dof_index_i, dof_index_i, dofI, dofI) = model.mCustomJoints[kI]->S.transpose() * F_Nd;

                unsigned int j = i;
                unsigned int dof_index_j = dof_index_i;

                while(model.lambda[j] != 0)
                {
                    F_Nd = model.bodyFrames[j]->getTransformFromParent().toMatrixTranspose() * (F_Nd);
                    j = model.lambda[j];
                    dof_index_j = model.mJoints[j].q_index;

                    if(model.mJoints[j].mJointType != JointTypeCustom)
                    {
                        if(model.mJoints[j].mDoFCount == 1)
                        {
                            MatrixNd H_temp2 = F_Nd.transpose() * (model.S[j]);
                            H.block(dof_index_i, dof_index_j, H_temp2.rows(), H_temp2.cols()) = H_temp2;
                            H.block(dof_index_j, dof_index_i, H_temp2.cols(), H_temp2.rows()) = H_temp2.transpose();
                        }
                        else if(model.mJoints[j].mDoFCount == 3)
                        {
                            MatrixNd H_temp2 = F_Nd.transpose() * (model.multdof3_S[j]);
                            H.block(dof_index_i, dof_index_j, H_temp2.rows(), H_temp2.cols()) = H_temp2;
                            H.block(dof_index_j, dof_index_i, H_temp2.cols(), H_temp2.rows()) = H_temp2.transpose();
                        }
                    }
                    else if(model.mJoints[j].mJointType == JointTypeCustom)
                    {
                        unsigned int k = model.mJoints[j].custom_joint_index;
                        unsigned int dof = model.mCustomJoints[k]->mDoFCount;

                        MatrixNd H_temp2 = F_Nd.transpose() * (model.mCustomJoints[k]->S);

                        H.block(dof_index_i, dof_index_j, 3, dof) = H_temp2;
                        H.block(dof_index_j, dof_index_i, dof, 3) = H_temp2.transpose();
                    }
                }
            }
        }
    }

    RDL_DLLAPI void forwardDynamics(Model &model, const VectorNd &Q, const VectorNd &QDot, const VectorNd &Tau,
            VectorNd &QDDot, std::vector<ForceVector> *f_ext)
    {
        unsigned int i = 0;

        // Reset the velocity of the root body
        model.v[0].setZero();

        for(i = 1; i < model.mBodies.size(); i++)
        {
            unsigned int lambda = model.lambda[i];
            ReferenceFrame *bodyFrame = model.bodyFrames[i].get();

            jcalc(model, i, Q, QDot);

            model.v[i].set(model.v[lambda].transform_copy(bodyFrame->getTransformFromParent()) + model.v_J[i]);

            model.c[i] = model.c_J[i] + SpatialVector(model.v[i] % model.v_J[i]);
            model.I[i].setSpatialMatrix(model.IA[i]);

            model.pA[i] = model.v[i].cross(model.I[i] * model.v[i]);

            if(f_ext != NULL && (*f_ext)[i] != SpatialVectorZero)
            {
                model.pA[i] -= (*f_ext)[i].transform_copy(bodyFrame->getInverseTransformToRoot());
            }
        }

        for(unsigned int i = 0; i<model.fixedBodyFrames.size();i++)
        {
            model.fixedBodyFrames[i]->update();
        }

        for(i = model.mBodies.size() - 1; i > 0; i--)
        {
            unsigned int q_index = model.mJoints[i].q_index;
            ReferenceFrame *bodyFrame = model.bodyFrames[i].get();

            if(model.mJoints[i].mDoFCount == 1 && model.mJoints[i].mJointType != JointTypeCustom)
            {
                model.U[i] = model.IA[i] * model.S[i];
                model.d[i] = model.S[i].dot(model.U[i]);
                model.u[i] = Tau[q_index] - model.S[i].dot(model.pA[i]);

                unsigned int lambda = model.lambda[i];
                if(lambda != 0)
                {
                    SpatialMatrix Ia = model.IA[i] - model.U[i] * (model.U[i] / model.d[i]).transpose();
                    ForceVector pa = model.pA[i] + Ia * model.c[i] + model.U[i] * model.u[i] / model.d[i];

                    model.IA[lambda].noalias() += bodyFrame->getTransformFromParent().toMatrixTranspose() * Ia * bodyFrame->getTransformFromParent().toMatrix();
                    model.pA[lambda].noalias() += pa.transformTranspose_copy(bodyFrame->getTransformFromParent());
                }
            }
            else if(model.mJoints[i].mDoFCount == 3 && model.mJoints[i].mJointType != JointTypeCustom)
            {
                model.multdof3_U[i] = model.IA[i]*model.multdof3_S[i];

                model.multdof3_Dinv[i] = (model.multdof3_S[i].transpose() * model.multdof3_U[i]).inverse().eval();

                Vector3d tau_temp(Tau[q_index], Tau[q_index + 1], Tau[q_index + 2]);
                model.multdof3_u[i] = tau_temp - model.multdof3_S[i].transpose() * model.pA[i];

                unsigned int lambda = model.lambda[i];
                if(lambda != 0)
                {
                    SpatialMatrix Ia = model.IA[i]- model.multdof3_U[i] * model.multdof3_Dinv[i] * model.multdof3_U[i].transpose();
                    ForceVector pa = model.pA[i] + Ia * model.c[i] + model.multdof3_U[i] * model.multdof3_Dinv[i] * model.multdof3_u[i];

                    model.IA[lambda].noalias() += bodyFrame->getTransformFromParent().toMatrixTranspose() * Ia * bodyFrame->getTransformFromParent().toMatrix();
                    model.pA[lambda].noalias() += pa.transformTranspose_copy(bodyFrame->getTransformFromParent());
                }
            }
            else if(model.mJoints[i].mJointType == JointTypeCustom)
            {
                unsigned int kI = model.mJoints[i].custom_joint_index;
                model.mCustomJoints[kI]->U = model.IA[i] * model.mCustomJoints[kI]->S;

                model.mCustomJoints[kI]->Dinv = (model.mCustomJoints[kI]->S.transpose() * model.mCustomJoints[kI]->U).inverse().eval();

                for(unsigned int z = 0; z < model.mCustomJoints[kI]->mDoFCount; ++z)
                {
                    model.mCustomJoints[kI]->tmp_ndof_vec(z) = Tau[q_index + z];
                }

                model.mCustomJoints[kI]->u = model.mCustomJoints[kI]->tmp_ndof_vec - model.mCustomJoints[kI]->S.transpose() * model.pA[i];

                unsigned int lambda = model.lambda[i];
                if(lambda != 0)
                {
                    SpatialMatrix Ia = model.IA[i] - (model.mCustomJoints[kI]->U * model.mCustomJoints[kI]->Dinv * model.mCustomJoints[kI]->U.transpose());
                    ForceVector pa = model.pA[i] + Ia * model.c[i] + (model.mCustomJoints[kI]->U * model.mCustomJoints[kI]->Dinv * model.mCustomJoints[kI]->u);

                    model.IA[lambda].noalias() += bodyFrame->getTransformFromParent().toMatrixTranspose() * Ia * bodyFrame->getTransformFromParent().toMatrix();
                    model.pA[lambda].noalias() += pa.transformTranspose_copy(bodyFrame->getTransformFromParent());
                }
            }
        }

         model.a[0].set(model.gravity * -1.);

         for(i = 1; i < model.mBodies.size(); i++)
         {
             unsigned int q_index = model.mJoints[i].q_index;

             model.a[i].set(model.a[model.lambda[i]].transform_copy(model.bodyFrames[i]->getTransformFromParent()) + model.c[i]);

             if(model.mJoints[i].mDoFCount == 1 && model.mJoints[i].mJointType != JointTypeCustom)
             {
                 QDDot[q_index] = (1. / model.d[i]) * (model.u[i] - model.U[i].dot(model.a[i]));
                 model.a[i].set(model.a[i] + model.S[i] * QDDot[q_index]);
             }
             else if(model.mJoints[i].mDoFCount == 3 && model.mJoints[i].mJointType != JointTypeCustom)
             {
                 Vector3d qdd_temp = model.multdof3_Dinv[i] * (model.multdof3_u[i] - model.multdof3_U[i].transpose() * model.a[i]);
                 QDDot[q_index] = qdd_temp[0];
                 QDDot[q_index + 1] = qdd_temp[1];
                 QDDot[q_index + 2] = qdd_temp[2];
                 model.a[i].set(model.a[i] + model.multdof3_S[i] * qdd_temp);
             }
             else if(model.mJoints[i].mJointType == JointTypeCustom)
             {
                 unsigned int kI = model.mJoints[i].custom_joint_index;
                 unsigned int dofI = model.mCustomJoints[kI]->mDoFCount;

                 VectorNd qdd_temp = model.mCustomJoints[kI]->Dinv * (model.mCustomJoints[kI]->u - model.mCustomJoints[kI]->U.transpose() * model.a[i]);

                 for(unsigned int z = 0; z < dofI; ++z)
                 {
                     QDDot[q_index + z] = qdd_temp[z];
                 }

                 model.a[i].set(model.a[i] + model.mCustomJoints[kI]->S * qdd_temp);
             }
         }
    }

    RDL_DLLAPI void forwardDynamicsLagrangian(Model &model, const VectorNd &Q, const VectorNd &QDot,
            const VectorNd &Tau, VectorNd &QDDot, Math::LinearSolver linear_solver, std::vector<ForceVector> *f_ext,
            Math::MatrixNd *H, Math::VectorNd *C)
    {
        bool free_H = false;
        bool free_C = false;

        if(!H)
        {
            model.tmp_ndof_mat.setZero();
            H = &model.tmp_ndof_mat;
            free_H = true;
        }

        if(!C)
        {
            model.tmp_ndof_vec.setZero();
            C = &model.tmp_ndof_vec;
            free_C = true;
        }

        // we set QDDot to zero to compute C properly with the InverseDynamics
        // method.
        QDDot.setZero();

        inverseDynamics(model, Q, QDot, QDDot, (*C), f_ext);
        compositeRigidBodyAlgorithm(model, Q, *H, false);

        bool solve_successful = linSolveGaussElimPivot (*H, *C * -1. + Tau, QDDot);
        assert (solve_successful);

        if(free_H)
        {
            // cppcheck-suppress uselessAssignmentPtrArg
            H = nullptr;
        }

        if(free_C)
        {
            // cppcheck-suppress uselessAssignmentPtrArg
            C = nullptr;
        }
    }

    RDL_DLLAPI void calcMInvTimesTau(Model &model, const VectorNd &Q, const VectorNd &Tau, VectorNd &QDDot,
            bool update_kinematics)
    {
        model.v[0].setZero();
        model.a[0].setZero();

        if(update_kinematics)
        {
            for(unsigned int i = 1; i < model.mBodies.size(); i++)
            {
                jcalc_X_lambda_S(model, model.mJointUpdateOrder[i], Q);

                model.v_J[i].setZero();
                model.v[i].setZero();
                model.c_J[i].setZero();
                model.pA[i].setZero();
                model.I[i].setSpatialMatrix(model.IA[i]);
            }

            for(unsigned int i = 0; i<model.fixedBodyFrames.size();i++)
            {
                model.fixedBodyFrames[i]->update();
            }
        }

        for(unsigned int i = 1; i < model.mBodies.size(); i++)
        {
            model.pA[i].setZero();
        }

        if(update_kinematics)
        {
            // Compute Articulate Body Inertias
            for(unsigned int i = model.mBodies.size() - 1; i > 0; i--)
            {
                if(model.mJoints[i].mDoFCount == 1 && model.mJoints[i].mJointType != JointTypeCustom)
                {
                    model.U[i] = model.IA[i] * model.S[i];
                    model.d[i] = model.S[i].dot(model.U[i]);
                    unsigned int lambda = model.lambda[i];

                    if(lambda != 0)
                    {
                        SpatialMatrix Ia = model.IA[i] - model.U[i] * (model.U[i] / model.d[i]).transpose();
                        model.IA[lambda].noalias() += model.bodyFrames[i]->getTransformFromParent().toMatrixTranspose() * Ia * model.bodyFrames[i]->getTransformFromParent().toMatrix();
                    }
                }
                else if(model.mJoints[i].mDoFCount == 3 && model.mJoints[i].mJointType != JointTypeCustom)
                {

                    model.multdof3_U[i] = model.IA[i] * model.multdof3_S[i];
                    model.multdof3_Dinv[i] = (model.multdof3_S[i].transpose() * model.multdof3_U[i]).inverse().eval();

                    unsigned int lambda = model.lambda[i];

                    if(lambda != 0)
                    {
                        SpatialMatrix Ia = model.IA[i] - (model.multdof3_U[i] * model.multdof3_Dinv[i] * model.multdof3_U[i].transpose());
                        model.IA[lambda].noalias() += model.bodyFrames[i]->getTransformFromParent().toMatrixTranspose() * Ia * model.bodyFrames[i]->getTransformFromParent().toMatrix();
                    }
                }
                else if(model.mJoints[i].mJointType == JointTypeCustom)
                {
                    unsigned int kI = model.mJoints[i].custom_joint_index;
                    model.mCustomJoints[kI]->U = model.IA[i] * model.mCustomJoints[kI]->S;
                    model.mCustomJoints[kI]->Dinv = (model.mCustomJoints[kI]->S.transpose() * model.mCustomJoints[kI]->U).inverse().eval();

                    unsigned int lambda = model.lambda[i];

                    if(lambda != 0)
                    {
                        SpatialMatrix Ia = model.IA[i] - (model.mCustomJoints[kI]->U * model.mCustomJoints[kI]->Dinv * model.mCustomJoints[kI]->U.transpose());
                        model.IA[lambda].noalias() += model.bodyFrames[i]->getTransformFromParent().toMatrixTranspose() * Ia * model.bodyFrames[i]->getTransformFromParent().toMatrix();
                    }
                }
            }
        }

        // compute articulated bias forces
        for(unsigned int i = model.mBodies.size() - 1; i > 0; i--)
        {
            unsigned int q_index = model.mJoints[i].q_index;

            if(model.mJoints[i].mDoFCount == 1 && model.mJoints[i].mJointType != JointTypeCustom)
            {

                model.u[i] = Tau[q_index] - model.S[i].dot(model.pA[i]);
                unsigned int lambda = model.lambda[i];
                if(lambda != 0)
                {
                    SpatialVector pa = model.pA[i] + model.U[i] * model.u[i] / model.d[i];
                    model.pA[lambda].noalias() += model.bodyFrames[i]->getTransformFromParent().applyTranspose(pa);
                }
            }
            else if(model.mJoints[i].mDoFCount == 3 && model.mJoints[i].mJointType != JointTypeCustom)
            {
                Vector3d tau_temp(Tau[q_index], Tau[q_index + 1], Tau[q_index + 2]);
                model.multdof3_u[i] = tau_temp - model.multdof3_S[i].transpose() * model.pA[i];
                unsigned int lambda = model.lambda[i];

                if(lambda != 0)
                {
                    SpatialVector pa = model.pA[i] + model.multdof3_U[i] * model.multdof3_Dinv[i] * model.multdof3_u[i];
                    model.pA[lambda].noalias() += model.bodyFrames[i]->getTransformFromParent().applyTranspose(pa);
                }
            }
            else if(model.mJoints[i].mJointType == JointTypeCustom)
            {
                unsigned int kI = model.mJoints[i].custom_joint_index;
                unsigned int dofI = model.mCustomJoints[kI]->mDoFCount;
                VectorNd tau_temp = model.mCustomJoints[kI]->ndof0_vec;

                for(unsigned int z = 0; z < dofI; ++z)
                {
                    tau_temp(z) = Tau[q_index + z];
                }
                model.mCustomJoints[kI]->u = tau_temp - (model.mCustomJoints[kI]->S.transpose() * model.pA[i]);

                unsigned int lambda = model.lambda[i];

                if(lambda != 0)
                {
                    SpatialVector pa = model.pA[i] + (model.mCustomJoints[kI]->U * model.mCustomJoints[kI]->Dinv * model.mCustomJoints[kI]->u);
                    model.pA[lambda].noalias() += model.bodyFrames[i]->getTransformFromParent().applyTranspose(pa);
                }
            }
        }

        for(unsigned int i = 1; i < model.mBodies.size(); i++)
        {
            unsigned int q_index = model.mJoints[i].q_index;
            unsigned int lambda = model.lambda[i];

            model.a[i].set(model.bodyFrames[i]->getTransformFromParent().apply(model.a[lambda]) + model.c[i]);

            if(model.mJoints[i].mDoFCount == 1 && model.mJoints[i].mJointType != JointTypeCustom)
            {
                QDDot[q_index] = (1. / model.d[i]) * (model.u[i] - model.U[i].dot(model.a[i]));
                model.a[i].set(model.a[i] + model.S[i] * QDDot[q_index]);
            }
            else if(model.mJoints[i].mDoFCount == 3 && model.mJoints[i].mJointType != JointTypeCustom)
            {
                Vector3d qdd_temp = model.multdof3_Dinv[i] * (model.multdof3_u[i] - model.multdof3_U[i].transpose() * model.a[i]);

                QDDot[q_index] = qdd_temp[0];
                QDDot[q_index + 1] = qdd_temp[1];
                QDDot[q_index + 2] = qdd_temp[2];
                model.a[i].set(model.a[i] + model.multdof3_S[i] * qdd_temp);
            }
            else if(model.mJoints[i].mJointType == JointTypeCustom)
            {
                unsigned int kI = model.mJoints[i].custom_joint_index;
                unsigned int dofI = model.mCustomJoints[kI]->mDoFCount;

                VectorNd qdd_temp = model.mCustomJoints[kI]->Dinv * (model.mCustomJoints[kI]->u - model.mCustomJoints[kI]->U.transpose() * model.a[i]);

                for(unsigned int z = 0; z < dofI; ++z)
                {
                    QDDot[q_index + z] = qdd_temp[z];
                }

                model.a[i].set(model.a[i] + model.mCustomJoints[kI]->S * qdd_temp);
            }
        }
    }

} /* namespace RobotDynamics */
