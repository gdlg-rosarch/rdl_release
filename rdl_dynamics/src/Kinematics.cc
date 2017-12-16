/*
 * Original Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 *
 * RDL - Robot Dynamics Library
 * Modifications Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <iostream>
#include <limits>
#include <cstring>
#include <assert.h>

#include "rdl_dynamics/rdl_mathutils.h"
#include "rdl_dynamics/Kinematics.h"

namespace RobotDynamics
{
    using namespace Math;

    RDL_DLLAPI void updateKinematics(Model &model, const VectorNd &Q, const VectorNd &QDot, const VectorNd &QDDot)
    {
        model.a[0].setZero();

        for(unsigned int i = 1; i < model.mBodies.size(); i++)
        {
            unsigned int q_index = model.mJoints[i].q_index;
            ReferenceFrame *bodyFrame = model.bodyFrames[i].get();

            unsigned int lambda = model.lambda[i];

            jcalc(model, i, Q, QDot);

            if(lambda != 0)
            {
                model.v[i].set(model.v[lambda].transform_copy(bodyFrame->getTransformFromParent()) + model.v_J[i]);
            }
            else
            {
                model.v[i] = model.v_J[i];
            }

            model.c[i] = model.c_J[i] + model.v[i] % model.v_J[i];
            model.a[i].set(model.a[lambda].transform_copy(bodyFrame->getTransformFromParent())+ model.c[i]);

            if(model.mJoints[i].mJointType != JointTypeCustom)
            {
                if(model.mJoints[i].mDoFCount == 1)
                {
                    model.a[i].set(model.a[i] + model.S[i] * QDDot[q_index]);
                }
                else if(model.mJoints[i].mDoFCount == 3)
                {
                    model.a[i].set(model.a[i] + model.multdof3_S[i] * Vector3d(QDDot[q_index], QDDot[q_index + 1], QDDot[q_index + 2]));
                }
            }
            else
            {
                unsigned int custom_index = model.mJoints[i].custom_joint_index;
                const CustomJoint *custom_joint = model.mCustomJoints[custom_index];
                unsigned int joint_dof_count = custom_joint->mDoFCount;

                model.a[i].set(model.a[i] + (model.mCustomJoints[custom_index]->S * QDDot.block(q_index, 0, joint_dof_count, 1)));
            }
        }

        for(auto const &fixedBodyFrame : model.fixedBodyFrames)
        {
            fixedBodyFrame->update();
        }
    }

    RDL_DLLAPI void updateKinematicsCustom(Model &model, const VectorNd *Q, const VectorNd *QDot,
            const VectorNd *QDDot)
    {
        unsigned int i;

        if(Q)
        {
            for(i = 1; i < model.mBodies.size(); i++)
            {
                if(!QDot)
                {
                    jcalc(model, i, (*Q), model.q0_vec);
                }
                else
                {
                    jcalc(model, i, (*Q), (*QDot));
                }
            }

            for(auto const &fixedBodyFrame : model.fixedBodyFrames)
            {
                fixedBodyFrame->update();
            }
        }

        if(QDot)
        {
            for(i = 1; i < model.mBodies.size(); i++)
            {
                ReferenceFrame *bodyFrame = model.bodyFrames[i].get();
                unsigned int lambda = model.lambda[i];

                if(lambda != 0)
                {
                    model.v[i].set(model.v[lambda].transform_copy(bodyFrame->getTransformFromParent()) + model.v_J[i]);
                    model.c[i] = model.c_J[i] + model.v[i] % model.v_J[i];
                }
                else
                {
                    model.v[i].set(model.v_J[i]);
                    model.c[i] = model.c_J[i] + model.v[i] % model.v_J[i];
                }
            }
        }

        if(QDDot)
        {
            for(i = 1; i < model.mBodies.size(); i++)
            {
                unsigned int q_index = model.mJoints[i].q_index;
                ReferenceFrame *bodyFrame = model.bodyFrames[i].get();
                unsigned int lambda = model.lambda[i];

                if(lambda != 0)
                {
                    model.a[i].set(model.a[lambda].transform_copy(bodyFrame->getTransformFromParent()) + model.c[i]);
                }
                else
                {
                    model.a[i].set(model.c[i]);
                }

                if(model.mJoints[i].mJointType != JointTypeCustom)
                {
                    if(model.mJoints[i].mDoFCount == 1)
                    {
                        model.a[i].set(model.a[i] + model.S[i] * (*QDDot)[q_index]);
                    }
                    else if(model.mJoints[i].mDoFCount == 3)
                    {
                        Vector3d omegadot_temp((*QDDot)[q_index], (*QDDot)[q_index + 1], (*QDDot)[q_index + 2]);
                        model.a[i].set(model.a[i] + model.multdof3_S[i] * omegadot_temp);
                    }
                }
                else
                {
                    unsigned int custom_index = model.mJoints[i].custom_joint_index;
                    const CustomJoint *custom_joint = model.mCustomJoints[custom_index];
                    unsigned int joint_dof_count = custom_joint->mDoFCount;

                    model.a[i].set(model.a[i] + (model.mCustomJoints[custom_index]->S * (QDDot->block(q_index, 0, joint_dof_count, 1))));
                }
            }
        }
    }

    void calcRelativeBodySpatialJacobian(Model & model, const Math::VectorNd& Q, Math::MatrixNd &G, ReferenceFrame* baseFrame, ReferenceFrame* relativeFrame, ReferenceFrame* expressedInFrame, bool update_kinematics)
    {
        ReferenceFrame* expressedInFrameRef = expressedInFrame==nullptr ? baseFrame : expressedInFrame;

        if(update_kinematics)
        {
            updateKinematicsCustom(model, &Q, NULL, NULL);
        }

        unsigned int common_parent_id = model.getCommonMovableParentId(baseFrame->getMovableBodyId(),relativeFrame->getMovableBodyId());
        unsigned int j = baseFrame->getMovableBodyId();

        while(j > common_parent_id)
        {
            if(model.mJoints[j].mJointType != JointTypeCustom)
            {
                if(model.mJoints[j].mDoFCount == 1)
                {
                    G.col(model.mJoints[j].q_index) = model.S[j].transform_copy(model.bodyFrames[j]->getTransformToDesiredFrame(expressedInFrameRef));
                }
                else if(model.mJoints[j].mDoFCount == 3)
                {
                    for(int k = 0;k<3;k++)
                    {
                        G.col(model.mJoints[j].q_index+k) = MotionVector(model.multdof3_S[j].col(k)).transform_copy(model.bodyFrames[j]->getTransformToDesiredFrame(expressedInFrameRef));
                    }
                }
            }
            else
            {
                unsigned int k = model.mJoints[j].custom_joint_index;

                G.block(0, model.mJoints[j].q_index, 6, model.mCustomJoints[k]->mDoFCount) = model.bodyFrames[j]->getTransformToDesiredFrame(expressedInFrameRef).toMatrix() * model.mCustomJoints[k]->S;
            }

            j = model.lambda[j];
        }

        j = relativeFrame->getMovableBodyId();

        while(j > common_parent_id)
        {
            if(model.mJoints[j].mJointType != JointTypeCustom)
            {
                if(model.mJoints[j].mDoFCount == 1)
                {
                    G.col(model.mJoints[j].q_index) -= model.S[j].transform_copy(model.bodyFrames[j]->getTransformToDesiredFrame(expressedInFrameRef));
                }
                else if(model.mJoints[j].mDoFCount == 3)
                {
                    for(int k = 0;k<3;k++)
                    {
                        G.col(model.mJoints[j].q_index+k) -= MotionVector(model.multdof3_S[j].col(k)).transform_copy(model.bodyFrames[j]->getTransformToDesiredFrame(expressedInFrameRef));
                    }
                }
            }
            else
            {
                unsigned int k = model.mJoints[j].custom_joint_index;

                G.block(0, model.mJoints[j].q_index, 6, model.mCustomJoints[k]->mDoFCount) -= model.bodyFrames[j]->getTransformToDesiredFrame(expressedInFrameRef).toMatrix() * model.mCustomJoints[k]->S;
            }

            j = model.lambda[j];
        }
    }

    void calcRelativeBodySpatialJacobianDot(Model & model, const Math::VectorNd& Q, const Math::VectorNd& QDot, Math::MatrixNd &G, ReferenceFrame* baseFrame, ReferenceFrame* relativeFrame, ReferenceFrame* expressedInFrame, bool update_kinematics)
    {
        ReferenceFrame* expressedInFrameRef = expressedInFrame==nullptr ? baseFrame : expressedInFrame;

        if(update_kinematics)
        {
            updateKinematicsCustom(model, &Q, &QDot, NULL);
        }

        unsigned int common_parent_id = model.getCommonMovableParentId(baseFrame->getMovableBodyId(),relativeFrame->getMovableBodyId());
        unsigned int j = baseFrame->getMovableBodyId();

        MotionVector v_cur_wrt_new_expr_in_new = calcSpatialVelocity(model,Q,QDot,baseFrame,expressedInFrameRef,expressedInFrameRef,false);
        while(j != 0)
        {
            SpatialTransform X_j_to_expressedInFrame = model.bodyFrames[j]->getTransformToDesiredFrame(expressedInFrameRef);
            if(model.mJoints[j].mJointType != JointTypeCustom)
            {
                if(model.mJoints[j].mDoFCount == 1)
                {
                    G.col(model.mJoints[j].q_index) = v_cur_wrt_new_expr_in_new%model.S[j].transform_copy(X_j_to_expressedInFrame);
                    // The reason this check is here is bc it's a waste to add this bit on to then subtract it off below, which is
                    // what would be happening. We could add this bit in always, but it would result in a bunch of redundant operations.
                    // Same goes with the loop below, you will see the same check being done.
                    if(j>common_parent_id)
                    {
                        G.col(model.mJoints[j].q_index) += MotionVector(model.S_o[j] + model.v[j]%model.S[j]).transform_copy(X_j_to_expressedInFrame);
                    }
                }
                else if(model.mJoints[j].mDoFCount == 3)
                {
                    for(int k = 0;k<3;k++)
                    {
                        G.col(model.mJoints[j].q_index+k) = v_cur_wrt_new_expr_in_new.crossm()*MotionVector(model.multdof3_S[j].col(k)).transform_copy(X_j_to_expressedInFrame);

                        if(j>common_parent_id)
                        {
                            G.col(model.mJoints[j].q_index+k) += MotionVector(model.multdof3_S_o[j].col(k) + model.v[j].crossm()*model.multdof3_S[j].col(k)).transform_copy(X_j_to_expressedInFrame);
                        }
                    }
                }
            }
            else
            {
                unsigned int k = model.mJoints[j].custom_joint_index;

                G.block(0, model.mJoints[j].q_index, 6, model.mCustomJoints[k]->mDoFCount) = v_cur_wrt_new_expr_in_new.crossm()*X_j_to_expressedInFrame.toMatrix()*model.mCustomJoints[k]->S;

                if(j>common_parent_id)
                {
                    G.block(0, model.mJoints[j].q_index, 6, model.mCustomJoints[k]->mDoFCount) += X_j_to_expressedInFrame.toMatrix() * (model.mCustomJoints[k]->S_o + model.v[j].crossm()*model.mCustomJoints[k]->S);
                }
            }

            j = model.lambda[j];
        }

        j = relativeFrame->getMovableBodyId();
        v_cur_wrt_new_expr_in_new = calcSpatialVelocity(model,Q,QDot,relativeFrame,expressedInFrameRef,expressedInFrameRef,false);
        while(j != 0)
        {
            SpatialTransform X_j_to_expressedInFrame = model.bodyFrames[j]->getTransformToDesiredFrame(expressedInFrameRef);
            if(model.mJoints[j].mJointType != JointTypeCustom)
            {
                if(model.mJoints[j].mDoFCount == 1)
                {
                    G.col(model.mJoints[j].q_index) -= v_cur_wrt_new_expr_in_new%model.S[j].transform_copy(X_j_to_expressedInFrame);
                    if(j>common_parent_id)
                    {
                        G.col(model.mJoints[j].q_index) -= MotionVector(model.S_o[j] + model.v[j]%model.S[j]).transform_copy(X_j_to_expressedInFrame);
                    }
                }
                else if(model.mJoints[j].mDoFCount == 3)
                {
                    for(int k = 0;k<3;k++)
                    {
                        G.col(model.mJoints[j].q_index+k) -= v_cur_wrt_new_expr_in_new.crossm()*MotionVector(model.multdof3_S[j].col(k)).transform_copy(X_j_to_expressedInFrame);
                        if(j>common_parent_id)
                        {
                            G.col(model.mJoints[j].q_index+k) -= MotionVector(model.multdof3_S_o[j].col(k) + model.v[j].crossm()*model.multdof3_S[j].col(k)).transform_copy(X_j_to_expressedInFrame);
                        }
                    }
                }
            }
            else
            {
                unsigned int k = model.mJoints[j].custom_joint_index;

                G.block(0, model.mJoints[j].q_index, 6, model.mCustomJoints[k]->mDoFCount) -= v_cur_wrt_new_expr_in_new.crossm()*X_j_to_expressedInFrame.toMatrix()*model.mCustomJoints[k]->S;
                if(j>common_parent_id)
                {
                    G.block(0, model.mJoints[j].q_index, 6, model.mCustomJoints[k]->mDoFCount) -= X_j_to_expressedInFrame.toMatrix() * (model.mCustomJoints[k]->S_o + model.v[j].crossm()*model.mCustomJoints[k]->S);
                }
            }

            j = model.lambda[j];
        }
    }

    RDL_DLLAPI void calcPointJacobian(Model &model, const VectorNd &Q, unsigned int body_id,
            const Vector3d &point_position, MatrixNd &G, bool update_kinematics)
    {
        // update the Kinematics if necessary
        if(update_kinematics)
        {
            updateKinematicsCustom(model, &Q, NULL, NULL);
        }

        unsigned int reference_body_id = body_id;

        if(model.IsFixedBodyId(body_id))
        {
            unsigned int fbody_id = body_id - model.fixed_body_discriminator;
            reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;

            model.tmpFramePoint.setIncludingFrame(point_position,model.fixedBodyFrames[fbody_id].get());
        }
        else
        {
            model.tmpFramePoint.setIncludingFrame(point_position,model.bodyFrames[body_id].get());
        }

        model.tmpFramePoint.changeFrame(model.worldFrame.get());
        SpatialTransform point_trans = SpatialTransform(Matrix3d::Identity(3, 3), model.tmpFramePoint.vec());

        assert (G.rows() == 3 && G.cols() == model.qdot_size);

        unsigned int j = reference_body_id;

        // e[j] is set to 1 if joint j contributes to the jacobian that we are
        // computing. For all other joints the column will be zero.
        while(j != 0)
        {
            if(model.mJoints[j].mJointType != JointTypeCustom)
            {
                if(model.mJoints[j].mDoFCount == 1)
                {
                    G.col(model.mJoints[j].q_index) = model.S[j].transform_copy(point_trans*model.bodyFrames[j]->getTransformToRoot()).getLinearPart();
                }
                else if(model.mJoints[j].mDoFCount == 3)
                {
                    G.block(0, model.mJoints[j].q_index, 3, 3) = ((point_trans * model.bodyFrames[j]->getTransformToRoot()).toMatrix() * model.multdof3_S[j]).block(3, 0, 3, 3);
                }
            }
            else
            {
                unsigned int k = model.mJoints[j].custom_joint_index;

                G.block(0, model.mJoints[j].q_index, 3, model.mCustomJoints[k]->mDoFCount) = ((point_trans * model.bodyFrames[j]->getTransformToRoot()).toMatrix() * model.mCustomJoints[k]->S).block(3, 0, 3, model.mCustomJoints[k]->mDoFCount);
            }

            j = model.lambda[j];
        }
    }

    RDL_DLLAPI void calcPointJacobianDot(Model &model, const VectorNd &Q, const VectorNd &QDot, unsigned int body_id,
            const Vector3d &point_position, MatrixNd &G, bool update_kinematics)
    {
        assert (G.rows() == 3 && G.cols() == model.qdot_size);

        MatrixNd GDot6D = MatrixNd::Constant(6,model.qdot_size,0.);
        calcPointJacobianDot6D(model,Q,QDot,body_id,point_position,GDot6D,update_kinematics);

        G = GDot6D.block(3,0,3,model.qdot_size);
    }

    RDL_DLLAPI void calcPointJacobian6D(Model &model, const VectorNd &Q, unsigned int body_id,
            const Vector3d &point_position, MatrixNd &G, bool update_kinematics)
    {
        // update the Kinematics if necessary
        if(update_kinematics)
        {
            updateKinematicsCustom(model, &Q, nullptr, nullptr);
        }

        Math::FramePointd p;

        unsigned int reference_body_id = body_id;

        if(model.IsFixedBodyId(body_id))
        {
            unsigned int fbody_id = body_id - model.fixed_body_discriminator;
            reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;

            p.setIncludingFrame(point_position,model.fixedBodyFrames[fbody_id].get());
        }
        else
        {
            p.setIncludingFrame(point_position,model.bodyFrames[body_id].get());
        }

        p.changeFrame(ReferenceFrame::getWorldFrame().get());
        SpatialTransform point_trans(Matrix3dIdentity,p.vec());

        assert (G.rows() == 6 && G.cols() == model.qdot_size);

        unsigned int j = reference_body_id;

        while(j != 0)
        {
            if(model.mJoints[j].mJointType != JointTypeCustom)
            {
                if(model.mJoints[j].mDoFCount == 1)
                {
                    G.col(model.mJoints[j].q_index) = model.S[j].transform_copy(model.bodyFrames[j]->getTransformToRoot()).transform_copy(point_trans);
                }
                else if(model.mJoints[j].mDoFCount == 3)
                {
                    ReferenceFrame* frame = model.bodyFrames[j].get();
                    Matrix63 m = model.multdof3_S[j];
                    for(int k = 0;k<3;k++)
                    {
                        model.tmpMotionVector = m.col(k);
                        model.tmpMotionVector.transform(frame->getTransformToRoot());
                        G.col(model.mJoints[j].q_index+k) = model.tmpMotionVector.transform_copy(point_trans);
                    }
                }
            }
            else
            {
                unsigned int k = model.mJoints[j].custom_joint_index;

                G.block(0, model.mJoints[j].q_index, 6, model.mCustomJoints[k]->mDoFCount) = ((point_trans * model.bodyFrames[j]->getTransformToRoot()).toMatrix() * model.mCustomJoints[k]->S).block(0, 0, 6, model.mCustomJoints[k]->mDoFCount);
            }

            j = model.lambda[j];
        }
    }

    RDL_DLLAPI void calcPointJacobianDot6D(Model &model, const VectorNd &Q, const VectorNd &QDot, unsigned int body_id,
            const Vector3d &point_position, MatrixNd &G, bool update_kinematics)
    {
        // update the Kinematics if necessary
        if(update_kinematics)
        {
            updateKinematicsCustom(model, &Q, &QDot, nullptr);
        }

        Math::FramePointd p;

        unsigned int reference_body_id = body_id;

        if(model.IsFixedBodyId(body_id))
        {
            unsigned int fbody_id = body_id - model.fixed_body_discriminator;
            reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;

            p.setIncludingFrame(point_position,model.fixedBodyFrames[fbody_id].get());
        }
        else
        {
            p.setIncludingFrame(point_position,model.bodyFrames[body_id].get());
        }

        p.changeFrame(ReferenceFrame::getWorldFrame().get());
        SpatialTransform point_trans(Matrix3dIdentity,p.vec());
        MotionVector v_p = -calcPointVelocity6D(model,Q,QDot,body_id,point_position,false);
        MotionVector vw_x_S;

        assert (G.rows() == 6 && G.cols() == model.qdot_size);

        unsigned int j = reference_body_id;

        while(j != 0)
        {
            if(model.mJoints[j].mJointType != JointTypeCustom)
            {
                if(model.mJoints[j].mDoFCount == 1)
                {
                    vw_x_S.setLinearPart(model.S[j].transform_copy(model.bodyFrames[j]->getTransformToRoot()).transform_copy(point_trans).getLinearPart());
                    vw_x_S.setLinearPart(-v_p.getAngularPart().cross(vw_x_S.getLinearPart()));

                    G.col(model.mJoints[j].q_index) = MotionVector(model.S_o[j] + model.v[j]%model.S[j]).transform_copy(model.bodyFrames[j]->getTransformToRoot()).transform_copy(point_trans) +
                            vw_x_S;
                }
                else if(model.mJoints[j].mDoFCount == 3)
                {
                    for(int k = 0;k<3;k++)
                    {
                        vw_x_S.setLinearPart(MotionVector(model.multdof3_S[j].col(k)).transform_copy(model.bodyFrames[j]->getTransformToRoot()).transform_copy(point_trans).getLinearPart());
                        vw_x_S.setLinearPart(-v_p.getAngularPart().cross(vw_x_S.getLinearPart()));
                        G.col(model.mJoints[j].q_index+k) = MotionVector(model.multdof3_S_o[j].col(k) + model.v[j].crossm()*model.multdof3_S[j].col(k)).transform_copy(model.bodyFrames[j]->getTransformToRoot()).transform_copy(point_trans) +
                        vw_x_S;
                    }
                }
            }
            else
            {
                unsigned int k = model.mJoints[j].custom_joint_index;

                for(unsigned int l = 0; l<model.mCustomJoints[k]->mDoFCount;l++)
                {
                    vw_x_S.setLinearPart(MotionVector(model.mCustomJoints[k]->S.col(l)).transform_copy(model.bodyFrames[j]->getTransformToRoot()).transform_copy(point_trans).getLinearPart());
                    vw_x_S.setLinearPart(-v_p.getAngularPart().cross(vw_x_S.getLinearPart()));
                    G.col(model.mJoints[j].q_index+l) = (point_trans * model.bodyFrames[j]->getTransformToRoot()).toMatrix() * (model.mCustomJoints[k]->S_o.col(l) + model.v[j].crossm()*model.mCustomJoints[k]->S.col(l)) +
                            vw_x_S;
                }
            }

            j = model.lambda[j];
        }
    }

    RDL_DLLAPI void calcBodySpatialJacobian(Model &model, const VectorNd &Q, unsigned int body_id, MatrixNd &G,
            bool update_kinematics)
    {
        // update the Kinematics if necessary
        if(update_kinematics)
        {
            updateKinematicsCustom(model, &Q, NULL, NULL);
        }

        assert (G.rows() == 6 && G.cols() == model.qdot_size);

        unsigned int reference_body_id = body_id;

        SpatialTransform base_to_body;
        ReferenceFrame* bodyFrame;

        if(model.IsFixedBodyId(body_id))
        {
            unsigned int fbody_id = body_id - model.fixed_body_discriminator;

            reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
            bodyFrame = model.fixedBodyFrames[fbody_id].get();
        }
        else
        {
            bodyFrame = model.bodyFrames[reference_body_id].get();
        }

        unsigned int j = reference_body_id;

        while(j != 0)
        {
            if(model.mJoints[j].mJointType != JointTypeCustom)
            {
                if(model.mJoints[j].mDoFCount == 1)
                {
                    G.col(model.mJoints[j].q_index) = model.S[j].transform_copy(model.bodyFrames[j]->getTransformToDesiredFrame(bodyFrame));
                }
                else if(model.mJoints[j].mDoFCount == 3)
                {
                    for(int k = 0;k<3;k++)
                    {
                        G.col(model.mJoints[j].q_index+k) = MotionVector(model.multdof3_S[j].col(k)).transform_copy(model.bodyFrames[j]->getTransformToDesiredFrame(bodyFrame));
                    }
                }
            }
            else if(model.mJoints[j].mJointType == JointTypeCustom)
            {
                unsigned int k = model.mJoints[j].custom_joint_index;

                G.block(0, model.mJoints[j].q_index, 6, model.mCustomJoints[k]->mDoFCount) = model.bodyFrames[j]->getTransformToDesiredFrame(bodyFrame).toMatrix() * model.mCustomJoints[k]->S;
            }

            j = model.lambda[j];
        }
    }

    RDL_DLLAPI void calcBodySpatialJacobianDot(Model &model, const Math::VectorNd &Q, const Math::VectorNd QDot, unsigned int body_id, Math::MatrixNd &G, const bool update_kinematics)
    {
        // update the Kinematics if necessary
        if(update_kinematics)
        {
            updateKinematicsCustom(model, &Q, &QDot, NULL);
        }

        assert (G.rows() == 6 && G.cols() == model.qdot_size);

        unsigned int reference_body_id = body_id;

        ReferenceFrame* bodyFrame;

        if(model.IsFixedBodyId(body_id))
        {
            unsigned int fbody_id = body_id - model.fixed_body_discriminator;

            reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
            bodyFrame = model.fixedBodyFrames[fbody_id].get();
        }
        else
        {
            bodyFrame = model.bodyFrames[reference_body_id].get();
        }

        unsigned int j = reference_body_id;

        while(j != 0)
        {
            if(model.mJoints[j].mJointType != JointTypeCustom)
            {
                if(model.mJoints[j].mDoFCount == 1)
                {
                    G.col(model.mJoints[j].q_index) = MotionVector(model.S_o[j] + model.v[j]%model.S[j]).transform_copy(model.bodyFrames[j]->getTransformToDesiredFrame(bodyFrame));
                }
                else if(model.mJoints[j].mDoFCount == 3)
                {
                    for(int k = 0;k<3;k++)
                    {
                        G.col(model.mJoints[j].q_index+k) = MotionVector(model.multdof3_S_o[j].col(k) + model.v[j].crossm()*model.multdof3_S[j].col(k)).transform_copy(model.bodyFrames[j]->getTransformToDesiredFrame(bodyFrame));
                    }
                }
            }
            else if(model.mJoints[j].mJointType == JointTypeCustom)
            {
                unsigned int k = model.mJoints[j].custom_joint_index;

                G.block(0, model.mJoints[j].q_index, 6, model.mCustomJoints[k]->mDoFCount) = model.bodyFrames[j]->getTransformToDesiredFrame(bodyFrame).toMatrix() * (model.mCustomJoints[k]->S_o + model.v[j].crossm()*model.mCustomJoints[k]->S);
            }

            j = model.lambda[j];
        }
    }

    RDL_DLLAPI Math::FrameVector calcPointVelocity(Model &model, const VectorNd &Q, const VectorNd &QDot, unsigned int body_id,
            const Vector3d &point_position, bool update_kinematics)
    {
        assert (model.IsBodyId(body_id));
        assert (model.q_size == Q.size());
        assert (model.qdot_size == QDot.size());

        // Reset the velocity of the root body
        model.v[0].setZero();

        if(update_kinematics)
        {
            updateKinematicsCustom(model, &Q, &QDot, nullptr);
        }

        unsigned int reference_body_id = body_id;
        Math::FramePointd p;

        if(model.IsFixedBodyId(body_id))
        {
            unsigned int fbody_id = body_id - model.fixed_body_discriminator;
            reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;

            p.setIncludingFrame(point_position,model.fixedBodyFrames[fbody_id].get());

            p.changeFrame(model.bodyFrames[reference_body_id].get());
        }
        else
        {
            p.setIncludingFrame(point_position,model.bodyFrames[body_id].get());
        }

        // E*(v-rxw)
        SpatialVector v = model.v[reference_body_id];
        return Math::FrameVector(model.worldFrame.get(),model.bodyFrames[reference_body_id]->getTransformToRoot().E*(v.getLinearPart()-p.vec().cross(v.getAngularPart())));
    }

    RDL_DLLAPI Math::MotionVector calcPointVelocity6D(Model &model, const Math::VectorNd &Q,
            const Math::VectorNd &QDot, unsigned int body_id, const Math::Vector3d &point_position,
            bool update_kinematics)
    {
        assert (model.IsBodyId(body_id));
        assert (model.q_size == Q.size());
        assert (model.qdot_size == QDot.size());

        // Reset the velocity of the root body
        model.v[0].setZero();

        if(update_kinematics)
        {
            updateKinematicsCustom(model, &Q, &QDot, NULL);
        }

        unsigned int reference_body_id = body_id;
        Math::FramePointd p;

        if(model.IsFixedBodyId(body_id))
        {
            unsigned int fbody_id = body_id - model.fixed_body_discriminator;
            reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
            p.setIncludingFrame(point_position,model.fixedBodyFrames[fbody_id].get());
            p.changeFrame(model.bodyFrames[reference_body_id].get());
        }
        else
        {
            p.setIncludingFrame(point_position,model.bodyFrames[body_id].get());
        }

        return MotionVector(SpatialTransform(model.bodyFrames[reference_body_id]->getTransformToRoot().E,p.vec()).apply(model.v[reference_body_id]));
    }

    RDL_DLLAPI Math::FrameVector calcPointAcceleration(Model &model, const VectorNd &Q, const VectorNd &QDot,
            const VectorNd &QDDot, unsigned int body_id, const Vector3d &point_position, bool update_kinematics)
    {

        // Reset the velocity of the root body
        model.v[0].setZero();
        model.a[0].setZero();

        if(update_kinematics)
        {
            updateKinematics(model, Q, QDot, QDDot);
        }

        unsigned int reference_body_id = body_id;
        Math::FramePointd p;

        if(model.IsFixedBodyId(body_id))
        {
            unsigned int fbody_id = body_id - model.fixed_body_discriminator;
            reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;

            p.setIncludingFrame(point_position,model.fixedBodyFrames[fbody_id].get());
            p.changeFrame(model.bodyFrames[reference_body_id].get());
        }
        else
        {
            p.setIncludingFrame(point_position,model.bodyFrames[body_id].get());
        }

        SpatialTransform p_X_i(model.bodyFrames[reference_body_id]->getTransformToRoot().E, p.vec());

        MotionVector p_v_i = model.v[reference_body_id].transform_copy(p_X_i);
        Vector3d a_dash = Vector3d(p_v_i[0], p_v_i[1], p_v_i[2]).cross(Vector3d(p_v_i[3], p_v_i[4], p_v_i[5]));
        MotionVector p_a_i = model.a[reference_body_id].transform_copy(p_X_i);

        return Math::FrameVector(model.worldFrame.get(),p_a_i[3] + a_dash[0], p_a_i[4] + a_dash[1], p_a_i[5] + a_dash[2]);
    }

    RDL_DLLAPI Math::MotionVector calcPointAcceleration6D(Model &model, const VectorNd &Q, const VectorNd &QDot,
            const VectorNd &QDDot, unsigned int body_id, const Vector3d &point_position, bool update_kinematics)
    {
        // Reset the velocity of the root body
        model.v[0].setZero();
        model.a[0].setZero();

        if(update_kinematics)
        {
            updateKinematics(model, Q, QDot, QDDot);
        }

        unsigned int reference_body_id = body_id;
        Math::FramePointd p;

        if(model.IsFixedBodyId(body_id))
        {
            unsigned int fbody_id = body_id - model.fixed_body_discriminator;
            reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;

            p.setIncludingFrame(point_position,model.fixedBodyFrames[fbody_id].get());
            p.changeFrame(model.bodyFrames[reference_body_id]);
        }
        else
        {
            p.setIncludingFrame(point_position,model.bodyFrames[body_id].get());
        }

        SpatialTransform p_X_i(model.bodyFrames[reference_body_id]->getTransformToRoot().E, p.vec());

        MotionVector p_v_i = model.v[reference_body_id].transform_copy(p_X_i);
        Vector3d a_dash = Vector3d(p_v_i[0], p_v_i[1], p_v_i[2]).cross(Vector3d(p_v_i[3], p_v_i[4], p_v_i[5]));
        return MotionVector(model.a[reference_body_id].transform_copy(p_X_i) + MotionVector(0, 0, 0, a_dash[0], a_dash[1], a_dash[2]));
    }

    Math::SpatialMotion calcSpatialVelocity(Model               & model,
                                            const Math::VectorNd& Q,
                                            const Math::VectorNd& QDot,
                                            ReferenceFrame*       body_frame,
                                            ReferenceFrame*       relative_body_frame,
                                            ReferenceFrame*       expressedInFrame,
                                            const bool            update_kinematics)
    {

        ReferenceFrame* expressedInFrameRef;
        expressedInFrameRef = expressedInFrame==nullptr ? body_frame : expressedInFrame;

        if(body_frame==relative_body_frame)
        {
            Math::SpatialMotion m(body_frame,relative_body_frame,body_frame,Math::SpatialVectorZero);
            m.changeFrame(expressedInFrameRef);
            return m;
        }

        if(update_kinematics)
        {
            updateKinematicsCustom(model,&Q,&QDot,nullptr);
        }

        Math::SpatialMotion v_body = model.v[body_frame->getMovableBodyId()];
        if(!body_frame->getIsBodyFrame())
        {
            // Change the frame if the reference frame is a fixed frame
            v_body.changeFrame(body_frame);
            v_body.setBodyFrame(body_frame);
        }

        Math::SpatialMotion v_relative_body=model.v[relative_body_frame->getMovableBodyId()];
        if(!relative_body_frame->getIsBodyFrame())
        {
            // Change the frame if the reference frame is a fixed frame
            v_relative_body.changeFrame(relative_body_frame);
            v_relative_body.setBodyFrame(relative_body_frame);
        }

        v_relative_body.changeFrame(expressedInFrameRef);
        v_body.changeFrame(expressedInFrameRef);

        return (v_body-v_relative_body);
    }

    Math::SpatialMotion calcSpatialVelocity(Model &model, const Math::VectorNd &Q, const Math::VectorNd &QDot, const unsigned int body_id, const unsigned int relative_body_id,
                        ReferenceFrame* expressedInFrame, const bool update_kinematics)
    {
        ReferenceFrame* body_frame;
        ReferenceFrame* relative_body_frame;
        if(model.IsFixedBodyId(body_id))
        {
            body_frame = model.fixedBodyFrames[body_id-model.fixed_body_discriminator].get();
        }
        else
        {
            body_frame = model.bodyFrames[body_id].get();
        }

        ReferenceFrame* expressedInFrameRef;
        expressedInFrameRef = expressedInFrame==nullptr ? body_frame : expressedInFrame;

        if(model.IsFixedBodyId(relative_body_id))
        {
            relative_body_frame = model.fixedBodyFrames[relative_body_id-model.fixed_body_discriminator].get();
        }
        else
        {
            relative_body_frame = model.bodyFrames[relative_body_id].get();
        }

        if(body_id==relative_body_id)
        {
            return Math::SpatialMotion(body_frame,relative_body_frame,body_frame,Math::SpatialVectorZero);
        }

        if(update_kinematics)
        {
            updateKinematicsCustom(model,&Q,&QDot,nullptr);
        }

        Math::SpatialMotion v_body;
        if(model.IsFixedBodyId(body_id))
        {
            v_body = model.v[model.mFixedBodies[body_id-model.fixed_body_discriminator].mMovableParent];
            v_body.changeFrame(model.fixedBodyFrames[body_id-model.fixed_body_discriminator]);
            v_body.setBodyFrame(body_frame);
        }
        else
        {
            v_body = model.v[body_id];
        }

        Math::SpatialMotion v_relative_body;
        if(model.IsFixedBodyId(relative_body_id))
        {
            v_relative_body = model.v[model.mFixedBodies[relative_body_id-model.fixed_body_discriminator].mMovableParent];
            v_relative_body.changeFrame(model.fixedBodyFrames[relative_body_id-model.fixed_body_discriminator]);
            v_relative_body.setBodyFrame(relative_body_frame);
        }
        else
        {
            v_relative_body = model.v[relative_body_id];
        }

        v_relative_body.changeFrame(expressedInFrameRef);
        v_body.changeFrame(expressedInFrameRef);

        return (v_body-v_relative_body);
    }

    Math::SpatialAcceleration calcSpatialAcceleration(Model     & model,
                                                      const Math::VectorNd& Q,
                                                      const Math::VectorNd& QDot,
                                                      const Math::VectorNd& QDDot,
                                                      ReferenceFrame*       body_frame,
                                                      ReferenceFrame*       relative_body_frame,
                                                      ReferenceFrame*       expressedInFrame,
                                                      const bool            update_kinematics)
    {
        if(body_frame==relative_body_frame)
        {
            return Math::SpatialAcceleration(body_frame,relative_body_frame,body_frame,Math::SpatialVectorZero);
        }

        if(update_kinematics)
        {
            updateKinematics(model,Q,QDot,QDDot);
        }

        Math::SpatialAcceleration a_body = model.a[body_frame->getMovableBodyId()];
        if(!body_frame->getIsBodyFrame())
        {
            // There is no relative motion btw fixed body and its moveable parent, so this is kosher
            a_body.changeFrame(body_frame);
            a_body.setBodyFrame(body_frame);
        }

        Math::SpatialAcceleration a_relative_body = model.a[relative_body_frame->getMovableBodyId()];
        if(!relative_body_frame->getIsBodyFrame())
        {
            // @todo - Can add in the velocity thing here to save cycles
            // There is no relative motion btw fixed body and its moveable parent, so this is kosher
            a_relative_body.changeFrame(relative_body_frame);
            a_relative_body.setBodyFrame(relative_body_frame);
        }

        SpatialMotion twistOfCurrentFrameWithRespectToNewFrame =
                calcSpatialVelocity(model,Q,QDot,a_relative_body.getReferenceFrame(),a_body.getReferenceFrame(),a_relative_body.getReferenceFrame(),false);
        SpatialMotion twistOfBodyWrtBaseExpressedInCurrent =
                calcSpatialVelocity(model,Q,QDot,a_relative_body.getBodyFrame(),a_relative_body.getBaseFrame(),a_relative_body.getReferenceFrame(),false);

        a_relative_body.changeFrameWithRelativeMotion(a_body.getReferenceFrame(),twistOfCurrentFrameWithRespectToNewFrame,twistOfBodyWrtBaseExpressedInCurrent);

        a_body-=a_relative_body;

        if(expressedInFrame==nullptr || expressedInFrame==a_body.getReferenceFrame())
        {
            /**
             * If the desired expressedInFrame of the result is body_frame, then only a_relative_body needs to have its frame change.
             */
            return a_body;
        }
        else
        {
            /**
             * User wants the result in a different frame. So change the frame of a_body to the desired frame then
             * return it
             */
            twistOfCurrentFrameWithRespectToNewFrame =
                    calcSpatialVelocity(model,Q,QDot,a_body.getReferenceFrame(),expressedInFrame,a_body.getReferenceFrame(),false);
            twistOfBodyWrtBaseExpressedInCurrent =
                    calcSpatialVelocity(model,Q,QDot,a_body.getBodyFrame(),a_body.getBaseFrame(),a_body.getReferenceFrame(),false);

            a_body.changeFrameWithRelativeMotion(expressedInFrame,twistOfCurrentFrameWithRespectToNewFrame,twistOfBodyWrtBaseExpressedInCurrent);
            return a_body;
        }
    }

    Math::SpatialAcceleration calcSpatialAcceleration(Model &model, const Math::VectorNd &Q, const Math::VectorNd &QDot, const Math::VectorNd &QDDot, const unsigned int body_id, const unsigned int relative_body_id,
                              ReferenceFrame* expressedInFrame, const bool update_kinematics)
    {
        ReferenceFrame* body_frame;
        ReferenceFrame* relative_body_frame;
        if(model.IsFixedBodyId(body_id))
        {
            body_frame = model.fixedBodyFrames[body_id-model.fixed_body_discriminator].get();
        }
        else
        {
            body_frame = model.bodyFrames[body_id].get();
        }

        if(model.IsFixedBodyId(relative_body_id))
        {
            relative_body_frame = model.fixedBodyFrames[relative_body_id-model.fixed_body_discriminator].get();
        }
        else
        {
            relative_body_frame = model.bodyFrames[relative_body_id].get();
        }

        return calcSpatialAcceleration(model,Q,QDot,QDDot,body_frame,relative_body_frame,expressedInFrame,update_kinematics);
    }

    RDL_DLLAPI bool inverseKinematics(Model &model, const VectorNd &Qinit, const std::vector<unsigned int> &body_id,
            const std::vector<Vector3d> &body_point, const std::vector<Vector3d> &target_pos, VectorNd &Qres,
            double step_tol, double lambda, unsigned int max_iter)
    {
        assert (Qinit.size() == model.q_size);
        assert (body_id.size() == body_point.size());
        assert (body_id.size() == target_pos.size());

        MatrixNd J = MatrixNd::Zero(3 * body_id.size(), model.qdot_size);
        VectorNd e = VectorNd::Zero(3 * body_id.size());
        MatrixNd G = model.three_x_qd0;
        Math::FramePointd p;
        Vector3d tmp_vec;

        Qres = Qinit;

        for(unsigned int ik_iter = 0; ik_iter < max_iter; ik_iter++)
        {
            updateKinematicsCustom(model, &Qres, NULL, NULL);

            for(unsigned int k = 0; k < body_id.size(); k++)
            {
                G.setZero();
                calcPointJacobian(model, Qres, body_id[k], body_point[k], G, false);
                p.setIncludingFrame(body_point[k],model.bodyFrames[body_id[k]].get());
                p.changeFrame(model.worldFrame);
                tmp_vec = p.vec();

                for(unsigned int i = 0; i < 3; i++)
                {
                    for(unsigned int j = 0; j < model.qdot_size; j++)
                    {
                        unsigned int row = k * 3 + i;
                        J(row, j) = G(i, j);
                    }

                    e[k * 3 + i] = target_pos[k][i] - tmp_vec[i];
                }
            }

            // abort if we are getting "close"
            if(e.norm() < step_tol)
            {
                return true;
            }

            MatrixNd JJTe_lambda2_I = J * J.transpose() + lambda * lambda * MatrixNd::Identity(e.size(), e.size());

            VectorNd z(body_id.size() * 3);
            bool solve_successful = linSolveGaussElimPivot (JJTe_lambda2_I, e, z);
            assert (solve_successful);

            VectorNd delta_theta = J.transpose() * z;

            Qres = Qres + delta_theta;

            if(delta_theta.norm() < step_tol)
            {
                return true;
            }

            VectorNd test_1(z.size());
            VectorNd test_res(z.size());

            test_1.setZero();

            for(unsigned int i = 0; i < z.size(); i++)
            {
                test_1[i] = 1.;

                VectorNd test_delta = J.transpose() * test_1;

                test_res[i] = test_delta.squaredNorm();

                test_1[i] = 0.;
            }
        }

        return false;
    }

}
