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
#include <assert.h>

#include "rdl_dynamics/rdl_mathutils.h"
#include "rdl_dynamics/Model.h"
#include "rdl_dynamics/Contacts.h"
#include "rdl_dynamics/Dynamics.h"

namespace RobotDynamics
{

    using namespace Math;

    unsigned int ConstraintSet::addConstraint(unsigned int body_id, const Vector3d &body_point,
            const Vector3d &world_normal, const char *constraint_name, double normal_acceleration)
    {
        assert (bound == false);

        std::string name_str;
        if(constraint_name != NULL)
        {
            name_str = constraint_name;
        }

        name.push_back(name_str);
        body.push_back(body_id);
        point.push_back(body_point);
        normal.push_back(world_normal);

        unsigned int n_constr = acceleration.size() + 1;

        acceleration.conservativeResize(n_constr);
        acceleration[n_constr - 1] = normal_acceleration;

        force.conservativeResize(n_constr);
        force[n_constr - 1] = 0.;

        impulse.conservativeResize(n_constr);
        impulse[n_constr - 1] = 0.;

        v_plus.conservativeResize(n_constr);
        v_plus[n_constr - 1] = 0.;

        d_multdof3_u = std::vector<Math::Vector3d>(n_constr, Math::Vector3d::Zero());

        return n_constr - 1;
    }

    bool ConstraintSet::bind(const Model &model)
    {
        assert (bound == false);

        if(bound)
        {
            std::cerr << "Error: binding an already bound constraint set!" << std::endl;
            abort();
        }
        unsigned int n_constr = size();

        H = model.ndof0_mat;

        C = model.ndof0_vec;

        gamma.conservativeResize(n_constr);
        gamma.setZero();
        G.conservativeResize(n_constr, model.dof_count);
        G.setZero();
        A.conservativeResize(model.dof_count + n_constr, model.dof_count + n_constr);
        A.setZero();
        b.conservativeResize(model.dof_count + n_constr);
        b.setZero();
        x.conservativeResize(model.dof_count + n_constr);
        x.setZero();

        GT_qr = Eigen::HouseholderQR<Math::MatrixNd>(G.transpose());
        GT_qr_Q = model.ndof0_mat;
        Y = MatrixNd::Zero(model.dof_count, G.rows());
        Z = MatrixNd::Zero(model.dof_count, model.dof_count - G.rows());
        qddot_y = model.ndof0_vec;
        qddot_z = model.ndof0_vec;

        K.conservativeResize(n_constr, n_constr);
        K.setZero();
        a.conservativeResize(n_constr);
        a.setZero();
        QDDot_t = model.ndof0_vec;
        QDDot_0 = model.ndof0_vec;

        f_t.resize(n_constr, SpatialVectorZero);
        f_ext_constraints.resize(model.mBodies.size(), SpatialVectorZero);
        point_accel_0.resize(n_constr, Vector3d::Zero());

        d_pA = std::vector<SpatialVector>(model.mBodies.size(), SpatialVectorZero);
        d_a = std::vector<SpatialVector>(model.mBodies.size(), SpatialVectorZero);
        d_u = model.nbodies0_vec;

        d_IA = std::vector<SpatialMatrix>(model.mBodies.size(), SpatialMatrixIdentity);
        d_U = std::vector<SpatialVector>(model.mBodies.size(), SpatialVectorZero);
        d_d = model.nbodies0_vec;

        d_multdof3_u = std::vector<Math::Vector3d>(model.mBodies.size(), Math::Vector3d::Zero());

        bound = true;

        return bound;
    }

    void ConstraintSet::clear()
    {
        acceleration.setZero();
        force.setZero();
        impulse.setZero();

        H.setZero();
        C.setZero();
        gamma.setZero();
        G.setZero();
        A.setZero();
        b.setZero();
        x.setZero();

        K.setZero();
        a.setZero();
        QDDot_t.setZero();
        QDDot_0.setZero();

        unsigned int i;
        for(i = 0; i < f_t.size(); i++)
        {
            f_t[i].setZero();
        }

        for(i = 0; i < f_ext_constraints.size(); i++)
        {
            f_ext_constraints[i].setZero();
        }

        for(i = 0; i < point_accel_0.size(); i++)
        {
            point_accel_0[i].setZero();
        }

        for(i = 0; i < d_pA.size(); i++)
        {
            d_pA[i].setZero();
        }

        for(i = 0; i < d_a.size(); i++)
        {
            d_a[i].setZero();
        }

        d_u.setZero();
    }

    RDL_DLLAPI void solveContactSystemDirect(Math::MatrixNd &H, const Math::MatrixNd &G, const Math::VectorNd &c,
            const Math::VectorNd &gamma, Math::VectorNd &qddot, Math::VectorNd &lambda, Math::MatrixNd &A,
            Math::VectorNd &b, Math::VectorNd &x, Math::LinearSolver &linear_solver)
    {
        // Build the system: Copy H
        A.block(0, 0, c.rows(), c.rows()) = H;

        // Copy G and G^T
        A.block(0, c.rows(), c.rows(), gamma.rows()) = G.transpose();
        A.block(c.rows(), 0, gamma.rows(), c.rows()) = G;

        // Build the system: Copy -C + \tau
        b.block(0, 0, c.rows(), 1) = c;
        b.block(c.rows(), 0, gamma.rows(), 1) = gamma;

        switch(linear_solver)
        {
            case (LinearSolverPartialPivLU) :
                x = A.partialPivLu().solve(b);
                break;
            case (LinearSolverColPivHouseholderQR) :
                x = A.colPivHouseholderQr().solve(b);
                break;
            case (LinearSolverHouseholderQR) :
                x = A.householderQr().solve(b);
                break;
            default:
                assert (0);
                break;
        }
    }

    RDL_DLLAPI void solveContactSystemRangeSpaceSparse(Model &model, Math::MatrixNd &H, const Math::MatrixNd &G,
            const Math::VectorNd &c, const Math::VectorNd &gamma, Math::VectorNd &qddot, Math::VectorNd &lambda,
            Math::MatrixNd &K, Math::VectorNd &a, Math::LinearSolver linear_solver)
    {
        SparseFactorizeLTL(model, H);

        MatrixNd Y(G.transpose());

        for(unsigned int i = 0; i < Y.cols(); i++)
        {
            VectorNd Y_col = Y.block(0, i, Y.rows(), 1);
            SparseSolveLTx(model, H, Y_col);
            Y.block(0, i, Y.rows(), 1) = Y_col;
        }

        VectorNd z(c);
        SparseSolveLTx(model, H, z);

        K = Y.transpose() * Y;

        a = gamma - Y.transpose() * z;

        lambda = K.llt().solve(a);

        qddot = c + G.transpose() * lambda;
        SparseSolveLTx(model, H, qddot);
        SparseSolveLx(model, H, qddot);
    }

    RDL_DLLAPI
    void solveContactSystemNullSpace(Math::MatrixNd &H, const Math::MatrixNd &G, const Math::VectorNd &c,
            const Math::VectorNd &gamma, Math::VectorNd &qddot, Math::VectorNd &lambda, Math::MatrixNd &Y,
            Math::MatrixNd &Z, Math::VectorNd &qddot_y, Math::VectorNd &qddot_z, Math::LinearSolver &linear_solver)
    {
        switch(linear_solver)
        {
            case (LinearSolverPartialPivLU) :
                qddot_y = (G * Y).partialPivLu().solve(gamma);
                break;
            case (LinearSolverColPivHouseholderQR) :
                qddot_y = (G * Y).colPivHouseholderQr().solve(gamma);
                break;
            case (LinearSolverHouseholderQR) :
                qddot_y = (G * Y).householderQr().solve(gamma);
                break;
            default:
                assert (0);
                break;
        }

        qddot_z = (Z.transpose() * H * Z).llt().solve(Z.transpose() * (c - H * Y * qddot_y));

        qddot = Y * qddot_y + Z * qddot_z;

        switch(linear_solver)
        {
            case (LinearSolverPartialPivLU) :
                lambda = (G * Y).partialPivLu().solve(Y.transpose() * (H * qddot - c));
                break;
            case (LinearSolverColPivHouseholderQR) :
                lambda = (G * Y).colPivHouseholderQr().solve(Y.transpose() * (H * qddot - c));
                break;
            case (LinearSolverHouseholderQR) :
                lambda = (G * Y).householderQr().solve(Y.transpose() * (H * qddot - c));
                break;
            default:
                assert (0);
                break;
        }
    }

    RDL_DLLAPI
    void calcContactJacobian(Model &model, const Math::VectorNd &Q, const ConstraintSet &CS, Math::MatrixNd &G,
            bool update_kinematics)
    {
        if(update_kinematics)
        {
            updateKinematicsCustom(model, &Q, NULL, NULL);
        }

        unsigned int i, j;

        // variables to check whether we need to recompute G
        unsigned int prev_body_id = 0;
        Vector3d prev_body_point = Vector3d::Zero();
        MatrixNd Gi(3, model.dof_count);

        for(i = 0; i < CS.size(); i++)
        {
            // only compute the matrix Gi if actually needed
            if(prev_body_id != CS.body[i] || prev_body_point != CS.point[i])
            {
                Gi.setZero();
                calcPointJacobian(model, Q, CS.body[i], CS.point[i], Gi, false);
                prev_body_id = CS.body[i];
                prev_body_point = CS.point[i];
            }

            for(j = 0; j < model.dof_count; j++)
            {
                Vector3d gaxis(Gi(0, j), Gi(1, j), Gi(2, j));
                G(i, j) = gaxis.transpose() * CS.normal[i];
            }
        }
    }

    RDL_DLLAPI
    void calcContactSystemVariables(Model &model, const Math::VectorNd &Q, const Math::VectorNd &QDot,
            const Math::VectorNd &Tau, ConstraintSet &CS)
    {
        // Compute C
        nonlinearEffects(model, Q, QDot, CS.C);
        assert (CS.H.cols() == model.dof_count && CS.H.rows() == model.dof_count);

        // Compute H
        compositeRigidBodyAlgorithm(model, Q, CS.H, false);

        calcContactJacobian(model, Q, CS, CS.G, false);

        // Compute gamma
        unsigned int prev_body_id = 0;
        Vector3d prev_body_point = Vector3d::Zero();
        Vector3d gamma_i = Vector3d::Zero();

        CS.QDDot_0.setZero();
        updateKinematicsCustom(model, NULL, NULL, &CS.QDDot_0);

        for(unsigned int i = 0; i < CS.size(); i++)
        {
            // only compute point accelerations when necessary
            if(prev_body_id != CS.body[i] || prev_body_point != CS.point[i])
            {
                gamma_i = calcPointAcceleration(model, Q, QDot, CS.QDDot_0, CS.body[i], CS.point[i], false);
                prev_body_id = CS.body[i];
                prev_body_point = CS.point[i];
            }

            // we also substract ContactData[i].acceleration such that the contact
            // point will have the desired acceleration
            CS.gamma[i] = CS.acceleration[i] - CS.normal[i].dot(gamma_i);
        }
    }

    RDL_DLLAPI
    void forwardDynamicsContactsDirect(Model &model, const VectorNd &Q, const VectorNd &QDot, const VectorNd &Tau,
            ConstraintSet &CS, VectorNd &QDDot)
    {
        calcContactSystemVariables(model, Q, QDot, Tau, CS);

        solveContactSystemDirect(CS.H, CS.G, Tau - CS.C, CS.gamma, QDDot, CS.force, CS.A, CS.b, CS.x, CS.linear_solver);

        // Copy back QDDot
        for(unsigned int i = 0; i < model.dof_count; i++)
        {
            QDDot[i] = CS.x[i];
        }

        // Copy back contact forces
        for(unsigned int i = 0; i < CS.size(); i++)
        {
            CS.force[i] = -CS.x[model.dof_count + i];
        }
    }

    RDL_DLLAPI
    void forwardDynamicsContactsRangeSpaceSparse(Model &model, const Math::VectorNd &Q, const Math::VectorNd &QDot,
            const Math::VectorNd &Tau, ConstraintSet &CS, Math::VectorNd &QDDot)
    {
        calcContactSystemVariables(model, Q, QDot, Tau, CS);

        solveContactSystemRangeSpaceSparse(model, CS.H, CS.G, Tau - CS.C, CS.gamma, QDDot, CS.force, CS.K, CS.a, CS.linear_solver);
    }

    RDL_DLLAPI
    void forwardDynamicsContactsNullSpace(Model &model, const VectorNd &Q, const VectorNd &QDot, const VectorNd &Tau,
            ConstraintSet &CS, VectorNd &QDDot)
    {
        calcContactSystemVariables(model, Q, QDot, Tau, CS);

        CS.GT_qr.compute(CS.G.transpose());
        CS.GT_qr.householderQ().evalTo(CS.GT_qr_Q);

        CS.Y = CS.GT_qr_Q.block(0, 0, QDot.rows(), CS.G.rows());
        CS.Z = CS.GT_qr_Q.block(0, CS.G.rows(), QDot.rows(), QDot.rows() - CS.G.rows());

        solveContactSystemNullSpace(CS.H, CS.G, Tau - CS.C, CS.gamma, QDDot, CS.force, CS.Y, CS.Z, CS.qddot_y, CS.qddot_z, CS.linear_solver);
    }

    RDL_DLLAPI
    void computeContactImpulsesDirect(Model &model, const Math::VectorNd &Q, const Math::VectorNd &QDotMinus,
            ConstraintSet &CS, Math::VectorNd &QDotPlus)
    {
        // Compute H
        updateKinematicsCustom(model, &Q, NULL, NULL);
        compositeRigidBodyAlgorithm(model, Q, CS.H, false);

        // Compute G
        calcContactJacobian(model, Q, CS, CS.G, false);

        solveContactSystemDirect(CS.H, CS.G, CS.H * QDotMinus, CS.v_plus, QDotPlus, CS.impulse, CS.A, CS.b, CS.x, CS.linear_solver);

        // Copy back QDotPlus
        for(unsigned int i = 0; i < model.dof_count; i++)
        {
            QDotPlus[i] = CS.x[i];
        }

        // Copy back constraint impulses
        for(unsigned int i = 0; i < CS.size(); i++)
        {
            CS.impulse[i] = CS.x[model.dof_count + i];
        }
    }

    RDL_DLLAPI
    void computeContactImpulsesRangeSpaceSparse(Model &model, const Math::VectorNd &Q, const Math::VectorNd &QDotMinus,
            ConstraintSet &CS, Math::VectorNd &QDotPlus)
    {
        // Compute H
        updateKinematicsCustom(model, &Q, NULL, NULL);
        compositeRigidBodyAlgorithm(model, Q, CS.H, false);

        // Compute G
        calcContactJacobian(model, Q, CS, CS.G, false);

        solveContactSystemRangeSpaceSparse(model, CS.H, CS.G, CS.H * QDotMinus, CS.v_plus, QDotPlus, CS.impulse, CS.K, CS.a, CS.linear_solver);
    }

    RDL_DLLAPI
    void computeContactImpulsesNullSpace(Model &model, const Math::VectorNd &Q, const Math::VectorNd &QDotMinus,
            ConstraintSet &CS, Math::VectorNd &QDotPlus)
    {
        // Compute H
        updateKinematicsCustom(model, &Q, NULL, NULL);
        compositeRigidBodyAlgorithm(model, Q, CS.H, false);

        // Compute G
        calcContactJacobian(model, Q, CS, CS.G, false);

        CS.GT_qr.compute(CS.G.transpose());
        CS.GT_qr_Q = CS.GT_qr.householderQ();

        CS.Y = CS.GT_qr_Q.block(0, 0, QDotMinus.rows(), CS.G.rows());
        CS.Z = CS.GT_qr_Q.block(0, CS.G.rows(), QDotMinus.rows(), QDotMinus.rows() - CS.G.rows());

        solveContactSystemNullSpace(CS.H, CS.G, CS.H * QDotMinus, CS.v_plus, QDotPlus, CS.impulse, CS.Y, CS.Z, CS.qddot_y, CS.qddot_z, CS.linear_solver);
    }

    /** \brief Compute only the effects of external forces on the generalized accelerations
     *
     * This function is a reduced version of ForwardDynamics() which only
     * computes the effects of the external forces on the generalized
     * accelerations.
     *
     */
    RDL_DLLAPI
    void forwardDynamicsApplyConstraintForces(Model &model, const VectorNd &Tau, ConstraintSet &CS, VectorNd &QDDot)
    {
        assert (QDDot.size() == model.dof_count);

        unsigned int i = 0;

        for(i = 1; i < model.mBodies.size(); i++)
        {
            model.I[i].setSpatialMatrix(model.IA[i]);
            model.pA[i] = model.v[i]%Momentum(model.I[i],model.v[i]);

            if(CS.f_ext_constraints[i] != SpatialVectorZero)
            {
                model.pA[i] -= model.bodyFrames[i]->getInverseTransformToRoot().toMatrixAdjoint() * CS.f_ext_constraints[i];
            }
        }

        for(i = model.mBodies.size() - 1; i > 0; i--)
        {
            unsigned int q_index = model.mJoints[i].q_index;
            ReferenceFrame* bodyFrame = model.bodyFrames[i].get();

            if(model.mJoints[i].mDoFCount == 3 && model.mJoints[i].mJointType != JointTypeCustom)
            {
                unsigned int lambda = model.lambda[i];
                model.multdof3_u[i] = Vector3d(Tau[q_index], Tau[q_index + 1], Tau[q_index + 2]) - model.multdof3_S[i].transpose() * model.pA[i];

                if(lambda != 0)
                {
                    SpatialMatrix Ia = model.IA[i] - (model.multdof3_U[i] * model.multdof3_Dinv[i] * model.multdof3_U[i].transpose());

                    SpatialVector pa = model.pA[i] + Ia * model.c[i] + model.multdof3_U[i] * model.multdof3_Dinv[i] * model.multdof3_u[i];

                    model.IA[lambda].noalias() += (bodyFrame->getTransformFromParent().toMatrixTranspose() * Ia * bodyFrame->getTransformFromParent().toMatrix());
                    model.pA[lambda].noalias() += bodyFrame->getTransformFromParent().applyTranspose(pa);
                }
            }
            else if(model.mJoints[i].mDoFCount == 1 && model.mJoints[i].mJointType != JointTypeCustom)
            {
                model.u[i] = Tau[q_index] - model.S[i].dot(model.pA[i]);

                unsigned int lambda = model.lambda[i];
                if(lambda != 0)
                {
                    SpatialMatrix Ia = model.IA[i] - model.U[i] * (model.U[i] / model.d[i]).transpose();
                    SpatialVector pa = model.pA[i] + Ia * model.c[i] + model.U[i] * model.u[i] / model.d[i];

                    model.IA[lambda].noalias() += (bodyFrame->getTransformFromParent().toMatrixTranspose() * Ia * bodyFrame->getTransformFromParent().toMatrix());
                    model.pA[lambda].noalias() += bodyFrame->getTransformFromParent().applyTranspose(pa);
                }
            }
            else if(model.mJoints[i].mJointType == JointTypeCustom)
            {

                unsigned int kI = model.mJoints[i].custom_joint_index;
                unsigned int dofI = model.mCustomJoints[kI]->mDoFCount;
                unsigned int lambda = model.lambda[i];
                VectorNd tau_temp = VectorNd::Zero(dofI);

                for(unsigned int z = 0; z < dofI; ++z)
                {
                    tau_temp[z] = Tau[q_index + z];
                }

                model.mCustomJoints[kI]->u = tau_temp - (model.mCustomJoints[kI]->S.transpose() * model.pA[i]);

                if(lambda != 0)
                {
                    SpatialMatrix Ia = model.IA[i] - (model.mCustomJoints[kI]->U * model.mCustomJoints[kI]->Dinv * model.mCustomJoints[kI]->U.transpose());
                    SpatialVector pa = model.pA[i] + Ia * model.c[i] + (model.mCustomJoints[kI]->U * model.mCustomJoints[kI]->Dinv * model.mCustomJoints[kI]->u);

                    model.IA[lambda].noalias() += bodyFrame->getTransformFromParent().toMatrixTranspose() * Ia * bodyFrame->getTransformFromParent().toMatrix();
                    model.pA[lambda].noalias() += bodyFrame->getTransformFromParent().applyTranspose(pa);
                }
            }
        }

        model.a[0].set(SpatialVector(-model.gravity));

        for(i = 1; i < model.mBodies.size(); i++)
        {
            unsigned int q_index = model.mJoints[i].q_index;
            unsigned int lambda = model.lambda[i];

            model.a[i].set(model.a[lambda].transform_copy(model.bodyFrames[i]->getTransformFromParent()) + model.c[i]);

            if(model.mJoints[i].mDoFCount == 3 && model.mJoints[i].mJointType != JointTypeCustom)
            {
                Vector3d qdd_temp = model.multdof3_Dinv[i] * (model.multdof3_u[i] - model.multdof3_U[i].transpose() * model.a[i]);

                QDDot[q_index] = qdd_temp[0];
                QDDot[q_index + 1] = qdd_temp[1];
                QDDot[q_index + 2] = qdd_temp[2];
                model.a[i].set(model.a[i] + model.multdof3_S[i] * qdd_temp);
            }
            else if(model.mJoints[i].mDoFCount == 1 && model.mJoints[i].mJointType != JointTypeCustom)
            {
                QDDot[q_index] = (1. / model.d[i]) * (model.u[i] - model.U[i].dot(model.a[i]));
                model.a[i].set(model.a[i] + model.S[i] * QDDot[q_index]);
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

                model.a[i].set(model.a[i] + (model.mCustomJoints[kI]->S * qdd_temp));
            }
        }
    }

    /** \brief Computes the effect of external forces on the generalized accelerations.
     *
     * This function is essentially similar to ForwardDynamics() except that it
     * tries to only perform computations of variables that change due to
     * external forces defined in f_t.
     */
    RDL_DLLAPI
    void forwardDynamicsAccelerationDeltas(Model &model, ConstraintSet &CS, VectorNd &QDDot_t,
            const unsigned int body_id, const std::vector<SpatialVector> &f_t)
    {
        assert (CS.d_pA.size() == model.mBodies.size());
        assert (CS.d_a.size() == model.mBodies.size());
        assert (CS.d_u.size() == model.mBodies.size());

        // TODO reset all values (debug)
        for(unsigned int i = 0; i < model.mBodies.size(); i++)
        {
            CS.d_pA[i].setZero();
            CS.d_a[i].setZero();
            CS.d_u[i] = 0.;
            CS.d_multdof3_u[i].setZero();
        }
        for(unsigned int i = 0; i < model.mCustomJoints.size(); i++)
        {
            model.mCustomJoints[i]->d_u.setZero();
        }

        for(unsigned int i = body_id; i > 0; i--)
        {
            if(i == body_id)
            {
                CS.d_pA[i] = -model.bodyFrames[i]->getInverseTransformToRoot().applyAdjoint(f_t[i]);
            }

            if(model.mJoints[i].mDoFCount == 3 && model.mJoints[i].mJointType != JointTypeCustom)
            {
                CS.d_multdof3_u[i] = -model.multdof3_S[i].transpose() * (CS.d_pA[i]);

                unsigned int lambda = model.lambda[i];
                if(lambda != 0)
                {
                    CS.d_pA[lambda] = CS.d_pA[lambda] + model.bodyFrames[i]->getTransformFromParent().applyTranspose(CS.d_pA[i] + (model.multdof3_U[i] * model.multdof3_Dinv[i] * CS.d_multdof3_u[i]));
                }
            }
            else if(model.mJoints[i].mDoFCount == 1 && model.mJoints[i].mJointType != JointTypeCustom)
            {
                CS.d_u[i] = -model.S[i].dot(CS.d_pA[i]);
                unsigned int lambda = model.lambda[i];

                if(lambda != 0)
                {
                    CS.d_pA[lambda] = CS.d_pA[lambda] + model.bodyFrames[i]->getTransformFromParent().applyTranspose(CS.d_pA[i] + model.U[i] * CS.d_u[i] / model.d[i]);
                }
            }
            else if(model.mJoints[i].mJointType == JointTypeCustom)
            {

                unsigned int kI = model.mJoints[i].custom_joint_index;
                //CS.
                model.mCustomJoints[kI]->d_u = -model.mCustomJoints[kI]->S.transpose() * (CS.d_pA[i]);
                unsigned int lambda = model.lambda[i];
                if(lambda != 0)
                {
                    CS.d_pA[lambda] = CS.d_pA[lambda] + model.bodyFrames[i]->getTransformFromParent().applyTranspose(CS.d_pA[i] + (model.mCustomJoints[kI]->U * model.mCustomJoints[kI]->Dinv * model.mCustomJoints[kI]->d_u));
                }
            }
        }

        QDDot_t[0] = 0.;
        CS.d_a[0] = model.a[0];

        for(unsigned int i = 1; i < model.mBodies.size(); i++)
        {
            unsigned int q_index = model.mJoints[i].q_index;
            unsigned int lambda = model.lambda[i];

            SpatialVector Xa = model.bodyFrames[i]->getTransformFromParent().apply(CS.d_a[lambda]);

            if(model.mJoints[i].mDoFCount == 3 && model.mJoints[i].mJointType != JointTypeCustom)
            {
                Vector3d qdd_temp = model.multdof3_Dinv[i] * (CS.d_multdof3_u[i] - model.multdof3_U[i].transpose() * Xa);

                QDDot_t[q_index] = qdd_temp[0];
                QDDot_t[q_index + 1] = qdd_temp[1];
                QDDot_t[q_index + 2] = qdd_temp[2];
                model.a[i].set(model.a[i] + model.multdof3_S[i] * qdd_temp);
                CS.d_a[i] = Xa + model.multdof3_S[i] * qdd_temp;
            }
            else if(model.mJoints[i].mDoFCount == 1 && model.mJoints[i].mJointType != JointTypeCustom)
            {

                QDDot_t[q_index] = (CS.d_u[i] - model.U[i].dot(Xa)) / model.d[i];
                CS.d_a[i] = Xa + model.S[i] * QDDot_t[q_index];
            }
            else if(model.mJoints[i].mJointType == JointTypeCustom)
            {
                unsigned int kI = model.mJoints[i].custom_joint_index;
                unsigned int dofI = model.mCustomJoints[kI]->mDoFCount;

                VectorNd qdd_temp = model.mCustomJoints[kI]->Dinv * (model.mCustomJoints[kI]->d_u - model.mCustomJoints[kI]->U.transpose() * Xa);

                for(unsigned int z = 0; z < dofI; ++z)
                {
                    QDDot_t[q_index + z] = qdd_temp[z];
                }

                model.a[i].set(model.a[i] + model.mCustomJoints[kI]->S * qdd_temp);
                CS.d_a[i] = Xa + model.mCustomJoints[kI]->S * qdd_temp;
            }
        }
    }

    inline void set_zero(std::vector<SpatialVector> &spatial_values)
    {
        for(unsigned int i = 0; i < spatial_values.size(); i++)
        {
            spatial_values[i].setZero();
        }
    }

    RDL_DLLAPI
    void forwardDynamicsContactsKokkevis(Model &model, const VectorNd &Q, const VectorNd &QDot, const VectorNd &Tau,
            ConstraintSet &CS, VectorNd &QDDot)
    {
        assert (CS.f_ext_constraints.size() == model.mBodies.size());
        assert (CS.QDDot_0.size() == model.dof_count);
        assert (CS.QDDot_t.size() == model.dof_count);
        assert (CS.f_t.size() == CS.size());
        assert (CS.point_accel_0.size() == CS.size());
        assert (CS.K.rows() == CS.size());
        assert (CS.K.cols() == CS.size());
        assert (CS.force.size() == CS.size());
        assert (CS.a.size() == CS.size());

        Vector3d point_accel_t;

        unsigned int ci = 0;

        // The default acceleration only needs to be computed once
        forwardDynamics(model, Q, QDot, Tau, CS.QDDot_0);

        // we have to compute the standard accelerations first as we use them to
        // compute the effects of each test force
        for(ci = 0; ci < CS.size(); ci++)
        {
            unsigned int body_id = CS.body[ci];
            Vector3d point = CS.point[ci];
            Vector3d normal = CS.normal[ci];
            double acceleration = CS.acceleration[ci];

            updateKinematicsCustom(model, NULL, NULL, &CS.QDDot_0);
            CS.point_accel_0[ci] = calcPointAcceleration(model, Q, QDot, CS.QDDot_0, body_id, point, false);

            CS.a[ci] = -acceleration + normal.dot(CS.point_accel_0[ci]);
        }

        // Now we can compute and apply the test forces and use their net effect
        // to compute the inverse articlated inertia to fill K.
        Math::FramePointd p;
        for(ci = 0; ci < CS.size(); ci++)
        {
            unsigned int body_id = CS.body[ci];
            Vector3d point = CS.point[ci];
            Vector3d normal = CS.normal[ci];

            unsigned int movable_body_id = body_id;
            Vector3d point_global;
            if(model.IsFixedBodyId(body_id))
            {
                unsigned int fbody_id = body_id - model.fixed_body_discriminator;
                movable_body_id = model.mFixedBodies[fbody_id].mMovableParent;
                p.setIncludingFrame(point,model.fixedBodyFrames[fbody_id].get());
                p.changeFrame(model.worldFrame.get());
                point_global = p.vec();
            }
            else
            {
                p.setIncludingFrame(point,model.bodyFrames[body_id].get());
                p.changeFrame(model.worldFrame.get());
                point_global = p.vec();
            }

            CS.f_t[ci] = SpatialTransform(Matrix3d::Identity(), -point_global).applyAdjoint(SpatialVector(0., 0., 0., -normal[0], -normal[1], -normal[2]));
            CS.f_ext_constraints[movable_body_id] = CS.f_t[ci];

            forwardDynamicsAccelerationDeltas(model, CS, CS.QDDot_t, movable_body_id, CS.f_ext_constraints);

            CS.f_ext_constraints[movable_body_id].setZero();

            CS.QDDot_t += CS.QDDot_0;

            // compute the resulting acceleration
            updateKinematicsCustom(model, NULL, NULL, &CS.QDDot_t);

            for(unsigned int cj = 0; cj < CS.size(); cj++)
            {
                point_accel_t = calcPointAcceleration(model, Q, QDot, CS.QDDot_t, CS.body[cj], CS.point[cj], false);

                CS.K(ci, cj) = CS.normal[cj].dot(point_accel_t - CS.point_accel_0[cj]);
            }
        }

        bool solve_successful = linSolveGaussElimPivot (CS.K, CS.a, CS.force);
        assert (solve_successful);

        for(ci = 0; ci < CS.size(); ci++)
        {
            unsigned int body_id = CS.body[ci];
            unsigned int movable_body_id = body_id;

            if(model.IsFixedBodyId(body_id))
            {
                unsigned int fbody_id = body_id - model.fixed_body_discriminator;
                movable_body_id = model.mFixedBodies[fbody_id].mMovableParent;
            }

            CS.f_ext_constraints[movable_body_id] -= CS.f_t[ci] * CS.force[ci];
        }

        forwardDynamicsApplyConstraintForces(model, Tau, CS, QDDot);
    }

} /* namespace RobotDynamics */
