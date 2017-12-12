/*
 * Original Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 *
 * RDL - Robot Dynamics Library
 * Modifications Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef RDL_DYNAMICS_H
#define RDL_DYNAMICS_H

#include <assert.h>
#include <iostream>

#include "rdl_dynamics/rdl_mathutils.h"
#include "rdl_dynamics/Kinematics.h"

namespace RobotDynamics
{
    struct Model;

    /** \page dynamics_page Dynamics
     *
     * All functions related to kinematics are specified in the \ref
     * dynamics_group "Dynamics Module".
     *
     * \defgroup dynamics_group Dynamics
     * @{
     */

    /** \@brief Computes inverse dynamics with the Newton-Euler Algorithm
     *
     * This function computes the generalized forces from given generalized
     * states, velocities, and accelerations:
     *   \f$ \tau = M(q) \ddot{q} + N(q, \dot{q}) \f$
     *
     * @param model rigid body model
     * @param Q     state vector of the internal joints
     * @param QDot  velocity vector of the internal joints
     * @param QDDot accelerations of the internals joints
     * @param Tau   actuations of the internal joints (output)
     * @param f_ext External forces acting on the body in base coordinates (optional, defaults to NULL)
     */
    RDL_DLLAPI void inverseDynamics(
        Model                         & model,
        const Math::VectorNd          & Q,
        const Math::VectorNd          & QDot,
        const Math::VectorNd          & QDDot,
        Math::VectorNd                & Tau,
        std::vector<Math::ForceVector> *f_ext = NULL);

    /** @brief Computes the coriolis forces
     *
     * This function computes the generalized forces from given generalized
     * states, velocities, and accelerations:
     *   \f$ \tau = M(q) \ddot{q} + N(q, \dot{q}) \f$
     *
     * @param model rigid body model
     * @param Q     state vector of the internal joints
     * @param QDot  velocity vector of the internal joints
     * @param Tau   actuations of the internal joints (output)
     */
    RDL_DLLAPI void nonlinearEffects(
        Model& model, const Math::VectorNd& Q, const Math::VectorNd& QDot, Math::VectorNd& Tau);

    /** @brief Computes the joint space inertia matrix by using the Composite Rigid Body Algorithm
     *
     * This function computes the joint space inertia matrix from a given model and
     * the generalized state vector:
     *   \f$ M(q) \f$
     *
     * @param model rigid body model
     * @param Q     state vector of the model
     * @param H     a matrix where the result will be stored in
     * @param update_kinematics  whether the kinematics should be updated (safer, but at a higher computational cost!)
     *
     * @note This function only evaluates the entries of H that are non-zero. One
     * Before calling this function one has to ensure that all other values
     * have been set to zero, e.g. by calling H.setZero().
     */
    RDL_DLLAPI void compositeRigidBodyAlgorithm(
        Model& model, const Math::VectorNd& Q, Math::MatrixNd& H, bool update_kinematics = true);

    /** @brief Computes forward dynamics with the Articulated Body Algorithm
     *
     * This function computes the generalized accelerations from given
     * generalized states, velocities and forces:
     *   \f$ \ddot{q} = M(q)^{-1} ( -N(q, \dot{q}) + \tau)\f$
     * It does this by using the recursive Articulated Body Algorithm that runs
     * in \f$O(n_{dof})\f$ with \f$n_{dof}\f$ being the number of joints.
     *
     * @param model rigid body model
     * @param Q     state vector of the internal joints
     * @param QDot  velocity vector of the internal joints
     * @param Tau   actuations of the internal joints
     * @param QDDot accelerations of the internal joints (output)
     * @param f_ext External forces acting on the body in base coordinates (optional, defaults to NULL)
     */
    RDL_DLLAPI void forwardDynamics(
        Model                         & model,
        const Math::VectorNd          & Q,
        const Math::VectorNd          & QDot,
        const Math::VectorNd          & Tau,
        Math::VectorNd                & QDDot,
        std::vector<Math::ForceVector> *f_ext = NULL);

    /** @brief Computes forward dynamics by building and solving the full Lagrangian equation
     *
     * This method builds and solves the linear system
     * \f[      H \ddot{q} = -C + \tau	\f]
     * for \f$\ddot{q}\f$ where \f$H\f$ is the joint space inertia matrix
     * computed with the CompositeRigidBodyAlgorithm(), \f$C\f$ the bias
     * force (sometimes called "non-linear effects").
     *
     * @param model rigid body model
     * @param Q     state vector of the internal joints
     * @param QDot  velocity vector of the internal joints
     * @param Tau   actuations of the internal joints
     * @param QDDot accelerations of the internal joints (output)
     * @param linear_solver specification which method should be used for solving the linear system
     * @param f_ext External forces acting on the body in base coordinates (optional, defaults to NULL)
     * @param H     preallocated workspace area for the joint space inertia matrix of size dof_count x dof_count (optional,
     * defaults to NULL and allocates temporary matrix)
     * @param C     preallocated workspace area for the right hand side vector of size dof_count x 1 (optional, defaults to NULL
     * and allocates temporary vector)
     */
    RDL_DLLAPI void forwardDynamicsLagrangian(
        Model                         & model,
        const Math::VectorNd          & Q,
        const Math::VectorNd          & QDot,
        const Math::VectorNd          & Tau,
        Math::VectorNd                & QDDot,
        Math::LinearSolver              linear_solver = Math::LinearSolverColPivHouseholderQR,
        std::vector<Math::ForceVector> *f_ext = nullptr,
        Math::MatrixNd                 *H = nullptr,
        Math::VectorNd                 *C = nullptr);

    /** @brief Computes the effect of multiplying the inverse of the joint
     * space inertia matrix with a vector in linear time.
     *
     * @param model rigid body model
     * @param Q     state vector of the generalized positions
     * @param Tau   the vector that should be multiplied with the inverse of
     *              the joint space inertia matrix
     * @param QDDot vector where the result will be stored
     * @param update_kinematics whether the kinematics should be updated (safer, but at a higher computational cost)
     *
     * This function uses a reduced version of the Articulated %Body Algorithm
     * to compute
     *
     *   \f$ \ddot{q} = M(q)^{-1}*\tau \f$
     *
     * in \f$O(n_{\textit{dof}}\f$) time.
     *
     * \note The first time calling this function you MUST supply update_kinematics=true.
     * Subsequent calls within the same kinematics tic though should supply update_kinematics=false
     * to avoids expensive recomputations of transformations and articulated body inertias.
     */
    RDL_DLLAPI void calcMInvTimesTau(
        Model& model, const Math::VectorNd& Q, const Math::VectorNd& Tau, Math::VectorNd& QDDot, bool update_kinematics = true);

    /** @} */
}

/* RDL_DYNAMICS_H */
#endif // ifndef RDL_DYNAMICS_H
