/*
 * Original Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 *
 * RDL - Robot Dynamics Library
 * Modifications Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef RDL_KINEMATICS_H
#define RDL_KINEMATICS_H

/**
 * @file Kinematics.h
 */

#include <assert.h>
#include <iostream>
#include "rdl_dynamics/Model.h"

namespace RobotDynamics
{
    /** \page kinematics_page Kinematics
     * All functions related to kinematics are specified in the \ref
     * kinematics_group "Kinematics Module".
     *
     * \note Please note that in the Robot Dynamics Library all angles
     * are specified in radians.
     *
     * \defgroup kinematics_group Kinematics
     * @{
     *
     * \note Please note that in the Robot Dynamics Library all angles
     * are specified in radians.
     */

    /** @brief Updates and computes velocities and accelerations of the bodies
     *
     * This function updates the kinematic variables such as body velocities, accelerations,
     * and reference frames in the model to reflect the variables passed to this function.
     *
     * @param model the model
     * @param Q     the positional variables of the model
     * @param QDot  the generalized velocities of the joints
     * @param QDDot the generalized accelerations of the joints
     */
    RDL_DLLAPI void updateKinematics(Model               & model,
                                     const Math::VectorNd& Q,
                                     const Math::VectorNd& QDot,
                                     const Math::VectorNd& QDDot);

    /** @brief Selectively updates model internal states of body positions, velocities and/or accelerations.
     *
     * This function updates the kinematic variables such as body velocities,
     * accelerations, and reference frames in the model to reflect the variables passed to this function.
     *
     * In contrast to RobotDynamics::updateKinematics() this function allows to update the model
     * state with values one is interested and thus reduce computations (e.g. only
     * positions, only positions + accelerations, only velocities, etc.).

     * @param model the model
     * @param Q     the positional variables of the model
     * @param QDot  the generalized velocities of the joints
     * @param QDDot the generalized accelerations of the joints
     */
    RDL_DLLAPI void updateKinematicsCustom(Model               & model,
                                           const Math::VectorNd *Q,
                                           const Math::VectorNd *QDot,
                                           const Math::VectorNd *QDDot);

    /** @brief Computes the point jacobian for a point on a body
     *
     * If a position of a point is computed by a function \f$g(q(t))\f$ for which its
     * time derivative is \f$\frac{d}{dt} g(q(t)) = G(q)\dot{q}\f$ then this
     * function computes the jacobian matrix \f$G(q)\f$.
     *
     * @param model   rigid body model
     * @param Q       state vector of the internal joints
     * @param body_id the id of the body
     * @param point_position the position of the point in body-local data
     * @param G       a matrix of dimensions 3 x \#qdot_size where the result will be stored in
     * @param update_kinematics whether UpdateKinematics() should be called or not (default: true)
     *
     * The result will be returned via the G argument.
     *
     * @note This function only evaluates the entries of G that are non-zero. One
     * Before calling this function one has to ensure that all other values
     * have been set to zero, e.g. by calling G.setZero().
     *
     */
    RDL_DLLAPI void calcPointJacobian(Model               & model,
                                      const Math::VectorNd& Q,
                                      unsigned int          body_id,
                                      const Math::Vector3d& point_position,
                                      Math::MatrixNd      & G,
                                      bool                  update_kinematics = true);

    /**
     * @brief Computes the jacobian of a frame \f$J_{i}^{k,j}\f$ with \f$i\f$ being the "base" frame, \f$j\f$ being the
     *"relative" frame, and \f$k\f$ being
     * the "expressed in" frame. Multiplying this jacobian by a vector of joint velocities will result in the spatial motion of
     *the baseFrame w.r.t
     * relativeFrame expressed in expressedInFrame, a.k.a \f$v_{i}^{k,j} = J_{i}^{k,j}\dot{q}\f$
     *
     * @param model
     * @param Q
     * @param G
     * @param baseFrame
     * @param relativeFrame
     * @param expressedInFrame
     * @param update_kinematics
     *
     * @note This function only modifies the elements of the jacobian that are nonzero, so be sure the other elements are
     * not nonzero because those elements will not be zero'd. Best practice would be to call G.setZero() before calling
     * this function
     *
     * @note <b>If expressedInFrame=nullptr then the expressedInFrame will be set to baseFrame</b>
     */
    void calcRelativeBodySpatialJacobian(Model               & model,
                                         const Math::VectorNd& Q,
                                         Math::MatrixNd      & G,
                                         ReferenceFrame       *baseFrame,
                                         ReferenceFrame       *relativeFrame,
                                         ReferenceFrame       *expressedInFrame = nullptr,
                                         bool                  update_kinematics = true);

    /**
     * @brief Computes the rate of change of the jacobian, \f$\dot{J}_{i}^{k,j}\f$, with \f$i\f$ being the "base" frame,
     * \f$j\f$ being the relative" frame, and \f$k\f$ being the "expressed in" frame. This jacobian is such that the
     * following is true, \f$a_{i}^{k,j} = J_{i}^{k,j}\ddot{q} + \dot{J}_{i}^{k,j}\dot{q}\f$ where
     * \f$a_{i}^{k,j}\f$ is the spatial acceleration of frame \f$i\f$ w.r.t frame \f$j\f$, expressed in
     * frame \f$k\f$. Additionally the jacobian \f$J_{i}^{k,j}\f$ can be computed by
     * RobotDynamics::calcRelativeBodySpatialJacobian.
     * @param model
     * @param Q
     * @param QDot
     * @param G
     * @param baseFrame
     * @param relativeFrame
     * @param expressedInFrame
     * @param update_kinematics
     *
     * * @note <b>If expressedInFrame=nullptr then the expressedInFrame will be set to baseFrame</b>
     */
    void calcRelativeBodySpatialJacobianDot(Model & model,
                                            const Math::VectorNd& Q,
                                            const Math::VectorNd& QDot,
                                            Math::MatrixNd &G, ReferenceFrame* baseFrame,
                                            ReferenceFrame* relativeFrame,
                                            ReferenceFrame* expressedInFrame=nullptr,
                                            bool update_kinematics=true);


    /**
     * @brief Computes the time derivative of the linear components the a point jacobian on a body.
     *
     * For a point \f$p_{i}\f$ on body \f$i\f$ expressed in body \f$i\f$'s frame, this function
     * computes \f$ \dot{J} \f$ such that \f$a_{i} = J(q)\ddot{q} + \dot{J}(q,\dot{q})\dot{q}\f$, where
     * \f$a_{i}\f$ is the 3D linear acceleration of point \f$p\f$, and \f$a_{i}\f$ is expressed in
     * world frame.
     *
     * @param model Rigid body model
     * @param Q Joint positions
     * @param QDot Joint velocities
     * @param body_id Id of the body
     * @param point_position 3D position on body
     * @param G Matrix where the result is stored
     * @param update_kinematics Defaults to true. If true, kinematics will be calculated. If kinematics have already
     * been calculated, setting this to false will save time
     *
     * @note This function only evaluates the entries of G that are non-zero. One
     * Before calling this function one has to ensure that all other values
     * have been set to zero, e.g. by calling G.setZero().
     */
    RDL_DLLAPI void calcPointJacobianDot(Model               & model,
                                         const Math::VectorNd& Q,
                                         const Math::VectorNd& QDot,
                                         unsigned int          body_id,
                                         const Math::Vector3d& point_position,
                                         Math::MatrixNd      & G,
                                         bool                  update_kinematics = true);

    /** @brief Computes a 6-D Jacobian for a point on a body
     *
     * Computes the 6-D Jacobian \f$G(q)\f$ that when multiplied with
     * \f$\dot{q}\f$ gives a 6-D vector that has the angular velocity as the
     * first three entries and the linear velocity as the last three entries.
     *
     * @param model   rigid body model
     * @param Q       state vector of the internal joints
     * @param body_id the id of the body
     * @param point_position the position of the point in body-local data
     * @param G       a matrix of dimensions 6 x \#qdot_size where the result will be stored in
     * @param update_kinematics whether UpdateKinematics() should be called or not (default: true)
     *
     * The result will be returned via the G argument.
     *
     * @note This function only evaluates the entries of G that are non-zero. One
     * Before calling this function one has to ensure that all other values
     * have been set to zero, e.g. by calling G.setZero().
     *
     */
    RDL_DLLAPI void calcPointJacobian6D(Model               & model,
                                        const Math::VectorNd& Q,
                                        unsigned int          body_id,
                                        const Math::Vector3d& point_position,
                                        Math::MatrixNd      & G,
                                        bool                  update_kinematics = true);

    /**
     * @brief Computes the 6D time derivative of a point jacobian on a body.
     *
     * For a point \f$ p_{i} \f$ on body \f$i\f$ expressed in body \f$ i \f$'s frame, this function
     * computes \f$ \dot{J} \f$ such that \f$ a_{i} = J(q)\ddot{q} + \dot{J}(q,\dot{q})\dot{q}\f$, where
     * \f$a_{i}\f$ is the 6D acceleration of point \f$p\f$, and \f$a_{i}\f$ is expressed in
     * world frame.
     *
     * Essentially, this method computes the jacobian of a reference frame located at point \f$p_{i}\f$, but aligned
     * with the world frame.
     *
     * @param model Rigid body model
     * @param Q Joint positions
     * @param QDot Joint velocities
     * @param body_id Id of the body
     * @param point_position 3D position on body
     * @param G Matrix where the result is stored
     * @param update_kinematics Defaults to true. If true, kinematics will be calculated. If kinematics have already
     * been calculated, setting this to false will save time
     *
     * @note This function only evaluates the entries of G that are non-zero. One
     * Before calling this function one has to ensure that all other values
     * have been set to zero, e.g. by calling G.setZero().
     */
    RDL_DLLAPI void calcPointJacobianDot6D(Model               & model,
                                           const Math::VectorNd& Q,
                                           const Math::VectorNd& QDot,
                                           unsigned int          body_id,
                                           const Math::Vector3d& point_position,
                                           Math::MatrixNd      & G,
                                           bool                  update_kinematics = true);

    /** @brief Computes the spatial jacobian for a body.
     * The result will be returned via the G argument and represents the
     * body Jacobian expressed at the origin of the body. The corresponding
     * spatial velocity of the body w.r.t the world frame expressed in
     * body frame can be calculated, for the \f$ ith \f$ body, as
     * \f$ v^{i,0}_{i} = G(q) \dot{q} \f$.
     *
     * @param model   rigid body model
     * @param Q       state vector of the internal joints
     * @param body_id the id of the body
     * @param G       a matrix of size 6 x \#qdot_size where the result will be stored in
     * @param update_kinematics whether UpdateKinematics() should be called or not (default: true)
     *
     * @note This function only evaluates the entries of G that are non-zero. One
     * Before calling this function one has to ensure that all other values
     * have been set to zero, e.g. by calling G.setZero().
     */
    RDL_DLLAPI void calcBodySpatialJacobian(
        Model& model, const Math::VectorNd& Q, unsigned int body_id, Math::MatrixNd& G, bool update_kinematics = true);

    /** @brief Computes the time derivative of the spatial jacobian for a body.
     * The result will be returned via the G argument and represents the time
     * derivative of the body Jacobian expressed at the origin of the body. The corresponding
     * spatial acceleration of the body w.r.t the world frame expressed in
     * body frame can be calculated, for the \f$i\f$th body, as
     * \f$ a^{i,0}_{i} = G(q) \ddot{q} + \dot{G}(q)\dot{q} \f$ where \f$G(q)\f$
     * is the body jacobian of body \f$i\f$.
     *
     * @param model   rigid body model
     * @param Q       joint positions
     * @param QDot    joint velocities
     * @param body_id the id of the body
     * @param G       a matrix of size 6 x \#qdot_size where the result will be stored in
     * @param update_kinematics whether UpdateKinematics() should be called or not (default: true)
     *
     * @note This function only evaluates the entries of G that are non-zero. One
     * Before calling this function one has to ensure that all other values
     * have been set to zero, e.g. by calling G.setZero().
     */
    RDL_DLLAPI void calcBodySpatialJacobianDot(Model               & model,
                                               const Math::VectorNd& Q,
                                               const Math::VectorNd  QDot,
                                               unsigned int          body_id,
                                               Math::MatrixNd      & G,
                                               const bool            update_kinematics = true);

    /** @brief Computes the velocity of a point on a body
     *
     * @param model   rigid body model
     * @param Q       state vector of the internal joints
     * @param QDot    velocity vector of the internal joints
     * @param body_id the id of the body
     * @param point_position the position of the point in body-local data
     * @param update_kinematics whether UpdateKinematics() should be called or not (default: true)
     *
     * @returns A FrameVector representing the points velocity in world frame.
     */
    RDL_DLLAPI Math::FrameVector calcPointVelocity(
        Model               & model,
        const Math::VectorNd& Q,
        const Math::VectorNd& QDot,
        unsigned int          body_id,
        const Math::Vector3d& point_position,
        bool                  update_kinematics = true);

    /** @brief Computes angular and linear velocity of a point that is fixed on a body
     *
     * @param model   rigid body model
     * @param Q       state vector of the internal joints
     * @param QDot    velocity vector of the internal joints
     * @param body_id the id of the body
     * @param point_position the position of the point in body-local data
     * @param update_kinematics whether UpdateKinematics() should be called or not (default: true)
     *
     * @returns The a 6-D vector for which the first three elements are the
     * angular velocity and the last three elements the linear velocity in the
     * global reference system.
     */
    RDL_DLLAPI
    Math::MotionVector calcPointVelocity6D(
        Model               & model,
        const Math::VectorNd& Q,
        const Math::VectorNd& QDot,
        unsigned int          body_id,
        const Math::Vector3d& point_position,
        bool                  update_kinematics = true);

    /** @brief Computes the linear acceleration of a point on a body
     *
     * @param model   rigid body model
     * @param Q       state vector of the internal joints
     * @param QDot    velocity vector of the internal joints
     * @param QDDot    velocity vector of the internal joints
     * @param body_id the id of the body
     * @param point_position the position of the point in body-local data
     * @param update_kinematics whether UpdateKinematics() should be called or not (default: true)
     *
     * @returns The cartesian acceleration of the point in world frame (output)
     *
     * The kinematic state of the model has to be updated before valid
     * values can be obtained. This can either be done by calling
     * UpdateKinematics() or setting the last parameter update_kinematics to
     * true (default).
     *
     * @warning  If this function is called after ForwardDynamics() without
     * an update of the kinematic state one has to add the gravity
     * acceleration has to be added to the result.
     */
    RDL_DLLAPI
    Math::FrameVector calcPointAcceleration(
        Model               & model,
        const Math::VectorNd& Q,
        const Math::VectorNd& QDot,
        const Math::VectorNd& QDDot,
        unsigned int          body_id,
        const Math::Vector3d& point_position,
        bool                  update_kinematics = true);

    /**
     * @brief Compute the spatial velocity of any frame with respect to any other frame, expressed in an arbirtary third frame. The returned
     **RobotDynamicS::Math::SpatialMotion
     * is expressed in body_frame unless th expressedInFrame is provided. Each time RobotDynamics::updateKinematics is called, the spatial velocity
     **of
     * each body with respect to the world, and expressed in body frame is calculated. For body \f$ i \f$ this is written
     * as \f$ v^{i,0}_{i} \f$. Given another body, body \f$ j \f$, the velocity of body \f$ i \f$ relative to body \f$ j \f$ and
     * expressed in body frame \f$ i \f$ is computed by,
     * \f[
     *  v^{i,j}_{i} = v^{i,0}_{i} - ^{i}X_{j} v^{j,0}_{j}
     * \f]
     *
     * @param model A RDL robot model
     * @param Q Vector of joint positions
     * @param QDot Vector of joint velocities
     * @param body_frame The primary frame
     * @param relative_body_frame The frame the result will be w.r.t
     * @param expressedInFrame Frame the result will be expressed in. Default value=nullptr in which case the result will be expressed in the frame corresponding to body_id
     * @param update_kinematics
     * @return
     *
     * @note If no expressedInFrame is given then the result will be expressed in the frame corresponding to the arg body_id.
     */
    Math::SpatialMotion calcSpatialVelocity(Model               & model,
                                            const Math::VectorNd& Q,
                                            const Math::VectorNd& QDot,
                                            ReferenceFrame*       body_frame,
                                            ReferenceFrame*       relative_body_frame,
                                            ReferenceFrame*       expressedInFrame=nullptr,
                                            const bool            update_kinematics = true);

    /**
     * @brief Compute the spatial velocity of any body with respect to any other body. The returned
     **RobotDynamicS::Math::SpatialMotion
     * is expressed in the body frame of body body_id. Each time RobotDynamics::updateKinematics is called, the spatial velocity
     **of
     * each body with respect to the world, and expressed in body frame is calculated. For body \f$ i \f$ this is written
     * as \f$ v^{i,0}_{i} \f$. Given another body, body \f$ j \f$, the velocity of body \f$ i \f$ relative to body \f$ j \f$ and
     * expressed in body frame \f$ i \f$ is computed by,
     * \f[
     *  v^{i,j}_{i} = v^{i,0}_{i} - ^{i}X_{j} v^{j,0}_{j}
     * \f]
     *
     * @param model A RDL robot model
     * @param Q Vector of joint positions
     * @param QDot Vector of joint velocities
     * @param body_id The
     * @param relative_body_id
     * @param expressedInFrame Frame the result will be expressed in. Default value=nullptr in which case the result will be expressed in the frame corresponding to body_id
     * @param update_kinematics
     * @return
     *
     * @note If no expressedInFrame is given then the result will be expressed in the frame corresponding to the arg body_id.
     */
    Math::SpatialMotion calcSpatialVelocity(Model               & model,
                                            const Math::VectorNd& Q,
                                            const Math::VectorNd& QDot,
                                            const unsigned int    body_id,
                                            const unsigned int    relative_body_id,
                                            ReferenceFrame*       expressedInFrame=nullptr,
                                            const bool            update_kinematics = true);

    /**
     * @brief Compute the spatial acceleration of any body with respect to any other body and express the
     * result in an arbitrary reference frame. The returned RobotDynamicS::Math::SpatialAcceleration
     * can be expressed in the reference frame of choice by supplying the expressedInFrame argument. If left to
     * the default value of nullptr, the result will be expressed in the reference frame corresponding to body_id.
     *
     * @note If the expressedInFrame arg is nullptr, the resulting acceleration will be expressed in the referenceFrame corresponding to body_id
     *
     * @note that in this notation \f$\widetilde{v}\f$ is the same as the cross operator,\f$v\times \f$, commonly used in RBD
     * by Featherstone.
     *
     * @param model A RDL robot model
     * @param Q Vector of joint positions
     * @param QDot Vector of joint velocities
     * @param QDDot Vector of joint accelerations
     * @param body_id The
     * @param relative_body_id
     * @param update_kinematics
     * @return
     */
    Math::SpatialAcceleration calcSpatialAcceleration(Model     & model,
                                            const Math::VectorNd& Q,
                                            const Math::VectorNd& QDot,
                                            const Math::VectorNd& QDDot,
                                            const unsigned int    body_id,
                                            const unsigned int    relative_body_id,
                                            ReferenceFrame*       expressedInFrame=nullptr,
                                            const bool            update_kinematics = true);

    /**
     * @brief Compute the spatial acceleration of any frame with respect to any other frame and express the
     * result in an arbitrary third frame. The returned RobotDynamicS::Math::SpatialAcceleration
     * can be expressed in the reference frame of choice by supplying the expressedInFrame argument. If left to
     * the default value of nullptr, the result will be expressed in body_frame.
     *
     * @note If the expressedInFrame arg is nullptr, the resulting acceleration will be expressed in body_frame
     *
     * @note that in this notation \f$\widetilde{v}\f$ is the same as the cross operator,\f$v\times \f$, commonly used in RBD
     * by Featherstone.
     *
     * @param model A RDL robot model
     * @param Q Vector of joint positions
     * @param QDot Vector of joint velocities
     * @param QDDot Vector of joint accelerations
     * @param body_id The
     * @param relative_body_id
     * @param update_kinematics
     * @return
     */
    Math::SpatialAcceleration calcSpatialAcceleration(Model     & model,
                                                      const Math::VectorNd& Q,
                                                      const Math::VectorNd& QDot,
                                                      const Math::VectorNd& QDDot,
                                                      ReferenceFrame*       body_frame,
                                                      ReferenceFrame*       relative_body_frame,
                                                      ReferenceFrame*       expressedInFrame=nullptr,
                                                      const bool            update_kinematics = true);

    /** @brief Computes linear and angular acceleration of a point on a body
     *
     * @param model   rigid body model
     * @param Q       state vector of the internal joints
     * @param QDot    velocity vector of the internal joints
     * @param QDDot    velocity vector of the internal joints
     * @param body_id the id of the body
     * @param point_position the position of the point in body-local data
     * @param update_kinematics whether UpdateKinematics() should be called or not (default: true)
     *
     * @returns A SpatialAcceleration of the desired point with respect to world expressed in world frame
     *
     * The kinematic state of the model has to be updated before valid
     * values can be obtained. This can either be done by calling
     * updateKinematics() or setting the last parameter update_kinematics to
     * true (default).
     *
     * @warning  If this function is called after ForwardDynamics() without
     * an update of the kinematic state one has to add the gravity
     * acceleration has to be added to the result.
     */
    RDL_DLLAPI
    Math::MotionVector calcPointAcceleration6D(
        Model               & model,
        const Math::VectorNd& Q,
        const Math::VectorNd& QDot,
        const Math::VectorNd& QDDot,
        unsigned int          body_id,
        const Math::Vector3d& point_position,
        bool                  update_kinematics = true);

    /** \brief Computes the inverse kinematics iteratively using a damped Levenberg-Marquardt method (also known as Damped Least
       Squares method)
     *
     * @param model rigid body model
     * @param Qinit initial guess for the state
     * @param body_id a vector of all bodies for which we we have kinematic target positions
     * @param body_point a vector of points in body local coordinates that are
     * to be matched to target positions
     * @param target_pos a vector of target positions
     * @param Qres output of the computed inverse kinematics
     * @param step_tol tolerance used for convergence detection
     * @param lambda damping factor for the least squares function
     * @param max_iter maximum number of steps that should be performed
     * @returns true on success, false otherwise
     *
     * This function repeatedly computes
     *   \f[ Qres = Qres + \Delta \theta\f]
     *   \f[ \Delta \theta = G^T (G^T G + \lambda^2 I)^{-1} e \f]
     * where \f$G = G(q) = \frac{d}{dt} g(q(t))\f$ and \f$e\f$ is the
     * correction of the body points so that they coincide with the target
     * positions. The function returns true when \f$||\Delta \theta||_2 \le\f$
     * step_tol or if the error between body points and target gets smaller
     * than step_tol. Otherwise it returns false.
     *
     * The parameter \f$\lambda\f$ is the damping factor that has to
     * be chosen carefully. In case of unreachable positions higher values (e.g
     * 0.9) can be helpful. Otherwise values of 0.0001, 0.001, 0.01, 0.1 might
     * yield good results. See the literature for best practices.
     *
     * @warning The actual accuracy might be rather low (~1.0e-2)! Use this function with a
     * grain of suspicion.
     */
    RDL_DLLAPI
    bool inverseKinematics(
        Model                            & model,
        const Math::VectorNd             & Qinit,
        const std::vector<unsigned int>  & body_id,
        const std::vector<Math::Vector3d>& body_point,
        const std::vector<Math::Vector3d>& target_pos,
        Math::VectorNd                   & Qres,
        double                             step_tol = 1.0e-12,
        double                             lambda = 0.01,
        unsigned int                       max_iter = 50);

    /** @} */
}

/* RDL_KINEMATICS_H */
#endif // ifndef RDL_KINEMATICS_H
