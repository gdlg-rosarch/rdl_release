/*
 * Original Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 *
 * RDL - Robot Dynamics Library
 * Modifications Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef __RDL_UTILS_H__
#define __RDL_UTILS_H__

/**
 * @file rdl_utils.h
 */

#include <string>
#include <rdl_dynamics/rdl_config.h>
#include <rdl_dynamics/rdl_math.hpp>

namespace RobotDynamics
{
    /** \page utils_page Utility Functions
     *
     * All utility functions are specified in the \ref utils_group "Utilities Module".
     *
     * \defgroup utils_group Utilities
     * @{
     *
     * Utility functions are those not necessarily required for kinematics/dynamics, but
     * provide utility by giving model information, calculating useful quantities such as
     * the position/velocity of the center of mass, etc.
     *
     */

    struct Model;

    /** \brief Namespace that contains optional helper functions */
    namespace Utils
    {
        /** @brief get string abbreviation for dof name from spatial vector. */
        std::string getDofName(const Math::SpatialVector& joint_dof);

        /** @brief get body name, returns empty string if bodyid is virtual and has multiple child bodies */
        std::string getBodyName(const RobotDynamics::Model& model, unsigned int body_id);

        /** \brief Creates a human readable overview of the model. */
        std::string getModelHierarchy(const Model& model);

        /** \brief Creates a human readable overview of the Degrees of Freedom. */
        std::string getModelDOFOverview(const Model& model);

        /** \brief Creates a human readable overview of the locations of all bodies that have names. */
        std::string getNamedBodyOriginsOverview(Model& model);

        /** @brief Computes the Center of Mass (COM) position.
         *
         * @param model The model for which we want to compute the COM
         * @param q The current joint positions
         * @param com (output) location of the Center of Mass of the model in world frame
         * @param update_kinematics (optional input) whether the kinematics should be updated (defaults to true)
         */
        void        calcCenterOfMass(Model               & model,
                                     const Math::VectorNd& q,
                                     Math::Vector3d      & com,
                                     bool                  update_kinematics = true);

        /** @brief Computes the Center of Mass (COM) position.
         *
         * @param model The model for which we want to compute the COM
         * @param q The current joint positions
         * @param com (output) location of the Center of Mass of the model in world frame
         * @param update_kinematics (optional input) whether the kinematics should be updated (defaults to true)
         */
        void calcCenterOfMass(Model               & model,
                              const Math::VectorNd& q,
                              Math::FramePointd   & com,
                              bool                  update_kinematics = true);

        /** @brief Computes the Center of Mass (COM) and optionally its linear velocity.
         *
         * When only interested in computing the location of the COM you can use
         * nullptr as value for com_velocity.
         *
         * @param model The model for which we want to compute the COM
         * @param q The current joint positions
         * @param qdot The current joint velocities
         * @param mass (output) total mass of the model
         * @param com (output) location of the Center of Mass of the model in world frame
         * @param com_velocity (optional output) linear velocity of the COM in world frame
         * @param angular_momentum (optional output) angular momentum of the model at the COM in world frame
         * @param update_kinematics (optional input) whether the kinematics should be updated (defaults to true)
         */
        void calcCenterOfMass(Model               & model,
                              const Math::VectorNd& q,
                              const Math::VectorNd& qdot,
                              double              & mass,
                              Math::Vector3d      & com,
                              Math::Vector3d       *com_velocity = NULL,
                              Math::Vector3d       *angular_momentum = NULL,
                              bool                  update_kinematics = true);

        /** @brief Computes the Center of Mass (COM) and optionally its linear velocity and/or angular momentum.
         *
         * When only interested in computing the location of the COM you can use
         * nullptr as value for com_velocity/angular_momentum.
         *
         * @ingroup reference_frame
         *
         * @param model The model for which we want to compute the COM
         * @param q The current joint positions
         * @param qdot The current joint velocities
         * @param mass (output) total mass of the model
         * @param com (output) location of the Center of Mass of the model in world frame
         * @param com_velocity (optional output) linear velocity of the COM in world frame
         * @param angular_momentum (optional output) angular momentum of the model at the COM in a reference frame aligned with
         * the world frame, but located at the center of mass
         * @param update_kinematics (optional input) whether the kinematics should be updated (defaults to true)
         */
        void calcCenterOfMass(Model               & model,
                              const Math::VectorNd& q,
                              const Math::VectorNd& qdot,
                              double              & mass,
                              Math::FramePointd   & com,
                              Math::FrameVector    *com_velocity = nullptr,
                              Math::FrameVector    *angular_momentum = nullptr,
                              bool                  update_kinematics = true);

        /** @brief Computes the Center of Mass (COM) and optionally its linear velocity.
         *
         * When only interested in computing the location of the COM you can use
         * nullptr as value for com_velocity.
         *
         * @ingroup reference_frame
         *
         * @param model The model for which we want to compute the COM
         * @param q The current joint positions
         * @param qdot The current joint velocities
         * @param com (output) location of the Center of Mass of the model in world frame
         * @param com_velocity (optional output) linear velocity of the COM in world frame
         * @param update_kinematics (optional input) whether the kinematics should be updated (defaults to true)
         */
        void calcCenterOfMass(Model               & model,
                              const Math::VectorNd& q,
                              const Math::VectorNd& qdot,
                              Math::FramePointd   & com,
                              Math::FrameVector    *com_velocity = nullptr,
                              bool                  update_kinematics = true);

        /** @brief Computes the potential energy of the full model.
         *
         * @deprecated Use Utils::calcPotentialEnergy instead
         * @param model
         * @param q
         * @param update_kinematics (optional) Defaults to true
         *
         * @return Potential energy
         */
        double calcPotentialEnergy(Model& model, const Math::VectorNd& q, bool update_kinematics = true);

        /** @brief Computes the kinetic energy of the full model.
         *
         * @param model
         * @param q
         * @param qdot
         * @param update_kinematics (optional) Defaults to true
         *
         * @return Kinetic energy
         */
        double calcKineticEnergy(Model               & model,
                                 const Math::VectorNd& q,
                                 const Math::VectorNd& qdot,
                                 bool                  update_kinematics = true);

        /**
         * @brief Computes the matrix \f$J_{com}\f$ such that \f$v_{com} = J_{com}  \dot{q} \f$
         *
         * @param model The model for which the COM jacobian will be computed for
         * @param q The current joint positions
         * @param jCom A 3 x model.qdot_size matrix where the jacobian will be stored.
         * @param update_kinematics If true, kinematic variables will be computed. Default = true.
         */
        void calcCenterOfMassJacobian(Model               & model,
                                      const Math::VectorNd& q,
                                      Math::MatrixNd      & jCom,
                                      bool                  update_kinematics = true);

        /**
         * @brief Calculates the center of mass of a subtree starting with the body with ID bodyId and scales it by the total
         **mass of the subtree.
         *
         * @param model The model used for the calculation
         * @param bodyId The ID of the first body in the desired subtree
         * @param q The current joint positions
         * @param updateKinematics If true, kinematic variables will be computed. Default = true
         */
        Math::FramePointd calcSubtreeCenterOfMassScaledByMass(Model               & model,
                                                              const unsigned int    bodyId,
                                                              const Math::VectorNd& q,
                                                              bool                  updateKinematics = true);

        /**
         * @brief Calculates the total mass of the subtree beginning with body bodyId and traversing outwards from there
         *
         * @param model The model used for the calculation
         * @param bodyId The ID of the first body in the desired subtree
         *
         * @return Mass of subtree
         */
        double calcSubtreeMass(Model& model, const unsigned int bodyId);

        /**
         * @brief Calculates the centroidal momentum matrix, \f$ A(q) \f$ for a model. The centroidal momentum
         * matrix is a \f$ 6 \times N \f$ matrix such that the 6dof centroidal momentum vector is computed by,
         * \f[
         *  h = A(q) \dot{q}
         * \f]
         *
         * @note It is crucial that the \f$ A \f$ matrix be all zeros when this function is called, as the elements will
         * be added to. To be sure, call A.setZero() before calling this function.
         *
         * @param model RDL model
         * @param q Vector of joint positions
         * @param A 6 by N matrix where the result will be stored
         * @param update_kinematics If true, calculates kinematic parameters
         */
        void   calcCentroidalMomentumMatrix(Model               & model,
                                            const Math::VectorNd& q,
                                            Math::MatrixNd      & A,
                                            bool                  update_kinematics = true);

        /**
         * @brief Calculates the time derivative of the centroidal momentum matrix, i.e. the matrix computed by
         * RobotDynamics::Utils::calcCentroidalMomentumMatrix and stores it in the Adot argument
         *
         * @param model
         * @param q
         * @param qdot
         * @param Adot
         * @param update_kinematics
         */
        void calcCentroidalMomentumMatrixDot(Model               & model,
                                             const Math::VectorNd& q,
                                             const Math::VectorNd& qdot,
                                             Math::MatrixNd      & Adot,
                                             bool                  update_kinematics = true);
    }
}

/**
 * @}
 */

#endif // ifndef __RDL_UTILS_H__
