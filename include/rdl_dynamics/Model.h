/*
 * Original Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 *
 * RDL - Robot Dynamics Library
 * Modifications Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef RDL_MODEL_H
#define RDL_MODEL_H

/**
 * @file Model.h
 */

#include "rdl_dynamics/rdl_math.hpp"
#include <map>
#include <list>
#include <assert.h>
#include <iostream>
#include <limits>
#include <cstring>

#include "rdl_dynamics/Joint.h"
#include "rdl_dynamics/Body.h"
#include "rdl_dynamics/SpatialAlgebraOperators.h"
#include "rdl_dynamics/Quaternion.h"
#include "rdl_dynamics/RdlExceptions.hpp"

#include <type_traits>

// std::vectors containing any objects that have Eigen matrices or vectors
// as members need to have a special allocater. This can be achieved with
// the following macro.

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RobotDynamics::Joint);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RobotDynamics::Body);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RobotDynamics::FixedBody);

/** @brief Namespace for all structures of the RobotDynamics library
 */
namespace RobotDynamics
{
    /** \page modeling_page Model
     *
     * \section model_structure Model Structure
     *
     * RDL stores the model internally in the \link RobotDynamics::Model
     * Model Structure\endlink. For each \link RobotDynamics::Body Body
     * \endlink it contains spatial velocities, accelerations and other
     * variables that describe the state of the rigid body system. Furthermore
     * it contains variables that are used as temporary variables in the
     * algorithms.
     *
     * The construction of \link RobotDynamics::Model Model Structures
     * \endlink makes use of carefully designed constructors of the classes
     * \link RobotDynamics::Body Body \endlink and \link
     * RobotDynamics::Joint Joint \endlink to ease the process of
     * creating articulated models.
     *
     * \link RobotDynamics::Body Bodies \endlink are created by calling one
     * of its constructors. Usually they are created by specifying the mass,
     * center of mass and the inertia at the center of mass.
     * \link RobotDynamics::Joint Joints \endlink are similarly created and is
     * described in detail in \ref joint_description.
     *
     * Adding bodies to the model is done by specifying the
     * parent body by its id, the transformation from the parent origin to the
     * joint origin, the joint specification as an object, and the body itself.
     * These parameters are then fed to the function
     * RobotDynamics::Model::addBody() or
     * RobotDynamics::Model::appendBody().
     *
     * To create a model with a floating base (a.k.a a model with a free-flyer
     * joint) it is recommended to use \link
     * RobotDynamics::Model::SetFloatingBaseBody
     * Model::SetFloatingBaseBody(...)\endlink.
     *
     * Once this is done, the model structure can be used with the functions of \ref
     * kinematics_group, \ref dynamics_group, \ref contacts_page, to perform
     * computations.
     *
     * After creating a model, access the body frames for moving bodies in the
     * RobotDynamics::Model:bodyFrames vector and for fixed bodies in
     * RobotDynamics::Model::fixedBodyFrames.
     *
     * \section modeling_urdf Using URDF
     *
     * For this see the documentation see \link
     * RobotDynamics::Urdf::urdfReadFromFile \endlink.
     */

    /** @brief Contains all information about the rigid body model
     *
     * This class contains all information required to perform the forward
     * dynamics calculation. The variables in this class are also used for
     * storage of temporary values. It is designed for use of the Articulated
     * Rigid Body Algorithm (which is implemented in ForwardDynamics()) and
     * follows the numbering as described in Featherstones book.
     *
     * Please note that body 0 is the root body and the moving bodies start at
     * index 1. This numbering scheme is very beneficial in terms of
     * readability of the code as the resulting code is very similar to the
     * pseudo-code in the RBDA book. The generalized variables q, qdot, qddot
     * and tau however start at 0 such that the first entry (e.g. q[0]) always
     * specifies the value for the first moving body.
     *
     *
     * @todo - Need to implement proper copy constructor, so I can use raw pointers and can clean them up properly
     *
     *
     * \note To query the number of degrees of freedom use Model::dof_count.
     */
    struct RDL_DLLAPI Model
    {
        Model();

        ~Model()
        {}

        // Structural information

        std::shared_ptr<ReferenceFrame>worldFrame;                       /**< Pointer to world frame */

        std::vector<std::shared_ptr<ReferenceFrame> >bodyFrames;         /**< Reference frames for each body. Frame names are
                                                                            the same as the body name */

        std::vector<std::shared_ptr<ReferenceFrame> >bodyCenteredFrames; /**< Body centered frames are aligned with the body
                                                                            frame, but located at a body's center of mass */

        std::vector<std::shared_ptr<ReferenceFrame> >fixedBodyFrames;

        std::vector<unsigned int>lambda;                     /**< The id of the parents body */
        std::vector<std::vector<unsigned int>> lambda_chain; /**< Each body's chain of parents */

        std::vector<unsigned int>lambda_q;         /**< The index of the parent degree of freedom that is directly influencing
                                                      the current one */

        std::vector<std::vector<unsigned int> >mu; /**< Contains the ids of all the children of a given body */

        /** @brief number of degrees of freedoms of the model
         *
         * This value contains the number of entries in the generalized state (q)
         * velocity (qdot), acceleration (qddot), and force (tau) vector.
         */
        unsigned int dof_count;

        /** @brief The size of the \f$\mathbf{q}\f$-vector.
         * For models without spherical joints the value is the same as
         * Model::dof_count, otherwise additional values for the w-component of the
         * Quaternion is stored at the end of \f$\mathbf{q}\f$.
         *
         * \sa \ref joint_description for more details.
         */
        unsigned int q_size;

        /** @brief The size of the
         *
         * (\f$\mathbf{\dot{q}}, \mathbf{\ddot{q}}\f$,
         * and \f$\mathbf{\tau}\f$-vector.
         *
         * \sa \ref joint_description for more details.
         */
        unsigned int qdot_size;

        /// @brief Id of the previously added body, required for Model::appendBody()
        unsigned int previously_added_body_id;

        /// @brief the cartesian vector of the gravity
        Math::MotionVector gravity;

        // State information
        /// @brief The spatial velocity of the bodies
        std::vector<Math::SpatialAcceleration>a; /**< The spatial acceleration of the bodies expressed in body frame */

        std::vector<Math::SpatialMotion>v;       /**< The spatial velocity of the bodies expressed in body frame */
        Math::SpatialMotion             tmpSpatialMotion;
        Math::MotionVector              tmpMotionVector;

        Math::FramePointd tmpFramePoint;

        Math::SpatialAcceleration tmpSpatialAcceleration; /**< A temporary spatial acceleration */

        /**
         * Temporary transform for use in intermediate calculations
         */
        Math::SpatialTransform X_tmp;

        ////////////////////////////////////
        // Joints

        /// @brief All joints

        std::vector<Joint>mJoints;

        std::vector<Math::MotionVector>S;   /**< Motion subspace for 1 dof joints */
        std::vector<Math::MotionVector>S_o; /**< The ring derivative of S, a.k.a the rate of growth */

        // Joint state variables
        std::vector<Math::SpatialTransform>X_J;

        std::vector<Math::SpatialVector>c_J; /**< Apparent derivative of v_J in a coordinate moving with the successor body */

        std::vector<Math::SpatialMotion>v_J; /**< Joints velocity w.r.t its parent joint, expressed in body frame */

        std::vector<unsigned int>mJointUpdateOrder;

        /// @brief Transformations from the parent body to the frame of the joint.
        // It is expressed in the coordinate frame of the parent.
        std::vector<Math::SpatialTransform>X_T;

        /// @brief The number of fixed joints that have been declared before
        ///  each joint.
        std::vector<unsigned int>mFixedJointCount;

        ////////////////////////////////////
        // Special variables for joints with 3 degrees of freedom
        /// @brief Motion subspace for joints with 3 degrees of freedom

        std::vector<Math::Matrix63>multdof3_S;
        std::vector<Math::Matrix63>multdof3_S_o;
        std::vector<Math::Matrix63>multdof3_U;
        std::vector<Math::Matrix3d>multdof3_Dinv;
        std::vector<Math::Vector3d>multdof3_u;
        std::vector<unsigned int>  multdof3_w_index;

        std::vector<CustomJoint *>mCustomJoints;

        ////////////////////////////////////
        // Dynamics variables

        /// @brief The velocity dependent spatial acceleration
        std::vector<Math::SpatialVector>c;

        /// @brief The spatial inertia of the bodies
        std::vector<Math::SpatialMatrix>IA;

        /// @brief The spatial bias force
        std::vector<Math::SpatialVector>pA;

        /// @brief Temporary variable U_i (RBDA p. 130)
        std::vector<Math::SpatialVector>U;

        /// @brief Temporary variable D_i (RBDA p. 130)
        Math::VectorNd d;

        /// @brief Temporary variable u (RBDA p. 130)
        Math::VectorNd u;

        /// @brief Internal forces on the body (used only InverseDynamics())
        std::vector<Math::SpatialVector>f;
        std::vector<Math::ForceVector>  f_b;

        std::vector<Math::SpatialInertia>  I;  /**< Body inertia expressed in body frame */
        std::vector<Math::RigidBodyInertia>Ic; /**< A body inertia. See RobotDynamics::compositeRigidBodyAlgorithm for the uses
                                                  of this */

        std::vector<Math::SpatialInertia>Ib_c; /**< Ib_c is the inertia of (non-fixed) bodies expressed in the body's center of
                                                  mass frame. */
        std::vector<Math::SpatialVector>hc;

        /**
         * @note It's CRITICAL that the elements of the zero matrices
         * never be modified! It will break absolutely everything.
         */
        Math::MatrixNd ndof0_mat;    /**< ndof x ndof matrix of zeros */
        Math::MatrixNd three_x_qd0;  /**< 3 x qdot_size matrix of zeros */
        Math::VectorNd q0_vec;       /**< q_size x 1 vector of zeros */
        Math::VectorNd ndof0_vec;    /**< ndof x 1 vector of zeros */
        Math::VectorNd nbodies0_vec; /**< mBodies.size() x 1 vector of zeros */

        Math::MatrixNd tmp_ndof_mat; /**< ndof x ndof matrix for intermediate calculations */
        Math::VectorNd tmp_ndof_vec; /**< ndof x 1 vector for intermediate calculations */

        ////////////////////////////////////
        // Bodies

        /** @brief Transformation from the parent body to the current body
         * \f[
         *	X_{\lambda(i)} = {}^{i} X_{\lambda(i)}
         * \f]
         */
        std::vector<Math::SpatialTransform>X_lambda;

        /// @brief All bodies that are attached to a body via a fixed joint.
        std::vector<FixedBody>mFixedBodies;

        /** @brief Value that is used to discriminate between fixed and movable
         * bodies.
         *
         * Bodies with id 1 .. (fixed_body_discriminator - 1) are moving bodies
         * while bodies with id fixed_body_discriminator .. max (unsigned int)
         * are fixed to a moving body. The value of max(unsigned int) is
         * determined via std::numeric_limits<unsigned int>::max() and the
         * default value of fixed_body_discriminator is max (unsigned int) / 2.
         *
         * On normal systems max (unsigned int) is 4294967294 which means there
         * could be a total of 2147483646 movable and / or fixed bodies.
         */
        unsigned int fixed_body_discriminator;

        /** @brief All bodies 0 ... N_B, including the base
         *
         * mBodies[0] - base body <br>
         * mBodies[1] - 1st moveable body <br>
         * ... <br>
         * mBodies[N_B] - N_Bth moveable body <br>
         */
        std::vector<Body>mBodies;

        /// @brief Human readable names for the bodies
        std::map<std::string, unsigned int>mBodyNameMap;

        std::map<std::string, ReferenceFrame *>bodyFrameMap;

        /** @brief Connects a given body to the model
         *
         * When adding a body there are basically informations required:
         * - what kind of body will be added?
         * - where is the new body to be added?
         * - by what kind of joint should the body be added?
         *
         * The first information "what kind of body will be added" is contained
         * in the Body class that is given as a parameter.
         *
         * The question "where is the new body to be added?" is split up in two
         * parts: first the parent (or successor) body to which it is added and
         * second the transformation to the origin of the joint that connects the
         * two bodies. With these two informations one specifies the relative
         * positions of the bodies when the joint is in neutral position.gk
         *
         * The last question "by what kind of joint should the body be added?" is
         * again simply contained in the Joint class.
         *
         * \param parent_id   id of the parent body
         * \param joint_frame the transformation from the parent frame to the origin
         *                    of the joint frame (represents X_T in RBDA)
         * \param joint       specification for the joint that describes the
         *                    connection
         * \param body        specification of the body itself
         *
         * \param body_name   human readable name for the body (can be used to
         *                    retrieve its id with GetBodyId())
         *
         * \returns id of the added body
         */
        unsigned int addBody(const unsigned int            parent_id,
                             const Math::SpatialTransform& joint_frame,
                             const Joint                 & joint,
                             const Body                  & body,
                             std::string                   body_name = "");

        unsigned int addBodySphericalJoint(const unsigned int            parent_id,
                                           const Math::SpatialTransform& joint_frame,
                                           const Joint                 & joint,
                                           const Body                  & body,
                                           std::string                   body_name = "");

        /** @brief Adds a Body to the model such that the previously added Body
         * is the Parent.
         *
         * This function is basically the same as Model::addBody() however the
         * most recently added body (or body 0) is taken as parent.
         */
        unsigned int appendBody(const Math::SpatialTransform& joint_frame,
                                const Joint                 & joint,
                                const Body                  & body,
                                std::string                   body_name = "");

        unsigned int addBodyCustomJoint(const unsigned int            parent_id,
                                        const Math::SpatialTransform& joint_frame,
                                        CustomJoint                  *custom_joint,
                                        const Body                  & body,
                                        std::string                   body_name = "");

        /** @brief Specifies the dynamical parameters of the first body and
         *  @brief assigns it a 6 DoF joint.
         *
         * The 6 DoF joint is simulated by adding 5 massless bodies at the base
         * which are connected with joints. The body that is specified as a
         * parameter of this function is then added by a 6th joint to the model.
         *
         * The floating base has the following order of degrees of freedom:
         *
         * \li translation X
         * \li translation Y
         * \li translation Z
         * \li rotation Z
         * \li rotation Y
         * \li rotation X
         *
         * To specify a different ordering, it is recommended to create a 6 DoF
         * joint. See \link RobotDynamics::Joint Joint\endlink for more
         * information.
         *
         * \param body Properties of the floating base body.
         *
         *  \returns id of the body with 6 DoF
         */
        unsigned int SetFloatingBaseBody(const Body& body);

        /** @brief Returns the id of a body that was passed to addBody()
         *
         * Bodies can be given a human readable name. This function allows to
         * resolve its name to the numeric id.
         *
         * \note Instead of querying this function repeatedly, it might be
         * advisable to query it once and reuse the returned id.
         *
         * \returns the id of the body or @code{.cpp} std::numeric_limits<unsigned int>::max() @endcode if the id was not found.
         */
        unsigned int GetBodyId(const char *body_name) const
        {
            if (mBodyNameMap.count(body_name) == 0)
            {
                return std::numeric_limits<unsigned int>::max();
            }

            return mBodyNameMap.find(body_name)->second;
        }

        ReferenceFrame* getBodyFrame(const std::string& bodyName) const
        {
            if (bodyFrameMap.count(bodyName) == 0)
            {
                return nullptr;
            }

            return bodyFrameMap.find(bodyName)->second;
        }

        /** @brief Returns the name of a body for a given body id */
        std::string GetBodyName(unsigned int body_id) const
        {
            std::map<std::string, unsigned int>::const_iterator iter = mBodyNameMap.begin();

            while (iter != mBodyNameMap.end())
            {
                if (iter->second == body_id)
                {
                    return iter->first;
                }

                iter++;
            }

            return "";
        }

        /** @brief Checks whether the body is rigidly attached to another body.
         */
        bool IsFixedBodyId(unsigned int body_id) const
        {
            if ((body_id >= fixed_body_discriminator) && (body_id < std::numeric_limits<unsigned int>::max()) &&
                (body_id - fixed_body_discriminator < mFixedBodies.size()))
            {
                return true;
            }
            return false;
        }

        bool IsBodyId(unsigned int id) const
        {
            if ((id > 0) && (id < mBodies.size()))
            {
                return true;
            }

            if ((id >= fixed_body_discriminator) && (id < std::numeric_limits<unsigned int>::max()))
            {
                if (id - fixed_body_discriminator < mFixedBodies.size())
                {
                    return true;
                }
            }
            return false;
        }

        /** @brief Determines id the actual parent body.
         *
         * @note When adding bodies using joints with multiple degrees of
         * freedom, additional virtual bodies are added for each degree of
         * freedom. This function returns the id of the actual
         * non-virtual parent body.
         */
        unsigned int GetParentBodyId(unsigned int id)
        {
            if (id >= fixed_body_discriminator)
            {
                return mFixedBodies[id - fixed_body_discriminator].mMovableParent;
            }

            unsigned int parent_id = lambda[id];

            while (mBodies[parent_id].mIsVirtual)
            {
                parent_id = lambda[parent_id];
            }

            return parent_id;
        }

        /**
         * @brief Returns the joint frame transformtion, i.e. the second argument to
           Model::addBody().
         */
        Math::SpatialTransform GetJointFrame(unsigned int id)
        {
            if (id >= fixed_body_discriminator)
            {
                return mFixedBodies[id - fixed_body_discriminator].mParentTransform;
            }

            unsigned int child_id  = id;
            unsigned int parent_id = lambda[id];

            if (mBodies[parent_id].mIsVirtual)
            {
                while (mBodies[parent_id].mIsVirtual)
                {
                    child_id  = parent_id;
                    parent_id = lambda[child_id];
                }
                return X_T[child_id];
            }
            else
            {
                return X_T[id];
            }
        }

        /**
         * @brief Sets the joint frame transformtion, i.e. the second argument to Model::addBody().
         */
        void SetJointFrame(unsigned int id, const Math::SpatialTransform& transform)
        {
            if (id >= fixed_body_discriminator)
            {
                throw RdlException("Error: setting of parent transform not supported for fixed bodies!");
            }

            unsigned int child_id  = id;
            unsigned int parent_id = lambda[id];

            if (mBodies[parent_id].mIsVirtual)
            {
                while (mBodies[parent_id].mIsVirtual)
                {
                    child_id  = parent_id;
                    parent_id = lambda[child_id];
                }
                X_T[child_id] = transform;
            }
            else if (id > 0)
            {
                X_T[id] = transform;
            }
        }

        /**
         * @brief Gets the quaternion for body i (only valid if body i is connected by
         * a JointTypeSpherical joint)
         *
         * See \ref joint_singularities for details.
         */
        Math::Quaternion GetQuaternion(unsigned int i, const Math::VectorNd& Q) const
        {
            assert(mJoints[i].mJointType == JointTypeSpherical);
            unsigned int q_index = mJoints[i].q_index;
            return Math::Quaternion(Q[q_index], Q[q_index + 1], Q[q_index + 2], Q[multdof3_w_index[i]]);
        }

        /**
         * @brief Sets the quaternion for body i (only valid if body i is connected by
         * a JointTypeSpherical joint)
         *
         * See \ref joint_singularities for details.
         */
        void SetQuaternion(unsigned int i, const Math::Quaternion& quat, Math::VectorNd& Q) const
        {
            assert(mJoints[i].mJointType == JointTypeSpherical);
            unsigned int q_index = mJoints[i].q_index;

            Q[q_index]             = quat[0];
            Q[q_index + 1]         = quat[1];
            Q[q_index + 2]         = quat[2];
            Q[multdof3_w_index[i]] = quat[3];
        }

        unsigned int getCommonMovableParentId(unsigned int id_1, unsigned int id_2) const
        {
            if(IsFixedBodyId(id_1))
            {
                id_1 = mFixedBodies[id_1-fixed_body_discriminator].mMovableParent;
            }

            if(IsFixedBodyId(id_2))
            {
                id_2 = mFixedBodies[id_2-fixed_body_discriminator].mMovableParent;
            }

            if(id_1==id_2)
            {
                return id_1;
            }

            if(id_1==0 || id_2==0)
            {
                return 0;
            }

            unsigned int chain_1_size = lambda_chain[id_1].size();
            unsigned int chain_2_size = lambda_chain[id_2].size();

            if(chain_1_size<=chain_2_size)
            {
                // Can start at 1 bc every list should start with world body id which is zero
                for(unsigned int i = 1; i<chain_1_size; i++)
                {
                    if(lambda_chain[id_1][i]!=lambda_chain[id_2][i])
                    {
                        return lambda_chain[id_1][i-1];
                    }
                }

                return lambda_chain[id_1][chain_1_size-1];
            }
            else
            {
                // Can start at 1 bc every list should start with world body id which is zero
                for(unsigned int i = 1; i<chain_2_size; i++)
                {
                    if(lambda_chain[id_1][i]!=lambda_chain[id_2][i])
                    {
                        return lambda_chain[id_2][i-1];
                    }
                }

                return lambda_chain[id_2][chain_2_size-1];
            }
        }
    };

    /** @} */
}

/* _MODEL_H */
#endif // ifndef RDL_MODEL_H
