/*
 * RDL - Robot Dynamics Library
 * Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef __RDL_SPATIAL_RIGID_BODY_INERTIA_HPP__
#define __RDL_SPATIAL_RIGID_BODY_INERTIA_HPP__

/**
 * @file SpatialRigidBodyInertia.hpp
 */

#include "rdl_dynamics/FrameObject.hpp"
#include "rdl_dynamics/RigidBodyInertia.hpp"

namespace RobotDynamics
{
    namespace Math
    {
        /**
         * @class SpatialInertia
         * @ingroup reference_frame
         * @brief A Math::SpatialInertia is a RigidBodyInertia explicitly expressed in a RobotDynamics::ReferenceFrame. The
         **frame
         * a Math::SpatialInertia is expressed in can be changed by calling RobotDynamics::FrameObject::changeFrame.
         */
        class SpatialInertia : public RigidBodyInertia, public FrameObject
        {
public:

            /**
             * @brief Empty constructor. Initializes FrameObject::referenceFrame to nullptr
             */
            SpatialInertia() : RigidBodyInertia(), FrameObject(nullptr)
            {}

            /**
             * @brief Constructor
             * @param referenceFrame
             */
            SpatialInertia(ReferenceFrame *referenceFrame) : RigidBodyInertia(), FrameObject(referenceFrame)
            {}

            /**
             * @brief Constructor
             * @param referenceFrame
             * @param mass
             * @param com_mass Mass scaled center of mass vector
             * @param inertia 3x3 inertia tensor
             */
            SpatialInertia(ReferenceFrame *referenceFrame, double mass, const Vector3d& com_mass,
                           const Matrix3d& inertia) : RigidBodyInertia(mass, com_mass, inertia), FrameObject(referenceFrame)
            {}

            /**
             * @brief Constructor
             * @param referenceFrame
             * @param m
             * @param h Mass scaled center of mass
             * @param Ixx
             * @param Iyx
             * @param Iyy
             * @param Izx
             * @param Izy
             * @param Izz
             */
            SpatialInertia(ReferenceFrame *referenceFrame,
                           double          m,
                           const Vector3d& h,
                           const double    Ixx,
                           const double    Iyx,
                           const double    Iyy,
                           const double    Izx,
                           const double    Izy,
                           const double    Izz) : RigidBodyInertia(m, h, Ixx, Iyx, Iyy, Izx, Izy,
                                                                                   Izz),
                                                                  FrameObject(referenceFrame)
            {}

            /**
             * @brief Constructor
             * @param referenceFrame
             * @param inertia
             */
            SpatialInertia(ReferenceFrame         *referenceFrame,
                           const RigidBodyInertia& inertia) : RigidBodyInertia(inertia), FrameObject(referenceFrame)
            {}

            /**
             * @brief Copy constructor
             * @param inertia
             */
            SpatialInertia(
                const SpatialInertia& inertia) : RigidBodyInertia(inertia), FrameObject(inertia.referenceFrame)
            {}

            TransformableGeometricObject* getTransformableGeometricObject()
            {
                return this;
            }

            /**
             * @brief Get a copy of a Math::SpatialInertia as type Math::RigidBodyInertia
             * @return A copy as type Math::RigidBodyInertia
             */
            RigidBodyInertia toRigidBodyInertia() const
            {
                return RigidBodyInertia(this->m, this->h, this->Ixx, this->Iyx, this->Iyy, this->Izx, this->Izy, this->Izz);
            }

            /**
             * @brief Overloaded += operator. Performs frame checks
             * @param spatialInertia
             */
            void      operator+=(const SpatialInertia& spatialInertia)
            {
                checkReferenceFramesMatch(&spatialInertia);
                this->RigidBodyInertia::operator+=(spatialInertia);
            }
        };

        /**
         * @brief Adds two Math::SpatialInertia together. Performs frame checks.
         * @param inertia_a
         * @param inertia_b
         * @return The addition of the two arguments
         */
        inline SpatialInertia operator+(SpatialInertia inertia_a, const SpatialInertia& inertia_b)
        {
            inertia_a += inertia_b;
            return inertia_a;
        }
    }
}

#endif //__SPATIAL_RIGID_BODY_INERTIA_HPP__
