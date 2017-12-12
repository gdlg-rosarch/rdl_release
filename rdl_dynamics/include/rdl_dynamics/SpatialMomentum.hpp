/*
 * RDL - Robot Dynamics Library
 * Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef __RDL_SPATIAL_MOMENTUM_HP__
#define __RDL_SPATIAL_MOMENTUM_HP__

/**
 * @file SpatialMomentum.hpp
 */

#include "rdl_dynamics/SpatialRigidBodyInertia.hpp"
#include "rdl_dynamics/Momentum.hpp"
#include "rdl_dynamics/SpatialMotion.hpp"

namespace RobotDynamics
{
    namespace Math
    {
        /**
         * @class SpatialMomentum
         * @ingroup reference_frame
         * @brief A SpatialMomentum is a Momentum expressed in a RobotDynamics::ReferenceFrame. The angular portion of the
         * vector is referred to as \f$k\f$ and the linear portion as \f$l\f$
         */
        class SpatialMomentum : public Momentum, public FrameObject
        {
public:

            /**
             * @brief Constructor. RobotDynamics::ReferenceFrame is initialized to nullptr
             */
            SpatialMomentum() : Momentum(), FrameObject(nullptr)
            {}

            /**
             * @brief Costructor
             * @param referenceFrame RobotDynamics::ReferenceFrame the momentum is expressed in
             * @param kx x-Angular
             * @param ky y-Angular
             * @param kz z-Angular
             * @param lx x-Linear
             * @param ly y-Linear
             * @param lz z-Linear
             */
            SpatialMomentum(ReferenceFrame *referenceFrame,
                            const double    kx,
                            const double    ky,
                            const double    kz,
                            const double    lx,
                            const double    ly,
                            const double    lz) : Momentum(kx, ky, kz, lx, ly, lz), FrameObject(referenceFrame)
            {}

            /**
             * @brief Constructor
             * @param referenceFrame RobotDynamics::ReferenceFrame the momentum is expressed in
             * @param k Angular part
             * @param l Linear part
             */
            SpatialMomentum(ReferenceFrame *referenceFrame, const Vector3d& k, const Vector3d l) : Momentum(k.x(), k.y(), k.z(),
                                                                                                            l.x(), l.y(),
                                                                                                            l.z()),
                                                                                                   FrameObject(referenceFrame)
            {}

            /**
             * @brief Copy constructor
             * @param SpatialMomentum
             */
            SpatialMomentum(
                const SpatialMomentum& SpatialMomentum) : Momentum(SpatialMomentum), FrameObject(SpatialMomentum.referenceFrame)
            {}

            /**
             * @brief Constructor
             * @param spatialForce
             */
            SpatialMomentum(
                const SpatialForce& spatialForce) : Momentum(spatialForce), FrameObject(spatialForce.getReferenceFrame())
            {}

            /**
             * @brief Constructor
             * @param referenceFrame RobotDynamics::ReferenceFrame the momentum is expressed in
             * @param forceVector
             */
            SpatialMomentum(ReferenceFrame    *referenceFrame,
                            const ForceVector& forceVector) : Momentum(forceVector), FrameObject(referenceFrame)
            {}

            /**
             * @brief Constructor
             * @param referenceFrame RobotDynamics::ReferenceFrame the momentum is expressed in
             */
            SpatialMomentum(ReferenceFrame *referenceFrame) : Momentum(), FrameObject(referenceFrame)
            {}

            /**
             * @brief Constructor. Momentum is computed from SpatialInertia and SpatialMotion by \f$ p=I*m \f$. Frame checks
             * are performed
             * @param inertia
             * @param vector
             */
            SpatialMomentum(const SpatialInertia& inertia,
                            const SpatialMotion & vector) : Momentum(), FrameObject(inertia.getReferenceFrame())
            {
                set(inertia, vector);
            }

            /**
             * @brief Set a momentums value by computing a momentum from the supplied SpatialInertia and SpatialMotion. Frame
             * checks will be performed.
             * @param inertia
             * @param vector
             */
            void set(const SpatialInertia& inertia, const SpatialMotion& vector)
            {
                inertia.getReferenceFrame()->checkReferenceFramesMatch(vector.getReferenceFrame());
                this->referenceFrame = inertia.getReferenceFrame();
                computeMomentum(inertia, vector);
            }

            /**
             * @brief Set a momentum by multiplying the supplied motion vector by the supplied inertia. Set the stored
             * RobotDynamics::FrameObject::referenceFrame to the supplied ReferenceFrame pointer. No frame checks will be
             **performed
             * @param referenceFrame
             * @param inertia
             * @param vector
             */
            void set(ReferenceFrame *referenceFrame, const RigidBodyInertia& inertia, const MotionVector& vector)
            {
                this->referenceFrame = referenceFrame;
                computeMomentum(inertia, vector);
            }

            /**
             * @brief Constructor
             * @param referenceFrame
             * @param f
             */
            void set(ReferenceFrame *referenceFrame, const ForceVector& f)
            {
                this->referenceFrame = referenceFrame;
                ForceVector::set(f);
            }

            /**
             * @brief Operator that returns the kinetic energy. Given a momentum, \f$m\f$ and motion vector \f$v\f$,
             * the resulting kinetic energy is computed by \f$ E_k = 0.5 m \cdot v \f$. Frame checks will be performed.
             * @return double The kinetic energy
             */
            inline double operator*(const SpatialMotion& vector)
            {
                this->checkReferenceFramesMatch(&vector);
                return    Momentum::operator*(vector);
            }

protected:

            TransformableGeometricObject* getTransformableGeometricObject()
            {
                return this;
            }
        };

        /**
         * Overloaded * operator for computing a SpatialMomentum from a SpatialInertia and SpatialMotion. Frame checks will
         * be performed.
         * @param inertia
         * @param vector
         * @return The resulting SpatialMomentum
         */
        inline SpatialMomentum operator*(const SpatialInertia& inertia, const SpatialMotion& vector)
        {
            inertia.getReferenceFrame()->checkReferenceFramesMatch(vector.getReferenceFrame());
            SpatialMomentum m(inertia, vector);

            return m;
        }
    }
}

#endif //__RDL_SPATIAL_MOMENTUM_HP__
