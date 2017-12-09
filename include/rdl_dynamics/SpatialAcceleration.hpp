/*
 * RDL - Robot Dynamics Library
 * Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef __RDL_SPATIAL_ACCELERATION_HPP__
#define __RDL_SPATIAL_ACCELERATION_HPP__

#include "rdl_dynamics/SpatialMotion.hpp"

/**
 * @file SpatialAcceleration.hpp
 */

namespace RobotDynamics
{
    namespace Math
    {
        /**
         * @class SpatialAcceleration
         * @ingroup reference_frame
         * @brief SpatialAcceleration. For clarity, the ReferenceFrames are stated as follows. A spatial acceleration is
         * the acceleration of the SpatialAcceleration::bodyFrame with respect to the SpatialAcceleration::baseFrame
         * expressed in the SpatialAcceleration::expressedInFrame
         */
        class SpatialAcceleration : public SpatialMotion
        {
public:

            /**
             * @brief Constructor. ReferenceFrame is initialized to nullptr and vector elements to 0.
             */
            SpatialAcceleration() : SpatialMotion()
            {}

            /**
             * @brief Constructor
             * @param bodyFrame The acceleration is of this frame
             * @param baseFrame The acceleration is relative to this frame
             * @param expressedInFrame The acceleration is expressed in this frame
             * @param wx Angular x-coordinate
             * @param wy Angular y-coordinate
             * @param wz Angular z-coordinate
             * @param vx Linear x-coordinate
             * @param vy Linear y-coordinate
             * @param vz Linear z-coordinate
             */
            SpatialAcceleration(ReferenceFrame *bodyFrame,
                                ReferenceFrame *baseFrame,
                                ReferenceFrame *expressedInFrame,
                                const double    wx,
                                const double    wy,
                                const double    wz,
                                const double    vx,
                                const double    vy,
                                const double    vz) : SpatialMotion(bodyFrame, baseFrame, expressedInFrame, wx, wy, wz, vx, vy,
                                                                                         vz)
            {}

            /**
             * @brief Constructor
             * @param bodyFrame The acceleration is of this frame
             * @param baseFrame The acceleration is relative to this frame
             * @param expressedInFrame The acceleration is expressed in this frame
             * @param w Vector containig the angular component of the acceleration
             * @param v Vector containing the linear component of the acceleration
             */
            SpatialAcceleration(ReferenceFrame *bodyFrame,
                                ReferenceFrame *baseFrame,
                                ReferenceFrame *expressedInFrame,
                                const Vector3d& w,
                                const Vector3d  v) : SpatialMotion(bodyFrame, baseFrame, expressedInFrame, w.x(), w.y(), w.z(),
                                                                                        v
                                                                                        .x(), v.y(), v.z())
            {}

            /**
             * @brief Constructor
             * @param bodyFrame The acceleration is of this frame
             * @param baseFrame The acceleration is relative to this frame
             * @param expressedInFrame The acceleration is expressed in this frame
             * @param v
             */
            SpatialAcceleration(ReferenceFrame      *bodyFrame,
                                ReferenceFrame      *baseFrame,
                                ReferenceFrame      *expressedInFrame,
                                const SpatialVector& v) : SpatialMotion(bodyFrame, baseFrame, expressedInFrame, v)
            {}

            /**
             * @brief Constructor
             * @param spatialAcceleration
             */
            SpatialAcceleration(
                const SpatialAcceleration& spatialAcceleration) : SpatialMotion(spatialAcceleration.getBodyFrame(),
                                                                                spatialAcceleration.getBaseFrame(),
                                                                                spatialAcceleration.getReferenceFrame(),
                                                                                spatialAcceleration)
            {}

            /**
             * @brief Use this method to change the ReferenceFrame::expressedInFrame of a SpatialAcceleration if there is
             * relative acceleration between the current frame and the desired frame. See V. Duindam 2.8(c) for how to
             * transform a twist(i.e. SpatialMotion), and take the derivative.
             * @param newFrame The new RobotDynamics::ReferenceFrame
             * @param twistOfCurrentFrameWithRespectToNewFrame
             * @param twistOfBodyWrtBaseExpressedInCurrent
             */
            void changeFrameWithRelativeMotion(
                ReferenceFrame *newFrame, SpatialMotion twistOfCurrentFrameWithRespectToNewFrame, const SpatialMotion
                &twistOfBodyWrtBaseExpressedInCurrent);

            /**
             * @brief Add a SpatialAcceleration, \f$ a_sp=a_sp+vector \f$. Frame checks are performed to ensure frame
             * rules are adhered to
             * @param vector
             */
            void operator+=(const SpatialAcceleration& vector);

            /**
             * @brief Subtract a SpatialAcceleration, \f$ a_sp=a_sp-vector \f$. Frame checks are performed to ensure frame
             * rules are adhered to
             * @param vector
             */
            void operator-=(const SpatialAcceleration& vector);
        };

        /**
         * @brief Add two SpatialAccelerations. Frame check will be performed to ensure frame rules are abided by
         * @param v1
         * @param v2
         * @return A new SpatialAcceleration that is the addition of the two arguments.
         */
        inline SpatialAcceleration operator+(SpatialAcceleration v1, const SpatialAcceleration& v2)
        {
            v1 += v2;
            return v1;
        }

        /**
         * @brief Subtract two SpatialAccelerations. Frame check will be performed to ensure frame rules are abided by
         * @param v1
         * @param v2
         * @return A new SpatialAcceleration that is the subtraction of the second argument from the first.
         */
        inline SpatialAcceleration operator-(SpatialAcceleration v1, const SpatialAcceleration& v2)
        {
            v1 -= v2;
            return v1;
        }
    }
}

#endif //__RDL_SPATIAL_ACCELERATION_HPP__
