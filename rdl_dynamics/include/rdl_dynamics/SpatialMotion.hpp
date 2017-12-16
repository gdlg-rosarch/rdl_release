/*
 * RDL - Robot Dynamics Library
 * Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef __RDL_SPATIAL_MOTION_HPP__
#define __RDL_SPATIAL_MOTION_HPP__

/**
 * @file SpatialMotion.hpp
 */

#include "rdl_dynamics/FrameObject.hpp"
#include "rdl_dynamics/FrameVector.hpp"
#include "rdl_dynamics/TransformableGeometricObject.hpp"

namespace RobotDynamics
{
    namespace Math
    {
        /**
         * @class SpatialMotion
         * @ingroup reference_frame
         * @brief A SpatialMotion vector is a MotionVector with a RobotDynamics::ReferenceFrame it is expressed in. This allows
         * for runtime checks that frame rules are obeyed and makes it easy to change the frame the metion vector is
         * expressed in. As with a SpatialAcceleration, a SpatialMotion vector is the spatial velocity of a
         **SpatialMotion::bodyFrame
         * relative to a SpatialMotion::baseFrame and is expressed in RobotDynamics::FrameObject::referenceFrame
         */
        class SpatialMotion : public MotionVector, public FrameObject
        {
public:

            /**
             * @brief Constructor. SpatialMotion::bodyFrame, SpatialMotion::baseFrame, and
             **RobotDynamics::FrameObject::referenceFrame initialized to nullptr
             */
            SpatialMotion() : MotionVector(), FrameObject(nullptr)
            {
                bodyFrame = nullptr;
                baseFrame = nullptr;
            }

            /**
             * @brief Constructor
             * @param bodyFrame RobotDynamics::ReferenceFrame in motion
             * @param baseFrame RobotDynamics::ReferenceFrame the motion is relative to
             * @param expressedInFrame RobotDynamics::ReferenceFrame the motion vector is expressed in
             * @param wx x-Angular part
             * @param wy y-Angular part
             * @param wz z-Angular part
             * @param vx x-Linear part
             * @param vy y-Linear part
             * @param vz z-Linear part
             */
            SpatialMotion(ReferenceFrame *bodyFrame,
                          ReferenceFrame *baseFrame,
                          ReferenceFrame *expressedInFrame,
                          const double    wx,
                          const double    wy,
                          const double    wz,
                          const double    vx,
                          const double    vy,
                          const double    vz) : MotionVector(wx, wy, wz, vx, vy, vz), FrameObject(expressedInFrame)
            {
                this->bodyFrame = bodyFrame;
                this->baseFrame = baseFrame;
            }

            /**
             * @brief Constructor
             * @param bodyFrame RobotDynamics::ReferenceFrame in motion
             * @param baseFrame RobotDynamics::ReferenceFrame the motion is relative to
             * @param expressedInFrame RobotDynamics::ReferenceFrame the motion vector is expressed in
             * @param w Angular part
             * @param v Linear part
             */
            SpatialMotion(ReferenceFrame *bodyFrame,
                          ReferenceFrame *baseFrame,
                          ReferenceFrame *expressedInFrame,
                          const Vector3d& w,
                          const Vector3d  v) : MotionVector(w.x(), w.y(), w.z(), v.x(), v.y(),
                                                                           v
                                                                           .z()), FrameObject(expressedInFrame)
            {
                this->bodyFrame = bodyFrame;
                this->baseFrame = baseFrame;
            }

            /**
             * @brief Constructor
             * @param bodyFrame RobotDynamics::ReferenceFrame in motion
             * @param baseFrame RobotDynamics::ReferenceFrame the motion is relative to
             * @param expressedInFrame RobotDynamics::ReferenceFrame the motion vector is expressed in
             * @param v
             */
            SpatialMotion(ReferenceFrame      *bodyFrame,
                          ReferenceFrame      *baseFrame,
                          ReferenceFrame      *expressedInFrame,
                          const SpatialVector& v) : MotionVector(v), FrameObject(expressedInFrame)
            {
                this->bodyFrame = bodyFrame;
                this->baseFrame = baseFrame;
            }

            /**
             * @brief Copy constructor
             * @param spatialMotion
             */
            SpatialMotion(
                const SpatialMotion& spatialMotion) : MotionVector(spatialMotion), FrameObject(spatialMotion.referenceFrame)
            {
                this->bodyFrame = spatialMotion.bodyFrame;
                this->baseFrame = spatialMotion.baseFrame;
            }

            /**
             * @brief Get linear part of spatial motion as a frame vector
             * @return FrameVector consisting of the reference frame and the linear portion
             */
            FrameVector getFramedLinearPart() const
            {
                return FrameVector(this->referenceFrame, this->getLinearPart());
            }

            /**
             * @brief Get angular part of spatial motion as a frame vector
             * @return FrameVector consisting of the reference frame and the angular portion
             */
            FrameVector getFramedAngularPart() const
            {
                return FrameVector(this->referenceFrame, this->getAngularPart());
            }

            Math::TransformableGeometricObject* getTransformableGeometricObject()
            {
                return this;
            }

           /**
            * @brief Copy and change frame
            * @param referenceFrame
            * @return Copied spatial transform with frame changed
            */
            SpatialMotion changeFrameAndCopy(std::shared_ptr<ReferenceFrame> referenceFrame) const
            {
                return changeFrameAndCopy(referenceFrame.get());
            }

            /**
             * @brief Copy and change frame
             * @param referenceFrame
             * @return Copied spatial transform with frame changed
             */
            SpatialMotion changeFrameAndCopy(ReferenceFrame* referenceFrame) const
            {
                SpatialMotion ret = *this;
                ret.changeFrame(referenceFrame);
                return ret;
            }

            /**
             * @brief Copy this into spatial motion and change frame
             * @param referenceFrame
             * @param spatialMotion Storage for copied result
             */
            void changeFrameAndCopy(ReferenceFrame* referenceFrame, SpatialMotion &spatialMotion) const
            {
                spatialMotion = *this;
                spatialMotion.changeFrame(referenceFrame);
            }

           /**
            * @brief Copy this into spatial motion and change frame
            * @param referenceFrame
            * @param spatialMotion Storage for copied result
            */
            void changeFrameAndCopy(std::shared_ptr<ReferenceFrame> referenceFrame, SpatialMotion &spatialMotion) const
            {
                this->changeFrameAndCopy(referenceFrame.get(),spatialMotion);
            }

            /**
             * @brief Get a SpatialMotions SpatialMotion::bodyFrame
             * @return The body frame, i.e. the principal moving frame
             */
            inline ReferenceFrame* getBodyFrame() const
            {
                return bodyFrame;
            }

            /**
             * @brief Get a SpatialMotions SpatialMotion::baseFrame
             * @return The base frame, i.e. the frame the SpatialMotion is relative to
             */
            inline ReferenceFrame* getBaseFrame() const
            {
                return baseFrame;
            }

            /**
             * @brief Get a copy of this SpatialMotion as type MotionVector
             * @return MotionVector that is a copy of a SpatialMotion
             */
            inline MotionVector toMotionVector() const
            {
                return *this;
            }

            void setIncludingFrame(ReferenceFrame *referenceFrame, const SpatialVector& v)
            {
                this->set(v);
                this->referenceFrame = referenceFrame;
            }

            void setIncludingFrame(ReferenceFrame *referenceFrame,
                                   double          wx,
                                   double          wy,
                                   double          wz,
                                   double          vx,
                                   double          vy,
                                   double          vz)
            {
                this->SpatialVector::set(wx, wy, wz, vx, vy, vz);
                this->referenceFrame = referenceFrame;
            }

            /**
             * @brief Overloaded += operator. Frame checks are performed.
             * @param v
             */
            void operator+=(const SpatialMotion& v);

            /**
             * @brief Overloaded -= operator. Frame checks are performed
             * @param v
             */
            void operator-=(const SpatialMotion& v);

            /**
             * @brief This is an operator for performing what is referred to in Featherstone as the spatial vector
             **cross(\f$\times\f$) operator times a MotionVector. In
             * VDuindam the spatial vector cross operator is referred to as the adjoint of a twist, denoted \f$ad_T\f$. This
             **method performs,
             * \f[
             * m =(v\times) m = ad_{v}m =
             *  \begin{bmatrix}
             * \omega\times & \mathbf{0} \\
             * v\times & \omega\times
             *  \end{bmatrix} m
             * \f]
             * The 3d vector \f$\times\f$ operator is equivalent to the \f$\sim\f$ operator. See Math::toTildeForm.
             * @param v
             */
            void operator%=(const SpatialMotion& v);

            /**
             * @brief Sets the body frame of a spatial motion
             *
             * @param referenceFrame
             */
            void setBodyFrame(ReferenceFrame *referenceFrame)
            {
                this->bodyFrame = referenceFrame;
            }

            /**
             * @brief Sets the base frame of a spatial motion
             * @param referenceFrame
             */
            void setBaseFrame(ReferenceFrame *referenceFrame)
            {
                this->baseFrame = referenceFrame;
            }

protected:

            ReferenceFrame *bodyFrame, *baseFrame;
        };

        EIGEN_STRONG_INLINE SpatialMotion operator+(SpatialMotion v1, const SpatialMotion& v2)
        {
            v1 += v2;
            return v1;
        }

        EIGEN_STRONG_INLINE SpatialMotion operator-(SpatialMotion v1, const SpatialMotion& v2)
        {
            v1 -= v2;
            return v1;
        }

        EIGEN_STRONG_INLINE SpatialMotion operator%(SpatialMotion v1, const SpatialMotion& v2)
        {
            v1.SpatialMotion::operator%=(v2);
            return v1;
        }
    }
}

#endif //__RDL_SPATIAL_MOTION_HPP__
