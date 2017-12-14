/*
 * RDL - Robot Dynamics Library
 * Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef __RDL_DYNAMICS_FRAME_ORIENTATION_HPP__
#define __RDL_DYNAMICS_FRAME_ORIENTATION_HPP__

/**
 * @file FrameOrientation.hpp
 */

#include "rdl_dynamics/FrameObject.hpp"
#include "rdl_dynamics/Quaternion.h"
#include "rdl_dynamics/TransformableGeometricObject.hpp"

namespace RobotDynamics
{
    namespace Math
    {
        /**
         * @class FrameOrientation
         * @ingroup reference_frame
         * @brief A Frame object that represents an orientation(quaternion) relative to a reference frame
         */
        class FrameOrientation : public FrameObject, TransformableGeometricObject
        {
public:

            FrameOrientation() : FrameObject(nullptr), q(0., 0., 0., 1.)
            {}

            FrameOrientation(ReferenceFrame *referenceFrame) : FrameObject(referenceFrame), q(0., 0., 0., 1.)
            {}

            FrameOrientation(ReferenceFrame *referenceFrame, Quaternion quat) : FrameObject(referenceFrame), q(quat)
            {}

            /**
             *
             * @param referenceFrame
             * @param rotation
             *
             * @note assumes rotation is a valid, orthogonal rotation matrix
             */
            FrameOrientation(ReferenceFrame *referenceFrame,
                             const Matrix3d& rotation) : FrameObject(referenceFrame), q(Quaternion::fromMatrix(rotation))
            {}

            Quaternion getOrientation() const
            {
                return q;
            }

            Quaternion* getOrientationPtr()
            {
                return &q;
            }

            void setOrientation(const Matrix3d& rotation)
            {
                this->q = Quaternion::fromMatrix(rotation);
            }

            void setOrientation(const Quaternion& rotation)
            {
                this->q = rotation;
            }

            void setIncludingFrame(const Quaternion& q, ReferenceFrame *referenceFrame)
            {
                this->referenceFrame = referenceFrame;
                setOrientation(q);
            }

            void transform(const RobotDynamics::Math::SpatialTransform& X)
            {
                q *= Quaternion::fromMatrix(X.E.transpose());
            }

            FrameOrientation changeFrameAndCopy(ReferenceFrame* referenceFrame) const
            {
                FrameOrientation ret = *this;
                ret.changeFrame(referenceFrame);
                return ret;
            }

            void changeFrameAndCopy(ReferenceFrame* referenceFrame, FrameOrientation &frameOrientation) const
            {
                frameOrientation = *this;
                frameOrientation.changeFrame(referenceFrame);
            }

            FrameOrientation changeFrameAndCopy(std::shared_ptr<ReferenceFrame> referenceFrame) const
            {
                return changeFrameAndCopy(referenceFrame.get());
            }

            void changeFrameAndCopy(std::shared_ptr<ReferenceFrame> referenceFrame, FrameOrientation &frameOrientation) const
            {
                changeFrameAndCopy(referenceFrame.get(),frameOrientation);
            }

            Math::TransformableGeometricObject* getTransformableGeometricObject()
            {
                return this;
            }

protected:

            Quaternion q;
        };
    }
}

#endif //__RDL_DYNAMICS_FRAME_ORIENTATION_HPP__
