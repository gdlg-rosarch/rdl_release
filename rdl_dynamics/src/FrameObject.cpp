/*
 * RDL - Robot Dynamics Library
 * Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include "rdl_dynamics/FrameObject.hpp"

namespace RobotDynamics
{
    void FrameObject::changeFrame(ReferenceFrame *desiredFrame)
    {
        if(desiredFrame == nullptr || this->referenceFrame == nullptr)
        {
            throw ReferenceFrameException("Either this reference frame or desired reference frame is nullptr!");
        }

        if(desiredFrame != this->referenceFrame)
        {
            this->referenceFrame->verifyFramesHaveSameRoot(desiredFrame);

            this->getTransformableGeometricObject()->transform(desiredFrame->getInverseTransformToRoot() * this->referenceFrame->getTransformToRoot());
            this->referenceFrame = desiredFrame;
        }
    }
}