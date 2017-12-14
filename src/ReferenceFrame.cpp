/*
 * RDL - Robot Dynamics Library
 * Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include "rdl_dynamics/ReferenceFrame.hpp"

namespace RobotDynamics
{
    std::shared_ptr<ReferenceFrame> ReferenceFrame::worldFrame = ReferenceFrame::createAWorldFrame("World");

    std::vector<ReferenceFrame *> ReferenceFrame::constructVectorOfFramesStartingWithRootEndingWithThis(
            ReferenceFrame *thisFrame)
    {
        ReferenceFrame *parentFrame = thisFrame->getParentFrame();

        if(parentFrame == nullptr)
        {
            // referenceFrame is the root frame.
            std::vector<ReferenceFrame *> vector(1);
            vector[0] = thisFrame;

            return vector;
        }

        // Need to add referenceFrame to the chain.
        int nElements = parentFrame->framesStartingWithRootEndingWithThis.size() + 1;
        std::vector<ReferenceFrame *> vector(nElements);

        for(int i = 0; i < (nElements - 1); i++)
        {
            vector[i] = parentFrame->framesStartingWithRootEndingWithThis[i];
        }

        vector[nElements - 1] = thisFrame;

        return vector;
    }

    RobotDynamics::Math::SpatialTransform ReferenceFrame::getTransformToDesiredFrame(ReferenceFrame *desiredFrame)
    {
        verifyFramesHaveSameRoot(desiredFrame);

        return desiredFrame->inverseTransformToRoot * this->transformToRoot;
    }

    void ReferenceFrame::update()
    {
        if(parentFrame == nullptr)
        {
            return;
        }

        inverseTransformToRoot = transformFromParent * parentFrame->inverseTransformToRoot;
        transformToRoot = inverseTransformToRoot.inverse();
    }

    void ReferenceFrame::checkReferenceFramesMatch(ReferenceFrame *referenceFrame) const
    {
        if(referenceFrame == nullptr)
        {
            throw ReferenceFrameException("Reference frame is nullptr!");
        }

        if(referenceFrame != this)
        {
            throw ReferenceFrameException("Reference frames do not match!");
        }
    }

    void ReferenceFrame::verifyFramesHaveSameRoot(ReferenceFrame *frame)
    {
        if(!(frame->getRootFrame() == this->getRootFrame()))
        {
            std::string msg = "Frames " + frame->getName() + " and " + this->getName() + " have mismatched roots!";
            throw ReferenceFrameException(msg);
        }
    }
}