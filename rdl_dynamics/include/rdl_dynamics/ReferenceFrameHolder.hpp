/*
 * RDL - Robot Dynamics Library
 * Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef __RDL_REFERENCE_FRAME_HOLDER_HPP__
#define __RDL_REFERENCE_FRAME_HOLDER_HPP__

#include "rdl_dynamics/ReferenceFrame.hpp"

namespace RobotDynamics
{
    /**
     * @class ReferenceFrameHolder
     * @ingroup reference_frame
     * @brief This class and its
     * implementation are an adaptation of ReferenceFrameHolder.java by <a href="http://robots.ihmc.us/">Jerry Pratt and the
     **IHMC Robotics Group</a>.
     */
    class RDL_DLLAPI ReferenceFrameHolder
    {
public:

        /**
         * @brief Objects inheriting from this class must implement this method. This makes it possible to
         * access any ReferenceFrameHolders referenceFrame to do runtime frame checks on reference frames
         * @return A pointer to a ReferenceFrameHolder's ReferenceFrame
         */
        virtual ReferenceFrame* getReferenceFrame() const = 0;

        /**
         * @brief Check if two ReferenceFrames are the same
         * @param referenceFrame
         */
        void                    checkReferenceFramesMatch(ReferenceFrame *referenceFrame) const
        {
            getReferenceFrame()->checkReferenceFramesMatch(referenceFrame);
        }

        /**
         * @brief Check if two ReferenceFrameHolders hold the same ReferenceFrame
         * @param referenceFrameHolder
         */
        void checkReferenceFramesMatch(const ReferenceFrameHolder *referenceFrameHolder) const
        {
            getReferenceFrame()->checkReferenceFramesMatch(referenceFrameHolder->getReferenceFrame());
        }

        /**
         * @brief Check if two ReferenceFrameHolders hold the same ReferenceFrame
         * @param referenceFrameHolder
         */
        void checkReferenceFramesMatch(ReferenceFrameHolder *referenceFrameHolder) const
        {
            getReferenceFrame()->checkReferenceFramesMatch(referenceFrameHolder->getReferenceFrame());
        }
    };
}

#endif // ifndef __RDL_REFERENCE_FRAME_HOLDER_HPP__
