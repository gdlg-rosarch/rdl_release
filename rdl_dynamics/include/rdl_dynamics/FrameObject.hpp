/*
 * RDL - Robot Dynamics Library
 * Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef __RDL_FRAME_OBJECT_HPP__
#define __RDL_FRAME_OBJECT_HPP__

#include "rdl_dynamics/TransformableGeometricObject.hpp"
#include "rdl_dynamics/ReferenceFrameHolder.hpp"

namespace RobotDynamics
{
    /**
     * @class FrameObject
     * @ingroup reference_frame
     * @brief An interface that objects with a ReferenceFrame extend to inherit the FrameObject::changeFrame method
     */
    class FrameObject : public ReferenceFrameHolder
    {
public:

        EIGEN_STRONG_INLINE FrameObject(ReferenceFrame *referenceFrame)
        {
            this->referenceFrame = referenceFrame;
        }

        /**
         * @brief Destructor
         */
        virtual ~FrameObject()
        {}

        /**
         * @brief Change the ReferenceFrame this FrameObject is expressed in.
         * @param desiredFrame A pointer to the ReferenceFrame this FrameObject is to be transformed to
         */
        virtual void changeFrame(ReferenceFrame *desiredFrame);

        /**
         * @brief Change the ReferenceFrame this FrameObject is expressed in
         * @param desiredFrame A shared pointer to the ReferenceFrame this FrameObject is to be transformed to
         */
        inline void  changeFrame(std::shared_ptr<ReferenceFrame>desiredFrame)
        {
            this->changeFrame(desiredFrame.get());
        }

        /**
         * @brief Get a pointer to the reference frame this FrameObject is expressed in
         * @return Pointer to the ReferenceFrame this FrameObject is expressed in
         */
        inline ReferenceFrame* getReferenceFrame() const
        {
            return referenceFrame;
        }

protected:

        ReferenceFrame *referenceFrame /**< Pointer to a ReferenceFrame*/;

        /**
         * @brief Pure virtual method that FrameObjects are required to implement so the FrameObject::changeFrame method is able
         * to access the TransformableGeometricObject which is required to implement the TransformableGeometricObject::transform
         **method
         * @return
         */
        virtual Math::TransformableGeometricObject* getTransformableGeometricObject() = 0;
    };
}
#endif // ifndef __RDL_FRAME_OBJECT_HPP__
