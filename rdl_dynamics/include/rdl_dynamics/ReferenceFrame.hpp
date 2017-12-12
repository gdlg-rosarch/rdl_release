/*
 * RDL - Robot Dynamics Library
 * Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef __RDL_REFERENCE_FRAME_HPP__
#define __RDL_REFERENCE_FRAME_HPP__

/**
 * @file ReferenceFrame.hpp
 */

#include <memory>
#include "rdl_dynamics/SpatialAlgebraOperators.h"
#include <string>
#include <vector>
#include <climits>
#include "rdl_dynamics/FrameExceptions.hpp"

namespace RobotDynamics
{
    /**
     * \page reference_frame_page ReferenceFrame
     * Detailed information on reference frames can be found in the @ref reference_frame "Reference Frame Module"
     *
     * @defgroup reference_frame Reference Frame
     * @{
     *
     * The ReferenceFrame object is the foundation of the way kinematics are handled. Each
     * time a body is added via RobotDynamics::Model::addBody(), a pointer to a ReferenceFrame
     * is created for that body and added to either RobotDynamics::Model::bodyFrames or RobotDynamics::Model::fixedBodyFrames
     * depending on the type of body that is added. Additionally a reference frame is placed on the center of mass of each
     * body and is stored in RobotDynamics::Model::bodyCenteredFrames.
     *
     * These reference frames can then be used to explicitly express geometric entities(SpatialMotion, SpatialMomentum,
     **SpatialForce, etc)
     * in a reference frame. To query the name of a reference frame a geometric object is expressed in, you may call the
     **ReferenceFrame::getName()
     * method on the objects ReferenceFrame. Furthermore, to transform a framed geometric entity(anything that inheritcs from
     **the
     * RobotDynamics::FrameObject) into a different reference frame, you may simply call the
     **RobotDynamics::FrameObject::changeFrame
     * method and supply it a pointer to the frame you want the geometric object to be transformed to.
     *
     * \note If you create your own ReferenceFrame outside of those stored in the RobotDynamics::Model, you are required to call
     * RobotDynamics::ReferenceFrame::update every time the created frames' ReferenceFrame::transformFromParent changes.
     * If it is a static frame that doesn't have a joint between it and it's parent frame, then that transform doesn't change,
     **so
     * you don't have to call RobotDynamics::frames::ReferenceFrame::update after construction.
     */

    /**
     * @class ReferenceFrame
     * @brief ReferenceFrame object used to tell what frame objects are expressed in. Every ReferenceFrame has a pointer to its
     **parent
     * ReferenceFrame. This parent frame is NOT allowed to be nullptr. The ONLY ReferenceFrame that is allowed to have
     **parentFrame=nullptr is the world frame. There is
     * only one world frame and it can be accessed by the static method ReferenceFrame::getWorldFrame() which will return a
     **shared_ptr to the world frame.
     * This class and its implementation are an adaptation of ReferenceFrame.java by <a href="http://robots.ihmc.us/">Jerry
     **Pratt and the IHMC Robotics Group</a>.
     */
    class RDL_DLLAPI ReferenceFrame
    {
public:

        /**
         * @brief Copy constructor
         * @param referenceFrameToCopy
         */
        ReferenceFrame(
            const ReferenceFrame& referenceFrameToCopy) : framesStartingWithRootEndingWithThis(referenceFrameToCopy.
                                                                                               framesStartingWithRootEndingWithThis),
                                                          frameName(referenceFrameToCopy.frameName),
                                                          parentFrame(referenceFrameToCopy.parentFrame),
                                                          transformFromParent(referenceFrameToCopy.transformFromParent),
                                                          transformToRoot(referenceFrameToCopy.transformToRoot),
                                                          inverseTransformToRoot(referenceFrameToCopy.inverseTransformToRoot),
                                                          isWorldFrame(referenceFrameToCopy.isWorldFrame),
                                                          isBodyFrame(referenceFrameToCopy.isBodyFrame),
                                                          movableBodyId(referenceFrameToCopy.movableBodyId)
        {}

        /**
         * @brief Constructor
         * @param frameName Name of frame
         * @param parentFrame Pointer to this frames parent ReferenceFrame
         * @param transformFromParent Spatial transform from the parent frames perspective to this frame
         * @param isBodyFrame Boolean true if this frame is a body frame(i.e. will be stored in Model::bodyFrames), false
         **otherwise
         * @param movableBodyId The ID of the movable body this frame is attached to. For frames attached to fixed bodies, this
         **should be FixedBody::mMovableParent
         */
        ReferenceFrame(const std::string                          & frameName,
                       ReferenceFrame                              *parentFrame,
                       const RobotDynamics::Math::SpatialTransform& transformFromParent,
                       bool                                         isBodyFrame,
                       unsigned int                                 movableBodyId) :
            frameName(frameName),
            parentFrame(parentFrame),
            transformFromParent(transformFromParent),
            isWorldFrame(false),
            isBodyFrame(isBodyFrame),
            movableBodyId(movableBodyId)
        {
            if (parentFrame == nullptr)
            {
                throw ReferenceFrameException(
                          "You are not allowed to create a frame with parentFrame=nullptr. Only a root frame and the world frame may have parentFrame=nullptr");
            }

            framesStartingWithRootEndingWithThis = constructVectorOfFramesStartingWithRootEndingWithThis(this);

            update();
        }

        /**
         * @brief Empty constructor. All contained ptrs will be initialize to nullptr
         */
        ReferenceFrame()
        {}

        /**
         * @brief Destructor
         */
        ~ReferenceFrame()
        {}

        /**
         * @brief Recalculates this frames ReferenceFrame::transformToRoot and ReferenceFrame::inverseTransformToRoot which are
         **used
         * by FrameObject::changeFrame to change the ReferenceFrame FrameObjects are expressed in. If
         * you create a ReferenceFrame you MUST call update every tic. If there is a joint between a frame and its' parent, you
         **need
         * to call ReferenceFrame::setTransformFromParent before calling update in order to update the framse
         **ReferenceFrame::transformFromParent.
         */
        void        update();

        /**
         * @brief Get the spatial transform from this frame to desiredFrame and store it in transformToPack
         * @param transformToPack Resulting transform to the desired frame will be stored here
         * @param desiredFrame The resulting transform will transform vectors into desiredFrame
         */
        inline void getTransformToDesiredFrame(RobotDynamics::Math::SpatialTransform& transformToPack,
                                               ReferenceFrame                        *desiredFrame)
        {
            transformToPack = getTransformToDesiredFrame(desiredFrame);
        }

        /**
         * @brief Get the spatial transform from this frame to desiredFrame and store it in transformToPack
         * @param desiredFrame The resulting transform will transform vectors into desiredFrame
         * @return Spatial transform that will transform vectors into the ReferenceFrame desiredFrame
         */
        RobotDynamics::Math::SpatialTransform getTransformToDesiredFrame(ReferenceFrame *desiredFrame);

        /**
         * @brief Check if two frames have the same roots
         * @throws ReferenceFrameException if the frames do not have the same roots
         * @param frame
         */
        void                                  verifyFramesHaveSameRoot(ReferenceFrame *frame);

        /**
         * @brief Set a frames ReferenceFrame::transformToParent. For frames connected by a joint, this needs to be
         * updated every tic BEFORE calling the ReferenceFrame::update method
         * @param transformFromParent The new ReferenceFrame::transformFromParent
         */
        inline void                           setTransformFromParent(
            const RobotDynamics::Math::SpatialTransform& transformFromParent)
        {
            this->transformFromParent = transformFromParent;
        }

        /**
         * @brief Check if the argument ReferenceFrame equals this
         * @param referenceFrame
         */
        void                                         checkReferenceFramesMatch(ReferenceFrame *referenceFrame) const;

        /**
         * @brief Get this frames ReferenceFrame::transformToRoot
         * @return ReferenceFrame::transformToRoot
         */
        inline RobotDynamics::Math::SpatialTransform getTransformToRoot()
        {
            return this->transformToRoot;
        }

        /**
         * @brief Get this frames ReferenceFrame::inverseTransformToRoot
         * @return ReferenceFrame::inverseTransformToRoot
         */
        inline RobotDynamics::Math::SpatialTransform getInverseTransformToRoot()
        {
            return this->inverseTransformToRoot;
        }

        /**
         * @brief Set this frames ReferenceFrame::transformToRoot. Make absolutely sure you know what you are doing
         * before manually setting this. A frames ReferenceFrame::transformToRoot is computed automatically every time
         * ReferenceFrame::update is called, so if you set it manually it will be overridden if you then call
         **ReferenceFrame::update
         * @param transformToRoot
         */
        inline void setTransformToRoot(const RobotDynamics::Math::SpatialTransform& transformToRoot)
        {
            this->transformToRoot = transformToRoot;
        }

        /**
         * @brief Get a pointer to this frames root frame
         * @return Pointer to this frames root ReferenceFrame
         */
        inline ReferenceFrame* getRootFrame()
        {
            return this->framesStartingWithRootEndingWithThis[0];
        }

        /**
         * @brief Get a vector containing all frames in the chain from this frames root to this frame
         * @return A vector of framse with element 0 correspordind to this frames root, and the final element coresponding
         * to this frame
         */
        inline std::vector<ReferenceFrame *>getFramesStartingWithRootEndingWithThis()
        {
            return framesStartingWithRootEndingWithThis;
        }

        /**
         * @brief get a pointer to this frames parent
         * @return ReferenceFrame::parentFrame
         */
        inline ReferenceFrame* getParentFrame()
        {
            return this->parentFrame;
        }

        /**
         * @brief Get the frame name
         * @return ReferenceFrame::framenName
         */
        inline std::string getName() const
        {
            return this->frameName;
        }

        /**
         * @brief Creates a root frame with ReferenceFrame::parentFrame=nullptr
         * @param frameName
         * @return pointer to the created root frame
         */
        static std::shared_ptr<ReferenceFrame>createARootFrame(const std::string& frameName)
        {
            return std::shared_ptr<ReferenceFrame>(new ReferenceFrame(frameName, false, 0, true));
        }

        /**
         * @brief Get a pointer to the world frame
         * @return Pointer to world frame
         */
        static std::shared_ptr<ReferenceFrame>getWorldFrame()
        {
            return worldFrame;
        }

        /**
         * @brief Get spatial transform from parent to this frame
         * @return SpatialTransform from parent to this frame
         */
        inline RobotDynamics::Math::SpatialTransform getTransformFromParent()
        {
            return this->transformFromParent;
        }

        /**
         * @brief Get spatial transform this frame to its parent
         * @return SpatialTransform from this frame to parent
         */
        inline RobotDynamics::Math::SpatialTransform getTransformToParent()
        {
            return this->transformFromParent.inverse();
        }

        /**
         * @brief Get a boolean telling if this frame is the world frame
         * @return Boolean true if this frame is world frame, false otherwise
         */
        inline bool getIsWorldFrame() const
        {
            return this->isWorldFrame;
        }

        /**
         * @brief Get the ID of the movable body this frame is attached to
         * @return Unsigned int corresponding to the movable body this frame is attached to
         */
        inline unsigned int getMovableBodyId() const
        {
            return this->movableBodyId;
        }

        /**
         * @brief Get boolean telling if this frame is a body frame or not. If it is a body frame, A pointer
         * to this frame would be stored in Model::bodyFrames vector
         * @return Boolean true if this frame is a body frame, false otherwise
         */
        inline bool getIsBodyFrame() const
        {
            return this->isBodyFrame;
        }

protected:

        /**
         * Constructor that creates a top level frame with parent=nullptr, and transforms to root being identity transform.
         * @param frameName
         * @param isWorldFrame True of creating the world frame, false if it's just a root frame
         * @param movableBodyId
         * @param isBodyFrame
         */
        ReferenceFrame(const std::string& frameName, bool isWorldFrame, unsigned int movableBodyId, bool isBodyFrame)
        {
            this->frameName     = frameName;
            this->isWorldFrame  = isWorldFrame;
            this->movableBodyId = movableBodyId;
            this->parentFrame   = nullptr;

            this->isBodyFrame                    = isBodyFrame;
            framesStartingWithRootEndingWithThis = constructVectorOfFramesStartingWithRootEndingWithThis(this);

            update();
        }

        /**
         * @brief Helper method to create a world frame
         * @param frameName
         * @return Pointer to the created world frame
         */
        static std::shared_ptr<ReferenceFrame>createAWorldFrame(const std::string& frameName)
        {
            return std::shared_ptr<ReferenceFrame>(new ReferenceFrame(frameName, true, 0, true));
        }

        /**
         * @brief Upon creation of a ReferenceFrame, this method is called to create the vector of frames each ReferenceFrame
         **holds
         * @param thisFrame
         * @return ReferenceFrame::framesStartingWithRootEndingWithThis
         */
        static std::vector<ReferenceFrame *>constructVectorOfFramesStartingWithRootEndingWithThis(
            ReferenceFrame *thisFrame);

        static std::shared_ptr<ReferenceFrame> worldFrame;                  /**< Static world frame pointer */
        std::vector<ReferenceFrame *> framesStartingWithRootEndingWithThis; /**< A vector of frames holding pointers to all
                                                                               frames in the chain from root to this frame */
        std::string frameName;                                              /**< A frames name */
        ReferenceFrame *parentFrame;                                        /**< Pointer to a frames parent frames */
        RobotDynamics::Math::SpatialTransform transformFromParent;          /**< SpatialTransform to a frame from its parent*/
        RobotDynamics::Math::SpatialTransform transformToRoot;              /**< SpatialTransform from a frame to the root frame
                                                                             */
        RobotDynamics::Math::SpatialTransform inverseTransformToRoot;       /**< SpatialTransform to a frame from the root frame
                                                                             */

        bool isWorldFrame;                                                  /**< True if a frame is the world frame, false
                                                                               otherwise */
        bool isBodyFrame;                                                   /**< True if a frame is a body frame, false
                                                                               otherwise. Body frame pointers are stored in
                                                                               Model::bodyFrames */
        unsigned int movableBodyId;                                         /**< The body ID of the movable body a frame is
                                                                               attached to */
    };
}

/**
 * @}
 */

#endif // ifndef __RDL_REFERENCE_FRAME_HPP__
