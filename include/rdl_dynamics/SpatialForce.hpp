/*
 * RDL - Robot Dynamics Library
 * Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef __RDL_SPATIAL_FORCE_HPP__
#define __RDL_SPATIAL_FORCE_HPP__

/**
 * @file SpatialForce.hpp
 */

#include "rdl_dynamics/FrameObject.hpp"
#include "rdl_dynamics/TransformableGeometricObject.hpp"

namespace RobotDynamics
{
    namespace Math
    {
        /**
         * @class SpatialForce
         * @ingroup reference_frame
         * @brief A SpatialForce is a spatial vector with the angular part being three moments and the linear part
         * being 3 linear forces
         */
        class SpatialForce : public ForceVector, public FrameObject
        {
public:

            /**
             * @brief Constructor. RobotDynamics::FrameObject::referenceFrame is initialized to nullptr
             */
            SpatialForce() : ForceVector(), FrameObject(nullptr)
            {}

            /**
             * @brief Constructor
             * @param referenceFrame ReferenceFrame this force vector is expressed in
             * @param mx x-Momentum
             * @param my y-Momentum
             * @param mz z-Momentum
             * @param fx x-Force
             * @param fy y-Force
             * @param fz z-Force
             */
            SpatialForce(ReferenceFrame *referenceFrame,
                         const double    mx,
                         const double    my,
                         const double    mz,
                         const double    fx,
                         const double    fy,
                         const double    fz) : ForceVector(mx, my, mz, fx, fy, fz), FrameObject(referenceFrame)
            {}

            /**
             * @brief Constructor
             * @param referenceFrame ReferenceFrame this force vector is expressed in
             * @param m Vector containing the angular component
             * @param f Vector containing the linear component
             */
            SpatialForce(ReferenceFrame *referenceFrame, const Vector3d& m, const Vector3d f) : ForceVector(m.x(), m.y(), m.z(),
                                                                                                            f.x(), f.y(),
                                                                                                            f.z()),
                                                                                                FrameObject(referenceFrame)
            {}

            Math::TransformableGeometricObject* getTransformableGeometricObject()
            {
                return this;
            }

          /**
          * @brief Copy and change frame
          * @param referenceFrame
          * @return Copied spatial transform with frame changed
          */
          SpatialForce changeFrameAndCopy(std::shared_ptr<ReferenceFrame> referenceFrame) const
          {
              return changeFrameAndCopy(referenceFrame.get());
          }

          /**
           * @brief Copy and change frame
           * @param referenceFrame
           * @return Copied spatial transform with frame changed
           */
          SpatialForce changeFrameAndCopy(ReferenceFrame* referenceFrame) const
          {
              SpatialForce ret = *this;
              ret.changeFrame(referenceFrame);
              return ret;
          }

          /**
           * @brief Copy this into spatial motion and change frame
           * @param referenceFrame
           * @param spatialForce Storage for copied result
           */
          void changeFrameAndCopy(ReferenceFrame* referenceFrame, SpatialForce &spatialForce) const
          {
              spatialForce = *this;
              spatialForce.changeFrame(referenceFrame);
          }

          /**
           * @brief Copy this into spatial motion and change frame
           * @param referenceFrame
           * @param spatialForce Storage for copied result
           */
          void changeFrameAndCopy(std::shared_ptr<ReferenceFrame> referenceFrame, SpatialForce &spatialForce) const
          {
              this->changeFrameAndCopy(referenceFrame.get(),spatialForce);
          }

            /**
             * @brief Get linear part of spatial force as a frame vector
             * @return FrameVector consisting of the reference frame and the linear portion
             */
            FrameVector getFramedLinearPart() const
            {
                return FrameVector(this->referenceFrame, this->getLinearPart());
            }

            /**
             * @brief Get angular part of spatial force as a frame vector
             * @return FrameVector consisting of the reference frame and the angular portion
             */
            FrameVector getFramedAngularPart() const
            {
                return FrameVector(this->referenceFrame, this->getAngularPart());
            }

            /**
             * @brief Copy constructor
             * @param spatialForce
             */
            SpatialForce(
                const SpatialForce& spatialForce) : ForceVector(spatialForce), FrameObject(spatialForce.referenceFrame)
            {}

            /**
             * @brief Constructor
             * @param referenceFrame ReferenceFrame this force vector is expressed in
             * @param spatialVector
             */
            SpatialForce(ReferenceFrame      *referenceFrame,
                         const SpatialVector& spatialVector) : ForceVector(spatialVector), FrameObject(referenceFrame)
            {}

            /**
             * @brief Get copy of this SpatialForce as type ForceVector
             * @return A SpatialForce as type ForceVector
             */
            inline ForceVector toForceVector() const
            {
                return ForceVector(this->mx(), this->my(), this->mz(), this->fx(), this->fy(), this->fz());
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
             * @brief Overloaded += operator. Frame checks are performed
             * @param f
             * @return \f$ f_sp=f_sp+f \f$
             */
            inline SpatialForce operator+=(const SpatialForce& f)
            {
                this->checkReferenceFramesMatch(&f);

                this->mx() += f.mx();
                this->my() += f.my();
                this->mz() += f.mz();

                this->fx() += f.fx();
                this->fy() += f.fy();
                this->fz() += f.fz();

                return *this;
            }

            /**
             * @brief Overloaded -= operator. Frame checks are performed
             * @param f
             * @return \f$ f_sp=f_sp-f \f$
             */
            inline SpatialForce operator-=(const SpatialForce& f)
            {
                this->checkReferenceFramesMatch(&f);

                this->mx() -= f.mx();
                this->my() -= f.my();
                this->mz() -= f.mz();

                this->fx() -= f.fx();
                this->fy() -= f.fy();
                this->fz() -= f.fz();

                return *this;
            }
        };

        /**
         * @brief Overloaded + operator to add two SpatialForce vectors and return the result in a new SpatialForce. Frame
         **checks are performed
         * @param f1
         * @param f2
         * @return \f$ f_sp=f_sp_1+f_sp_2 \f$
         */
        inline SpatialForce operator+(SpatialForce f1, const SpatialForce& f2)
        {
            return f1 += f2;
        }

        /**
         * @brief Overloaded - operator to add two SpatialForce vectors and return the result in a new SpatialForce. Frame
         **checks are performed
         * @param f1
         * @param f2
         * @return \f$ f_sp=f_sp_1-f_sp_2 \f$
         */
        inline SpatialForce operator-(SpatialForce f1, const SpatialForce& f2)
        {
            return f1 -= f2;
        }
    }
}

#endif // ifndef __RDL_SPATIAL_FORCE_HPP__
