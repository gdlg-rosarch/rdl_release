/*
 * RDL - Robot Dynamics Library
 * Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

/**
 * @file FramePoint.hpp
 * @brief File containing the FramePoint<T> object definition
 */

#ifndef __RDL_FRAME_POINT_HPP__
#define __RDL_FRAME_POINT_HPP__

#include "rdl_dynamics/FrameObject.hpp"
#include "rdl_dynamics/Point3.hpp"
#include "rdl_dynamics/FrameExceptions.hpp"

namespace RobotDynamics
{
    namespace Math
    {
        /**
         * @class FramePoint
         * @ingroup reference_frame
         * @brief A FramePoint is a 3D point that is expressed in a ReferenceFrame. To change the ReferenceFrame a
         * FramePoint is expressed in, you may call the inhereted FrameObject::changeFrame method and supply it a
         * pointer to the ReferenceFrame you wish to have the FramePoint expressed in. This class and its
         * implementation are an adaptation of FramePoint.java by <a href="http://robots.ihmc.us/">Jerry Pratt and the IHMC
         **Robotics Group</a>.
         */
        class RDL_DLLAPI FramePoint : public FrameObject, public Math::Point3d
        {
public:

            /**
             * @brief Constructor
             * @param referenceFrame A pointer to the ReferenceFrame the point will be expressed in
             * @param x The x-component of the point
             * @param y The y-component of the point
             * @param z The z-component of the point
             */
            FramePoint(ReferenceFrame *referenceFrame, const double x, const double y,
                       const double z) : FrameObject(referenceFrame), Math::Point3d(x, y, z)
            {}

            /**
             * @brief Constructor
             * @param referenceFrame A pointer to the ReferenceFrame the point will be expressed in
             * @param array An array that will be used to set the components of this FramePoint
             */
            FramePoint(ReferenceFrame *referenceFrame, double array[3]) : FrameObject(referenceFrame), Math::Point3d(array)
            {}

            /**
             * @brief Constructor
             * @param referenceFrame A pointer to the ReferenceFrame the point will be expressed in
             * @param v A Vector3d that will be used to set the components of this FramePoint
             */
            FramePoint(ReferenceFrame *referenceFrame,
                       Math::Vector3d  v) : FrameObject(referenceFrame), Math::Point3d(v[0], v[1], v[2])
            {}

            /**
             * @brief Constructor
             * @param referenceFrame A pointer to the ReferenceFrame the point will be expressed in
             * @param vector A templated stl vector that will be used to set the components of this FramePoint
             */
            FramePoint(ReferenceFrame    *referenceFrame,
                       std::vector<double>vector) : FrameObject(referenceFrame), Math::Point3d(vector)
            {}

            /**
             * @brief Constructor
             * @param referenceFrame A pointer to the ReferenceFrame the point will be expressed in
             * @param point A Math::Point3 that will be used to set the components of this FramePoint
             */
            FramePoint(ReferenceFrame      *referenceFrame,
                       const Math::Point3d& point) : FrameObject(referenceFrame), Math::Point3d(point)
            {}

            /**
             * @brief Copy constructor
             * @param framePoint A FramePoint to copy
             */
            FramePoint(
                const FramePoint& framePoint) : FrameObject(framePoint.getReferenceFrame()),
                                                Math::Point3d(framePoint.x(), framePoint.y(), framePoint.z())
            {}

            /**
             * @brief Constructor that initializes to (x,y,z) = (0,0,0)
             * @param referenceFrame A pointer to the ReferenceFrame the point will be expressed in
             */
            FramePoint(ReferenceFrame *referenceFrame) : FrameObject(referenceFrame), Math::Point3d()
            {}

            /**
             * @brief Empty constructor that creates a point with ReferencFrame=nullptr and (x,y,z)=(0,0,0)
             */
            FramePoint() : FrameObject(nullptr), Math::Point3d()
            {}

            /**
             * @brief Destructor
             */
            ~FramePoint()
            {}

            /**
             * @brief Return a pointer to this as base class type Math::TransformableGeometricObject. See
             **FrameObject::changeFrame for how this method is used
             * @return Pointer to this object as type Math::TransformableGeometricObject
             */
            Math::TransformableGeometricObject* getTransformableGeometricObject()
            {
                return this;
            }

            /**
             * @brief copy into new frame point and change the frame of that
             * @param referenceFrame
             * @return
             */
            FramePoint changeFrameAndCopy(ReferenceFrame* referenceFrame) const
            {
                FramePoint p = *this;
                p.changeFrame(referenceFrame);
                return p;
            }

            /**
             * @brief copy into new frame point and change the frame of that
             * @param referenceFrame
             * @return
             */
            FramePoint changeFrameAndCopy(std::shared_ptr<ReferenceFrame> referenceFrame) const
            {
                return changeFrameAndCopy(referenceFrame.get());
            }

            /**
             * @brief Copy *this into p and change its frame
             * @param referenceFrame
             * @param p Modified
             * @return
             */
            void changeFrameAndCopy(ReferenceFrame* referenceFrame, FramePoint &p) const
            {
                p = *this;
                p.changeFrame(referenceFrame);
            }

            /**
             * @brief Copy *this into p and change its frame
             * @param referenceFrame
             * @param p Modified
             * @return
             */
            void changeFrameAndCopy(std::shared_ptr<ReferenceFrame> referenceFrame, FramePoint &p) const
            {
                changeFrameAndCopy(referenceFrame.get(),p);
            }

            /**
             * @brief Set both the ReferenceFrame this object is expressed in as well as the (x,y,z) coordinates of the point
             * @param v Vector3d that this point will be set to
             * @param referenceFrame Pointer to the ReferenceFrame this object will be expressed in
             */
            EIGEN_STRONG_INLINE void setIncludingFrame(const Math::Vector3d& v, ReferenceFrame *referenceFrame)
            {
                setIncludingFrame(v(0), v(1), v(2), referenceFrame);
            }

            /**
             * @brief Set both the ReferenceFrame the point is expressed in as well as the (x,y,z) coordinates
             * @param x The x coordinate
             * @param y The y coordinate
             * @param z The z coordinate
             * @param referenceFrame The ReferenceFrame this point is to be expressed in
             */
            void setIncludingFrame(const double x, const double y, const double z, ReferenceFrame *referenceFrame)
            {
                if (!referenceFrame)
                {
                    throw ReferenceFrameException("Reference frame is nullptr!");
                }

                this->set(x, y, z);
                this->referenceFrame = referenceFrame;
            }

            /**
             * @brief Set both the ReferenceFrame the point is expressed in as well as the (x,y,z) coordinates
             * @param point Math::Point3d to set this point to
             * @param referenceFrame Pointer to ReferenceFrame this point will be expressed in
             */
            void setIncludingFrame(const Math::Point3d& point, ReferenceFrame *referenceFrame)
            {
                if (!referenceFrame)
                {
                    throw ReferenceFrameException("Reference frame cannot be nullptr!");
                }

                this->x()            = point.x();
                this->y()            = point.y();
                this->z()            = point.z();
                this->referenceFrame = referenceFrame;
            }

            /**
             * @brief Calculate the distance squared between two FramePoints.  \f$\Delta_x^2+\Delta_y^2+\Delta_z^2\f$
             * @throws ReferenceFrameException If both points are not expressed in the same ReferenceFrame
             * @param point FramePoint to calculate squared distance to
             * @return Distance squared
             */
            double distanceSquared(const FramePoint& point) const
            {
                checkReferenceFramesMatch(&point);

                double dx = this->x() - point.x();
                double dy = this->y() - point.y();
                double dz = this->z() - point.z();
                return dx * dx + dy * dy + dz * dz;
            }

            /**
             * @brief Calculate the distance between two FramePoints. \f$\sqrt{\Delta_x^2+\Delta_y^2+\Delta_z^2}\f$
             * @throws ReferenceFrameException If both points are not expressed in the same ReferenceFrame
             * @param point FramePoint to calculate distance to
             * @return Distance between points as template type T
             */
            double distance(const FramePoint& point) const
            {
                checkReferenceFramesMatch(&point);

                return sqrt(distanceSquared(point));
            }

            /**
             * @brief Calculate the L1 distance between two FramePoints by \f$|\Delta_x| + |\Delta_y| + |\Delta_z|\f$
             * @throws ReferenceFrameException If both points are not expressed in the same ReferenceFrame
             * @param point FramePoint to calculate distance to
             * @return Distance between points as template type T
             */
            double distanceL1(const FramePoint& point) const
            {
                checkReferenceFramesMatch(&point);

                return fabs(this->x() - point.x()) + fabs(this->y() - point.y()) + fabs(this->z() - point.z());
            }

            /**
             * @brief Calculate the LInfinity distance between two FramePoints by \f$max(|\Delta_x|,|\Delta_y|,|\Delta_z|)\f$
             * @throws ReferenceFrameException If both points are not expressed in the same ReferenceFrame
             * @param point FramePoint to calculate distance to
             * @return Distance between points as template type T
             */
            double distanceLinf(const FramePoint& point) const
            {
                checkReferenceFramesMatch(&point);

                double dx = this->x() - point.x();
                double dy = this->y() - point.y();
                double dz = this->z() - point.z();

                double tmp = fabs(dx) > fabs(dy) ? fabs(dx) : fabs(dy);

                return tmp > fabs(dz) ? tmp : fabs(dz);
            }

            /**
             * @brief Perform addition of FramePoint argument
             * @param point The FramePoint to be added
             * @throws ReferenceFrameException If both points are not expressed in the same ReferenceFrame
             */
            void add(const FramePoint& point)
            {
                checkReferenceFramesMatch(&point);
                this->x() += point.x();
                this->y() += point.y();
                this->z() += point.z();
            }

            /**
             * @brief Perform subtraction of FramePoint argument
             * @param point The FramePoint to be subtracted
             * @throws ReferenceFrameException If both points are not expressed in the same ReferenceFrame
             */
            void subtract(const FramePoint& point)
            {
                checkReferenceFramesMatch(&point);
                this->x() -= point.x();
                this->y() -= point.y();
                this->z() -= point.z();
            }

            /**
             * @brief Return true FramePoint argument is within epsilon of this, false otherwise
             * @param point The FramePoint to be compared
             * @param epsilon The tolerance of the comparison check
             * @throws ReferenceFrameException If both points are not expressed in the same ReferenceFrame
             */
            bool epsilonEquals(const FramePoint& point, const double epsilon) const
            {
                checkReferenceFramesMatch(&point);
                return fabs(this->x() - point.x()) < epsilon && fabs(this->y() - point.y()) < epsilon &&
                       fabs(this->z() - point.z()) < epsilon;
            }

            /**
             * @brief Overloaded += operator, performs this = this + point
             * @param point FramePoint to be added
             * @return FramePoint representing this = this + point
             */
            void operator+=(const FramePoint& point)
            {
                checkReferenceFramesMatch(&point);
                this->x() += point.x();
                this->y() += point.y();
                this->z() += point.z();
            }

            /**
             * @brief Overloaded += operator, performs this = this + vec
             * @param vec Vector3d to be added
             * @return FramePoint representing this = this + vec
             * @note No frame checking/changing occers
             */
            void operator+=(const Vector3d& vec)
            {
                this->x() += vec.x();
                this->y() += vec.y();
                this->z() += vec.z();
            }

            /**
             * @brief Overloaded -= operator, performs this = this - point
             * @param point FramePoint to be subtracted
             * @return FramePoint representing this = this - point
             */
            void operator-=(const FramePoint& point)
            {
                checkReferenceFramesMatch(&point);
                this->x() -= point.x();
                this->y() -= point.y();
                this->z() -= point.z();
            }

            /**
             * @brief Overloaded -= operator, performs this = this - vec
             * @param vec Vector3d to be subtracted
             * @return FramePoint representing this = this - vec
             * @note No frame checking/changing occers
             */
            void operator-=(const Vector3d& vec)
            {
                this->x() -= vec.x();
                this->y() -= vec.y();
                this->z() -= vec.z();
            }

            /**
             * @brief Overloaded *= operator, performs this = this*scala
             * @param scale Scalar to scale each element of this FramePoint by
             * @return FramePoint representing this = this*scale
             */
            void operator*=(const double scale)
            {
                this->x() *= scale;
                this->y() *= scale;
                this->z() *= scale;
            }

            /**
             * @brief Overloaded /= operator, performs this = this*scala
             * @param scale Scalar to divide each element of this FramePoint by
             * @return FramePoint representing this = this/scale
             */
            void operator/=(const double scale)
            {
                this->x() /= scale;
                this->y() /= scale;
                this->z() /= scale;
            }

protected:
        };

        /**
         * @brief Check if two FramePoints are equal
         * @param lhs
         * @param rhs
         * @throws ReferenceFrameException If the two FramePoint arguments are not expressed in the same ReferenceFrame
         * @return bool true if equal, false if not equal
         */
        inline bool operator==(const FramePoint& lhs, const FramePoint& rhs)
        {
            lhs.checkReferenceFramesMatch(rhs.getReferenceFrame());

            if (lhs.x() != rhs.x())
            {
                return false;
            }

            if (lhs.y() != rhs.y())
            {
                return false;
            }

            if (lhs.z() != rhs.z())
            {
                return false;
            }

            return true;
        }

        /**
         * @brief Add two FramePoints and return result in newly created FramePoint
         * @param p1
         * @param p2
         * @throws ReferenceFrameException If the two FramePoint arguments are not expressed in the same ReferenceFrame
         * @return A FramePoint that is the addition of the two argument FramePoints, i.e. p1+=p2
         */
        inline FramePoint operator+(FramePoint p1, const FramePoint& p2)
        {
            p1 += p2;
            return p1;
        }

        /**
         * @brief Subtract two FramePoints and return result in newly created FramePoint
         * @param p1
         * @param p2
         * @throws ReferenceFrameException If the two FramePoint arguments are not expressed in the same ReferenceFrame
         * @return A FramePoint that is the difference of the two argument FramePoints, i.e. p1-=p2
         */
        inline FramePoint operator-(FramePoint p1, const FramePoint& p2)
        {
            p1 -= p2;
            return p1;
        }

        /**
         * @brief Check if two FramePoints are not equal
         * @param lhs
         * @param rhs
         * @throws ReferenceFrameException If the two FramePoint arguments are not expressed in the same ReferenceFrame
         * @return bool false if equal, true if not equal
         */
        inline bool operator!=(const FramePoint& lhs, const FramePoint& rhs)
        {
            return !operator==(lhs, rhs);
        }

        inline std::ostream& operator<<(std::ostream& output, const FramePoint& framePoint)
        {
            output << "ReferenceFrame = " << framePoint.getReferenceFrame()->getName() << std::endl;
            output << "x = " << framePoint.x() << " y = " << framePoint.y() << " z = " << framePoint.z() << std::endl;
            return output;
        }

        /**
         * @typedef FramePoint<double> FramePointd
         * @brief A type definition for a FramePoint<double>
         */
        typedef FramePoint FramePointd;
    }
}
#endif // ifndef __RDL_FRAME_POINT_HPP__
