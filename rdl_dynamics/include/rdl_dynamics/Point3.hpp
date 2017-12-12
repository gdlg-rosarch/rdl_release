/*
 * RDL - Robot Dynamics Library
 * Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef __RDL_POINT_3_HPP__
#define __RDL_POINT_3_HPP__

/**
 * @file Point3.hpp
 */

#include <math.h>
#include <stdexcept>
#include "rdl_dynamics/TransformableGeometricObject.hpp"

namespace RobotDynamics
{
    namespace Math
    {
        /**
         * @class Point3
         * @brief A generic 3D point
         */
        class RDL_DLLAPI Point3d : public Math::TransformableGeometricObject
        {
public:

            Point3d(const double x, const double y, const double z)
            {
                point[0] = x;
                point[1] = y;
                point[2] = z;
            }

            Point3d(const std::vector<double>& vector)
            {
                x() = vector[0];
                y() = vector[1];
                z() = vector[2];
            }

            EIGEN_STRONG_INLINE Point3d(const Point3d& point)
            {
                x() = point.x();
                y() = point.y();
                z() = point.z();
            }

            EIGEN_STRONG_INLINE Point3d(const double array[3])
            {
                x() = array[0];
                y() = array[1];
                z() = array[2];
            }

            EIGEN_STRONG_INLINE Point3d()
            {
                setToZero();
            }

            virtual ~Point3d()
            {}

            /**
             * @brief Create a skew symmetric matrix, m, from a 3d vector such that, given two vectors \f$v_1\f$ and \f$v_2\f$,
             * a 3rd vector which is the cross product of the first two is given by, \f$v_3=\tilde{v_1}v_2\f$. The \f$\sim\f$
             * operator is referred to in Featherstones RBDA as the 3d vector cross(\f$\times\f$) operator.
             * @return A skew symmetric matrix
             */
            Matrix3d toTildeForm()
            {
                return Matrix3d(0., -point[2], point[1], point[2], 0., -point[0], -point[1], point[0], 0.);
            }

            /**
             * @brief Performs in place point transform. Given a point, \f$p\f$, this performs \f$ p = -X.E X.r + X.E p \f$
             * @param X
             */
            inline void transform(const Math::SpatialTransform& X)
            {
                set(-X.E * X.r + X.E * this->vec());
            }

            EIGEN_STRONG_INLINE void set(const std::vector<double>& vector)
            {
                set(vector[0], vector[1], vector[2]);
            }

            EIGEN_STRONG_INLINE void set(const Point3d& point)
            {
                set(point.x(), point.y(), point.z());
            }

            void set(const Math::Vector3d& v)
            {
                set(v[0], v[1], v[2]);
            }

            void set(const double x, const double y, const double z)
            {
                point[0] = x;
                point[1] = y;
                point[2] = z;
            }

            EIGEN_STRONG_INLINE void setToZero()
            {
                set(0., 0., 0.);
            }

            /**
             * Compares if two points are within epsilon of each other
             * @param point
             * @param epsilon
             * @return true if they are within epsilon, false otherwise
             */
            EIGEN_STRONG_INLINE bool epsilonEquals(const Point3d& point, const double epsilon) const
            {
                return fabs(this->x() - point.x()) < epsilon && fabs(this->y() - point.y()) < epsilon && fabs(
                           this->z() - point.z()) < epsilon;
            }

            /**
             * @brief clamp any values that are less than min to min
             * @param min
             */
            void clampMin(const double min)
            {
                if (x() < min)
                {
                    x() = min;
                }

                if (y() < min)
                {
                    y() = min;
                }

                if (z() < min)
                {
                    z() = min;
                }
            }

            /**
             * @brief clamp any values that are greater than make to max
             * @param max
             */
            void clampMax(const double max)
            {
                if (x() > max)
                {
                    x() = max;
                }

                if (y() > max)
                {
                    y() = max;
                }

                if (z() > max)
                {
                    z() = max;
                }
            }

            /**
             * @brief clamp any values greater than max to max, and any value less than min to min
             * @param min
             * @param max
             */
            void clampMinMax(const double min, const double max)
            {
                this->clampMin(min);
                this->clampMax(max);
            }

            /**
             * @brief Set each element to the absolute value
             */
            void absoluteValue()
            {
                this->x() = fabs(this->x());
                this->y() = fabs(this->y());
                this->z() = fabs(this->z());
            }

            /**
             * @brief Square of the distance between two ponts, \f$ x^2 + y^2 + z^2 \f$
             * @param point
             * @return Distance squared
             */
            double distanceSquared(const Point3d& point) const
            {
                double dx = x() - point.x();
                double dy = y() - point.y();
                double dz = z() - point.z();

                return dx * dx + dy * dy + dz * dz;
            }

            /**
             * brief Distance between two points, \f$ \sqrt{x^2 + y^2 + z^2} \f$
             * @param point
             * @return Distance
             */
            double distance(const Point3d& point) const
            {
                return sqrt(distanceSquared(point));
            }

            /**
             * @brief L1 norm of two points
             * @param point
             * @return L1 norm
             */
            double distanceL1(const Point3d& point) const
            {
                return fabs(x() - point.x()) + fabs(y() - point.y()) + fabs(z() - point.z());
            }

            /**
             * @brief Cross product between a point and vector
             * @param v
             * @return
             */
            Vector3d cross(const Vector3d& v)
            {
                return Vector3d(point[1] * v[2] - point[2] * v[1],
                                point[2] * v[0] - point[0] * v[2],
                                point[0] * v[1] - point[1] * v[0]);
            }

            /**
             * L-infinity norm
             * @param point
             * @return
             */
            double distanceLinf(const Point3d& point) const
            {
                double dx = x() - point.x();
                double dy = y() - point.y();
                double dz = z() - point.z();

                double tmp = fabs(dx) > fabs(dy) ? fabs(dx) : fabs(dy);

                return tmp > fabs(dz) ? tmp : fabs(dz);
            }

            EIGEN_STRONG_INLINE double& x()
            {
                return point[0];
            }

            EIGEN_STRONG_INLINE double x() const
            {
                return point[0];
            }

            EIGEN_STRONG_INLINE double& y()
            {
                return point[1];
            }

            EIGEN_STRONG_INLINE double y() const
            {
                return point[1];
            }

            EIGEN_STRONG_INLINE double& z()
            {
                return point[2];
            }

            EIGEN_STRONG_INLINE double z() const
            {
                return point[2];
            }

            EIGEN_STRONG_INLINE double* data()
            {
                return point;
            }

            EIGEN_STRONG_INLINE Math::Vector3d vec() const
            {
                return Math::Vector3d(point[0], point[1], point[2]);
            }

            void operator+=(const Point3d& point)
            {
                this->x() += point.x();
                this->y() += point.y();
                this->z() += point.z();
            }

            void operator-=(const Point3d& point)
            {
                this->x() -= point.x();
                this->y() -= point.y();
                this->z() -= point.z();
            }

            void operator*=(const double scale)
            {
                this->x() *= scale;
                this->y() *= scale;
                this->z() *= scale;
            }

            void operator/=(const double scale)
            {
                this->x() /= scale;
                this->y() /= scale;
                this->z() /= scale;
            }

            bool operator==(const Point3d& rhs)
            {
                if ((this->x() != rhs.x()) || (this->y() != rhs.y()) || (this->z() != rhs.z()))
                {
                    return false;
                }

                return true;
            }

            bool operator!=(const Point3d& rhs)
            {
                return !this->operator==(rhs);
            }

protected:

            double point[3];
        };

        inline Point3d operator+(Point3d leftHandSide, const Point3d& point)
        {
            leftHandSide += point;
            return leftHandSide;
        }

        inline Point3d operator-(Point3d leftHandSide, const Point3d& point)
        {
            leftHandSide -= point;
            return leftHandSide;
        }

        inline Point3d operator*(Point3d leftHandSide, const double scale)
        {
            leftHandSide *= scale;
            return leftHandSide;
        }

        inline std::ostream& operator<<(std::ostream& os, const Point3d& point)
        {
            os << "x: " << point.x() << '\n' << "y: " << point.y() << '\n' << "z: " << point.z() << "\n";
            return os;
        }
    }
}
#endif // ifndef __RDL_POINT_3_HPP__
