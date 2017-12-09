/*
 * Original Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 *
 * RDL - Robot Dynamics Library
 * Modifications Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef __RDL_QUATERNION_H__
#define __RDL_QUATERNION_H__

#include <cmath>
#include "rdl_dynamics/rdl_eigenmath.h"

namespace RobotDynamics
{
    namespace Math
    {
        /** \brief Quaternion that are used for \ref joint_singularities "singularity free" joints.
         *
         * order: x,y,z,w
         */
        class Quaternion : public Vector4d
        {
public:

            /**
             * @brief Constructor
             */
            Quaternion() : Vector4d(0., 0., 0., 1.)
            {}

            Quaternion(const Eigen::Quaterniond& q) : Quaternion(q.x(), q.y(), q.z(), q.w())
            {}

            /**
             * @brief Constructor
             * @param vec4
             */
            Quaternion(const Vector4d& vec4) : Vector4d(vec4)
            {}

            /**
             * @brief Constructor
             * @param x
             * @param y
             * @param z
             * @param w
             */
            Quaternion(double x, double y, double z, double w) : Vector4d(x, y, z, w)
            {}

            EIGEN_STRONG_INLINE double& x()
            {
                return this->data()[0];
            }

            EIGEN_STRONG_INLINE double x() const
            {
                return this->data()[0];
            }

            EIGEN_STRONG_INLINE double& y()
            {
                return this->data()[1];
            }

            EIGEN_STRONG_INLINE double y() const
            {
                return this->data()[1];
            }

            EIGEN_STRONG_INLINE double& z()
            {
                return this->data()[2];
            }

            EIGEN_STRONG_INLINE double z() const
            {
                return this->data()[2];
            }

            EIGEN_STRONG_INLINE double& w()
            {
                return this->data()[3];
            }

            EIGEN_STRONG_INLINE double w() const
            {
                return this->data()[3];
            }

            void set(double x, double y, double z, double w)
            {
                this->data()[0] = x;
                this->data()[1] = y;
                this->data()[2] = z;
                this->data()[3] = w;
            }

            /**
             * @brief Get vector part
             * @return Vector part
             */
            EIGEN_STRONG_INLINE Vector3d getVectorPart() const
            {
                return Vector3d(this->data()[0], this->data()[1], this->data()[2]);
            }

            /**
             * @brief Get scalar part
             * @return Vector part
             */
            EIGEN_STRONG_INLINE double getScalarPart() const
            {
                return this->data()[3];
            }

            void operator=(const Eigen::Quaterniond& q)
            {
                set(q.x(), q.y(), q.z(), q.w());
            }

            /**
             * @brief Method to scale the elements of a quaternion by a constant. Normalization is NOT performed
             * @param s
             * @return Scaled quaternion
             */
            Quaternion operator*(const double& s) const
            {
                return Quaternion((*this)[0] * s, (*this)[1] * s, (*this)[2] * s, (*this)[3] * s);
            }

            /**
             * @brief Quaternion multiplication
             * @param q Quaternion to multiply by
             * @return New multiplied quaternion result
             */
            Quaternion operator*(const Quaternion& q) const
            {
                return Quaternion((*this)[3] * q[0] + (*this)[0] * q[3] + (*this)[1] * q[2] - (*this)[2] * q[1],
                                  (*this)[3] * q[1] + (*this)[1] * q[3] + (*this)[2] * q[0] - (*this)[0] * q[2],
                                  (*this)[3] * q[2] + (*this)[2] * q[3] + (*this)[0] * q[1] - (*this)[1] * q[0],
                                  (*this)[3] * q[3] - (*this)[0] * q[0] - (*this)[1] * q[1] - (*this)[2] * q[2]);
            }

            /**
             * @brief Overloaded *= operator for quaternion multiplication
             * @param q Quaternion to multiply by
             * @return Modified result of the multiplication
             */
            Quaternion& operator*=(const Quaternion& q)
            {
                set((*this)[3] * q[0] + (*this)[0] * q[3] + (*this)[1] * q[2] - (*this)[2] * q[1],
                    (*this)[3] * q[1] +
                    (*this)[1] * q[3] + (*this)[2] * q[0] - (*this)[0] * q[2],
                    (*this)[3] * q[2] + (*this)[2] * q[3] + (*this)[0] * q[1] - (*this)[1] * q[0],
                    (*this)[3] * q[3] - (*this)[0] * q[0] - (*this)[1] * q[1] - (*this)[2] * q[2]);
                return *this;
            }

            /**
             * @brief From Wikipedia: In computer graphics, Slerp is shorthand for spherical linear interpolation,
             * introduced by Ken Shoemake in the context of quaternion interpolation for the
             * purpose of animating 3D rotation. It refers to constant-speed motion along a unit-radius
             * great circle arc, given the ends and an interpolation parameter between 0 and 1
             * @note Only unit quaternions are valid rotations, so make sure to normalize
             * @param alpha Interpolation parameter. Should be between 0 and 1
             * @param quat Quaternion to interpolate between
             * @return Interpolated quaternion
             */
            Quaternion slerp(double alpha, const Quaternion& quat) const
            {
                // check whether one of the two has 0 length
                double s = std::sqrt(squaredNorm() * quat.squaredNorm());

                // division by 0.f is unhealthy!
                assert(s != 0.);

                double angle = acos(dot(quat) / s);

                if ((angle == 0.) || std::isnan(angle))
                {
                    return *this;
                }
                assert(!std::isnan(angle));

                double d  = 1. / std::sin(angle);
                double p0 = std::sin((1. - alpha) * angle);
                double p1 = std::sin(alpha * angle);

                if (dot(quat) < 0.)
                {
                    return Quaternion(((*this) * p0 - quat * p1) * d);
                }
                return Quaternion(((*this) * p0 + quat * p1) * d);
            }

            static Quaternion fromAxisAngle(const Vector3d& axis, double angle_rad)
            {
                double d  = axis.norm();
                double s2 = std::sin(angle_rad * 0.5) / d;

                return Quaternion(axis[0] * s2, axis[1] * s2, axis[2] * s2, std::cos(angle_rad * 0.5));
            }

            /**
             * @brief Creates a quaternion from a rotation matrix.
             * @note This method assumes the argument is an orthogonal matrix
             * @param mat Orothogonal rotation matrix
             * @return
             */
            static inline Quaternion fromMatrix(const Matrix3d& mat)
            {
                Quaternion q;

                fromMatrix(mat, q);
                return q;
            }

            static void fromMatrix(const Matrix3d& mat, Quaternion& quat)
            {
                double trace = mat.trace();

                if (trace > 0.)
                {
                    double s = 2. * std::sqrt(trace + 1.);
                    quat.set((mat(1, 2) - mat(2, 1)) / s, (mat(2, 0) - mat(0, 2)) / s, (mat(0, 1) - mat(1,
                                                                                                        0)) / s,
                             0.25 * s);
                }
                else if ((mat(0, 0) > mat(1, 1)) && (mat(0, 0) > mat(2, 2)))
                {
                    double s = 2. * std::sqrt(1. + mat(0, 0) - mat(1, 1) - mat(2, 2));
                    quat.set(-0.25 * s,
                             (-mat(0, 1) - mat(1, 0)) / s, (-mat(0, 2) - mat(2, 0)) / s, (mat(2, 1) - mat(1, 2)) / s);
                }
                else if (mat(1, 1) > mat(2, 2))
                {
                    double s = 2. * std::sqrt(1. + mat(1, 1) - mat(0, 0) - mat(2, 2));
                    quat.set((-mat(0,
                                   1) - mat(1, 0)) / s, -0.25 * s, (-mat(1, 2) - mat(2, 1)) / s, (mat(0, 2) - mat(2,
                                                                                                                  0))
                             / s);
                }
                else
                {
                    double s = 2. * std::sqrt(1. + mat(2, 2) - mat(0, 0) - mat(1, 1));
                    quat.set((-mat(0,
                                   2) - mat(2, 0)) / s, (-mat(1, 2) - mat(2, 1)) / s, -0.25 * s, (mat(1, 0) - mat(0,
                                                                                                                  1))
                             / s);
                }
            }

            static Quaternion fromZYXAngles(const Vector3d& zyx_angles)
            {
                return Quaternion::fromAxisAngle(Vector3d(0., 0.,
                                                          1.), zyx_angles[0]) * Quaternion::fromAxisAngle(Vector3d(0.,
                                                                                                                   1.,
                                                                                                                   0.),
                                                                                                          zyx_angles[1]) *
                       Quaternion::fromAxisAngle(
                           Vector3d(1., 0., 0.),
                           zyx_angles[2]);
            }

            static Quaternion fromYXZAngles(const Vector3d& yxz_angles)
            {
                return Quaternion::fromAxisAngle(Vector3d(0., 1.,
                                                          0.), yxz_angles[0]) * Quaternion::fromAxisAngle(Vector3d(1.,
                                                                                                                   0.,
                                                                                                                   0.),
                                                                                                          yxz_angles[1]) *
                       Quaternion::fromAxisAngle(
                           Vector3d(0., 0., 1.),
                           yxz_angles[2]);
            }

            static Quaternion fromXYZAngles(const Vector3d& xyz_angles)
            {
                return Quaternion::fromAxisAngle(Vector3d(1., 0.,
                                                          0.), xyz_angles[0]) * Quaternion::fromAxisAngle(Vector3d(0.,
                                                                                                                   1.,
                                                                                                                   0.),
                                                                                                          xyz_angles[1]) *
                       Quaternion::fromAxisAngle(
                           Vector3d(0., 0., 1.),
                           xyz_angles[2]);
            }

            Matrix3d toMatrix() const
            {
                double x = (*this)[0];
                double y = (*this)[1];
                double z = (*this)[2];
                double w = (*this)[3];

                return Matrix3d(1. - 2. * (y * y + z * z), 2. * x * y + 2. * w * z, 2. * x * z - 2. * w * y,
                                2. * x * y - 2. * w * z, 1. - 2. * (x * x + z * z), 2. * y * z + 2. * w * x,
                                2. * x * z + 2. * w * y, 2. * y * z - 2. * w * x, 1. - 2. * (x * x + y * y)
                                );
            }

            Quaternion conjugate() const
            {
                return Quaternion(-(*this)[0], -(*this)[1], -(*this)[2], (*this)[3]);
            }

            Quaternion timeStep(const Vector3d& omega, double dt)
            {
                double omega_norm = omega.norm();

                return Quaternion::fromAxisAngle(omega / omega_norm, dt * omega_norm) * (*this);
            }

            Vector3d rotate(const Vector3d& vec) const
            {
                Vector3d   vn(vec);
                Quaternion vec_quat(vn[0], vn[1], vn[2], 0.f), res_quat;

                res_quat = vec_quat * (*this);
                res_quat = conjugate() * res_quat;

                return Vector3d(res_quat[0], res_quat[1], res_quat[2]);
            }
        };
    }
}

/* __RDL_QUATERNION_H__ */
#endif // ifndef __RDL_QUATERNION_H__
