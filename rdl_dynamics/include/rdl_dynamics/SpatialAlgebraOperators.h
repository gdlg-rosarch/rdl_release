/*
 * Original Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 *
 * RDL - Robot Dynamics Library
 * Modifications Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

/**
 * @file SpatialAlgebraOperators.h
 */

#ifndef __RDL_SPATIALALGEBRAOPERATORS_H__
#define __RDL_SPATIALALGEBRAOPERATORS_H__

#include <iostream>
#include <cmath>
#include "rdl_dynamics/rdl_eigenmath.h"

namespace RobotDynamics
{
    namespace Math
    {
        /**
         * @brief Create a skew symmetric matrix, m, from a 3d vector such that, given two vectors \f$v_1\f$ and \f$v_2\f$,
         * a 3rd vector which is the cross product of the first two is given by, \f$v_3=\tilde{v_1}v_2\f$. The \f$\sim\f$
         * operator is referred to in Featherstones RBDA as the 3d vector cross(\f$\times\f$) operator.
         * @param vector
         * @return A skew symmetric matrix
         */
        inline Matrix3d toTildeForm(const Vector3d& vector)
        {
            return Matrix3d(0., -vector[2], vector[1], vector[2], 0., -vector[0], -vector[1], vector[0], 0.);
        }

        /** \brief Compact representation of spatial transformations.
         *
         * Instead of using a verbose 6x6 matrix, this structure only stores a 3x3
         * matrix and a 3-d vector to store spatial transformations. It also
         * encapsulates efficient operations such as concatenations and
         * transformation of spatial vectors.
         */
        struct RDL_DLLAPI SpatialTransform
        {
            /**
             * @brief Constructor
             */
            SpatialTransform() : E(Matrix3d::Identity(3, 3)), r(Vector3d::Zero(3, 1))
            {}

            /**
             * @brief Constructor
             * @param rotation Orthogonal rotation matrix
             * @param x X component
             * @param y Y component
             * @param z Z component
             */
            SpatialTransform(const Matrix3d& rotation, const double x, const double y, const double z) : E(rotation), r(x, y, z)
            {}

            /**
             * @brief Constructor
             * @param rotation Orthogonal rotation matrix
             * @param translation 3D translational component
             */
            SpatialTransform(const Matrix3d& rotation, const Vector3d& translation) : E(rotation), r(translation)
            {}

            /**
             * @brief Transform a spatial vector. Same as \f$ X * v \f$
             * @param v_sp Spatial motion vector to be copied/transformed
             *
             * @returns Transformed spatial vector. \f$ \begin{bmatrix} E * w \\ - E * (r \times w) + E * v \end{bmatrix} \f$
             */
            SpatialVector apply(const SpatialVector& v_sp) const
            {
                Vector3d v_rxw(v_sp[3] - r[1] * v_sp[2] + r[2] * v_sp[1],
                               v_sp[4] - r[2] * v_sp[0] + r[0] * v_sp[2],
                               v_sp[5] - r[0] * v_sp[1] + r[1] * v_sp[0]);

                return SpatialVector(E(0, 0) * v_sp[0] + E(0, 1) * v_sp[1] + E(0, 2) * v_sp[2], E(1, 0) * v_sp[0] + E(1,
                                                                                                                      1) *
                                     v_sp[1] +
                                     E(1, 2) * v_sp[2], E(2, 0) * v_sp[0] + E(2, 1) * v_sp[1] + E(2, 2) * v_sp[2], E(0,
                                                                                                                     0)
                                     * v_rxw[0] + E(0, 1) * v_rxw[1] + E(0, 2) * v_rxw[2], E(1, 0) * v_rxw[0] + E(1,
                                                                                                                  1) *
                                     v_rxw[1] + E(1,
                                                  2)
                                     * v_rxw[2], E(2, 0) * v_rxw[0] + E(2, 1) * v_rxw[1] + E(2, 2) * v_rxw[2]);
            }

            /**
             * @brief Applies \f$ X^T * f \f$
             * @param f_sp Spatial force
             *
             * @returns \f$ \begin{bmatrix} E^T * n + r \times * E^T * f \\ E^T * f \end{bmatrix} \f$
             */
            SpatialVector applyTranspose(const SpatialVector& f_sp) const
            {
                Vector3d E_T_f(E(0, 0) * f_sp[3] + E(1, 0) * f_sp[4] + E(2, 0) * f_sp[5], E(0, 1) * f_sp[3] + E(1,
                                                                                                                1) * f_sp[4] +
                               E(2, 1) * f_sp[5], E(0, 2) * f_sp[3] + E(1, 2) * f_sp[4] + E(2, 2) * f_sp[5]);

                return SpatialVector(E(0, 0) * f_sp[0] + E(1, 0) * f_sp[1] + E(2,
                                                                               0) * f_sp[2] - r[2] * E_T_f[1] + r[1] * E_T_f[2],
                                     E(0, 1) * f_sp[0] + E(1, 1) * f_sp[1] + E(2,
                                                                               1) * f_sp[2] + r[2] * E_T_f[0] - r[0] * E_T_f[2],
                                     E(0,
                                       2)
                                     * f_sp[0] + E(1, 2) * f_sp[1] + E(2,
                                                                       2)
                                     * f_sp[2] - r[1] * E_T_f[0] + r[0] * E_T_f[1], E_T_f[0], E_T_f[1], E_T_f[2]);
            }

            /**
             * @brief Applies \f$ X * f \f$ where \f$ f \f$ is a spatial force
             * @param f_sp Spatial force vector
             *
             * @return \f$ \begin{bmatrix} E * n - E * (r \times f) \\ E * f \end{bmatrix} \f$
             */
            SpatialVector applyAdjoint(const SpatialVector& f_sp) const
            {
                Vector3d En_rxf = E * (Vector3d(f_sp[0], f_sp[1], f_sp[2]) - r.cross(Vector3d(f_sp[3], f_sp[4], f_sp[5])));

                return SpatialVector(En_rxf[0], En_rxf[1], En_rxf[2], E(0, 0) * f_sp[3] + E(0, 1) * f_sp[4] + E(0,
                                                                                                                2) * f_sp[5],
                                     E(1, 0) * f_sp[3] + E(1, 1) * f_sp[4] + E(1, 2) * f_sp[5], E(2, 0) * f_sp[3] + E(2,
                                                                                                                      1)
                                     * f_sp[4] + E(2, 2) * f_sp[5]);
            }

            /**
             * @brief Return transform as 6x6 spatial matrix
             *
             * @return \f$ \begin{bmatrix} E & \mathbf{0} \\ -E * r\times & E \end{bmatrix} \f$
             */
            SpatialMatrix toMatrix() const
            {
                Matrix3d _Erx = E * Matrix3d(0., -r[2], r[1], r[2], 0., -r[0], -r[1], r[0], 0.);
                SpatialMatrix result;

                result.block<3, 3>(0, 0) = E;
                result.block<3, 3>(0, 3) = Matrix3d::Zero(3, 3);
                result.block<3, 3>(3, 0) = -_Erx;
                result.block<3, 3>(3, 3) = E;

                return result;
            }

            /**
             * @brief Returns Spatial transform that transforms spatial force vectors
             *
             * @return \f$ \begin{bmatrix} E & -E * r\times \\ \mathbf{0} & E \end{bmatrix} \f$
             */
            SpatialMatrix toMatrixAdjoint() const
            {
                Matrix3d _Erx = E * Matrix3d(0., -r[2], r[1], r[2], 0., -r[0], -r[1], r[0], 0.);
                SpatialMatrix result;

                result.block<3, 3>(0, 0) = E;
                result.block<3, 3>(0, 3) = -_Erx;
                result.block<3, 3>(3, 0) = Matrix3d::Zero(3, 3);
                result.block<3, 3>(3, 3) = E;

                return result;
            }

            /**
             * @brief Returns spatial force transform transposed
             *
             * @return \f$ \begin{bmatrix} E^{T} & (-E r\times)^{T} \\ \mathbf{0} & E^{T} \end{bmatrix} \f$
             */
            SpatialMatrix toMatrixTranspose() const
            {
                Matrix3d _Erx = E * Matrix3d(0., -r[2], r[1], r[2], 0., -r[0], -r[1], r[0], 0.);
                SpatialMatrix result;

                result.block<3, 3>(0, 0) = E.transpose();
                result.block<3, 3>(0, 3) = -_Erx.transpose();
                result.block<3, 3>(3, 0) = Matrix3d::Zero(3, 3);
                result.block<3, 3>(3, 3) = E.transpose();

                return result;
            }

            /**
             * @brief Returns inverse of transform
             *
             * @return \f$ X^{-1} \f$
             */
            SpatialTransform inverse() const
            {
                return SpatialTransform(E.transpose(), -E * r);
            }

            /**
             * @brief Inverts in place. \f$ this = this^{-1} \f$
             */
            void invert()
            {
                r = -E * r;
                E.transposeInPlace();
            }

            /**
             * @brief Overloaded * operator for combining transforms
             * @param XT
             * @return Combined rotation
             */
            SpatialTransform operator*(const SpatialTransform& XT) const
            {
                return SpatialTransform(E * XT.E, XT.r + XT.E.transpose() * r);
            }

            void operator*=(const SpatialTransform& XT)
            {
                r  = XT.r + XT.E.transpose() * r;
                E *= XT.E;
            }

            Matrix3d E;
            Vector3d r;
        };

        inline std::ostream& operator<<(std::ostream& output, const SpatialTransform& X)
        {
            output << "X.E = " << std::endl << X.E << std::endl;
            output << "X.r = " << X.r.transpose();
            return output;
        }

        /**
         * @brief Get spatial transform from angle and axis
         * @param angle_rad angle magnitude
         * @param axis normalized 3d vector
         * @return Spatial transform
         */
        inline SpatialTransform Xrot(double angle_rad, const Vector3d& axis)
        {
            double s, c;

            s = sin(angle_rad);
            c = cos(angle_rad);

            return SpatialTransform(Matrix3d(axis[0] * axis[0] * (1.0f - c) + c, axis[1] * axis[0] * (1.0f - c) + axis[2] * s,
                                             axis[0] * axis[2] * (1.0f - c) - axis[1] * s,

                                             axis[0] * axis[1] * (1.0f - c) - axis[2] * s, axis[1] * axis[1] * (1.0f - c) + c,
                                             axis[1] * axis[2] * (1.0f - c) + axis[0] * s,

                                             axis[0] * axis[2] *
                                             (1.0f -
                                              c) + axis[1] * s, axis[1] * axis[2] *
                                             (1.0f - c) - axis[0] * s, axis[2] * axis[2] * (1.0f - c) + c

                                             ), Vector3d(0., 0., 0.));
        }

        /**
         * @brief Get transform with zero translation and pure rotation about x axis
         * @param xrot
         * @return Transform with zero translation and x-rotation
         */
        inline SpatialTransform Xrotx(const double& xrot)
        {
            double s, c;

            s = sin(xrot);
            c = cos(xrot);
            return SpatialTransform(Matrix3d(1., 0., 0., 0., c, s, 0., -s, c), Vector3d(0., 0., 0.));
        }

        /**
         * @brief Get transform with zero translation and pure rotation about y axis
         * @param yrot
         * @return Transform with zero translation and y-rotation
         */
        inline SpatialTransform Xroty(const double& yrot)
        {
            double s, c;

            s = sin(yrot);
            c = cos(yrot);
            return SpatialTransform(Matrix3d(c, 0., -s, 0., 1., 0., s, 0., c), Vector3d(0., 0., 0.));
        }

        /**
         * @brief Get transform with zero translation and pure rotation about z axis
         * @param zrot
         * @return Transform with zero translation and z-rotation
         */
        inline SpatialTransform Xrotz(const double& zrot)
        {
            double s, c;

            s = sin(zrot);
            c = cos(zrot);
            return SpatialTransform(Matrix3d(c, s, 0., -s, c, 0., 0., 0., 1.), Vector3d(0., 0., 0.));
        }

        /**
         * @brief Get pure translation transform
         * @param r
         * @return Transform with identity rotation and translation \f$ r \f$
         */
        inline SpatialTransform Xtrans(const Vector3d& r)
        {
            return SpatialTransform(Matrix3d::Identity(3, 3), r);
        }

        /**
         * @brief Get the spatial motion cross matrix
         * @param v
         * @return \f$ v\times \f$
         */
        inline SpatialMatrix crossm(const SpatialVector& v)
        {
            return SpatialMatrix(0,
                                 -v[2],
                                 v[1],
                                 0,
                                 0,
                                 0,
                                 v[2],
                                 0,
                                 -v[0],
                                 0,
                                 0,
                                 0,
                                 -v[1],
                                 v[0],
                                 0,
                                 0,
                                 0,
                                 0,
                                 0,
                                 -v[5],
                                 v[4],
                                 0,
                                 -v[2],
                                 v[1],
                                 v[5],
                                 0,
                                 -v[3],
                                 v[2],
                                 0,
                                 -v[0],
                                 -v[4],
                                 v[3],
                                 0,
                                 -v[1],
                                 v[0],
                                 0);
        }

        /**
         * @brief Spatial motion cross times spatial motion
         * @param v1
         * @param v2
         * @return \f$ v1\times v2 \f$
         */
        inline SpatialVector crossm(const SpatialVector& v1, const SpatialVector& v2)
        {
            return SpatialVector(-v1[2] * v2[1] + v1[1] * v2[2],
                                 v1[2] * v2[0] - v1[0] * v2[2],
                                 -v1[1] * v2[0] + v1[0] * v2[1],
                                 -v1[5] * v2[1] + v1[4] * v2[2] - v1[2] * v2[4] + v1[1] * v2[5],
                                 v1[5] * v2[0] - v1[3] * v2[2] + v1[2] * v2[3] - v1[0] * v2[5],
                                 -v1[4] * v2[0] + v1[3] * v2[1] - v1[1] * v2[3] + v1[0] * v2[4]);
        }

        /**
         * @brief Get the spatial force cross matrix
         * @param v
         * @return \f$ v\times* \f$
         */
        inline SpatialMatrix crossf(const SpatialVector& v)
        {
            return SpatialMatrix(0,
                                 -v[2],
                                 v[1],
                                 0,
                                 -v[5],
                                 v[4],
                                 v[2],
                                 0,
                                 -v[0],
                                 v[5],
                                 0,
                                 -v[3],
                                 -v[1],
                                 v[0],
                                 0,
                                 -v[4],
                                 v[3],
                                 0,
                                 0,
                                 0,
                                 0,
                                 0,
                                 -v[2],
                                 v[1],
                                 0,
                                 0,
                                 0,
                                 v[2],
                                 0,
                                 -v[0],
                                 0,
                                 0,
                                 0,
                                 -v[1],
                                 v[0],
                                 0);
        }

        /**
         * @brief Spatial motion cross spatial force
         * @param v1 Spatial motion
         * @param v2 Spatial force
         * @return \f$ v1\times* v2 \f$
         */
        inline SpatialVector crossf(const SpatialVector& v1, const SpatialVector& v2)
        {
            return SpatialVector(-v1[2] * v2[1] + v1[1] * v2[2] - v1[5] * v2[4] + v1[4] * v2[5],
                                 v1[2] * v2[0] - v1[0] * v2[2] + v1[5] * v2[3] - v1[3] * v2[5],
                                 -v1[1] * v2[0] + v1[0] * v2[1] - v1[4] * v2[3] + v1[3] * v2[4],
                                 -v1[2] * v2[4] + v1[1] * v2[5],
                                 +v1[2] * v2[3] - v1[0] * v2[5],
                                 -v1[1] * v2[3] + v1[0] * v2[4]);
        }

        /**
         * @brief Get the rotated linear portion of the spatial vector
         * @param v Spatial vector
         * @param X Spatial transform
         * @return Rotated linear portion of the spatial vector argument
         */
        inline Vector3d getLinearPartTransformed(const SpatialVector& v, const SpatialTransform& X)
        {
            double v_rxw_x = v[3] - X.r[1] * v[2] + X.r[2] * v[1];
            double v_rxw_y = v[4] - X.r[2] * v[0] + X.r[0] * v[2];
            double v_rxw_z = v[5] - X.r[0] * v[1] + X.r[1] * v[0];

            return Vector3d(X.E(0, 0) * v_rxw_x + X.E(0, 1) * v_rxw_y + X.E(0, 2) * v_rxw_z,
                            X.E(1, 0) * v_rxw_x + X.E(1, 1) * v_rxw_y + X.E(1, 2) * v_rxw_z,
                            X.E(2, 0) * v_rxw_x + X.E(2, 1) * v_rxw_y + X.E(2, 2) * v_rxw_z);
        }
    } /* Math */
}     /* RobotDynamics */

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RobotDynamics::Math::SpatialTransform)

/* __RDL_SPATIALALGEBRAOPERATORS_H__*/
#endif // ifndef __RDL_SPATIALALGEBRAOPERATORS_H__
