/*
 * Original Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 *
 * RDL - Robot Dynamics Library
 * Modifications Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef __RDL_MATHUTILS_H__
#define __RDL_MATHUTILS_H__

/**
 * @file rdl_mathutils.h
 */

#include <assert.h>
#include <cmath>

#include "rdl_dynamics/rdl_eigenmath.h"

namespace RobotDynamics
{
    struct Model;

    namespace Math
    {
        /**
         * @enum LinearSolver
         * @brief Available solver methods for the linear systems.
         */
        enum RDL_DLLAPI LinearSolver
        {
            LinearSolverUnknown = 0,
            LinearSolverPartialPivLU,
            LinearSolverColPivHouseholderQR,
            LinearSolverHouseholderQR,
            LinearSolverLLT,
            LinearSolverLast,
        };

        extern RDL_DLLAPI Vector3d Vector3dZero;
        extern RDL_DLLAPI Matrix3d Matrix3dIdentity;
        extern RDL_DLLAPI Matrix3d Matrix3dZero;

        extern RDL_DLLAPI SpatialVector SpatialVectorZero;
        extern RDL_DLLAPI SpatialMatrix SpatialMatrixIdentity;
        extern RDL_DLLAPI SpatialMatrix SpatialMatrixZero;

        RDL_DLLAPI inline VectorNd vectorFromPtr(unsigned int n, double *ptr)
        {
            // TODO: use memory mapping operators for Eigen
            VectorNd result(n);

            for (unsigned int i = 0; i < n; i++)
            {
                result[i] = ptr[i];
            }

            return result;
        }

        RDL_DLLAPI inline MatrixNd matrixFromPtr(unsigned int rows, unsigned int cols, double *ptr, bool row_major = true)
        {
            MatrixNd result(rows, cols);

            if (row_major)
            {
                for (unsigned int i = 0; i < rows; i++)
                {
                    for (unsigned int j = 0; j < cols; j++)
                    {
                        result(i, j) = ptr[i * cols + j];
                    }
                }
            }
            else
            {
                for (unsigned int i = 0; i < rows; i++)
                {
                    for (unsigned int j = 0; j < cols; j++)
                    {
                        result(i, j) = ptr[i + j * rows];
                    }
                }
            }

            return result;
        }

        /// \brief Solves a linear system using gaussian elimination with pivoting
        RDL_DLLAPI bool            linSolveGaussElimPivot(MatrixNd A, VectorNd b, VectorNd& x);

        /** \brief Translates the inertia matrix to a new center. */
        RDL_DLLAPI Matrix3d        parallel_axis(const Matrix3d& inertia, double mass, const Vector3d& com);

        /** \brief Creates a transformation of a linear displacement
         *
         * This can be used to specify the translation to the joint center when
         * adding a body to a model. See also section 2.8 in RBDA.
         *
         * \note The transformation returned is for motions. For a transformation for forces
         * \note one has to conjugate the matrix.
         *
         * \param displacement The displacement as a 3D vector
         */
        RDL_DLLAPI SpatialMatrix   Xtrans_mat(const Vector3d& displacement);

        /** \brief Creates a rotational transformation around the Z-axis
         *
         * Creates a rotation around the current Z-axis by the given angle
         * (specified in radians).
         *
         * \param zrot Rotation angle in radians.
         */
        RDL_DLLAPI SpatialMatrix   Xrotz_mat(const double& zrot);

        /** \brief Creates a rotational transformation around the Y-axis
         *
         * Creates a rotation around the current Y-axis by the given angle
         * (specified in radians).
         *
         * \param yrot Rotation angle in radians.
         */
        RDL_DLLAPI SpatialMatrix   Xroty_mat(const double& yrot);

        /** \brief Creates a rotational transformation around the X-axis
         *
         * Creates a rotation around the current X-axis by the given angle
         * (specified in radians).
         *
         * \param xrot Rotation angle in radians.
         */
        RDL_DLLAPI SpatialMatrix   Xrotx_mat(const double& xrot);

        /** \brief Creates a spatial transformation for given parameters
         *
         * Creates a transformation to a coordinate system that is first rotated
         * and then translated.
         *
         * \param displacement The displacement to the new origin
         * \param zyx_euler The orientation of the new coordinate system, specifyed
         * by ZYX-Euler angles.
         */
        RDL_DLLAPI SpatialMatrix   XtransRotZYXEuler(const Vector3d& displacement, const Vector3d& zyx_euler);

        RDL_DLLAPI inline Matrix3d rotx(const double& xrot)
        {
            double s, c;

            s = sin(xrot);
            c = cos(xrot);
            return Matrix3d(
                       1., 0., 0.,
                       0., c, s,
                       0., -s, c
                       );
        }

        RDL_DLLAPI inline Matrix3d roty(const double& yrot)
        {
            double s, c;

            s = sin(yrot);
            c = cos(yrot);
            return Matrix3d(
                       c, 0., -s,
                       0., 1., 0.,
                       s, 0., c
                       );
        }

        RDL_DLLAPI inline Matrix3d rotz(const double& zrot)
        {
            double s, c;

            s = sin(zrot);
            c = cos(zrot);
            return Matrix3d(
                       c, s, 0.,
                       -s, c, 0.,
                       0., 0., 1.
                       );
        }

        RDL_DLLAPI inline Matrix3d rotxdot(const double& x, const double& xdot)
        {
            double s, c;

            s = sin(x);
            c = cos(x);
            return Matrix3d(
                       0., 0., 0.,
                       0., -s * xdot, c * xdot,
                       0., -c * xdot, -s * xdot
                       );
        }

        RDL_DLLAPI inline Matrix3d rotydot(const double& y, const double& ydot)
        {
            double s, c;

            s = sin(y);
            c = cos(y);
            return Matrix3d(
                       -s * ydot, 0., -c * ydot,
                       0., 0., 0.,
                       c * ydot, 0., -s * ydot
                       );
        }

        RDL_DLLAPI inline Matrix3d rotzdot(const double& z, const double& zdot)
        {
            double s, c;

            s = sin(z);
            c = cos(z);
            return Matrix3d(
                       -s * zdot, c * zdot, 0.,
                       -c * zdot, -s * zdot, 0.,
                       0., 0., 0.
                       );
        }

        RDL_DLLAPI inline Vector3d angular_velocity_from_angle_rates(const Vector3d& zyx_angles,
                                                                     const Vector3d& zyx_angle_rates)
        {
            double sy = sin(zyx_angles[1]);
            double cy = cos(zyx_angles[1]);
            double sx = sin(zyx_angles[2]);
            double cx = cos(zyx_angles[2]);

            return Vector3d(
                       zyx_angle_rates[2] - sy * zyx_angle_rates[0],
                       cx * zyx_angle_rates[1] + sx * cy * zyx_angle_rates[0],
                       -sx * zyx_angle_rates[1] + cx * cy * zyx_angle_rates[0]
                       );
        }

        RDL_DLLAPI inline Vector3d global_angular_velocity_from_rates(const Vector3d& zyx_angles, const Vector3d& zyx_rates)
        {
            Matrix3d RzT = rotz(zyx_angles[0]).transpose();
            Matrix3d RyT = roty(zyx_angles[1]).transpose();

            return Vector3d(
                       Vector3d(0., 0., zyx_rates[0])
                       + RzT * Vector3d(0., zyx_rates[1], 0.)
                       + RzT * RyT * Vector3d(zyx_rates[2], 0., 0.)
                       );
        }

        RDL_DLLAPI inline Vector3d angular_acceleration_from_angle_rates(const Vector3d& zyx_angles,
                                                                         const Vector3d& zyx_angle_rates,
                                                                         const Vector3d& zyx_angle_rates_dot)
        {
            double sy    = sin(zyx_angles[1]);
            double cy    = cos(zyx_angles[1]);
            double sx    = sin(zyx_angles[2]);
            double cx    = cos(zyx_angles[2]);
            double xdot  = zyx_angle_rates[2];
            double ydot  = zyx_angle_rates[1];
            double zdot  = zyx_angle_rates[0];
            double xddot = zyx_angle_rates_dot[2];
            double yddot = zyx_angle_rates_dot[1];
            double zddot = zyx_angle_rates_dot[0];

            return Vector3d(
                       xddot - (cy * ydot * zdot + sy * zddot),
                       -sx * xdot * ydot + cx * yddot + cx * xdot * cy * zdot + sx * (-sy * ydot * zdot + cy * zddot),
                       -cx * xdot * ydot - sx * yddot - sx * xdot * cy * zdot + cx * (-sy * ydot * zdot + cy * zddot)
                       );
        }

        RDL_DLLAPI
        void SparseFactorizeLTL(Model& model, Math::MatrixNd& H);

        RDL_DLLAPI
        void SparseMultiplyHx(Model& model, Math::MatrixNd& L);

        RDL_DLLAPI
        void SparseMultiplyLx(Model& model, Math::MatrixNd& L);

        RDL_DLLAPI
        void SparseMultiplyLTx(Model& model, Math::MatrixNd& L);

        RDL_DLLAPI
        void SparseSolveLx(Model& model, Math::MatrixNd& L, Math::VectorNd& x);

        RDL_DLLAPI
        void SparseSolveLTx(Model& model, Math::MatrixNd& L, Math::VectorNd& x);
    } /* Math */
}     /* RobotDynamics */

/* __RDL_MATHUTILS_H__ */
#endif // ifndef __RDL_MATHUTILS_H__
