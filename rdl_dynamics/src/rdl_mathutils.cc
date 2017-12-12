/*
 * Original Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 *
 * RDL - Robot Dynamics Library
 * Modifications Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <cmath>
#include <limits>

#include <iostream>
#include <assert.h>

#include "rdl_dynamics/rdl_mathutils.h"
#include <rdl_dynamics/Model.h>

namespace RobotDynamics
{
    namespace Math
    {

        RDL_DLLAPI Vector3d Vector3dZero(0., 0., 0.);
        RDL_DLLAPI Matrix3d Matrix3dIdentity(1., 0., 0., 0., 1., 0., 0., 0., 1);
        RDL_DLLAPI Matrix3d Matrix3dZero(0., 0., 0., 0., 0., 0., 0., 0., 0.);

        RDL_DLLAPI SpatialVector SpatialVectorZero(0., 0., 0., 0., 0., 0.);

        RDL_DLLAPI SpatialMatrix SpatialMatrixIdentity(1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 1.);

        RDL_DLLAPI SpatialMatrix SpatialMatrixZero(0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);

        RDL_DLLAPI bool linSolveGaussElimPivot(MatrixNd A, VectorNd b, VectorNd &x)
        {
            x.setZero();

            // We can only solve quadratic systems
            assert (A.rows() == A.cols());

            unsigned int n = A.rows();
            // cppcheck-suppress variableScope
            unsigned int pi;

            // the pivots
            size_t *pivot = new size_t[n];

            // temporary result vector which contains the pivoted result
            VectorNd px(x);

            unsigned int i, j, k;

            for(i = 0; i < n; i++)
            {
                pivot[i] = i;
            }

            for(j = 0; j < n; j++)
            {
                pi = j;
                double pv = fabs(A(j, pivot[j]));

                // find the pivot
                for(k = j; k < n; k++)
                {
                    double pt = fabs(A(j, pivot[k]));
                    if(pt > pv)
                    {
                        pv = pt;
                        pi = k;
                        unsigned int p_swap = pivot[j];
                        pivot[j] = pivot[pi];
                        pivot[pi] = p_swap;
                    }
                }

                for(i = j + 1; i < n; i++)
                {
                    if(fabs(A(j, pivot[j])) <= std::numeric_limits<double>::epsilon())
                    {
                        std::cerr << "Error: pivoting failed for matrix A = " << std::endl;
                        std::cerr << "A = " << std::endl << A << std::endl;
                        std::cerr << "b = " << b << std::endl;
                    }
                    //		assert (fabs(A(j,pivot[j])) > std::numeric_limits<double>::epsilon());
                    double d = A(i, pivot[j]) / A(j, pivot[j]);

                    b[i] -= b[j] * d;

                    for(k = j; k < n; k++)
                    {
                        A(i, pivot[k]) -= A(j, pivot[k]) * d;
                    }
                }
            }

            // warning: i is an unsigned int, therefore a for loop of the
            // form "for (i = n - 1; i >= 0; i--)" might end up in getting an invalid
            // value for i!
            i = n;
            do
            {
                i--;

                for(j = i + 1; j < n; j++)
                {
                    px[i] += A(i, pivot[j]) * px[j];
                }
                px[i] = (b[i] - px[i]) / A(i, pivot[i]);

            } while(i > 0);

            // Unswapping
            for(i = 0; i < n; i++)
            {
                x[pivot[i]] = px[i];
            }

            delete[] pivot;

            return true;
        }

        RDL_DLLAPI Matrix3d parallel_axis(const Matrix3d &inertia, double mass, const Vector3d &com)
        {
            Matrix3d com_cross = toTildeForm(com);

            return inertia + mass * com_cross * com_cross.transpose();
        }

        RDL_DLLAPI SpatialMatrix Xtrans_mat(const Vector3d &r)
        {
            return SpatialMatrix(1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., r[2], -r[1], 1., 0., 0., -r[2], 0., r[0], 0., 1., 0., r[1], -r[0], 0., 0., 0., 1.);
        }

        RDL_DLLAPI SpatialMatrix Xrotx_mat(const double &xrot)
        {
            double s, c;
            s = sin(xrot);
            c = cos(xrot);

            return SpatialMatrix(1., 0., 0., 0., 0., 0., 0., c, s, 0., 0., 0., 0., -s, c, 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., c, s, 0., 0., 0., 0., -s, c);
        }

        RDL_DLLAPI SpatialMatrix Xroty_mat(const double &yrot)
        {
            double s, c;
            s = sin(yrot);
            c = cos(yrot);

            return SpatialMatrix(c, 0., -s, 0., 0., 0., 0., 1., 0., 0., 0., 0., s, 0., c, 0., 0., 0., 0., 0., 0., c, 0., -s, 0., 0., 0., 0., 1., 0., 0., 0., 0., s, 0., c);
        }

        RDL_DLLAPI SpatialMatrix Xrotz_mat(const double &zrot)
        {
            double s, c;
            s = sin(zrot);
            c = cos(zrot);

            return SpatialMatrix(c, s, 0., 0., 0., 0., -s, c, 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., c, s, 0., 0., 0., 0., -s, c, 0., 0., 0., 0., 0., 0., 1.);
        }

        RDL_DLLAPI SpatialMatrix XtransRotZYXEuler(const Vector3d &displacement, const Vector3d &zyx_euler)
        {
            return Xrotz_mat(zyx_euler[0]) * Xroty_mat(zyx_euler[1]) * Xrotx_mat(zyx_euler[2]) * Xtrans_mat(displacement);
        }

        RDL_DLLAPI void SparseFactorizeLTL(Model &model, Math::MatrixNd &H)
        {
            for(unsigned int i = 0; i < model.qdot_size; i++)
            {
                for(unsigned int j = i + 1; j < model.qdot_size; j++)
                {
                    H(i, j) = 0.;
                }
            }

            for(unsigned int k = model.qdot_size; k > 0; k--)
            {
                H(k - 1, k - 1) = sqrt(H(k - 1, k - 1));
                unsigned int i = model.lambda_q[k];
                while(i != 0)
                {
                    H(k - 1, i - 1) = H(k - 1, i - 1) / H(k - 1, k - 1);
                    i = model.lambda_q[i];
                }

                i = model.lambda_q[k];
                while(i != 0)
                {
                    unsigned int j = i;
                    while(j != 0)
                    {
                        H(i - 1, j - 1) = H(i - 1, j - 1) - H(k - 1, i - 1) * H(k - 1, j - 1);
                        j = model.lambda_q[j];
                    }
                    i = model.lambda_q[i];
                }
            }
        }

        RDL_DLLAPI void SparseMultiplyHx(Model &model, Math::MatrixNd &L)
        {
            assert (0 && !"Not yet implemented!");
        }

        RDL_DLLAPI void SparseMultiplyLx(Model &model, Math::MatrixNd &L)
        {
            assert (0 && !"Not yet implemented!");
        }

        RDL_DLLAPI void SparseMultiplyLTx(Model &model, Math::MatrixNd &L)
        {
            assert (0 && !"Not yet implemented!");
        }

        RDL_DLLAPI void SparseSolveLx(Model &model, Math::MatrixNd &L, Math::VectorNd &x)
        {
            for(unsigned int i = 1; i <= model.qdot_size; i++)
            {
                unsigned int j = model.lambda_q[i];
                while(j != 0)
                {
                    x[i - 1] = x[i - 1] - L(i - 1, j - 1) * x[j - 1];
                    j = model.lambda_q[j];
                }
                x[i - 1] = x[i - 1] / L(i - 1, i - 1);
            }
        }

        RDL_DLLAPI void SparseSolveLTx(Model &model, Math::MatrixNd &L, Math::VectorNd &x)
        {
            for(int i = model.qdot_size; i > 0; i--)
            {
                x[i - 1] = x[i - 1] / L(i - 1, i - 1);
                unsigned int j = model.lambda_q[i];
                while(j != 0)
                {
                    x[j - 1] = x[j - 1] - L(i - 1, j - 1) * x[i - 1];
                    j = model.lambda_q[j];
                }
            }
        }
    } /* Math */
} /* RobotDynamics */
