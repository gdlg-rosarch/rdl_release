/*
 * RDL - Robot Dynamics Library
 * Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include "rdl_dynamics/RigidBodyInertia.hpp"

namespace RobotDynamics
{
    namespace Math
    {
        void RigidBodyInertia::operator+=(const RigidBodyInertia &rbi)
        {
            this->m+=rbi.m;
            this->h+=rbi.h;
            this->Ixx+=rbi.Ixx;
            this->Iyx+=rbi.Iyx;
            this->Iyy+=rbi.Iyy;
            this->Izx+=rbi.Izx;
            this->Izy+=rbi.Izy;
            this->Izz+=rbi.Izz;
        }

        void RigidBodyInertia::transform_transpose_slow(const SpatialTransform &X)
        {
            Matrix3d E = X.E;
            Vector3d r = X.r;

            Vector3d E_T_mr = E.transpose() * h + m * r;

            Matrix3d I = E.transpose() * Matrix3d(Ixx, Iyx, Izx, Iyx, Iyy, Izy, Izx, Izy, Izz) * E - toTildeForm(r) * toTildeForm(E.transpose() * h) - toTildeForm(E_T_mr) * toTildeForm(r);

            // Apparently this assign chunk takes forever.
            Ixx = I(0,0);
            Iyx = I(0,1);
            Izx = I(0,2);
            Iyy = I(1,1);
            Izy = I(1,2);
            Izz = I(2,2);

            h = E_T_mr;
        }

        /** Same as X^T I X
         * This method is optimized for speed. Indexing is super expensive
        */
        RigidBodyInertia RigidBodyInertia::transform_transpose_copy(const SpatialTransform &X) const
        {
            double E00 = X.E(0,0);
            double E01 = X.E(0,1);
            double E02 = X.E(0,2);
            double E10 = X.E(1,0);
            double E11 = X.E(1,1);
            double E12 = X.E(1,2);
            double E20 = X.E(2,0);
            double E21 = X.E(2,1);
            double E22 = X.E(2,2);

            double h0 = h(0);
            double h1 = h(1);
            double h2 = h(2);

            double r0 = X.r(0);
            double r1 = X.r(1);
            double r2 = X.r(2);

            // E_T_h -> alpha
            double alpha_x = E00*h0 + E10*h1 + E20*h2;
            double alpha_y = E01*h0 + E11*h1 + E21*h2;
            double alpha_z = E02*h0 + E12*h1 + E22*h2;

            double d = r2*alpha_z;
            double e = r1*alpha_y;
            double f = r0*alpha_x;

            double pi00 = E00*Ixx + E10*Iyx + E20*Izx;
            double pi01 = E00*Iyx + E10*Iyy + E20*Izy;
            double pi02 = E00*Izx + E10*Izy + E20*Izz;
            double pi10 = E01*Ixx + E11*Iyx + E21*Izx;
            double pi11 = E01*Iyx + E11*Iyy + E21*Izy;
            double pi12 = E01*Izx + E11*Izy + E21*Izz;
            double pi20 = E02*Ixx + E12*Iyx + E22*Izx;
            double pi21 = E02*Iyx + E12*Iyy + E22*Izy;
            double pi22 = E02*Izx + E12*Izy + E22*Izz;

            double eps_x = alpha_x + m*r0;
            double eps_y = alpha_y + m*r1;
            double eps_z = alpha_z + m*r2;

            double a = eps_z*r2;
            double b = eps_y*r1;
            double c = eps_x*r0;

            double phi00 = -a - b;
            double phi01 = eps_y*r0;
            double phi02 = eps_z*r0;
            double phi11 = -a - c;
            double phi12 = eps_z*r1;
            double phi22 = -b - c;

            return RigidBodyInertia(m,Vector3d(eps_x,eps_y,eps_z),
                pi00*E00 + pi01*E10 + pi02*E20 - (-d - e) - phi00, 
                pi00*E01 + pi01*E11 + pi02*E21 - (r1*alpha_x) - phi01,
                pi10*E01 + pi11*E11 + pi12*E21 - (-d - f) - phi11,
                pi00*E02 + pi01*E12 + pi02*E22 - (r2*alpha_x) - phi02,
                pi10*E02 + pi11*E12 + pi12*E22 - (r2*alpha_y) - phi12,
                pi20*E02 + pi21*E12 + pi22*E22 - (-e - f) - phi22);
        }

        RigidBodyInertia RigidBodyInertia::transform_copy(const SpatialTransform &X) const
        {
            double hx = h[0];
            double hy = h[1];
            double hz = h[2];

            double rx = X.r[0];
            double ry = X.r[1];
            double rz = X.r[2];

            double E00 = X.E(0,0);
            double E01 = X.E(0,1);
            double E02 = X.E(0,2);
            double E10 = X.E(1,0);
            double E11 = X.E(1,1);
            double E12 = X.E(1,2);
            double E20 = X.E(2,0);
            double E21 = X.E(2,1);
            double E22 = X.E(2,2);

            double alpha = 2.*hx*rx;
            double beta = 2.*hy*ry;
            double gamma = 2.*hz*rz;
            double epsilon = m*rx*rx;
            double phi = m*ry*ry;
            double psi = m*rz*rz;
            double pi = m*rx*ry;
            double a = m*rx*rz;
            double b = m*ry*rz;
            double c = hx*ry;
            double d = hx*rz;
            double e = hy*rz;
            double f = hy*rx;
            double g = hz*rx;
            double j = hz*ry;

            double b00 = Ixx - gamma - beta + psi + phi;
            double b01 = Iyx + f + c - pi;
            double b02 = Izx + g + d - a;
            double b11 = Iyy - gamma - alpha + psi + epsilon;
            double b12 = Izy + j + e - b;
            double b22 = Izz - beta - alpha + phi + epsilon;

            double w00 = E00*b00 + E01*b01 + E02*b02;
            double w01 = E00*b01 + E01*b11 + E02*b12;
            double w02 = E00*b02 + E01*b12 + E02*b22;

            double w10 = E10*b00 + E11*b01 + E12*b02;
            double w11 = E10*b01 + E11*b11 + E12*b12;
            double w12 = E10*b02 + E11*b12 + E12*b22;

            double w20 = E20*b00 + E21*b01 + E22*b02;
            double w21 = E20*b01 + E21*b11 + E22*b12;
            double w22 = E20*b02 + E21*b12 + E22*b22;

            return RigidBodyInertia(m,X.E * (h - m * X.r),
                w00*E00 + w01*E01 + w02*E02,
                w00*E10 + w01*E11 + w02*E12,
                w10*E10 + w11*E11 + w12*E12,
                w00*E20 + w01*E21 + w02*E22,
                w10*E20 + w11*E21 + w12*E22,
                w20*E20 + w21*E21 + w22*E22);
        }

        void RigidBodyInertia::transform_slow(const SpatialTransform &X)
        {
            Vector3d newH = X.E * (h - m * X.r);
            Matrix3d newI = X.E * (Matrix3d(Ixx, Iyx, Izx, Iyx, Iyy, Izy, Izx, Izy, Izz) + toTildeForm(X.r) * toTildeForm(h) + (toTildeForm(h - m * X.r) * toTildeForm(X.r))) * X.E.transpose();
            h = newH;
            Ixx = newI(0,0);
            Iyx = newI(0,1);
            Izx = newI(0,2);
            Iyy = newI(1,1);
            Izy = newI(1,2);
            Izz = newI(2,2);
        }

        void RigidBodyInertia::createFromMatrix(const SpatialMatrix &Ic)
        {
            m = Ic(3, 3);
            h.set(-Ic(1, 5), Ic(0, 5), -Ic(0, 4));
            Ixx = Ic(0, 0);
            Iyx = Ic(1, 0);
            Iyy = Ic(1, 1);
            Izx = Ic(2, 0);
            Izy = Ic(2, 1);
            Izz = Ic(2, 2);
        }

        SpatialMatrix RigidBodyInertia::toMatrix() const
        {
            SpatialMatrix result;
            result(0, 0) = Ixx;
            result(0, 1) = Iyx;
            result(0, 2) = Izx;
            result(1, 0) = Iyx;
            result(1, 1) = Iyy;
            result(1, 2) = Izy;
            result(2, 0) = Izx;
            result(2, 1) = Izy;
            result(2, 2) = Izz;

            result.block<3, 3>(0, 3) = toTildeForm(h);
            result.block<3, 3>(3, 0) = -toTildeForm(h);
            result.block<3, 3>(3, 3) = Matrix3d::Identity(3, 3) * m;

            return result;
        }

        void RigidBodyInertia::setSpatialMatrix(SpatialMatrix &mat) const
        {
            mat(0, 0) = Ixx;
            mat(0, 1) = Iyx;
            mat(0, 2) = Izx;
            mat(1, 0) = Iyx;
            mat(1, 1) = Iyy;
            mat(1, 2) = Izy;
            mat(2, 0) = Izx;
            mat(2, 1) = Izy;
            mat(2, 2) = Izz;

            mat(3, 0) = 0.;
            mat(3, 1) = h[2];
            mat(3, 2) = -h[1];
            mat(4, 0) = -h[2];
            mat(4, 1) = 0.;
            mat(4, 2) = h[0];
            mat(5, 0) = h[1];
            mat(5, 1) = -h[0];
            mat(5, 2) = 0.;

            mat(0, 3) = 0.;
            mat(0, 4) = -h[2];
            mat(0, 5) = h[1];
            mat(1, 3) = h[2];
            mat(1, 4) = 0.;
            mat(1, 5) = -h[0];
            mat(2, 3) = -h[1];
            mat(2, 4) = h[0];
            mat(2, 5) = 0.;

            mat(3, 3) = m;
            mat(3, 4) = 0.;
            mat(3, 5) = 0.;
            mat(4, 3) = 0.;
            mat(4, 4) = m;
            mat(4, 5) = 0.;
            mat(5, 3) = 0.;
            mat(5, 4) = 0.;
            mat(5, 5) = m;
        }
    }
}