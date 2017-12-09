/*
 * RDL - Robot Dynamics Library
 * Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef __RDL_MOMENTUM_HPP__
#define __RDL_MOMENTUM_HPP__

/**
 * @file Momentum.hpp
 */

#include "rdl_dynamics/RigidBodyInertia.hpp"
#include "rdl_dynamics/TransformableGeometricObject.hpp"

namespace RobotDynamics
{
    namespace Math
    {
        /**
         * @class Momentum
         * @brief Momentum is mass/inertia multiplied by velocity
         */
        class Momentum : public ForceVector
        {
public:

            Momentum() : ForceVector()
            {}

            Momentum(const double kx, const double ky, const double kz, const double lx, const double ly,
                     const double lz) : ForceVector(kx, ky, kz, lx, ly, lz)
            {}

            Momentum(const Vector3d& k, const Vector3d l) : ForceVector(k.x(), k.y(), k.z(), l.x(), l.y(), l.z())
            {}

            Momentum(const Momentum& momentum) : ForceVector(momentum)
            {}

            Momentum(const ForceVector& forceVector) : ForceVector(forceVector)
            {}

            Momentum(const RigidBodyInertia& inertia, const MotionVector& vector) : ForceVector()
            {
                computeMomentum(inertia, vector);
            }

            /**
             * @brief Set momentum by computing it from the inertia and velocity
             * @param inertia
             * @param vector
             */
            inline void set(const RigidBodyInertia& inertia, const MotionVector& vector)
            {
                computeMomentum(inertia, vector);
            }

            EIGEN_STRONG_INLINE double& kx()
            {
                return this->data()[0];
            }

            EIGEN_STRONG_INLINE double& ky()
            {
                return this->data()[1];
            }

            EIGEN_STRONG_INLINE double& kz()
            {
                return this->data()[2];
            }

            EIGEN_STRONG_INLINE double kx() const
            {
                return this->data()[0];
            }

            EIGEN_STRONG_INLINE double ky() const
            {
                return this->data()[1];
            }

            EIGEN_STRONG_INLINE double kz() const
            {
                return this->data()[2];
            }

            EIGEN_STRONG_INLINE double& lx()
            {
                return this->data()[3];
            }

            EIGEN_STRONG_INLINE double& ly()
            {
                return this->data()[4];
            }

            EIGEN_STRONG_INLINE double& lz()
            {
                return this->data()[5];
            }

            EIGEN_STRONG_INLINE double lx() const
            {
                return this->data()[3];
            }

            EIGEN_STRONG_INLINE double ly() const
            {
                return this->data()[4];
            }

            EIGEN_STRONG_INLINE double lz() const
            {
                return this->data()[5];
            }

            /**
             * @brief Operator for computing kinetic energy. With momentum, \f$ m \f$ and Math::MotionVector, \f$ v \f$ this
             **performs performs \f$ \frac{m\cdot v}{2} \f$
             * @param vector
             * @return kinetic energy
             */
            EIGEN_STRONG_INLINE double operator*(const MotionVector& vector)
            {
                return this->dot(vector) * 0.5;
            }

protected:

            /**
             * @brief Computes momentum from inertia and velocity
             * @param I
             * @param v
             */
            inline void computeMomentum(const RigidBodyInertia& I, const MotionVector& v)
            {
                (*this) << I.Ixx * v[0] + I.Iyx * v[1] + I.Izx * v[2] + I.h[1] * v[5] - I.h[2] * v[4],
                I.Iyx * v[0] + I.Iyy * v[1] + I.Izy * v[2] - I.h[0] * v[5] + I.h[2] * v[3],
                I.Izx * v[0] + I.Izy * v[1] + I.Izz * v[2] + I.h[0] * v[4] - I.h[1] * v[3],
                -I.h[1] * v[2] + I.h[2] * v[1] + I.m * v[3],
                I.h[0] * v[2] - I.h[2] * v[0] + I.m * v[4],
                -I.h[0] * v[1] + I.h[1] * v[0] + I.m * v[5];
            }
        };

        /**
         * @brief Multiplying a Math::RigidBodyInertia by a Math::MotionVector returns a Momentum
         * @param I
         * @param v
         * @return
         */
        inline Momentum operator*(const RigidBodyInertia& I, const MotionVector& v)
        {
            return Momentum(I, v);
        }
    }
}

#endif // ifndef __RDL_MOMENTUM_HPP__
