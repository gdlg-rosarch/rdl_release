/*
 * RDL - Robot Dynamics Library
 * Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef __RDL_TRANSFORMABLE_GEOMETRIC_OBJECT_HPP__
#define __RDL_TRANSFORMABLE_GEOMETRIC_OBJECT_HPP__

/**
 * @file TransformableGeometricObject.hpp
 * @brief Contains various geometric objects that have methods for transforming themselves into different frames of reference
 */

#include "rdl_dynamics/rdl_eigenmath.h"
#include "rdl_dynamics/SpatialAlgebraOperators.h"

namespace RobotDynamics
{
    namespace Math
    {
        /**
         * @class TransformableGeometricObject
         * @brief The TransformableGeometricObject class is an essential interface because it forces all geometric objects
         * to implement a method that tells how to transform them. This makes in possible for frame transformations of any
         * TransformableGeometricObject can be done via the FrameObject::changeFrame method.
         */
        class RDL_DLLAPI TransformableGeometricObject
        {
public:

            /**
             * @brief Pure virtual object. This object forces objects that inherit from it to have a method that tells
             * how that object is transformed.
             * @param X SpatialTransform
             */
            virtual void transform(const RobotDynamics::Math::SpatialTransform& X) = 0;
        };

        /**
         * @class Vector3d_X
         * @brief This object is here because a Vector3d type needs to implement the transform method so a frames::FrameVector
         **can
         * inherit from it, and since Math::SpatialTransform already depends on Vector3d, a new class had to be created that
         **could
         * implement the transform mehod AND depend on Math::SpatialTransform.
         */
        class RDL_DLLAPI Vector3d_X : public Vector3d, public TransformableGeometricObject
        {
public:

            /**
             * @brief Constructor
             * @param other
             */
            template<typename OtherDerived>
            Vector3d_X(const Eigen::MatrixBase<OtherDerived>& other)
                : Vector3d(other)
            {}

            /**
             * @brief Overloaded equal operator
             * @param other
             * @return Copied Vector3d_X
             */
            template<typename OtherDerived>
            Vector3d_X& operator=(const Eigen::MatrixBase<OtherDerived>& other)
            {
                Vector3d::operator=(other);
                return *this;
            }

            /**
             * @brief Empty constructor
             */
            EIGEN_STRONG_INLINE Vector3d_X()
            {}

            /**
             * @brief Constructor
             * @param v0 x
             * @param v1 y
             * @param v2 z
             */
            EIGEN_STRONG_INLINE Vector3d_X(const double v0, const double v1, const double v2)
            {
                Base::_check_template_params();

                (*this) << v0, v1, v2;
            }

            /**
             * @brief Rotates a 3D vector by the given SpatialTransform. Performs \f$ v = E v \f$
             * @param X SpatialTransform used to rotate the 3D vector
             */
            inline void transform(const RobotDynamics::Math::SpatialTransform& X)
            {
                this->set(this->transform_copy(X));
            }

            void operator+=(const Vector3d_X& v)
            {
                this->x() += v.x();
                this->y() += v.y();
                this->z() += v.z();
            }

            void operator-=(const Vector3d_X& v)
            {
                this->x() -= v.x();
                this->y() -= v.y();
                this->z() -= v.z();
            }

            void operator*=(double scale)
            {
                this->x() *= scale;
                this->y() *= scale;
                this->z() *= scale;
            }

            /**
             * @brief Rotates a 3D vector by the given SpatialTransform and returns the result in a new Vector3d_X. Performs \f$
             **v_2 = E v_1 \f$
             * @param X SpatialTransform used to rotate the 3D vector
             * @return Returns a copied, transformed 3D vector
             */
            inline Vector3d_X transform_copy(const RobotDynamics::Math::SpatialTransform& X)
            {
                return X.E * (*this);
            }
        };

        inline Vector3d_X operator+(Vector3d_X v1, const Vector3d_X& v2)
        {
            v1 += v2;
            return v1;
        }

        inline Vector3d_X operator-(Vector3d_X v1, const Vector3d_X& v2)
        {
            v1 -= v2;
            return v1;
        }

        inline Vector3d_X operator*(Vector3d_X v, double scale)
        {
            v *= scale;
            return v;
        }

        /**
         * @brief Operator that performs a 3D vector transform. Does not alter the Vector3d_X argument.
         */
        EIGEN_STRONG_INLINE Vector3d_X operator*(const SpatialTransform& X, Vector3d_X v)
        {
            v.transform(X);
            return v;
        }

        /**
         * @class ForceVector
         * @brief A ForceVector is a SpatialVector containing 3 moments and 3 linear forces.
         */
        class RDL_DLLAPI ForceVector : public SpatialVector, public TransformableGeometricObject
        {
public:

            /**
             * @brief Constructor
             * @param other
             */
            template<typename OtherDerived>
            ForceVector(const Eigen::MatrixBase<OtherDerived>& other)
                : SpatialVector(other)
            {}

            /**
             * Overloaded equal operator
             * @param other
             * @return Copied ForceVector
             */
            template<typename OtherDerived>
            ForceVector& operator=(const Eigen::MatrixBase<OtherDerived>& other)
            {
                SpatialVector::operator=(other);
                return *this;
            }

            /**
             * @brief Empty constructor
             */
            EIGEN_STRONG_INLINE ForceVector()
            {}

            /**
             * @brief Get a copy of a ForceVector as type SpatialVector
             */
            EIGEN_STRONG_INLINE SpatialVector toSpatialVector() const
            {
                return *this;
            }

            /**
             * @brief Constructor
             * @param mx x-Moment
             * @param my y-Moment
             * @param mz z-Moment
             * @param fx x-Linear Force
             * @param fy y-Linear Force
             * @param fz z-Linear Force
             */
            ForceVector(const double mx, const double my, const double mz, const double fx, const double fy, const double fz)
            {
                Base::_check_template_params();

                (*this) << mx, my, mz, fx, fy, fz;
            }

            /**
             * @brief Setter
             * @param f
             */
            EIGEN_STRONG_INLINE void set(const ForceVector& f)
            {
                (*this) << f.data()[0], f.data()[1], f.data()[2], f.data()[3], f.data()[4], f.data()[5];
            }

            /**
             * @brief Get reference to x-angular component
             * @return Reference to x-angular component
             */
            EIGEN_STRONG_INLINE double& mx()
            {
                return this->data()[0];
            }

            /**
             * @brief Get reference to y-angular component
             * @return Reference to y-angular component
             */
            EIGEN_STRONG_INLINE double& my()
            {
                return this->data()[1];
            }

            /**
             * @brief Get reference to z-angular component
             * @return Reference to z-angular component
             */
            EIGEN_STRONG_INLINE double& mz()
            {
                return this->data()[2];
            }

            /**
             * @brief Get copy of x-angular component
             * @return Copy of x-angular component
             */
            EIGEN_STRONG_INLINE double mx() const
            {
                return this->data()[0];
            }

            /**
             * @brief Get copy of y-angular component
             * @return Copy of y-angular component
             */
            EIGEN_STRONG_INLINE double my() const
            {
                return this->data()[1];
            }

            /**
             * @brief Get copy of z-angular component
             * @return Copy of z-angular component
             */
            EIGEN_STRONG_INLINE double mz() const
            {
                return this->data()[2];
            }

            /**
             * @brief Get reference to x-linear component
             * @return Reference to x-linear component
             */
            EIGEN_STRONG_INLINE double& fx()
            {
                return this->data()[3];
            }

            /**
             * @brief Get reference to y-linear component
             * @return Reference to y-linear component
             */
            EIGEN_STRONG_INLINE double& fy()
            {
                return this->data()[4];
            }

            /**
             * @brief Get reference to z-linear component
             * @return Reference to z-linear component
             */
            EIGEN_STRONG_INLINE double& fz()
            {
                return this->data()[5];
            }

            /**
             * @brief Get copy of x-linear component
             * @return Copy of x-linear component
             */
            EIGEN_STRONG_INLINE double fx() const
            {
                return this->data()[3];
            }

            /**
             * @brief Get copy of y-linear component
             * @return Copy of y-linear component
             */
            EIGEN_STRONG_INLINE double fy() const
            {
                return this->data()[4];
            }

            /**
             * @brief Get copy of z-linear component
             * @return Copy of z-linear component
             */
            EIGEN_STRONG_INLINE double fz() const
            {
                return this->data()[5];
            }

            /**
             * @brief Transform a force vector. Performs \f$ f= X^T*f \f$
             * @param X SpatialTransform
             */
            inline void transformTranspose(const SpatialTransform& X)
            {
                this->set(this->transformTranspose_copy(X));
            }

            /**
             * @brief Copies then transforms a ForceVector
             * @param X
             * @return Returns a copied, transform ForceVector
             */
            ForceVector transformTranspose_copy(const SpatialTransform& X) const
            {
                Vector3d E_T_f(X.E(0, 0) * this->data()[3] + X.E(1, 0) * this->data()[4] + X.E(2, 0) * this->data()[5], X.E(0,
                                                                                                                            1)
                               * this->data()[3] + X.E(1, 1) * this->data()[4] + X.E(2, 1) * this->data()[5], X.E(0,
                                                                                                                  2)
                               * this->data()[3] + X.E(1, 2) * this->data()[4] + X.E(2, 2) * this->data()[5]);

                return ForceVector(X.E(0, 0) * this->data()[0] + X.E(1, 0) * this->data()[1] + X.E(2,
                                                                                                   0) * this->data()[2] -
                                   X.r[2] * E_T_f[1] + X.r[1] * E_T_f[2], X.E(0, 1) * this->data()[0] + X.E(1,
                                                                                                            1)
                                   * this->data()[1] +
                                   X.E(2, 1) * this->data()[2] + X.r[2] * E_T_f[0] - X.r[0] * E_T_f[2], X.E(0,
                                                                                                            2)
                                   * this->data()[0] + X.E(1, 2) * this->data()[1] + X.E(2,
                                                                                         2)
                                   * this->data()[2] - X.r[1] * E_T_f[0] + X.r[0] * E_T_f[1], E_T_f[0], E_T_f[1], E_T_f[2]);
            }

            /**
             * @brief Performs the following in place transform \f[
             * f =
             * \begin{bmatrix}
             * X.E & -X.E(X.r\times) \\
             * \mathbf{0} & X.E
             * \end{bmatrix}
             * f
             * \f]
             * @param X
             */
            void transform(const SpatialTransform& X)
            {
                this->set(this->transform_copy(X));
            }

            /**
             * @brief Copy then transform a ForceVector by \f[
             * f_2 =
             * \begin{bmatrix}
             * X.E & -X.E(X.r\times) \\
             * \mathbf{0} & X.E
             * \end{bmatrix}
             * f_1
             * \f]
             * @param X
             * @return Copied, transformed ForceVector
             */
            ForceVector transform_copy(const SpatialTransform& X) const
            {
                Vector3d En_rxf = X.E * (this->getAngularPart() - X.r.cross(this->getLinearPart()));

                return ForceVector(En_rxf[0], En_rxf[1], En_rxf[2], X.E(0, 0) * this->data()[3] + X.E(0,
                                                                                                      1) * this->data()[4] +
                                   X.E(0, 2) * this->data()[5], X.E(1, 0) * this->data()[3] + X.E(1,
                                                                                                  1)
                                   * this->data()[4] + X.E(1, 2) * this->data()[5], X.E(2, 0) * this->data()[3] + X.E(2,
                                                                                                                      1)
                                   * this->data()[4] + X.E(2, 2) * this->data()[5]);
            }

            /**
             * @brief Overloaded plus-equals operator
             * @return \f$ f_2=f_2+f_1 \f$
             */
            inline ForceVector operator+=(const ForceVector& v)
            {
                (*this) << (this->mx() += v.mx()), (this->my() += v.my()), (this->mz() += v.mz()), (this->fx() += v.fx()),
                (this->fy() += v.fy()), (this->fz() += v.fz());
                return *this;
            }
        };

        /**
         * @brief Operator for transforming a ForceVector. Calls the ForceVector::transform method.
         * @param X SpatialTransform
         * @param f ForceVector to be transformed
         * @return Transformed ForceVector
         */
        inline ForceVector operator*(const SpatialTransform& X, ForceVector f)
        {
            f.transform(X);
            return f;
        }

        /**
         * @class MotionVector
         * @brief A MotionVector is a SpatialVector with 3 angular velocities and 3 linear velocities.
         */
        class RDL_DLLAPI MotionVector : public SpatialVector, public TransformableGeometricObject
        {
public:

            /**
             * @brief Constructor
             * @param other
             */
            template<typename OtherDerived>
            MotionVector(const Eigen::MatrixBase<OtherDerived>& other)
                : SpatialVector(other)
            {}

            /**
             * @brief Overload equal operator
             * @param other
             * @return A motion vector
             */
            template<typename OtherDerived>
            MotionVector& operator=(const Eigen::MatrixBase<OtherDerived>& other)
            {
                SpatialVector::operator=(other);
                return *this;
            }

            /**
             * @brief Empty constructor
             */
            EIGEN_STRONG_INLINE MotionVector()
            {}

            /**
             * @brief Constructor
             * @param v0 x-angular
             * @param v1 y-angular
             * @param v2 z-angular
             * @param v3 x-linear
             * @param v4 y-linear
             * @param v5 z-linear
             */
            MotionVector(const double v0, const double v1, const double v2, const double v3, const double v4, const double v5)
            {
                Base::_check_template_params();

                (*this) << v0, v1, v2, v3, v4, v5;
            }

            /**
             * @brief Get a copy of a MotionVector as a SpatialVector
             * @return SpatialVector copy of a MotionVectors
             */
            EIGEN_STRONG_INLINE SpatialVector toSpatialVector() const
            {
                return *this;
            }

            /**
             * @brief Setter
             * @param v Sets the values equal to those stored in v
             */
            EIGEN_STRONG_INLINE void set(const MotionVector& v)
            {
                (*this) << v.data()[0], v.data()[1], v.data()[2], v.data()[3], v.data()[4], v.data()[5];
            }

            /**
             * @brief Get a reference to the angular-x component
             * @return A reference to the angular x-component
             */
            EIGEN_STRONG_INLINE double& wx()
            {
                return this->data()[0];
            }

            /**
             * @brief Get a reference to the angular-y component
             * @return A reference to the angular y-component
             */
            EIGEN_STRONG_INLINE double& wy()
            {
                return this->data()[1];
            }

            /**
             * @brief Get a reference to the angular-z component
             * @return A copy reference to the angular z-component
             */
            EIGEN_STRONG_INLINE double& wz()
            {
                return this->data()[2];
            }

            /**
             * @brief Get a copy of the angular-x component
             * @return A copy of the angular x-component
             */
            EIGEN_STRONG_INLINE double wx() const
            {
                return this->data()[0];
            }

            /**
             * @brief Get a copy of the angular-y component
             * @return A copy of the angular y-component
             */
            EIGEN_STRONG_INLINE double wy() const
            {
                return this->data()[1];
            }

            /**
             * @brief Get a copy of the angular-z component
             * @return A copy of the angular z-component
             */
            EIGEN_STRONG_INLINE double wz() const
            {
                return this->data()[2];
            }

            /**
             * @brief Get a reference to the linear-x component
             * @return A reference to the linear x-component
             */
            EIGEN_STRONG_INLINE double& vx()
            {
                return this->data()[3];
            }

            /**
             * @brief Get a reference to the linear-y component
             * @return A reference to the linear y-component
             */
            EIGEN_STRONG_INLINE double& vy()
            {
                return this->data()[4];
            }

            /**
             * @brief Get a reference to the linear-z component
             * @return A reference to the linear z-component
             */
            EIGEN_STRONG_INLINE double& vz()
            {
                return this->data()[5];
            }

            /**
             * @brief Get a copy of the linear-x component
             * @return A copy of the linear x-component
             */
            EIGEN_STRONG_INLINE double vx() const
            {
                return this->data()[3];
            }

            /**
             * @brief Get a copy of the linear-y component
             * @return A copy of the linear y-component
             */
            EIGEN_STRONG_INLINE double vy() const
            {
                return this->data()[4];
            }

            /**
             * @brief Get a copy of the linear-z component
             * @return A copy of the linear z-component
             */
            EIGEN_STRONG_INLINE double vz() const
            {
                return this->data()[5];
            }

            /**
             * @brief Transforms a motion vector. Performs \f$ v= X v \f$
             * @param X
             */
            inline void transform(const SpatialTransform& X)
            {
                this->set(this->transform_copy(X));
            }

            /**
             * @brief Copies, transforms, and returns a MotionVector. Performs \f$ v_2=X v_1 \f$
             * @param X
             * @return Copied, transformed MotionVector
             */
            MotionVector transform_copy(const SpatialTransform& X) const
            {
                Vector3d v_rxw(this->data()[3] - X.r[1] * this->data()[2] + X.r[2] * this->data()[1],
                               this->data()[4] - X.r[2] * this->data()[0] + X.r[0] * this->data()[2],
                               this->data()[5] - X.r[0] * this->data()[1] + X.r[1] * this->data()[0]);

                return MotionVector(X.E(0, 0) * this->data()[0] + X.E(0, 1) * this->data()[1] + X.E(0,
                                                                                                    2) * this->data()[2],
                                    X.E(1, 0) * this->data()[0] + X.E(1, 1) * this->data()[1] + X.E(1,
                                                                                                    2)
                                    * this->data()[2], X.E(2, 0) * this->data()[0] + X.E(2, 1) * this->data()[1] + X.E(2,
                                                                                                                       2)
                                    * this->data()[2],
                                    X.E(0, 0) * v_rxw[0] + X.E(0, 1) * v_rxw[1] + X.E(0, 2) * v_rxw[2], X.E(1,
                                                                                                            0)
                                    * v_rxw[0] +
                                    X.E(1, 1) * v_rxw[1] + X.E(1, 2) * v_rxw[2], X.E(2, 0) * v_rxw[0] + X.E(2,
                                                                                                            1) * v_rxw[1] +
                                    X.E(2, 2) * v_rxw[2]);
            }

            /**
             * @brief See V. Duindum thesis p.25 for an explanation of what \f$ad_T\f$ operator is.
             * It is also in Featherstone p. 25 eq. 2.31 & 2.32. For featherstone notation,
             * it is essentially the \f$\times\f$ operator for spatial vectors. Given two SpatialMotion vectors, \f$v_1\f$ and
             * \f$v_2\f$, this method returns \f$ v_3 = (v_1\times) v_2 \f$. Expanded, it looks like,
             ** \f[
             * v_3 = (v_1 \times) v_2 = ad_{v_1} v_2 =
             *  \begin{bmatrix}
             * \omega_{v_1} \times & \mathbf{0} \\
             * v_{v_1}\times & \omega_{v_1}\times
             *  \end{bmatrix} v_2
             * \f]
             * The 3d vector \f$\times\f$ operator is equivalent to the \f$\sim\f$ operator. See Math::toTildeForm.
             */
            MotionVector         cross(const MotionVector& v);

            /**
             * @brief See Featherstone p. 25 eq. 2.31 & 2.32. For featherstone notation,
             * it is essentially the \f$\times *\f$ operator for spatial vectors. Given a SpatialMotion vector, \f$v_1\f$, and
             * SpatialForceVector, \f$f_1\f$, this method returns \f$ f_2 = (v_1\times *) f_1 \f$. Expanded, it looks like,
             * \f[
             *  =(v_1\times *) f_1 =
             *  \begin{bmatrix}
             * \omega_{v_1}\times & v_{v_1}\times \\
             * \mathbf{0} & \omega_{v_1}\times
             *  \end{bmatrix} f
             * \f]
             * The 3d vector \f$\times\f$ operator is equivalent to the \f$\sim\f$ operator. See Math::toTildeForm.
             */
            ForceVector          cross(const ForceVector& v);

            /**
             * @brief Get the spatial motion cross matrix,
             * \f[
             *  m\times =
             *  \begin{bmatrix}
             *   \omega\times & \mathbf{0} \\
             *   v\times & \omega\times
             *  \end{bmatrix}
             * \f]
             * @return \f$ v\times \f$
             */
            inline SpatialMatrix crossm()
            {
                return SpatialMatrix(0, -this->data()[2], this->data()[1], 0, 0, 0, this->data()[2], 0,
                                     -this->data()[0], 0, 0, 0, -this->data()[1], this->data()[0], 0, 0, 0, 0, 0,
                                     -this->data()[5],
                                     this->data()[4], 0, -this->data()[2], this->data()[1], this->data()[5], 0,
                                     -this->data()[3],
                                     this->data()[2], 0, -this->data()[0], -this->data()[4], this->data()[3], 0,
                                     -this->data()[1], this->data()[0], 0);
            }

            /**
             * @brief Get the spatial force cross matrix
             * @return \f$ v\times* \f$
             */
            inline SpatialMatrix crossf()
            {
                return SpatialMatrix(0, -this->data()[2], this->data()[1], 0, -this->data()[5], this->data()[4],
                                     this->data()[2], 0, -this->data()[0], this->data()[5], 0, -this->data()[3],
                                     -this->data()[1],
                                     this->data()[0], 0, -this->data()[4], this->data()[3], 0, 0, 0, 0, 0, -this->data()[2],
                                     this->data()[1], 0, 0, 0, this->data()[2], 0, -this->data()[0], 0, 0, 0,
                                     -this->data()[1], this->data()[0], 0);
            }

            /**
             * @brief Operator for performing the RBDA \f$\times\f$ operator for two motion vectors, i.e.
             **\f$m_2=(m_1\times)m_2\f$.
             * @param v
             * @return A MotionVector
             */
            inline MotionVector operator%=(const MotionVector& v)
            {
                return this->cross(v);
            }

            /**
             * @brief Operator for performing the RBDA \f$ \times * \f$ operator for a motion vector and a force vector, i.e.
             **\f$ f=(m\times *)f \f$.
             * @param v
             * @return A ForceVector
             */
            inline ForceVector operator%=(const ForceVector& v)
            {
                return this->cross(v);
            }

            /**
             * @brief Overloaded += operator for a MotionVector
             * @return MotionVector
             */
            inline MotionVector operator+=(const MotionVector& v)
            {
                (*this) << (this->wx() += v.wx()), (this->wy() += v.wy()), (this->wz() += v.wz()), (this->vx() += v.vx()),
                (this->vy() += v.vy()), (this->vz() += v.vz());
                return *this;
            }
        };

        /**
         * @brief Operator for transforming a MotionVector
         * @param X SpatialTransform to the desired frame
         * @param v
         * @return Transformed MotionVector
         */
        inline MotionVector operator*(const SpatialTransform& X, MotionVector v)
        {
            v.transform(X);
            return v;
        }

        /**
         * Operator for performing the spatial vector \f$\times\f$ operator.
         * @param v
         * @param v2
         * @return Returns \f$ ret = (v \times)v2 \f$
         */
        inline MotionVector operator%(MotionVector v, const MotionVector& v2)
        {
            return v.MotionVector::operator%=(v2);
        }

        /**
         * Operator for performing the spatial vector \f$\times*\f$ operator
         * @param v
         * @param v2
         * @return Returns \f$ ret = (v\times*)v2 \f$
         */
        inline ForceVector operator%(MotionVector v, const ForceVector& v2)
        {
            return v.MotionVector::operator%=(v2);
        }
    }
}

#endif //__TRANSFORMABLE_GEOMETRIC_OBJECT_HPP__
