/*
 * Original Copyright (c) 2011-2016 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 *
 * RDL - Robot Dynamics Library
 * Modifications Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef __RDL_EIGENMATH_H__
#define __RDL_EIGENMATH_H__

#include "rdl_dynamics/rdl_config.h"

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/QR>
#include <eigen3/Eigen/Eigen>

namespace RobotDynamics
{
    /** \brief Math types such as vectors and matrices and utility functions. */
    namespace Math
    {
        typedef Eigen::Matrix<double, 6, 3>Matrix63;
        typedef Eigen::VectorXd            VectorNd;
        typedef Eigen::MatrixXd            MatrixNd;
    } /* Math */
}     /* RobotDynamics */

namespace RobotDynamics
{
    namespace Math
    {
        class RDL_DLLAPI Vector3d : public Eigen::Vector3d
        {
public:

            typedef Eigen::Vector3d Base;

            template<typename OtherDerived>
            Vector3d(const Eigen::MatrixBase<OtherDerived>& other)
                : Eigen::Vector3d(other)
            {}

            template<typename OtherDerived>
            Vector3d& operator=(const Eigen::MatrixBase<OtherDerived>& other)
            {
                this->Base::operator=(other);
                return *this;
            }

            EIGEN_STRONG_INLINE Vector3d()
            {}

            EIGEN_STRONG_INLINE Vector3d(const double& v0, const double& v1, const double& v2)
            {
                Base::_check_template_params();

                (*this) << v0, v1, v2;
            }

            void set(const Eigen::Vector3d& v)
            {
                Base::_check_template_params();

                set(v[0], v[1], v[2]);
            }

            void set(const double& v0, const double& v1, const double& v2)
            {
                Base::_check_template_params();

                (*this) << v0, v1, v2;
            }
        };

        class RDL_DLLAPI Matrix3d : public Eigen::Matrix3d
        {
public:

            typedef Eigen::Matrix3d Base;

            template<typename OtherDerived>
            Matrix3d(const Eigen::MatrixBase<OtherDerived>& other)
                : Eigen::Matrix3d(other)
            {}

            template<typename OtherDerived>
            Matrix3d& operator=(const Eigen::MatrixBase<OtherDerived>& other)
            {
                this->Base::operator=(other);
                return *this;
            }

            EIGEN_STRONG_INLINE Matrix3d()
            {}

            EIGEN_STRONG_INLINE Matrix3d(const double& m00,
                                         const double& m01,
                                         const double& m02,
                                         const double& m10,
                                         const double& m11,
                                         const double& m12,
                                         const double& m20,
                                         const double& m21,
                                         const double& m22)
            {
                Base::_check_template_params();

                (*this) << m00, m01, m02, m10, m11, m12, m20, m21, m22;
            }
        };

        class RDL_DLLAPI Vector4d : public Eigen::Vector4d
        {
public:

            typedef Eigen::Vector4d Base;

            template<typename OtherDerived>
            Vector4d(const Eigen::MatrixBase<OtherDerived>& other)
                : Eigen::Vector4d(other)
            {}

            template<typename OtherDerived>
            Vector4d& operator=(const Eigen::MatrixBase<OtherDerived>& other)
            {
                this->Base::operator=(other);
                return *this;
            }

            EIGEN_STRONG_INLINE Vector4d()
            {}

            EIGEN_STRONG_INLINE Vector4d(const double& v0, const double& v1, const double& v2, const double& v3)
            {
                Base::_check_template_params();

                (*this) << v0, v1, v2, v3;
            }

            void set(const double& v0, const double& v1, const double& v2, const double& v3)
            {
                Base::_check_template_params();

                (*this) << v0, v1, v2, v3;
            }
        };

        class RDL_DLLAPI SpatialVector : public Eigen::Matrix<double, 6, 1>
        {
public:

            typedef Eigen::Matrix<double, 6, 1>Base;

            template<typename OtherDerived>
            SpatialVector(const Eigen::MatrixBase<OtherDerived>& other)
                : Eigen::Matrix<double, 6, 1>(other)
            {}

            template<typename OtherDerived>
            SpatialVector& operator=(const Eigen::MatrixBase<OtherDerived>& other)
            {
                this->Base::operator=(other);
                return *this;
            }

            EIGEN_STRONG_INLINE SpatialVector()
            {
                (*this) << 0., 0., 0., 0., 0., 0.;
            }

            EIGEN_STRONG_INLINE SpatialVector(const double& v0,
                                              const double& v1,
                                              const double& v2,
                                              const double& v3,
                                              const double& v4,
                                              const double& v5)
            {
                Base::_check_template_params();

                (*this) << v0, v1, v2, v3, v4, v5;
            }

            EIGEN_STRONG_INLINE void set(const double& v0,
                                         const double& v1,
                                         const double& v2,
                                         const double& v3,
                                         const double& v4,
                                         const double& v5)
            {
                Base::_check_template_params();

                (*this) << v0, v1, v2, v3, v4, v5;
            }

            EIGEN_STRONG_INLINE Vector3d getAngularPart() const
            {
                return Vector3d(this->data()[0], this->data()[1], this->data()[2]);
            }

            EIGEN_STRONG_INLINE Vector3d getLinearPart() const
            {
                return Vector3d(this->data()[3], this->data()[4], this->data()[5]);
            }

            inline void setAngularPart(const Vector3d& v)
            {
                this->data()[0] = v(0);
                this->data()[1] = v(1);
                this->data()[2] = v(2);
            }

            inline void setLinearPart(const Vector3d& v)
            {
                this->data()[3] = v(0);
                this->data()[4] = v(1);
                this->data()[5] = v(2);
            }

            EIGEN_STRONG_INLINE void set(const Vector3d& angularPart, const Vector3d& linearPart)
            {
                Base::_check_template_params();

                (*this) << angularPart[0], angularPart[1], angularPart[2], linearPart[0], linearPart[1], linearPart[2];
            }
        };

        class RDL_DLLAPI Matrix4d : public Eigen::Matrix<double, 4, 4>
        {
public:

            typedef Eigen::Matrix<double, 4, 4>Base;

            template<typename OtherDerived>
            Matrix4d(const Eigen::MatrixBase<OtherDerived>& other)
                : Eigen::Matrix<double, 4, 4>(other)
            {}

            template<typename OtherDerived>
            Matrix4d& operator=(const Eigen::MatrixBase<OtherDerived>& other)
            {
                this->Base::operator=(other);
                return *this;
            }

            EIGEN_STRONG_INLINE Matrix4d()
            {}

            EIGEN_STRONG_INLINE Matrix4d(const Scalar& m00,
                                         const Scalar& m01,
                                         const Scalar& m02,
                                         const Scalar& m03,
                                         const Scalar& m10,
                                         const Scalar& m11,
                                         const Scalar& m12,
                                         const Scalar& m13,
                                         const Scalar& m20,
                                         const Scalar& m21,
                                         const Scalar& m22,
                                         const Scalar& m23,
                                         const Scalar& m30,
                                         const Scalar& m31,
                                         const Scalar& m32,
                                         const Scalar& m33)
            {
                Base::_check_template_params();

                (*this) << m00, m01, m02, m03, m10, m11, m12, m13, m20, m21, m22, m23, m30, m31, m32, m33;
            }

            void set(const Scalar& m00,
                     const Scalar& m01,
                     const Scalar& m02,
                     const Scalar& m03,
                     const Scalar& m10,
                     const Scalar& m11,
                     const Scalar& m12,
                     const Scalar& m13,
                     const Scalar& m20,
                     const Scalar& m21,
                     const Scalar& m22,
                     const Scalar& m23,
                     const Scalar& m30,
                     const Scalar& m31,
                     const Scalar& m32,
                     const Scalar& m33)
            {
                Base::_check_template_params();

                (*this) << m00, m01, m02, m03, m10, m11, m12, m13, m20, m21, m22, m23, m30, m31, m32, m33;
            }
        };

        class RDL_DLLAPI SpatialMatrix : public Eigen::Matrix<double, 6, 6>
        {
public:

            typedef Eigen::Matrix<double, 6, 6>Base;

            template<typename OtherDerived>
            SpatialMatrix(const Eigen::MatrixBase<OtherDerived>& other)
                : Eigen::Matrix<double, 6, 6>(other)
            {}

            template<typename OtherDerived>
            SpatialMatrix& operator=(const Eigen::MatrixBase<OtherDerived>& other)
            {
                this->Base::operator=(other);
                return *this;
            }

            EIGEN_STRONG_INLINE SpatialMatrix()
            {
                Base::_check_template_params();

                (*this) << 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.;
            }

            EIGEN_STRONG_INLINE SpatialMatrix(const Scalar& m00,
                                              const Scalar& m01,
                                              const Scalar& m02,
                                              const Scalar& m03,
                                              const Scalar& m04,
                                              const Scalar& m05,
                                              const Scalar& m10,
                                              const Scalar& m11,
                                              const Scalar& m12,
                                              const Scalar& m13,
                                              const Scalar& m14,
                                              const Scalar& m15,
                                              const Scalar& m20,
                                              const Scalar& m21,
                                              const Scalar& m22,
                                              const Scalar& m23,
                                              const Scalar& m24,
                                              const Scalar& m25,
                                              const Scalar& m30,
                                              const Scalar& m31,
                                              const Scalar& m32,
                                              const Scalar& m33,
                                              const Scalar& m34,
                                              const Scalar& m35,
                                              const Scalar& m40,
                                              const Scalar& m41,
                                              const Scalar& m42,
                                              const Scalar& m43,
                                              const Scalar& m44,
                                              const Scalar& m45,
                                              const Scalar& m50,
                                              const Scalar& m51,
                                              const Scalar& m52,
                                              const Scalar& m53,
                                              const Scalar& m54,
                                              const Scalar& m55)
            {
                Base::_check_template_params();

                (*this) << m00, m01, m02, m03, m04, m05, m10, m11, m12, m13, m14, m15, m20, m21, m22, m23, m24, m25, m30, m31,
                m32, m33, m34, m35, m40, m41, m42, m43, m44, m45, m50, m51, m52, m53, m54, m55;
            }

            void set(const Scalar& m00,
                     const Scalar& m01,
                     const Scalar& m02,
                     const Scalar& m03,
                     const Scalar& m04,
                     const Scalar& m05,
                     const Scalar& m10,
                     const Scalar& m11,
                     const Scalar& m12,
                     const Scalar& m13,
                     const Scalar& m14,
                     const Scalar& m15,
                     const Scalar& m20,
                     const Scalar& m21,
                     const Scalar& m22,
                     const Scalar& m23,
                     const Scalar& m24,
                     const Scalar& m25,
                     const Scalar& m30,
                     const Scalar& m31,
                     const Scalar& m32,
                     const Scalar& m33,
                     const Scalar& m34,
                     const Scalar& m35,
                     const Scalar& m40,
                     const Scalar& m41,
                     const Scalar& m42,
                     const Scalar& m43,
                     const Scalar& m44,
                     const Scalar& m45,
                     const Scalar& m50,
                     const Scalar& m51,
                     const Scalar& m52,
                     const Scalar& m53,
                     const Scalar& m54,
                     const Scalar& m55)
            {
                Base::_check_template_params();

                (*this) << m00, m01, m02, m03, m04, m05, m10, m11, m12, m13, m14, m15, m20, m21, m22, m23, m24, m25, m30, m31,
                m32, m33, m34, m35, m40, m41, m42, m43, m44, m45, m50, m51, m52, m53, m54, m55;
            }
        };
    }
}

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RobotDynamics::Math::SpatialVector)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RobotDynamics::Math::SpatialMatrix)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RobotDynamics::Math::Matrix63)

/* ___RDL_EIGENMATH_H__ */
#endif // ifndef __RDL_EIGENMATH_H__
