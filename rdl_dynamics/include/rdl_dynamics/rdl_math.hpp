/*
 * RDL - Robot Dynamics Library
 * Copyright (c) 2017 Jordan Lack <jlack1987@gmail.com>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef __RDL_MATH_HPP__
#define __RDL_MATH_HPP__

#include "rdl_dynamics/rdl_eigenmath.h"
#include "rdl_dynamics/FramePoint.hpp"
#include "rdl_dynamics/RigidBodyInertia.hpp"
#include "rdl_dynamics/ReferenceFrame.hpp"
#include "rdl_dynamics/FrameVector.hpp"
#include "rdl_dynamics/SpatialMotion.hpp"
#include "rdl_dynamics/SpatialForce.hpp"
#include "rdl_dynamics/SpatialMomentum.hpp"
#include "rdl_dynamics/SpatialRigidBodyInertia.hpp"
#include "rdl_dynamics/SpatialAcceleration.hpp"

// If we use Eigen3 we have to create specializations of the STL
// std::vector such that the alignment is done properly.
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RobotDynamics::Math::RigidBodyInertia);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RobotDynamics::Math::SpatialAcceleration);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RobotDynamics::Math::SpatialInertia);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RobotDynamics::Math::Vector3d_X);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RobotDynamics::Math::MotionVector);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RobotDynamics::Math::ForceVector);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RobotDynamics::Math::SpatialForce);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RobotDynamics::Math::SpatialMotion);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RobotDynamics::Math::SpatialMomentum);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RobotDynamics::ReferenceFrame);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(RobotDynamics::Math::FrameVector);

/* __RDL_MATH_HPP__ */
#endif // ifndef __RDL_MATH_HPP__
