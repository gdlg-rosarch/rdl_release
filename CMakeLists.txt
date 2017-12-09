CMAKE_MINIMUM_REQUIRED(VERSION 2.8.3)
PROJECT (rdl_dynamics)

# check c++11 / c++0x
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "-std=c++0x")
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

find_package(catkin REQUIRED rdl_cmake)
find_package(Eigen3 REQUIRED )

SET(INCLUDES_DIR ${PROJECT_SOURCE_DIR}/include)
# Source files for RDL
SET ( RDL_SOURCES
    src/rdl_mathutils.cc
    src/rdl_utils.cc
	src/Contacts.cc
	src/Dynamics.cc
	src/Joint.cc
	src/Model.cc
	src/Kinematics.cc
    src/ReferenceFrame.cpp
    src/SpatialMotion.cpp
    src/SpatialAcceleration.cpp
    src/TransformableGeometricObject.cpp
    src/FrameObject.cpp
    src/RigidBodyInertia.cpp
	)

include_directories( ${INCLUDES_DIR} ${EIGEN3_INCLUDE_DIR})

catkin_package(
 INCLUDE_DIRS ${INCLUDES_DIR} ${EIGEN3_INCLUDE_DIR}
 LIBRARIES ${PROJECT_NAME}
# CATKIN_DEPENDS 
 DEPENDS Eigen3
)

include_directories( ${INCLUDES_DIR} ${EIGEN3_INCLUDE_DIR} ${catkin_INCLUDE_DIRS})

add_library(${PROJECT_NAME} ${RDL_SOURCES})
target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -Wno-sign-compare -Wno-unused-function -Wno-unused-variable")

set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS})

install(
  TARGETS 
    ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
  FILES_MATCHING PATTERN "*.h*"
)

include(Format)
include(Lint)

lint_library(INCLUDE_DIRS include/${PROJECT_NAME} SRCS ${RDL_SOURCES})
format_code_target()

# Add gtest based cpp test target and link libraries
if (CATKIN_ENABLE_TESTING)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

    catkin_add_gtest(ReferenceFrameTests tests/ReferenceFrameTest.cpp)
    catkin_add_gtest(Point3Tests tests/Point3Test.cpp)
    catkin_add_gtest(FramePointTests tests/FramePointTest.cpp)
    catkin_add_gtest(FrameVectorTests tests/FrameVectorTest.cpp)
    catkin_add_gtest(ModelFrameTests tests/ModelFrameTests.cpp)
    catkin_add_gtest(RigidBodyInertiaTests tests/RigidBodyInertiaTests.cpp)
    catkin_add_gtest(MotionVectorTests tests/MotionVectorTests.cpp)
    catkin_add_gtest(ForceVectorTests tests/ForceVectorTests.cpp)
    catkin_add_gtest(SpatialMotionTests tests/SpatialMotionTests.cpp)
    catkin_add_gtest(SpatialForceTests tests/SpatialForceTests.cpp)
    catkin_add_gtest(SpatialMomentumTests tests/SpatialMomentumTests.cpp)
    catkin_add_gtest(SpatialRigidBodyInertiaTests tests/SpatialRigidBodyInertiaTests.cpp)
    catkin_add_gtest(MomentumTests tests/MomentumTests.cpp)
    catkin_add_gtest(SpatialAccelerationTests tests/SpatialAccelerationTests.cpp)
    catkin_add_gtest(RdlCalcAccelerationsTests tests/RdlCalcAccelerationTests.cpp)
    catkin_add_gtest(RdlBodyTests tests/RdlBodyTests.cpp)
    catkin_add_gtest(RdlKinematicsTests tests/RdlKinematicsTests.cpp)
    catkin_add_gtest(RdlCalcVelocitiesTests tests/RdlCalcVelocitiesTests.cpp)
    catkin_add_gtest(RdlDynamicsTests tests/RdlDynamicsTests.cpp)
    catkin_add_gtest(RdlContactsTests tests/RdlContactsTests.cpp)
    catkin_add_gtest(RdlFloatingBaseTests tests/RdlFloatingBaseTests.cpp)
    catkin_add_gtest(RdlMultiDofTests tests/RdlMultiDofTests.cpp)
    catkin_add_gtest(RdlCompositeRigidBodyTests tests/RdlCompositeRigidBodyTests.cpp)
    catkin_add_gtest(RdlSparseFactorizationTests tests/RdlSparseFactorizationTests.cpp)
    catkin_add_gtest(RdlInverseDynamicsTests tests/RdlInverseDynamicsTests.cpp)
    catkin_add_gtest(RdlTwoLegModelTests tests/RdlTwoLegModelTests.cpp)
    catkin_add_gtest(RdlImpulsesTests tests/RdlImpulsesTests.cpp)
    catkin_add_gtest(MathTests tests/MathTests.cc)
    catkin_add_gtest(RdlModelTests tests/RdlModelTests.cpp)
    catkin_add_gtest(SpatialAlgebraTests tests/SpatialAlgebraTests.cc)
    catkin_add_gtest(RdlUtilsTests tests/RdlUtilsTests.cpp)
    catkin_add_gtest(UnitTestUtilsTests tests/UnitTestUtilsTests.cpp)
    catkin_add_gtest(JointTests tests/JointTests.cpp)
    catkin_add_gtest(RdlCustomJointTests tests/RdlCustomJointTests.cpp)
    catkin_add_gtest(RdlCustomJointSingleBodyTests tests/RdlCustomJointSingleBodyTests.cpp)
    catkin_add_gtest(RdlCustomJointMultiBodyTests tests/RdlCustomJointMultiBodyTests.cpp)
    catkin_add_gtest(QuaternionTests tests/QuaternionTests.cpp)
    catkin_add_gtest(FrameOrientationTests tests/FrameOrientationTests.cpp)
	#
    
    target_link_libraries(ReferenceFrameTests ${PROJECT_NAME})
    target_link_libraries(Point3Tests ${PROJECT_NAME})
    target_link_libraries(FramePointTests ${PROJECT_NAME})
    target_link_libraries(FrameVectorTests ${PROJECT_NAME})
    target_link_libraries(ModelFrameTests ${PROJECT_NAME})
    target_link_libraries(RigidBodyInertiaTests ${PROJECT_NAME})
    target_link_libraries(MotionVectorTests ${PROJECT_NAME})
    target_link_libraries(SpatialMotionTests ${PROJECT_NAME})
    target_link_libraries(ForceVectorTests ${PROJECT_NAME})
    target_link_libraries(SpatialForceTests ${PROJECT_NAME})
    target_link_libraries(SpatialMomentumTests ${PROJECT_NAME})
    target_link_libraries(SpatialRigidBodyInertiaTests ${PROJECT_NAME})
    target_link_libraries(MomentumTests ${PROJECT_NAME})
    target_link_libraries(SpatialAccelerationTests ${PROJECT_NAME})
    target_link_libraries(RdlKinematicsTests ${PROJECT_NAME})
    target_link_libraries(RdlCalcVelocitiesTests ${PROJECT_NAME})
    target_link_libraries(RdlBodyTests ${PROJECT_NAME})
    target_link_libraries(RdlCalcAccelerationsTests ${PROJECT_NAME})
    target_link_libraries(RdlDynamicsTests ${PROJECT_NAME})
    target_link_libraries(RdlContactsTests ${PROJECT_NAME})
    target_link_libraries(RdlFloatingBaseTests ${PROJECT_NAME})
    target_link_libraries(RdlMultiDofTests ${PROJECT_NAME})
    target_link_libraries(RdlCompositeRigidBodyTests ${PROJECT_NAME})
    target_link_libraries(RdlSparseFactorizationTests ${PROJECT_NAME})
    target_link_libraries(RdlInverseDynamicsTests ${PROJECT_NAME})
    target_link_libraries(RdlTwoLegModelTests ${PROJECT_NAME})
    target_link_libraries(RdlImpulsesTests ${PROJECT_NAME})
    target_link_libraries(MathTests ${PROJECT_NAME})
    target_link_libraries(RdlModelTests ${PROJECT_NAME})
    target_link_libraries(SpatialAlgebraTests ${PROJECT_NAME})
    target_link_libraries(RdlUtilsTests ${PROJECT_NAME})
    target_link_libraries(UnitTestUtilsTests ${PROJECT_NAME})
    target_link_libraries(JointTests ${PROJECT_NAME})
    target_link_libraries(RdlCustomJointTests ${PROJECT_NAME})
    target_link_libraries(RdlCustomJointSingleBodyTests ${PROJECT_NAME})
    target_link_libraries(RdlCustomJointMultiBodyTests ${PROJECT_NAME})
    target_link_libraries(QuaternionTests ${PROJECT_NAME})
    target_link_libraries(FrameOrientationTests ${PROJECT_NAME})

    include(CodeCoverage)
    coverage_add_target(run_tests_${PROJECT_NAME})

endif(CATKIN_ENABLE_TESTING)