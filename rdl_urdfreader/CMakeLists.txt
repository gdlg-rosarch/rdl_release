CMAKE_MINIMUM_REQUIRED(VERSION 2.8.3)
PROJECT(rdl_urdfreader)

find_package(catkin REQUIRED COMPONENTS roscpp urdf rdl_dynamics rdl_cmake roslib rdl_cmake)

# check c++11 / c++0x
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if (COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "-std=c++0x")
else ()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif ()

include_directories(include ${catkin_INCLUDE_DIRS} )

catkin_package(
        INCLUDE_DIRS include
        LIBRARIES ${PROJECT_NAME}
        CATKIN_DEPENDS roscpp urdf rdl_dynamics roslib
#        DEPENDS
)

SET(URDFREADER_SOURCES
        src/urdfreader.cc
        )

add_library(${PROJECT_NAME} ${URDFREADER_SOURCES})

target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -Wno-sign-compare -Wno-unused-function -Wno-unused-variable")
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS})

add_executable(${PROJECT_NAME}_util src/rdl_urdfreader_util.cc)

target_link_libraries(${PROJECT_NAME}_util ${PROJECT_NAME})

install(DIRECTORY include/${PROJECT_NAME}/
        DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
        FILES_MATCHING PATTERN "*.h"
        )

install(TARGETS ${PROJECT_NAME} ${PROJECT_NAME}_util
        ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
        LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
        RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
        )

include(Format)
include(Lint)

set(SRC_FILES src/urdfreader.cc)

lint_library(INCLUDE_DIRS include/${PROJECT_NAME} SRCS ${SRC_FILES})
format_code_target()

if (CATKIN_ENABLE_TESTING)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

    catkin_add_gtest(UrdfReaderTests tests/UrdfReaderTests.cpp)
    target_link_libraries(UrdfReaderTests ${PROJECT_NAME})

    set_target_properties(UrdfReaderTests PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS})

    install(DIRECTORY tests/urdf/
            DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}/tests/urdf
            )

    include(CodeCoverage)
    coverage_add_target(run_tests_${PROJECT_NAME})
endif ()