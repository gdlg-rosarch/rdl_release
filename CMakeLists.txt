CMAKE_MINIMUM_REQUIRED(VERSION 2.8.3)
PROJECT (rdl_benchmark)

find_package(catkin REQUIRED COMPONENTS rdl_cmake rdl_dynamics rdl_urdfreader)

# check c++11 / c++0x
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "-std=c++0x")
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

set(INCLUDES_DIR include)

# Options
set( BENCHMARK_SOURCES 
  src/model_generator.cc
  src/Human36Model.cc
  src/benchmark.cc
  )

catkin_package(
 INCLUDE_DIRS ${INCLUDES_DIR} 
 # LIBRARIES ${PROJECT_NAME}
CATKIN_DEPENDS rdl_dynamics rdl_urdfreader
 # DEPENDS 
)

include_directories(${INCLUDES_DIR} ${catkin_INCLUDE_DIRS})

add_executable( ${PROJECT_NAME} ${BENCHMARK_SOURCES} )

target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES})

install(TARGETS ${PROJECT_NAME} 
    ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  	LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  	RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
    )