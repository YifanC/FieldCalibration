#----------------------------------------------------------------------------
# Setup the project
#
project(FieldCalibration C CXX)

# An option to enable testing, enable it with "cmake blabla -DTESTING=On"
option(TESTING "Testing" OFF)

#----------------------------------------------------------------------------
# cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
#
cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

#----------------------------------------------------------------------------
# Use C++11 extention and O3 optimization (-founding-math because of CGAL)
#
if (APPLE)
  message(STATUS "No optimization for mac :-)")
else()
  set(CMAKE_CXX_FLAGS "-O3 -frounding-math")
endif()
add_definitions("-std=c++11")

#----------------------------------------------------------------------------
# Find CGAL library
#
find_package(CGAL QUIET)
if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )
else()
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED)
add_definitions(${ROOT_CXX_FLAGS})

#----------------------------------------------------------------------------
# Find package for multi-threading
#
find_package (Threads)

find_package(OpenMP)
if (OPENMP_FOUND)
  message(STATUS "Using OpenMP")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#----------------------------------------------------------------------------
# Add sources and headers
#
include_directories(./include)

set(SOURCE_FILES
        source/Laser.cpp
        source/LaserTrack.cpp
        source/Matrix3x3.cpp
        #source/ThreeVector.cpp
        source/TPCVolumeHandler.cpp
        source/Interpolation3D.cpp
        source/Utilities.cpp
        source/DriftVelocity.cpp
        source/EfieldCalculation.cpp
        source/LaserDisplacementIteration.cpp
        source/WeightAverage.cpp
       # source/VoxelMesh.cpp
        )

#set(HEADER_FILES include/Laser.hpp include/LaserTrack.hpp include/Laser.hpp include/ThreeVector.hpp include/Matrix3x3.hpp include/TPCVolumeHandler include/Interpolation3D.hpp source/Utilities.cpp include/Utilities.hpp)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
if (APPLE)
  message(STATUS "So Yifan is trying somethings.")
  execute_process(COMMAND "root-config" "--incdir" OUTPUT_VARIABLE rootinc RESULT_VARIABLE rs)
  include_directories(${rootinc})
else()
  include_directories(${ROOT_INCLUDE_DIR})
endif()


#----------------------------------------------------------------------------
# Add the executable, and link it to the necessary libraries
#
if(TESTING)
  add_library(fieldcal SHARED ${SOURCE_FILES} )
endif()

add_executable(FieldCal FieldCal.cpp ${SOURCE_FILES} )
target_link_libraries(FieldCal ${ROOT_LIBRARIES} ${CGAL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

#----------------------------------------------------------------------------
# Copy all required files to build directory.
#
#set(FieldCal_Files Field.root)

#foreach(_File ${FieldCal_Files})
#  configure_file(
#    ${PROJECT_SOURCE_DIR}/${_File}
#    ${PROJECT_BINARY_DIR}/${_File}
#    COPYONLY
#    )
#endforeach()


# Testing stuff
if(TESTING)
  message(STATUS "Testing enabled.")
  find_package(GTest REQUIRED)
  if( GTEST_FOUND )
    enable_testing()
    add_subdirectory(tests)
  else()
    message(STATUS "Cannot find GTest, which is required if testing is enabled. Either install GTest or disable
            Testing again.")
  endif()
endif()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS FieldCal DESTINATION ${CMAKE_SOURCE_DIR}/bin)
#install(TARGETS fieldcal DESTINATION ${CMAKE_SOURCE_DIR}/lib )

