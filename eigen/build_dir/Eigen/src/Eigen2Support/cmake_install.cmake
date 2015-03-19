# Install script for directory: /home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigen2Support

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/home/chakro23/exa2ct/shark/eigen/build_dir")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigen2Support/CwiseOperators.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigen2Support/MathFunctions.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigen2Support/Block.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigen2Support/TriangularSolver.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigen2Support/QR.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigen2Support/Lazy.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigen2Support/Cwise.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigen2Support/LU.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigen2Support/SVD.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigen2Support/Memory.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigen2Support/Macros.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigen2Support/VectorBlock.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigen2Support/Meta.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigen2Support/LeastSquares.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigen2Support/Minor.h")
FILE(INSTALL DESTINATION "/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigen2Support" TYPE FILE FILES
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigen2Support/CwiseOperators.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigen2Support/MathFunctions.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigen2Support/Block.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigen2Support/TriangularSolver.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigen2Support/QR.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigen2Support/Lazy.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigen2Support/Cwise.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigen2Support/LU.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigen2Support/SVD.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigen2Support/Memory.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigen2Support/Macros.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigen2Support/VectorBlock.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigen2Support/Meta.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigen2Support/LeastSquares.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigen2Support/Minor.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigen2Support/Geometry/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

