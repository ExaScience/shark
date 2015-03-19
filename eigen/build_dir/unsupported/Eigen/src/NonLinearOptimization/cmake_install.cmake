# Install script for directory: /home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/src/NonLinearOptimization

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
   "/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/src/NonLinearOptimization/fdjac1.h;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/src/NonLinearOptimization/HybridNonLinearSolver.h;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/src/NonLinearOptimization/rwupdt.h;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/src/NonLinearOptimization/qrsolv.h;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/src/NonLinearOptimization/dogleg.h;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/src/NonLinearOptimization/r1mpyq.h;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/src/NonLinearOptimization/chkder.h;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/src/NonLinearOptimization/lmpar.h;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/src/NonLinearOptimization/covar.h;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/src/NonLinearOptimization/LevenbergMarquardt.h;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/src/NonLinearOptimization/r1updt.h")
FILE(INSTALL DESTINATION "/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/src/NonLinearOptimization" TYPE FILE FILES
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/src/NonLinearOptimization/fdjac1.h"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/src/NonLinearOptimization/HybridNonLinearSolver.h"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/src/NonLinearOptimization/rwupdt.h"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/src/NonLinearOptimization/qrsolv.h"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/src/NonLinearOptimization/dogleg.h"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/src/NonLinearOptimization/r1mpyq.h"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/src/NonLinearOptimization/chkder.h"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/src/NonLinearOptimization/lmpar.h"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/src/NonLinearOptimization/covar.h"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/src/NonLinearOptimization/LevenbergMarquardt.h"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/src/NonLinearOptimization/r1updt.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")

