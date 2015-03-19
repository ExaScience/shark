# Install script for directory: /home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigenvalues

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
   "/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigenvalues/GeneralizedEigenSolver.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigenvalues/Tridiagonalization.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigenvalues/RealSchur_MKL.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigenvalues/HessenbergDecomposition.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigenvalues/ComplexSchur_MKL.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigenvalues/ComplexEigenSolver.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigenvalues/MatrixBaseEigenvalues.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigenvalues/SelfAdjointEigenSolver_MKL.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigenvalues/RealQZ.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigenvalues/GeneralizedSelfAdjointEigenSolver.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigenvalues/RealSchur.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigenvalues/EigenSolver.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigenvalues/ComplexSchur.h")
FILE(INSTALL DESTINATION "/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Eigenvalues" TYPE FILE FILES
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigenvalues/GeneralizedEigenSolver.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigenvalues/Tridiagonalization.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigenvalues/RealSchur_MKL.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigenvalues/HessenbergDecomposition.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigenvalues/ComplexSchur_MKL.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigenvalues/ComplexEigenSolver.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigenvalues/MatrixBaseEigenvalues.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigenvalues/SelfAdjointEigenSolver_MKL.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigenvalues/RealQZ.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigenvalues/GeneralizedSelfAdjointEigenSolver.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigenvalues/RealSchur.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigenvalues/EigenSolver.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Eigenvalues/ComplexSchur.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")

