# Install script for directory: /home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products

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
   "/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/CoeffBasedProduct.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/TriangularMatrixVector.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/TriangularSolverMatrix_MKL.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/Parallelizer.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/SelfadjointMatrixMatrix_MKL.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/SelfadjointMatrixVector.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/TriangularSolverMatrix.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/TriangularSolverVector.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/TriangularMatrixMatrix_MKL.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/GeneralMatrixMatrixTriangular_MKL.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/GeneralBlockPanelKernel.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/SelfadjointProduct.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/GeneralMatrixMatrix.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/SelfadjointRank2Update.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/TriangularMatrixMatrix.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/TriangularMatrixVector_MKL.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/GeneralMatrixMatrixTriangular.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/GeneralMatrixVector.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/SelfadjointMatrixVector_MKL.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/GeneralMatrixVector_MKL.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/GeneralMatrixMatrix_MKL.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products/SelfadjointMatrixMatrix.h")
FILE(INSTALL DESTINATION "/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/Core/products" TYPE FILE FILES
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/CoeffBasedProduct.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/TriangularMatrixVector.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/TriangularSolverMatrix_MKL.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/Parallelizer.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/SelfadjointMatrixMatrix_MKL.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/SelfadjointMatrixVector.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/TriangularSolverMatrix.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/TriangularSolverVector.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/TriangularMatrixMatrix_MKL.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/GeneralMatrixMatrixTriangular_MKL.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/GeneralBlockPanelKernel.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/SelfadjointProduct.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/GeneralMatrixMatrix.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/SelfadjointRank2Update.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/TriangularMatrixMatrix.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/TriangularMatrixVector_MKL.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/GeneralMatrixMatrixTriangular.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/GeneralMatrixVector.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/SelfadjointMatrixVector_MKL.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/GeneralMatrixVector_MKL.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/GeneralMatrixMatrix_MKL.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/Core/products/SelfadjointMatrixMatrix.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")

