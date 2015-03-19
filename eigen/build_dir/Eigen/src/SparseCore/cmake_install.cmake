# Install script for directory: /home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore

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
   "/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/SparseFuzzy.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/SparseProduct.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/SparseVector.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/TriangularSolver.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/SparseDenseProduct.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/CompressedStorage.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/SparseSelfAdjointView.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/SparseTriangularView.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/SparseDot.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/MappedSparseMatrix.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/SparseView.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/SparseCwiseBinaryOp.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/SparseDiagonalProduct.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/SparseUtil.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/SparseMatrix.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/SparsePermutation.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/ConservativeSparseSparseProduct.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/SparseRedux.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/AmbiVector.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/SparseSparseProductWithPruning.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/SparseColEtree.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/SparseTranspose.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/SparseCwiseUnaryOp.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/SparseBlock.h;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore/SparseMatrixBase.h")
FILE(INSTALL DESTINATION "/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/SparseCore" TYPE FILE FILES
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/SparseFuzzy.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/SparseProduct.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/SparseVector.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/TriangularSolver.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/SparseDenseProduct.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/CompressedStorage.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/SparseSelfAdjointView.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/SparseTriangularView.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/SparseDot.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/MappedSparseMatrix.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/SparseView.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/SparseCwiseBinaryOp.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/SparseDiagonalProduct.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/SparseUtil.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/SparseMatrix.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/SparsePermutation.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/ConservativeSparseSparseProduct.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/SparseRedux.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/AmbiVector.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/SparseSparseProductWithPruning.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/SparseColEtree.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/SparseTranspose.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/SparseCwiseUnaryOp.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/SparseBlock.h"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/src/SparseCore/SparseMatrixBase.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")

