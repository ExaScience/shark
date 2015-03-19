# Install script for directory: /home/chakro23/exa2ct/shark/eigen/Eigen

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
   "/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/UmfPackSupport;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/StdList;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/PardisoSupport;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/LU;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/StdDeque;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/Dense;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/MetisSupport;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/SparseCore;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/SparseQR;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/Core;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/LeastSquares;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/OrderingMethods;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/SVD;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/IterativeLinearSolvers;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/SuperLUSupport;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/SparseCholesky;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/SPQRSupport;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/Array;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/Sparse;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/Jacobi;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/Eigen;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/Cholesky;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/Geometry;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/QR;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/QtAlignedMalloc;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/CholmodSupport;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/StdVector;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/Householder;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/Eigenvalues;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/Eigen2Support;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/PaStiXSupport;/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/SparseLU")
FILE(INSTALL DESTINATION "/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen" TYPE FILE FILES
    "/home/chakro23/exa2ct/shark/eigen/Eigen/UmfPackSupport"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/StdList"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/PardisoSupport"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/LU"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/StdDeque"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/Dense"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/MetisSupport"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/SparseCore"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/SparseQR"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/Core"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/LeastSquares"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/OrderingMethods"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/SVD"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/IterativeLinearSolvers"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/SuperLUSupport"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/SparseCholesky"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/SPQRSupport"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/Array"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/Sparse"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/Jacobi"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/Eigen"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/Cholesky"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/Geometry"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/QR"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/QtAlignedMalloc"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/CholmodSupport"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/StdVector"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/Householder"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/Eigenvalues"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/Eigen2Support"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/PaStiXSupport"
    "/home/chakro23/exa2ct/shark/eigen/Eigen/SparseLU"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/home/chakro23/exa2ct/shark/eigen/build_dir/Eigen/src/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

