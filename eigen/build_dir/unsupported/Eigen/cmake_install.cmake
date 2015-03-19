# Install script for directory: /home/chakro23/exa2ct/shark/eigen/unsupported/Eigen

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
   "/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/AdolcForward;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/BVH;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/IterativeSolvers;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/MatrixFunctions;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/MoreVectorization;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/AutoDiff;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/AlignedVector3;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/Polynomials;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/FFT;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/NonLinearOptimization;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/SparseExtra;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/IterativeSolvers;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/NumericalDiff;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/Skyline;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/MPRealSupport;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/OpenGLSupport;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/KroneckerProduct;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/Splines;/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/LevenbergMarquardt")
FILE(INSTALL DESTINATION "/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen" TYPE FILE FILES
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/AdolcForward"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/BVH"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/IterativeSolvers"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/MatrixFunctions"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/MoreVectorization"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/AutoDiff"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/AlignedVector3"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/Polynomials"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/FFT"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/NonLinearOptimization"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/SparseExtra"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/IterativeSolvers"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/NumericalDiff"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/Skyline"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/MPRealSupport"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/OpenGLSupport"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/KroneckerProduct"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/Splines"
    "/home/chakro23/exa2ct/shark/eigen/unsupported/Eigen/LevenbergMarquardt"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/home/chakro23/exa2ct/shark/eigen/build_dir/unsupported/Eigen/src/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

