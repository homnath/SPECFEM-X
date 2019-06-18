# cmake file for the cross-platform configure of the
# software package SPECFEMX
# REVISION:
#   HNG, Mar 24,2017; HNG, Sep 06,2011;

# required cmake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.4 FATAL_ERROR)

# define project
PROJECT(semgeotech Fortran C)

#SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules")

# Disallow in-source build
IF("${semgeotech_SOURCE_DIR}"  STREQUAL "${semgeotech_BINARY_DIR}")
  MESSAGE(FATAL_ERROR "SPECFEMX requires an out-of-source build. Please create a separate build directory and run CMake there.")
ENDIF("${semgeotech_SOURCE_DIR}"  STREQUAL "${semgeotech_BINARY_DIR}")

# Set a default build type if none was specified                                
IF(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)                      
  MESSAGE(STATUS "Setting build type to 'Debug' as none was specified.")        
  SET(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)    
  # Set the possible values of build type for cmake-gui                         
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"        
  "MinSizeRel" "RelWithDebInfo")                                              
ENDIF()

# Test Fortran compiler
TRY_COMPILE(TEST_F90 ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/src/testf90.f90)
IF(NOT TEST_F90)
  MESSAGE(FATAL_ERROR "Fortran compiler: ${CMAKE_Fortran_COMPILER} doesn't work! Please set working CMAKE_Fortran_COMPILER!")
ENDIF(NOT TEST_F90)

# build directories
SET(BDIR ${CMAKE_BINARY_DIR}/bin)
SET(DDIR ${CMAKE_BINARY_DIR}/doc)
SET(IDIR ${CMAKE_BINARY_DIR}/input)
SET(LDIR ${CMAKE_BINARY_DIR}/lib)
SET(ODIR ${CMAKE_BINARY_DIR}/output)
SET(PDIR ${CMAKE_BINARY_DIR}/partition)
SET(PDIR ${CMAKE_BINARY_DIR}/tmp)
SET(SDIR ${CMAKE_SOURCE_DIR})
FILE(MAKE_DIRECTORY ${BDIR})
#FILE(MAKE_DIRECTORY ${DDIR})
FILE(MAKE_DIRECTORY ${IDIR}) 
FILE(MAKE_DIRECTORY ${ODIR}) 
FILE(MAKE_DIRECTORY ${PDIR}) 

# set paths to executable and libraries
SET(EXECUTABLE_OUTPUT_PATH ${BDIR})
SET(LIBRARY_OUTPUT_PATH ${LDIR})

# set internal default logical variables
SET(ISDOCUMENTATION OFF)
SET(ISPARTMESH OFF)
SET(ISPARALLEL OFF)

# build documentation option
SET(BUILD_DOCUMENTATION OFF CACHE BOOL "Build documentation for
SPECFEM3D_GEOTECH")
IF(BUILD_DOCUMENTATION)
  SET(ISDOCUMENTATION ON)
  FIND_PACKAGE(LATEX REQUIRED)
  IF(NOT PDFLATEX_COMPILER OR NOT BIBTEX_COMPILER)
    MESSAGE(WARNING "LaTeX compiler/s not found! DOCUMENTATION cannot be built!"    )
  SET(ISDOCUMENTATION OFF)
  ENDIF(NOT PDFLATEX_COMPILER OR NOT BIBTEX_COMPILER)
ENDIF(BUILD_DOCUMENTATION)

#SET(BUILD_SLOPE ON CACHE BOOL "Build slope stability program")
#SET(BUILD_EXCAVATION ON CACHE BOOL "Build excavation program")

# build partmesh option
SET(BUILD_PARTMESH OFF CACHE BOOL "Build mesh partitioning program")
IF(BUILD_PARTMESH)
  SET(ISPARTMESH ON)
  # find library
  FIND_PATH(SCOTCH_LIBRARY_PATH NAMES libscotch.a libscotcherr.a
  PATHS /usr/lib /usr/local/lib ENV LD_LIBRARY_PATH)
  IF(NOT SCOTCH_LIBRARY_PATH)
    MESSAGE(WARNING "SCOTCH library not found! Please set correct SCOTCH_LIBRARY_PATH manually!")
    SET(ISPARTMESH OFF)
  ENDIF(NOT SCOTCH_LIBRARY_PATH)
  # find include
  FIND_PATH(SCOTCH_INCLUDE_PATH NAMES scotch.h ptscotch.h
  PATHS /usr/include /usr/local/include ENV LD_INCLUDE_PATH)
  IF(NOT SCOTCH_INCLUDE_PATH)
    MESSAGE(WARNING "SCOTCH include not found! Please set correct SCOTCH_INCLUDE_PATH manually!")
    SET(ISPARTMESH OFF)
  ENDIF(NOT SCOTCH_INCLUDE_PATH)
ENDIF(BUILD_PARTMESH)
SET(BUILD_UTIL_EXODUSNEW2SEM OFF CACHE BOOL "Build exodusnew2sem program")
SET(BUILD_UTIL_EXODUSOLD2SEM OFF CACHE BOOL "Build exodusold2sem program")
SET(BUILD_UTIL_WRITE_SOS OFF CACHE BOOL "Build write_sos program")

# MPI for parallel processing
SET(ENABLE_MPI ON CACHE BOOL "Build semgeotech with MPI support")
IF(ENABLE_MPI)
  SET(ISPARALLEL ON)
  FIND_PACKAGE(MPI COMPONENTS Fortran REQUIRED)
  IF(NOT MPI_FOUND)
    MESSAGE(WARNING "MPI Fortran compiler not found! Parallel semgeotech cannot be built!")
    SET(ISPARALLEL OFF)
  ENDIF(NOT MPI_FOUND)
  
  # Test MPI Fortran compiler
  TRY_COMPILE(TEST_MPI ${CMAKE_BINARY_DIR}/CMakeTmp ${CMAKE_SOURCE_DIR}/src/testmpi.f90
  #COMPILE_DEFINITIONS USE_MPI
  #CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${MPI_Fortran_INCLUDE_PATH}"
  LINK_LIBRARIES ${MPI_Fortran_LIBRARIES}
  OUTPUT_VARIABLE test)
  IF(NOT TEST_MPI)
    MESSAGE(WARNING "MPI Fortran compiler: ${MPI_Fortran_COMPILER} doesn't work! Parallel semgeotech cannot be built! ${test}")
    SET(ISPARALLEL OFF)
  ENDIF(NOT TEST_MPI) 
ENDIF(ENABLE_MPI)

# to be used with cmake-modules
#FIND_PACKAGE(PETSc COMPONENTS Fortran REQUIRED)
#IF(NOT MPI_FOUND)
#  MESSAGE(WARNING "MPI Fortran compiler not found! Parallel semgeotech cannot be built!")
#  SET(ISPARALLEL OFF)
#ENDIF(NOT MPI_FOUND)
# PETSC library for solver
SET(ENABLE_PETSC ON CACHE BOOL "Build with PETSC support")
IF(ENABLE_PETSC)
  SET(ISPETSC ON)
  IF(EXISTS $ENV{PETSC_DIR})
    SET(PETSC_DIR $ENV{PETSC_DIR} CACHE PATH "PETSc directory.")
  ELSE()
    MESSAGE(WARNING "Please set correct PETSC_DIR environement variable to use PETSc!")
    SET(ISPETSC OFF)
  ENDIF()
  IF(EXISTS $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
    SET(PETSC_ARCH $ENV{PETSC_DIR}/$ENV{PETSC_ARCH} CACHE PATH "PETSC_ARCH directory.")
  ELSE()
    MESSAGE(WARNING "Please set correct PETSC_ARCH environement variable to use PETSc!")
    SET(ISPETSC OFF)
  ENDIF()
ENDIF(ENABLE_PETSC)

# process CMakeLists.txt in other subdirectories
ADD_SUBDIRECTORY(doc)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(util)

# notes
#SET_SOURCE_FILES_PROPERTIES(${SRC} PROPERTIES LANGUAGE "C")

