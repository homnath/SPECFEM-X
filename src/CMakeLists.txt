# build main programs
# REVISION:
#   HNG, Mar 21,2017; HNG, Sep 06,2011;

# build sem library
SET(srcsem
  global.f90
  string_library.f90
  math_library.f90
  shape_library.f90
  element.f90
  degrees_of_freedom.f90
  elastic.f90
  infinite_layer.f90
  recreate_fault.f90
  read_input.f90
  gll_library.f90
  mesh_spec.f90
  integration.f90
  infinite_element.f90
  utmgeo.f90
  visual.f90
  fault.f90
  weakform.f90
  prepare_gravity.f90
  apply_bc.f90
  preprocess.f90
  gindex.f90
  apply_traction.f90
  apply_mtraction.f90
  map_location.f90
  free_surface.f90
  compute_pressure.f90
  viscoelastic.f90
  plastic_library.f90
  okada_solution.f90
  postprocess.f90)
ADD_LIBRARY(sem ${srcsem})

# build parallel library
IF(ISPARALLEL)
  SET(srcparallel
     set_precision_mpi.f90
     mpi_library.f90
     math_library_mpi.f90
     ghost_library_mpi.f90
     prepare_sparse.f90
     parsolver.f90
     parsolver_petsc.F90)
  ADD_LIBRARY(parallel ${srcparallel})
  ADD_DEPENDENCIES(parallel sem)
  TARGET_INCLUDE_DIRECTORIES(parallel PUBLIC "${PETSC_DIR}/include" "${PETSC_ARCH}/include")
  TARGET_LINK_LIBRARIES(parallel ${PETSC_ARCH}/lib/libpetsc.so)
# build serial library
ELSE(ISPARALLEL)
  SET(srcserial
     serial_library.f90
     math_library_serial.f90
     prepare_sparse_serial.f90
     solver.f90
     solver_petsc.f90)
  ADD_LIBRARY(serial ${srcserial})
  ADD_DEPENDENCIES(serial sem)
ENDIF(ISPARALLEL)


SET(srcboth
   model.F90
   cmtsolution.F90
   earthquake_load.F90
   station.F90
   benchmark.F90)
#compiler preprocessor directives
#ADD_DEFINITIONS(CMAKE_Fortran_FLAGS "-DUSE_MPI")
IF(ISPARALLEL)
  SET_PROPERTY(SOURCE ${srcboth} PROPERTY COMPILE_DEFINITIONS USE_MPI)
  #SET_PROPERTY(SOURCE ${srcgeotech} PROPERTY INCLUDE_DIRECTORIES ${MPI_Fortran_INCLUDE_PATH})
ENDIF(ISPARALLEL)
ADD_LIBRARY(both ${srcboth})
IF(ISPARALLEL)
  ADD_DEPENDENCIES(both sem parallel)
ELSE(ISPARALLEL)
  ADD_DEPENDENCIES(both sem serial)
ENDIF(ISPARALLEL)

# build geotech library
#SET(srcgeotech "")
#message("${srcgeotech}")
#IF(BUILD_EXCAVATION)
#  LIST(APPEND srcgeotech semexcav3d.f90)
#ENDIF(BUILD_EXCAVATION)
#IF(BUILD_SLOPE)
#  LIST(APPEND srcgeotech semslope3d.f90)
#ENDIF(BUILD_SLOPE)
SET(srcgeotech
   specfem3d.F90)
#compiler preprocessor directives
#ADD_DEFINITIONS(CMAKE_Fortran_FLAGS "-DUSE_MPI")
IF(ISPARALLEL)
  SET_PROPERTY(SOURCE ${srcgeotech} PROPERTY COMPILE_DEFINITIONS USE_MPI)
  #SET_PROPERTY(SOURCE ${srcgeotech} PROPERTY INCLUDE_DIRECTORIES ${MPI_Fortran_INCLUDE_PATH})
ENDIF(ISPARALLEL)
ADD_LIBRARY(geotech ${srcgeotech})
IF(ISPARALLEL)
  ADD_DEPENDENCIES(geotech sem parallel both)
ELSE(ISPARALLEL)
  ADD_DEPENDENCIES(geotech sem serial both)
ENDIF(ISPARALLEL)

# build specfem-x
#TARGET_LINK_LIBRARIES(geotech ${semall})
#SET_TARGET_PROPERTIES(geotech PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
IF(ISPARALLEL)
  SET_PROPERTY(SOURCE specfem-x.F90 PROPERTY COMPILE_DEFINITIONS USE_MPI)
  #SET_PROPERTY(SOURCE specfem-x.F90 PROPERTY INCLUDE_DIRECTORIES ${MPI_Fortran_INCLUDE_PATH})
  SET(LIBS geotech both parallel sem ${MPI_Fortran_LIBRARIES})
  SET(mainprog pspecfem-x)
ELSE(ISPARALLEL)
  SET(LIBS geotech both serial sem)
  SET(mainprog specfem-x)
ENDIF(ISPARALLEL)
ADD_EXECUTABLE(${mainprog}
              specfem-x.F90)
TARGET_LINK_LIBRARIES(${mainprog} ${LIBS})
#IF(ISPARALLEL)
#  SET(semall geotech parallel sem)
#ELSE(ISPARALLEL)
#  SET(semall geotech serial sem)
#ENDIF(ISPARALLEL)
#TARGET_LINK_LIBRARIES(specfem-x ${semall} ${MPI_Fortran_LIBRARIES})

# Build partmesh program
# partmesh program partitions the mesh for the parallel processing
IF(ISPARTMESH)
  SET(srcpartmesh
     global.f90
     string_library.f90
     serial_library.f90
     math_library.f90
     element.f90
     infinite_layer.f90
     partmesh_library.f90
     recreate_fault.f90
     read_input.f90
     partmesh_scotch.f90
     partmesh.f90
  )
  ADD_LIBRARY(partscotch ${srcpartmesh})
  TARGET_INCLUDE_DIRECTORIES(partscotch PUBLIC
${SCOTCH_INCLUDE_PATH})
  ADD_EXECUTABLE(partmesh
              partmesh.f90)
  TARGET_LINK_LIBRARIES(partmesh partscotch
${SCOTCH_LIBRARY_PATH}/libscotch.a ${SCOTCH_LIBRARY_PATH}/libscotcherr.a)
ENDIF(ISPARTMESH)

