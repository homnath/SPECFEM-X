# Makefile for SPECFEM-X
# REVISION
#   HNG, Jul 12,2011; HNG, Apr 09,2010
# binary directory
BIN = ../bin
#include ${PETSC_DIR}/lib/petsc/conf/variables

FORTRAN_COMPILER = INTEL

ifeq ($(FORTRAN_COMPILER), INTEL)
	FC = ifort
	MPIFC = mpiifort #mpif90
	FFLAGS = -std08 -O3 -axAVX -xHost -mcmodel medium -shared-intel -check all -check bounds -check pointers -check stack -check uninit -fp-stack-check -traceback -ftz -funroll-loops -unroll5 -align sequence -qopt-report0 -qopt-report-phase=vec -implicitnone -warn truncated_source -warn unused -warn declarations -warn alignments -warn ignore_loc -warn interface -warn usage -heap-arrays

	CC = icc
	CFLAGS = -O3 -axAVX -vec-report0 -xHost
else
	FC = gfortran
	MPIFC = mpif90
	FFLAGS = -std=f2008 -O3 -mcmodel=medium -fimplicit-none -frange-check -fbounds-check -fmax-errors=10 -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -Warray-temporaries -Wimplicit-interface -fbacktrace -fdump-core -Wunused-parameter
	CC = gcc
	CFLAGS = -O3
endif


NVCC = nvcc
GENCODE_60 = -gencode=arch=compute_60,code=\"sm_60,compute_60\"
NVCC_FLAGS_BASE = $(CUDA_FLAGS) $(CUDA_INC) $(MPI_INCLUDES) $(COND_MPI_CPPFLAGS)
NVCC_FLAGS = $(NVCC_FLAGS_BASE) -DCUDA $(GENCODE_60) -ftz=true  -L/usr/local/cuda-10.1/lib64
CUDA_LINK =  -L/usr/local/cuda-10.1/lib64 -lcudart -lstdc++ -lcudadevrt -lcublas

## serial fortran compiler
##FC = gfortran
##FFLAGS = -g -std=f2003 -fbounds-check -fimplicit-none -O5 -pedantic -pedantic-errors -Wline-truncation -fno-trapping-math -fno-backslash
#
## parallel fortran compiler
#FC = mpif90 #~/openmpi-1.4.3/bin/mpif90 # OR path to appropriate MPI forrtan compiler
##FFLAGS = -g -std=f2008 -fbounds-check -fimplicit-none -O5 -pedantic -pedantic-errors -Wline-truncation -fno-trapping-math -fno-backslash
##FFLAGS = -std08 -xHost -O3 -mcmodel medium -shared-intel -check all -check bounds -fp-stack-check -C -g -traceback -ftz -funroll-loops -unroll5 -align sequence -implicitnone -warn truncated_source -warn declarations -warn alignments -warn ignore_loc -warn usage
#
#FFLAGS = -std08 -xHost -O3 -mcmodel medium -shared-intel -check all -check bounds -fp-stack-check -traceback -ftz -funroll-loops -unroll5 -align sequence -implicitnone -warn truncated_source -warn declarations -warn alignments -warn ignore_loc -warn usage -warn unused -gen-interfaces -warn interfaces
#debug: FFLAGS += -g
#FC_COMPILE = $(FC) $(FFLAGS)
#
## C compiler
#CC = icc
#CFLAGS = -O3 -Werror -Wall -pedantic-errors
#debug: CFLAGS += -g
#CC_COMPILE = $(CC) $(CFLAGS)

debug: FFLAGS += -g
FC_COMPILE = $(MPIFC) $(FFLAGS)
debug: CFLAGS += -g
CC_COMPILE = $(CC) $(CFLAGS)
LIB =
LIBMPI = -lmpi #-L/~/openmpi-1.4.3/lib/lmpi # OR path to appropriate MPI library
#LIBSCOTCH = -L/home/hgharti/lsoft/scotch_5.1.12/lib -lscotch -lscotcherr
LIBSCOTCH = -L${SCOTCH_DIR}/lib -lscotch -lscotcherr
INCSCOTCH = -I${SCOTCH_DIR}/include

#INCPETSC = -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include

# object tree
OBJ = $(BIN)/global.o            \
      $(BIN)/package_version.o   \
      $(BIN)/string_library.o    \
      $(BIN)/math_library.o      \
      $(BIN)/shape_library.o     \
      $(BIN)/element.o           \
      $(BIN)/degrees_of_freedom.o\
      $(BIN)/elastic.o           \
      $(BIN)/infinite_layer.o    \
      $(BIN)/recreate_fault.o    \
      $(BIN)/read_input.o        \
      $(BIN)/gll_library.o       \
      $(BIN)/mesh_spec.o         \
      $(BIN)/integration.o       \
      $(BIN)/infinite_element.o  \
      $(BIN)/utmgeo.o            \
      $(BIN)/visual.o            \
      $(BIN)/fault.o             \
      $(BIN)/weakform.o          \
      $(BIN)/prepare_gravity.o   \
      $(BIN)/viscoelastic.o      \
      $(BIN)/preprocess.o        \
      $(BIN)/gindex.o		         \
      $(BIN)/map_location.o      \
      $(BIN)/free_surface.o      \
      $(BIN)/apply_bc.o          \
      $(BIN)/apply_traction.o    \
      $(BIN)/apply_mtraction.o   \
      $(BIN)/plastic_library.o   \
      $(BIN)/postprocess.o       \
      $(BIN)/okada_solution.o    \
      $(BIN)/GPU_broken_loop.cuda.o          \
      $(BIN)/initialize_cuda.cuda.o  \
      $(BIN)/helper_functions.cuda.o  \
      $(BIN)/prepare_GPU.cuda.o  \
	    $(EMPTY_MACRO)

# object tree for SEM serial
OBJ_SEM = $(BIN)/serial_library.o \
          $(BIN)/math_library_serial.o   \
	        $(BIN)/prepare_sparse_serial.o	\
          $(BIN)/solver.o         \
          $(BIN)/solver_petsc.o   \
	        $(EMPTY_MACRO)

# object tree for SEMMPI parallel
OBJ_SEMMPI = $(BIN)/set_precision_mpi.o \
             $(BIN)/mpi_library.o       \
             $(BIN)/math_library_mpi.o  \
             $(BIN)/ghost_library_mpi.o \
	           $(BIN)/prepare_sparse.o	\
             $(BIN)/parsolver.o	\
	    			 $(BIN)/parsolver_petsc.o \
	           $(EMPTY_MACRO)

# object for both serial and MPI
OBJ_BOTH_SERIAL = $(BIN)/model.serial.o           \
                  $(BIN)/cmtsolution.serial.o     \
                  $(BIN)/earthquake_load.serial.o \
                  $(BIN)/station.serial.o         \
                  $(BIN)/benchmark.serial.o       \
	                $(EMPTY_MACRO)
OBJ_BOTH_PARALLEL = $(BIN)/model.parallel.o            \
                    $(BIN)/cmtsolution.parallel.o     \
                    $(BIN)/earthquake_load.parallel.o \
                    $(BIN)/station.parallel.o \
                    $(BIN)/benchmark.parallel.o       \
	                  $(EMPTY_MACRO)

# object tree for PARTMESH
OBJ_PARTMESH = $(BIN)/global.o           \
               $(BIN)/string_library.o   \
							 $(BIN)/serial_library.o   \
							 $(BIN)/math_library.o     \
               $(BIN)/element.o          \
               $(BIN)/infinite_layer.o \
               $(BIN)/partmesh_library.o \
							 $(BIN)/recreate_fault.o   \
							 $(BIN)/read_input.o       \
	             $(BIN)/partmesh_scotch.o  \
               $(BIN)/partmesh.o         \
	             $(EMPTY_MACRO)

OBJ_specfem-x =  $(OBJ)               \
                $(OBJ_SEM)           \
							  $(OBJ_BOTH_SERIAL)   \
							  $(BIN)/specfem3d.o  \
							  $(BIN)/specfemx.o    \
	              $(EMPTY_MACRO)
OBJ_pspecfem-x = $(OBJ)               \
                $(OBJ_SEMMPI)        \
								$(OBJ_BOTH_PARALLEL) \
								$(BIN)/pspecfem3d.o \
								$(BIN)/pspecfemx.o   \
	              $(EMPTY_MACRO)

# main executable for spectral-element method
specfem-x = $(BIN)/specfem-x
pspecfem-x = $(BIN)/pspecfem-x
PARTMESH = $(BIN)/partmesh

default: specfem-x pspecfem-x

clean:
	rm -f *.mod $(BIN)/*.o core src/*.*.swp
test:
	echo ${PETSC_FC_INCLUDES}
all: clean partmesh specfem-x pspecfem-x

debug: clean partmesh specfem-x pspecfem-x 

specfem-x: $(OBJ_specfem-x)
	$(FC_COMPILE) -o $(specfem-x) $(OBJ_specfem-x) $(LIB)  $(CUDA_LINK)

pspecfem-x: $(OBJ_pspecfem-x)
	$(FC_COMPILE) -o $(pspecfem-x) $(OBJ_pspecfem-x) $(LIBMPI) ${PETSC_KSP_LIB}

partmesh: $(OBJ_PARTMESH)
	$(FC_COMPILE) -o ${PARTMESH} $(OBJ_PARTMESH) $(LIBSCOTCH)

# compile file which can be compiled as both serial and parllel versions
# serial
$(BIN)/%.serial.o: %.F90
	${FC_COMPILE} -c -o $@ $<
# parallel
$(BIN)/%.parallel.o: %.F90
	${FC_COMPILE} -DUSE_MPI -c -o $@ $<

# compile instruction for each source
$(BIN)/set_precision_mpi.o: set_precision_mpi.f90
	$(FC_COMPILE) -c -o $(BIN)/set_precision_mpi.o set_precision_mpi.f90

$(BIN)/global.o: global.f90
	$(FC_COMPILE) -c -o $(BIN)/global.o global.f90

$(BIN)/package_version.o: package_version.f90
	$(FC_COMPILE) -c -o $(BIN)/package_version.o package_version.f90

$(BIN)/string_library.o: string_library.f90
	$(FC_COMPILE) -c -o $(BIN)/string_library.o string_library.f90

$(BIN)/math_library.o: math_library.f90
	$(FC_COMPILE) -c -o $(BIN)/math_library.o math_library.f90

$(BIN)/mpi_library.o: mpi_library.f90
	$(FC_COMPILE) -c -o $(BIN)/mpi_library.o mpi_library.f90

$(BIN)/serial_library.o: serial_library.f90
	$(FC_COMPILE) -c -o $(BIN)/serial_library.o serial_library.f90

$(BIN)/math_library_serial.o: math_library_serial.f90
	$(FC_COMPILE) -c -o $(BIN)/math_library_serial.o math_library_serial.f90

$(BIN)/math_library_mpi.o: math_library_mpi.f90
	$(FC_COMPILE) -c -o $(BIN)/math_library_mpi.o math_library_mpi.f90

$(BIN)/recreate_fault.o: recreate_fault.f90
	$(FC_COMPILE) -c -o $(BIN)/recreate_fault.o recreate_fault.f90

$(BIN)/read_input.o: read_input.f90
	$(FC_COMPILE) -c -o $(BIN)/read_input.o read_input.f90

$(BIN)/gll_library.o: gll_library.f90
	$(FC_COMPILE) -c -o $(BIN)/gll_library.o gll_library.f90

$(BIN)/mesh_spec.o: mesh_spec.f90
	$(FC_COMPILE) -c -o $(BIN)/mesh_spec.o mesh_spec.f90

$(BIN)/integration.o: integration.f90
	$(FC_COMPILE) -c -o $(BIN)/integration.o integration.f90

$(BIN)/infinite_layer.o: infinite_layer.f90
	$(FC_COMPILE) -c -o $(BIN)/infinite_layer.o infinite_layer.f90

$(BIN)/infinite_element.o: infinite_element.f90
	$(FC_COMPILE) -c -o $(BIN)/infinite_element.o infinite_element.f90

$(BIN)/utmgeo.o: utmgeo.f90
	$(FC_COMPILE) -c -o $(BIN)/utmgeo.o utmgeo.f90

$(BIN)/shape_library.o: shape_library.f90
	$(FC_COMPILE) -c -o $(BIN)/shape_library.o shape_library.f90

$(BIN)/element.o: element.f90
	$(FC_COMPILE) -c -o $(BIN)/element.o element.f90

$(BIN)/degrees_of_freedom.o: degrees_of_freedom.f90
	$(FC_COMPILE) -c -o $(BIN)/degrees_of_freedom.o degrees_of_freedom.f90

$(BIN)/elastic.o: elastic.f90
	$(FC_COMPILE) -c -o $(BIN)/elastic.o elastic.f90

$(BIN)/gauss_library.o: gauss_library.f90
	$(FC_COMPILE) -c -o $(BIN)/gauss_library.o gauss_library.f90

$(BIN)/fault.o: fault.f90
	$(FC_COMPILE) -c -o $(BIN)/fault.o fault.f90

$(BIN)/apply_bc.o: apply_bc.f90
	$(FC_COMPILE) -c -o $(BIN)/apply_bc.o apply_bc.f90

$(BIN)/prepare_gravity.o: prepare_gravity.f90
	$(FC_COMPILE) -c -o $(BIN)/prepare_gravity.o prepare_gravity.f90

$(BIN)/viscoelastic.o: viscoelastic.f90
	$(FC_COMPILE) -c -o $(BIN)/viscoelastic.o viscoelastic.f90

$(BIN)/preprocess.o: preprocess.f90
	$(FC_COMPILE) -c -o $(BIN)/preprocess.o preprocess.f90

$(BIN)/weakform.o: weakform.f90
	$(FC_COMPILE) -c -o $(BIN)/weakform.o weakform.f90

$(BIN)/map_location.o: map_location.f90
	$(FC_COMPILE) -c -o $(BIN)/map_location.o map_location.f90

$(BIN)/free_surface.o: free_surface.f90
	$(FC_COMPILE) -c -o $(BIN)/free_surface.o free_surface.f90

$(BIN)/apply_traction.o: apply_traction.f90
	$(FC_COMPILE) -c -o $(BIN)/apply_traction.o apply_traction.f90

$(BIN)/apply_mtraction.o: apply_mtraction.f90
	$(FC_COMPILE) -c -o $(BIN)/apply_mtraction.o apply_mtraction.f90

$(BIN)/model.o: model.F90
	$(FC_COMPILE) -DUSE_MPI -c -o $(BIN)/model.o model.F90

$(BIN)/cmtsolution.o: cmtsolution.F90
	$(FC_COMPILE) -DUSE_MPI -c -o $(BIN)/cmtsolution.o cmtsolution.F90

$(BIN)/earthquake_load.o: earthquake_load.F90
	$(FC_COMPILE) -DUSE_MPI -c -o $(BIN)/earthquake_load.o earthquake_load.F90

$(BIN)/station.o: station.F90
	$(FC_COMPILE) -DUSE_MPI -c -o $(BIN)/station.o station.F90

$(BIN)/prepare_sparse_serial.o: prepare_sparse_serial.f90
	$(FC_COMPILE) -c -o $(BIN)/prepare_sparse_serial.o prepare_sparse_serial.f90

$(BIN)/prepare_sparse.o: prepare_sparse.f90
	$(FC_COMPILE) -c -o $(BIN)/prepare_sparse.o prepare_sparse.f90

$(BIN)/gindex.o: gindex.f90
	$(FC_COMPILE) -c -o $(BIN)/gindex.o gindex.f90

$(BIN)/ghost_library_mpi.o: ghost_library_mpi.f90
	$(FC_COMPILE) -c -o $(BIN)/ghost_library_mpi.o ghost_library_mpi.f90

$(BIN)/plastic_library.o: plastic_library.f90
	$(FC_COMPILE) -c -o $(BIN)/plastic_library.o plastic_library.f90

$(BIN)/okada_solution.o: okada_solution.f90
	$(FC_COMPILE) -c -o $(BIN)/okada_solution.o okada_solution.f90

$(BIN)/visual.o: visual.f90
	$(FC_COMPILE) -c -o $(BIN)/visual.o visual.f90

$(BIN)/postprocess.o: postprocess.f90
	$(FC_COMPILE) -c -o $(BIN)/postprocess.o postprocess.f90

$(BIN)/benchmark.o: benchmark.F90
	$(FC_COMPILE) -DUSE_MPI -c -o $(BIN)/benchmark.o benchmark.F90

$(BIN)/solver.o: solver.f90
	$(FC_COMPILE) -c -o $(BIN)/solver.o solver.f90

$(BIN)/solver_petsc.o: solver_petsc.f90
	$(FC_COMPILE) -c -o $(BIN)/solver_petsc.o solver_petsc.f90

$(BIN)/parsolver.o: parsolver.f90
	$(FC_COMPILE) -c -o $(BIN)/parsolver.o parsolver.f90

#$(BIN)/parsolver_petsc.o: parsolver_petsc.F90
#	$(FC_COMPILE) ${INCPETSC} -cpp -c -o $(BIN)/parsolver_petsc.o parsolver_petsc.F90 -L${PETSC_DIR}/lib/petsc/conf/variables

$(BIN)/parsolver_petsc.o: parsolver_petsc.F90
	$(FC_COMPILE) ${INCPETSC} -cpp -c -o $(BIN)/parsolver_petsc.o parsolver_petsc.F90

#	$(FC_COMPILE) ${PETSC_FC_INCLUDES} -cpp -c -o $(BIN)/parsolver_petsc.o parsolver_petsc.f90 ${PETSC_KSP_LIB}
# both of the above two work fine
#
# following objective files are compiled from the same file but with
# different compiler preprocessing directives options.
# -DUSE_MPI provides the flag to use the MPI library for the compilation.
# the flag -DUSE_MPI is equivalent to #define USE_MPI 1.
$(BIN)/specfem3d.o: specfem3d.F90
	$(FC_COMPILE) -c -o $(BIN)/specfem3d.o specfem3d.F90

#$(BIN)/pspecfem3d.o: specfem3d.F90
#	$(FC_COMPILE) -DUSE_MPI -c -o $(BIN)/pspecfem3d.o specfem3d.F90

$(BIN)/pspecfem3d.o: specfem3d.F90
	$(FC_COMPILE) ${INCPETSC} -DUSE_MPI -c -o $(BIN)/pspecfem3d.o specfem3d.F90

# following objective files are compiled from the same file but with
# different compiler preprocessing directives options.
$(BIN)/specfemx.o: specfemx.F90
	$(FC_COMPILE) -c -o $(BIN)/specfemx.o specfemx.F90

$(BIN)/pspecfemx.o: specfemx.F90
	$(FC_COMPILE) -DUSE_MPI -c -o $(BIN)/pspecfemx.o specfemx.F90

# mesh partition
$(BIN)/partmesh_library.o: partmesh_library.f90
	$(FC_COMPILE) -c -o $(BIN)/partmesh_library.o partmesh_library.f90

$(BIN)/partmesh_scotch.o: partmesh_scotch.f90 partmesh_library.f90
	$(FC_COMPILE) -c -o $(BIN)/partmesh_scotch.o partmesh_scotch.f90 $(INCSCOTCH)

$(BIN)/partmesh.o: partmesh.f90
	$(FC_COMPILE) -c -o $(BIN)/partmesh.o partmesh.f90

$(BIN)/%.cuda.o: %.cu GPU.h
	${NVCC} -c $< -o $@ $(NVCC_FLAGS)


