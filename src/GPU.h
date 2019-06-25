
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <cuda.h>
//#include <cuda_runtime.h>
#include <cublas_v2.h>
//
//#include "config.h"

#define FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define FC_FUNC_(name,NAME) name ## _

typedef double realw;

typedef struct mesh_ {

int nelmt;
int neq;
int nedof;

realw * K;
int * gdof_elmt ;  
realw* kp;
realw * u;
realw * f;
realw* dprecon; 
realw* r;
realw * p;


// Scalar values
realw * KSP_rtol;
realw * pkp;
realw * beta;
realw * rz;

} Mesh;

