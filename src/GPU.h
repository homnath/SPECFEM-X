#ifndef GPU_H
#define GPU_H


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

void start_timing_cuda(cudaEvent_t* start,cudaEvent_t* stop);
void stop_timing_cuda(cudaEvent_t* start,cudaEvent_t* stop,const char* info_str);
void stop_timing_cuda(cudaEvent_t* start,cudaEvent_t* stop,const char* info_str,float* t);
void get_free_memory(double* free_db, double* used_db, double* total_db);

typedef struct mesh_ {

int nelmt;
int neq;
int nedof;
int NPROC;

// MPI variables
realw * MPI_send_recv_buffer;
int max_point;
int ngpart;

// Main GPU arrays
realw * K;
int * gdof_elmt ;  
realw * kp;
realw * u;
realw * f;
realw * dprecon; 
realw* r;
realw * p;
realw* z;

realw * p_loc;
realw * kp_loc;

// Scalar values
realw  KSP_rtol;
realw * pkp;
realw * beta;
realw * rz;
realw * rz_new;
realw * alpha;

//Cublas
cublasHandle_t cublas_handle;

} Mesh;

#endif
