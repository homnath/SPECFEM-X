
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

realw * p_loc;
realw * kp_loc;

// Scalar values
realw * KSP_rtol;
realw * pkp;
realw * beta;
realw * rz;
//Cublas
cublasHandle_t cublas_handle;

} Mesh;

void start_timing_cuda(cudaEvent_t* start,cudaEvent_t* stop){
     cudaEventCreate(start);
       cudaEventCreate(stop);
         cudaEventRecord( *start, 0 );
         }
  
         /* ----------------------------------------------------------------------------------------------- */
  
         void stop_timing_cuda(cudaEvent_t* start,cudaEvent_t* stop,const char* info_str){
           float time;
  //           // stops events
               cudaEventRecord( *stop, 0 );
                 cudaEventSynchronize( *stop );
                   cudaEventElapsedTime( &time, *start, *stop );
                     cudaEventDestroy( *start );
                       cudaEventDestroy( *stop );
                         // user output
                           printf("%s: Execution Time = %f ms\n",info_str,time);
                           }
  
                           /* ----------------------------------------------------------------------------------------------- */
  
                           void stop_timing_cuda(cudaEvent_t* start,cudaEvent_t* stop,const char* info_str,float* t){
                             float time;
                               // stops events
                                 cudaEventRecord( *stop, 0 );
                                   cudaEventSynchronize( *stop );
                                    cudaEventElapsedTime( &time, *start, *stop );
                                       cudaEventDestroy( *start );
                                         cudaEventDestroy( *stop );
                                           // user output
                                             printf("%s: Execution Time = %f ms\n",info_str,time);
  
                                               // returns time
                                                 *t = time;
                                                 }
