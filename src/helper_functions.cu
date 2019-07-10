#include "GPU.h"

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

void get_free_memory(double* free_db, double* used_db, double* total_db) {

  // gets memory usage in byte
  size_t free_byte ;
  size_t total_byte ;
  cudaError_t cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;
  if (cudaSuccess != cuda_status) {
    printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
  }

  *free_db = (double)free_byte ;
  *total_db = (double)total_byte ;
  *used_db = *total_db - *free_db ;
  return;
}

