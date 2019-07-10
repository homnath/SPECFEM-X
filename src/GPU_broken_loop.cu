#include "GPU.h"


// daxpy like routines
__global__ void vecdaxpy(realw *v1, realw *v2, realw * c, int N) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    // vec1 = vec1 + c * vec2
    if (index < N) v1[index] += *c * v2[index];
}

__global__ void vecdaxmy(realw *v1, realw *v2, realw * c, int N) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    // vec1 = vec1 - c * vec2
    if (index < N) v1[index] -= *c * v2[index];
}


__global__ void vecdiv(realw *v1, realw *v2, realw *c, int N) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    // vec1 = vec1 - c * vec2
    if (index < N) c[index] = v1[index] / v2[index];
}

__global__ void vecmult(realw *v1, realw *v2, realw *v3, int N) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    // vec1 = vec2  * vec3
    if (index < N) v3[index] = v1[index] * v2[index];
}

__global__ void vecadd2(realw *v1, realw *v2, realw *c, int N) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    // vec1 = c * vec1 + vec2
    if (index < N) v1[index] = *c * v1[index] + v2[index] ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void get_p_loc_vector(int *gdof_elmt,realw *p,realw * p_loc){
  int tid =threadIdx.x;
  int ivec = blockIdx.x;
  int nedof = blockDim.x;
  p_loc[ivec * nedof + tid] = p[gdof_elmt[ivec * nedof + tid]];
}

__global__ void assemble_kp_vector(int *gdof_elmt,realw *kp,realw * kp_loc){
  int tid =threadIdx.x;
  int ivec = blockIdx.x;
  int nedof = blockDim.x;
  atomicAdd(&kp[gdof_elmt[ivec * nedof + tid]],kp_loc[ivec * nedof + tid]);
}

// prepares a device array with with all inter-element edge-nodes -- this
// is followed by a memcpy and MPI operations
__global__ void prepare_buffer_on_device(realw* kp,
                                         realw* send_buffer,
                                         const int npart,
                                         const int max_point,
                                         const int* mpi_count,
                                         const int* mpi_gdof) {

  int id = threadIdx.x + blockIdx.x*blockDim.x; int ientry,idof;

  for(int i_part=0; i_part < npart; i_part++) {
    if (id< mpi_count[i_part]) {
      // entry in interface array
      ientry = id + max_point * i_part;
      // global index in wavefield
      idof = mpi_gdof[ientry] - 1;
      send_buffer[ientry] = kp[idof];
    }
  }
}

__global__ void add_neighbors_contrib_on_device(realw* kp,
                                                realw* recv_buffer,
                                                const int npart,
                                                const int max_point,
                                                const int* mpi_count,
                                                const int* mpi_gdof) {

  int id = threadIdx.x + blockIdx.x*blockDim.x ; int ientry,idof;
  for( int i_part=0; i_part < npart; i_part++) {
    if (id< mpi_count[i_part]) {
      // entry in interface array
      ientry = id + max_point * i_part;
      // global index in wavefield
      idof = mpi_gdof[ientry] - 1;
      atomicAdd(&kp[idof],recv_buffer[ientry]);
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

extern "C"
void FC_FUNC_(gpu_loop1,
              GPU_LOOP1)(long* gpu_pointer,realw * h_send_buffer){

  cudaEvent_t start,stop,start1,stop1;
  Mesh* mp = (Mesh*)(*gpu_pointer); //get mesh pointer out of fortran integer container

  float time;
  start_timing_cuda(&start,&stop);
  start_timing_cuda(&start1,&stop1);
  int N = mp->nelmt ;

  cudaMemset(mp->kp,0,(mp->neq + 1)*sizeof(realw));
  int nthreads = 256;
  int nblocks = ceil((mp->neq + 1)/nthreads) + 1 ;
  stop_timing_cuda(&start,&stop,"after memset",&time);
  start_timing_cuda(&start,&stop);
   nthreads = mp->nedof;
  int nblock = N;

  const double beta = 0.0;
  const double alpha = 1.0; 
  
  int &m = mp->nedof;
  int n = 1;
  int &k = mp->nedof;
  int &lda = m;
  int &ldb = k;
  int &ldc = m;
  realw * A = mp->K;
  double * &B = mp->p_loc;
  double * &C = mp->kp_loc;
  stop_timing_cuda(&start,&stop,"in between",&time);
  // Cuda timing
  start_timing_cuda(&start,&stop);
  
  get_p_loc_vector<<<nblock,nthreads>>>(mp->gdof_elmt, mp->p, B);
  stop_timing_cuda(&start,&stop,"p_loc",&time);
 start_timing_cuda(&start,&stop);

  cublasDgemmStridedBatched(mp->cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, A, lda, mp->nedof * mp->nedof, B, ldb, mp->nedof, &beta, C, ldc, mp->nedof, mp->nelmt);

  stop_timing_cuda(&start,&stop,"cublas",&time);
 start_timing_cuda(&start,&stop);

  assemble_kp_vector<<<nblock,nthreads>>>(mp->gdof_elmt, mp->kp, C);

  if ( mp->NPROC > 1 ){
    nthreads = 32;
    nblocks = ceil((mp->max_point)/nthreads) + 1 ;

    prepare_buffer_on_device<<<nblocks,nthreads>>>(mp->kp,
                                                   mp->MPI_send_recv_buffer,
                                                   mp->ngpart,
                                                   mp->max_point,
                                                   mp->mpi_count,
                                                   mp->mpi_gdof);

    cudaMemcpy(h_send_buffer,mp->MPI_send_recv_buffer,
               mp->ngpart * mp->max_point *sizeof(realw),cudaMemcpyDeviceToHost);
  } 

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

extern "C"
void FC_FUNC_(gpu_loop2,
              GPU_LOOP2)(long* gpu_pointer,realw * h_max_p, realw* h_max_u, realw* h_alpha, realw * h_recv_buffer){

  Mesh* mp = (Mesh*)(*gpu_pointer); //get mesh pointer out of fortran integer container

  int nthreads, nblocks;

  if ( mp->NPROC > 1 ){

    cudaMemcpy(mp->MPI_send_recv_buffer,h_recv_buffer,
               mp->ngpart * mp->max_point *sizeof(realw),cudaMemcpyHostToDevice);

    nthreads = 32;
    nblocks = ceil((mp->max_point)/nthreads) + 1 ;

    add_neighbors_contrib_on_device<<<nblocks,nthreads>>>(mp->kp,
                                                          mp->MPI_send_recv_buffer,
                                                          mp->ngpart,
                                                          mp->max_point,
                                                          mp->mpi_count,
                                                          mp->mpi_gdof);
  }

  cublasDdot_v2(mp->cublas_handle, mp->neq+1, mp->r, 1, mp->z, 1, mp->rz);
  cublasDdot_v2(mp->cublas_handle, mp->neq+1, mp->p, 1, mp->kp, 1, mp->pkp);

  vecdiv<<<1,1>>>( mp->rz,mp->pkp,mp->alpha,1);
  nthreads =128;
  nblocks = ceil((mp->neq+1)/nthreads ) + 1;

  vecdaxpy<<<nblocks,nthreads>>>(mp->u,mp->p,mp->alpha,mp->neq+1);

  int pidx = 0;
  int uidx = 0;
  cublasIdamax(mp->cublas_handle, mp->neq + 1, mp->p, 1, &pidx);
  cublasIdamax(mp->cublas_handle, mp->neq + 1, mp->u, 1, &uidx);
  cudaMemcpy(h_max_p, mp->p + pidx - 1, sizeof(realw), cudaMemcpyDeviceToHost);
  cudaMemcpy(h_max_u, mp->u + uidx - 1, sizeof(realw), cudaMemcpyDeviceToHost);
  *h_max_p = abs(*h_max_p);
  *h_max_u = abs(*h_max_u);

  cudaMemcpy(h_alpha,mp->alpha,sizeof(realw),cudaMemcpyDeviceToHost);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

extern "C"
void FC_FUNC_(gpu_loop3,
              GPU_LOOP3)(long* gpu_pointer){

  Mesh* mp = (Mesh*)(*gpu_pointer); //get mesh pointer out of fortran integer container
  int nthreads =128;
  int nblocks = ceil((mp->neq+1)/nthreads ) + 1;

  vecdaxmy<<<nblocks,nthreads>>>(mp->r,mp->kp,mp->alpha,mp->neq+1);
  vecmult<<<nblocks,nthreads>>>(mp->dprecon,mp->r,mp->z,mp->neq+1);

  cublasDdot_v2(mp->cublas_handle, mp->neq+1, mp->r, 1, mp->z, 1, mp->rz_new);

  vecdiv<<<1,1>>>( mp->rz_new,mp->rz,mp->beta,1);
  vecadd2<<<nblocks,nthreads>>>(mp->p, mp->z,mp->beta,mp->neq+1);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

extern "C"
void FC_FUNC_(gpu_loop4,
              GPU_LOOP4)(long* gpu_pointer,realw * h_u){
  Mesh* mp = (Mesh*)(*gpu_pointer); //get mesh pointer out of fortran integer container
  cudaMemcpy(h_u,mp->u,(mp->neq+1)*sizeof(realw),cudaMemcpyDeviceToHost);
}
