#include "GPU.h"

extern "C"
void FC_FUNC_(prepare_gpu,
              PREPARE_GPU)(long* gpu_pointer, realw * h_K, int * nedof, int * nelmt, int * h_gdof_elmt, int * neq,realw * f, realw* dprecon, realw* u, realw * r, realw * p,realw*z, realw * KSP_rtol, int * myrank, int * NPROC, int * max_point, int * ngpart, int * mpi_count,int * mpi_gdof){

  Mesh* mp = (Mesh*) malloc( sizeof(Mesh) );
  *gpu_pointer = (long)mp;

  mp->nelmt = *nelmt ;
  mp->neq = *neq ;
  mp->nedof = *nedof;
  mp->NPROC = *NPROC;

  cublasCreate(&mp->cublas_handle);

  cudaMalloc((void**) &mp->K,(*nedof)*(*nedof)*(*nelmt)*sizeof(realw));
  cudaMemcpy(mp->K,h_K,sizeof(realw)*(*nelmt)*(*nedof)*(*nedof),cudaMemcpyHostToDevice);

  cudaMalloc((void**) &mp->gdof_elmt,(*nedof)*(*nelmt)*sizeof(int));
  cudaMemcpy(mp->gdof_elmt,h_gdof_elmt,sizeof(int)*(*nelmt)*(*nedof),cudaMemcpyHostToDevice);

  cudaMalloc((void**) &mp->kp,(*neq + 1)*sizeof(realw));

  int nthreads = 256;
  int nblocks = ceil((mp->neq + 1)/nthreads) + 1 ;
  cudaMemset(mp->kp,0,(mp->neq + 1)*sizeof(realw));

  cudaMalloc((void**) &mp->u,(*neq + 1)*sizeof(realw));
  cudaMemcpy(mp->u,u,sizeof(realw)*(*neq+1),cudaMemcpyHostToDevice);

  cudaMalloc((void**) &mp->f,(*neq + 1)*sizeof(realw));
  cudaMemcpy(mp->f,f,sizeof(realw)*(*neq+1),cudaMemcpyHostToDevice);

  cudaMalloc((void**) &mp->dprecon,(*neq + 1)*sizeof(realw));
  cudaMemcpy(mp->dprecon,dprecon,sizeof(realw)*(*neq+1),cudaMemcpyHostToDevice);

  cudaMalloc((void**) &mp->r,(*neq + 1)*sizeof(realw));
  cudaMemcpy(mp->r,r,sizeof(realw)*(*neq+1),cudaMemcpyHostToDevice);

  cudaMalloc((void**) &mp->z,(*neq + 1)*sizeof(realw));
  cudaMemcpy(mp->z,z,sizeof(realw)*(*neq+1),cudaMemcpyHostToDevice);

  cudaMalloc((void**) &mp->p,(*neq + 1)*sizeof(realw));
  cudaMemcpy(mp->p,p,sizeof(realw)*(*neq+1),cudaMemcpyHostToDevice);

  cudaMalloc((void**) &mp->kp,(*neq + 1)*sizeof(realw));
  cudaMalloc((void**) &mp->pkp,sizeof(realw));

  cudaMalloc((void**) &mp->p_loc,(*nedof*(*nelmt))*sizeof(realw));
  cudaMalloc((void**) &mp->kp_loc,(*nedof*(*nelmt))*sizeof(realw));

  mp->KSP_rtol = *KSP_rtol ;

  cudaMalloc((void**) &mp->rz,sizeof(realw));
  cudaMalloc((void**) &mp->rz_new,sizeof(realw));
  cudaMalloc((void**) &mp->beta,sizeof(realw));
  cudaMalloc((void**) &mp->pkp,sizeof(realw));
  cudaMalloc((void**) &mp->alpha,sizeof(realw));

  if (*NPROC > 1) {
    mp->ngpart = *ngpart; 
    mp->max_point = *max_point; 
    cudaMalloc((void**) &mp->MPI_send_recv_buffer,(*ngpart)*(*max_point)*sizeof(realw));

    cudaMalloc((void**) &mp->mpi_count,(*ngpart)*sizeof(int));
    cudaMemcpy(mp->mpi_count,mpi_count,(*ngpart)*sizeof(int),cudaMemcpyHostToDevice);

    cudaMalloc((void**) &mp->mpi_gdof,(*ngpart)*(*max_point)*sizeof(int));
    cudaMemcpy(mp->mpi_gdof,mpi_gdof,(*ngpart)*(*max_point)*sizeof(int),cudaMemcpyHostToDevice);
  }

  double free_db,used_db,total_db;
  get_free_memory(&free_db,&used_db,&total_db);
  printf("memory usage:\n");
  printf("  rank %d: GPU memory usage: used = %f MB, free = %f MB, total = %f MB\n",*myrank,
              used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);

}
