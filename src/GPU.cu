
#include "GPU.h"


__global__ void plus_reduce(realw *v1,realw*v2,int N,realw *total){

int tid =threadIdx.x;
int i =blockIdx.x*blockDim.x +threadIdx.x;

//Eachblockloadsitselementsintosharedmemory
__shared__ realw x[128];


//x[tid]=v1[i]*v2[i];
x[tid]=(i<N)?v1[i]:0.0;
x[tid]*=(i<N)?v2[i]:0.0;


__syncthreads();

//Buildsummationtreeoverelements.
for (unsigned int s=1; s<128 ; s *= 2) {
   if (tid % (2*s) == 0) x[tid] += x[tid + s];
   __syncthreads();
    }

//Thread0addsthepartialsumtothetotalsum
if(tid ==0)atomicAdd(total,x[tid]);
}



extern "C"
void FC_FUNC_(gpu_dot_product,
              GPU_DOT_PRODUCT)(long* gpu_pointer, realw * v1, realw * v2, int * size,realw * product) {
  realw * d_v1;
  cudaMalloc((void**) &d_v1, *size*sizeof(realw)); 
  cudaMemcpy(d_v1,v1,sizeof(realw)*(*size),cudaMemcpyHostToDevice);
  realw * d_v2;
  cudaMalloc((void**) &d_v2, *size*sizeof(realw));
  cudaMemcpy(d_v2,v2,sizeof(realw)*(*size),cudaMemcpyHostToDevice);
  
  realw * d_product;
  cudaMalloc((void**) &d_product,sizeof(realw));
  cudaMemset(d_product,0,sizeof(realw));

  cublasHandle_t cublas_handle = NULL;
  cublasCreate(&cublas_handle);
  cublasDdot_v2(cublas_handle, *size, d_v1, 1, d_v2, 1, d_product);
  cublasDestroy(cublas_handle);
  // int nthreads =128;
  // int nblocks = ceil(*size/nthreads ) + 1;
  // plus_reduce<<<nblocks,nthreads>>>(d_v1,d_v2,*size,d_product); 
  cudaMemcpy(product,d_product,sizeof(realw),cudaMemcpyDeviceToHost);

   cudaFree(d_v1);
   cudaFree(d_v2);
}


extern "C"
void FC_FUNC_(prepare_gpu,
              PREPARE_GPU)(long* gpu_pointer, realw * h_K, int * nedof, int * nelmt, int * h_gdof_elmt, int * neq,realw * f, realw* dprecon, realw* u, realw * r, realw * p,realw * KSP_rtol){

  Mesh* mp = (Mesh*) malloc( sizeof(Mesh) );
  *gpu_pointer = (long)mp;

  mp->nelmt = *nelmt ;
  mp->neq = *neq ;
  mp->nedof = *nedof;

  cudaMalloc((void**) &mp->K,(*nedof)*(*nedof)*(*nelmt)*sizeof(realw));
  cudaMemcpy(mp->K,h_K,sizeof(realw)*(*nelmt)*(*nedof)*(*nedof),cudaMemcpyHostToDevice);

  cudaMalloc((void**) &mp->gdof_elmt,(*nedof)*(*nelmt)*sizeof(int));
  cudaMemcpy(mp->gdof_elmt,h_gdof_elmt,sizeof(int)*(*nelmt)*(*nedof),cudaMemcpyHostToDevice);

  cudaMalloc((void**) &mp->kp,(*neq + 1)*sizeof(realw));
  cudaMemset(&mp->kp,0,(*neq + 1)*sizeof(realw));

  cudaMalloc((void**) &mp->u,(*neq + 1)*sizeof(realw));
  cudaMemcpy(mp->u,u,sizeof(realw)*(*neq+1),cudaMemcpyHostToDevice);

  cudaMalloc((void**) &mp->f,(*neq + 1)*sizeof(realw));
  cudaMemcpy(mp->f,f,sizeof(realw)*(*neq+1),cudaMemcpyHostToDevice);

  cudaMalloc((void**) &mp->dprecon,(*neq + 1)*sizeof(realw));
  cudaMemcpy(mp->dprecon,dprecon,sizeof(realw)*(*neq+1),cudaMemcpyHostToDevice);

  cudaMalloc((void**) &mp->r,(*neq + 1)*sizeof(realw));
  cudaMemcpy(mp->r,r,sizeof(realw)*(*neq+1),cudaMemcpyHostToDevice);

  cudaMalloc((void**) &mp->p,(*neq + 1)*sizeof(realw));
  cudaMemcpy(mp->p,p,sizeof(realw)*(*neq+1),cudaMemcpyHostToDevice);

  cudaMalloc((void**) &mp->kp,(*neq + 1)*sizeof(realw));
  cudaMalloc((void**) &mp->pkp,sizeof(realw));

  cudaMalloc((void**) &mp->KSP_rtol,sizeof(realw));
  cudaMemcpy(mp->KSP_rtol,KSP_rtol,sizeof(realw),cudaMemcpyHostToDevice);

  cudaMalloc((void**) &mp->rz,sizeof(realw));
  cudaMalloc((void**) &mp->beta,sizeof(realw));
  cudaMalloc((void**) &mp->pkp,sizeof(realw));

}
