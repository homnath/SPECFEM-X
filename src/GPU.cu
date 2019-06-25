
#include "GPU.h"


__global__ void plus_reduce(realw *v1,realw*v2,int N,realw *total){

int tid =threadIdx.x;
int i =blockIdx.x*blockDim.x +threadIdx.x;

//Eachblockloadsitselementsintosharedmemory
__shared__ realw x[128];


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

  int nthreads =128;
  int nblocks = ceil(*size/nthreads ) + 1;
  plus_reduce<<<nblocks,nthreads>>>(d_v1,d_v2,*size,d_product); 
  cudaMemcpy(product,d_product,sizeof(realw),cudaMemcpyDeviceToHost);

   cudaFree(d_v1);
   cudaFree(d_v2);
}


extern "C"
void FC_FUNC_(prepare_gpu,
              PREPARE_GPU)(long* gpu_pointer, realw * h_K, int * NGLL3, int * size){

  Mesh* mp = (Mesh*) malloc( sizeof(Mesh) );
  *gpu_pointer = (long)mp;

  cudaMalloc((void**) &mp->K,(*NGLL3)*(*NGLL3)*(*size)*sizeof(realw));
  cudaMemcpy(mp->K,h_K,sizeof(realw)*(*size)*(*NGLL3)*(*NGLL3),cudaMemcpyHostToDevice);


}
