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


// daxpy like routines
__global__ void vecadd(realw *v1, realw *v2, realw c, int N) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    // vec1 = vec1 + c * vec2
    if (index < N) v1[index] = v1[index] + c * v2[index];
}

__global__ void vecsub(realw *v1, realw *v2, realw c, int N) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    // vec1 = vec1 - c * vec2
    if (index < N) v1[index] = v1[index] - c * v2[index];
}

__global__ void vecmult(realw *v1, realw *v2, realw *v3, int N) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    // vec1 = vec2  * vec3
    if (index < N) v3[index] = v1[index] + v2[index];
}

__global__ void vecadd2(realw *v1, realw *v2, realw c, int N) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    // vec1 = c * vec1 + vec2
    if (index < N) v1[index] = c * v1[index] + v2[index] ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////


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
void FC_FUNC_(gpu_daxpy_1,
              GPU_DAXPY_1)(long* gpu_pointer, realw * v1, realw * v2, realw * scalar, int * size) {
  
  // Takes in two vectors and a scalar and recomputes the first vector:
  // v1 = v1 + v2*scalar 
  
  realw * d_v1;
  cudaMalloc((void**) &d_v1, *size*sizeof(realw)); 
  cudaMemcpy(d_v1,v1,sizeof(realw)*(*size),cudaMemcpyHostToDevice);
  realw * d_v2;
  cudaMalloc((void**) &d_v2, *size*sizeof(realw));
  cudaMemcpy(d_v2,v2,sizeof(realw)*(*size),cudaMemcpyHostToDevice);
  realw * d_scalar;
  cudaMalloc((void**) &d_scalar, sizeof(realw));
  cudaMemcpy(d_scalar,scalar,sizeof(realw),cudaMemcpyHostToDevice);

  int nthreads =128;
  int nblocks = ceil(*size/nthreads ) + 1;
  vecadd<<<nblocks,nthreads>>>(d_v1,d_v2, *d_scalar, *size); 
  cudaMemcpy(v1,d_v1,sizeof(realw),cudaMemcpyDeviceToHost);

   cudaFree(d_v1);
   cudaFree(d_v2);
   cudaFree(d_scalar);
}




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

//kp[gdof_elmt[tid]] += kp_loc[tid];
atomicAdd(&kp[gdof_elmt[ivec * nedof + tid]],kp_loc[ivec * nedof + tid]);

}

__global__ void set_vec_zero(realw *kp, int N) {
  // 
  int index = threadIdx.x + blockIdx.x * blockDim.x; 
  if (index < N ){
    kp[index] = 0.0;
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

extern "C"
void FC_FUNC_(compute_matvec_prod,
              COMPUTE_MATVEC_PROD)(long* gpu_pointer, realw * h_p, realw * h_kp){

  Mesh* mp = (Mesh*)(*gpu_pointer); //get mesh pointer out of fortran integer container

  cudaMemcpy(mp->p,h_p,sizeof(realw)*(mp->neq+1),cudaMemcpyHostToDevice);

  int N = mp->nelmt ;





  realw * p_loc, * kp_loc;
  cudaMalloc((void**)&p_loc,N*(mp->nedof)*sizeof(realw));
  cudaMalloc((void**)&kp_loc,N*(mp->nedof)*sizeof(realw));
  //cudaMemset(&mp->kp,0,(mp->neq + 1)*sizeof(realw));
  int nthreads = 256;
  int nblocks = ceil((mp->neq + 1)/nthreads) + 1 ;
  set_vec_zero<<<nblocks, nthreads>>>(mp->kp, (mp->neq + 1));

  cublasHandle_t cublas_handle=NULL;
  cublasCreate(&cublas_handle);

  for (int ielm = 0 ; ielm < mp->nelmt ; ielm+=N){
 
   printf("in loop %d\n",ielm);

   int nthreads = mp->nedof;
   int nblock = N;

   //get_p_loc_vector<<<nblock,nthreads>>>(mp->gdof_elmt + ielm * mp->nedof, mp->p,p_loc);
   const double beta = 0.0;
   const double alpha = 1.0; 
   
   int &m = mp->nedof;
   int &n = N;
   int &k = mp->nedof;
   int &lda = m;
   int &ldb = k;
   int &ldc = m;
   realw * A = mp->K + ielm * mp->nedof * mp->nedof;
   double * &B = p_loc;
   double * &C = kp_loc;

   get_p_loc_vector<<<nblock,nthreads>>>(mp->gdof_elmt + ielm * mp->nedof, mp->p, B);
   cudaDeviceSynchronize();
   cublasDgemm_v2(cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
   cudaDeviceSynchronize();
   assemble_kp_vector<<<nblock,nthreads>>>(mp->gdof_elmt + ielm * mp->nedof, mp->kp, C);

   }

   printf("finished loop\n");

   cublasDestroy(cublas_handle);
   cudaFree(p_loc);
   cudaFree(kp_loc);

   cudaMemcpy(h_kp,mp->kp,sizeof(realw)*(mp->neq+1),cudaMemcpyDeviceToHost);

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

  int nthreads = 256;
  int nblocks = ceil((mp->neq + 1)/nthreads) + 1 ;
  set_vec_zero<<<nblocks, nthreads>>>(mp->kp, (mp->neq + 1));



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
