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

  int nthreads =128;
  int nblocks = ceil(*size/nthreads ) + 1;
 // vecdaxpy<<<nblocks,nthreads>>>(d_v1,d_v2,* scalar,*size); 
  cudaMemcpy(v1,d_v1,sizeof(realw)*(*size),cudaMemcpyDeviceToHost);

  cudaFree(d_v1);
  cudaFree(d_v2);
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


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

extern "C"
void FC_FUNC_(gpu_superloop,
              GPU_SUPERLOOP)(long* gpu_pointer, int * MAX_ITER, realw * h_u, int * h_iter, int * h_errcode){
  cudaEvent_t start,stop,start1,stop1;
  Mesh* mp = (Mesh*)(*gpu_pointer); //get mesh pointer out of fortran integer container

  int errcode = 1;
  int iter = 1 ;

  while(errcode != 0 && iter < *MAX_ITER){


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

  stop_timing_cuda(&start,&stop,"assemble",&time);
  start_timing_cuda(&start,&stop);

cublasDdot_v2(mp->cublas_handle, mp->neq+1, mp->r, 1, mp->z, 1, mp->rz);

cublasDdot_v2(mp->cublas_handle, mp->neq+1, mp->p, 1, mp->kp, 1, mp->pkp);


  stop_timing_cuda(&start,&stop,"dot products",&time);
  start_timing_cuda(&start,&stop);

vecdiv<<<1,1>>>( mp->rz,mp->pkp,mp->alpha,1);
nthreads =128;
nblocks = ceil((mp->neq+1)/nthreads ) + 1;

vecdaxpy<<<nblocks,nthreads>>>(mp->u,mp->p,mp->alpha,mp->neq+1);
cudaDeviceSynchronize();


  stop_timing_cuda(&start,&stop,"daxpy",&time);
  start_timing_cuda(&start,&stop);

realw alpha2;
realw max_p = 0;
realw max_u = 0;

int pidx = 0;
int uidx = 0;
cublasIdamax(mp->cublas_handle, mp->neq + 1, mp->p, 1, &pidx);
cublasIdamax(mp->cublas_handle, mp->neq + 1, mp->u, 1, &uidx);
cudaMemcpy(&max_p, mp->p + pidx - 1, sizeof(realw), cudaMemcpyDeviceToHost);
cudaMemcpy(&max_u, mp->u + uidx - 1, sizeof(realw), cudaMemcpyDeviceToHost);
max_p = abs(max_p);
max_u = abs(max_u);
//printf("max_p=%e %e, max_u=%lf %lf\n", max_p, max_p2, max_u, max_u2);


cudaMemcpy(&alpha2,mp->alpha,sizeof(realw),cudaMemcpyDeviceToHost);

if (abs(alpha2)*max_p/max_u <= mp->KSP_rtol) {
    errcode=0;
    break;
}


  stop_timing_cuda(&start,&stop,"max",&time);
  start_timing_cuda(&start,&stop);


vecdaxmy<<<nblocks,nthreads>>>(mp->r,mp->kp,mp->alpha,mp->neq+1);
vecmult<<<nblocks,nthreads>>>(mp->dprecon,mp->r,mp->z,mp->neq+1);

cublasDdot_v2(mp->cublas_handle, mp->neq+1, mp->r, 1, mp->z, 1, mp->rz_new);

vecdiv<<<1,1>>>( mp->rz_new,mp->rz,mp->beta,1);
vecadd2<<<nblocks,nthreads>>>(mp->p, mp->z,mp->beta,mp->neq+1);

  stop_timing_cuda(&start,&stop,"last loop part",&time);
  start_timing_cuda(&start,&stop);

  stop_timing_cuda(&start1,&stop1,"full loop",&time);
printf("");
//  stop_timing_cuda(&start,&stop,"memcpy",&time);
iter += 1 ; 
}

if ( errcode == 0 ) cudaMemcpy(h_u,mp->u,(mp->neq+1)*sizeof(realw),cudaMemcpyDeviceToHost);
*h_iter = iter;
*h_errcode = errcode;

}


extern "C"
void FC_FUNC_(prepare_gpu,
              PREPARE_GPU)(long* gpu_pointer, realw * h_K, int * nedof, int * nelmt, int * h_gdof_elmt, int * neq,realw * f, realw* dprecon, realw* u, realw * r, realw * p,realw*z, realw * KSP_rtol){

  Mesh* mp = (Mesh*) malloc( sizeof(Mesh) );
  *gpu_pointer = (long)mp;

  mp->nelmt = *nelmt ;
  mp->neq = *neq ;
  mp->nedof = *nedof;


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
}
