#include "GPU.h"

/* ----------------------------------------------------------------------------------------------- */

// GPU initialization

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(initialize_cuda_device,
              INITIALIZE_CUDA_DEVICE)(int* myrank_f,int* ncuda_devices) {

  int device;
  int device_count;

  // Gets rank number of MPI process
  int myrank = *myrank_f;

  // note: from here on we use the runtime API  ...

  // Gets number of GPU devices
  device_count = 0;
  cudaGetDeviceCount(&device_count);
  // Do not check if command failed with `exit_on_cuda_error` since it calls cudaDevice()/ThreadSynchronize():
  // If multiple MPI tasks access multiple GPUs per node, they will try to synchronize
  // GPU 0 and depending on the order of the calls, an error will be raised
  // when setting the device number. If MPS is enabled, some GPUs will silently not be used.
  //
  // being verbose and catches error from first call to CUDA runtime function, without synchronize call
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess){
    fprintf(stderr,"Error after cudaGetDeviceCount: %s\n", cudaGetErrorString(err));
  }

  // returns device count to fortran
  if (device_count == 0) printf("CUDA runtime error: there is no device supporting CUDA\n");
  *ncuda_devices = device_count;

  // Sets the active device
  if (device_count >= 1) {
    // generalized for more GPUs per node
    // note: without previous context release, cudaSetDevice will complain with the cuda error
    //         "setting the device when a process is active is not allowed"

    // releases previous contexts
    cudaDeviceReset();

    //printf("rank %d: cuda device count = %d sets device = %d \n",myrank,device_count,myrank % device_count);
    //MPI_Barrier(MPI_COMM_WORLD);

    // sets active device
    // device changes for different mpi processes according to number of device per node
    // (assumes that number of devices per node is the same for different compute nodes)
    device = myrank % device_count;

    cudaSetDevice( device );
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess){
      fprintf(stderr,"Error %s\n", cudaGetErrorString(err));
    }

    // double check that device was  properly selected
    cudaGetDevice(&device);
    if (device != (myrank % device_count)) {
       printf("error rank: %d devices: %d \n",myrank,device_count);
       printf("  cudaSetDevice()=%d\n  cudaGetDevice()=%d\n",myrank%device_count,device);
    }
  }

  // returns a handle to the active device
  cudaGetDevice(&device);

  // get device properties
  struct cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp,device);

  // exit if the machine has no CUDA-enabled device
  if (deviceProp.major == 9999 && deviceProp.minor == 9999){
    fprintf(stderr,"No CUDA-enabled device found, exiting...\n\n");
  }

  // outputs device infos to file
  char filename[BUFSIZ];
  FILE* fp;
  int do_output_info;

  // by default, only master process outputs device infos to avoid file cluttering
  do_output_info = 0;
  if (myrank == 0) {
    do_output_info = 1;
    sprintf(filename,"./gpu_device_info.txt");
  }

  // output to file
  if (do_output_info) {
    fp = fopen(filename,"w");
    if (fp != NULL){
      // display device properties
      fprintf(fp,"Device Name = %s\n",deviceProp.name);
      fprintf(fp,"memory:\n");
      fprintf(fp,"  totalGlobalMem (in MB): %f\n",(unsigned long) deviceProp.totalGlobalMem / (1024.f * 1024.f));
      fprintf(fp,"  totalGlobalMem (in GB): %f\n",(unsigned long) deviceProp.totalGlobalMem / (1024.f * 1024.f * 1024.f));
      fprintf(fp,"  totalConstMem (in bytes): %lu\n",(unsigned long) deviceProp.totalConstMem);
      fprintf(fp,"  Maximum 1D texture size (in bytes): %lu\n",(unsigned long) deviceProp.maxTexture1D);
      fprintf(fp,"  sharedMemPerBlock (in bytes): %lu\n",(unsigned long) deviceProp.sharedMemPerBlock);
      fprintf(fp,"  regsPerBlock (in bytes): %lu\n",(unsigned long) deviceProp.regsPerBlock);
      fprintf(fp,"blocks:\n");
      fprintf(fp,"  Maximum number of threads per block: %d\n",deviceProp.maxThreadsPerBlock);
      fprintf(fp,"  Maximum size of each dimension of a block: %d x %d x %d\n",
              deviceProp.maxThreadsDim[0],deviceProp.maxThreadsDim[1],deviceProp.maxThreadsDim[2]);
      fprintf(fp,"  Maximum sizes of each dimension of a grid: %d x %d x %d\n",
              deviceProp.maxGridSize[0],deviceProp.maxGridSize[1],deviceProp.maxGridSize[2]);
      fprintf(fp,"features:\n");
      fprintf(fp,"  Compute capability of the device = %d.%d\n", deviceProp.major, deviceProp.minor);
      fprintf(fp,"  multiProcessorCount: %d\n",deviceProp.multiProcessorCount);
      if (deviceProp.canMapHostMemory){
        fprintf(fp,"  canMapHostMemory: TRUE\n");
      }else{
        fprintf(fp,"  canMapHostMemory: FALSE\n");
      }
      if (deviceProp.deviceOverlap){
        fprintf(fp,"  deviceOverlap: TRUE\n");
      }else{
        fprintf(fp,"  deviceOverlap: FALSE\n");
      }
      if (deviceProp.concurrentKernels){
        fprintf(fp,"  concurrentKernels: TRUE\n");
      }else{
        fprintf(fp,"  concurrentKernels: FALSE\n");
      }
      // outputs initial memory infos via cudaMemGetInfo()
      double free_db,used_db,total_db;
      get_free_memory(&free_db,&used_db,&total_db);
      fprintf(fp,"memory usage:\n");
      fprintf(fp,"  rank %d: GPU memory usage: used = %f MB, free = %f MB, total = %f MB\n",myrank,
              used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);

      // closes output file
      fclose(fp);
    }
  }


}
