/*
// Fast bootstap for rCI using a naive matrix based approach.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <curand.h>
#include <R.h>
#include <Rinternals.h>
#include "cudaHelpers.h"


const int numThreads = 64;


// Code to create indicies properly from the uniform random numbers. 
__global__
void truncate_to_index(double *randomDoubles, uint64_t *randomInt, uint64_t N, uint64_t maxI){
   
  uint64_t i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i >= maxI){
   return;
  }
  randomInt[i] = (uint64_t)ceil(randomDoubles[i] * N) - 1;

}



__global__
void runBootOnDevice(double *rcimat, double *outVec, uint64_t *permVector, uint64_t N, uint64_t R){
  
  uint64_t i = blockIdx.x*blockDim.x + threadIdx.x;
  
  if(i >= R){
    return;
  }

  uint64_t *permIdx;

  double currCI;
  double curVal;
  double RS_numerator, RS_denominator;
  RS_numerator = 0;
  RS_denominator = 0;
  permIdx = permVector + i*N;
  
  for(uint64_t j = 0; j < N; j++){

      for(uint64_t k = 0; k < N; k++){
        curVal = rcimat[permIdx[j]*N + permIdx[k]];

        RS_numerator += (curVal * (double)(curVal > 0));
        RS_denominator += (double)(curVal != 0) * 2;
      }

    }

    currCI = (RS_numerator)/(RS_denominator);
    
    outVec[i] = currCI;

}


void bootOnCuda(double *rcimat, double *outVec, uint64_t R, uint64_t N, int xties, int yties, uint64_t *state){


  double *devrcimat, *devOutVec, *devRandomNumbers;

  curandGenerator_t gen;
  size_t free_mem, total_mem;
  size_t mem_needed_per_R;

  uint64_t *permVector;

  uint64_t RperLoop, Roffset, loopI;

  gpuErrchkRand(curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT));

  cudaDeviceSynchronize();
  gpuErrchkRand(curandSetPseudoRandomGeneratorSeed(gen, *state));
  cudaDeviceSynchronize();

  // Need to calculate here to make sure we don't run out of memory while using the GPU. 
  // Will do the computation in batches over R. 

  gpuErrchk(cudaMalloc(&devrcimat, N*N*sizeof(double))); 
  gpuErrchk(cudaMemcpy(devrcimat, rcimat, N*N*sizeof(double), cudaMemcpyHostToDevice));


  gpuErrchk(cudaMemGetInfo(&free_mem, &total_mem));



  mem_needed_per_R = (2*(N+1))*sizeof(double);

  free_mem = free_mem - 50*sizeof(double); // keeping some extra buffer space of 50 doubles for variables allocated in kernels

  RperLoop = min(free_mem / mem_needed_per_R, R);

  //printf("R per loop: %lld\n", RperLoop);

  for(loopI = 0; loopI < (R/RperLoop)+1; loopI++){
    //printf(" loopI: %lld", loopI);
    //printf(" R/RperLoop + 1: %lld", (R/RperLoop)+1);
    Roffset = loopI * RperLoop;
    //printf(" R offset: %lld \n", Roffset);

    if((RperLoop + Roffset) > R){
      RperLoop = min(R - Roffset, 0);
    }
    if(RperLoop == 0){
      break;
    }

    gpuErrchk(cudaMalloc(&devOutVec, RperLoop*sizeof(double)));

    gpuErrchk(cudaMalloc(&devRandomNumbers, RperLoop*N*sizeof(double)));
    gpuErrchk(cudaMalloc(&permVector, RperLoop*N*sizeof(uint64_t)));






    gpuErrchkRand(curandGenerateUniformDouble(gen, devRandomNumbers, RperLoop*N));
    cudaDeviceSynchronize();
  // Creating permutation indicies from uniform doubles
    truncate_to_index<<<(RperLoop*N+(numThreads-1))/numThreads, numThreads>>>(devRandomNumbers, permVector, N, RperLoop*N);
    cudaDeviceSynchronize();

    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess) {
            // print the CUDA error message and exit
      printf("CUDA error: %s\n", cudaGetErrorString(error));
      exit(-1);
    }

    // Running one bootstrap instance per thread.  
    runBootOnDevice<<<(RperLoop+(numThreads-1))/numThreads, numThreads>>>(devrcimat, devOutVec, permVector, N, RperLoop);
    cudaDeviceSynchronize();

    error = cudaGetLastError();
    if(error != cudaSuccess) {
            // print the CUDA error message and exit
      printf("CUDA error: %s\n", cudaGetErrorString(error));
      exit(-1);
    }
    //Copying back results
    gpuErrchk(cudaMemcpy(outVec+Roffset, devOutVec, RperLoop*sizeof(double), cudaMemcpyDeviceToHost));
    cudaDeviceSynchronize();
    
    gpuErrchk(cudaFree(permVector));
    gpuErrchk(cudaFree(devRandomNumbers));
    gpuErrchk(cudaFree(devOutVec));
  }

  

  // Freeing Memory
  
  gpuErrchk(cudaFree(devrcimat));
  curandDestroyGenerator(gen);

}

extern "C"
SEXP bootCUDA(SEXP prcimat,
             // SEXP pobsCI,
             SEXP pR,
             // SEXP pB,
             SEXP pn,
             SEXP pxties,
             SEXP pyties,
             SEXP pseed){
  
  double Ndouble = *REAL(pn);
  
  // int discard_x_ties = *INTEGER(pdiscard_x_ties);
  // int discard_y_ties = *INTEGER(pdiscard_y_ties);
  
  // double obsCI = *REAL(pobsCI);

  double Rdouble = *REAL(pR);  
  // double Bdouble = *REAL(pB);

  uint64_t N = (uint64_t) Ndouble;
  uint64_t R = (uint64_t) Rdouble;
  // uint64_t B = (uint64_t) Bdouble;

  int xties = *INTEGER(pxties);
  int yties = *INTEGER(pyties);
  
  SEXP pout = PROTECT(allocVector(REALSXP,R));
  
  // double *out = REAL(pout);
  
  double *seed = REAL(pseed);
  uint64_t *state = (uint64_t*) seed;

  // double *rcimat2 = malloc(N * N * sizeof(double));
  bootOnCuda(REAL(prcimat), REAL(pout), R, N, xties, yties, state);
  // printf("%f\n", out[0]);
  // rciBootWithCopy(REAL(prcimat), rcimat2, REAL(pout), R, N, xties, yties, state);


  // free(rcimat2);
  UNPROTECT(1);
  
  return pout;
  
}


