/*
// Fast permutations for rCI using a naive matrix based approach.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <curand.h>
#include <R.h>
#include <Rinternals.h>

// #include "xoroshiro128+.h"


// #define x) do { if((x)!=cudaSuccess) { \
//     printf("Error at %s:%d\n",__FILE__,__LINE__);\
//     return EXIT_FAILURE;}} while(0)
// #define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
//     printf("Error at %s:%d\n",__FILE__,__LINE__);\
//     return EXIT_FAILURE;}} while(0)


const int numThreads = 256;


// Code to create indicies properly from the uniform random numbers. 
__global__
void truncate_to_index(double *randomDoubles, uint64_t *randomInt, uint64_t N, uint64_t maxI){
   
  uint64_t i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i >= maxI){
   return;
  }
  randomInt[i] = (uint64_t)ceil(randomDoubles[i] * N);

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


  uint64_t *permVector;


  cudaMalloc(&devrcimat, N*N*sizeof(double)); 
  cudaMalloc(&devOutVec, R*sizeof(double));

  cudaMalloc(&devRandomNumbers, R*N*sizeof(double));
  cudaMalloc(&permVector, R*N*sizeof(uint64_t));


  cudaMemcpy(devrcimat, rcimat, N*N*sizeof(double), cudaMemcpyHostToDevice);

  curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
  curandSetPseudoRandomGeneratorSeed(gen, *state);

  curandGenerateUniformDouble(gen, devRandomNumbers, R*N);

  // Creating permutation indicies from uniform doubles
  truncate_to_index<<<(R*N+(numThreads-1))/numThreads, numThreads>>>(devRandomNumbers, permVector, N, R*N);

  // Running one bootstrap instance per thread.  
  runBootOnDevice<<<(R+(numThreads-1))/numThreads, numThreads>>>(devrcimat, devOutVec, permVector, N, R);


  //Copying back results
  cudaMemcpy(outVec, devOutVec, R*sizeof(double), cudaMemcpyDeviceToHost));

  // Freeing Memory
  cudaFree(permVector);
  cudaFree(devRandomNumbers);
  cudaFree(devOutVec);
  cudaFree(devrcimat);

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


