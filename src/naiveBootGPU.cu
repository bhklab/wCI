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

static const char *curandGetErrorString(curandStatus_t error)
{
	switch (error)
	{
		case CURAND_STATUS_SUCCESS:
			return "CURAND_STATUS_SUCCESS";

		case CURAND_STATUS_VERSION_MISMATCH:
			return "CURAND_STATUS_VERSION_MISMATCH";

		case CURAND_STATUS_NOT_INITIALIZED:
			return "CURAND_STATUS_NOT_INITIALIZED";

		case CURAND_STATUS_ALLOCATION_FAILED:
			return "CURAND_STATUS_ALLOCATION_FAILED";

		case CURAND_STATUS_TYPE_ERROR:
			return "CURAND_STATUS_TYPE_ERROR";

		case CURAND_STATUS_OUT_OF_RANGE:
			return "CURAND_STATUS_OUT_OF_RANGE";

		case CURAND_STATUS_LENGTH_NOT_MULTIPLE:
			return "CURAND_STATUS_LENGTH_NOT_MULTIPLE";

		case CURAND_STATUS_DOUBLE_PRECISION_REQUIRED:
			return "CURAND_STATUS_DOUBLE_PRECISION_REQUIRED";

		case CURAND_STATUS_LAUNCH_FAILURE:
			return "CURAND_STATUS_LAUNCH_FAILURE";

		case CURAND_STATUS_PREEXISTING_FAILURE:
			return "CURAND_STATUS_PREEXISTING_FAILURE";

		case CURAND_STATUS_INITIALIZATION_FAILED:
			return "CURAND_STATUS_INITIALIZATION_FAILED";

		case CURAND_STATUS_ARCH_MISMATCH:
			return "CURAND_STATUS_ARCH_MISMATCH";

		case CURAND_STATUS_INTERNAL_ERROR:
			return "CURAND_STATUS_INTERNAL_ERROR";
	}

	return "<unknown>";
}


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
	   if (code != cudaSuccess) 
		      {
			            fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
				          if (abort) exit(code);
					     }
}

#define gpuErrchkRand(ans) { gpuAssertRand((ans), __FILE__, __LINE__); }
inline void gpuAssertRand(curandStatus_t code, const char *file, int line, bool abort=true)
{
	if (code != CURAND_STATUS_SUCCESS)
		{
		        fprintf(stderr,"GPUassert: %s %s %d\n", curandGetErrorString(code), file, line);
			if (abort) exit(code);
	       	}
}

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


  gpuErrchk(cudaMalloc(&devrcimat, N*N*sizeof(double))); 
  gpuErrchk(cudaMalloc(&devOutVec, R*sizeof(double)));

  gpuErrchk(cudaMalloc(&devRandomNumbers, R*N*sizeof(double)));
  gpuErrchk(cudaMalloc(&permVector, R*N*sizeof(uint64_t)));


  gpuErrchk(cudaMemcpy(devrcimat, rcimat, N*N*sizeof(double), cudaMemcpyHostToDevice));

  gpuErrchkRand(curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT));

  cudaDeviceSynchronize();
  gpuErrchkRand(curandSetPseudoRandomGeneratorSeed(gen, *state));
  cudaDeviceSynchronize();

  gpuErrchkRand(curandGenerateUniformDouble(gen, devRandomNumbers, R*N));
  cudaDeviceSynchronize();
  // Creating permutation indicies from uniform doubles
  truncate_to_index<<<(R*N+(numThreads-1))/numThreads, numThreads>>>(devRandomNumbers, permVector, N, R*N);
  cudaDeviceSynchronize();

  cudaError_t error = cudaGetLastError();
  if(error != cudaSuccess) {
		        // print the CUDA error message and exit
	printf("CUDA error: %s\n", cudaGetErrorString(error));
	 exit(-1);
  }

  // Running one bootstrap instance per thread.  
  runBootOnDevice<<<(R+(numThreads-1))/numThreads, numThreads>>>(devrcimat, devOutVec, permVector, N, R);
  cudaDeviceSynchronize();

  error = cudaGetLastError();
  if(error != cudaSuccess) {
		        // print the CUDA error message and exit
	printf("CUDA error: %s\n", cudaGetErrorString(error));
	 exit(-1);
  }
  //Copying back results
  gpuErrchk(cudaMemcpy(outVec, devOutVec, R*sizeof(double), cudaMemcpyDeviceToHost));
  cudaDeviceSynchronize();

  // Freeing Memory
  gpuErrchk(cudaFree(permVector));
  curandDestroyGenerator(gen);
  gpuErrchk(cudaFree(devRandomNumbers));
  gpuErrchk(cudaFree(devOutVec));
  gpuErrchk(cudaFree(devrcimat));

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

