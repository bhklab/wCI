
// Fast permutations for rCI using a naive matrix based approach.




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <cuda.h>
#include <curand.h>
#include <R.h>
#include <Rinternals.h>
#include "cudaHelpers.h"

const int numThreads = 16; 
//The way its set up now, since each thread needs to do more computation as the size of N increases, this number is optimal to be smaller 
// as N gets larger... can I work around that? 

// Code to create indicies properly from the uniform random numbers. 
// This could be optimized by operating on more units of memory per kernel invocation (striding)
__global__
void truncate_to_index_for_fisher_yates(float *randomDoubles, uint64_t *randomInt, uint64_t N, uint64_t maxI){
   
  uint64_t i = blockIdx.x*blockDim.x + threadIdx.x;
  uint64_t N_cur = (i % N) + 1;

  if(i >= maxI){
   return;
  }

  randomInt[i] = (uint64_t)floor(randomDoubles[i] * N_cur);


  // #if __CUDA_ARCH__>=200
  //   printf("Thread: %lld , Random Int: %lld \n", i, randomInt[i]);

  // #endif 


}

__global__
void permuteYVecAndSave(float *devYvec, float *devYPerms, uint64_t *permVector, uint64_t N, uint64_t R){


  uint64_t i_start = (blockIdx.x*blockDim.x + threadIdx.x)*N;
  

  if(i_start >= R*N){
    return;
  }

  // Finishing an inside out Fisher-Yates shuffle, the indicies are already in [j,N-1], 
  // so I just need to do the shuffle  
  // In this case permVector[i_start + i] is "j", and everything is offset by i_start except devYvec
  for(uint64_t i = 0; i < N; i++){

    if(permVector[i_start + i] != i){
      devYPerms[i_start + i] = devYPerms[i_start + permVector[i_start + i]];
    }

    devYPerms[i_start + permVector[i_start + i]] = devYvec[i];

  }



  // #if __CUDA_ARCH__>=200
  // for(uint64_t i = 0; i < N; i++){

  //   printf("Shuffled Y element %lld is %lf \n", i, devYPerms[i_start + i]);
  // }
  // #endif 

}


// create a binary bit vector to record the inversion or not of each pair in input vector
// Note that the GPU needs to do this in 32 bytes to be efficient 

// Implement this on GPU by: calculating first which R we are in, then how many of the pairs offset we are

// First start off with a DEV version
__device__
void createInvVecDev(float *ValueVec, unsigned int *invVec, uint64_t N){

  uint64_t i = 0;
  uint64_t j = 0;
  uint64_t totIndex = 0;
  uint64_t curInt = 0;

  uint64_t curBit = 0;

  for(i = 0; i < N; i++){

    for(j = 0; j < i; j++){

      curInt = totIndex/(sizeof(unsigned int) * CHAR_BIT);
      curBit = totIndex%(sizeof(unsigned int) * CHAR_BIT);

      if(ValueVec[j] > ValueVec[i]){
        invVec[curInt] |= (1 << (((sizeof(unsigned int) * CHAR_BIT) - 1 ) - curBit)); // Sending 1 to appropriate place to signify inverted
      }
      totIndex += 1;
    }
  }
}

__device__
void createValidVecDev(float *ValueVec, unsigned int *validVec, float Delta, uint64_t N){

  uint64_t i = 0;
  uint64_t j = 0;
  uint64_t totIndex = 0;
  uint64_t curInt = 0;

  uint64_t curBit = 0;

  for(i = 0; i < N; i++){

    for(j = 0; j < i; j++){

      curInt = totIndex/(sizeof(unsigned int) * CHAR_BIT);
      curBit = totIndex%(sizeof(unsigned int) * CHAR_BIT);

      if(fabs(ValueVec[j] - ValueVec[i]) >= Delta){
        validVec[curInt] |= (1 << (((sizeof(unsigned int) * CHAR_BIT) - 1 ) - curBit)); // Sending 1 to appropriate place to signify inverted
      }
      totIndex += 1;
    }
  }
}


__global__
void createYinvKernel(float *devYPerms, unsigned int *devYinvVec,uint64_t N, uint64_t R){

  uint64_t i = blockIdx.x*blockDim.x + threadIdx.x; // which permutation are we working on
  uint64_t numIntsNeeded = ((N * N - N)/(2*sizeof(unsigned int)*CHAR_BIT) + 1);


  if(i >= R){ // out of permutations to work on
    return;
  }

  unsigned int *curYinvVec = devYinvVec + i * numIntsNeeded;

  createInvVecDev(devYPerms + i * N, curYinvVec, N);

}


__global__
void createYvalidKernel(float *devYPerms, unsigned int *devYvalidVec, uint64_t N, float DeltaY, uint64_t R){

  uint64_t i = blockIdx.x*blockDim.x + threadIdx.x; // which permutation are we working on
  uint64_t numIntsNeeded = ((N * N - N)/(2*sizeof(unsigned int)*CHAR_BIT) + 1);


  if(i >= R){ // out of permutations to work on
    return;
  }

  unsigned int *curYValidVec = devYvalidVec + i * numIntsNeeded;

  createValidVecDev(devYPerms + i * N, curYValidVec, DeltaY, N);

}


// Split this up into a strided computation and a reduction?
__global__
void calculateCI(unsigned int *devXinvVec, unsigned int *devXvalidVec, 
                 unsigned int *devYinvVec, unsigned int *devYvalidVec,
                 float *outVec, uint64_t N, uint64_t R){
  
  uint64_t i = blockIdx.x*blockDim.x + threadIdx.x; // which permutation are we working on
  uint64_t numIntsNeeded = ((N * N - N)/(2*sizeof(unsigned int)*CHAR_BIT) + 1);

  uint64_t totalDisconcordant = 0;
  uint64_t totalValid = 0;

  unsigned int *curYinvVec = devYinvVec + i * numIntsNeeded;
  unsigned int *curYValidVec = devYvalidVec + i * numIntsNeeded;

  if(i >= R){ // out of permutations to work on
    return;
  }
  // printf("curYinvVec: %u\n", *curYinvVec);

  // printf("devYPerms 0: %f\n", *devYPerms);
  // printf("devYPerms 1: %f\n", *(devYPerms+1));

  // createInvVecDev(devYPerms + i * N, curYinvVec, N);
  // createValidVecDev(devYPerms + i* N, curYValidVec, DeltaY, N);

  // printf("curYinvVec: %u\n", *curYinvVec);
  // printf("curYValidVec: %u\n", *curYValidVec);

  // printf("devXinvVec: %u\n", *devXinvVec);
  // printf("devXvalidVec: %u\n", *devXvalidVec);



  // printf("numIntsNeeded: %u\n", numIntsNeeded);

  for(uint64_t j = 0; j < numIntsNeeded; j++){

    totalDisconcordant += __popc((curYinvVec[j]^*devXinvVec)&(curYValidVec[j]&*devXvalidVec));
    totalValid += __popc((curYValidVec[j]&*devXvalidVec));


  }
  // currCI = (RS_numerator)/(RS_denominator);
    // printf("totalDisconcordant: %lld\n", totalDisconcordant);
    // printf("totalValid: %lld\n", totalValid);

  outVec[i] = 1-(float)totalDisconcordant/(float)totalValid;
  // outVec[i] = 5.0;
}


/*
nvcc -O3 -arch=sm_86 -G -I/usr/local/cuda/include -I/usr/share/R/include -L/usr/lib/R/lib -lR -L/usr/local/cuda/lib64 -lcurand --shared -Xcompiler -fPIC -o naivePermGPU.so naivePermGPU.cu

dyn.load("naivePermGPU.so")
t <- .Call("permCUDA", as.numeric(1:100), runif(100), 1, 100, 0, 0, runif(1));
*/


void permOnCuda(unsigned int *XinvVec, unsigned int *XvalidVec, float *yvec, 
                float *outVec, uint64_t R, uint64_t N, float DeltaY, uint64_t *state){


  float *devYvec, *devOutVec, *devRandomNumbers; 
  // TODO: convert the double precision vectors to single for efficiency

  curandGenerator_t gen;
  cudaStream_t stream1, stream2;

  gpuErrchk(cudaStreamCreate(&stream1));
  gpuErrchk(cudaStreamCreate(&stream2));

  size_t free_mem, total_mem;
  size_t mem_needed_per_R;

  uint64_t *permVector;

  uint64_t RperLoop, Roffset, loopI;

  unsigned int *devXinvVec, *devXvalidVec, *devYinvVec, *devYvalidVec;

  uint64_t numIntsNeeded = ((N * N - N)/(2*sizeof(unsigned int)*CHAR_BIT) + 1);


  gpuErrchkRand(curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT));

  cudaDeviceSynchronize();
  gpuErrchkRand(curandSetPseudoRandomGeneratorSeed(gen, *state));
  cudaDeviceSynchronize();


  gpuErrchk(cudaMalloc(&devYvec, N*sizeof(float))); 
  gpuErrchk(cudaMemcpy(devYvec, yvec, N*sizeof(float), cudaMemcpyHostToDevice));

  gpuErrchk(cudaMalloc(&devXinvVec, numIntsNeeded * sizeof(unsigned int))); 
  gpuErrchk(cudaMalloc(&devXvalidVec, numIntsNeeded * sizeof(unsigned int))); 

  gpuErrchk(cudaMemcpy(devXinvVec, XinvVec,numIntsNeeded * sizeof(unsigned int), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(devXvalidVec, XvalidVec,numIntsNeeded * sizeof(unsigned int), cudaMemcpyHostToDevice));


  gpuErrchk(cudaMemGetInfo(&free_mem, &total_mem));




  // Need to calculate here to make sure we don't run out of memory while using the GPU. 
  // Will do the computation in batches over R. 

  mem_needed_per_R = (N+1)*sizeof(float) + N*sizeof(uint64_t) + 2 * (size_t)numIntsNeeded * sizeof(unsigned int);

  free_mem = free_mem - 1000*sizeof(double); // keeping some extra buffer space of 50 doubles for variables allocated in kernels

  RperLoop = (uint64_t)min(free_mem*0.95 / mem_needed_per_R, R);

  // printf("R per loop: %lld \n", RperLoop);
  // printf("R: %lld \n", R);
  // printf("Stopping condition: %lld \n", (R/RperLoop)+1);
  // printf("mem_needed_per_R: %lld \n", mem_needed_per_R);

  uint64_t stoppingCond = (R/RperLoop)+1;
  
/*
dyn.load("naivePermGPU.so")
system.time(t <- .Call("permCUDA", as.numeric(1:1000), runif(1000), 1e4, 1000, 0, 0, runif(1)))

system.time(t <- .Call("permCUDA", as.numeric(1:100), runif(100), 1e7, 100, 0, 0.2, runif(1)))
*/

  for(loopI = 0; loopI < stoppingCond; loopI++){ 

    // printf("Loop I: %lld \n", loopI);
    // printf("Stopping condition: %lld", (R/RperLoop)+1);

    Roffset = loopI * RperLoop;
    // printf("Roffset: %lld \n ", Roffset);

    if((RperLoop + Roffset) > R){
      RperLoop = R - Roffset; // this should never underflow because the loop should stop before that happens
    }
    // printf("Roffset: %lld \n", Roffset);

    if(RperLoop == 0){
      break;
    }

    gpuErrchk(cudaMalloc(&devOutVec, RperLoop*sizeof(float)));


    // TODO:: I can reuse the devRandomNumbers vector to store my permuted values!
    // TODO:: Make sure I am freeing everything I allocate 

    gpuErrchk(cudaMalloc(&devRandomNumbers, RperLoop*N*sizeof(float)));
    // This holds randomly sampled indicies for the Fisher Yates shuffle (inside out)
    gpuErrchk(cudaMalloc(&permVector, RperLoop*N*sizeof(uint64_t))); 

    // This will hold a permuted copy of Y for each vector
    // gpuErrchk(cudaMalloc(&devRandomNumbers, RperLoop*N*sizeof(float))); 

    // devRandomNumbers = devRandomNumbers;

    // Allocating the y inv and valid vectors for each R
    gpuErrchk(cudaMalloc(&devYinvVec, numIntsNeeded * sizeof(unsigned int) * RperLoop)); 
    gpuErrchk(cudaMalloc(&devYvalidVec, numIntsNeeded * sizeof(unsigned int) * RperLoop)); 

    // Initialize to 0 (required as only bits which need to be flipped will be flipped)

    gpuErrchk(cudaMemset(devYinvVec, 0, numIntsNeeded * sizeof(unsigned int) * RperLoop)); 
    gpuErrchk(cudaMemset(devYvalidVec, 0, numIntsNeeded * sizeof(unsigned int) * RperLoop));     

    gpuErrchkRand(curandGenerateUniform(gen, devRandomNumbers, RperLoop*N));

    // Creating permutation indicies from uniform doubles
    truncate_to_index_for_fisher_yates<<<(RperLoop*N+(numThreads-1))/numThreads, numThreads>>>(devRandomNumbers, permVector, N, RperLoop*N);

    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess) {
            // print the CUDA error message and exit
      printf("CUDA error: %s\n", cudaGetErrorString(error));
      exit(-1);
    }


    permuteYVecAndSave<<<(RperLoop+(numThreads-1))/numThreads, numThreads>>>(devYvec, devRandomNumbers, permVector, N, RperLoop);

    error = cudaGetLastError();
    if(error != cudaSuccess) {
            // print the CUDA error message and exit
      printf("CUDA error: %s\n", cudaGetErrorString(error));
      exit(-1);
    }


    // Think about running this in streams?
    createYinvKernel<<<(RperLoop+(numThreads-1))/numThreads, numThreads, 0, stream1>>>(devRandomNumbers, devYinvVec, N, RperLoop);
    createYvalidKernel<<<(RperLoop+(numThreads-1))/numThreads, numThreads, 0, stream2>>>(devRandomNumbers, devYvalidVec, N, DeltaY, RperLoop);

    // Running one CI calculation per thread.
    cudaDeviceSynchronize();  
    calculateCI<<<(RperLoop+(numThreads-1))/numThreads, numThreads>>>(devXinvVec, devXvalidVec, devYinvVec,
                                                                      devYvalidVec, devOutVec, N, RperLoop);

    error = cudaGetLastError();
    if(error != cudaSuccess) {
            // print the CUDA error message and exit
      printf("CUDA error: %s\n", cudaGetErrorString(error));
      exit(-1);
    }


    gpuErrchk(cudaMemcpy(outVec+(Roffset), devOutVec, RperLoop*sizeof(float), cudaMemcpyDeviceToHost));
    






    // error = cudaGetLastError();
    // if(error != cudaSuccess) {
    //         // print the CUDA error message and exit
    //   printf("CUDA error: %s\n", cudaGetErrorString(error));
    //   exit(-1);
    // }
    // //Copying back results
    // gpuErrchk(cudaMemcpy(outVec+Roffset, devOutVec, RperLoop*sizeof(double), cudaMemcpyDeviceToHost));
    // cudaDeviceSynchronize();
    
    gpuErrchk(cudaFree(permVector));
    gpuErrchk(cudaFree(devRandomNumbers));
    gpuErrchk(cudaFree(devOutVec));
    // gpuErrchk(cudaFree(devYPerms));
    gpuErrchk(cudaFree(devYinvVec));
    gpuErrchk(cudaFree(devYvalidVec));

  }

  

  // Freeing Memory
  
  gpuErrchk(cudaFree(devYvec));
  gpuErrchk(cudaFree(devXinvVec));
  gpuErrchk(cudaFree(devXvalidVec));

  curandDestroyGenerator(gen);

  gpuErrchk(cudaStreamDestroy(stream1));
  gpuErrchk(cudaStreamDestroy(stream2));


}
__host__
void createInvVec(double *ValueVec, unsigned int *invVec, uint64_t N){

  uint64_t i = 0;
  uint64_t j = 0;
  uint64_t totIndex = 0;
  uint64_t curInt = 0;

  uint64_t curBit = 0;

  for(i = 0; i < N; i++){

    for(j = 0; j < i; j++){

      curInt = totIndex/(sizeof(unsigned int) * CHAR_BIT);
      curBit = totIndex%(sizeof(unsigned int) * CHAR_BIT);

      if(ValueVec[j] > ValueVec[i]){
        invVec[curInt] |= (1 << (((sizeof(unsigned int) * CHAR_BIT) - 1 ) - curBit)); // Sending 1 to appropriate place to signify inverted
      }
      totIndex += 1;
    }
  }
}

__host__
void createValidVec(double *ValueVec, unsigned int *validVec, double Delta, uint64_t N){

  uint64_t i = 0;
  uint64_t j = 0;
  uint64_t totIndex = 0;
  uint64_t curInt = 0;

  uint64_t curBit = 0;

  for(i = 0; i < N; i++){

    for(j = 0; j < i; j++){

      curInt = totIndex/(sizeof(unsigned int) * CHAR_BIT);
      curBit = totIndex%(sizeof(unsigned int) * CHAR_BIT);

      if(fabs(ValueVec[j] - ValueVec[i]) >= Delta){
        validVec[curInt] |= (1 << (((sizeof(unsigned int) * CHAR_BIT) - 1 ) - curBit)); // Sending 1 to appropriate place to signify inverted
      }
      totIndex += 1;
    }
  }
}

//https://stackoverflow.com/questions/111928/is-there-a-printf-converter-to-print-in-binary-format
#define BYTE_TO_BINARY_PATTERN "%c%c%c%c%c%c%c%c"
#define BYTE_TO_BINARY(byte)  \
  (byte & 0x80 ? '1' : '0'), \
  (byte & 0x40 ? '1' : '0'), \
  (byte & 0x20 ? '1' : '0'), \
  (byte & 0x10 ? '1' : '0'), \
  (byte & 0x08 ? '1' : '0'), \
  (byte & 0x04 ? '1' : '0'), \
  (byte & 0x02 ? '1' : '0'), \
  (byte & 0x01 ? '1' : '0') 


void convDoubleToSingle(double *doubleVec, float *singleVec, uint64_t N){

  for(uint64_t i = 0; i < N; i++){
    singleVec[i] = (float)doubleVec[i];
  }

}

void convSingleToDouble( float *singleVec, double *doubleVec, uint64_t N){

  for(uint64_t i = 0; i < N; i++){
    doubleVec[i] = (double)singleVec[i];
  }

}
extern "C"
SEXP permCUDA(SEXP pxvec,
              SEXP pyvec,
             // SEXP pobsCI,
             SEXP pR,
             // SEXP pB,
             SEXP pn,
             SEXP pDeltaX,
             SEXP pDeltaY,
             SEXP pseed){
  
  double Ndouble = *REAL(pn);  
  // double obsCI = *REAL(pobsCI);
  double Rdouble = *REAL(pR);  

  double DeltaX = *REAL(pDeltaX);
  float DeltaY = (float) (*REAL(pDeltaY));

  uint64_t N = (uint64_t) Ndouble;
  uint64_t R = (uint64_t) Rdouble;
  // uint64_t B = (uint64_t) Bdouble;

  double *seed = REAL(pseed);
  uint64_t *state = (uint64_t*) seed;

  uint64_t numIntsNeededX = ((N * N - N)/(2*sizeof(unsigned int)*CHAR_BIT) + 1);
  // I am allocating 2 arrays with 1 bit for each pair in X. One array to hold whether a pair is valid, and one to hold if its inverted. 
  //consider aligning this to 32 bytes for the GPU to be happy? Alhough googling suggests this is not necessary 
  unsigned int *XinvVec = (unsigned int *)malloc(numIntsNeededX * sizeof(unsigned int)); //Sometimes allocating an extra int, but 32 bits shouldn't 
  memset(XinvVec, 0, numIntsNeededX * sizeof(unsigned int));
 
  unsigned int *XvalidVec = (unsigned int *)malloc(numIntsNeededX * sizeof(unsigned int)); //Sometimes allocating an extra int, but 32 bits shouldn't 
  memset(XvalidVec, 0, numIntsNeededX * sizeof(unsigned int));
  

  SEXP pout = PROTECT(allocVector(REALSXP,R));
  

  float *singleY = (float *)malloc(N *sizeof(float));
  float *singleOut = (float *)malloc(R * sizeof(float));

  convDoubleToSingle(REAL(pyvec), singleY, N);


  // double *out = REAL(pout);
  


  createInvVec(REAL(pxvec), XinvVec,  N);
  // printf("\nFirst Int of Inversion Vec: "BYTE_TO_BINARY_PATTERN" "BYTE_TO_BINARY_PATTERN" "BYTE_TO_BINARY_PATTERN" "BYTE_TO_BINARY_PATTERN"\n", 
  //       BYTE_TO_BINARY(XinvVec[0]>>24), BYTE_TO_BINARY(XinvVec[0]>>16), BYTE_TO_BINARY(XinvVec[0]>>8), BYTE_TO_BINARY(XinvVec[0]>>0));


  createValidVec(REAL(pxvec), XvalidVec, DeltaX, N);
  // printf("\nFirst Int of Valid Vec: "BYTE_TO_BINARY_PATTERN" "BYTE_TO_BINARY_PATTERN" "BYTE_TO_BINARY_PATTERN" "BYTE_TO_BINARY_PATTERN"\n", 
  //       BYTE_TO_BINARY(XvalidVec[0]>>24), BYTE_TO_BINARY(XvalidVec[0]>>16), BYTE_TO_BINARY(XvalidVec[0]>>8), BYTE_TO_BINARY(XvalidVec[0]>>0));


  // double *rcimat2 = malloc(N * N * sizeof(double));
  permOnCuda(XinvVec, XvalidVec, singleY, singleOut, R, N, DeltaY, state);
  // printf("%f\n", out[0]);
  // rciBootWithCopy(REAL(prcimat), rcimat2, REAL(pout), R, N, xties, yties, state);
  convSingleToDouble(singleOut, REAL(pout), R);

  // free(rcimat2);
  free(XinvVec);
  free(XvalidVec);

  free(singleY);
  free(singleOut);

  UNPROTECT(1);
  
  return pout;
  
}



/*
nvcc -O3 -arch=sm_86 -G -I/usr/local/cuda/include -I/usr/share/R/include -L/usr/lib/R/lib -lR -L/usr/local/cuda/lib64 -lcurand --shared -Xcompiler -fPIC -o naivePermGPU.so naivePermGPU.cu

dyn.load("naivePermGPU.so")
system.time(t <- .Call("permCUDA", as.numeric(1:100), runif(100), 1e4, 100, 0, 0, runif(1)))

*/
