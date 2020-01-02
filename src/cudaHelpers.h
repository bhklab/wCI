#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <curand.h>



static const char *curandGetErrorString(curandStatus_t error);

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true){
	   if (code != cudaSuccess) 
		      {
			            fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
				          if (abort) exit(code);
					     }
}

#define gpuErrchkRand(ans) { gpuAssertRand((ans), __FILE__, __LINE__); }
inline void gpuAssertRand(curandStatus_t code, const char *file, int line, bool abort=true){
	if (code != CURAND_STATUS_SUCCESS)
		{
		        fprintf(stderr,"GPUassert: %s %s %d\n", curandGetErrorString(code), file, line);
			if (abort) exit(code);
	       	}
}