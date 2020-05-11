/*
// Fast permutations for rCI using a naive matrix based approach.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>

#include "xoroshiro128+.h"


struct runningStat2{
    double numerator;
    double denominator;
};

typedef struct runningStat2 runningStat2;


// struct returnRes2{
//     double pval;
//     int earlyStop;
// };

// typedef struct returnRes2 returnRes2;

// void printVec(uint64_t *list, uint64_t N){
//   for(uint64_t i = 0; i < N; i ++){
//     printf("%lld, ", list[i]);
//   }
// }

// Using "inside out" Fisher Yates https://en.wikipedia.org/wiki/Fisher–Yates_shuffle
void sampleIdxBoot(uint64_t N, uint64_t *permPointer, uint64_t *state){
  // uint64_t j;
  // permPointer[0] = 0;

  for(uint64_t i = 0; i < N; i++){
      // j = generate_random_index(state, i);
      // if(j != i){
      //   permPointer[i] = permPointer[j];
      // }
      permPointer[i] = generate_random_index(state, N - 1);
  }

}
  
  


// runningStat2 returnCurrentIncrement2(uint64_t xVal, uint64_t yVal, int xties, int yties){

//   runningStat2 res;

//   if(xties == yties & yties == 0){
//     // printf("Ignore ties\n. ");
//       res.numerator = (double) ((xVal == yVal) & (xVal != 0) & (yVal != 0)); // Datasets agree and are non-zero
//       res.denominator = (double) ((xVal != 0) & (yVal != 0)); // Both are non-zero
//     } else if (xties == yties & yties == 1){

//       res.numerator = ((double) xVal == yVal)*0.5 + (double) ((xVal != 0) & (yVal != 0)) * 0.5 ; // Datasets agree (0.5) non zero (0.5)
//       res.denominator = 1; // All pairs are counted. 

//     } else {
//       if (xties == 1 & yties == 0){
//         // printf("OUTX=FALSE\n");
//         res.numerator = (((double) ((xVal == yVal) & (yVal != 0))) + (double) ((yVal != 0) & (xVal == 0)) * 0.5); // Believe it or not, this is actually faster than using short circuting if statements!!
//         res.denominator = (double) (yVal != 0);
//       }
//       if (xties == 0 & yties == 1){
//         res.numerator = (((double) ((xVal == yVal) & (yVal != 0))) + (double) ((yVal != 0) & (xVal == 0)) * 0.5);
//         res.denominator = (double) (xVal != 0);
//       }
//     }
//     return res;
// }


uint64_t returnLinearIdx2(uint64_t row, uint64_t column, uint64_t N){
  return(column*N + row);
}



void rciBoot(double *rcimat, double *outVec, uint64_t R, uint64_t N, int xties, int yties, uint64_t *state){

  double currCI;
  uint64_t i = 0;
  // uint64_t j = 0;
  uint64_t *permIdx = malloc(N * sizeof(uint64_t));


  double curVal;

  runningStat2 RS;
  RS.numerator = 0;
  RS.denominator = 0;
  // runningStat2 CS;

  while(i < R){

    sampleIdxBoot(N, permIdx, state);
    // printVec(permIdx, N);
    // printf("\n");

    for(uint64_t j = 0; j < N; j++){

      for(uint64_t k = 0; k < N; k++){

        curVal = rcimat[returnLinearIdx2(permIdx[k],permIdx[j],N)];
        // yVal = ymat[returnLinearIdx2(permIdx[k],permIdx[j],N)];
        // printf("Cur Val: %f.", curVal);

        // CS = returnCurrentIncrement2(xVal, yVal, xties, yties);
        // if(xVal == 0 & yVal != 0){
        //   printf("CurNum: %f. \n CurDenom: %f.\n", CS.numerator, CS.denominator);
        // }
        // printf("CurNum: %f. \n CurDenom: %f.\n", CS.numerator, CS.denominator);

        RS.numerator += (curVal * (double)(curVal > 0));
        // RS.denominator += (double)(curVal != 0);
        RS.denominator += (double)(curVal != 0) * 2;
        // if(curVal == 0){
        //   printf("We have a 0");
        // }

      }

    }

    currCI = ((double)RS.numerator)/((double)RS.denominator);
    // printf("Obs Num: %f. \n Obs Denom: %f.\n", RS.numerator, RS.denominator);

    // printf("Current I: %d. \n Curr CI: %f.\n", i, currCI);
    outVec[i] = currCI;

    RS.numerator = 0;
    RS.denominator = 0;
    i++;
  }

  free(permIdx);

}



void rciBootWithCopy(double *rcimat, double *rcimat2, double *outVec, uint64_t R, uint64_t N, int xties, int yties, uint64_t *state){

  double currCI;
  uint64_t i = 0;
  // uint64_t j = 0;
  uint64_t *permIdx = malloc(N * sizeof(uint64_t));


  double curVal;

  runningStat2 RS;
  RS.numerator = 0;
  RS.denominator = 0;
  // runningStat2 CS;

  while(i < R){

    sampleIdxBoot(N, permIdx, state);
    // printVec(permIdx, N);
    // printf("\n");

    for(uint64_t j = 0; j < N; j++){

      for(uint64_t k = 0; k < N; k++){

        rcimat2[returnLinearIdx2(k,j,N)] = rcimat[returnLinearIdx2(permIdx[k],permIdx[j],N)];
        // yVal = ymat[returnLinearIdx2(permIdx[k],permIdx[j],N)];
        // printf("Cur Val: %f.", curVal);

        // CS = returnCurrentIncrement2(xVal, yVal, xties, yties);
        // if(xVal == 0 & yVal != 0){
        //   printf("CurNum: %f. \n CurDenom: %f.\n", CS.numerator, CS.denominator);
        // }
        // printf("CurNum: %f. \n CurDenom: %f.\n", CS.numerator, CS.denominator);

       
        // if(curVal == 0){
        //   printf("We have a 0");
        // }

      }

    }
    for(uint64_t j = 0; j < N; j++){
      for(uint64_t k = 0; k < N; k++){
        curVal = rcimat2[returnLinearIdx2(k,j,N)];
        RS.numerator += (curVal * (double)(curVal > 0));
        RS.denominator += (double)(curVal != 0);
      }
    }  
    currCI = ((double)RS.numerator)/((double)RS.denominator);
    // printf("Obs Num: %f. \n Obs Denom: %f.\n", RS.numerator, RS.denominator);

    // printf("Current I: %d. \n Curr CI: %f.\n", i, currCI);
    outVec[i] = currCI;

    RS.numerator = 0;
    RS.denominator = 0;
    i++;
  }

  free(permIdx);

}

SEXP bootC(SEXP prcimat,
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
  rciBoot(REAL(prcimat), REAL(pout), R, N, xties, yties, state);
  // printf("%f\n", out[0]);
  // rciBootWithCopy(REAL(prcimat), rcimat2, REAL(pout), R, N, xties, yties, state);


  // free(rcimat2);
  UNPROTECT(1);
  
  return pout;
  
}


