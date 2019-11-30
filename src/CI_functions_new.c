#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>

#include "xoroshiro128+.h"



struct retRes{
    double C;
    double D;
    double CC;
    double CD;
    double DD;
    uint64_t N;
};

typedef struct retRes retRes;



void computePair(double x1, double x2, double y1, double y2, double deltaX, double deltaY, int xties, int yties, 
  double *C1, double *C2, double *D1, double *D2, int logic){


  if(xties == yties & yties == 0){

    if(logic){
      if((fabs(x1 - x2) > deltaX) & (fabs(y1 - y2) > deltaY)){
        if((x1 - x2 > 0  & y1 - y2 > 0) | (x1 - x2 < 0  & y1 - y2 < 0)){
          *C1 += 1;
          *C2 += 1;
        } else {
          *D1 += 1;
          *D2 += 1;
        }
      }
    } else {
      if((fabs(x1 - x2) > deltaX) | (fabs(y1 - y2) > deltaY)){
        if((x1 - x2 > 0  & y1 - y2 > 0) | (x1 - x2 < 0  & y1 - y2 < 0)){
          *C1 += 1;
          *C2 += 1;
        } else {
          *D1 += 1;
          *D2 += 1;
        }
      }
    }


  } else if (xties == yties & yties == 1){

    if(logic){
      if((fabs(x1 - x2) > deltaX) & (fabs(y1 - y2) > deltaY)){
        if((x1 - x2 > 0  & y1 - y2 > 0) | (x1 - x2 < 0  & y1 - y2 < 0)){
          *C1 += 1;
          *C2 += 1;
        } else {
          *D1 += 1;
          *D2 += 1;
        }
      } else {
        *C1 += 0.5;
        *C2 += 0.5;
        *D1 += 0.5;
        *D2 += 0.5;
      }
    } else {
      if((fabs(x1 - x2) > deltaX) | (fabs(y1 - y2) > deltaY)){
        if((x1 - x2 > 0  & y1 - y2 > 0) | (x1 - x2 < 0  & y1 - y2 < 0)){
          *C1 += 1;
          *C2 += 1;
        } else {
          *D1 += 1;
          *D2 += 1;
        } 
      } else {
        *C1 += 0.5;
        *C2 += 0.5;
        *D1 += 0.5;
        *D2 += 0.5;
      }
    }


  } else if (xties == 1 & yties == 0){

    if(logic){
      if((fabs(x1 - x2) > deltaX) & (fabs(y1 - y2) > deltaY)){
        if((x1 - x2 > 0  & y1 - y2 > 0) | (x1 - x2 < 0  & y1 - y2 < 0)){
          *C1 += 1;
          *C2 += 1;
        } else {
          *D1 += 1;
          *D2 += 1;
        }
      } else if(fabs(y1 - y2) > deltaY){
        *C1 += 0.5;
        *C2 += 0.5;
        *D1 += 0.5;
        *D2 += 0.5;
      }
    } else {
      if((fabs(x1 - x2) > deltaX) | (fabs(y1 - y2) > deltaY)){
        if((x1 - x2 > 0  & y1 - y2 > 0) | (x1 - x2 < 0  & y1 - y2 < 0)){
          *C1 += 1;
          *C2 += 1;
        } else {
          *D1 += 1;
          *D2 += 1;
        } 
      } 
    }

  } else if (xties == 0 & yties == 1){

    if(logic){
      if((fabs(x1 - x2) > deltaX) & (fabs(y1 - y2) > deltaY)){
        if((x1 - x2 > 0  & y1 - y2 > 0) | (x1 - x2 < 0  & y1 - y2 < 0)){
          *C1 += 1;
          *C2 += 1;
        } else {
          *D1 += 1;
          *D2 += 1;
        }
      } else if((fabs(x1 - x2) > deltaX)){
        *C1 += 0.5;
        *C2 += 0.5;
        *D1 += 0.5;
        *D2 += 0.5;
      }
    } else {
      if((fabs(x1 - x2) > deltaX) | (fabs(y1 - y2) > deltaY)){
        if((x1 - x2 > 0  & y1 - y2 > 0) | (x1 - x2 < 0  & y1 - y2 < 0)){
          *C1 += 1;
          *C2 += 1;
        } else {
          *D1 += 1;
          *D2 += 1;
        } 
      } 
    }

  } else {
    // double error = 0/0; // Need to think about how to do this gracefully?
  }
  return;
}


retRes modified_concordance_index(double *x, double *y, uint64_t N, double deltaX, double deltaY, int xties, int yties, int logic) {

  double c[N];
  double d[N];

  // std::list<bool> cdseq;

  for (int i = 0; i < N; ++i) {
    c[i] = 0;
    d[i] = 0;
  }


  for (int i = 0; i < N - 1; ++i) {
    for (int j = i + 1; j < N; ++j) {

      computePair(x[i], x[j], y[i], y[j], deltaX, deltaY, xties, yties, &c[i], &c[j], &d[i], &d[j], logic);

    }
  }

  double C = 0.0;
  double D = 0.0;
  double CC = 0.0;
  double CD = 0.0;
  double DD = 0.0;



  for (int i = 0; i < N; ++i) {
    C += c[i];
    D += d[i];
    CC += c[i] * (c[i] - 1);
    DD += d[i] * (d[i] - 1);
    CD += c[i] * d[i];
  }


  retRes ret;
  ret.C = C;
  ret.D = D;
  ret.CC = CC;
  ret.DD = DD;
  ret.CD = CD;
  ret.N = N;
  // reteq.cdseq = cdseq;
  return ret;


}



SEXP newPairedConcIndex(SEXP pin_x,
             SEXP pin_y,
             SEXP pn,
             SEXP pxties,
             SEXP pyties,
             SEXP pdeltaX,
             SEXP pdeltaY,
             SEXP plogic){
  
  double Ndouble = *REAL(pn);
  
  // int discard_x_ties = *INTEGER(pdiscard_x_ties);
  // int discard_y_ties = *INTEGER(pdiscard_y_ties);

  uint64_t N = (uint64_t) Ndouble;
  
  
  SEXP pout = PROTECT(allocVector(REALSXP,6));
  
  double *out = REAL(pout);
  
  
  retRes res = modified_concordance_index(REAL(pin_x), REAL(pin_y), N, *REAL(pdeltaX), *REAL(pdeltaY), *INTEGER(pxties), *INTEGER(pyties), *INTEGER(plogic));
  

  out[0] = res.C;
  out[1] = res.D; 
  out[2] = res.CC;
  out[3] = res.DD;
  out[4] = res.CD;
  out[5] = res.N;

  // printf("%f\n", out[0]);
  
  UNPROTECT(1);
  
  return pout;
  
}
