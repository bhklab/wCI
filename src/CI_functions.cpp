#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>
#include <random>

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
//#include <omp.h>



// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
double concordanceIndex_modified(NumericVector x, NumericVector y, double deltaX, double deltaY, double alpha, int outx) {
  double N = (double)x.size();
  double ch = 0;
  double dh = 0;
  double uh = 0;
  int pos = 0;
  for(int i=0; i < (x.size()-1); i++){
    for(int j=i+1; j < x.size(); j++){
      pos++;
      int firstVectorPair = fabs(x[i] - x[j]) > deltaX;
      int secondVectorPair = fabs(y[i] - y[j]) > deltaY;
      if(firstVectorPair | secondVectorPair){
        NumericVector pp(2);
        if(x[i] < x[j]){
          pp[0] = 1;
          pp[1] = 2;

        } else if(x[i] > x[j]){
          pp[0] = 2;
          pp[1] = 1;

        } else{
          pp[0] = 1;
          pp[1] = 1;
        }
        NumericVector oo(2);
        if(y[i] < y[j]){
          oo[0] = 1;
          oo[1] = 2;
        } else if(y[i] > y[j]){
          oo[0] = 2;
          oo[1] = 1;
        } else{
          oo[0] = 1;
          oo[1] = 1;
        }
        if((pp[0]==pp[1])||(oo[0]==oo[1])){
          if(outx){
            uh = uh + 1;
          }else{
            dh = dh + 1;
          }
        }else{
          if((pp[0]==oo[0]) && (pp[1]==oo[1])){
            ch = ch + 1;
          }else{
            dh = dh + 1;
          }
        }
        //		free(pp);
        //		free(oo);
      }
    }
  }
  double pc = (1/N * (N-1)) * ch;
  double pd = (1/N * (N-1)) * dh;
  double cindex = pc / (pc + pd);
  double dci = (cindex - 0.5)*2.0;
  ////return(cindex);
  return(dci);
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
NumericVector shuffle(NumericVector array) {

  int n =  array.size();
  NumericVector array1(n);

  for(int i=0; i< n;i++){
    array1[i] = array[i];
  }

  srand((unsigned)time(NULL));
  for (int i = 0; i < n - 1; i++) {
    size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
    double t = array1[j];
    array1[j] = array1[i];
    array1[i] = t;
  }

  return(array1);
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
NumericVector fisher_yates(NumericVector array) {

  int n =  array.size();
  NumericVector array1(n);

  for(int i=0; i< n;i++){
    array1[i] = array[i];
  }
  srand (time(NULL));

  for(int i=n-1; i>0; i--) {
    int idx = (rand() % (i+1));
    std::swap(array1[idx], array1[i]);
  }

  return(array1);
}



// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
NumericVector shuffle2(NumericVector array) {

  int n =  array.size();
  NumericVector array1(n);
  for(int i=0; i< n;i++){
    array1[i] = array[i];
  }
  srand (time(NULL));
  std::random_device rd;
  std::mt19937_64 g(rd());
  std::shuffle(array1.begin(),array1.end(),g);
  return(array1);
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
NumericVector permute_concordanceIndex_modified(NumericVector x, NumericVector y, double deltaX, double deltaY, double alpha, int outx, int permutations) {

  NumericVector randomPermut(permutations);

  NumericVector xShuffled(x.size());

  //	#pragma omp parallel for private(xShuffled)
  for(int i=0; i < permutations;i++){
    xShuffled = shuffle2(x);
    //	std::shuffle
    //xShuffled = shuffle(x);
    //NumericVector yShuffled = shuffle(y,sizeY);
    randomPermut[i] = concordanceIndex_modified(xShuffled,y,deltaX,deltaY,alpha,outx);
  }
  return(randomPermut);
}

