#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>
#include <random>

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>



// [[Rcpp::plugins(openmp)]]

/* function to compute the concordance index with constraints on which pairs to compare
 * input:
 * x,y: vectors that will be compared. Has to have the same length and no missing values.
 * deltaX, deltaY: the delta range which pairs within are not compared
 * outx = True: do not include pairs that are equal
 * 
 * output:
 * c-index [-1,1]
 */

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
double concordanceIndex_modified(std::vector<double> x, std::vector<double> y, double deltaX, double deltaY, double alpha, int outx) {
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
        std::vector<double> pp(2);
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
        std::vector<double> oo(2);
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



/* function to shuffle an array.
 * input:
 * array:  a vector of at least length 2
 * 
 * output:
 * shuffled array
 */
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
std::vector<double> shuffle(std::vector<double> array) {

  int n =  array.size();
  std::vector<double> array1(n);
  for(int i=0; i< n;i++){
    array1[i] = array[i];
  }
  srand (time(NULL));
  std::random_device rd;
  std::mt19937_64 g(rd());
  std::shuffle(array1.begin(),array1.end(),g);
  return(array1);
}



/* function to compute the pvalue based on the modified concordance index
 * input:
 * x,y: vectors that will be compared. Has to have the same length and no missing values.
 * deltaX, deltaY: the delta range which pairs within are not compared
 * outx = True: do not include pairs that are equal
 * permutations: number of permutations to do on the shuffled array
 * nThreads: number of threads to use.
 * 
 * output:
 * an array of length(permutation) with cindecies of the shuffled vectors
 */

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
std::vector<double> permute_concordanceIndex_modified(std::vector<double> x, std::vector<double> y, double deltaX, double deltaY, double alpha, int outx, int permutations,int nThreads) {
  
  std::vector<double> randomPermut(permutations);
  
  omp_set_dynamic(0);     // Explicitly disable dynamic teams
  omp_set_num_threads(nThreads); // Use n threads for all consecutive parallel regions
  
  #pragma omp parallel for
  for(int i=0; i < permutations;i++){
    
    std::vector<double> xShuffled(x.size());
    xShuffled = shuffle(x);
    randomPermut[i] = concordanceIndex_modified(xShuffled,y,deltaX,deltaY,alpha,outx);
    
  }
  
  
  return(randomPermut);
  
  
}
