#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>
#include <cmath>
#include <omp.h>
#include <random>
#include <Rcpp.h>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <vector>

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

/* function tests whether a pair is usable for concordance index calculation purposes. */
bool usable (std::vector<double> x, double cutoff, double delta) {
  if ((x[0] >= cutoff || x[1] >= cutoff) && std::abs(x[0] - x[1]) >= delta) {
    return true;
  } else {
    return false;
  }
}

/* function returns sign of argument */
int sgn(double x) {
  if (x > 0) {
    return(1);
  } else if (x == 0) {
    return(0);
  } else if (x < 0) {
    return(-1);
  } else {
    throw std::invalid_argument ("erfinv received an argument that was not a real number.");
  }
}

/* function computes inverse of error function for normal quantile calculation */
double erfinv(double x){
  int sgnx = sgn(x);
  x = 1 - x * x;
  double lnx = log(x);
  double a = 4.3307467508 + lnx / 2; // magic number for Winitzki approximation
  double b = 6.80272108844 * lnx; // other magic number for approximation
  return(sgnx * sqrt(sqrt(a * a - b) - a));
}

/* function calculates modified concordance index.
   Input: predictions x, observations y, cutoffs for x and y, deltas for x and y, confidence level alpha, flag outx, string alternative*/
// [[Rcpp::export]]
double concordanceIndex_modified(std::vector<double> x, std::vector<double> y, std::vector<double> cutoff, std::vector<double> delta, double alpha, bool outx, std::string alternative) {
  
  int N = static_cast<int>(x.size());
  std::vector<int> c(N);
  std::vector<int> d(N);
  
  for (int i = 0; i < N; ++i) {
    c[i] = 0;
    d[i] = 0;
  }
  
  for (int i = 0; i < N - 1; ++i) {
    for (int j = i + 1; j < N; ++j) {
      if (usable(x, cutoff[0], delta[0]) && usable(y, cutoff[1], delta[1])) {
        if (outx == false && (x[i] == x[j] || y[i] == y[j])) { // should this be an xor?
          ++d[i];
          ++d[j];
        } else {
          if ((x[0] > x[1] && y[0] > y[1]) || (x[0] < x[1] && y[0] < y[1])) {
            ++c[i];
            ++c[j];
          } else {
            ++d[i];
            ++d[j];
          }
        }
      }
    }
  }
  
  int C = 0;
  int D = 0;
  int CC = 0;
  int CD = 0;
  int DD = 0;
  
  for (int i = 0; i < N; ++i) {
    C += c[i];
    D += d[i];
    CC += c[i] * (c[i] - 1);
    DD += d[i] * (d[i] - 1);
    CD += c[i] * d[i];
  }
  
  if (C == 0 && D == 0) {
    throw std::invalid_argument ("All pairs were thrown out. Consider changing cutoff and/or delta.");
  }
  
  double cindex = static_cast<double>(C) / (C + D); //static_cast prevents integer division from occurring
  double varp = 4 * N * (N - 1) / (N - 2) * (pow(D, 2.0) * CC - 2 * C * D * CD + pow(C, 2.0) * DD) / pow(C + D, 4.0);
  if (varp / N >= 0) {
    double sterr = sqrt(varp / N);
    double ci = sterr * 1.41421356237 * erfinv(2 * alpha); // magic number is sqrt(2)
    double p = (1 + erf((cindex - 0.5) / sterr / 1.41421356237)) / 2;
    if (std::string::compare(alternative, "less") == 0) {
    } else if (std::string::compare(alternative, "greater") == 0) {
      p <- 1 - p;
    } else if (std::string::compare(alternative, "two.sided") == 0) {
      p <- 2 * std::min(p, 1 - p);
    } else {
      throw std::invalid_argument ("'alternative' must be set to one of 'less', 'greater', or 'two.sided'.");
    }
    std::cout << "Calculated Concordance Index: " << cindex << "\n";
    std::cout << (1 - alpha) * 100 << "% confidence interval: [" << std::max(cindex - ci, 0.0) << ", " << std::min(cindex + ci, 1.0) << "]\n";
    std::cout << "p-value: " << p << "\n";
    return cindex;
  } else {
    throw std::invalid_argument ("CI calculation failed.");
  }
}

/* function to shuffle an array.
 * input:
 * array:  a vector of at least length 2
 * 
 * output:
 * shuffled array
 */

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

// [[Rcpp::export]]
std::vector<double> permute_concordanceIndex_modified(std::vector<double> x, std::vector<double> y, double deltaX, double deltaY, double alpha, int outx, int permutations,int nThreads) {
  
  std::vector<double> randomPermut(permutations);
  
  omp_set_dynamic(0);     // Explicitly disable dynamic teams
  omp_set_num_threads(nThreads); // Use n threads for all consecutive parallel regions
  
  #pragma omp parallel
  {
    std::vector<double> xShuffled(x.size());
    #pragma omp for
  for(int i=0; i < permutations;i++){
    
    
    xShuffled = shuffle(x);
    randomPermut[i] = concordanceIndex_modified(xShuffled,y,deltaX,deltaY,alpha,outx);
    
  }
  
  }
  return(randomPermut);
  
  
}
