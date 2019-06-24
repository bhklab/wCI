#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>
#include <cmath>
#include <random>
#include <Rcpp.h>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <vector>
#include <string>

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
// [[Rcpp::export]]
bool usable (double x1, double x2, double delta) {
  if (fabs(x1 - x2) >= delta) {
    return true;
  } else {
    return false;
  }
}

/* function tests whether a pair is usable for concordance index calculation purposes. */
// [[Rcpp::export]]
bool usableHard (double x1, double x2, double delta) {
  if (fabs(x1 - x2) > delta) {
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


//function to do and, or
double logicOpF(bool x, bool y, std::string logicOp){

  if (logicOp.compare("and") == 0){
    return(x && y);
  }else if(logicOp.compare("or") == 0){
    return(x || y);
  }

  throw std::invalid_argument ("CI calculation failed.");

}


/* function calculates modified concordance index.
 Input: predictions x, observations y, cutoffs for x and y, deltas for x and y, confidence level alpha, flag outx, string alternative*/
// [[Rcpp::export]]
List concordanceIndex_modified_helper(std::vector<double> x, std::vector<double> y, double deltaX, double deltaY, double alpha, bool outx, std::string alternative, std::string logicOp) {

  int N = static_cast<int>(x.size());
  std::vector<int> c(N);
  std::vector<int> d(N);


  std::list<bool> cdseq;

  for (int i = 0; i < N; ++i) {
    c[i] = 0;
    d[i] = 0;
  }

  double numOfPairs = 0;

  for (int i = 0; i < N - 1; ++i) {
    for (int j = i + 1; j < N; ++j) {

      if (logicOpF(usable(x[i],x[j], deltaX), usable(y[i],y[j], deltaY), logicOp)) {
        if(usableHard(y[i],y[j], deltaY) && logicOp.compare("and") == 0){
          ++numOfPairs;

          if (outx == false && (!usableHard(x[i],x[j], deltaX))) {
            d[i] = d[i] + 1;
            c[i] = c[i] + 1;
            cdseq.push_back(true);
            cdseq.push_back(false);
          } else {

            if ((x[i] > x[j] & y[i] > y[j]) || (x[i] < x[j] & y[i] < y[j])) {
              c[i] = c[i] + 1;
              c[j] = c[j] + 1;
              cdseq.push_back(true);
              cdseq.push_back(true);
            } else {
              if(outx == true && (!usableHard(x[i],x[j], deltaX))){
                --numOfPairs;
              }else{

                d[i] = d[i] + 1;
                d[j] = d[j] + 1;
                cdseq.push_back(false);
                cdseq.push_back(false);
              }
            }
          }
        }else if (logicOp.compare("or") == 0){

          if(logicOpF(usable(x[i],x[j], deltaX), usableHard(y[i],y[j], deltaY), "and")){
            ++numOfPairs;

            if (outx == false && (!usableHard(x[i],x[j], deltaX))) {
              d[i] = d[i] + 1;
              c[i] = c[i] + 1;
              cdseq.push_back(true);
              cdseq.push_back(false);
            } else {

              if ((x[i] > x[j] & y[i] > y[j]) || (x[i] < x[j] & y[i] < y[j])) {
                c[i] = c[i] + 1;
                c[j] = c[j] + 1;
                cdseq.push_back(true);
                cdseq.push_back(true);
              } else {
                if(outx == true && (!usableHard(x[i],x[j], deltaX))){
                  --numOfPairs;
                }else{

                  d[i] = d[i] + 1;
                  d[j] = d[j] + 1;
                  cdseq.push_back(false);
                  cdseq.push_back(false);
                }
              }
            }
          }else if(usable(x[i],x[j], deltaX)){

            if(y[i]!=y[j]){
              ++numOfPairs;
              if (outx == false && (!usableHard(x[i],x[j], deltaX))) {
                d[i] = d[i] + 1;
                c[i] = c[i] + 1;
                cdseq.push_back(true);
                cdseq.push_back(false);

              } else {

                if ((x[i] > x[j] & y[i] > y[j]) || (x[i] < x[j] & y[i] < y[j])) {
                  c[i] = c[i] + 1;
                  c[j] = c[j] + 1;
                  cdseq.push_back(true);
                  cdseq.push_back(true);
                } else {
                  if(outx == true && (!usableHard(x[i],x[j], deltaX))){
                    --numOfPairs;
                  }else{
                    d[i] = d[i] + 1;
                    d[j] = d[j] + 1;
                    cdseq.push_back(false);
                    cdseq.push_back(false);
                  }
                }
              }
            }
          }else if(usableHard(y[i],y[j], deltaY)){
            ++numOfPairs;

            if (outx == false && (x[i]==x[j])) {
              d[i] = d[i] + 1;
              c[i] = c[i] + 1;
              cdseq.push_back(true);
              cdseq.push_back(false);
            } else {

              if ((x[i] > x[j] & y[i] > y[j]) || (x[i] < x[j] & y[i] < y[j])) {
                c[i] = c[i] + 1;
                c[j] = c[j] + 1;
                cdseq.push_back(true);
                cdseq.push_back(true);
              } else {
                if(outx == true && (x[i]==x[j])){
                  --numOfPairs;
                }else{

                  d[i] = d[i] + 1;
                  d[j] = d[j] + 1;
                  cdseq.push_back(false);
                  cdseq.push_back(false);
                }
              }
            }
          }


        }
      }
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

    List ret;
    ret["C"] = C;
    ret["D"] = D;
    ret["CC"] = CC;
    ret["DD"] = DD;
    ret["CD"] = CD;
    ret["N"] = N;
    ret["cdseq"] = cdseq;
    return ret;


}

