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


/* function tests whether a pair is usable for concordance index calculation purposes [Soft threshold]. */
// [[Rcpp::export]]
bool usable (double x1, double x2, double delta) {
  if (fabs(x1 - x2) >= delta) {
    return true;
  } else {
    return false;
  }
}

/* function tests whether a pair is usable for concordance index calculation purposes [Strict threshold]. */
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
 Input: predictions x, observations y, deltaValues for x, deltaValues for y, confidence level alpha, flag outx representing if ties in x are counted, string alternative["two.sides", "less","greater"], string logicOp ["and","or"]
 Output: stats needed to calculate CI, confidence interval, and p-value
 */
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

            if ((x[i] > x[j] && y[i] > y[j]) || (x[i] < x[j] && y[i] < y[j])) {
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

              if ((x[i] > x[j]&&y[i] > y[j]) || (x[i] < x[j] && y[i] < y[j])) {
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

                if ((x[i] > x[j] && y[i] > y[j]) || (x[i] < x[j] && y[i] < y[j])) {
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

              if ((x[i] > x[j] && y[i] > y[j]) || (x[i] < x[j] && y[i] < y[j])) {
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





/* functions to calculate kernels for weights of deltas instead of fixed values
 */
// [[Rcpp::export]]
double kernel_gaussian_C( double x, double m, double s){

  const double PI_cons = 3.141592653589793238463;

  return((1/sqrt(2*PI_cons*pow(s,2)))*exp(-pow((x-m),2)/(2*pow(s,2))));

}

/* functions to calculate kernels for weights of deltas instead of fixed values
 */
// [[Rcpp::export]]
double kernel_laplace_C( double x, double m, double b){

  return((1/(2*b)) * exp(-fabs(x - m) / b));

}

/* function calculates modified concordance index.
Input: predictions x, observations y, deltaValues for x, deltaValues for y, kernelFunction fo x deltas, kernelFunction fo y deltas, confidence level alpha, flag outx representing if ties in x are counted, string alternative["two.sides", "less","greater"], string logicOp ["and","or"]
Output: stats needed to calculate CI, confidence interval, and p-value
*/
// [[Rcpp::export]]
List concordanceIndex_modified_helper_weighted(std::vector<double> x, std::vector<double> y, std::vector<double> deltaX, std::vector<double> deltaY,std::string weightingFun_pred,std::string weightingFun_obs,
                                               double alpha, bool outx, std::string alternative, std::string logicOp,double max_weight, double max_weight_obs) {

  int N = static_cast<int>(x.size());
  std::vector<double> c(N);
  std::vector<double> d(N);

  NumericVector w_order(N);

  for (int i = 0; i < N; ++i) {
    w_order[i] = i;
  }

  std::list<bool> cdseq;

  for (int i = 0; i < N; ++i) {
    c[i] = 0;
    d[i] = 0;
  }

  double numOfPairs = 0;
  double w = 0;
  double obs_w = 0;
  double pred_w = 0;

  for (int i = 0; i < N - 1; ++i) {
    for (int j = i + 1; j < N; ++j) {

      if((weightingFun_obs.compare("kernel_gaussian") == 0 && weightingFun_pred.compare("kernel_gaussian") == 0) | (weightingFun_obs.compare("kernel_laplace") == 0 && weightingFun_pred.compare("kernel_laplace") == 0)){

        if(weightingFun_obs.compare("kernel_gaussian") == 0){
          obs_w = fabs(log10(kernel_gaussian_C(y[w_order[i]] - y[w_order[j]],0.0002001131,0.0939948369)));

          pred_w = fabs(log10(kernel_gaussian_C(x[w_order[i]] - x[w_order[j]],0.0002001131,0.0939948369)));
          w = 1/max_weight * std::max(obs_w, pred_w);
        }else if(weightingFun_obs.compare("kernel_laplace") == 0){
          obs_w = fabs(log10(kernel_laplace_C(y[w_order[i]] - y[w_order[j]],-0.001785626,0.061982848)));
          pred_w = fabs(log10(kernel_laplace_C(x[w_order[i]] - x[w_order[j]],-0.001785626,0.061982848)));
          w = 1/max_weight * std::max(obs_w, pred_w);
        }

      }else if((weightingFun_obs.compare("kernel_gaussian") == 0) | (weightingFun_obs.compare("kernel_laplace") == 0)){
        if(weightingFun_obs.compare("kernel_gaussian") == 0){
          obs_w = fabs(log10(kernel_gaussian_C(y[w_order[i]] - y[w_order[j]],0.0002001131,0.0939948369)));
        }else if(weightingFun_obs.compare("kernel_laplace") == 0){
          obs_w = fabs(log10(kernel_laplace_C(y[w_order[i]] - y[w_order[j]],-0.001785626,0.061982848)));
        }

        w = 1/max_weight * obs_w;
      }else{
        w = 1;
      }


      NumericVector dX(2);
      NumericVector dY(2);

      dX[0] = deltaX[i];
      dX[1] = deltaX[j];


      dY[0] = deltaY[i];
      dY[1] = deltaY[j];

      double deltaXF = Rcpp::sample(dX,2,false,R_NilValue)[1];
      double deltaYF = Rcpp::sample(dY,2,false,R_NilValue)[1];


      if (logicOpF(usable(x[i],x[j], deltaXF), usable(y[i],y[j], deltaYF), logicOp)) {
        if(usableHard(y[i],y[j], deltaYF) && logicOp.compare("and") == 0){
          ++numOfPairs;

          if (outx == false && (!usableHard(x[i],x[j], deltaXF))) {
            d[i] = d[i] + w;
            c[i] = c[i] + w;
            cdseq.push_back(true);
            cdseq.push_back(false);
          } else {

            if ((x[i] > x[j] && y[i] > y[j]) || (x[i] < x[j] && y[i] < y[j])) {
              c[i] = c[i] + w;
              c[j] = c[j] + w;
              cdseq.push_back(true);
              cdseq.push_back(true);
            } else {
              if(outx == true && (!usableHard(x[i],x[j], deltaXF))){
                --numOfPairs;
              }else{

                d[i] = d[i] + w;
                d[j] = d[j] + w;
                cdseq.push_back(false);
                cdseq.push_back(false);
              }
            }
          }
        }else if (logicOp.compare("or") == 0){

          if(logicOpF(usable(x[i],x[j], deltaXF), usableHard(y[i],y[j], deltaYF), "and")){
            ++numOfPairs;

            if (outx == false && (!usableHard(x[i],x[j], deltaXF))) {
              d[i] = d[i] + w;
              c[i] = c[i] + w;
              cdseq.push_back(true);
              cdseq.push_back(false);
            } else {

              if ((x[i] > x[j] && y[i] > y[j]) || (x[i] < x[j] && y[i] < y[j])) {
                c[i] = c[i] + w;
                c[j] = c[j] + w;
                cdseq.push_back(true);
                cdseq.push_back(true);
              } else {
                if(outx == true && (!usableHard(x[i],x[j], deltaXF))){
                  --numOfPairs;
                }else{

                  d[i] = d[i] + w;
                  d[j] = d[j] + w;
                  cdseq.push_back(false);
                  cdseq.push_back(false);
                }
              }
            }
          }else if(usable(x[i],x[j], deltaXF)){

            if(y[i]!=y[j]){
              ++numOfPairs;
              if (outx == false && (!usableHard(x[i],x[j], deltaXF))) {
                d[i] = d[i] + w;
                c[i] = c[i] + w;
                cdseq.push_back(true);
                cdseq.push_back(false);

              } else {

                if ((x[i] > x[j] && y[i] > y[j]) || (x[i] < x[j] && y[i] < y[j])) {
                  c[i] = c[i] + w;
                  c[j] = c[j] + w;
                  cdseq.push_back(true);
                  cdseq.push_back(true);
                } else {
                  if(outx == true && (!usableHard(x[i],x[j], deltaXF))){
                    --numOfPairs;
                  }else{
                    d[i] = d[i] + w;
                    d[j] = d[j] + w;
                    cdseq.push_back(false);
                    cdseq.push_back(false);
                  }
                }
              }
            }
          }else if(usableHard(y[i],y[j], deltaYF)){
            ++numOfPairs;

            if (outx == false && (x[i]==x[j])) {
              d[i] = d[i] + w;
              c[i] = c[i] + w;
              cdseq.push_back(true);
              cdseq.push_back(false);
            } else {

              if ((x[i] > x[j] && y[i] > y[j]) || (x[i] < x[j] && y[i] < y[j])) {
                c[i] = c[i] + w;
                c[j] = c[j] + w;
                cdseq.push_back(true);
                cdseq.push_back(true);
              } else {
                if(outx == true && (x[i]==x[j])){
                  --numOfPairs;
                }else{

                  d[i] = d[i] + w;
                  d[j] = d[j] + w;
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


