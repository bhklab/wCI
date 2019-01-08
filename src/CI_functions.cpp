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
bool usable (double x1, double x2, double delta) {
  if (std::abs(x1 - x2) >= delta) {
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
        if(y[i]!=y[j]){
          ++numOfPairs;
          if (outx == false && (x[i] == x[j])) {
            d[i] = d[i] + 1;
            //d[j] = d[j] + 0.25;
            c[i] = c[i] + 1;
            //c[j] = c[j] + 0.25;
            cdseq.push_back(true);
            cdseq.push_back(false);
          } else {


            if ((x[i] > x[j] & y[i] > y[j]) || (x[i] < x[j] & y[i] < y[j])) {
              ++c[i];
              ++c[j];
              cdseq.push_back(true);
              cdseq.push_back(true);
            } else {
              if(outx == true && (x[i] == x[j])){
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





/* functions to calculate kernels for weights
 */
// [[Rcpp::export]]
double kernel_gaussian_C( double x, double m, double s){

  const double PI_cons = 3.141592653589793238463;

  return((1/sqrt(2*PI_cons*pow(s,2)))*exp(-pow((x-m),2)/(2*pow(s,2))));

}

// [[Rcpp::export]]
double kernel_laplace_C( double x, double m, double b){

  return((1/(2*b)) * exp(-fabs(x - m) / b));

}
/* function calculates modified concordance index.
 Input: predictions x, observations y, cutoffs for x and y, deltas for x and y, confidence level alpha, flag outx, string alternative*/
// [[Rcpp::export]]
List concordanceIndex_modified_helper_weighted(std::vector<double> x, std::vector<double> y, double deltaX, double deltaY,std::string weightingFun_pred,std::string weightingFun_obs,
                                               double alpha, bool outx, std::string alternative, std::string logicOp,double max_weight, double max_weight_obs, bool permute_weights) {

  int N = static_cast<int>(x.size());
  std::vector<double> c(N);
  std::vector<double> d(N);

  NumericVector w_order(N);

  for (int i = 0; i < N; ++i) {
    w_order[i] = i;
  }
  if (permute_weights == true){
    w_order = Rcpp::sample(w_order,N,false,R_NilValue);
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

      if((weightingFun_obs.compare("kernel_gaussian") == 0 & weightingFun_pred.compare("kernel_gaussian") == 0) | (weightingFun_obs.compare("kernel_laplace") == 0 & weightingFun_pred.compare("kernel_laplace") == 0)){

        if(weightingFun_obs.compare("kernel_gaussian") == 0){
          //w = fabs(log10(kernel_gaussian_C(y[i] - y[j],0.0002037366,0.0919937995))) * fabs(log10(kernel_gaussian_C(x[i] - x[j],0.0002037366,0.0919937995)));
          obs_w = fabs(log10(kernel_gaussian_C(y[w_order[i]] - y[w_order[j]],0.0002001131,0.0939948369)));
          //if(obs_w < 0){
          //  obs_w = 0;
          //}
          pred_w = fabs(log10(kernel_gaussian_C(x[w_order[i]] - x[w_order[j]],0.0002001131,0.0939948369)));
          //if(pred_w < 0){
          //  pred_w = 0;
          //}
          w = 1/max_weight * std::max(obs_w, pred_w);
        }else if(weightingFun_obs.compare("kernel_laplace") == 0){
          obs_w = fabs(log10(kernel_laplace_C(y[w_order[i]] - y[w_order[j]],-0.001785626,0.061982848)));
          //if(obs_w < 0){
          //  obs_w = 0;
          //}
          pred_w = fabs(log10(kernel_laplace_C(x[w_order[i]] - x[w_order[j]],-0.001785626,0.061982848)));
          //if(pred_w < 0){
          //  pred_w = 0;
          //}
          w = 1/max_weight * std::max(obs_w, pred_w);
        }

        // w <- abs(log10(weightingFun_obs(observations[i] - observations[j]))) * abs(log10(weightingFun_obs(predictions[i] - predictions[j])))
       // w = 1;
      }else if((weightingFun_obs.compare("kernel_gaussian") == 0) | (weightingFun_obs.compare("kernel_laplace") == 0)){
        if(weightingFun_obs.compare("kernel_gaussian") == 0){
          obs_w = 1/max_weight_obs * fabs(log10(kernel_gaussian_C(y[w_order[i]] - y[w_order[j]],0.0002001131,0.0939948369)));
          //if(obs_w < 0){
          //  obs_w = 0;
          //}
        }else if(weightingFun_obs.compare("kernel_laplace") == 0){
          obs_w = 1/max_weight_obs * fabs(log10(kernel_laplace_C(y[w_order[i]] - y[w_order[j]],-0.001785626,0.061982848)));
          //if(obs_w < 0){
          //  obs_w = 0;
          //}
        }
        // w <- abs(log10(weightingFun_obs(observations[i] - observations[j])))
        //w = 1;
      }else{
        w = 1;
      }


      if (logicOpF(usable(x[i],x[j], deltaX), usable(y[i],y[j], deltaY), logicOp)) {
        if(y[i]!=y[j]){
          ++numOfPairs;
          if (outx == false && (x[i] == x[j])) {
            d[i] = d[i] + w;
            c[i] = c[i] + w;
            cdseq.push_back(true);
            cdseq.push_back(false);
          } else {


            if ((x[i] > x[j] & y[i] > y[j]) || (x[i] < x[j] & y[i] < y[j])) {
              c[i] = c[i] + w;
              c[j] = c[j] + w;
              cdseq.push_back(true);
              cdseq.push_back(true);
            } else {
              if(outx == true && (x[i] == x[j])){
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



/* function calculates modified concordance index.
 Input: predictions x, observations y, cutoffs for x and y, deltas for x and y, confidence level alpha, flag outx, string alternative*/
// [[Rcpp::export]]
List concordanceIndex_modified_helper_parallel(std::vector<double> x, std::vector<double> y, double deltaX, double deltaY, double alpha, bool outx, std::string alternative, std::string logicOp) {

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
        if(y[i]!=y[j]){
          ++numOfPairs;
          if (outx == false && (x[i] == x[j])) { // should this be an xor?
            ++d[i];
            ++d[j];
            cdseq.push_back(false);
            cdseq.push_back(false);
          } else {


            if ((x[i] > x[j] & y[i] > y[j]) || (x[i] < x[j] & y[i] < y[j])) {
              ++c[i];
              ++c[j];
              cdseq.push_back(true);
              cdseq.push_back(true);
            } else {
              if(outx == true && (x[i] == x[j])){
                --numOfPairs;
              }else{
                ++d[i];
                ++d[j];
                cdseq.push_back(false);
                cdseq.push_back(false);
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



/* function calculates modified concordance index.
   Input: predictions x, observations y, cutoffs for x and y, deltas for x and y, confidence level alpha, flag outx, string alternative*/
// [[Rcpp::export]]
List concordanceIndex_modified_AllinC(std::vector<double> x, std::vector<double> y, double deltaX, double deltaY, double alpha, bool outx, std::string alternative, std::string logicOp) {

  int N = static_cast<int>(x.size());
  std::vector<int> c(N);
  std::vector<int> d(N);

  for (int i = 0; i < N; ++i) {
    c[i] = 0;
    d[i] = 0;
  }

  logicOp = "|";

  double numOfPairs = 0;

  for (int i = 0; i < N - 1; ++i) {
    for (int j = i + 1; j < N; ++j) {
      if (usable(x[i],x[j], deltaX) && usable(y[i],y[j], deltaY)) {
        if(y[i]!=y[j]){
          ++numOfPairs;
        if (outx == false && (x[i] == x[j])) { // should this be an xor?
          ++d[i];
          ++d[j];
        } else {


          if ((x[i] > x[j] & y[i] > y[j]) || (x[i] < x[j] & y[i] < y[j])) {
            ++c[i];
            ++c[j];
          } else {
            if(outx == true && (x[i] == x[j])){
              --numOfPairs;
            }else{
              ++d[i];
              ++d[j];
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


  std::cout << "C: " << C << ", D:" << D << "\n";
  if (C == 0 && D == 0) {
    throw std::invalid_argument ("All pairs were thrown out. Consider changing cutoff and/or delta.");
  }
  double cindex = static_cast<double>(C) / (C + D); //static_cast prevents integer division from occurring

  double varp = 4 * (( (pow(D, 2.0) * CC ) - (2 * C * D * CD) + (pow(C, 2.0) * DD)) / pow(C + D, 4.0)) * N * (N - 1) / (N - 2);

  if (varp / N >= 0) {
    double sterr = sqrt(varp / N);



    double ci = sterr * 1.41421356237 * erfinv(2 * alpha); // magic number is sqrt(2)
    double p = (1 + erf((cindex - 0.5) / sterr / 1.41421356237)) / 2;
    if (alternative.compare("less") == 0) {
    } else if (alternative.compare("greater") == 0) {
      p = 1 - p;
    } else if (alternative.compare("two.sided") == 0) {
      p = 2 * std::min(p, 1 - p);
    } else {
      throw std::invalid_argument ("'alternative' must be set to one of 'less', 'greater', or 'two.sided'.");
    }
    std::cout << "Calculated Concordance Index: " << cindex << "\n";
    std::cout << (1 - alpha) * 100 << "% confidence interval: [" << std::max(cindex - ci, 0.0) << ", " << std::min(cindex + ci, 1.0) << "]\n";
    std::cout << "p-value: " << p << "\n";

    std::cout << "ci: " << ci << "\n";
    std::cout << "sterr: " << sterr << "\n";

    List ret;
    ret["cindex"] = cindex;
    ret["p.value"] = p;
    ret["lower"] = std::max(cindex - ci, 0.0);
    ret["upper"] = std::min(cindex + ci, 1.0);
    ret["relevant.pairs.no"] = numOfPairs;
    ret["p.concordant.pairs"] = C;
    return ret;

//    return cindex;
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


// [[Rcpp::export]]
std::vector<int> shuffleInt(std::vector<int> array) {

  int n =  array.size();
  std::vector<int> array1(n);
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

/*  omp_set_dynamic(0);     // Explicitly disable dynamic teams
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
  */
  return(randomPermut);


}

