#' Takes two numerical vectors and computes the concordance index between them
#' by comparing the order of values for two pairs of data each time
#'
#' This function returns the concordance index and its p-value
#' along with the lower and upper confidence intervals of said p-value.
#'
#'
#' @examples
#' data(PLX4720_data)
#' pci_PLX4720 <- paired.concordance.index(predictions = PLX4720_data[ ,"AAC_CTRPv2"],
#' observations = PLX4720_data[ ,"AAC_GDSC"], delta.pred = 0, delta.obs = 0,
#' outx = TRUE)
#' pci_PLX4720$cindex
#'
#' @param predictions {numeric} A vector of predicted drug responces which could
#' be either continuous or discrete
#' @param observations {numeric} A vector of observed continuous drug responces
#' @param delta.pred {numeric} The minimunm reliable difference between two
#' values in the predictions vector to be considered as significantly various
#' values.
#' @param delta.obs {numeric} The minimunm reliable difference between two
#' values in the observations vector to be considered as significantly various
#' values. In drug sensitivity , default value for delta.pred is picked by
#' looking into delta auc values (drug response metric) between biological
#' replicates across three large pharmacogenomic studies,
#' CTRPv2 (370 drugs over ~15-20 cells), GDSC (1 drug over ~600 cells),
#' GRAY (85 drugs over ~10-50 cells)
#' @param alpha {numeric} alpha level to compute confidence interval
#' @param outx {boolean} set to TRUE to not count pairs of predictions that are
#' tied as a relevant pair. This results in a Goodman-Kruskal gamma type rank
#' correlation.
#' @param alternative {character} What is the alternative hypothesis? Must be
#' one of "two.sides", "less", and "greater" and defaults to two.sides".
#' @param logic.operator {character} determines how strict should the test be to
#' remove noisy pairs. Must be one of "and" or "or" and defaults to "and".
#' @param CPP {boolean} Whether to use the C version of the code for faster
#' execution
#' @param p_method {character} Either "Permutation", or "Asymptotic", picks a
#' method to use for calculating p-values. If Permutation, then "alpha"/"num_hypothesis"
#' is used to determine the effective alpha used for estimating number of required 
#' permutations. 
#' @param conf_int_method {character} Either "Bootstrap" or "Asymptotic", picks a method for
#' estimating the confidence interval corresponding to 1-"alpha". 
#' @param num_hypothesis {numeric} Total number of hypothesis being tested in analysis. Used
#' for adjusting number of required permutations when using the permutation method of computing
#' p values. Default 1. Ignored if using asymptotic p value.
#' @param perm_p_confidence {numeric} Maximum permited 1 SD confidence interval of our estimated 
#' permutation p value around the true p value, as a fraction of "alpha"/"num_hypothesis". Ignored
#' if using asymptotic p value, no guarantee on correctness exists. 
#' @param boot_num {numeric} number of samples to use for bootstrap. Default 5000. Ignored
#' if using asymptotic confidence interval. 
#' @param comppairs {numeric} minimum number of pairs to calculate a valid CI.
#' @importFrom stats complete.cases qnorm pnorm 
#' @importFrom boot boot.ci
#' @import Rcpp
#' @useDynLib wCI _wCI_concordanceIndex_modified_helper
#' @return [list] ! list of concordance index and its pvalue
#' along with the lower and upper confidence intervals
#' @export
#'
paired.concordance.index <- function(predictions, observations, delta.pred=0,
                                     delta.obs=0, alpha = 0.05, outx=FALSE,
                                     alternative = c("two.sided", "less", "greater"),
                                     logic.operator=c("and", "or"),
                                     CPP=TRUE, 
                                     p_method = c("Permutation", "Asymptotic"),
                                     conf_int_method = c("Bootstrap", "Asymptotic"),
                                     num_hypothesis = 1,
                                     perm_p_confidence = 0.2,
                                     boot_num  = 5000,
                                     comppairs=10) {
  alternative <- match.arg(alternative)
  logic.operator <- match.arg(logic.operator)
  predictions[which(is.nan(predictions))] <- NA
  observations[which(is.nan(observations))] <- NA
  cc.ix <- complete.cases(predictions, observations)
  predictions <- predictions[which(cc.ix)]
  observations <- observations[which(cc.ix)]
  
  p_method <- match.arg(p_method)
  conf_int_method <- match.arg(conf_int_method)
  
  
  if(!CPP){
    logic.operator <- ifelse(logic.operator=="or", "|", "&")
    N <- length(which(cc.ix))
    if(length(delta.pred) == 1){
      delta.pred <- rep(delta.pred, N)
    }else{
      delta.pred <- delta.pred[which(cc.ix)]
    }
    if(length(delta.obs) == 1){
      delta.obs <- rep(delta.obs, N)
    }else{
      delta.obs <- delta.obs[which(cc.ix)]
    }
    c <- d <- u <- matrix(0, nrow = 1, ncol = N)
    c.d.seq <- NULL
    for (i in seq(from = 1, to = N - 1)) {
      for (j in seq(from = i + 1, to = N)) {
        pair <- c(i, j)
        iff <- as.logical(outer(abs(predictions[i] - predictions[j]) >
                                  max(delta.pred[i], delta.pred[j]),
                                abs(observations[i] - observations[j]) >
                                  max(delta.obs[i], delta.obs[j]), logic.operator))
        if(logic.operator == "&"){
          ife <- abs(predictions[i] - predictions[j]) == max(delta.pred[i],
                                                             delta.pred[j])
        }else{
          ife <- !iff
        }
        if(iff){ #add flag to replace 'or' behaviour with 'xor' behaviour
          pp <- (predictions[i] < predictions[j])
          oo <- (observations[i] < observations[j])
          if (pp == oo) {
            c[pair] <- c[pair] + 1
            c.d.seq <- c(c.d.seq, TRUE)
            c.d.seq <- c(c.d.seq, TRUE)
          } else {
            d[pair] <- d[pair] + 1
            c.d.seq <- c(c.d.seq, FALSE)
            c.d.seq <- c(c.d.seq, FALSE)
          }
        }else if (ife){
          if(outx | abs(observations[i] - observations[j]) <= max(delta.obs[i],
                                                                  delta.obs[j])){
            u[pair] <- u[pair] + 1
          }else{
            d[pair] <- d[pair] + 0.5
            c[pair] <- c[pair] + 0.5
            c.d.seq <- c(c.d.seq, TRUE)
            c.d.seq <- c(c.d.seq, FALSE)
          }
        }
      }
    }
    C <- sum(c)
    D <- sum(d)
    CC <- sum(c * (c - 1))
    DD <- sum(d * (d - 1))
    CD <- sum(c * d)
  }else{
    values <- concordanceIndex_modified_helper(x=predictions, y=observations,
                                               deltaX=delta.pred, deltaY=delta.obs,
                                               alpha=alpha, outx=outx,
                                               alternative=alternative,
                                               logicOp=logic.operator)
    C <- values$C
    D <- values$D
    CC <- values$CC
    DD <- values$DD
    CD <- values$CD
    N <- values$N
    c.d.seq <- values$cdseq
  }
  
  if (N < 3 || (C == 0 && D == 0)) {
    return(list("cindex"=NA, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA,
                "relevant.pairs.no"=0))
  }
  
  returnList <- list("cindex"=NA, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA,
                     "relevant.pairs.no"= (C + D) / 2, "concordant.pairs"=c.d.seq)
  if(C!=0 & D==0){
    return(list("cindex"=1, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA,
                "relevant.pairs.no"=(C + D) / 2, "concordant.pairs"=c.d.seq))
  }

  if(C==0 & D!=0){
    return(list("cindex"=0, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA,
                "relevant.pairs.no"=(C + D) / 2, "concordant.pairs"=c.d.seq))
  }
  
  
  cindex <- returnList$cindex <- C / (C + D)
  varp <- 4 * ((D ^ 2 * CC - 2 * C * D * CD + C ^ 2 * DD) / (C + D) ^ 4) * N *
    (N - 1) / (N - 2)
  if (varp >= 0) {
    sterr <- sqrt(varp / N)
    ci <- qnorm(p = alpha / 2, lower.tail = FALSE) * sterr
    p <- pnorm((cindex - 0.5) / sterr)
  }

  if(p_method == "Asymptotic"){
    if(C==0 || D==0 || C * (C - 1)==0 || D * (D - 1)==0 || C * D==0 || (C + D) <
       comppairs){
      return(list("cindex"=NA, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA,
                  "relevant.pairs.no"=(C + D) / 2, "concordant.pairs"=c.d.seq))
    }
    returnList$p.value <- switch(alternative, less=p, greater=1 - p, two.sided=2 *
                                   min(p, 1 - p))
  } else if(p_method == "Permutation"){
    if(alternative != "two.sided") {warning("Only 2 sided p value currently implemented for permutation.")}
    returnList$p.value <- naiveRCIPerm(x = predictions, y = observations, delta_x = delta.pred, 
                 delta_y = delta.obs, tie.method.x = ifelse(outx, "ignore", "half"), 
                 required_alpha = alpha/num_hypothesis, p_confidence = perm_p_confidence, C=CPP)
  }
  
 if(conf_int_method == "Asymptotic"){
   sterr <- sqrt(varp / N)
   ci <- qnorm(p = alpha / 2, lower.tail = FALSE) * sterr
   returnList$lower <- max(cindex - ci, 0)
   returnList$upper <- min(cindex + ci, 1)
   returnList$sterr <- sterr
 } else if (conf_int_method == "Bootstrap"){
   boot.out <- naiveRCIBoot(x = predictions, y = observations, delta_x = delta.pred, 
                            delta_y = delta.obs, tie.method.x = ifelse(outx, "ignore", "half"), R=boot_num)
   ci.obj <- boot.ci(boot.out, type="bca")
   returnList$lower <- max(ci.obj$bca[4], 0)
   returnList$upper <- min(ci.obj$bca[5], 1)
 }
  

  
  return(returnList)
}