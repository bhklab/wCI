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
#' @param outy {boolean} set to TRUE to not count pairs of predictions that are
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
#' @importFrom sn selm coef.selm psn
#' @import Rcpp
#' @useDynLib wCI _wCI_newPCI
#' @return [list] ! list of concordance index and its pvalue
#' along with the lower and upper confidence intervals
#' @export
#'
paired.concordance.index.new <- function(predictions, observations, delta.pred=0,
                                     delta.obs=0, alpha = 0.05, outx=FALSE, outy=FALSE,
                                     alternative = c("two.sided", "less", "greater"),
                                     logic.operator=c("and", "or"),
                                     CPP=TRUE, 
                                     p_method = c("Permutation", "Asymptotic", "SkewNormal"),
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
  N <- length(predictions)

  if(N < 3){
    return(list("cindex"=NA, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA,
                "relevant.pairs.no"=0))
  }

  p_method <- match.arg(p_method)
  conf_int_method <- match.arg(conf_int_method)
  
  
    values <- newPCI(pin_x=predictions, pin_y=observations, as.numeric(length(predictions)),
                                               pdeltaX=delta.pred, pdeltaY=delta.obs,
                                               pxties=ifelse(outx, 0L, 1L), pyties = ifelse(outx, 0L, 1L),
                                               plogic=ifelse(logic.operator == "and", 1L, 0L))
    C <- values[1]
    D <- values[2]
    CC <- values[3]
    DD <- values[4]
    CD <- values[5]
    N <- values[6]
    # c.d.seq <- values$cdseq
  
  if ((C == 0 && D == 0)) {
    return(list("cindex"=NA, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA,
                "relevant.pairs.no"=0))
  }
  
  returnList <- list("cindex"=NA, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA,
                     "relevant.pairs.no"= (C + D) / 2)
  if(C!=0 & D==0){
    return(list("cindex"=1, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA,
                "relevant.pairs.no"=(C + D) / 2))
  }
  
  if(C==0 & D!=0){
    return(list("cindex"=0, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA,
                "relevant.pairs.no"=(C + D) / 2))
  }
  
  
  cindex <- returnList$cindex <- C / (C + D)
  varp <- 4 * ((D ^ 2 * CC - 2 * C * D * CD + C ^ 2 * DD) / (C + D) ^ 4) * N *
    (N - 1) / (N - 2)
  if (varp >= 0) {
    sterr <- sqrt(varp / N)
    ci <- qnorm(p = alpha / 2, lower.tail = FALSE) * sterr
    p <- pnorm((cindex - 0.5) / sterr)
  } else {
    sterr <- CI <- p <- NA_real_
  }
  if(conf_int_method == "Asymptotic"){
    sterr <- sqrt(varp / N)
    ci <- qnorm(p = alpha / 2, lower.tail = FALSE) * sterr
    returnList$lower <- max(cindex - ci, 0)
    returnList$upper <- min(cindex + ci, 1)
    returnList$sterr <- sterr
  } else if (conf_int_method == "Bootstrap"){
    boot.out <- naiveRCIBoot(x = predictions, y = observations, delta_x = delta.pred, 
                             delta_y = delta.obs, tie.method.x = ifelse(outx, "ignore", "half"), R=boot_num, C=CPP )
    ci.obj <- boot.ci(boot.out, type="bca")
    returnList$lower <- max(ci.obj$bca[4], 0)
    returnList$upper <- min(ci.obj$bca[5], 1)
    returnList$sterr <- sd(boot.out$t[,1])
    returnList$boot.out <- boot.out
  }
  
  
  if(p_method == "Asymptotic"){
    if(C==0 || D==0 || C * (C - 1)==0 || D * (D - 1)==0 || C * D==0 || (C + D) <
       comppairs){
      return(list("cindex"=NA, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA,
                  "relevant.pairs.no"=(C + D) / 2))
    }
    returnList$p.value <- switch(alternative, less=p, greater=1 - p, two.sided=2 *
                                   min(p, 1 - p))
  } else if(p_method == "Permutation"){
    # if(alternative != "two.sided") {warning("Only 2 sided p value currently implemented for permutation.")}
    returnList$p.value <- naiveRCIPerm(x = predictions, y = observations, delta_x = delta.pred, 
                                       delta_y = delta.obs, tie.method.x = ifelse(outx, "ignore", "half"), 
                                       required_alpha = alpha/num_hypothesis/2, p_confidence = perm_p_confidence, C=CPP, alternative = alternative)
  } else if(p_method == "SkewNormal"){
    if(boot_num < 10000){
      warning("At least 10,000 bootstrap samples are recommended to ensure good estimation of sample distribution Skew")
    }
    if(!conf_int_method == "Bootstrap") {# check if bootstrap already exists
      boot.out <- naiveRCIBoot(x = predictions, y = observations, delta_x = delta.pred, 
                             delta_y = delta.obs, tie.method.x = ifelse(outx, "ignore", "half"), R=boot_num )
      returnList$boot.out <- boot.out
    }
    colnames(boot.out$t) <- "obs"
    sn.fit <- selm(obs~1, data=as.data.frame(boot.out$t))
    p <- psn(0.5, dp=coef(sn.fit, "dp"))
    returnList$p.value <- switch(alternative, less=p, greater=1 - p, two.sided=2 *
                                 min(p, 1 - p))

  }
  
  
  
  return(returnList)
}