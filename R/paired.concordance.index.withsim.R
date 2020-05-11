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
#' @param comppairs {numeric} minimum number of pairs to calculate a valid CI
#' @importFrom stats complete.cases qnorm pnorm
#' @import Rcpp
#' @import parallel
#' @useDynLib wCI _wCI_concordanceIndex_modified_helper
#' @return [list] ! list of concordance index and its pvalue
#' along with the lower and upper confidence intervals
#' @export
#'
paired.concordance.index.withsim <- function(predictions, observations, delta.pred=0,
                                     delta.obs=0, alpha = 0.05, outx=FALSE,
                                     alternative = c("two.sided", "less", "greater"),
                                     logic.operator=c("and", "or"),
                                     CPP=TRUE, comppairs=10,
                                     simulate.p.value = FALSE, ncores=2) {
  alternative <- match.arg(alternative)
  logic.operator <- match.arg(logic.operator)
  predictions[which(is.nan(predictions))] <- NA
  observations[which(is.nan(observations))] <- NA
  cc.ix <- complete.cases(predictions, observations)
  predictions <- predictions[which(cc.ix)]
  observations <- observations[which(cc.ix)]
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
  if(C!=0 & D==0){
    return(list("cindex"=1, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA,
                "relevant.pairs.no"=(C + D) / 2, "concordant.pairs"=c.d.seq))
  }
  if(C==0 || D==0 || C * (C - 1)==0 || D * (D - 1)==0 || C * D==0 || (C + D) <
     comppairs){
    return(list("cindex"=NA, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA,
                "relevant.pairs.no"=(C + D) / 2, "concordant.pairs"=c.d.seq))
  }

  cindex <- C / (C + D)

  if(simulate.p.value != FALSE)
  {
    if(simulate.p.value <100)
    {
      warning("simulate.p.value to low. Setting to 100")
      simulate.p.value = 100
    }

    #sci = unlist(lapply(1:simulate.p.value, function(x)
    #  {
    #  v=paired.concordance.index(sample(predictions), sample(observations),
    #                             delta.pred, delta.obs, alpha, outx, alternative,
    #                             logic.operator, CPP, comppairs, simulate.p.value=FALSE)
    #  return(v$cindex)
    #  }))

    sci <- unlist(parallel::mclapply(1:simulate.p.value, function(x)
    {
      v=paired.concordance.index(sample(predictions), sample(observations),
                                 delta.pred, delta.obs, alpha, outx, alternative,
                                 logic.operator, CPP, comppairs, simulate.p.value=FALSE)
      return(v$cindex)
    }, mc.cores = ncores) )


    if(cindex>0.5)
    {
      left.ci <- sum(sci < (1-cindex))
      rigt.ci <- sum(sci > cindex)
    } else {
      left.ci <- sum(sci < cindex)
      rigt.ci <- sum(sci > (1-cindex))
    }
    pl <- 2*left.ci/length(sci); pr <- 2*rigt.ci/length(sci)
    p <- min(pl, pr)
    if(p==0){p=1/(length(sci)+1)}

    return(list("cindex"=cindex,
                "p.value"=p,
                "sterr"=NA, "lower"=NA, "upper"=NA,
                "relevant.pairs.no"=(C + D) / 2,
                "concordant.pairs"=c.d.seq #, sci=sci
                ))
  }
  varp <- 4 * ((D ^ 2 * CC - 2 * C * D * CD + C ^ 2 * DD) / (C + D) ^ 4) * N *
    (N - 1) / (N - 2)
  if (varp >= 0) {
    sterr <- sqrt(varp / N)
    ci <- qnorm(p = alpha / 2, lower.tail = FALSE) * sterr
    p <- pnorm((cindex - 0.5) / sterr)
  } else {
    return(list("cindex"=cindex,
                "p.value"=1,
                "sterr"=NA,
                "lower"=0,
                "upper"=0,
                "relevant.pairs.no"=(C + D) / 2,
                "concordant.pairs"=c.d.seq))
  }
  return(list("cindex"=cindex,
              "p.value"=switch(alternative, less=p, greater=1 - p, two.sided=2 *
                                min(p, 1 - p)),
              "sterr"=sterr,
              "lower"=max(cindex - ci, 0),
              "upper"=min(cindex + ci, 1),
              "relevant.pairs.no"=(C + D) / 2,
              "concordant.pairs"=c.d.seq))
}