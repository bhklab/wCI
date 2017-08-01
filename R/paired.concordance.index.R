#' Takes two numerical vectors and computes the concordance index between them
#' by comparing the order of values for two pairs of data each time
#'
#' This function return the concordance index and its p-value
#' along with the lower and upper confidence intervals of said p-value.
#'
#'
#' @examples
#' gdsc.auc <- summarizeSensitivityProfiles(GDSC, sensitivity.measure=‘auc_published’)
#' xx <- summarizeSensitivityProfiles(predicted.by.model, gdsc.auc$Erlotinib)
#'
#' @param predictions {numeric} A vector of predicted drug responces which could be either continuous or discrete
#' @param observations {numeric} A vector of observed continuous drug responces
#' @param cutoff {numeric} A drug responce threshold which is used to classify cells to sensitive vs resistant to drug.
#' @param delta {numeric} The minimunm reliable difference between two drug sensitivity values to be considered as significantly various responses.
#' default value for delta is picked by looking into delta auc values between biological replicates across three
#' large pharmacogenomic studies, CTRPv2(370 drugs over ~15-20 cells) , GDSC(1 drug over ~600 cells), GRAY (85 drugs over ~10-50)
#' @param alpha {numeric} alpha level to compute confidence interval
#' @param outx {boolean} set to TRUE to not count pairs of observations tied on x as a relevant pair.
#' This results in a Goodman-Kruskal gamma type rank correlation.
#' @param alternative {character} what is the alternative hypothesis? Must be one of "two.sides", "less", and "greater".
#' @return [list] ! list of concordance index and its pvalue
#' along with the lower and upper confidence intervals
#' @export

paired.concordance.index <- function(predictions, observations, delta.pred=0.2, delta.obs=0.2, alpha = 0.05, outx = TRUE, alternative = c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)
  predictions[which(is.nan(predictions))] <- NA
  observations[which(is.nan(observations))] <- NA
  cc.ix <- complete.cases(predictions, observations)
  predictions <- predictions[which(cc.ix)]
  observations <- observations[which(cc.ix)]
  N <- length(which(cc.ix))
  c <- d <- u <- matrix(0, nrow = 1, ncol = N)
  for (i in seq(from = 1, to = N - 1)) {
    for (j in seq(from = i + 1, to = N)) {
      pair <- c(i, j)
      #if ((any(predictions[pair] >= cutoff) && abs(predictions[i] - predictions[j]) >= delta) ||
         #(any(observations[pair] >= cutoff) && abs(observations[i] - observations[j]) >= delta)) {
        #if ((predictions[i] == predictions[j] || observations[i] == observations[j])||
        if(abs(observations[i] - observations[j]) >= delta.obs || abs(predictions[i] - predictions[j]) >= delta.pred) { #add flag to replace 'or' behaviour with 'xor' behaviour
          pp <- (predictions[i] < predictions[j])
          oo <- (observations[i] < observations[j])
          if (pp == oo) {
            c[pair] <- c[pair] + 1
          } else {
            d[pair] <- d[pair] + 1
          }

        } else {
          if(outx){
            u[pair] <- u[pair] + 1
          }else{
            d[pair] <- d[pair] + 1
          }
        }
    }
  }

  C <- sum(c)
  D <- sum(d)

  if (N < 3 || (C == 0 && D == 0)) {
    return(list("cindex"=NA, "p.value"=NA, "lower"=NA, "upper"=NA, "relevant.pairs.no"=0))
  }

  cindex <- C / (C + D)

  CC <- sum(c * (c - 1))
  DD <- sum(d * (d - 1))
  CD <- sum(c * d)
  varp <- 4 * (D ^ 2 * CC - 2 * C * D * CD + C ^ 2 * DD) / (C + D) ^ 4 * N * (N - 1) / (N - 2)
  if (varp >= 0) {
    sterr <- sqrt(varp / N)
    ci <- qnorm(p = alpha / 2, lower.tail = FALSE) * sterr
    p <- pnorm((cindex - 0.5) / sterr)
  } else {
    return(list("cindex"=cindex, "p.value"=1, "lower"=0, "upper"=0, "relevant.pairs.no"=(C + D) / 2))
  }
  return(list("cindex"=cindex, 
              "p.value"=switch(alternative, less=p, greater=1 - p, two.sided=2 * min(p, 1 - p)), 
              "lower"=max(cindex - ci, 0), 
              "upper"=min(cindex + ci, 1), 
              "relevant.pairs.no"=(C + D) / 2))
}