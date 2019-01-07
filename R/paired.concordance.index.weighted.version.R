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

paired.concordance.index.weighted.version <- function(predictions, observations,
                                                      delta.pred=.2, delta.obs=.2,
                                                      weightingFun_pred, weightingFun_obs,
                                                      alpha=0.05, outx=FALSE, alternative=c("two.sided", "less", "greater"), logic.operator=c("and", "or"), CPP=TRUE, comppairs=10, permute=FALSE) {
  alternative <- match.arg(alternative)
  logic.operator <- match.arg(logic.operator)
  predictions[which(is.nan(predictions))] <- NA
  observations[which(is.nan(observations))] <- NA
  cc.ix <- complete.cases(predictions, observations)
  predictions <- predictions[which(cc.ix)]
  observations <- observations[which(cc.ix)]
  max_weight <- 1
  obs_dist <- outer(predictions, predictions, FUN="-")
  pred_dist <- outer(observations, observations, FUN="-")
  if(!missing(weightingFun_obs)){
    obs_weights <- abs(log10(weightingFun_obs(obs_dist)))
    if(permute){
      w_order <- sample(1:length(observations))
    }else{
      w_order <- 1:length(observations)
    }
    #obs_weights[which(obs_weights < 0)] <- 0
    if(sum(obs_weights)!=0){
      max_weight <- sum(obs_weights)
    }
  }
  if(!missing(weightingFun_obs) & !missing(weightingFun_pred)){
    pred_weights <- abs(log10(weightingFun_pred(pred_dist)))
    #pred_weights[which(pred_weights < 0)] <- 0
    obs_weights <- obs_weights/sum(obs_weights)
    pred_weights <- pred_weights/sum(pred_weights)
    jj <- sapply(1:length(obs_weights), function(i){max(obs_weights[i], pred_weights[i])})
    if(sum(jj)!=0){
      max_weight <- sum(jj)
    }
  }

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
        if(!missing(weightingFun_obs) & !missing(weightingFun_pred)){
          #w <- sqrt(abs(log(weightingFun_obs(observations[i] - observations[j]))) * abs(log(weightingFun_obs(predictions[i] - predictions[j]))))
          obs_w <- abs(log10(weightingFun_obs(observations[w_order[i]] - observations[w_order[j]])))
          #obs_w <- ifelse(obs_w < 0, 0, obs_w)
          pred_w <- abs(log10(weightingFun_pred(predictions[w_order[i]] - predictions[w_order[j]])))
          #pred_w <- ifelse(pred_w < 0, 0, pred_w)
          w <- 1/max_weight *  max(obs_w, pred_w)
        }else if(!missing(weightingFun_obs)){
          obs_w <- abs(log10(weightingFun_obs(observations[w_order[i]] - observations[w_order[j]])))
          #obs_w <- ifelse(obs_w < 0, 0, obs_w)
          w <- 1/max_weight *  obs_w
        }else{
          w <- 1
        }
        iff <- as.logical(outer(abs(predictions[i] - predictions[j]) > max(delta.pred[i], delta.pred[j]), abs(observations[i] - observations[j]) > max(delta.obs[i], delta.obs[j]), logic.operator))
        if(logic.operator == "&"){
          ife <- abs(predictions[i] - predictions[j]) == max(delta.pred[i], delta.pred[j])
        }else{
          ife <- !iff
        }
        if(iff | !missing(weightingFun_obs)){ #add flag to replace 'or' behaviour with 'xor' behaviour
          pp <- (predictions[i] < predictions[j])
          oo <- (observations[i] < observations[j])
          if (pp == oo) {
            c[pair] <- c[pair] + w
            c.d.seq <- c(c.d.seq, TRUE)
            c.d.seq <- c(c.d.seq, TRUE)
          } else {
            d[pair] <- d[pair] + w
            c.d.seq <- c(c.d.seq, FALSE)
            c.d.seq <- c(c.d.seq, FALSE)
          }
        }else if (ife){
          if(outx){
            u[pair] <- u[pair] + w
          }else{
            d[pair] <- d[pair] + w/2
            c[pair] <- c[pair] + w/2
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

    if(missing(weightingFun_obs)){
      f_obs <- "ignore"
    }else{
      f_obs <- find.original.name(weightingFun_obs)
    }

    if(missing(weightingFun_pred)){
      f_pred <- "ignore"
    }else{
      f_pred <- find.original.name(weightingFun_pred)
    }




    values <- concordanceIndex_modified_helper_weighted(x=predictions, y=observations,
                                               deltaX=delta.pred, deltaY=delta.obs, weightingFun_pred=f_pred, weightingFun_obs=f_obs,
                                               alpha=alpha, outx=outx, alternative=alternative, logicOp=logic.operator, max_weight, max_weight)
    C <- values$C
    D <- values$D
    CC <- values$CC
    DD <- values$DD
    CD <- values$CD
    N <- values$N
    c.d.seq <- values$cdseq
  }

  if (N < 3 || (C == 0 && D == 0)) {
    return(list("cindex"=NA, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA, "relevant.pairs.no"=0))
  }
  if(C!=0 & D==0){
    return(list("cindex"=1, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA, "relevant.pairs.no"=(C + D) / 2, "concordant.pairs"=c.d.seq))
  }
  if(C==0 || D==0 || C * (C - 1)==0 || D * (D - 1)==0 || C * D==0 || (C + D) < comppairs){
    return(list("cindex"=NA, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA, "relevant.pairs.no"=(C + D) / 2, "concordant.pairs"=c.d.seq))
  }
  cindex <- C / (C + D)
  varp <- 4 * ((D ^ 2 * CC - 2 * C * D * CD + C ^ 2 * DD) / (C + D) ^ 4) * N * (N - 1) / (N - 2)
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
              "p.value"=switch(alternative, less=p, greater=1 - p, two.sided=2 * min(p, 1 - p)),
              "sterr"=sterr,
              "lower"=max(cindex - ci, 0),
              "upper"=min(cindex + ci, 1),
              "relevant.pairs.no"=(C + D) / 2,
              "concordant.pairs"=c.d.seq))
}
