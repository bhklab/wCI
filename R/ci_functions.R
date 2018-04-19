
paired.concordance.index.usingC <-function(x,y, deltaX=0, deltaY=0, alpha =0, outx = 1, alternative = c("two.sided", "less", "greater"), logic.operator=c("and", "or"))
{
  values <- concordanceIndex_modified_helper(x = x, y = y,
                                      deltaX = deltaX, deltaY = deltaY,
                                      alpha = alpha, outx = outx,alternative=alternative,logicOp = logic.operator)

  C <- values$C
  D <- values$D
  CC <- values$CC
  DD <- values$DD
  CD <- values$CD
  N <- values$N
  c.d.seq <- values$cdseq
  cindex <- C / (C + D)
  varp <- 4 * ((D ^ 2 * CC - 2 * C * D * CD + C ^ 2 * DD) / (C + D) ^ 4) * N * (N - 1) / (N - 2)
  if (varp >= 0) {
    sterr <- sqrt(varp / N)
    ci <- qnorm(p = alpha / 2, lower.tail = FALSE) * sterr
    p <- pnorm((cindex - 0.5) / sterr)
  } else {
    return(list("cindex"=cindex, "p.value"=1, "lower"=0, "upper"=0, "relevant.pairs.no"=(C + D) / 2, "concordant.pairs"=c.d.seq))
  }
  return(list("cindex" = cindex, "p.value" = switch(alternative, less = p, greater = 1 - p, two.sided = 2 * min(p, 1 - p)), "lower" = max(cindex - ci, 0), "upper" = min(cindex + ci, 1), "relevant.pairs.no" = (C + D) / 2, "concordant.pairs"=c.d.seq))
}




.computePval <- function(cindex, permutCI)
{
  pv <- NA
  if(cindex >= 0)
  {
    pv <- sum(permutCI[ permutCI >= cindex ])/ length(permutCI)
  } else
  {
    pv <- abs(sum(permutCI[ permutCI < cindex ])/ length(permutCI))
  }
  return(pv)
}



CIinC <-function(x,y, deltaX=0, deltaY=0, alpha =0, outx = 1, npermut=1000000, ncpu=1)
{
  cindex <- concordanceIndex_modified(x = x, y = y,
                                      deltaX = deltaX, deltaY = deltaY,
                                      alpha = alpha, outx = outx)

#  permutCI <- permute_concordanceIndex_modified(x = x, y = y,
#                                                deltaX = deltaX, deltaY = deltaY,
#                                                alpha = alpha, outx = outx,
#                                                permutations = npermut, nThread=ncpu)
  pv <- .computePval(cindex, permutCI)
  return(list(ci=cindex, p.value=pv))
}


#' concordance index function
#'
#' To compute concordance index and p-value for given vector
#' @param x first vector
#' @param y second vector
#' @param deltaX
#' @param deltaY
#' @param alpha
#' @param outx
#' @param npermut default 1000000, number of permutations to run to calculate p-value
#' @param ncpu number of cpu to use
#' @return a list containg concordance index and p-value
#'
#' @examples
#' library(CI)
#' x <- 1:30; y <- x*x
#' CI(x=x, y=y, deltaX=0, deltaY=0, alpha =0, outx = 1, npermut=10000)
#' @importFrom foreach foreach %dopar%
#' @export
CI <-function(x,y, deltaX=0, deltaY=0, alpha=0, outx = 1, npermut=1000000,
              ncpu=1)
{
  cindex <- concordanceIndex_modified(x = x, y = y,
                                      deltaX = deltaX, deltaY = deltaY,
                                      alpha = alpha, outx = outx)
  if(ncpu>1)
  {
    cl <- parallel::makeCluster(ncpu)
    doParallel::registerDoParallel(cl)
    permutCI <-foreach(i=1:npermut) %dopar% {
      xs <- sample(x); ys <- sample(y)
      concordanceIndex_modified(x = xs, y = ys,
                                deltaX = deltaX, deltaY = deltaY,
                                alpha = alpha, outx = outx)
    }
    parallel::stopCluster(cl)
  } else
  {
    permutCI <- sapply(1:npermut, function(i){
                       xs <- sample(x); ys <- sample(y)
                       concordanceIndex_modified(x = xs, y = ys,
                                                 deltaX = deltaX, deltaY = deltaY,
                                                 alpha = alpha, outx = outx)
      })

  }
  permutCI <- unlist(permutCI)
  pv <- .computePval(cindex, permutCI)
  return(list(ci=cindex, p.value=pv))
}


testCI <- function()
{
  library(CI)
  library(MASS)
  cr <- -0.298
  out <- mvrnorm(50, mu = c(0,0), Sigma = matrix(c(1,cr,cr,1), ncol = 2),
                 empirical = TRUE)

  cor.test(out[,1], out[,2])
  cor.test(out[,1], out[,2], method = "pearson")
  cor.test(out[,1], out[,2], method = "spearman")
  #x=out[,1]; y=out[,2];deltaX=0;deltaY=0;alpha =0;outx=1;npermut=10000
  CI(x=out[,1], y=out[,2], deltaX=0, deltaY=0, alpha =0, outx = 1, npermut=10000,
     ncpu = 2)

  CIinC(out[,1], out[,2], deltaX=0, deltaY=0, alpha =0, outx = 1, npermut=10000,ncpu = 2)
}



