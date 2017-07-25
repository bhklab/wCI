
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



CIinC <-function(x,y, deltaX=0, deltaY=0, alpha =0, outx = 1, npermut=1000000)
{
  cindex <- concordanceIndex_modified(x = x, y = y,
                                      deltaX = deltaX, deltaY = deltaY,
                                      alpha = alpha, outx = outx)

  permutCI <- permute_concordanceIndex_modified(x = x, y = y,
                                                deltaX = deltaX, deltaY = deltaY,
                                                alpha = alpha, outx = outx,
                                                permutations = npermut)
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
    stopCluster(cl)
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

  CIinC(out[,1], out[,2], deltaX=0, deltaY=0, alpha =0, outx = 1, npermut=10000)
}



