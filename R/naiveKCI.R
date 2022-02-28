#' Calculation of the kernelized concordance index (kCI)
#' 
#' This function takes as input two numeric vectors and calculates the kernelized 
#' concordance index. 
#' 
#' @param x, y {numeric} Vectors of real numbers
#' @param mykernel {function} A univariate kernel function on the distance between two points, defaults to 
#' 1/(1+exp(-27.5512 * abs(x - 0.08))), as fit on replicates from pharmacogenomic data from PharmacoGx
#' @return [list] A list of the kCI and a number of fields output by other concordance index versions. 
#' @export

naiveKCI <- function(x, y, 
                      mykernel = function(x) (1/(1+exp(-27.5512 * (abs(x) - 0.0800)))), 
                      compute.p = c(TRUE, FALSE), 
                      alternative=c("two.sided", "greater", "less"),
                      p.method=c(),
                      alpha=0.05, 
                      interval=c("confidence", "prediction"), 
                      returnAll=c(FALSE, TRUE)){  

  myCompleteCases <- complete.cases(x,y)
  x <- x[myCompleteCases]
  y <- y[myCompleteCases]
  
  xmat <- matrix(rep(x, length(x)), ncol=length(x))
  ymat <- matrix(rep(y, length(y)), ncol=length(y))
 
  # Deltamat_ij = x_j - x_i 
  xdeltamat <- xmat - t(xmat)
  ydeltamat <- ymat - t(ymat)
  
  # RCI: sign(xdeltamat) * (abs(xdeltamat) > threshold)
  if (altkmat){
    kmat <- altkmat
  } else {
    kmat <- mykernel(abs(xdeltamat)) * mykernel(abs(ydeltamat))
  }
  kcimat <- sign(xdeltamat) * sign(ydeltamat) * kmat
  cindex <- (1 + sum(kcimat)/sum(abs(kcimat)))/2

  return(list(cindex=cindex, 
                p.value=1, 
                sterr=NA, 
                lower=NA, 
                upper=NA, 
                relevant.pairs.no=sum(abs(kcimat))/2))
}
