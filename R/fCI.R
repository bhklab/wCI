################################################
## These functions are for computing modified concordance index and simulated
## p-value. This function is optimized for performance in R. An implementation
## in C/C++ might increase speed.
#################################################

getDiffMat <- function(x, dx=0)
  {
    sortx <- names(sort(x))
    dmat <- matrix(data = 0, nrow = length(sortx), ncol = length(sortx))
    rownames(dmat) <- colnames(dmat) <- sortx
    lsor <- length(sortx)
    for(i  in 1:(lsor-1))
    {
      xsi = x[sortx[i]]
      for (j in (i+1):lsor)
      {
        sj = sortx[j]
        if((x[sj] - xsi) > dx)
        {
          dmat[i, j:lsor] = -1
          dmat[j:lsor, i] = 1
          break()
        }
      }
    }
    return(dmat)
  }

#' Takes two numerical vectors and computes the concordance index between them
#' by comparing the order of values for two pairs of data each time
#'
#' This function returns the concordance index and simulated p-value.
#'
#' @param x a numeric vector
#' @param y a numeric vector
#' @param dx The minimunm reliable difference between two values in
#' the x vector to be considered as significantly
#' @param dy The minimunm reliable difference between two values in
#' the y vector to be considered as significantly
#' @param n.simulate number of simulation for p value computation, default 1000.
#' If set to 0 only CI will be return
#' @param mc.cores number of cores for parallel processing, default 1
#'
#' @return a list of concordance index and its pvalue
#'
#' @examples
#' data("iris")
#' fwci(x=iris$Sepal.Length, y=iris$Sepal.Width, dx=0.1, dy=0.5)
#'
#' @author{ Arvind Mer }
#' @import dqrng
#' @import parallel
#' @export
#'
fwci <- function(x, y, dx=0, dy=0, n.simulate=1000, mc.cores=1)
{
  if(n.simulate > 0 & n.simulate<1000)
  {
    warning("n.simulate is low, p value might not be reliable")
  }
  names(x) = namesX = paste0("x", 1:length(x))
  names(y) = namesY = paste0("y", 1:length(y))

  dfx = getDiffMat(x, dx)
  dfy = getDiffMat(y, dy)
  dz = dfx[namesX, namesX] * dfy[namesY, namesY]
  conp= sum(dz==  1)/2
  disp= sum(dz== -1)/2
  cindex = conp/(conp+disp)

  p=NA
  if(n.simulate>0)
  {
    ns <- ceiling(n.simulate/mc.cores)

    sz=nrow(dfy)
    #profvis::profvis({
    sci <- unlist(parallel::mclapply(1:mc.cores, function(i)
    {
      #profvis::profvis({
      lapply(1:ns, function(si)
      {
        #rx=dqsample.int(sz,size=sz);
        ry=dqsample.int(sz,size=sz)
        dz = dfx * dfy[ry,ry]
        conp= sum(dz==  1)/2
        disp= sum(dz== -1)/2
        conp/(conp+disp)
      })
      #})
    }, mc.cores = mc.cores))
    #})

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
  }

  return(list("cindex"=cindex,
              "p.value"=p,
              "total.relevant.pairs"= (conp+disp),
              "concordant.pairs"=conp))
}


