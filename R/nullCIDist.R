#' Computes the distribution of CI under the assumption that all permutations are equally likely.
#"
#' Returns either a probability density distribution or a cumulative distribution (if cumulative=1)
#' over the number of discordant pairs k, with a vector of length choose(n,2)+1 corresponding to 0 
#' through choose(n,2) discordant pairs.  The function assumes that every permutation is equally likely.  
#'
#' Parameters: 
#' @param n {numeric} positive integer, length of the sequence, i.e. the number of elements
#' @param multvect {numeric} a vector of multiplicites of tied elements. If all elements are unique, 
#' this is a vector of ones of length n. Only the counts of tied elements need to be specified, e.g.
#' c(2,3,4) for sets of two, three, and four tied elements
#' @param makeplot {logical} whether to output a plot of the cdf (default=0)
#' @param cumulative {logical} whether to return a probability distribution (FALSE) or cumulative distribution (TRUE) 
#'
#' @export

nullCIDist <- function(n, multvect=c(1), makeplot=0, cumulative=1, return_polylist=0, force_sym=0){
  # Initialize to zeros; by symmetry, it's sufficient to count up to
  # choose(n,2)/2, but this is for thoroughness.  My dist is the probability
  # distribution on the number of inversions [0, choose(n,2)] on n elements.  
  
  ## Example:
  ##   mynull100 <- nullCIDist(n=100)
  if (sum(multvect) < n) {
    multvect <- c(multvect, rep(1, n - sum(multvect)))
  }
  
  # Given the multiplicities in v
  if (sum(multvect == 1) > 0){
    polylist <- list(getSimplePolyProduct(elts=sum(multvect == 1), range=sum(multvect == 1), norm=1))
  } else {
    polylist <- list(c(1))
  }

  multcount <- sum(multvect == 1)
  mvect <- multvect[multvect > 1]
  for (ii in mvect){
    polylist <- c(polylist, getSmartMultiplicity(multcount + ii, ii, norm=1))
    multcount <- multcount + ii
  }
  
  #mydist <- as.numeric(reduce(polylist, mult.p))
  mydist <- as.numeric(mult.plist(polylist))
  # There is a symmetry problem with mydist
  if (force_sym == 1){
    mydist[length(mydist):(length(mydist)/2 + 1)] = mydist[1:(length(mydist)/2)]
  }
  
  if (cumulative == 1){
    mydist <- cumsum(mydist)
  }
  
  if (return_polylist == 1){
    return(polylist) }  #for debugging
  else{
    return(as.numeric(mydist)) }
}

