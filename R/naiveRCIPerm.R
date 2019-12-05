# count_xties and count_yties are booleans that indicate whether ties should count as 0.5 correctly ordered 
# and 0.5 incorrectly ordered pairs.  if count_*ties = 0, ties are ignored. 

# This function has not yet been tested.

### TODO:: Need to fix this to only do rCI. Now sure how to generalize to pearson/spearman, think about that!
### XXX: this doesn't do anything correctly right now!!


library(tictoc)


computePearson <- function(mat) return(cor(mat[,1], mat[,2]))


#' Naive RCI Permutations
#' 
#' @useDynLib wCI _wCI_rCIPermC
naiveRCIPerm <- function(x, y, 
                         delta_x = 0.2, 
                         delta_y = 0.2, 
                         valid.logic = c("and", "or"),
                         tie.method.x = c("ignore", "half"), 
                         tie.method.y = c("ignore", "half"),
                         alternative = c("two.sided", "greater", "less"),
                         required_alpha = 1e-6,
                         p_confidence = 0.01,
                         C=TRUE,
                         verbose=FALSE){
  
  valid.logic <- match.arg(valid.logic)
  tie.method.x = match.arg(tie.method.x)
  tie.method.y = match.arg(tie.method.y)
  alternative = match.arg(alternative)
  alternative <- switch (alternative,
    "two.sided" = 0L,
    "greater" = 1L,
    "less" = -1L
  )
  
  seed <- .Random.seed
   
  myCompleteCases <- complete.cases(x,y)
  x <- x[myCompleteCases]
  y <- y[myCompleteCases]
  
  N <- length(x)
  
  ### Calculating Permutation stopping parameters
  
  B <- choose_b(required_alpha, p_confidence)
  R <- choose_r_fast(required_alpha, p_confidence)
  
  
  xmat <- matrix(rep(x, length(x)), ncol=length(x))
  ymat <- matrix(rep(y, length(y)), ncol=length(y))
  
  # Deltamat_ij = x_j - x_i 
  xdeltamat <- xmat - t(xmat)
  ydeltamat <- ymat - t(ymat)
  
  # RCI: sign(xdeltamat) * (abs(xdeltamat) > threshold)
  if (valid.logic == "and"){
    xmat <- sign(xdeltamat) * (abs(xdeltamat) > delta_x) 
    ymat <- sign(ydeltamat) * (abs(ydeltamat) > delta_y)
  } else {
    stop("Not Implemented")
  }
  compCI <- function(rcimat){
    if(tie.method.x == tie.method.y & tie.method.y == "ignore"){
      t0 <- sum(rcimat>0)/sum(rcimat != 0)
    } else if (tie.method.x == tie.method.y & tie.method.y == "half"){
      t0 <- (sum(rcimat>0) + 0.5*sum(rcimat==0)-N/2)/(choose(N,2)*2)
    } else {
      if (tie.method.x == "half" & tie.method.y == "ignore"){
        tieNum <- sum(abs(xdeltamat) <= delta_x) - N
      }
      if (tie.method.x == "ignore" & tie.method.y == "half"){
        tieNum <- sum(abs(ydeltamat) <= delta_y) - N
      }
      t0 <- (sum(rcimat>0) + 0.5*tieNum)/(sum(rcimat!=0) + tieNum)
    }
    return(t0)
  }
  
  t0 <- compCI(xmat * ymat)
  if(C){
    if(valid.logic == "or"){
      stop("Not Implemented Yet in C. Please use C=FALSE, but warning, very slow!")
    }
    pval <- .Call("_wCI_rCIPermC", PACKAGE = "wCI",
                  as.integer(xmat),
                  as.integer(ymat),
                  as.numeric(t0),
                  as.numeric(R),
                  as.numeric(B),
                  as.numeric(N),
                  as.integer(tie.method.x == "half"),
                  as.integer(tie.method.y == "half"), 
                  as.integer(alternative),
                  as.numeric(runif(2))
                  )
    if((getOption("verbose") || verbose ) & pval[2]){
      message("Permutations reached early stopping condition")
    }
    return(pval[1])
  } else {
    totalSeenLarger <- 0
    i <- 1
    while(i <= B){
      smpl <- sample.int(N)
      ymat <- ymat[smpl, smpl]
      t <- compCI(xmat * ymat)
      totalSeenLarger <- totalSeenLarger + (abs(t-0.5) > abs(t0-0.5))
      if(totalSeenLarger == R){
        break
      }
      i = i + 1
    }
    if(any(totalSeenLarger == R)){
      if((getOption("verbose") || verbose )){
        message("Permutations reached early stopping condition")
      }
      return(totalSeenLarger/i)
    } else {
      return((totalSeenLarger+1)/(B + 1))
    }
  }
}
