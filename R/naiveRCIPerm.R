# count_xties and count_yties are booleans that indicate whether ties should count as 0.5 correctly ordered 
# and 0.5 incorrectly ordered pairs.  if count_*ties = 0, ties are ignored. 

# This function has not yet been tested.

### TODO:: Need to fix this to only do rCI. Now sure how to generalize to pearson/spearman, think about that!
### XXX: this doesn't do anything correctly right now!!


source("adaptive_permutation.R")



computePearson <- function(mat) return(cor(mat[,1], mat[,2]))

naiveRCIPerm <- function(x, y, 
                         delta_x = 0.2, 
                         delta_y = 0.2, 
                         valid.logic = c("and", "or"),
                         tie.method.x = c("ignore", "half"), 
                         tie.method.y = c("ignore", "half"),
                         required_alpha = 1e-6,
                         p_confidence = 0.01){
  
  valid.logic <- match.arg(valid.logic)
  tie.method.x = match.arg(tie.method.x)
  tie.method.y = match.arg(tie.method.y)
  N <- length(x)
  
  seed <- .Random.seed
   
  myCompleteCases <- complete.cases(x,y)
  x <- x[myCompleteCases]
  y <- y[myCompleteCases]
  
  ### Calculating Permutation stopping parameters
  
  B <- choose_b(required_alpha, p_confidence)
  R <- choose_r(required_alpha, p_confidence)
  
  
  
  
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
  
  t <- numeric(B)
  
  totalSeenLarger <- 0
  
  i <- 1
  while(i <= B){
    smpl <- sample.int(N)
    ymat <- ymat[smpl, smpl]
    t[i] <- compCI(xmat * ymat)
    totalSeenLarger <- totalSeenLarger + (abs(t[i]-0.5) > abs(t0-0.5))
    if(totalSeenLarger == R){
      break
    }
    i = i + 1
  }
  browser()
  if(any(totalSeenLarger == R)){
    return(totalSeenLarger/i)
  } else {
    return((totalSeenLarger+1)/(B + 1))
  }
  
  return(res)
}
