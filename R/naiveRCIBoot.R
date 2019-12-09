# count_xties and count_yties are booleans that indicate whether ties should count as 0.5 correctly ordered 
# and 0.5 incorrectly ordered pairs.  if count_*ties = 0, ties are ignored. 

# This function has not yet been tested.

naiveRCIBoot <- function(x, y, 
                     delta_x = 0.2, 
                     delta_y = 0.2, 
                     valid.logic = c("and"),
                     tie.method.x = c("ignore", "half"), 
                     tie.method.y = c("ignore", "half"),
                     R = 2000, C=T){
  
  valid.logic <- match.arg(valid.logic)
  tie.method.x = match.arg(tie.method.x)
  tie.method.y = match.arg(tie.method.y)
  
  seed <- .Random.seed
  
  myCompleteCases <- complete.cases(x,y)
  x <- x[myCompleteCases]
  y <- y[myCompleteCases]
  N <- length(x)
  
  xmat <- matrix(rep(x, length(x)), ncol=length(x))
  ymat <- matrix(rep(y, length(y)), ncol=length(y))
  
  # Deltamat_ij = x_j - x_i 
  xdeltamat <- xmat - t(xmat)
  ydeltamat <- ymat - t(ymat)
  
  # RCI: sign(xdeltamat) * (abs(xdeltamat) > threshold)
  if (valid.logic == "and"){
    rcimat <- sign(xdeltamat) * sign(ydeltamat) * ((abs(xdeltamat) > delta_x) & (abs(ydeltamat) > delta_y)) 
    if(tie.method.x == tie.method.y & tie.method.y == "ignore"){
    } else if (tie.method.x == tie.method.y & tie.method.y == "half"){
      rcimat <- rcimat + 0.5*rcimat==0
      diag(rcimat) <- 0
    } else {
      if (tie.method.x == "half" & tie.method.y == "ignore"){
          rcimat <- rcimat + (abs(xdeltamat) <= delta_x & abs(ydeltamat) > delta_y) * 0.5
      }
      if (tie.method.x == "ignore" & tie.method.y == "half"){
          rcimat <- rcimat + (abs(ydeltamat) <= delta_y & abs(xdeltamat) > delta_x) * 0.5
      }
    }
  } else {
    stop("Valid logic or not supported anymore. ")
  }
  
  compCI <- function(rcimat){
    t0 <- sum(rcimat[rcimat > 0])/(sum(rcimat != 0))
    return(t0)
  }
  
  t0 <- compCI(rcimat)
  
  if(!C){
    t <- numeric(R)
  
    for(i in seq_along(t)){
      smpl <- sample.int(N,N,replace=TRUE)
      rcimatSample <- rcimat[smpl, smpl]
      t[i] <- compCI(rcimatSample)
    }
  } else {

    t <- .Call("_wCI_rCIBootC", PACKAGE = "wCI",
                  as.numeric(rcimat),
                  as.numeric(R),
                  as.numeric(N),
                  as.integer(tie.method.x == "half"),
                  as.integer(tie.method.y == "half"), 
                  as.numeric(runif(2))
                  )

  }
  
  


  dim(t) <- c(R,1)
  res <- list(
    t0 = t0,
    t = t,
    R = R,
    data = cbind(x,y),
    seed = seed,
    sim = "ordinary",
    stype = "i",
    call = match.call(),
    strata = rep(1, N),
    weights = rep(1/N, N)
  )
  attr(res, "class") <- "boot"
  attr(res, "boot_type") <- "boot"
  return(res)
}
