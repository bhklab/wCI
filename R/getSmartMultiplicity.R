getSmartMultiplicity <- function(elements, multiplicity, norm=1) {
  if (multiplicity == 1) {
    return(getSimplePolyProduct(elements, elements, norm=norm))
  }
  
  numlist <- getSimplePolyList(elements, multiplicity, norm=norm)
  for (ii in multiplicity:2){
    x <- rep(1/ii, ii)
    # This is very slow because of the polynum modulo operation and is very ghetto
    # The operation doesn't check that there exists an element jj for which x is a factor.  
    # Fix this.  Moving the new element to the end 
    for (jj in 1:length(numlist)){
      if (sum(abs(as.numeric(as.polynomial(numlist[[jj]]) %% x))) < 10^-14){
        numlist <- c(numlist, list(as.numeric(divide.p(numlist[[jj]], x))))
        numlist[jj] <- c()
        break
      }
    }
  }

  if (sum(unlist(lapply(numlist, function(x) length(x) - 1))) != (elements - multiplicity) * multiplicity){
    
  }
  return(numlist) #(reduce(numlist,mult.p))
}
