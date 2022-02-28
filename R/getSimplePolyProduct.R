getSimplePolyProduct <- function(elts, range, norm=0){
  if (range > elts){
    return("Error: Range (or number of polynomials) is greater than the input number of elements")}
  if (range < 1 | range != round(range)){
    return("Error: Range must be a positive integer")}

  #initialize - initial length is sum of (elts-1):(elts-range) + 1
  #This initialization is a bit gross, but it's faster to start at the bottom of the range
  retpoly <- c(rep(1, (elts-range+1)), rep(0, range/2 * (2*elts-range-1) + 1 - (elts-range+1))) #* 1/(elts-range+1)
  count <- (elts - range + 1)
  if (range > 1){
    for (k in (elts-range+2):elts){
      retpoly[1:(count+k)] = (cumsum(c(retpoly[1:count], rep(0,k))) - cumsum(c(rep(0,k), retpoly[1:count])))# * 1/k
      count <- count + k-1
      # Normalize stepwise, floor out guys at 10^-200?
    }
  }
  if (norm != 0){
    retpoly <- retpoly / sum(retpoly)
  }
  return(retpoly[1:count])
}
