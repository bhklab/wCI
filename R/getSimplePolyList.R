getSimplePolyList <- function(elts, range, norm=0){
  if (range > elts) {
    return("Error: Range (or number of polynomials) is greater than the input number of elements")}
  if (range < 1 | range != round(range)) {
    return("Error: Range must be a positive integer")}
  a <- list()

  for (ii in (elts - range + 1):elts){
    k <- ifelse(norm != 0, 1/ii, 1)
    a <- c(a, list(rep(k,ii)))
  }
  return(a)
}
