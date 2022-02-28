mult.plist <-
function(plist, outOrder){
  if (missing(outOrder)) {
    outOrder <- sum(unlist(lapply(plist, length))) - length(plist) + 1
  }
  if (length(plist) == 1){
    return(plist[[1]])
  }
  
  product <- fft(pad(plist[[1]], outOrder))
  for (ii in 2:length(plist)){
    product <- product * fft(pad(plist[[ii]], outOrder))
  }
  
  return(fft(1/length(product) * product, inverse = TRUE))
}
