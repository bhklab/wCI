mult.p <-
function(p1, p2, outOrder){
  if(missing(outOrder)){
    outOrder <- length(p1) + length(p2) - 1
  } 
  p1 <- pad(p1, outOrder)
  p2 <- pad(p2, outOrder)
  
  p1.fft <- fft(p1)
  p2.fft <- fft(p2)
  p3.fft <- p1.fft*p2.fft
  
  return(fft(1/length(p3.fft) * p3.fft, inverse = TRUE))
}
