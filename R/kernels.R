kernel_gaussian <- function(x, m=0.0002037366, s=0.0919937995){
  #dnorm(x, m, s)
  (1/sqrt(2*pi*s^2))*exp(-(x-m)^2/(2*s^2))
}
kernel_laplace <- function(x, m=-0.001630207, b=0.060597464){
  (1/(2*b)) * exp(-abs(x - m) / b)
}


find.original.name <- function(fun) {
  objects <- ls(envir = environment(fun))
  for (i in objects) {
    if (identical(fun, get(i, envir = environment(fun)))) {
      return(i)
    }
  }
}
