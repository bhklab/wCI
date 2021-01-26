dyn.load("naivePermGPU.so")
system.time(t <- .Call("permCUDA", as.numeric(1:500), runif(500), 1e5, 500, 0, 0, runif(1)))
