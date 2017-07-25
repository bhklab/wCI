

# parallelizing using R
concordanceIndex.modified <- function(x,y,deltaX,deltaY,alpha=0.05,outx=1,permutations,numCores=4){
  require(parallel)
  
  cindex <- concordanceIndex_modified(x = x,y,deltaX,deltaY,alpha,outx)
  
  C <- unlist(mclapply(seq(1:permutations),function(x){
    tmp <- sample(x)
    concordanceIndex_modified(tmp,y,deltaX,deltaY,alpha,outx)
  },mc.cores = numCores))
  
  psedueE <- 1e-200
  pvalue <- (sum(C>cindex,na.rm = T)+psedueE)/permutations
  
  return(c("cIndex"=cindex,"pvalue"=pvalue))
}


# method in C
concordanceIndex.modified_allC <- function(x,y,deltaX,deltaY,alpha=0.05,outx=1,permutations){
  
  cindex <- concordanceIndex_modified(x = x,y,deltaX,deltaY,alpha,outx)
  
  C <- permute_concordanceIndex_modified(x = x,y = y,deltaX = deltaX,deltaY = 0.01,alpha = alpha,outx = outx,permutations = permutations)
  psedueE <- 1e-200
  pvalue <- (sum(C>cindex,na.rm = T)+psedueE)/permutations
  
  return(c("cIndex"=cindex,"pvalue"=pvalue))
}
