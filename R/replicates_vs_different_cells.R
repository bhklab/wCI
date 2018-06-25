library(PharmacoGx)
load("Psets/CTRPv2.RData")
load("results/aac_param_bs.RData", verbose=T)
pdf("results/ctrpv2_aac.pdf", height=7, width=10)
dev.off()

null_aac <- hist(CTRPv2@sensitivity$profiles$auc_recomputed, main="CTRPv2", xlab="AUC", breaks=100, col="gray", ylim=c(0, 10000))
hist(aac[[1]][[1]]$`17596`$`bootstrapped AUCs`, col="red", xlim=c(0, 1), add=T)



d <- density(CTRPv2@sensitivity$profiles$auc_recomputed)
drug_index <- 1
cell_index <- 1
drug <- names(aac)[drug_index]
cell <- names(aac[[drug_index]])[cell_index]

pdf("results/null_drug_cell_rep.pdf", height=7, width=7)
d_rep1 <-  density(aac[[drug_index]][[cell_index]][[1]]$`bootstrapped AUCs`)
d_rep2 <-  density(aac[[drug_index]][[cell_index]][[2]]$`bootstrapped AUCs`)

plot(range(d$x, d_rep1$x, d_rep2$x), range(d$y, d_rep1$y, d_rep2$y), type = "n", xlab = "AAC", ylab = "Density", main=sprintf("CTRPv2, biological  replicates of, %s, %s", drug, cell))
lines(d, col = "black")
lines(d_rep1, col = "red")
lines(d_rep2, col = "blue")
dev.off()
delta_aac_all <- sapply(CTRPv2@sensitivity$profiles$auc_recomputed, function(x){x-CTRPv2@sensitivity$profiles$auc_recomputed})
hist(CTRPv2@sensitivity$profiles$auc_published)
hist(CTRPv2@sensitivity$profiles$auc_recomputed)

 rep_delta_aac <- NULL
for(drug in names(aac)){
  for(cell in names(aac[[drug]])){
    exps <- aac[[drug]][[cell]]
    #print(sprintf("%s %s %s", drug, cell, length(exps)))
    if(length(exps) > 1){
      comb_2 <- combn(length(exps), 2)
      for(ii in 1:ncol(comb_2)){
        xx <- aac[[drug]][[cell]][[comb_2[1, ii]]]
        yy <- aac[[drug]][[cell]][[comb_2[2, ii]]]
        delta.aac <- xx[["bootstrapped AUCs"]] - yy[["bootstrapped AUCs"]]
        rep_delta_aac <- c(rep_delta_aac, delta.aac)
      }
    }
  }
}

pdf("results/bs_rep.pdf", height=5, width=5)
hist(abs(rep_delta_aac), 100, col="gray", main=sprintf("replicates\n#exps=%s * 100 bs", length(rep_delta_aac)/100), xlab="delta aac")
abline(v=quantile(abs(rep_delta_aac), 0.95))
text(x=quantile(abs(rep_delta_aac), 0.95)+.06, y=15000, labels=sprintf("95%%=%.2f", quantile(abs(rep_delta_aac), 0.95)))
dev.off()


exp <- sample(1:length(CTRPv2@sensitivity$profiles$auc_recomputed), sqrt(length(rep_delta_aac) * 2), replace = F)
x <- CTRPv2@sensitivity$profiles$auc_recomputed[exp]
d <- sapply(x, function(b){b-x})
d_both_res_rem <- d
for(b in x){
  d_both_res_rem[which(x == 0 & b == 0)]<- NA
}
length(which(is.na(d_both_res_rem)))


#both_res_cases <- sapply(x, function(b){b-x})
d_upper_tri <- d[upper.tri(d, diag=F)]
hist(d_upper_tri)
pdf("results/diff.pdf", height=5, width=5)
hist(abs(d_upper_tri), 100, col="gray", main=sprintf("A sample of null dist with length ~ 8e5\nOnly %s points with delta aac=0 from fully resistant pairs", length(which(is.na(d_both_res_rem)))), xlab="delta aac", cex.main=.9)
dev.off()

diff_dist <- abs(d_upper_tri)
same_dist <- abs(rep_delta_aac)
hist(diff_dist)
hist(same_dist)
mean(diff_dist)

pdf("results/replicates_vs_different.pdf", height=5, width=5)
hist(diff_dist, col=rgb(1,0,0,0.5),xlim=c(0,1), ylim=c(0,4e5), main=sprintf("The likelihood of pairs mark as different \nonly goes high if delta aac > 0.07"), xlab="delta aac", cex.main=.8)
hist(same_dist, col=rgb(0,0,1,0.5), add=T)
abline(v=mean(diff_dist), col=rgb(1,0,0))
text(x=mean(diff_dist)+.065, y=3e5, labels=sprintf("mu=%.2f\nsd=%.2f", mean(diff_dist), sd(diff_dist)), col=rgb(1,0,0))
abline(v=mean(same_dist), col=rgb(0,0,1))
text(x=mean(same_dist)+.065, y=35e4, labels=sprintf("mu=%.2f\nsd=%.2f", mean(same_dist), sd(same_dist)), col=rgb(0,0,1))
legend("topright", legend=c("different", "replicates"), col=c(rgb(1,0,0), rgb(0,0,1)), pch=18, bty="n")
dev.off()


d_same <- density(same_dist)
d_diff <- density(diff_dist)
v <- 0.07
p_rep <- approx(d_same$x, d_same$y, xout=v)$y * (d_same$x[2] - d_same$x[1])
p_diff <- approx(d_diff$x, d_diff$y, xout=v)$y * (d_diff$x[2] - d_diff$x[1])
ifelse(p_rep < p_diff, "diff", "rep")


df <- approxfun(d_same)
df(v)

-log(approx(d_same$x, d_same$y, xout=v)$y)
-log(approx(d_diff$x, d_diff$y, xout=v)$y)
p_rep <- sum(d_same$y[d_same$x < v]) * (d_same$x[2] - d_same$x[1])
p_diff <- sum(d_diff$y[d_diff$x <= v]) * (d_diff$x[2] - d_diff$x[1])
dnorm(0.1, density(rep_delta_aac), sd(rep_delta_aac))
x=rnorm(N, mean = 3, sd = 2)
LL <- function(mu, sigma) {
       R = dnorm(x, mu, sigma)
       #
     -sum(log(R))
}
mle(LL, start = list(mu = 1, sigma=1))
drug_index <- 1
cell_index <- 1
drug <- names(aac)[drug_index]
cell <- names(aac[[drug_index]])[cell_index]
d_rep1 <-  aac[[drug_index]][[cell_index]][[1]]$`bootstrapped AUCs`
d_rep2 <-  aac[[drug_index]][[cell_index]][[2]]$`bootstrapped AUCs`
d_rep_delta <- d_rep1 - d_rep2
hist(d_rep_delta, col="red", xlim=c(-1,1))


