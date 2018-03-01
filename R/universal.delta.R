if(!file.exists("data/PSets")){
  dir.create("data/PSets", recursive=T)
}
###delta auc for replicates in GDSC
if(!file.exists("data/PSets/CTRPv2.RData")){
  downloadPSet("GDSC", saveDir="data/PSets/")
}
load("data/PSets/GDSC.RData")
ids <- unlist(strsplit(GDSC@drug[which(rownames(GDSC@drug) == "AZD6482"), "drugid"], split="/"))

tt <- GDSC@sensitivity$profiles[grep(paste0("drugid_", ids[1]), rownames(GDSC@sensitivity$info)), "auc_recomputed"]
ss <- GDSC@sensitivity$profiles[grep(paste0("drugid_", ids[2]), rownames(GDSC@sensitivity$info)),"auc_recomputed"]

names(tt) <- GDSC@sensitivity$info[grep(paste0("drugid_", ids[1]), rownames(GDSC@sensitivity$info)), "cellid"]
names(ss) <- GDSC@sensitivity$info[grep(paste0("drugid_", ids[2]), rownames(GDSC@sensitivity$info)), "cellid"]

nn <- intersect(names(tt), names(ss))
gdsc.delta.auc <- abs(ss[nn]-tt[nn])
qq.gdsc <- quantile(gdsc.delta.auc, probs = seq(0, 1, 0.05))[c("90%","95%")]

pdf("AZD6482_auc_delta.pdf", height=5, width=10)
par(mfrow=c(1,2))
plot(gdsc.delta.auc, pch=20, xlab="cells", ylab="delta auc", main="AZD6482")
hist(gdsc.delta.auc, col="gray", breaks=100, xlab="delta auc", main="AZD6482")
legend("topright", legend=paste(c("90%","95%"), round(qq.gdsc, digits=2), sep=":"), bty="n")
dev.off()

###delta auc for replicates in CTRPv2
if(!file.exists("data/PSets/CTRPv2.RData")){
  downloadPSet("CTRPv2", saveDir="data/PSets/")
}
load("data/PSets/CTRPv2.RData")
dd <- apply(CTRPv2@sensitivity$n, 2, function(x){length(which(x > 1))})
dd.rep <- names(which(dd > 15))
qq.ctrp <- NULL
ctrp.delta.auc <- NULL
for(drug in dd.rep){
  cells <- names(which(CTRPv2@sensitivity$n[,drug] > 1))
  ii <- sensitivityInfo(CTRPv2)[which(sensitivityInfo(CTRPv2)$drugid==drug & sensitivityInfo(CTRPv2)$cellid %in% cells),]
  delta.auc <- NULL
  for(cell in cells){
    xx <- sensitivityProfiles(CTRPv2)[rownames(ii)[which(ii$cellid == cell)], "auc_recomputed", drop=T]
    tmp.delta <- NULL
    for(i in 1:(length(xx)-1)){
      for(j in (i+1):length(xx)){
        tmp.delta <- c(abs(xx[i]-xx[j]), tmp.delta)
      }
    }
    delta.auc <- c(delta.auc, mean(tmp.delta))                                
  }
  qq.ctrp <- rbind(qq.ctrp, quantile(delta.auc, probs = seq(0, 1, 0.05))[c("90%","95%")])
  ctrp.delta.auc <- c(ctrp.delta.auc, delta.auc)
}
qq.ctrp.avg <- apply(qq.ctrp, 2, mean)
qq.ctrp.max <- apply(qq.ctrp, 2, max)
ctrp.drug.max.delta.auc <- apply(qq.ctrp, 2, which.max)
ctrp.quantile.over.combined.drugs.delta.auc <- quantile(ctrp.delta.auc, probs = seq(0, 1, 0.05))[c("90%","95%")]

pdf("ctrp_auc_delta.pdf", height=5, width=5)
hist(qq.ctrp, col="gray", breaks=100, xlab="delta auc", main="CTRP")
legend("topright", legend=paste(c("90%","95%"), round(ctrp.quantile.over.combined.drugs.delta.auc, digits=2), sep=":"), bty="n")
dev.off()

###delta auc for replicates in GRAY
if(!file.exists("data/PSets/CTRPv2.RData")){
  downloadPSet("GRAY", saveDir="data/PSets/")
}
load("data/PSets/GRAY.RData")
dd <- apply(GRAY@sensitivity$n, 2, function(x){length(which(x > 1))})
dd.rep <- names(which(dd > 10))
qq.gray <- NULL
gray.delta.auc <- NULL
for(drug in dd.rep){
  cells <- names(which(GRAY@sensitivity$n[,drug] > 1))
  ii <- sensitivityInfo(GRAY)[which(sensitivityInfo(GRAY)$drugid==drug & sensitivityInfo(GRAY)$cellid %in% cells),]
  delta.auc <- NULL
  for(cell in cells){
    xx <- sensitivityProfiles(GRAY)[rownames(ii)[which(ii$cellid == cell)], "auc_recomputed", drop=T]
    tmp.delta <- NULL
    for(i in 1:(length(xx)-1)){
      for(j in (i+1):length(xx)){
        tmp.delta <- c(abs( xx[i]-xx[j]), tmp.delta)
      }
    }
    delta.auc <- c(delta.auc, mean(tmp.delta))                                
  }
  qq.gray <- rbind(qq, quantile(delta.auc, probs = seq(0, 1, 0.05))[c("90%","95%")])
  gray.delta.auc <- c(gray.delta.auc, delta.auc)
}
qq.gray.avg <- apply(qq, 2, mean)
qq.gray.max <- apply(qq, 2, max)
gray.drug.max.delta.auc <- apply(qq, 2, which.max)
gray.quantile.over.combined.drugs.delta.auc <- quantile(gray.delta.auc, probs = seq(0, 1, 0.05))[c("90%","95%")]

pdf("gray_auc_delta.pdf", height=5, width=5)
hist(qq.ctrp, col="gray", breaks=100, xlab="delta auc", main="GRAY")
legend("topright", legend=paste(c("90%","95%"), round(gray.quantile.over.combined.drugs.delta.auc, digits=2), sep=":"), bty="n")
dev.off()

################
univerasal.delta1 <- quantile(c(ctrp.delta.auc, gray.delta.auc, gdsc.delta.auc), probs = seq(0, 1, 0.05))[c("90%","95%")]
pdf("all_datasets_auc_delta.pdf", height=5, width=5)
hist(c(ctrp.delta.auc, gray.delta.auc, gdsc.delta.auc), col="gray", breaks=100, xlab="delta auc", main="CTRP, GDSC, GRAY")
legend("topright", legend=paste(c("90%","95%"), round(univerasal.delta1, digits=2), sep=":"), bty="n")
dev.off()
###############

####for test
univerasal.delta2 <- apply(rbind(qq.gdsc, qq.ctrp.avg, qq.gray.avg), 2 , mean)
univerasal.delta3 <- apply(rbind(qq.gdsc, ctrp.quantile.over.combined.drugs.delta.auc, gray.quantile.over.combined.drugs.delta.auc), 2 , mean)

univerasal.delta4 <- apply(rbind(univerasal.delta1, univerasal.delta2, univerasal.delta3), 2, mean)

