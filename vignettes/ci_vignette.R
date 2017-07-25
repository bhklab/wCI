## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----ci1-----------------------------------------------------------------
library(CI)
library(MASS)

## ----ci2-----------------------------------------------------------------
cr = 0.95
df <- mvrnorm(50, mu = c(0,0), Sigma = matrix(c(1,cr,cr,1), ncol = 2), empirical = TRUE)

## ----ci3-----------------------------------------------------------------
CI(x=df[,1], y=df[,2], deltaX=0, deltaY=0, alpha =0, outx = 1, npermut=10000)

## ----ci4-----------------------------------------------------------------
cor.test(df[,1], df[,2])
cor.test(df[,1], df[,2], method = "pearson")
cor.test(df[,1], df[,2], method = "spearman")

## ----plotFun-------------------------------------------------------------
plotData <- function(x,y,title="", xlab="x", ylab="y")
{
  par(pty="s")
  plot(x,y, main=title, xlab=xlab, ylab=ylab, pch=16, col="red")
}
  

## ----ci5-----------------------------------------------------------------

dt <- lapply(seq(-1,1,0.1), function(cr){
      df <- mvrnorm(50, mu = c(0,0), Sigma = matrix(c(1,cr,cr,1), ncol = 2), empirical = TRUE)
      ci <- CI(x=df[,1], y=df[,2], deltaX=0, deltaY=0, alpha =0, outx = 1, npermut=10000, ncpu = 1)
      pr <- cor.test(df[,1], df[,2], method = "pearson")
      sm <- cor.test(df[,1], df[,2], method = "spearman")
      data.frame(ci=ci$ci, ci.p=ci$p.value,
                 pr=pr$estimate, pr.p=pr$p.value,
                 sm=sm$estimate, sm.p=sm$p.value)
      })

dt <- do.call(rbind.data.frame, dt)

plotData(dt$ci, dt$pr,title="CI vs pearson correlation", xlab="CI", ylab="pearson")
plotData(dt$ci, dt$sm,title="CI vs spearman correlation", xlab="CI", ylab="spearman")

pOfset = 1
plotData(log(dt$ci.p+pOfset), log(dt$pr.p+pOfset), title="p-value (CI vs pearson correlation)",
         xlab="CI p-value", ylab="pearson p-value")

plotData(log(dt$ci.p+pOfset), log(dt$sm.p+pOfset), title="p-value (CI vs spearman correlation)",
         xlab="CI p-value", ylab="spearman p-value")


