paramBootStrap <- function(conc, viability, conc_as_log=TRUE, viability_as_pct=FALSE, N=100, nthread=1, trunc=TRUE, area.type="Fitted", degree_of_freedom=1, verbose=FALSE) {
  availcore <- parallel::detectCores()
  if ( nthread > availcore) {
    nthread <- availcore
  }
  cix <- complete.cases(conc, viability)
  conc <- conc[cix]
  viability <- viability[cix]
  pars <- unlist(logLogisticRegression(conc,
                                       viability,
                                       conc_as_log=conc_as_log,
                                       viability_as_pct=viability_as_pct,
                                       trunc=trunc,
                                       verbose=verbose))
  trueAUC <- computeAUC(conc,
                        viability,
                        conc_as_log=conc_as_log,
                        viability_as_pct=viability_as_pct,
                        trunc=trunc,
                        area.type=area.type,
                        verbose=verbose)
  trueViability <- PharmacoGx:::.Hill(conc, pars)
  errors <- viability - trueViability
  errors <- errors[which(!is.na(errors))]
  #ee <- rnorm(1000, mean=mean(errors), sd=sd(errors))
  #hist(ee)
  library(parallel)
  AUCs <- parallel::mclapply(seq_len(N),
                             function(x, conc, trueViability, errors, conc_as_log, viability_as_pct, trunc, verbose){
                               if(degree_of_freedom == 1){
                                 sd_errors <- sd(errors)
                               }else{
                                 sd_errors <- sqrt(sum((errors - mean(errors)) ^ 2)/(length(errors) - degree_of_freedom))
                               }
                               if(area.type=="Fitted"){
                                 auc <- computeAUC(conc,
                                                   Hill_fit=unlist(logLogisticRegression(conc,
                                                                                         #trueViability + sample(errors, length(trueViability), replace=TRUE),
                                                                                         trueViability + rnorm(length(trueViability), mean(errors), sd_errors),
                                                                                         conc_as_log=conc_as_log,
                                                                                         viability_as_pct=viability_as_pct,
                                                                                         trunc=trunc,
                                                                                         verbose=verbose)),
                                                   conc_as_log=conc_as_log,
                                                   viability_as_pct=viability_as_pct, trunc=trunc, verbose=verbose)
                               }else{
                                 auc <- computeAUC(conc,
                                                   viability=trueViability + rnorm(length(trueViability), mean(errors), sd_errors),
                                                   conc_as_log=conc_as_log,
                                                   viability_as_pct=viability_as_pct,
                                                   trunc=trunc,
                                                   area.type=area.type,
                                                   verbose=verbose)
                               }
                               return(auc)
                             }, conc=conc, trueViability=trueViability, errors=errors, conc_as_log=conc_as_log, viability_as_pct=viability_as_pct, trunc=trunc, verbose=verbose, mc.cores=parallel::detectCores())
  AUCs <- do.call(c, AUCs)
  return(list("AUC"=trueAUC, "bootstrapped AUCs"=AUCs, "errors"=errors))
}
