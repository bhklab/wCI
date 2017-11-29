resampleResiduals <- function(conc, viability, conc_as_log=TRUE, viability_as_pct=FALSE, N=100, nthread=1) {
  availcore <- parallel::detectCores()
  if ( nthread > availcore) {
    nthread <- availcore
  }
  pars <- unlist(logLogisticRegression(conc,
                                       viability,
                                       conc_as_log=conc_as_log,
                                       viability_as_pct=viability_as_pct))
  trueAUC <- computeAUC(conc,
                        viability,
                        conc_as_log=conc_as_log,
                        viability_as_pct=viability_as_pct)
  trueViability <- PharmacoGx:::.Hill(conc, pars)
  errors <- viability - trueViability
  errors <- errors[which(!is.na(errors))]
  library(parallel)
  AUCs <- parallel::mclapply(seq_len(N),
                                    function(x, conc, trueViability, errors, conc_as_log, viability_as_pct){
                                              auc <- computeAUC(conc,
                                                                Hill_fit=unlist(logLogisticRegression(conc,
                                                                                                      trueViability + sample(errors, length(trueViability), replace=TRUE),
                                                                                                      conc_as_log=conc_as_log,
                                                                                                      viability_as_pct=viability_as_pct)),
                                                               conc_as_log=conc_as_log,
                                                                viability_as_pct=viability_as_pct)
                                              return(auc)
    }, conc=conc, trueViability=trueViability, errors=errors, conc_as_log=conc_as_log, viability_as_pct=viability_as_pct, mc.cores=parallel::detectCores())
  AUCs <- do.call(c, AUCs)
  return(list("AUC"=trueAUC, "bootstrapped mean AUC"=mean(AUCs), "bootstrapped SD of AUCs"=sd(AUCs)))
}
