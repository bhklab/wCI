resampleResiduals <- function(conc, viability, viability_as_pct = TRUE, conc_as_log = FALSE, N = 100, SDs = 2) {
  
  AUCs <- matrix(NA, nrow = 1, ncol = N)
  pars <- unlist(logLogisticRegression(conc,
                                       viability,
                                       conc_as_log = conc_as_log,
                                       viability_as_pct = viability_as_pct))
  trueAUC <- computeAUC(conc,
                        Hill_fit = pars,
                        conc_as_log = conc_as_log,
                        viability_as_pct = viability_as_pct)
  trueViability <- PharmacoGx:::.Hill(conc, pars)
  errors <- viability - trueViability
  errors <- errors[which(!is.na(errors))]
  for (i in seq_len(N)) {
    pars <- 
      AUCs[i] <- computeAUC(conc,
                            Hill_fit = unlist(logLogisticRegression(conc,
                                                         trueViability + sample(errors, length(trueViability), replace = TRUE),
                                                         conc_as_log = conc_as_log,
                                                         viability_as_pct = viability_as_pct)),
                            conc_as_log = conc_as_log,
                            viability_as_pct = viability_as_pct)
  }
  
  return(list("AUC" = trueAUC, "lower" = trueAUC - SDs * sd(AUCs), "upper" = trueAUC + SDs * sd(AUCs)))
}