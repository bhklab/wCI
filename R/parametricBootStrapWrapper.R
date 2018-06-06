parametricBootStrapWrapper <- function(pSet, cell, drug, experiment_id, N=100, nthread=1, trunc=TRUE, area.type="Fitted", conc_as_log=TRUE, viability_as_pct=FALSE, degree_of_freedom=1, verbose=FALSE){
  library(PharmacoGx)
  if(missing(experiment_id)){
    if(missing(cell) | missing(drug)){
      return("An experiment id or name of the cell and drug should be passed to the function!")
    }else{
      if(typeof(cell) != "character" | length(cell) != 1 | typeof(drug) != "character" | length(drug) != 1 | class(pSet) != "PharmacoSet" | length(pSet) != 1){
        return("Cell, drug and pSet are expected to be the name of one cell, one drug and one pSet!")
      }
      experiment_id <- which(PharmacoGx::sensitivityInfo(pSet)[,"cellid"] == cell & PharmacoGx::sensitivityInfo(pSet)[,"drugid"] == drug)[1]
    }
  }
  return(paramBootStrap(conc=log10(as.numeric(pSet@sensitivity$raw[experiment_id, , "Dose"])),
                        viability=as.numeric(pSet@sensitivity$raw[experiment_id, , "Viability"])/100,
                        N=N,
                        nthread=nthread,
                        trunc=trunc,
                        area.type=area.type,
                        conc_as_log=conc_as_log,
                        viability_as_pct=viability_as_pct,
                        degree_of_freedom=degree_of_freedom,
                        verbose=verbose))
}
