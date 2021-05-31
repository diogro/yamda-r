#' Calculate z-tranform of module correlations using a correlation matrix
#'
#' calculates the z-tranform for a given correlation matrix and modularity hypothesis
#'
#' @param hypot a single modularity hypothesis
#' @param x a correlation matrix
#' @param nneg logical. If TRUE, module coeficients can only be positive
#' @export
#' @importFrom evolqg CreateHypotMatrix
#' @importFrom penalized penalized coef
#' @importFrom utils capture.output
calcZTransCoef = function(hypot, x, nneg){
  z.x = ztrans(x)
  if(ncol(hypot) == 1 && all(hypot == 0)){
    m1 = lm(lt(z.x) ~ 1)
    ztrans_coef = coef(m1)
    names(ztrans_coef) = c("background")
  } else{
    if(is.null(colnames(hypot)))
      colnames(hypot) = paste("module", 1:ncol(hypot), sep = "_")
    n_modules = ncol(hypot)
    mod_pred = t(laply(CreateHypotMatrix(hypot), lt))[,1:n_modules]
    capture.output({
      m1 = penalized(lt(z.x), ~ mod_pred, ~ 1, lambda1 = 0, lambda2 = 0, positive = nneg)
    })
    ztrans_coef = coef(m1,"all")
    names(ztrans_coef) = c("background", paste0("p", colnames(hypot)))
  }
  ztrans_coef
}

calcModuleCorrelations <- function(hypot, ztrans_coef) {
  module_correlations = numeric(length(ztrans_coef))
  module_correlations[1] <- inv_ztrans(ztrans_coef[1])
  if(length(module_correlations) > 1){
    for(k in 2:length(module_correlations)){
      module_correlations[k] <- inv_ztrans(ztrans_coef[1] + 
                                           ztrans_coef[k]) - 
                                inv_ztrans(ztrans_coef[1])
    }
    names(module_correlations) = c("background", colnames(hypot))
  } else{
    names(module_correlations) = c("background")
  }
  module_correlations
}
