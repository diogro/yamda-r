#' Yet another modularity detection algorith in R
#'
#' Uses hierarchical partitions and MLE to rank modularity hypothesis using individual measurements.
#'
#' @export
#' @param data individual measurements or residuals
#' @param hypot_list list of matrices describing modularity hypothesis. Each element in the list should be one matrix, and each column in these matrices represents a module, and each row a trait. If position [i,j] is one, the trait i belong to module j, and if it is zero, trait i does not belong to module j. Modules can be overlapping.
#' @param nneg If true, belonging to the same module can only increase correlations, not decrease them.
#' @importFrom plyr laply
#' @importFrom penalized penalized coef
#' @importFrom bbmle logLik AICc
#' @examples
#' library(mvtnorm)
#' modules = matrix(c(rep(c(1, 0, 0), each = 5),
#'                    rep(c(0, 1, 0), each = 5),
#'                    rep(c(0, 0, 1), each = 5)), 15)
#' modules_ztrans = c(0.3, 0.5, .3, 1.7)
#' mod.cor = calcExpectedMatrix(modules, modules_ztrans)
#' sds = runif(15, 1, 10)
#' mod_cov = outer(sds, sds) * mod.cor
#' pop = rmvnorm(50, sigma = mod_cov)
#' YamdaMLE(pop, list(modules, modules[,-1]), TRUE)[[2]]
YamdaMLE = function(data, hypot_list, nneg = TRUE){
  n_models = length(hypot_list)
  if(is.null(names(hypot_list)))
    names(hypot_list) = paste("Hypothesis", 1:n_models, sep = "_")
  hypot_list[[n_models+1]] = matrix(0, ncol(data), 1)
  names(hypot_list)[[n_models+1]] = "No Modularity"
  stats = data.frame(hypothesis = character(), LL = numeric(), param = numeric(), AICc = numeric())
  models = vector("list", n_models+1)
  expected_matrices = vector("list", n_models+1)
  module_correlations = vector("list", n_models+1)
  ztrans_coef = vector("list", n_models+1)
  # Calculating for each actual hypothesis
  for(i in seq_along(hypot_list)){
    current_hypot = as.matrix(hypot_list[[i]])
    n_modules = ncol(current_hypot)
    if(is.null(colnames(current_hypot)))
      colnames(current_hypot) = paste("module", 1:n_modules, sep = "_")
    models[[i]] = fitZTransCoefML(data, current_hypot, nneg)
    if(nneg)
      ztrans_coef[[i]] = c(coef(models[[i]])[1], exp(coef(models[[i]])[-1]))
    else
      ztrans_coef[[i]] = coef(models[[i]])
    expected_matrices[[i]] = calcExpectedMatrix(current_hypot, ztrans_coef[[i]])
    module_correlations[[i]] = calcModuleCorrelations(current_hypot, ztrans_coef[[i]])
    stats = rbind(stats, calcModelStatsMLE(models[[i]], names(hypot_list)[[i]]))
  }
  names(expected_matrices) = names(hypot_list)
  names(module_correlations) = names(hypot_list)
  names(ztrans_coef) = names(hypot_list)
  names(models) = names(hypot_list)
  stats$dAICc <- stats$AIC - min(stats$AICc)
  stats$ModelLogL <- exp(-0.5 * stats$dAICc)
  stats$AkaikeWeight <- stats$ModelLogL/sum(stats$ModelLogL)
  stats = stats[order(stats$dAICc),]
  return(list(stats = stats,
              ztrans_coef = ztrans_coef[order(stats$AICc)],
              module_correlations = module_correlations[order(stats$AICc)],
              expected_matrices = expected_matrices[order(stats$AICc)],
              models = models[order(stats$AICc)]))
}
