#' Yet another modularity detection algorith in R
#'
#' Uses hierarchical partitions and linear models to rank modularity hypothesis using correlation matrices.
#'
#' @export
#' @param x correlation matrix
#' @param hypot_list list of matrices describing modularity hypothesis. Each element in the list should be one matrix, and each column in these matrices represents a module, and each row a trait. If position [i,j] is one, the trait i belong to module j, and if it is zero, trait i does not belong to module j. Modules can be overlapping.
#' @param n number of individuals used in the correlation matrix estimation
#' @param nneg If true, belonging to the same module can only increase correlations, not decrease them.
#' @importFrom evolqg CreateHypotMatrix
#' @importFrom plyr laply
#' @examples
#' data(toadCor)
#' data(toadHypo)
#' Yamda(toadCor, toadHypo, 25, nneg = FALSE)
YamdaLM = function(x, hypot_list, n, nneg = FALSE){
  n_models = length(hypot_list)
  if(is.null(names(hypot_list)))
    names(hypot_list) = paste("Hypothesis", 1:n_models, sep = "_")
  hypot_list[[n_models+1]] = matrix(0, ncol(x), 1)
  names(hypot_list)[[n_models+1]] = "No Modularity"
  stats = data.frame(hypothesis = character(), LL = numeric(), param = numeric(), AICc = numeric())
  expected_matrices = vector("list", n_models+1)
  module_correlations = vector("list", n_models+1)
  ztrans_coef = vector("list", n_models+1)
  # Calculating for each actual hypothesis
  for(i in seq_along(hypot_list)){
    current_hypot = as.matrix(hypot_list[[i]])
    n_modules = ncol(current_hypot)
    if(is.null(colnames(current_hypot)))
      colnames(current_hypot) = paste("module", 1:n_modules, sep = "_")
    ztrans_coef[[i]] = calcZTransCoef(current_hypot, x, nneg)
    expected_matrices[[i]] = calcExpectedMatrix(current_hypot, ztrans_coef[[i]])
    module_correlations[[i]] = calcModuleCorrelations(current_hypot, ztrans_coef[[i]])
    stats = rbind(stats, calcModelStats(x, expected_matrices[[i]], n, ztrans_coef[[i]], names(hypot_list)[[i]]))
  }
  names(expected_matrices) = names(hypot_list)
  names(module_correlations) = names(hypot_list)
  names(ztrans_coef) = names(hypot_list)
  stats$dAICc <- stats$AIC - min(stats$AICc)
  stats$ModelLogL <- exp(-0.5 * stats$dAICc)
  stats$AkaikeWeight <- stats$ModelLogL/sum(stats$ModelLogL)
  stats = stats[order(stats$dAICc),]
  return(list(stats = stats,
              ztrans_coef = ztrans_coef[order(stats$AICc)],
              module_correlations = module_correlations[order(stats$AICc)],
              expected_matrices = expected_matrices[order(stats$AICc)]))
}

LogL = function(z_r, z_p, var) {-0.5 * log(var) - ((z_r - z_p)^2) / (2 * var)}

calcModelStats <- function(x, expected_matrix, n, ztrans_coef, current_hypot_name) {
  z.x = ztrans(x)
  LL = sum(LogL(lt(expected_matrix), lt(z.x), var = 1/(n - 3)))
  n_corr = length(lt(z.x))
  param = length(ztrans_coef)
  AICc = - 2 * LL + 2 * param + (2 * param * (param + 1)) / (n_corr - param - 1)
  data.frame(Hypothesis = current_hypot_name, LL, param, AICc)
}
