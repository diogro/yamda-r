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
 #' @importFrom penalized penalized coef
 #' @examples
 #' data(toadCor)
 #' data(toadHypo)
 #'
 #' Yamda(toadCor, toadHypo, 25, nneg = FALSE)
Yamda = function(x, hypot_list, n, nneg = FALSE){
  z.x = ztrans(x)
  n_models = length(hypot_list)
  if(is.null(names(hypot_list)))
    names(hypot_list) = paste("Hypothesis", 1:n_models, sep = "_")
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

  m1 <- lm(lt(z.x)~1)
  module_correlations[[length(module_correlations)]] <- inv_ztrans(coef(m1))
  names(module_correlations[[length(module_correlations)]]) <- "background"
  ztrans_coef[[length(ztrans_coef)]] <- coef(m1)
  names(ztrans_coef[[length(ztrans_coef)]]) <- "background"
  expected = matrix(coef(m1)[1],nrow=dim(z.x)[1],ncol=dim(z.x)[2])
  expected = inv_ztrans(expected)
  diag(expected) = 1
  expected_matrices[[length(expected_matrices)]] <- expected

  LL = sum(LogL(lt(expected), lt(z.x), var = 1/(n - 3)))
  n_corr = length(lt(z.x))
  param = length(coef(m1))
  AICc = - 2 * LL + 2 * param + (2 * param * (param + 1)) / (n_corr - param - 1)
  stats = rbind(stats, data.frame(Hypothesis = "No Modularity", LL, param, AICc))

  names(expected_matrices) = c(names(hypot_list),"No Modularity")
  names(module_correlations) = c(names(hypot_list),"No Modularity")
  names(ztrans_coef) = c(names(hypot_list),"No Modularity")
  stats$dAICc <- stats$AIC - min(stats$AICc)
  stats$ModelLogL <- exp(-0.5 * stats$dAICc)
  stats$AkaikeWeight <- stats$ModelLogL/sum(stats$ModelLogL)
  stats = stats[order(stats$dAICc),]
  return(list(stats = stats,
              ztrans_coef = ztrans_coef[order(stats$AICc)],
              module_correlations = module_correlations[order(stats$AICc)],
              expected_matrices = expected_matrices[order(stats$AICc)]))
}

calcZTransCoef = function(current_hypot, x, nneg){
  z.x = ztrans(x)
  n_modules = ncol(current_hypot)
  mod_pred = t(laply(CreateHypotMatrix(current_hypot), lt))[,1:n_modules]
  m1 = penalized(lt(z.x), ~ mod_pred, ~ 1, lambda1 = 0, lambda2 = 0, positive = nneg)
  ztrans_coef = coef(m1,"all")
  names(ztrans_coef) = c("background", colnames(current_hypot))
  ztrans_coef
}

calcExpectedMatrix <- function(current_hypot, ztrans_coef) {
  n_modules = ncol(current_hypot)
  expected = Reduce("+", Map("*", ztrans_coef[2:(n_modules+1)],
                             CreateHypotMatrix(current_hypot)[1:n_modules])) + ztrans_coef[1]
  expected = inv_ztrans(expected)
  diag(expected) = 1
  expected
}

calcModuleCorrelations <- function(current_hypot, ztrans_coef) {
  module_correlations = numeric(length(ztrans_coef))
  module_correlations[1] <- inv_ztrans(ztrans_coef[1])
  for(k in 2:length(module_correlations)){
    module_correlations[k] <- inv_ztrans(ztrans_coef[1] + ztrans_coef[k]) - inv_ztrans(ztrans_coef[1])
  }
  names(module_correlations) = c("background", colnames(current_hypot))
  module_correlations
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
