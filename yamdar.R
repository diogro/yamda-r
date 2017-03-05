devtools::install_github("diogro/evolqg")
install.packages("arm")
install.packages("superheat")
install.packages("cowplot")
install.packages("penalized")

library(evolqg)
library(arm)
library(superheat)
library(cowplot)
library(penalized)

lt = function(x) x[lower.tri(x)]
ztrans = function(x) 0.5 * log((1 + x) / (1 - x))
inv_ztrans = function(x) (exp(2*x) - 1) / (exp(2*x) + 1)
Yamda = function(mod.cor, hypot_list, n, nneg = FALSE){
  z.mod.cor = ztrans(mod.cor)
  n_models = length(hypot_list)
  if(is.null(names(hypot_list)))
    names(hypot_list) = paste("Hypothesis", 1:n_models, sep = "_")
  stats = data.frame(hypothesis = character(), LL = numeric(), param = numeric(), AICc = numeric())
  expected_matrices = vector("list", n_models+1)
  module_correlations = vector("list", n_models+1)
  # Calculating for each actual hypothesis
  for(i in seq_along(hypot_list)){
    current_hypot = as.matrix(hypot_list[[i]])
    n_modules = ncol(current_hypot)
    if(is.null(colnames(current_hypot)))
      colnames(current_hypot) = paste("module", 1:n_modules, sep = "_")
    mod_pred = t(laply(CreateHypotMatrix(current_hypot), function(x) x[lower.tri(x)]))[,1:n_modules]
    m1 = penalized(lt(z.mod.cor), ~ mod_pred, ~1, lambda1 = 0, lambda2 = 0, positive=nneg)
    module_correlations[[i]] = inv_ztrans(coef(m1,"all"))
    names(module_correlations[[i]]) = c("background", colnames(current_hypot))
    expected = Reduce("+", Map("*", coef(m1,"all")[2:(n_modules+1)],
                               CreateHypotMatrix(current_hypot)[1:n_modules])) + coef(m1,"all")[1]
    expected = inv_ztrans(expected)
        for(k in 2:length(module_correlations[[i]])){
      module_correlations[[i]][k] <- inv_ztrans(coef(m1,"all")[1]+coef(m1,"all")[k])-inv_ztrans(coef(m1,"all")[1])
    }
    diag(expected) = 1
    expected_matrices[[i]] = expected
    LogL = function(z_r, z_p) {-0.5 * log(var) - ((z_r - z_p)^2) / (2 * var)}
    var = 1/(n - 3)
    LL = sum(LogL(lt(expected), lt(z.mod.cor)))
    n_corr = length(lt(z.mod.cor))
    param = length(coef(m1,"all"))
    AICc = - 2 * LL + 2 * param + (2 * param * (param + 1)) / (n_corr - param - 1)
    stats = rbind(stats, data.frame(Hypothesis = names(hypot_list)[[i]], LL, param, AICc))
  }
  m1<-lm(lt(z.mod.cor)~1)
  module_correlations[[length(module_correlations)]] <- inv_ztrans(coef(m1))
  names(module_correlations[[length(module_correlations)]]) <- "background"
  expected = matrix(coef(m1)[1],nrow=dim(z.mod.cor)[1],ncol=dim(z.mod.cor)[2])
  expected = inv_ztrans(expected)
  diag(expected) = 1
  expected_matrices[[length(expected_matrices)]]<-expected
  LogL = function(z_r, z_p) {-0.5 * log(var) - ((z_r - z_p)^2) / (2 * var)}
  var = 1/(n - 3)
  LL = sum(LogL(lt(expected), lt(z.mod.cor)))
  n_corr = length(lt(z.mod.cor))
  param = length(coef(m1))
  AICc = - 2 * LL + 2 * param + (2 * param * (param + 1)) / (n_corr - param - 1)
  stats = rbind(stats, data.frame(Hypothesis = "No Modularity", LL, param, AICc))
  
  names(expected_matrices) = c(names(hypot_list),"No Modularity")
  names(module_correlations) = c(names(hypot_list),"No Modularity")
  stats$dAICc <- stats$AIC - min(stats$AICc)
  stats = stats[order(stats$dAICc),]
  return(list(stats = stats, 
              module_correlations = module_correlations[order(stats$AICc)], 
              expected_matrices = expected_matrices[order(stats$AICc)]))
}
                           
#####
## Exemplo com uma matriz conhecida
#####

true_hypots = matrix( c(
           0, 1, 0, 1, 0,
           0, 1, 0, 1, 0,
           1, 1, 0, 0, 1,
           0, 1, 1, 0, 1,
           1, 0, 1, 0, 0,
           0, 0, 1, 0, 0,
           0, 0, 1, 0, 0,
           0, 0, 1, 1, 0,
           0, 0, 1, 1, 0,
           0, 0, 1, 1, 0), ncol = 5, byrow = T)
hypots = matrix(
                c(
           0, 1, 0,
           0, 1, 0,
           0, 0, 1,
           1, 0, 1,
           1, 0, 0,
           1, 0, 0,
           1, 0, 0,
           1, 0, 0,
           1, 0, 0,
           1, 0, 0), ncol = 3, byrow = T)
# Gerando a matriz a partir de uma da primeira hipotese
corrs = c(0.1, 0.3, 0.5, 0.6, 0.4)
mod.cor = matrix(ztrans(0.2), 10, 10)
for(i in 1:4){
    cor.hypot = CreateHypotMatrix(true_hypots)[[i]]
    hypot.mask = matrix(as.logical(cor.hypot), 10, 10)
    mod.cor[ hypot.mask] = mod.cor[ hypot.mask] + rnorm(length(mod.cor[hypot.mask]),
                                                        ztrans(corrs[i]), 0.1) # within-modules
}
mod.cor = (mod.cor + t(mod.cor)) / 2
mod.cor = inv_ztrans(mod.cor)
mod.cor = mod.cor + 0.05 * RandomMatrix(10)
diag(mod.cor) = 1

                           
# Lista com todas as hipoteses a serem comparadas, no esquema de sempre. Os modulos podem ser sobrepostos                           
hypot_list = list(true_hypots, hypots)
                           
# Função em si. Ela retorna uma lista com a comparação de modelos, as estimativas dos coeficientes pra cada modulo e a matriz "esperada"                           
x = Yamda(mod.cor, hypot_list, 10)

# Plot da matriz original com a esperada. Nesse caso é ótimo...
plot_grid(superheat(mod.cor)[[2]], superheat(x[[3]][[1]])[[2]] )
x$stats
