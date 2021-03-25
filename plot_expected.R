library(cowplot)
library(yamdar)
library(superheat)
library(evolqg)
library(mvtnorm)

plot_expected = function(original_matrix, yamda_object, n = 3, 
                         nrow = NULL, ncol = NULL, ...){
  expected_matrices = yamda_object$expected_matrices
  n_hypot = length(expected_matrices) - 1
  expected_matrices = expected_matrices[yamda_object$stats$Hypothesis]
  plots = vector("list", n + 1)
  plots[[1]] = superheat(original_matrix, print.plot = FALSE, title = "Observed", ...)[[2]]
  for(i in 1:n){
    plots[[i+1]] = superheat(expected_matrices[[i]], print.plot = FALSE, 
                             title =  names(expected_matrices)[i] ,...)[[2]]
  }
  plot_grid(plotlist = plots, nrow = nrow, ncol = ncol)
}

modules = matrix(c(rep(c(1, 0, 0), each = 5),
                   rep(c(0, 1, 0), each = 5),
                   rep(c(0, 0, 1), each = 5)), 15)
cor.hypot = CreateHypotMatrix(modules)[[4]]
hypot.mask = matrix(as.logical(cor.hypot), 15, 15)
mod.cor = matrix(NA, 15, 15)
mod.cor[ hypot.mask] = runif(length(mod.cor[ hypot.mask]), 0.8, 0.9) # within-modules
mod.cor[!hypot.mask] = runif(length(mod.cor[!hypot.mask]), 0.3, 0.4) # between-modules
diag(mod.cor) = 1
mod.cor = (mod.cor + t(mod.cor))/2 # correlation matrices s

pop = rmvnorm(100, sigma = mod.cor)
modules_2 = matrix(c(rep(c(1, 0, 0), each = 5),
                     rep(c(0, 1, 0), each = 5)),    15)
modules_1 = matrix(c(rep(c(1, 0, 0), each = 5)),    15)
hypot = list(modules, modules_1, modules_2)
result = Yamda(pop, hypot, factors = FALSE)

plot_expected(original_matrix = mod.cor, yamda_object = result, n = 3)
