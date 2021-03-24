library(cowplot)
library(yamdar)
library(superheat)
library(evolqg)
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

data("toadCor")
data("toadHypo")
full = hclustHypot(toadCor)
class(toadHypo[[1]])
toadHypo[[5]] = full
names(toadHypo)[5] = "hclust"
x = YamdaFactorsMLE(toadCor, full, T)
x[[1]]
x = Yamda(toadCor, toadHypo, 25, T)
all = do.call(cbind, toadHypo)
x = Yamda(toadCor, list(all), 25, T)

toadHypo[[5]] = evolqg::LModularity(toadCor)[[2]]

plot_expected(toadCor, x)
x11()
superheat(x$expected_matrices[[14]], row.dendrogram = T, col.dendrogram = T)
superheat(toadCor, row.dendrogram = T, col.dendrogram = T)

y = Yamda(toadCor, list(full[[3]][,-5]), 40, T)
toMembership(full[[5]])
full[[5]]

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
library(mvtnorm)
toadCor
pop = rmvnorm(100, sigma = mod.cor)
modules_2 = matrix(c(rep(c(1, 0, 0), each = 5),
                   rep(c(0, 1, 0), each = 5)),    15)
modules_1 = matrix(c(rep(c(1, 0, 0), each = 5)),    15)
modules_15 = diag(15)
modules_15[1,2] = 1
hypot = list(modules, modules_1, modules_2, modules_15[,2])
result = Yamda(mod.cor, hypot, n = 50)

n_corr = (15 * 15 - 15)/2
LL  = c(sapply(result[[4]], function(x) sum(dmvnorm(pop, sigma = x, log = T))),
        sum(dmvnorm(pop, sigma = mod.cor, log = T)))
param = c(result[[1]]$param, (15 * 15 - 15)/2 + 15)
-2 * LL[[1]] + 2 * param + (2 * param * (param + 1)) / (100 - param - 1)


