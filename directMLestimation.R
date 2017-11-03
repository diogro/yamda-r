library(evolqg)
library(corrgram)
library(mvtnorm)
library(penalized)

modules = matrix(c(rep(c(1, 0, 0), each = 5),
                   rep(c(0, 1, 0), each = 5),
                   rep(c(0, 0, 1), each = 5)), 15)
cor.hypot = CreateHypotMatrix(modules)[[4]]
hypot.mask = matrix(as.logical(cor.hypot), 15, 15)
mod.cor = matrix(NA, 15, 15)
mod.cor[ hypot.mask] = runif(length(mod.cor[ hypot.mask]), 0.8, 0.9) # within-modules
mod.cor[!hypot.mask] = runif(length(mod.cor[!hypot.mask]), 0.3, 0.4) # between-modules
diag(mod.cor) = 1
mod.cor = (mod.cor + t(mod.cor))/2 # correlation matrices should be symmetric
sds = runif(15, 1, 10)
mod_cov = outer(sds, sds) * mod.cor

pop = rmvnorm(50, sigma = mod_cov)

random_hypot = matrix(sample(c(0, 1), 3 * 15, replace = TRUE), 15, 3)

x = pop
hypot = modules
pars = params
modular_logL = function(x, hypot, pars){
  x_sds = apply(x, 2, sd)
  sum(dmvnorm(x, colMeans(x), sigma = outer(x_sds, pop_sds) * calcExpectedMatrix(hypot, pars), log = TRUE))
}

modular_logL(pop, modules, params)
params = calcZTransCoef(modules, mod.cor, F)
params = calcZTransCoef(hypot, mod.cor, F)
params = numeric(ncol(hypot) + 1)

