library(yamdar)
library(mvtnorm)

true_factors = matrix(c(rep(c(1, 0), c(10, 5)),
                        rep(c(0, 1), c(5, 10))), 15)
modules_coef = c(0.3, 0.5, 0.5)
mod.cor = calcExpectedMatrixFactors(true_factors, modules_coef)
sds = runif(15, 1, 3)
mod_cov = outer(sds, sds) * mod.cor

modules = matrix(c(rep(c(1, 0), c(10, 5)),
                   rep(c(0, 1), c(5, 10)),
                   rep(c(1, 0, 0), each = 5),
                   rep(c(0, 0, 1), each = 5),
                   sample(c(0, 1), 15, replace = T)), 15)
pop = rmvnorm(100, sigma = mod_cov)
hypot = list(tudo = modules, sem_zuera = modules[,-5], true = modules[,-c(3, 4, 5)])
YamdaFactorsMLE(pop, hypot, TRUE)[1:2]
YamdaFactorsMLE(pop, hypot, FALSE)[1:2]
y_mle = YamdaMLE(pop, hypot, FALSE)

true_factors = matrix(c(rep(c(1, 0), c(4, 2)),
                        rep(c(0, 1), c(2, 4))), 6)
modules_ztrans = c(0.3, 0.5, 0.5)
mod.cor = calcExpectedMatrixFactors(true_factors, modules_ztrans)
sds = runif(6, 1, 3)
mod_cov = outer(sds, sds) * mod.cor

modules = matrix(c(rep(c(1, 0), c(4, 2)),
                   rep(c(0, 1), c(2, 4)),
                   rep(c(1, 0, 0), each = 2),
                   rep(c(0, 0, 1), each = 2)), 6)
pop = rmvnorm(100, sigma = mod_cov)
hypot = list(tudo = modules, true = modules[,-c(3, 4)])

