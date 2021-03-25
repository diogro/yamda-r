library(yamdar)
library(mvtnorm)

### Generating matrices using overlaping factors
## Factor model needs less parameters, non-factor model needs more to account for
## correlations in factor-overlaps

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
Yamda(pop, hypot, factors = TRUE, nneg = TRUE)[1:2]
Yamda(pop, hypot, factors = TRUE)[1:2]
Yamda(pop, hypot, factors = FALSE)[1:2]

### Generating matrices using non-overlaping factors
## Both descriptions are pretty much equivalent

modules = matrix(c(rep(c(1, 0, 0), each = 5),
                   rep(c(0, 1, 0), each = 5),
                   rep(c(0, 0, 1), each = 5)), 15)
modules_z_coef = c(0.3, 0.5, .3, 1.7)
mod.cor = calcExpectedMatrixFactors(modules, modules_z_coef)
sds = runif(15, 1, 10)
mod_cov = outer(sds, sds) * mod.cor
pop = rmvnorm(50, sigma = mod_cov)
Yamda(pop, list(modules, modules[,-1]), TRUE)[1:2]
Yamda(pop, list(modules, modules[,-1]), FALSE)[1:2]


### Generating matrices using z-transform scheme of module contributions
## non-factor model nails parameter estimates 

true_factors = matrix(c(rep(c(1, 0), c(4, 2)),
                        rep(c(0, 1), c(2, 4))), 6)
modules_ztrans = c(0.3, 0.5, 0.5)
mod.cor = calcExpectedMatrix(true_factors, modules_ztrans)
sds = runif(6, 1, 3)
mod_cov = outer(sds, sds) * mod.cor

modules = matrix(c(rep(c(1, 0), c(4, 2)),
                   rep(c(0, 1), c(2, 4)),
                   rep(c(1, 0, 0), each = 2),
                   rep(c(0, 0, 1), each = 2)), 6)
pop = rmvnorm(100, sigma = mod_cov)
hypot = list(tudo = modules, true = modules[,-c(3, 4)])
Yamda(pop, hypot, factors = TRUE)[1:2]
Yamda(pop, hypot, factors = FALSE)[1:2]
