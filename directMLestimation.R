library(evolqg)
library(corrgram)
library(mvtnorm)
library(penalized)
library(bbmle)
library(lucid)
library(yamdar)

modules = matrix(c(rep(c(1, 0, 0), each = 5),
                   rep(c(0, 1, 0), each = 5),
                   rep(c(0, 0, 1), each = 5)), 15)
modules_ztrans = c(0.3, 0.5, 0.6, 1.7)
mod.cor = calcExpectedMatrix(modules, modules_ztrans)

# correlation matrices should be symmetric
sds = runif(15, 1, 10)
mod_cov = outer(sds, sds) * mod.cor

pop = rmvnorm(50, sigma = mod_cov)

random_hypot = matrix(sample(c(0, 1), 3 * 15, replace = TRUE), 15, 3)

# https://stackoverflow.com/questions/12982528/how-to-create-an-r-function-programmatically
make_function <- function(args, body, env = parent.frame()) {
  f <- function() {}
  formals(f) <- args
  body(f) <- body
  environment(f) <- env
  f
}
make_alist <- function(args) {
  res <- replicate(length(args), substitute())
  names(res) <- args
  res
}
calcZTransCoefML <- function (data, hypot, formula = NULL){
  corr_matrix = cor(data)
  initial_params = calcZTransCoef(hypot, corr_matrix, F)
  mod_names = names(initial_params)[-1]
  args = make_alist(c("x", "hypot", "background", mod_names))
  body = quote({
    x_sds = apply(x, 2, sd); 
    pars = eval(parse(text = paste0("c(background, ", paste(mod_names, collapse = ", "), ")")))
    -sum(dmvnorm(x, colMeans(x), sigma = outer(x_sds, x_sds) * calcExpectedMatrix(hypot, pars), log = TRUE))
    })
  f = make_function(args, body)
  m1 = mle2(f, start = as.list(initial_params), data = list(x = data, hypot = hypot))
  m1
}
pop = rmvnorm(100, sigma = mod_cov)
m1 = calcZTransCoefML(data = pop, hypot = modules)
calcZTransCoef(modules, x = cor(pop), F)
AIC(m1)

