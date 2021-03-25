#' Calculate the coefficients of module correlations using individual values
#'
#' Calculates the z-transform for a given population and modularity hypothesis
#'
#'@param data individual measurements or residuals
#'@param hypot a single modularity hypothesis
#'@param nneg Logical. If TRUE, coeficients are non-negative?
#'@param factors Use factors or z-transform additive correlations
#'@export
#'@importFrom bbmle logLik AICc mle2
#'@importFrom car logit
#'@importFrom mvtnorm dmvnorm	
fitModuleCoef <- function(data, hypot, nneg, factors){
  
  corr_matrix = cor(data)
  initial_params = calcZTransCoef(hypot, corr_matrix, nneg)
  mod_names = names(initial_params)[-1]
  
  if(factors){
    ExpectedMatrix = calcExpectedMatrixFactors
    initial_params = inv_ztrans(initial_params)
  } else{
    ExpectedMatrix = calcExpectedMatrix
  }

  args = make_alist(c("x", "hypot", "background", mod_names))
  body = quote({
    x_sds = apply(x, 2, sd);
    if(length(mod_names) > 0){
      pars = eval(parse(text = paste0("c(background, ", paste(mod_names, collapse = ", "), ")")))
    } else {
      pars = eval(parse(text = "c(background)"))
    }
    if(nneg){
      nn_pars = pars
      nn_pars[-1] = exp(nn_pars[-1])
      -sum(dmvnorm(x, colMeans(x), sigma = outer(x_sds, x_sds) * ExpectedMatrix(hypot, nn_pars), log = TRUE))
    }
    else{
      -sum(dmvnorm(x, colMeans(x), sigma = outer(x_sds, x_sds) * ExpectedMatrix(hypot, pars), log = TRUE))
    }
  })
  f = make_function(args, body)
  mle2(f, start = as.list(initial_params), data = list(x = data, hypot = hypot, nneg = nneg))
}
