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
#'@importFrom evolqg ExtendMatrix
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
  
  args = make_alist(c("x", "hypot", "nneg", "background", mod_names))
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
      Sigma = outer(x_sds, x_sds) * ExpectedMatrix(hypot, nn_pars)
    }
    else{
      Sigma = outer(x_sds, x_sds) * ExpectedMatrix(hypot, pars)
    }
    minusLL = -sum(dmvnorm(x, colMeans(x), sigma = Sigma, log = TRUE))
    if(!is.finite(minusLL)){
      eVal = eigen(Sigma)$values
      if(any(eVal < 0)){
        #warning("Expected matrix is no positive definite, bending to correct.")
        Sigma_ext = ExtendMatrix(Sigma, ret.dim = which(eVal < 0)[1] - 1)$ExtMat
        minusLL = -sum(dmvnorm(x, colMeans(x), sigma = Sigma_ext, log = TRUE))
      }
    }
    return(minusLL)
  })
  f = make_function(args, body)
  tryCatch({
    mle2(f, start = as.list(initial_params), data = list(x = data, hypot = hypot, nneg = nneg))
  },
  error = function(e){
    message("Optimization failed, trying with Nelder-Mead.")
    mle2(f, start = as.list(initial_params), data = list(x = data, hypot = hypot, nneg = nneg), method = "Nelder-Mead")
  })
}
