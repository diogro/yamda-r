#' Calculate the z-tranform of module correlations using individual values
#'
#' Calculates the z-tranform for a given population and modularity hypothesis
#'
#'@param data individual measurements or residuals
#'@param hypot a single modularity hypothesis
#'@param nneg Logical. If TRUE, coeficients are non-negative?
#'@export
#'@importFrom bbmle logLik AICc mle2
#'@importFrom car logit
#'@importFrom mvtnorm dmvnorm	
fitZTransCoefML <- function(data, hypot, nneg){
  corr_matrix = cor(data)
  initial_params = calcZTransCoef(hypot, corr_matrix, nneg)
  mod_names = names(initial_params)[-1]
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
      -sum(dmvnorm(x, colMeans(x), sigma = outer(x_sds, x_sds) * calcExpectedMatrix(hypot, nn_pars), log = TRUE))
    }
    else
      -sum(dmvnorm(x, colMeans(x), sigma = outer(x_sds, x_sds) * calcExpectedMatrix(hypot, pars), log = TRUE))

  })
  f = make_function(args, body)
  mle2(f, start = as.list(initial_params), data = list(x = data, hypot = hypot, nneg = nneg))
}
