#' @importFrom bbmle logLik AICc
calcModelStatsMLE <- function(m1, current_hypot_name) {
  LL = as.numeric(logLik(m1))
  param = length(coef(m1))
  AICc = AICc(m1)
  data.frame(Hypothesis = current_hypot_name, LL, param, AICc)
}