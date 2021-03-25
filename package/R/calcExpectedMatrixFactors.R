#' Calculate expected matrix using modularity Factors
#'
#' Calculates the expected correlation matrix for a given modularity hypothesis 
#' using background and within-modules coefficients
#'
#' @param hypot a single modularity hypothesis
#' @param coef a vector of coefficients for each module, 
#' starting with a background correlation
#' @param c2c Logical. If TRUE, transform resulting matrix to correlation matrix.
#' @export
calcExpectedMatrixFactors <- function(hypot, coef, c2c = TRUE) {
  n_modules = length(coef) - 1
  p = nrow(hypot)
  if(n_modules > 0){
    expected = Reduce("+", Map("*",
                               coef[2:(n_modules+1)],
                               CreateHypotMatrix(hypot)[1:n_modules])) + 
                               coef[1] + diag(p)
    
  } else{
    expected = matrix(coef, p, p) + diag(p)
  }
  if(c2c) suppressWarnings({cov2cor(expected)})
  else expected
}
