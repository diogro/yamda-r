#' Calculate expected matrix
#'
#' Calculates the expected correlation matrix for a given modularity hypothesis 
#' using background and within-modules z-transformed correlations
#'
#' @param hypot a single modularity hypothesis
#' @param ztrans_coef a vector of z-transformed correlations for each module, 
#' starting with a background correlation
#' @export
calcExpectedMatrix <- function(hypot, ztrans_coef) {
  n_modules = length(ztrans_coef) - 1
  if(n_modules > 0){
    expected = Reduce("+", Map("*", ztrans_coef[2:(n_modules+1)],
                               CreateHypotMatrix(hypot)[1:n_modules])) + 
                               ztrans_coef[1]

  } else{
    expected = matrix(ztrans_coef, nrow(hypot), nrow(hypot))
  }
  expected = inv_ztrans(expected)
  diag(expected) = 1
  expected
}
