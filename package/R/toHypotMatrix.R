#' to Hypot matrix
#'
#' Transform membership vector to a matrix of modularity hypothesis
#'
#' @export
#' @param membership vector of module membership.
toHypotMatrix <- function(membership){
  modules = unique(membership)
  n_modules = length(modules)
  n_traits = length(membership)
  if(is.null(names(membership))) names(membership) = 1:n_traits
  hypot = matrix(0, n_traits, n_modules)
  for(i in 1:n_modules){
    hypot[membership == modules[i], i] = 1
  }
  rownames(hypot) = names(membership)
  colnames(hypot) = modules
  hypot = hypot[,!colSums(hypot) <= 1, drop = FALSE]
  hypot = hypot[,!colSums(hypot) == n_traits, drop = FALSE]
  hypot
}
