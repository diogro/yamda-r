#' to Membership vector
#'
#' Transform a matrix of modularity hypothesis into a mebership vector
#'
#' @export
#' @param hypot matrix of modularity hypothesis. If position [i,j] is one, the trait i belong to module j, and if it is zero, trait i does not belong to module j. Modules can be overlapping.
#' @param membershipAsFactor Logical. Transform membership vector to factor?
#' @param noModule String. Module to assing traits that don't belong to any modules. Leave empty to use NA
toMembership <- function(hypot, membershipAsFactor = FALSE, noModule = NULL){
  n_modules = ncol(hypot)
  n_traits = nrow(hypot)
  if(is.null(colnames(hypot))) colnames(hypot) = paste0("module_", 1:n_modules)
  if(is.null(rownames(hypot))) rownames(hypot) = 1:n_traits
  membership = rep(NA, n_traits)
  names(membership) = rownames(hypot)
  for(i in 1:n_modules){
    current_module = hypot[,i]
    membership[current_module == 1] = ifelse(is.na(membership[current_module == 1]),
                                             colnames(hypot)[i],
                                             paste(membership[current_module == 1], colnames(hypot)[i], sep = "-"))
  }
  if(!is.null(noModule)) membership[is.na(membership)] = noModule
  if(membershipAsFactor) factor(membership)
  else membership
}
