#' Hierarchical clustering with hclust
#'
#' Generate modularity hypothesis using hierarchical clustering via hclust
#' @export
#' @param x correlation matrix
#' @param ... aditional arguments to hclust
#' @importFrom plyr llply
#' @importFrom stats as.dist
hclustHypot <- function(x, ...){
  dist_x = as.dist(sqrt(2*(1-x)^2))
  cluster = hclust(dist_x, ...)
  hypot = NULL
  n = length(cluster$height)
  hypot_list = vector("list", n-1)
  for(i in 1:(n-1)){
    h = cluster$height[n-i]
    membership = cutree(cluster, h = h)
    current_hypot = toHypotMatrix(membership)
    if(is.null(hypot)) hypot = current_hypot
    else hypot = cbind(hypot, current_hypot)
    hypot = hypot[,!colSums(hypot) <= 1, drop = FALSE]
    hypot <- hypot[,!duplicated(t(hypot))]
    colnames(hypot) = NULL
    #rownames(hypot) = rownames(x)
    hypot_list[[i]] = hypot
  }
  llply(unique(hypot_list), as.matrix)
}
