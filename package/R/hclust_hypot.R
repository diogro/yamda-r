#' Hierarchical clustering with hclust
#'
#' Generate modularity hypothesis using hierarchical clustering via hclust
#' @export
#' @param x correlation matrix
#' @param ... aditional arguments to hclust
#' @importFrom plyr llply
hclustHypot <- function(x, ...){
  cluster = hclust(dist(x), ...)
  hypot = NULL
  hypot_list = vector("list", length(cluster$height))
  for(i in seq_along(cluster$height)){
    h = cluster$height[i]
    membership = cutree(cluster, h = h)
    current_hypot = toHypotMatrix(membership)
    if(is.null(hypot)) hypot = current_hypot
    else hypot = cbind(hypot, current_hypot)
    hypot = hypot[,!colSums(hypot) <= 1, drop = FALSE]
    hypot <- hypot[, !duplicated(t(hypot))]
    colnames(hypot) = NULL
    #rownames(hypot) = rownames(x)
    hypot_list[[i]] = hypot
  }
  llply(unique(hypot_list), as.matrix)
}
