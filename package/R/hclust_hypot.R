#' Hierarchical clustering with hclust
#'
#' Generate modularity hypothesis using hierarchical clustering via hclust
#' @export
#' @param x correlation matrix
#' @param ... aditional arguments to hclust
hclustHypot <- function(x, ...){
  cluster = hclust(dist(x), ...)
  hypot = NULL
  for(h in cluster$height){
    membership = cutree(cluster, h = h)
    if(is.null(hypot)) hypot = toHypotMatrix(membership)
    else hypot = cbind(hypot, toHypotMatrix(membership))
  }
  hypot = hypot[,!colSums(hypot) <= 1, drop = FALSE]
  duplicated.columns <- duplicated(t(hypot))
  hypot <- hypot[, !duplicated.columns]
  hypot
}
