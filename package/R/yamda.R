#' Yet another modularity detection algorith in R
#'
#' Uses hierarchical partitions and linear models to rank modularity hypothesis using correlation matrices.
#'
#' @export
#' @param x correlation matrix
#' @param hypot_list list of matrices describing modularity hypothesis. Each element in the list should be one matrix, and each column in these matrices represents a module, and each row a trait. If position [i,j] is one, the trait i belong to module j, and if it is zero, trait i does not belong to module j. Modules can be overlapping.
#' @param n number of individuals used in the correlation matrix estimation
#' @param nneg If true, belonging to the same module can only increase correlations, not decrease them.
#' @importFrom evolqg CreateHypotMatrix
#' @importFrom plyr laply
#' @importFrom penalized penalized coef
#' @examples
#' data(toadCor)
#' data(toadHypo)
#' Yamda(toadCor, toadHypo, 25, nneg = FALSE)
#'
Yamda = function(cor_matrix, data, hypot_list, n, nneg = FALSE){}
