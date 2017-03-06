#' yamdar
#'
#' @name yamdar
#' @docType package
#' @importFrom methods is
#' @importFrom stats cor cov cov2cor lm quantile reorder residuals rnorm runif sd var cutree dist hclust
#' @importFrom utils write.csv write.table
NULL

#' Example multivariate data stesi
#'
#' Simulated example of 4 continuous bone lengths from 5 species.
#'
#' \itemize{
#' \item humerus
#' \item ulna
#' \item femur
#' \item tibia
#' \item species
#' }
#'
#' @docType data
#' @keywords datasets
#' @name dentus
#' @usage data(dentus)
#' @format A data frame with 300 rows and 5 variables
NULL

#' Tree for dentus example species
#'
#' Hypothetical tree for the species in the dentus data set.
#'
#' @docType data
#' @keywords datasets
#' @name dentus.tree
#' @usage data(dentus.tree)
#' @format ape tree object
NULL


#' Example correlation matrix for a toad skull
#'
#' Correlation matrix from "High evolutionary constraints limited adaptive responses to past climate changes in toad skulls" Monique Nouailhetas Simon, Fabio Andrade Machado, Gabriel Marroig Proc. R. Soc. B 2016 283 20161783; DOI: 10.1098/rspb.2016.1783. Published 26 October 2016
#'
#' @docType data
#' @keywords datasets
#' @name toadCor
#' @usage data(toadCor)
#' @format Matrix
NULL

#' Example modularity hypothesis for the toad skull
#'
#' List of modularity hypothesis matrices for the toad skull
#'
#' @docType data
#' @keywords datasets
#' @name toadHypo
#' @usage data(toadHypo)
#' @format Matrix
NULL
