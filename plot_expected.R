library(cowplot)
library(yamdar)
library(superheat)
plot_expected = function(original_matrix, yamda_object){
  expected_matrices = yamda_object$expected_matrices
  n_hypot = length(expected_matrices) - 1
  if(is.null(names(expected_matrices))) names(expected_matrices) = 1:n_hypot
  plots = vector("list", n_hypot + 1)
  plots[[1]] = superheat(original_matrix, print.plot = FALSE)[[2]]
  for(i in 1:n_hypot){
    plots[[i+1]] = superheat(expected_matrices[[i]], print.plot = FALSE)[[2]]
  }
  do.call(plot_grid, plots)
}
data("toadCor")
data("toadHypo")
toadHypo[[5]] = hclustHypot(toadCor)
names(toadHypo)[5] = "hclust"
x = Yamda(toadCor, toadHypo, 25, F)
plot_expected(toadCor, x)
#superheat(x$expected_matrices$hclust, row.dendrogram = T, col.dendrogram = T)
