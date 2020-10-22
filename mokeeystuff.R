load("~/Desktop/mammals.RData")
load("~/Dropbox/labbio/data/cov_bayes_data/Rdatas/monkeys.RData")

ls()
homo = cov2cor(main.data$Homo_sapiens$ed.cov)
nrow(homo)
noFM = grep("FM", colnames(homo))
homo = homo[-noFM, -noFM]
nrow(main.data$Homo_sapiens$info)
homo.yamda = Yamda(x = homo, list(face.neuro = hipoteses_mammals_C[,7:8],
                                  oral.aboboda = hipoteses_mammals_C[,c(2, 6:8)],
                                  total = hipoteses_mammals_C[,-1]), 160)
homo.yamda[[1]]
plot_expected(homo, homo.yamda)
superheat(homo.yamda[[3]][[2]], row.dendrogram = TRUE, col.dendrogram = TRUE)
superheat(homo, row.dendrogram = TRUE, col.dendrogram = TRUE)
homo.yamda[[2]]

(homo.yamda[[2]]$total[-1] + homo.yamda[[2]]$total["background"]) /  homo.yamda[[2]]$total["background"]

papio = cov2cor(main.data$Papio_papio$ed.cov)
nrow(papio)
noFM = grep("FM", colnames(papio))
papio = papio[-noFM, -noFM]
nrow(main.data$papio_sapiens$info)
papio.yamda = Yamda(x = papio, list(face.neuro = hipoteses_mammals_C[,7:8],
                                  oral.aboboda = hipoteses_mammals_C[,c(2, 6:8)],
                                  total = hipoteses_mammals_C[,-1]), 35, T)
papio.yamda[[1]]
papio.yamda[[2]]
(papio.yamda[[2]]$total[-1] + papio.yamda[[2]]$total["background"]) /  papio.yamda[[2]]$total["background"]

pera = cov2cor(as.matrix(read.csv("~/Desktop/Pelamelimorphia_covariance_matrix.csv")[-1]))
rownames(pera) = colnames(pera) = rownames(homo)
pera.yamda = Yamda(x = pera, list(face.neuro = hipoteses_mammals_C[,7:8],
                                    oral.aboboda = hipoteses_mammals_C[,c(2, 6:8)],
                                    total = hipoteses_mammals_C[,-1]), 59, T)
pera.yamda[[1]]
pera.yamda[[2]]
(pera.yamda[[2]]$total[-1] + pera.yamda[[2]]$total["background"]) /  pera.yamda[[2]]$total["background"]
superheat(pera, row.dendrogram = TRUE, col.dendrogram = TRUE)
