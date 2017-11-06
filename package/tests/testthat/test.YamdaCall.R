test_that("Yamda calls is ok", {
  data("toadCor")
  data("toadHypo")
  full = hclustHypot(toadCor)
  toadHypo[[5]] = full
  names(toadHypo)[5] = "hclust"
  x = YamdaLM(toadCor, full, 40, T)
})
