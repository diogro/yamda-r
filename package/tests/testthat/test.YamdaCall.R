test_that("Yamda calls is ok", {
  data("toadCor")
  data("toadHypo")
  full = hclustHypot(toadCor)
  class(toadHypo[[1]])
  toadHypo[[5]] = full
  names(toadHypo)[5] = "hclust"
  x = Yamda(toadCor, full, 40, T)
  x[[3]]
})
