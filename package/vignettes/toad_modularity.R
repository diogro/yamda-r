## ---- fig.show='hold'----------------------------------------------------
library(superheat)
library(yamdar)
data("toadCor")
superheat(toadCor, row.dendrogram = T, col.dendrogram = T,
          bottom.label.text.angle = 90, 
          bottom.label.text.size = 3,
          left.label.text.size = 3)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

