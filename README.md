genesorteR: Feature Ranking for Clustered Single Cell Data
======

**genesorteR is an R package for single cell data analysis. It calculate a specificity score to rank all genes in each cluster in the single cell data. It can then use this ranking to find sets of marker genes or to find highly variable genes. Read more [in genesorteR's pre-print](coming soon...).** 

JAMM was developed at the RWTH Aachen University Hospital by Mahmoud M Ibrahim.

If you have questions or need help running genesorteR please email us at [this email](http://scr.im/jammpro), we will be happy to help you. For bugs or feature requests, please post [here](https://github.com/mahmoudibrahim/genesorteR/issues).



Install genesorteR
------
```R
install.packages("devtools") #install devtools package from CRAN
library(devtools)
devtools::install_github("mahmoudibrahim/genesorteR") #install genesorteR from the Github repository.
```

Quick Tutorial
------

```R
data(kidneyTabulaMuris) #three cell types from kidney (Tabula Muris data)
sg = sortGenes(kidneyTabulaMuris$exp, kidneyTabulaMuris$cellType)
head(sg$specScore) #specificity scores for each gene in each cluster

mm = getMarkers(sg, quant = 0.99) #define a small set of markers

#cluster genes and make a heatmap
pp = plotMarkerHeat(sg$inputMat, sg$inputClass, mm$markers, clusterGenes=TRUE, outs = TRUE)
pp$gene_class_info #cell clusters
```

Also check out the examples for each function. Vignettes coming soon.


Latest News and Updates
------

* **June 2 2019:** *genesorteR's first release, v0.1.4


genesorteR Documentation
------

Here is the PDF manual for genesorteR v0.1.4. You can of course also access the documentation of each function like so `?sortGenes`.





---

*Hint: some real sorting here...*
[![IMAGE ALT TEXT HERE](https://i.ytimg.com/vi/kPRA0W1kECg/hqdefault.jpg?sqp=-oaymwEjCPYBEIoBSFryq4qpAxUIARUAAAAAGAElAADIQj0AgKJDeAE=&rs=AOn4CLAPZKv07Ear1Z7rPMJaFiCPWL4_Dw)](https://www.youtube.com/watch?v=kPRA0W1kECg)
