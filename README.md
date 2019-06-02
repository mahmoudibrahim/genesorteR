genesorteR: Feature Ranking for Clustered Single Cell Data
======

**genesorteR is an R package for single cell data analysis. It calculate a specificity score to rank all genes in each cluster in the single cell data. It can then use this ranking to find sets of marker genes or to find highly variable genes. Read more [in genesorteR's pre-print](coming soon...).** 

JAMM was developed at the RWTH Aachen University Hospital.

If you have questions or need help running genesorteR please email us at [this email](http://scr.im/jammpro), we will be happy to help you. For bugs or feature requests, please post [here](https://github.com/mahmoudibrahim/genesorteR/issues).



Install genesorteR
------
```R
#install devtools package from CRAN
install.packages("devtools") 

library(devtools)

#install genesorteR from the Github repository
devtools::install_github("mahmoudibrahim/genesorteR") 
```

Quick Tutorial
------

```R
data(kidneyTabulaMuris) #three cell types from kidney (Tabula Muris data)

#get specificity scores for each cell type
sg = sortGenes(kidneyTabulaMuris$exp, kidneyTabulaMuris$cellType)

head(sg$specScore) #specificity scores for each gene in each cluster

#define a small set of markers
mm = getMarkers(sg, quant = 0.99)

#cluster genes and make a heatmap
pp = plotMarkerHeat(sg$inputMat, sg$inputClass, mm$markers, clusterGenes=TRUE, outs = TRUE)

pp$gene_class_info #gene clusters
```

Also check out the examples in each function's documentation. Vignettes coming soon.


Latest News and Updates
------

* **June 2 2019:** *genesorteR's first release, v0.1.4


genesorteR Documentation
------

Here is the [PDF manual for genesorteR](https://github.com/mahmoudibrahim/genesorteR/genesorteR.pdf). 

You can of course also access the documentation of each function like so `?sortGenes`.





---

*Hint: some real sorting here (click the picture!)...*


[![IMAGE ALT TEXT HERE](http://img.youtube.com/vi/kPRA0W1kECg/0.jpg)](https://www.youtube.com/watch?v=kPRA0W1kECg)
