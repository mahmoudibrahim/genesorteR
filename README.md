genesorteR: Feature Ranking for Single Cell Data
======

**genesorteR is an R package for single cell data analysis. It calculates a specificity score to rank all genes in each cell cluster. It can then use this ranking to find sets of marker genes or to find highly variable genes.** 

**genesorteR is relatively quick, just seconds for 100k cells, few minutes for millions of cells. Read more [in genesorteR's pre-print](https://www.biorxiv.org/content/10.1101/676379v1).** 

genesorteR was developed at the RWTH Aachen University Hospital.

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
library(genesorteR)

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


Fits with Seurat?
------

Yes it does!

```R
#if "seuratObject" is the Seurat object that contains your data, I think this should work:
gs = sortGenes(seuratObject@assays$RNA@data, Idents(seuratObject))
```



Latest News and Updates
------
* **June 25 2019:** genesorteR's Pre-print is now [available on BioArxiv](https://www.biorxiv.org/content/10.1101/676379v1)!
* **June 24 2019:** new release, v0.2.0. Please check release notes.
* **June 2 2019:** genesorteR's first release, v0.1.4.


genesorteR Documentation
------

Here is the [PDF manual for genesorteR](https://github.com/mahmoudibrahim/genesorteR/blob/master/genesorteR.pdf). 

You can of course also access the documentation of each function like so `?sortGenes`.





---

*Hint: some real sorting here (click the picture!)...*


[![IMAGE ALT TEXT HERE](http://img.youtube.com/vi/kPRA0W1kECg/0.jpg)](https://www.youtube.com/watch?v=kPRA0W1kECg)
