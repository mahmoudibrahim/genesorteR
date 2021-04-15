genesorteR: Feature Ranking for Single Cell Data
======

**genesorteR is an R package for single cell data analysis. It calculates a specificity score to rank all genes in each cell cluster. It can then use this ranking to find sets of marker genes or to find highly variable or differentially expressed genes. genesorteR is applicable to scRNA-Seq data as well as other sparse single cell data like scATAC-Seq** 

**genesorteR is relatively quick, just seconds for 100k cells, few minutes for millions of cells. Read more [in genesorteR's pre-print](https://www.biorxiv.org/content/10.1101/676379v2).** 

**If you have questions or need help running genesorteR please email us at [this email](http://scr.im/jammpro), we will be happy to help you. For bugs or feature requests, please post [here](https://github.com/mahmoudibrahim/genesorteR/issues).**

genesorteR was developed at the RWTH Aachen University Hospital.


What genesorteR Can Do
------
* Rank genes (scRNA-Seq) or open chromatin regions (snATAC-Seq) by "specificity" in cell clusters
* Find small sets of marker genes and differentially expressed / highly variable genes or accessible regions
* Cluster genes/open chromatin regions and make heatmap summaries of single cell data


Cited by...
------

* [Cook & Vanderhyden, bioRxiv 2021](https://www.biorxiv.org/content/10.1101/2021.03.05.434142v1.full)
* [Kuppe, Ibrahim et al., Nature 2021](https://www.nature.com/articles/s41586-020-2941-1)
* [Vargo & Gilbert, BMC Bioinformatics 2020](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03641-z)
* [De Boeck, Ahn et al., Nature Communications 2020](https://www.nature.com/articles/s41467-020-18569-4)
* [Lu, Shiau, Yi et al., Dev. Cell 2020](https://doi.org/10.1016/j.devcel.2020.04.009)
* [Berg et al. bioRxiv 2020](https://www.biorxiv.org/content/10.1101/2020.03.31.018820v1)



Install genesorteR
------
```R
#install devtools package from CRAN
install.packages("devtools") 

#install genesorteR from the Github repository
devtools::install_github("mahmoudibrahim/genesorteR") 
```


genesorteR Documentation
------

Here is the [PDF manual for genesorteR](https://github.com/mahmoudibrahim/genesorteR/blob/master/genesorteR.pdf). 

You can of course also access the documentation of each function like so `?sortGenes`.


Wiki (Tutorials & FAQs)
------

* [Comparing genesorteR Specificity Scores to TF-IDF on the Cell Cluster Level in Single Cell Data](https://github.com/mahmoudibrahim/genesorteR/wiki/Specificity-Scores-versus-TF-IDF) (19 Apr. 2020)

* [Finding Genre-specific Actors and Clustering Film Genres](https://github.com/mahmoudibrahim/genesorteR/wiki/Finding-Genre-specific-Actors-and-Clustering-Film-Genres-(data-from-IMDB)) (21 Jan. 2020)

* [Pathway and gene set enrichment analysis in large single cell RNA-Seq data](https://github.com/mahmoudibrahim/genesorteR/wiki/From-Cluster-to-Pathway-Enrichment-in-Large-scRNA-Seq-Data) (02 Jan. 2020)

* [CorrNet: Plot beautiful networks from single cell data clustering results using genesorteR & ggraph](https://github.com/mahmoudibrahim/genesorteR/wiki/Visualize-single-cell-data-in-R-using-genesorteR-&-ggraph) (18 Oct. 2019)

more Wiki pages coming soon...


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
Check [this post](https://github.com/mahmoudibrahim/genesorteR/issues/1) for more info. Also check out the examples in each function's documentation. Vignettes coming soon.

Note that genesorteR does not currently accept expression matrices with negative entries.


Fits with Seurat?
------

Yes it does!

```R
#if "seuratObject" is the Seurat object that contains your data, I think this should work:
gs = sortGenes(seuratObject@assays$RNA@data, Idents(seuratObject))
```


---

*Hint: some real sorting here (click the picture!)...*


[![IMAGE ALT TEXT HERE](http://img.youtube.com/vi/kPRA0W1kECg/0.jpg)](https://www.youtube.com/watch?v=kPRA0W1kECg)
