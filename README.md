
# scMetab

<!-- badges: start -->
<!-- badges: end -->

The goal of scMetab is to ...

## Installation

You can install the released version of scMetab using devtools:

First click `scMetab.Rproj`, then: 

``` r
library(devtools)
install()
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(scMetab)
## basic example code

# percents <- ClusterPercent(seu_obj, assay = "RNA", group.by = "celltype")
# max_percent <- sort(apply(percents, 1, max), decreasing = TRUE)
# summary(max_percent)

norm.data <- GetAssayData(seu_obj, assay = "RNA", slot = "data")
# norm.data <- norm.data[names(max_percent)[1:5000], ]
# dim(norm.data)

gaude.list <- clean_glist(rownames(norm.data), used.gsets = gaudeMetabDf, min.size = 5)
metab.list <- gaude.list
metab.gene <- intersect(unique(gaudeMetabDf$GeneSymbol), rownames(seu_obj))

# PAS ======
pas_tools <- c("aucell", "mean", "pagoda2")
AUCell::plotGeneCount(as.matrix(norm.data))
pas_list = sapply(X = pas_tools, FUN = cal_pas, data = norm.data, glist = metab.list,
                  n_cores = 4, seed.use = seed.use, aucMaxRankPercent = 'auto',
                  simplify = FALSE)
seu_obj[["AUCell"]] <- CreateAssayObject(counts = t(pas_list$aucell))
```

