## Filter Genes with High Drop-out Rate
#
#' Filter Genes with High Drop-out Rate
#'
#' @param seu_obj input seurat object
#' @param percent_min minimum proportion of cells expressing a gene, default is 5, means 5% of cells.
#' Genes expressed in less than ncol(seu_obj) * percent_min / 100 cells will be removed.
#'
#' @return a seurat object with high drop-out genes removed.
#' @export
#'
#' @examples filter_genes(pbmc_small)
filter_genes <- function(seu_obj,
                        percent_min = 5
) {
    message("Caution: If there is SCT assay in your seurat object, it will be removed! You should re-run
            Seurat::sctransform() due to the change of library size.")
    SeuratObject::DefaultAssay(seu_obj) <- "RNA"
    if (!is.null(seu_obj@assays$SCT)) seu_obj[["SCT"]] <- NULL

    counts <- SeuratObject::GetAssayData(seu_obj, slot = "counts", assay = "RNA") > 0
    nCells <- Matrix::rowSums(counts)
    seu_obj <- seu_obj[nCells >= ncol(seu_obj) * percent_min / 100, ]
    print(paste("There were", nrow(counts), "genes before filtering, and", nrow(seu_obj), "left after filtering."))
    gc()
    return(seu_obj)
}


## pagoda2 pathway
#' Title
#'
#' @param data
#' @param glist
#' @param n_cores
#'
#' @return
#' @export
#'
#' @examples
cal_pagoda2 <- function(data,   # normalized count matrix
                        glist,  # gene set list
                        n_cores # number of cores to used
) {
    message(paste("Begin calculating pagoda2 score at", Sys.time()))
    # pre-processing
    p2 <- pagoda2::Pagoda2$new(data, n.cores = n_cores, log.scale = FALSE)
    pagoda2::p2$adjustVariance(plot = F)
    pagoda2::p2$calculatePcaReduction(nPCs = 5, use.odgenes = FALSE, fastpath = FALSE)

    # get pathway pca
    kegg.env <- list2env(glist)
    pagoda2::p2$testPathwayOverdispersion(setenv = kegg.env, verbose = T,
                                 recalculate.pca = T, return.table = F,
                                 min.pathway.size = 1)

    # un-significant pathway were null, remove it
    path_scores <- sapply(names(glist), function(i) {p2$misc$pwpca[[i]]$xp$scores})
    path_scores <- path_scores[!sapply(path_scores, is.null)]
    path_scores <- sapply(path_scores, identity)
    rownames(path_scores) <- colnames(data)

    message(paste("Calculate pagoda2 score over at", Sys.time()))
    gc()
    return(as.data.frame(path_scores))
}


## pagoda
#' Title
#'
#' @param data
#' @param glist
#' @param n_cores
#'
#' @return
#' @export
#'
#' @examples
cal_pagoda <- function(data,   # normalized count matrix
                       glist,  # gene set list
                       n_cores # number of cores to used
) {
    message(paste("Begin calculating pagoda score at", Sys.time()))
    # pre-processing
    knn <- scde::knn.error.models(data, n.cores = n_cores,
                                  min.nonfailed = 5, min.count.threshold = 1,
                                  max.model.plots = 10)
    varinfo <- scde::pagoda.varnorm(knn, counts = data, n.cores = n_cores, plot = FALSE)
    varinfo <- scde::pagoda.subtract.aspect(varinfo, Matrix::colSums(data[, rownames(knn)]>0))

    # get pathway pca
    kegg.env <- list2env(glist)
    pwpca <- scde::pagoda.pathway.wPCA(varinfo, kegg.env, n.components = 1, n.cores = n_cores,
                                       min.pathway.size = 1, max.pathway.size = 5000)
    path_scores <- sapply(pwpca, function(x) x$xp$scores)
    rownames(path_scores) <- colnames(data)

    message(paste("Calculate pagoda score over at", Sys.time()))
    gc()
    return(as.data.frame(path_scores))
}


## AUCell
#' Title
#'
#' @param data
#' @param glist
#' @param n_cores
#'
#' @return
#' @export
#'
#' @examples
cal_AUCell <- function(data,   # normalized count matrix
                       glist,  # gene set list
                       n_cores # number of cores to used
) {
    message(paste("Begin calculating AUCell score at", Sys.time()))
    ac_ranking <- AUCell::AUCell_buildRankings(data, nCores = n_cores, plotStats = F, verbose = T)
    sc_AUC <- AUCell::AUCell_calcAUC(glist, ac_ranking, normAUC = T,
                             aucMaxRank=ceiling(0.05 * nrow(ac_ranking)), verbose = T)

    message(paste("Calculate AUCell score over at", Sys.time()))
    gc()
    return(as.data.frame(t(getAUC(sc_AUC))))
}


### ssgsea
#' Title
#'
#' @param data
#' @param glist
#' @param n_cores
#'
#' @return
#' @export
#'
#' @examples
cal_ssgsea <- function(data,   # normalized count matrix
                       glist,  # gene set list
                       n_cores # number of cores to used
) {
    message(paste("Begin calculating ssgsea score at", Sys.time()))
    scores <- GSVA::gsva(data, glist, method = "ssgsea", parallel.sz = n_cores)

    message(paste("Calculate ssgsea score over at", Sys.time()))
    gc()
    return(as.data.frame(t(scores)))
}


### plage
#' Title
#'
#' @param data
#' @param glist
#' @param n_cores
#'
#' @return
#' @export
#'
#' @examples
cal_plage <- function(data,   # normalized count matrix
                      glist,  # gene set list
                      n_cores # number of cores to used
) {
    message(paste("Begin calculating plage score at", Sys.time()))
    scores <- GSVA::gsva(data, glist, method = "plage", parallel.sz = n_cores)

    message(paste("Calculate plage score over at", Sys.time()))
    gc()
    return(as.data.frame(t(scores)))
}


### mean
#' Title
#'
#' @param data
#' @param glist
#' @param n_cores
#'
#' @return
#' @export
#'
#' @examples
cal_mean <- function(
    data,   # normalized count matrix
    glist,  # gene set list
    n_cores # number of cores to used
) {
    message(paste("Begin calculating mean score at", Sys.time()))
    scores <- sapply(glist, function(x) {colMeans(data[x, ])})

    message(paste("Calculate mean score over at", Sys.time()))
    return(scores)
}


### plot function
#' Title
#'
#' @param pas
#' @param meta.data
#' @param out.dir
#'
#' @return
#' @export
#'
#' @examples
plotPath_x_libsize <- function(pas,
                               meta.data,
                               out.dir) {
    if (is.null(pas.tools)) {
        pas.tools <- names(pas)
    }
    if (is.null(out.dir)) {
        errorCondition("You should specify output directory!")
        return()
    }
    for (tools in pas.tools) {
        pdf(paste0(out.dir, "/", normalization, "_", tools, ".pdf"))
        df <- merge(pas[[tools]], meta.data, by = "row.names")
        for (pathway in colnames(pas[[tools]])) {
            print(ggpubr::ggscatter(df, x = pathway, y = libsize_type,
                            add = "reg.line", conf.int = TRUE,
                            add.params = list(fill = "red3")) +
                      stat_cor(method = "spearman"))
        }
        dev.off()
    }
}


#' Title
#'
#' @param pas
#' @param out.dir
#'
#' @return
#' @export
#'
#' @examples
plotPath_x_Path <- function(pas,
                            out.dir = NULL) {

    for (i in 1:(length(pas)-1)) {
        for (j in (i+1):length(pas)) {
            pdf(paste0(out.dir, "/", names(pas)[i], "_x_", names(pas)[j], ".pdf"))
            con_path <- intersect(colnames(pas[[i]]), colnames(pas[[j]]))
            for (pathway in con_path) {
                df <- data.frame(tools1 = pas[[i]][, pathway], tools2 = pas[[j]][, pathway])
                print(ggpubr::ggscatter(df, x = "tools1", y = "tools2",
                                add = "reg.line", conf.int = TRUE,
                                add.params = list(fill = "red3"),
                                xlab = names(pas)[i], ylab = names(pas)[j]) +
                          stat_cor(method = "spearman") +
                          ggtitle(pathway))
            }
            dev.off()
        }

    }
}



