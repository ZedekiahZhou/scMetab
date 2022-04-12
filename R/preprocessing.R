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


#' Clean Gene List
#'
#' Filter gene list according to genes expressed (in your count matrix), you should filter genes in you matrix first!
#'
#' Gene list could be provided as data.frame (to used.gsets) or list (to glist).
#'
#' @param expressed.gene genes expressed in your count matrix, can be `rownames(norm.data)`
#' @param used.gsets A data.frame represents gene set to be cleaned, can be
#' keggMetabDf, keggMetabDfSuper or reacMetabDf (see ?keggMetabDf for details).
#'
#' You can also provide genes sets by yourself, make sure each row represents a genes, and columns "GeneSymbol"
#' (up-to-date HGNC symbols) and "Pathway" (Pathway names) are included.
#' @param glist provide gene sets as a list
#' @param min.size gene sets with genes less than min.size will be removed
#' @param suffix logical, whether to add gene set size as suffix to gene set name, default is TRUE
#'
#' @return a list object, each element represents a cleaned gene sets.
#' @export
#'
#' @examples
#' genes <- c("HSD17B1", "PLA2G4A", "AACS", "MGLL", "MGST1")
#' # set min.size = 0 only for example
#' tmp <- clean_glist(genes, keggMetabDf, min.size = 0)
clean_glist <- function(expressed.gene,
                        used.gsets = NULL,
                        glist = NULL,
                        min.size = 5,
                        suffix = TRUE
) {
    if (!is.null(used.gsets)) {
        glist <- split(used.gsets$GeneSymbol, f = factor(used.gsets$Pathway, levels = unique(used.gsets$Pathway)))
    } else if (is.null(glist)) {
        stop("You should provide gene sets to either used.gsets of glist!")
    }
    nGene_prefilter <- length(unique(unlist(glist)))

    glist.filtered <- sapply(glist, function(x) x[x %in% expressed.gene])

    gset_len <- sapply(glist.filtered, length)
    if (suffix) {
        names(glist.filtered) <- paste0(names(glist.filtered), " (", gset_len, ")")
    }
    glist.filtered <- glist.filtered[gset_len >= min.size]
    nGene_postfilter <- length(unique(unlist(glist.filtered)))

    print("Note:")
    print(paste("There were", nGene_prefilter, "metab genes and", length(expressed.gene), "expressed genes."))
    print(paste("After filtered,", nGene_postfilter, "expressed metab genes keeped!"))
    return(glist.filtered)
}



