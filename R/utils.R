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
#' Filter gene list according to genes expressed (in your count matrix), you should run filter_genes() first!
#'
#' @param expressed.gene genes expressed in your count matrix, can be `rownames(norm.data)`
#' @param used.gsets A data.frame, gene sets to be cleaned, can be keggMetabDf or reacMetabDf. You can also provide genes sets
#' by yourself, make sure each row represents a genes, and columns "GeneSymbol" (up-to-date HGNC symbols) and
#' "Pathway" (Pathway names) are included.
#' @param min.size gene sets with genes less than min.size will be removed
#' @param suffix logical, whether to add gene set size as suffix to gene set name, default is TRUE
#'
#' @return a list object, each element represents a cleaned gene sets.
#' @export
#'
clean_glist <- function(expressed.gene,
                        used.gsets,
                        min.size = 5,
                        suffix = TRUE
) {
    glist <- split(used.gsets$GeneSymbol, f = used.gsets$Pathway)
    glist.filtered <- sapply(glist, function(x) x[x %in% expressed.gene])

    gset_len <- sapply(glist.filtered, length)
    if (suffix) {
        names(glist.filtered) <- paste0(names(glist.filtered), " (", gset_len, ")")
    }
    glist.filtered <- glist.filtered[gset_len >= min.size]

    print("Note:")
    print(paste("There were", length(unique(used.gsets$GeneSymbol)), "metab genes and",
                length(expressed.gene), "expressed genes."))
    print(paste("After filtered,", length(unique(unlist(glist.filtered))), "expressed metab genes keeped!"))
    return(glist.filtered)
}



### plot function
#' Title
#'
#' @param pas a list containing different type of pathway activity score
#' @param meta.data a data.frame containing metadata of cells
#' @param out.dir output directory for plot file
#'
#' @importFrom grDevices dev.off pdf
#' @importFrom rlang .data
#' @export
#'
plotPath_x_libsize <- function(pas,
                               meta.data,
                               out.dir) {
    if (is.null(pas.tools)) {
        pas.tools <- names(pas)
    }
    if (is.null(out.dir)) {
        stop("You should specify output directory!")
    }
    for (tools in pas.tools) {
        pdf(paste0(out.dir, "/", tools, ".pdf"))
        df <- merge(pas[[tools]], meta.data, by = "row.names")
        for (pathway in colnames(pas[[tools]])) {
            print(ggpubr::ggscatter(df, x = pathway, y = .data$libsize_type,
                            add = "reg.line", conf.int = TRUE,
                            add.params = list(fill = "red3")) +
                      ggpubr::stat_cor(method = "spearman"))
        }
        dev.off()
    }
}


#' Pathway-Pathway Plot from Different PAS Methods
#'
#'
#'
#' @param pas a list containing different type of pathway activity score
#' @param out.dir output directory for plot file
#' @param group.by variable used to group data points
#'
#' @export
#'
plotPath_x_Path <- function(pas,
                            out.dir = NULL,
                            group.by = NULL
) {
    cellname <- rownames(pas[[1]])
    if (is.null(group.by)) {
        color = "black"
    } else {
        color = group.by
    }

    for (i in 1:(length(pas)-1)) {
        for (j in (i+1):length(pas)) {
            pdf(paste0(out.dir, "/", names(pas)[i], "_x_", names(pas)[j], ".pdf"))
            con_path <- intersect(colnames(pas[[i]]), colnames(pas[[j]]))
            for (pathway in con_path) {
                df <- data.frame(tools1 = pas[[i]][, pathway], tools2 = pas[[j]][, pathway], color = color)
                print(ggpubr::ggscatter(df, x = "tools1", y = "tools2",
                                        color = "color",
                                        add = "reg.line", conf.int = TRUE,
                                        add.params = list(fill = "red3"),
                                        xlab = names(pas)[i], ylab = names(pas)[j]) +
                          ggpubr::stat_cor(method = "spearman") +
                          ggplot2::ggtitle(pathway))
            }
            dev.off()
        }

    }
}



