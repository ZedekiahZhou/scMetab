### plot function

# new geom proto
GeomSplitViolin <- ggplot2::ggproto("GeomSplitViolin",
                                    ggplot2::GeomViolin,
                                    draw_group = function(self, data, ..., draw_quantiles = NULL
) {
    data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x'])
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
        stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
        quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
        aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
        aesthetics$alpha <- rep(1, nrow(quantiles))
        both <- cbind(quantiles, aesthetics)
        quantile_grob <- GeomPath$draw_panel(both, ...)
        ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    }
    else {
        ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
})


#' Split violin plot
#'
#' see ?geom_violin for detail
#'
#' @param mapping mapping
#' @param data data to displayed
#' @param stat stat
#' @param position positon
#' @param ... ...
#' @param draw_quantiles draw_quantiles
#' @param trim trim
#' @param scale scale
#' @param na.rm na.rm
#' @param show.legend show.legend
#' @param inherit.aes inherit.aes
#'
#' @return split violin plot
#' @export
geom_split_violin <- function (mapping = NULL,
                               data = NULL,
                               stat = "ydensity",
                               position = "identity",
                               ...,
                               draw_quantiles = NULL,
                               trim = TRUE,
                               scale = "area",
                               na.rm = FALSE,
                               show.legend = NA,
                               inherit.aes = TRUE) {
    ggplot2::layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
                   position = position, show.legend = show.legend, inherit.aes = inherit.aes,
                   params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


#' Correlation Heatmap for Gene Expression or Reduction
#'
#' This function is used to plot cell-to-cell correlation heatmap based on (metabolic) gene expression
#' or reduction.
#'
#'
#' @param seu_obj seurat object
#' @param reduction Reduction name to be used to calculate cell-to-cell correlation, must can be found
#' in seu_obj. Default is NULL and metabolic gene expression will be used.
#' @param metab.gene Metabolic gene symbols to be used (if reduction is null). If not specified, it will
#' be set to genes in gaudeMetabDf.
#' @param dir.used directory used for *.pdf plot
#' @param file.prefix file prefix for *.pdf plot
#' @param dims Which dimensions to use to calculate correlation (if reduction is not null), default is 1:30.
#' @param slot which slot to use to fetch gene expression (if reduction is null), default is counts.
#' @param vars vars to set as heatmap row/column annotations, must included in seu_obj
#' @param split.by vars to split the plot, each will be plotted individually.
#' @param cor_method correlation method used, default is "spearman".
#' @param n number of cells to subsample for heatmap.
#' @param seed.use random seed used for subsample.
#' @param only.red only use red colors, default is TRUE
#'
#' @return *.pdf plot file
#' @export
cor_heatmap <- function(seu_obj,
                        reduction = NULL,
                        metab.gene = NULL,
                        dir.used = NULL,
                        file.prefix = NULL,
                        dims = 1:30,
                        slot = "counts",
                        vars = c("celltype", "patientID", "nFeature_metab", "nCount_metab"),
                        split.by = "celltype",
                        cor_method = "spearman",
                        n = 500,
                        seed.use = NULL,
                        only.red = TRUE
) {
    metab.df <- gaudeMetabDf
    metab.gene <- intersect(unique(metab.df$GeneSymbol), rownames(seu_obj))

    if (is.null(reduction)) {
        metab.data <- Seurat::GetAssayData(seu_obj[metab.gene, ], slot = slot)
        tmpfile <- paste0("plot/", dir.used, file.prefix, "correlation_heatmap_", split.by, "_MetabGene.pdf")
    } else {
        metab.data <- t(Seurat::Embeddings(seu_obj, reduction = reduction))[dims, ]
        tmpfile <- paste0("plot/", dir.used, file.prefix, "correlation_heatmap_",split.by, "_", reduction, ".pdf")
    }

    if (is.null(seed.use)) {
        set.seed(NULL)
    } else {
        withr::local_seed(seed.use)
    }

    if (only.red) {
        mycolor = grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "Reds"))(100)
    } else {
        mycolor = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
    }

    pdf(tmpfile, width = 12, height = 12)
    for (used_split in unique(seu_obj[[split.by, drop = TRUE]])) {
        print(paste("Plot for", used_split))
        cand_cells <- colnames(seu_obj)[seu_obj[[split.by, drop = TRUE]] == used_split]
        cells <- sample(cand_cells, size = min(length(cand_cells), n))
        annotation_col <- Seurat::FetchData(seu_obj, vars = vars)[cells, ]
        cor.data <- stats::cor(as.matrix(metab.data[, cells]), method = cor_method)
        diag(cor.data) <- 0

        print(pheatmap::pheatmap(cor.data, annotation_col = annotation_col, annotation_row = annotation_col,
                       show_rownames = FALSE, show_colnames = FALSE,
                       color = mycolor,
                       main = used_split))
    }
    dev.off()
}



#' Scatter plot of pathway activity score - library size
#'
#' @param pas a list containing different type of pathway activity score
#' @param pas.tools names of pas tools to plot, if NULL then set as names(pas)
#' @param meta.data a data.frame containing metadata of cells
#' @param libsize_type variable names of library size, one of "nCount_RNA" of "nFeature_RNA"
#' @param out.dir output directory for plot file
#'
#' @importFrom grDevices dev.off pdf
#' @importFrom rlang .data
#' @export
#'
plotPath_x_libsize <- function(pas,
                               pas.tools = NULL,
                               meta.data,
                               libsize_type = "nCount_RNA",
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
            print(ggpubr::ggscatter(df, x = pathway, y = libsize_type,
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



#' Title
#'
#' @param sc a Seurat object constructed with PAS data
#' @param file_prefix file prefix of .pdf plot file
#' @param out_dir directory for output
#' @param reduction_type reduction type used
#' @param group_type variable to group cells
#' @param qc_field QC field to plot using Seurat::FeaturePlot()
#'
#' @export
pathway_seurat_plot <- function(sc, file_prefix = NULL, out_dir = "./",
                                reduction_type = c("tsne", "umap"),
                                group_type = c("celltype", "seurat_clusters", "patientID", "Phase"),
                                qc_field = c("nFeature_RNA", "nCount_RNA", "percent.mt")
) {
    pdf(paste0(out_dir, "/", file_prefix, "_PAS_clustering.pdf"), width = 8, height = 6)
    for (reduction in reduction_type) {
        for (group.by in group_type) {
            print(Seurat::DimPlot(sc, reduction = reduction, group.by = group.by, label = T, repel = T))
        }
    }
    print(Seurat::FeaturePlot(sc, features = qc_field))
    dev.off()
}


#' DoHeatmap of Pathwa Genes
#'
#' Re-scale seurat object using specified pathway_genes
#' @param seu_obj seurat object, must include RNA assay
#' @param pathway_genes gene names of a pathway
#' @param group.by variable names to group the heatmap by, default is "celltype"
#' @param cells cells to used
#' @param colors a vector of colors for different groups
#'
#' @return a heatmap
#' @export
#'
pathway_gene_heatmap <- function(seu_obj,
                                 pathway_genes,
                                 group.by = "celltype",
                                 cells = NULL,
                                 colors = NULL) {
    Seurat::DefaultAssay(seu_obj) <- "RNA"
    seu_obj <- Seurat::ScaleData(seu_obj, features = pathway_genes, verbose = FALSE)
    if (is.null(colors)) {
        print(Seurat::DoHeatmap(seu_obj, features = pathway_genes, group.by = group.by, cells = cells))
    } else {
        suppressMessages(p <- Seurat::DoHeatmap(seu_obj, group.by = "celltype_x_TorN",
                                                features = pathway_genes, cells = cells,
                                                group.colors = colors) +
                            ggplot2::scale_color_manual(values = colors))
        print(p)
    }

}
