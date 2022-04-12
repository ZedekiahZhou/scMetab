#' Sample Cell Pairs for Comparison
#'
#' Down sample cell pairs from inter-group and intra-group for comparison
#'
#' @param seu_obj seurat object
#' @param group.by variables to group cells, cell from inter-group and intra-group will be compared
#' @param split.by variables to split cells, cells will be split first and then grouped individually
#' @param sample_type weather to sample inter-group cells or intra-group cells, default is both
#' @param n number of cell pairs for each split
#' @param seed.use seed used for subsample
#'
#' @return a data.frame of select cell pairs
#' @export
#'
sample_cell_pairs <- function(seu_obj,
                         group.by = "patientID",
                         split.by = "celltype",
                         sample_type = c("inter", "intra"),
                         n = 1000,
                         seed.use = NULL
) {
    if (length(setdiff(sample_type, c("inter", "intra"))) > 0) {
        stop(paste0("sample_type '", sample_type, "' not support!"))
    }

    if (is.null(seed.use)) {
        set.seed(NULL)
    } else {
        withr::local_seed(seed.use)
    }

    cellinfo <- seu_obj[[c(group.by, split.by)]]
    res <- data.frame()

    for (used_split in unique(cellinfo[[split.by]])) {
        sub_cells <- rownames(cellinfo)[cellinfo[[split.by]] == used_split]
        sub_cellinfo <- cellinfo[sub_cells, group.by, drop = TRUE]
        names(sub_cellinfo) <- sub_cells

        inter_cells <- data.frame() -> intra_cells

        if ("inter" %in% sample_type) {
            for (i in 1:n) {
                cell1 <- sample(sub_cells, 1)
                heters <- sub_cells[sub_cellinfo != sub_cellinfo[cell1]]
                cell2 <- sample(heters, 1)
                inter_cells <- rbind(inter_cells, data.frame(cell1 = cell1, cell2 = cell2,
                                                             cell1.group = sub_cellinfo[cell1],
                                                             cell2.group = sub_cellinfo[cell2]))
            }
            inter_cells$pair <- "inter"
        }

        if ("intra" %in% sample_type) {
            sub_cells_rmUnique <- sub_cells[sub_cellinfo %in% sub_cellinfo[duplicated(sub_cellinfo)]]
            for (i in 1:n) {
                cell1 <- sample(sub_cells_rmUnique, 1)
                homos <- setdiff(sub_cells[sub_cellinfo == sub_cellinfo[cell1]], cell1)
                cell2 <- sample(homos, 1)
                intra_cells <- rbind(intra_cells, data.frame(cell1 = cell1, cell2 = cell2,
                                                             cell1.group = sub_cellinfo[cell1],
                                                             cell2.group = sub_cellinfo[cell2]))
            }
            intra_cells$pair <- "intra"
        }

        cells <- rbind(inter_cells, intra_cells)
        cells$split.by = used_split
        res <- rbind(res, cells)
    }
    return(res)
}



#' Feature Percent by Group
#'
#' Calculate percent of each cells group which expressed specific features
#'
#' @param seu_obj a seurat object
#' @param assay Which assays to use, default is `Seaurat::DefaultAssay(seu_obj)`
#' @param features features to calculate percent, default is all features in the assay
#' @param group.by Categories for grouping, default is not grouping, calculate the percent of whole dataset.
#' @param threshold threshold to detect expressed genes, default is 0
#'
#' @return a matrix of percents with row as features and column as groups
#' @export
#'
ClusterPercent <- function(seu_obj,
                           assay = NULL,
                           features = NULL,
                           group.by = NULL,
                           threshold = 0)
{
    if (!is.null(assay)) Seurat::DefaultAssay(seu_obj) <- assay

    if (is.null(features)) features = rownames(seu_obj)

    if (is.null(group.by)) {
        seu_obj$whole <- "Seurat"
        group.by = "whole"
    }

    # PercentAbove <- function(x, threshold) {
    #     return(length(x = x[x > threshold]) / length(x = x))
    # }

    data.features <- Seurat::GetAssayData(seu_obj[features, ], slot = "counts")

    data.plot <- sapply(X = unique(x = seu_obj[[group.by]][[1]]),
                        FUN = function(ident)
    {
        data.use <- data.features[, colnames(seu_obj)[seu_obj[[group.by]][[1]] == ident]]
        pct.exp <- Matrix::rowMeans(data.use > threshold)
            # apply(X = data.use, MARGIN = 1, FUN = PercentAbove,
            #              threshold = 0)
        return((pct.exp))
        suppressMessages(gc())
    })
}



#' Jaccard Distance of Two Gene Sets
#'
#' @param x gene sets 1 (gs1)
#' @param y gene sets 2 (gs2)
#'
#' @return distance in \[0, 1\]
#' @export
#'
#' @examples
#' jaccard_dist(c("PTEN", "TP53"), c("PTEN", "DUSP2"))
jaccard_dist <- function(x, y) {
    a <- intersect(x, y)
    b <- union(x, y)
    if (length(b) == 0) {
        return(1)
    } else {
        return(1 - length(a)/length(b))
    }
}



#' Pathway Correlation and P Value
#'
#' @param seu_obj a Seuat object with an assay stored pathway activity score (PAS)
#' @param assay assay name of PAS, default is "AUCell"
#' @param features pathways to calculate pair-wise correlation
#' @param cor_method one of "spearman" or "pearson", default is "spearman"
#' @param glist gene set list to correspond to PAS, used to calculate pathway distance
#'
#'
#' @return A list containing the following components:
#' r    the matrix of correlations
#' n    the matrix of number of observations used in analyzing each pair f variables
#' P    the asymptotic P-values
#' dist the matrix of jaccard distance
#' radj correlations adjusted by dist
#' padj adjusted P-values
#' @export
pathway_correlation <- function(seu_obj,
                                assay = "AUCell",
                                features = NULL,
                                cor_method = "spearman",
                                glist =NULL
) {
    Seurat::DefaultAssay(seu_obj) <- assay

    df <- Seurat::FetchData(seu_obj, vars = features)

    # pathway correlation and P-values
    res <- Hmisc::rcorr(as.matrix(df), type = cor_method)

    cell_meta <- setdiff(features, names(glist))

    if (!is.null(cell_meta)) {
        tmp <- vector("list", length = length(cell_meta))
        names(tmp) <- cell_meta
        glist <- c(glist, tmp)
    }

    # pathway dist
    res$dist <- sapply(glist, function(x) {
        sapply(glist, function(y) jaccard_dist(x, y))
    })

    # adjust R by dist
    res$radj <- res$r * res$dist

    # adjust p values
    padj <- stats::p.adjust(res$P, method = "BH")
    res$padj <- matrix(padj, ncol = ncol(res$P), dimnames = dimnames(res$P))

    return(res)
}
