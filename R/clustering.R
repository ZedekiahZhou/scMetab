# ==================================================================================================
# seurat on pas
# ==================================================================================================

#' Default Seurat Clustering Pipeline
#'
#' Packaging seurat default clustering pipeline in a single function.
#'
#' @param seu a Seurat object
#' @param features Vector of features names to scale/center. Default is variable features.
#' @param vars.to.regress variables to regress out. For example, nUMI, or precent.mt
#' @param n_pcs number of PCs used. Default is 30
#' @param useHarmony whether to use harmony to integrate cells from different source. Default is FALSE
#' @param harmony.group.var var to integrate, must be specified if useHarmony is TRUE.
#' @param max.iter.hamony maximum iterations for harmony
#' @param reduction reduction used for dimensional reduction and clustering. Default is "pca". If
#' useHarmony is TRUE, it will be set as "harmony"
#' @param clustering do clustering or not. Default is TRUE.
#' @param resolution resolution parameter for Seurat::FindClusters
#' @param seed.use random seed used. Default is 2021
#' @param useTSNE whether to do TSNE dimensional reduction. Default is TRUE
#' @param verbose control verbosity. Default is TRUE
#' @param fast_tsne_path for windows, specify the path of fastTSNE.exe
#'
#' @return A Seurat object with clustering done.
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData RunPCA
#' @importFrom Seurat RunUMAP RunTSNE FindNeighbors FindClusters
#' @importFrom harmony RunHarmony
#' @export
seurat_pipe <- function(seu,
                        features = NULL,
                        vars.to.regress = NULL,
                        n_pcs = 30,
                        useHarmony = FALSE,
                        harmony.group.var = "samID",
                        max.iter.hamony = 30,
                        reduction = "pca",
                        clustering = TRUE,
                        resolution = c(0.2, 0.5, 0.8, 1.2),
                        seed.use = 2021,
                        useTSNE = TRUE,
                        verbose = TRUE,
                        fast_tsne_path = NULL) {
    seu <- seu %>%
        NormalizeData(verbose = verbose) %>%
        FindVariableFeatures(verbose = verbose) %>%
        ScaleData(verbose = verbose, vars.to.regress = vars.to.regress, features = features) %>%
        RunPCA(npcs = n_pcs, verbose = verbose, seed.use = seed.use, features = features)

    withr::local_seed(seed = seed.use)
    if (useHarmony) {
        seu <- RunHarmony(seu, group.by.vars = harmony.group.var, max.iter.hamony = max.iter.hamony)
        reduction = "harmony"
    }

    seu <- RunUMAP(seu, reduction = reduction, dims = 1:n_pcs, verbose = verbose, seed.use = seed.use)

    if (clustering) {
        seu <- seu %>%
            FindNeighbors(reduction = reduction, dims = 1:n_pcs, verbose = verbose) %>%
            FindClusters(resolution = resolution, random.seed = seed.use, verbose = verbose)
    }


    if (useTSNE) {
        seu <- RunTSNE(seu, reduction = reduction, dims = 1:n_pcs, verbose = verbose,
                       seed.use = seed.use, tsne.method = "FIt-SNE", nthreads = 10,
                       fast_tsne_path = fast_tsne_path)
    }
    return(seu)
}


#' Reorder Cluster Results of Seurat
#'
#' The clusters results identified by Seurat::FindClusters were a factor with levels from 1 to n.
#' However, when singletons (which is cluster consist of only one cell) produced during cluster
#' identifications, the levels of the factor will sort as character (eg: just like "0", "1", "10", "2", "3", ...).
#' This function was used to reorder cluster results as integer (eg: just like "0", "1", "2", "3", ..., "9", "10", "11")
#'
#' @param seu_obj input seurat object
#' @param resolution resolutions of cluster results. The cluster results of Seurat were names as "{assay_name}_snn_res.{resolution}"
#' @param assay_name assay_name of cluster results
#'
#' @return a reorder seurat object
#' @importFrom Seurat AddMetaData
#' @export
reorder_cluster_res <- function(seu_obj,
                                resolution,
                                assay_name = "RNA"
) {
    cluster_res <- seu_obj[[paste0(assay_name, "_snn_res.", resolution)]]
    cluster_res <- lapply(cluster_res, function(x) {factor(as.integer(as.character(x)))})
    cluster_res <- data.frame(cluster_res, row.names = colnames(seu_obj))
    return(AddMetaData(seu_obj, cluster_res))
}




#' Assign Cell Type for Each Cluster According to Marker Expression
#'
#' @param markers cell type marker used for identification. A data frame with two columns, "Markers" & "Celltype"
#' @param seu_obj seurat object used to identify cell type
#' @param cell_score_method one of "binary" or "raw"
#' @param cluster_score_method one of "freq" or "mean". If cell_score_method is "binary", only "mean" supported.
#' If cell_score_method is "raw", use "freq" as default, see detail.
#' @param cluster cell clusters to annotated, default is "seurat_clusters", which means seu_obj$seurat_clusters
#'
#' @details
#' For cell_score_method "binary", use data slot of seu_obj. Cells were scored as 1 if it expresses a marker (eg. counts > 0),
#' otherwise it was scored as 0. Then mean score across cells within a cluster and markers of a cell type were calcalated as
#' celltype score of a cluster. Clusters were assigned to cell type with the highest score.
#'
#' For cell_score_method "raw", use scale.data slot of seu_obj, and cluster_score_method was set as "freq". Cells were scored
#' using the scaled expression of markers, and means of markers of each cell type were calculated as cell type scores.
#' Then, if method is *freq*, annotate each cell to the cell type with highest scores, and annotate each cluster to
#' the cell type with most number of cells. If method is *mean*, calculate the cluster scores as the means of cells of
#' the cluster, then annotate cluster to the cell type with highest scores.
#'
#' @return a data frame annotate each cluster to cell types
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
celltype_assign <- function(markers,
                            seu_obj,
                            cell_score_method = "binary",
                            cluster_score_method = "mean",
                            cluster = "seurat_clusters"
) {
    if (cell_score_method == "binary") {
        markers_exp <- data.frame(t(Seurat::FetchData(seu_obj, vars = markers$Markers, slot = "data") > 0))
    } else if (cell_score_method == "raw") {
        # seu_obj <- Seurat::ScaleData(seu_obj, features = markers$Markers)
        markers_exp <- data.frame(t(Seurat::FetchData(seu_obj, vars = markers$Markers, slot = "data")))
    } else {
        stop(paste("Cell score method", cell_score_method, "not supported!"))
    }

    cell_score <- sapply(split(markers_exp, f = markers$Celltype), Matrix::colMeans)

    if (cluster_score_method == "mean") {
        cluster_score <- tapply(1:nrow(cell_score), seu_obj[[cluster]], simplify = FALSE, function(idx) {
            unlist(Matrix::colMeans(cell_score[idx, ]))
        })
        cluster_score <- do.call(rbind, cluster_score)
        cluster_celltype <- apply(cluster_score, 1, function(x) colnames(cluster_score)[which.max(x)])
        cluster_celltype_df <- data.frame(cluster = names(cluster_celltype), celltype = cluster_celltype)
    } else if (cluster_score_method == "freq") {
        if (cell_score_method == "binary") {
            stop("Cluster score method 'freq' not supported for cell_score_mthod 'binary', otherwise ties will produced!")
        }
        cell_celltype <- apply(cell_score, 1, function(x) colnames(cell_score)[which.max(x)])
        cluster_celltype_table <- as.data.frame(table(unlist(seu_obj[[cluster]]), cell_celltype))
        cluster_celltype <- cluster_celltype_table %>%
            dplyr::group_by(.data$Var1) %>%
            dplyr::filter(.data$Freq == max(.data$Freq)) %>%
            dplyr::arrange(.data$Var1)
        cluster_celltype_df <- data.frame(cluster = cluster_celltype$Var1, celltype = cluster_celltype$cell_celltype)
    } else {
        stop(paste("Cluster score method", cluster_score_method, "not supported!"))
    }

    return(cluster_celltype_df)
}


#' Distinguish Malignant and Non-malignant Cells Iteratively
#'
#' @param seu_obj seurat object contain epithelial cells
#' @param init_list Initial list of malignant and non_malignant features. A list with two elements
#' named as "malignant" and "non_malignant", each contain 50 features.
#' @param out_pdf pdf file for clustering plot of each iteration.
#' @param max_iter maximum iterations, default is 20
#' @param len_feature_list length of feature list to calculate score. Default if 50
#' @param max.cells.per.ident max cells per ident for `Seurat::FindMarkers`, default is 2000
#' @param error_rate is the portion of misclassified cell compared to last iteration if less than
#' error_rate, then it is regarded as convergent.
#'
#' @details
#' First, assign an initial malignant score and non-malignant score for each epithelial cell in
#' `seu_obj` according to init_list via `Seurat::AddModuleScore`. Then, k-means clustering algorithm
#' will be applied on the two scores to define malignant cells and non-malignant cells. Next,
#' Differentially expressed genes between this putative malignant and non-malignant was identified
#' as new features of malignant and non-malignant (50 features each). Finally, repeat the scoring
#' and clustering iteratively until the portion of misclassified cells less than `error_rate`.
#'
#' @return a data.frame of clustering result for each iteration
#' @export
#'
ident_malig <- function(seu_obj,
                        init_list,
                        out_pdf = "ident_malignant_iteration.pdf",
                        max_iter = 20,
                        len_feature_list = 50,
                        max.cells.per.ident = 2000,
                        error_rate = 0.001
) {
    # initial feature list for malignant and non-malignant
    malig_list <- init_list
    malig_res <- data.frame(row.names = colnames(seu_obj))
    i <- 1

    pdf(out_pdf, width = 16, height = 9)
    while (i < max_iter) {
        message(paste("Iteration", i))
        seu_obj <- Seurat::AddModuleScore(seu_obj, features = malig_list, name = names(malig_list))
        malig_scores <- seu_obj[[c("malignant1", "non_malignant2")]]
        kres <- stats::kmeans(malig_scores, centers = 2)
        if (kres$centers[1, "malignant1"] > kres$centers[1, "non_malignant2"]) {
            seu_obj$malig_clu <- ifelse(kres$cluster[colnames(seu_obj)] == 1, "malignant", "non_malignant")
        } else {
            seu_obj$malig_clu <- ifelse(kres$cluster[colnames(seu_obj)] == 2, "malignant", "non_malignant")
        }

        malig_res <- cbind(malig_res, seu_obj$malig_clu)
        colnames(malig_res) <- c(colnames(malig_res)[-ncol(malig_res)], paste0("Iter", i))
        print(ggplot2::ggplot(seu_obj@meta.data, mapping = ggplot2::aes(x = .data$malignant1,
                                                                        y = .data$non_malignant2,
                                                                        col = .data$malig_clu)) +
                  ggplot2::geom_point() +
                  ggplot2::geom_abline(slope = 1, intercept = 0) +
                  ggplot2::ggtitle(paste("Iter", i)) +
                  ggplot2::theme_classic())

        # check convergence
        if (i > 1) {
            gain <- sum(malig_res[, paste0("Iter", i)] != malig_res[, paste0("Iter", i-1)])/nrow(malig_res)
            if (gain < error_rate) {
                message(paste("Get convergence after", i, "iteration!"))
                break
            }
        }

        Seurat::Idents(seu_obj) <- "malig_clu"
        malig_markers <- Seurat::FindMarkers(seu_obj, ident.1 = "malignant", ident.2 = "non_malignant",
                                     max.cells.per.ident = max.cells.per.ident)
        malig_markers <- malig_markers[order(abs(malig_markers$avg_log2FC), decreasing = T), ]
        malig_list <- split(rownames(malig_markers), f = malig_markers$avg_log2FC < 0)
        malig_list <- lapply(malig_list, function(x) x[1:min(length(x), len_feature_list)])
        names(malig_list) <- c("malignant", "non_malignant")

        i <- i+1
    }
    dev.off()
    return(malig_res)
}




#' Seurat Clustering for Pathway Activity Score (PAS)
#'
#' @param pas PAS matrix or data.frame, cells in rows and pathways in columns
#' @param meta_data cell meta info
#' @param meta_var_used fields assumed to exist in meta_data, including:
#'  \itemize{
#'   \item{nCount_RNA: }{total UMI counts of the cell}
#'   \item{nFeature_RNA: }{number of genes expressed in the cell}
#'   \item{percent.mt: }{percentage of counts from mitochondria genes in total UMI counts}
#'   \item{celltype: }{identified main cell types}
#'   \item{sub_celltype: }{identified cell subtypes, optional}
#'   \item{patientID: }{patient ID}
#'   \item{S.Score, G2M.Score, Phase: }{results of Seurat::CellCycleScoring()}
#' }
#' @param n_pcs number of PC used for clustering, default is 5
#' @param resolution resolution for Seurat::FindClusters(), default is 0.1
#' @param seed.use random seed used, default is 2021
#'
#' @return a clustered Seurat object with an assay named "PAS"
#' @export
pathway_seurat_clustering <- function(pas,
                                      meta_data,
                                      meta_var_used = c("nCount_RNA", "nFeature_RNA", "percent.mt",
                                                        "celltype", "sub_celltype", "patientID",
                                                        "S.Score", "G2M.Score", "Phase"),
                                      n_pcs = 5,
                                      resolution = 0.1,
                                      seed.use = 2021
) {
    start_time <- Sys.time()
    print(paste("Begin seurat analysis at", start_time))

    sc <- Seurat::CreateSeuratObject(t(pas), assay = "PAS")
    sc <- sc %>%
        Seurat::ScaleData() %>%
        Seurat::FindVariableFeatures()
    sc <- sc %>%
        Seurat::RunPCA(npcs = n_pcs, seed.use = seed.use, verbose = FALSE) %>%
        Seurat::RunUMAP(dims = 1:n_pcs, seed.use = seed.use) %>%
        Seurat::RunTSNE(dims = 1:n_pcs, seed.use = seed.use)
    print(sc[["pca"]], dims = 1:2, nfeatures = 5)
    sc <- Seurat::FindNeighbors(sc, dims = 1:n_pcs)
    sc <- Seurat::FindClusters(sc, resolution = resolution, random.seed = seed.use)

    # Add meta data
    meta_var_used <- intersect(meta_var_used, colnames(meta_data))
    sc <- Seurat::AddMetaData(sc, meta_data[, meta_var_used])
    if ("celltype" %in% colnames(sc@meta.data)) Seurat::Idents(sc) <- "celltype"

    stop_time <- Sys.time()
    time_used <- stop_time - start_time
    print(paste0("Analysis complete at ", stop_time,
                 ", using ", round(time_used, digits = 2), " ", units(time_used)))
    tmpgc <- gc()
    return(sc)
}


#' Title
#'
#' @param sc a Seurat object constructed with PAS data
#' @param group.by variables to group cell for comparing, default is "celltype"
#' @param only.pos if only find positive markers, default is TRUE
#'
#' @return cell markers of different cell types
#' @export
pathway_seurat_markers <- function(sc,
                                   group.by = "celltype",
                                   only.pos = TRUE) {
    Seurat::Idents(sc) <- "celltype"
    markers <- Seurat::FindAllMarkers(sc, assay = "PAS", slot = "scale.data",
                                      logfc.threshold = 0, min.pct = 0, only.pos = only.pos,
                                      pseudocount.use = 1e-4)

    return(markers)
}


#' Title
#'
#' @param pas PAS matrix or data.frame, cells in rows and pathways in columns
#' @param meta_data cell meta info
#' @param file_prefix file prefix for .pdf plot file
#' @param out_dir output directory, default is "./"
#' @param ... other parameters
#'
#' @return cell type markers
#' @export
pathway_seurat <- function(pas, meta_data,
                           file_prefix = NULL, out_dir = "./",
                           ...

) {
    sc <- pathway_seurat_clustering(pas = pas, meta_data = meta_data)
    pathway_seurat_plot(sc, file_prefix = file_prefix, out_dir = out_dir)
    markers <- pathway_seurat_markers(sc)
    return(markers)
}
