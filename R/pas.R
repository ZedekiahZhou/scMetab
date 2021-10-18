### pagoda2 pathway
#' Using PAGODA2 to Calculate Pathway Activity Scores
#'
#' @param data normalized count matrix
#' @param glist gene set list
#' @param n_cores number of cores to used
#'
#' @return A data frame of pathway activity scores, with cells in rows and pathways in columns.
#' @importFrom pagoda2 Pagoda2
#' @export
cal_pagoda2 <- function(data,   # normalized count matrix
                        glist,  # gene set list
                        n_cores # number of cores to used
) {
    start_time <- Sys.time()
    message(paste("Begin calculating pagoda2 score at", start_time))

    # pre-processing
    p2 <- Pagoda2$new(data, n.cores = n_cores, log.scale = FALSE)
    p2$adjustVariance(plot = F)
    p2$calculatePcaReduction(nPCs = 5, use.odgenes = FALSE, fastpath = FALSE)

    # get pathway pca
    kegg.env <- list2env(glist)
    p2$testPathwayOverdispersion(setenv = kegg.env, verbose = T,
                                 recalculate.pca = T, return.table = F,
                                 min.pathway.size = 1)

    # un-significant pathway were null, remove it
    path_scores <- sapply(names(glist), function(i) {p2$misc$pwpca[[i]]$xp$scores})
    path_scores <- path_scores[!sapply(path_scores, is.null)]
    path_scores <- sapply(path_scores, identity)
    rownames(path_scores) <- colnames(data)

    stop_time <- Sys.time()
    time_used <- stop_time - start_time
    message(paste0("Calculating pagoda2 score complete at ", stop_time,
                   ", using ", round(time_used, digits = 2), " ", units(time_used)))
    tmpgc <- gc()

    return(as.data.frame(path_scores))
}


## pagoda
#' Using PAGODA to Calculate Pathway Activity Scores
#'
#' @param data normalized count matrix
#' @param glist gene set list
#' @param n_cores number of cores to used
#'
#' @return A data frame of pathway activity scores, with cells in rows and pathways in columns.
#' @export
cal_pagoda <- function(data,   # normalized count matrix
                       glist,  # gene set list
                       n_cores # number of cores to used
) {
    start_time <- Sys.time()
    message(paste("Begin calculating pagoda score at", start_time))

    # pre-processing
    genv <- list2env(glist)
    matw <- matrix(1, nrow(data), ncol(data))
    rownames(matw) <- rownames(data)
    colnames(matw) <- colnames(data)
    var <- apply(data, 1, stats::var)
    varinfo <- list('mat' = as.matrix(data), 'matw' = matw, 'arv' = var)
    rm(matw)
    tmpgc <- gc()

    # calculate score
    pwpca <- scde::pagoda.pathway.wPCA(varinfo, genv, n.components = 1, batch.center = FALSE,
                                       min.pathway.size = 5, n.cores = n_cores, verbose = TRUE)
    path_scores <- sapply(pwpca, function(path) {path$xp$scores})
    rownames(path_scores) <- colnames(data)
    tmpgc = gc()

    stop_time <- Sys.time()
    time_used <- stop_time - start_time
    message(paste0("Calculating pagoda score complete at ", stop_time,
                   ", using ", round(time_used, digits = 2), " ", units(time_used)))
    return(as.data.frame(path_scores))
}


### AUCell
#' Using AUCell to Calculate Pathway Activity Scores
#'
#' @param data normalized count matrix
#' @param glist gene set list
#' @param n_cores number of cores to used
#'
#' @return A data frame of pathway activity scores, with cells in rows and pathways in columns.
#' @export
cal_AUCell <- function(data,
                       glist,
                       n_cores
) {
    start_time <- Sys.time()
    message(paste("Begin calculating AUCell score at", start_time))

    # calculate AUC
    ac_ranking <- AUCell::AUCell_buildRankings(data, nCores = n_cores, plotStats = F, verbose = T)
    tmpgc <- gc()
    sc_AUC <- AUCell::AUCell_calcAUC(glist, ac_ranking, normAUC = T, nCores = n_cores,
                                     aucMaxRank=ceiling(0.05 * nrow(ac_ranking)), verbose = T)

    stop_time <- Sys.time()
    time_used <- stop_time - start_time
    message(paste0("Calculating AUCell score complete at ", stop_time,
                   ", using ", round(time_used, digits = 2), " ", units(time_used)))
    tmpgc <- gc()
    return(as.data.frame(t(AUCell::getAUC(sc_AUC))))
}


### ssgsea
#' Using ssgsea to Calculate Pathway Activity Scores
#'
#' @param data normalized count matrix
#' @param glist gene set list
#' @param n_cores number of cores to used
#'
#' @return A data frame of pathway activity scores, with cells in rows and pathways in columns.
#' @export
cal_ssgsea <- function(data,   # normalized count matrix
                       glist,  # gene set list
                       n_cores # number of cores to used
) {
    start_time <- Sys.time()
    message(paste("Begin calculating ssgsea score at", start_time))

    scores <- GSVA::gsva(data, glist, method = "ssgsea", parallel.sz = n_cores)

    stop_time <- Sys.time()
    time_used <- stop_time - start_time
    message(paste0("Calculating ssgsea score complete at ", stop_time,
                   ", using ", round(time_used, digits = 2), " ", units(time_used)))
    tmpgc <- gc()

    return(as.data.frame(t(scores)))
}


### plage
#' Using PLAGE to Calculate Pathway Activity Scores
#'
#' @param data normalized count matrix
#' @param glist gene set list
#' @param n_cores number of cores to used
#'
#' @return A data frame of pathway activity scores, with cells in rows and pathways in columns.
#' @export
cal_plage <- function(data,   # normalized count matrix
                      glist,  # gene set list
                      n_cores # number of cores to used
) {
    start_time <- Sys.time()
    message(paste("Begin calculating plage score at", start_time))

    scores <- GSVA::gsva(data, glist, method = "plage", parallel.sz = n_cores)

    stop_time <- Sys.time()
    time_used <- stop_time - start_time
    message(paste0("Calculating plage score complete at ", stop_time,
                   ", using ", round(time_used, digits = 2), " ", units(time_used)))
    tmpgc <- gc()

    return(as.data.frame(t(scores)))
}


### mean
#' Using Mean as Pathway Activity Scores
#'
#' Using mean of expression of pathway genes as scores.
#'
#' @param data normalized count matrix
#' @param glist gene set list
#' @param n_cores number of cores to used
#'
#' @return A data frame of pathway activity scores, with cells in rows and pathways in columns.
#' @export
cal_mean <- function(
    data,   # normalized count matrix
    glist,  # gene set list
    n_cores # number of cores to used
) {
    start_time <- Sys.time()
    message(paste("Begin calculating mean score at", start_time))

    # calculate AUC
    scores <- sapply(glist, function(x) {Matrix::colMeans(data[x, ])})

    stop_time <- Sys.time()
    time_used <- stop_time - start_time
    message(paste0("Calculating mean score complete at ", stop_time,
                   ", using ", round(time_used, digits = 2), " ", units(time_used)))
    tmpgc <- gc()

    return(as.data.frame(scores))
}
