# =========================================================================
#### Several Pathway Activity Score (PAS) Methods
# =========================================================================
#' Functions using different pathway activity scoring methods.
#'
#' @param data normalized count matrix
#' @param glist gene set list
#' @param pas_method pathway activity scoring method used.
#'  \itemize{
#'   \item{pagoda2: }{default method, see `?cal_pagoda2`}
#'   \item{default_pagoda2: }{see `?cal_default_pagoda2`}
#'   \item{pagoda: }{see `?cal_pagoda`}
#'   \item{aucell: }{see `?cal_aucell`}
#'   \item{mean: }{see `?cal_mean`}
#'   \item{ssgsea: }{see `?cal_ssgsea`}
#'   \item{plage: }{see `?cal_plage`}
#' }
#' @param n_cores number of cores to used
#' @param seed.use only used for method "aucell"
#' @param aucMaxRankPercent percent of top genes to calculate AUC, only work for pas_method "aucell"
#'
#' @return A data frame of pathway activity scores, with cells in rows and pathways in columns.
#' @export

cal_pas <- function(data,   # normalized count matrix
                    glist,  # gene set list
                    pas_method = "pagoda2", # pas method used
                    n_cores = 4, # number of cores to used
                    seed.use = NULL,
                    aucMaxRankPercent = 0.05 # only for AUCell, see `?cal_aucell` for detail.
) {
    if (pas_method == "aucell") {
        return(cal_aucell(data = data, glist = glist, n_cores = n_cores,
                          seed.use = seed.use, aucMaxRankPercent = aucMaxRankPercent))
    } else {
        pas_function <- switch(
            EXPR = pas_method,
            "pagoda2" = cal_pagoda2,
            "mean" = cal_mean,
            "default_pagoda2" = cal_default_pagoda2,
            "pagoda" = cal_pagoda,
            "ssgsea" = cal_ssgsea,
            "plage" = cal_plage,
            stop("Unkown PAS method: ", pas_method)
        )
        return(pas_function(data = data, glist = glist, n_cores = n_cores, seed.use = seed.use))
    }


}


### modified pagoda2 pathway
#' Using Modified PAGODA2 to Calculate Pathway Activity Scores
#'
#' Note::Default Pagoda2 method testPathwayOverdispersion was modified to testPathwayOverdispersion2, so
#' all input pathway scores (not just significant pathways) could be calculated.
#'
#' @param data normalized count matrix
#' @param glist gene set list
#' @param n_cores number of cores to used
#' @param seed.use not used
#'
#' @return A data frame of pathway activity scores, with cells in rows and pathways in columns.
#' @export
cal_pagoda2 <- function(data, glist, n_cores = 4, seed.use = NULL) {
    start_time <- Sys.time()
    print(paste("Begin calculating pagoda2 score at", start_time))

    # pre-processing
    mod_p2 <- mod_Pagoda2$new(data, n.cores = n_cores, min.cells.per.gene = 0,
                              min.transcripts.per.cell = -Inf, log.scale = FALSE)
    mod_p2$adjustVariance(plot = F)
    # p2$calculatePcaReduction(nPCs = 5, use.odgenes = FALSE, fastpath = FALSE)

    # get pathway pca
    kegg.env <- list2env(glist)
    mod_p2$testPathwayOverdispersion2(setenv = kegg.env, verbose = T,
                                 recalculate.pca = T, return.table = F,
                                 min.pathway.size = 1, max.pathway.size = 1000)

    # un-significant pathway were null, remove it
    path_scores <- sapply(names(glist), function(i) {mod_p2$misc$pwpca[[i]]$xp$scores})
    # path_scores <- path_scores[!sapply(path_scores, is.null)]
    # path_scores <- sapply(path_scores, identity)
    rownames(path_scores) <- colnames(data)

    stop_time <- Sys.time()
    time_used <- stop_time - start_time
    print(paste0("Calculating pagoda2 score complete at ", stop_time,
                   ", using ", round(time_used, digits = 2), " ", units(time_used)))
    tmpgc <- gc()

    return(as.data.frame(path_scores))
}

### pagoda2 pathway
#' Using Default PAGODA2 to Calculate Pathway Activity Scores
#'
#' @param data normalized count matrix
#' @param glist gene set list
#' @param n_cores number of cores to used
#' @param seed.use not used
#'
#' @return A data frame of pathway activity scores, with cells in rows and pathways in columns.
#' @export
cal_default_pagoda2 <- function(data, glist, n_cores = 4, seed.use = NULL) {
    start_time <- Sys.time()
    print(paste("Begin calculating pagoda2 score at", start_time))

    # pre-processing
    p2 <- pagoda2::Pagoda2$new(data, n.cores = n_cores, log.scale = FALSE)
    p2$adjustVariance(plot = F)
    # p2$calculatePcaReduction(nPCs = 5, use.odgenes = FALSE, fastpath = FALSE)

    # get pathway pca
    kegg.env <- list2env(glist)
    p2$testPathwayOverdispersion(setenv = kegg.env, verbose = T,
                                 recalculate.pca = T, return.table = F,
                                 min.pathway.size = 1, max.pathway.size = 1000)

    # un-significant pathway were null, remove it
    path_scores <- sapply(names(glist), function(i) {p2$misc$pwpca[[i]]$xp$scores})
    path_scores <- path_scores[!sapply(path_scores, is.null)]
    path_scores <- sapply(path_scores, identity)
    rownames(path_scores) <- colnames(data)

    stop_time <- Sys.time()
    time_used <- stop_time - start_time
    print(paste0("Calculating pagoda2 score complete at ", stop_time,
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
#' @param seed.use not used
#'
#' @return A data frame of pathway activity scores, with cells in rows and pathways in columns.
#' @export
cal_pagoda <- function(data, glist, n_cores = 4, seed.use = NULL) {
    start_time <- Sys.time()
    print(paste("Begin calculating pagoda score at", start_time))

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
    print(paste0("Calculating pagoda score complete at ", stop_time,
                   ", using ", round(time_used, digits = 2), " ", units(time_used)))
    return(as.data.frame(path_scores))
}


### AUCell
#' Using AUCell to Calculate Pathway Activity Scores
#'
#' @param data normalized count matrix
#' @param glist gene set list
#' @param n_cores number of cores to used
#' @param seed.use seed for AUCell_buildRankings
#' @param aucMaxRankPercent percent of top genes to calculate AUC, common values may range from 0.01 to 0.2,
#' default is 0.05. It can also be set as 'auto', which will set aucMaxRank as the number of expressed genes
#' at least 99% of cells had.
#'
#' @return A data frame of pathway activity scores, with cells in rows and pathways in columns.
#' @export
cal_aucell <- function(data, glist, n_cores = 4, seed.use = NULL, aucMaxRankPercent = 0.05) {
    start_time <- Sys.time()
    print(paste("Begin calculating AUCell score at", start_time))

    # calculate AUC
    withr::local_seed(seed = seed.use)
    ac_ranking <- AUCell::AUCell_buildRankings(data, nCores = n_cores, plotStats = FALSE, verbose = T)
    tmpgc <- gc()
    if (is.numeric(aucMaxRankPercent)) {
        aucMaxRank <- ceiling(aucMaxRankPercent * nrow(ac_ranking))
    } else {
        if (aucMaxRankPercent == "auto") {
            aucMaxRank <- ac_ranking@nGenesDetected["1%"]
        } else {
            stop("If aucMaxRankPercent is not numeric, it can only be 'auto'!")
        }
    }

    sc_AUC <- AUCell::AUCell_calcAUC(glist, ac_ranking, normAUC = T, nCores = n_cores,
                                     aucMaxRank=aucMaxRank, verbose = T)

    stop_time <- Sys.time()
    time_used <- stop_time - start_time
    print(paste0("Calculating AUCell score complete at ", stop_time,
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
#' @param seed.use not used
#'
#' @return A data frame of pathway activity scores, with cells in rows and pathways in columns.
#' @export
cal_ssgsea <- function(data, glist, n_cores = 4, seed.use = NULL) {
    start_time <- Sys.time()
    print(paste("Begin calculating ssgsea score at", start_time))

    scores <- GSVA::gsva(data, glist, method = "ssgsea", parallel.sz = n_cores)

    stop_time <- Sys.time()
    time_used <- stop_time - start_time
    print(paste0("Calculating ssgsea score complete at ", stop_time,
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
#' @param seed.use not used
#'
#' @return A data frame of pathway activity scores, with cells in rows and pathways in columns.
#' @export
cal_plage <- function(data, glist, n_cores = 4, seed.use = NULL) {
    start_time <- Sys.time()
    print(paste("Begin calculating plage score at", start_time))

    scores <- GSVA::gsva(data, glist, method = "plage", parallel.sz = n_cores)

    stop_time <- Sys.time()
    time_used <- stop_time - start_time
    print(paste0("Calculating plage score complete at ", stop_time,
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
#' @param seed.use not used
#'
#' @return A data frame of pathway activity scores, with cells in rows and pathways in columns.
#' @export
cal_mean <- function(data, glist, n_cores = 4, seed.use = NULL) {
    start_time <- Sys.time()
    print(paste("Begin calculating mean score at", start_time))

    # calculate AUC
    scores <- sapply(glist, function(x) {Matrix::colMeans(data[x, ])})

    stop_time <- Sys.time()
    time_used <- stop_time - start_time
    print(paste0("Calculating mean score complete at ", stop_time,
                   ", using ", round(time_used, digits = 2), " ", units(time_used)))
    tmpgc <- gc()

    return(as.data.frame(scores))
}
