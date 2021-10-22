# =========================================================================
#### An Modified Version of pagoda2::Pagoda2 Class. (could be fold up)
# =========================================================================
######

#' @importFrom stats qnorm
#' @import Matrix
NULL

# The default pagoda2::Padoda2$testPathwayOverdispersion() only calculated pathway scores which are significant. So
# the statement `if (avar > 0.5) ...` is commented and all pathway scores could be calculated.
#
# Note: This code is from pagoda2 package with little modification.
mod_Pagoda2 <- R6::R6Class("mod_Pagoda2",
    inherit = pagoda2::Pagoda2,
    public = list(

        testPathwayOverdispersion2 = function(setenv, type='counts', max.pathway.size=1e3, min.pathway.size=10,
            n.randomizations=5, verbose=FALSE, n.cores=self$n.cores, score.alpha=0.05,
            plot=FALSE, cells=NULL, adjusted.pvalues=TRUE,
            z.score = qnorm(0.05/2, lower.tail = FALSE), use.oe.scale = FALSE,
            return.table=FALSE, name='pathwayPCA',
            correlation.distance.threshold=0.2, loading.distance.threshold=0.01,
            top.aspects=Inf, recalculate.pca=FALSE, save.pca=TRUE) {

            if (!requireNamespace("scde", quietly=TRUE)){
                stop("You need to install package 'scde' to be able to use testPathwayOverdispersion().")
            }

            nPcs <- 1
            if (type=='counts') {
                x <- self$counts
                # apply scaling if using raw counts
                x@x <- x@x*rep(self$misc[['varinfo']][colnames(x),'gsf'],diff(x@p))
            } else {
                if (!type %in% names(self$reductions)) { stop("Reduction ",type,' not found')}
                x <- self$reductions[[type]]
            }
            if (!is.null(cells)) {
                x <- x[cells,]
            }

            proper.gene.names <- colnames(x)

            if (is.null(self$misc[['pwpca']]) || recalculate.pca) {
                if (verbose) {
                    message("determining valid pathways")
                }

                # determine valid pathways
                gsl <- ls(envir = setenv)
                gsl.ng <- unlist(parallel::mclapply(pagoda2:::sn(gsl), function(go) sum(unique(get(go, envir = setenv)) %in% proper.gene.names),mc.cores=n.cores,mc.preschedule=TRUE))
                gsl <- gsl[gsl.ng >= min.pathway.size & gsl.ng<= max.pathway.size]
                names(gsl) <- gsl

                if (verbose) {
                    message("processing ", length(gsl), " valid pathways")
                }

                cm <- Matrix::colMeans(x)

                pwpca <- pagoda2:::papply(gsl, function(sn) {
                    lab <- proper.gene.names %in% get(sn, envir = setenv)
                    if (sum(lab)<1) {
                        return(NULL)
                    }
                    pcs <- irlba::irlba(x[,lab], nv=nPcs, nu=0, center=cm[lab])
                    pcs$d <- pcs$d/sqrt(nrow(x))
                    pcs$rotation <- pcs$v
                    pcs$v <- NULL

                    # get standard deviations for the random samples
                    ngenes <- sum(lab)
                    z <- do.call(rbind,lapply(seq_len(n.randomizations), function(i) {
                        si <- sample(ncol(x), ngenes)
                        pcs <- irlba::irlba(x[,si], nv=nPcs, nu=0, center=cm[si])$d
                    }))
                    z <- z/sqrt(nrow(x))

                    # local normalization of each component relative to sampled PC1 sd
                    avar <- pmax(0, (pcs$d^2-mean(z[, 1]^2))/sd(z[, 1]^2))

                    #if (avar>0.5) {
                    # flip orientations to roughly correspond with the means
                    pcs$scores <- as.matrix(t(x[,lab] %*% pcs$rotation) - as.numeric((cm[lab] %*% pcs$rotation)))
                    cs <- unlist(lapply(seq_len(nrow(pcs$scores)), function(i) sign(cor(pcs$scores[i,], Matrix::colMeans(t(x[, lab, drop = FALSE])*abs(pcs$rotation[, i]))))))
                    pcs$scores <- pcs$scores*cs
                    pcs$rotation <- pcs$rotation*cs
                    rownames(pcs$rotation) <- colnames(x)[lab]
                    #} # don't bother otherwise - it's not significant
                    return(list(xp=pcs,z=z,n=ngenes))
                }, n.cores = n.cores,mc.preschedule=TRUE)
                if (save.pca) {
                    self$misc[['pwpca']] <- pwpca
                }
            } else {
                if (verbose) {
                    message("reusing previous overdispersion calculations")
                    pwpca <- self$misc[['pwpca']]
                }
            }

            if (verbose) {
                message("scoring pathway od signifcance")
            }

            # score overdispersion
            true.n.cells <- nrow(x)

            pagoda.effective.cells <- function(pwpca, start = NULL) {
                n.genes <- unlist(lapply(pwpca, function(x) rep(x$n, nrow(x$z))))
                var <- unlist(lapply(pwpca, function(x) x$z[, 1]))
                if (is.null(start)) { start <- true.n.cells*2 } # start with a high value
                of <- function(p, v, sp) {
                    sn <- p[1]
                    vfit <- (sn+sp)^2/(sn*sn+1/2) -1.2065335745820*(sn+sp)*((1/sn + 1/sp)^(1/3))/(sn*sn+1/2)
                    residuals <- (v-vfit)^2
                    return(sum(residuals))
                }
                x <- stats::nlminb(objective = of, start = c(start), v = var, sp = sqrt(n.genes-1/2), lower = c(1), upper = c(true.n.cells))
                return((x$par)^2+1/2)
            }
            n.cells <- pagoda.effective.cells(pwpca)

            vdf <- data.frame(do.call(rbind, lapply(seq_along(pwpca), function(i) {
                vars <- as.numeric((pwpca[[i]]$xp$d))
                cbind(i = i, var = vars, n = pwpca[[i]]$n, npc = seq(1:ncol(pwpca[[i]]$xp$rotation)))
            })))

            # fix p-to-q mistake in qWishartSpike
            qWishartSpikeFixed <- function (q, spike, ndf = NA, pdim = NA, var = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)  {
                params <- RMTstat::WishartSpikePar(spike, ndf, pdim, var, beta)
                qnorm(q, mean = params$centering, sd = params$scaling, lower.tail, log.p)
            }

            # add right tail approximation to ptw, which gives up quite early
            pWishartMaxFixed <- function (q, ndf, pdim, var = 1, beta = 1, lower.tail = TRUE) {
                params <- RMTstat::WishartMaxPar(ndf, pdim, var, beta)
                q.tw <- (q - params$centering)/(params$scaling)
                p <- RMTstat::ptw(q.tw, beta, lower.tail, log.p = TRUE)
                p[p == -Inf] <- stats::pgamma((2/3)*q.tw[p == -Inf]^(3/2), 2/3, lower.tail = FALSE, log.p = TRUE) + lgamma(2/3) + log((2/3)^(1/3))
                p
            }

            vshift <- 0
            ev <- 0

            vdf$var <- vdf$var-(vshift-ev)*vdf$n
            basevar <- 1
            vdf$exp <- RMTstat::qWishartMax(0.5, n.cells, vdf$n, var = basevar, lower.tail = FALSE)
            #vdf$z <- qnorm(pWishartMax(vdf$var, n.cells, vdf$n, log.p = TRUE, lower.tail = FALSE, var = basevar), lower.tail = FALSE, log.p = TRUE)
            vdf$z <- qnorm(pWishartMaxFixed(vdf$var, n.cells, vdf$n, lower.tail = FALSE, var = basevar), lower.tail = FALSE, log.p = TRUE)
            vdf$cz <- qnorm(pagoda2:::bh.adjust(pnorm(as.numeric(vdf$z), lower.tail = FALSE, log.p = TRUE), log = TRUE), lower.tail = FALSE, log.p = TRUE)
            vdf$ub <- RMTstat::qWishartMax(score.alpha/2, n.cells, vdf$n, var = basevar, lower.tail = FALSE)
            vdf$ub.stringent <- RMTstat::qWishartMax(score.alpha/nrow(vdf)/2, n.cells, vdf$n, var = basevar, lower.tail = FALSE)

            if (plot) {
                test_pathway_par <- par(mfrow = c(1, 1), mar = c(3.5, 3.5, 1.0, 1.0), mgp = c(2, 0.65, 0))
                on.exit(par(test_pathway_par))
                un <- sort(unique(vdf$n))
                on <- order(vdf$n, decreasing = FALSE)
                pccol <- colorRampPalette(c("black", "grey70"), space = "Lab")(max(vdf$npc))
                plot(vdf$n, vdf$var/vdf$n, xlab = "gene set size", ylab = "PC1 var/n", ylim = c(0, max(vdf$var/vdf$n)), col = adjustcolor(pccol[vdf$npc],alpha=0.1),pch=19)
                lines(vdf$n[on], (vdf$exp/vdf$n)[on], col = 2, lty = 1)
                lines(vdf$n[on], (vdf$ub.stringent/vdf$n)[on], col = 2, lty = 2)
            }

            rs <- (vshift-ev)*vdf$n
            vdf$oe <- (vdf$var+rs)/(vdf$exp+rs)
            vdf$oec <- (vdf$var+rs)/(vdf$ub+rs)

            df <- data.frame(name = names(pwpca)[vdf$i], npc = vdf$npc, n = vdf$n, score = vdf$oe, z = vdf$z, adj.z = vdf$cz, stringsAsFactors = FALSE)
            if (adjusted.pvalues) {
                vdf$valid <- vdf$cz  >=  z.score
            } else {
                vdf$valid <- vdf$z  >=  z.score
            }

            if (!any(vdf$valid)) {
                stop("No significantly overdispersed pathways found at z.score threshold of ",z.score)
            }

            # apply additional filtering based on >0.5 sd above the local random estimate
            vdf$valid <- vdf$valid & unlist(lapply(pwpca,function(x) !is.null(x$xp$scores)))
            vdf$name <- names(pwpca)[vdf$i]

            if (return.table) {
                df <- df[vdf$valid, ]
                df <- df[order(df$score, decreasing = TRUE), ]
                return(df)
            }
            if (verbose) {
                message("compiling pathway reduction")
            }
            # calculate pathway reduction matrix

            # return scaled patterns
            xmv <- do.call(rbind, lapply(pwpca[vdf$valid], function(x) {
                xm <- x$xp$scores
            }))

            if (use.oe.scale) {
                xmv <- (xmv -rowMeans(xmv))* (as.numeric(vdf$oe[vdf$valid])/sqrt(apply(xmv, 1, var)))
                vdf$sd <- as.numeric(vdf$oe)
            } else {
                # chi-squared
                xmv <- (xmv-rowMeans(xmv)) * sqrt((qchisq(pnorm(vdf$z[vdf$valid], lower.tail = FALSE, log.p = TRUE), n.cells, lower.tail = FALSE, log.p = TRUE)/n.cells)/apply(xmv, 1, var))
                vdf$sd <- sqrt((qchisq(pnorm(vdf$z, lower.tail = FALSE, log.p = TRUE), n.cells, lower.tail = FALSE, log.p = TRUE)/n.cells))

            }
            rownames(xmv) <- paste("#PC", vdf$npc[vdf$valid], "# ", names(pwpca)[vdf$i[vdf$valid]], sep = "")
            rownames(vdf) <- paste("#PC", vdf$npc, "# ", vdf$name, sep = "")
            self$misc[['pathwayODInfo']] <- vdf

            # collapse gene loading
            if (verbose) {
                message("clustering aspects based on gene loading ... ",appendLF=FALSE)
            }
            tam2 <- pagoda2::pagoda.reduce.loading.redundancy(list(xv=xmv,xvw=matrix(1,ncol=ncol(xmv),nrow=nrow(xmv))),pwpca,NULL,plot=FALSE,distance.threshold=loading.distance.threshold,n.cores=n.cores)
            if (verbose) {
                message(nrow(tam2$xv)," aspects remaining")
            }
            if (verbose) {
                message("clustering aspects based on pattern similarity ... ",appendLF=FALSE)
            }
            tam3 <- pagoda2::pagoda.reduce.redundancy(tam2, distance.threshold=correlation.distance.threshold,top=top.aspects)
            if (verbose) {
                message(nrow(tam3$xv)," aspects remaining\n")
            }
            tam2$xvw <- tam3$xvw <- NULL # to save space
            tam3$env <- setenv

            # clean up aspect names, as GO ids are meaningless
            names(tam3$cnam) <- rownames(tam3$xv) <- paste0('aspect',1:nrow(tam3$xv))

            self$misc[['pathwayOD']] <- tam3
            self$reductions[[name]] <- tam3$xv
            invisible(tam3)
        }
    )
)


#### Codes above could be be fold up) ======
