#' Get Rid of Outdated Gene Symbols
#'
#' Update previous gene symbols to up-to-date approved gene symbols, return a translation table
#'
#' @details Gene symbols are updated according to "hgnc_20210601.txt" file download from HGNC database (custom downloads).
#' Previous symbols are updated, approved symbols are preserved, unmatched symbols are translate to "".
#'
#' Input symbols are matched against approved symbols first, unmatched ones are then matched against previous
#' symbols. Symbols left unmatched finally were translate to "". Alias symbols were not under consideration
#' to reduce multi-mapping which is confused.
#'
#' If a previous symbol match to multiple approved symbols, or multiple previous symbols match to an approved
#' symbol, the first match will be kept according to the alphabetical order of approved symbols.
#'
#' @param geneSymbols input gene symbols to update, duplicates and NAs are allowed
#'
#' @return a translation table with three variables:
#' - original: original symbols of the same length and order as input
#' - updated: updated symbols:
#' - type: symbols was matched to approved symbols, previous symbols, or deleted
#' @export
#'
#' @examples
#' getUpdatedSymbols(c("KRAS2", "DUSP2", NA))
getUpdatedSymbols <- function(geneSymbols) {
    dup_geneSymbols <- geneSymbols
    geneSymbols <- unique(geneSymbols)
    geneSymbols <- geneSymbols[!is.na(geneSymbols)]

    hgnc.previous <- hgnc %>% tidyr::separate_rows(Previous.symbols, sep = ", ")
    hgnc.previous <- hgnc.previous %>% dplyr::filter(Previous.symbols != "")

    approved <- geneSymbols[geneSymbols %in% hgnc$Approved.symbol] # approved symbol

    not.approved <- setdiff(geneSymbols, approved) # not approved symbol
    previous <- hgnc.previous %>%
        dplyr::filter(Previous.symbols %in% not.approved) %>% # update previous symbol in not approved symbols
        dplyr::filter(!(Approved.symbol %in% approved)) # remove symbols which mismatch to already approved ones
    previous <- previous[order(previous$Approved.symbol, previous$Previous.symbols), ]
    previous <- previous %>%
        dplyr::filter(!duplicated(Approved.symbol)) %>%
        dplyr::filter(!duplicated(Previous.symbols)) # filter multi-map

    # combine approved and previous symbol
    approved <- data.frame(original = approved, updated = approved,
                           type = rep("Approved", length(approved)))
    previous <- data.frame(original = previous$Previous.symbols, updated = previous$Approved.symbol,
                           type = rep("Previous", nrow(previous)))
    res <- rbind(approved, previous)
    # add deleted symbols
    deleted <- setdiff(geneSymbols, res$original)
    deleted <- data.frame(original = deleted, updated = rep("", length(deleted)),
                          type = rep("Deleted", length(deleted)))
    res <- rbind(res, deleted)

    # return data frame
    rownames(res) <- res$original
    trans.table <- res[match(dup_geneSymbols, res$original), ]
    return(trans.table)
}
