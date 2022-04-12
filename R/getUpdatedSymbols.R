#' Get Rid of Outdated Gene Symbols
#'
#' Update previous gene symbols to up-to-date approved gene symbols, return a translation table
#'
#' @details Gene symbols are updated according to "hgnc_20210601.txt" file download from HGNC database (custom downloads).
#' Previous symbols are updated, approved symbols are preserved, unmatched symbols are also preserved.
#'
#' Input symbols are matched against approved symbols first, unmatched ones are then matched against previous
#' symbols. Symbols left unmatched finally were also preserved, but matched status were recorded.
#' Alias symbols were not under consideration to reduce multi-mapping which is confused.
#'
#' If a previous symbol match to multiple approved symbols, or multiple previous symbols match to an approved
#' symbol, the first match will be kept according to the alphabetical order of approved symbols.
#'
#' @param gene_symbols input gene symbols to update, duplicates and NAs are allowed
#'
#' @return a translation table with three variables:
#' - original: original symbols of the same length and order as input
#' - updated: updated symbols:
#' - type: symbols was matched to approved symbols, previous symbols, or deleted
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom tidyr separate_rows
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' update_symbols(c("KRAS2", "DUSP2", NA))
update_symbols <- function(gene_symbols) {

    dup_gene_symbols <- gene_symbols
    gene_symbols <- unique(gene_symbols)
    gene_symbols <- gene_symbols[!is.na(gene_symbols)]

    hgnc.previous <- hgnc %>% separate_rows(.data$Previous.symbols, sep = ", ")
    hgnc.previous <- hgnc.previous %>% filter(.data$Previous.symbols != "")

    approved <- gene_symbols[gene_symbols %in% hgnc$Approved.symbol] # approved symbol

    not.approved <- setdiff(gene_symbols, approved) # not approved symbol
    previous <- hgnc.previous %>%
        filter(.data$Previous.symbols %in% not.approved) %>% # update previous symbol in not approved symbols
        filter(!(.data$Approved.symbol %in% approved)) # remove symbols which mismatch to already approved ones
    previous <- previous[order(previous$Approved.symbol, previous$Previous.symbols), ]
    previous <- previous %>%
        filter(!duplicated(.data$Approved.symbol)) %>%
        filter(!duplicated(.data$Previous.symbols)) # filter multi-map

    # combine approved and previous symbol
    approved <- data.frame(original = approved, updated = approved,
                           type = rep("Approved", length(approved)))
    previous <- data.frame(original = previous$Previous.symbols, updated = previous$Approved.symbol,
                           type = rep("Previous", nrow(previous)))
    res <- rbind(approved, previous)
    # add deleted symbols
    deleted <- setdiff(gene_symbols, res$original)
    deleted <- data.frame(original = deleted, updated = deleted,
                          type = rep("Deleted", length(deleted)))
    res <- rbind(res, deleted)

    # return data frame
    rownames(res) <- res$original
    trans.table <- res[match(dup_gene_symbols, res$original), ]
    return(trans.table)
}


#' Updated Gene Symbols for Seurat Object
#'
#' Seurat do not allowed renaming features, so get counts data out, rename it and reconstructed seurat object.
#'
#' @param seu_obj Input seurat object
#' @param out_dir output dir for gene symbol updated table
#'
#' @return seurat object with gene symbols updated.
#' @export
#'
update_symbols_seurat <- function(seu_obj, out_dir = "res/") {
    message("Note: only raw count and meta.data will be preserved!")

    # get counts and meta.data
    if (!dir.exists(out_dir)) dir.create(out_dir)

    counts <- SeuratObject::GetAssayData(seu_obj, slot = "counts", assay = "RNA")
    meta.data <- seu_obj@meta.data

    # update gene symbol
    trans.table <- update_symbols(rownames(counts))
    rownames(counts) <- trans.table$updated
    utils::write.table(trans.table, file = file.path(out_dir, "Gene_symbol_updated_table.tsv"),
                       quote = F, sep = "\t", row.names = F)

    # return re-constructed seurat object
    return(SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta.data))
}
