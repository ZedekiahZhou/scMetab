#' @title
#' Approved HGNC Gene Symbols
#'
#' @description
#' A dataset containing the HGNC gene symbol information, which is download from
#' "https://www.genenames.org/download/custom/" with "Select status" only include "Approved"
#' in Jun 1, 2021, saved as "data-raw/hgnc_20210601.txt"
#'
#' @format A data frame with 42571 rows and 9 variables:
#' \describe{
#'   \item{Approved.symbol}{The official gene symbol that has been approved by the HGNC and is publicly available.
#'   Symbols are approved based on specific HGNC nomenclature guidelines. }
#'   \item{Previous.symbols}{Symbols previously approved by the HGNC for this gene.
#'   This field can contain multiple values as a comma delimited list.}
#' }
#' @source \url{https://www.genenames.org/download/custom/}
"hgnc"



#' @title
#' Metabolic Genes and Pathways from KEGG
#'
#' @description
#' A dataset containing metabolic genes information from KEGG metabolic pathways, generated from
#' "data-raw/hsa00001.json" files downloaded from KEGG website, which include all kegg pathway info of homo sapiens.
#'
#' @details
#' Note: There **are** duplicated genes in different pathway.
#'
#' @format A data frame with 2958 rows and 10 variables:
#' \describe{
#'   \item{CateID}{KEGG category ID}
#'   \item{Category}{KEGG category name}
#'   \item{PathID}{KEGG pathway ID}
#'   \item{Pathway}{KEGG pathway name}
#'   \item{GeneID}{NCBI gene ID}
#'   \item{GeneSymbol}{HGNC gene symbol}
#'   \item{GeneDescription}{Gene description}
#'   \item{ProteinID}{KEGG protein ID}
#'   \item{ProteinName}{}
#'   \item{ProteinDescription}{}
#'   ...
#' }
"keggMetabDf"


#' @title
#' Metabolic Genes and Super Pathways from KEGG
#'
#' @description
#' A dataset containing metabolic genes information from KEGG metabolic pathways, generated from
#' "data-raw/hsa00001.json" files downloaded from KEGG website, which include all kegg pathway info of homo sapiens.
#'
#' Note: use categories as super pathways
#'
#' @details
#' Note: There **are** duplicated genes in different super pathway.
#'
#' @format A data frame with 2958 rows and 10 variables:
#' \describe{
#'   \item{Category}{KEGG categories}
#'   \item{Pathway}{KEGG super pathway name}
#'   \item{GeneSymbol}{HGNC gene symbol}
#' }
"keggMetabDfSuper"


#' @title
#' Metabolic Genes and Pathways from Reactome
#'
#' @description
#' Seven metabolic super-pathways based on the latest Reactome annotations (Fabregat et al., 2016) were
#' curated by [Peng et al., 2018](https://www.cell.com/cell-reports/fulltext/S2211-1247(18)30438-8).
#' Included were amino acid metabolism (348 genes), carbohydrate metabolism (286 genes),
#' integration of energy (110 genes), lipid metabolism (766 genes), nucleotide metabolism (90 genes),
#' tricarboxylic acid cycle (TCA cycle, 148 genes) and vitamin & cofactor metabolism (168 genes).
#' Generated using Table S2 (save as "data-raw/Reactome_metab_TableS2.xlsx") from this article.
#'
#' @details
#' Note: There **are** duplicated genes in different pathway. Out-of-date gene symbols were updated
#' by `update_symbols()`
#'
#' @format A data frame with 1916 rows and 2 variables:
#' \describe{
#'   \item{GeneSymbol}{Gene symbol}
#'   \item{Pathway}{belong to which one of 7 super pathways}
#' }
"reacMetabDf"


#' @title
#' Metabolic Genes and Pathways from Gaude et. al.
#'
#' @description
#' Manually curated metabolic pathways from Gaude et. al. (https://www.nature.com/articles/ncomms13041#Sec7)
#'
#' @details
#' Note: There **are** duplicated genes in different pathway. Out-of-date gene symbols were updated
#' by `update_symbols()`
#'
#' @format A data frame with 1932 rows and 2 variables:
#' \describe{
#'   \item{GeneSymbol}{Gene symbol}
#'   \item{Pathway}{Pathway name}
#' }
"gaudeMetabDf"


#' @title
#' Cancer Hallmark Genes and Pathways from MSigDB.
#'
#' @description
#' Cancer Hallmark Genes and Pathways from MSigDB.
#' `HallmarksDf` contains all MSigDB hallmark gene sets
#' `cancerHallmarksDf` contains 6 cancer hallmarks and mTORC1 signaling pathway
#'
#' @details
#' Note: There **are** duplicated genes in different pathway. Out-of-date gene symbols were updated
#' by `update_symbols()`
#'
#' @format
#' \describe{
#'   \item{GeneSymbol}{Gene symbol}
#'   \item{Pathway}{Pathway name}
#' }
"HallmarksDf"
"cancerHallmarksDf"

#' @title
#' Human Transcriptional Factor Coding Genes
#'
#' @description
#' 1839 genes, with symbols updated use upadted_symbols()
#'
"TFgenes"


#' @title
#' Immune Genes
#'
#' @description
#' Immune Genes
#'
#' @details
#' Note: There **are** duplicated genes in different pathway. Out-of-date gene symbols were updated
#' by `update_symbols()`
#'
#' @format
#' \describe{
#'   \item{GeneSymbol}{Gene symbol}
#'   \item{Pathway}{Pathway name}
#' }
"immuneDf"
