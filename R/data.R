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
#'   \item{HGNC.ID}{A unique ID provided by the HGNC}
#'   \item{Approved.symbol}{The official gene symbol that has been approved by the HGNC and is publicly available.
#'   Symbols are approved based on specific HGNC nomenclature guidelines. }
#'   \item{Approved.name}{The official gene name that has been approved by the HGNC and is publicly available.
#'   Names are approved based on specific HGNC nomenclature guidelines.}
#'   \item{Status}{Indicates whether the gene is classified as: Approved, Entry withdrawn or Symbol withdrawn.
#'   only Approved ones are included in this dataset.}
#'   \item{Previous.symbols}{Symbols previously approved by the HGNC for this gene.
#'   This field can contain multiple values as a comma delimited list.}
#'   \item{Alias.symbols}{Other symbols used to refer to this gene.
#'   This field can contain multiple values as a comma delimited list.}
#'   \item{Chromosome}{Indicates the location of the gene or region on the chromosome.}
#'   \item{Accession.numbers}{Accession numbers for each entry selected by the HGNC.
#'   This field can contain multiple values as a comma delimited list.}
#'   \item{RefSeq.IDs}{The Reference Sequence (RefSeq) identifier for that entry, provided by the NCBI.}
#'   ...
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
#' by `getUpdatedSymbols()`
#'
#' @format A data frame with 1916 rows and 2 variables:
#' \describe{
#'   \item{GeneSymbol}{Gene symbol}
#'   \item{Pathway}{belong to which one of 7 super pathways}
#' }
"reacMetabDf"
