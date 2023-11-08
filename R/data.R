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
#' Translate Gaude metabolic pathways to lower case (capitalized)

"gaude_trans"

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

#' @rdname HallmarksDf
"cancerHallmarksDf"

#' @title
#' Human Transcriptional Factor Coding Genes
#'
#' @description
#' 1839 genes, with symbols updated use upadted_symbols()
#'
"TFgenes"


#' @title
#' Colors for Different Cell Type
#'
#' @description
#' used for plot
#'
"celltype_color"

#' @title
#' Colors for Different datasets
#'
#' @description
#' used for plot
#'
"dataset_color"

#' @title
#' Colors for Tumor and Normal
#'
#' @description
#' used for plot
#'
"TvN_color"


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


#' @title
#' Signature Genes of T Cell Subtype
#'
#' @description
#' This data frame contains top 50 signature genes for each T cell subtype identified by Zheng et al. They identified 17 CD8+ and
#' 24 CD4+ metaclusters and corresponding signatures.
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
#'
#' @references
#' Zheng, L. et al. Pan-cancer single-cell landscape of tumor-infiltrating T cells. Science (2021) doi:10.1126/science.abe6474.
#'
"tcellSig"


#' @title
#' Signature Genes of Myeloid
#'
#' @description
#' This data frame contains signatures of myeloid subtype identified by Cheng et al.
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
#'
#' @references
#' Cheng, L. et al. A pan-cancer single-cell transcriptional atlas of tumor infiltrating myeloid cells. Cell (2021)
#' doi:10.1016/j.cell.2021.01.010
#'
"myeloidSig"
