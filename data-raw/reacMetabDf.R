## code to prepare `reacMetabDf` dataset goes here

# Seven metabolic super-pathways based on the latest Reactome annotations (Fabregat et al., 2016) were
# curated by [Peng et al., 2018](https://www.cell.com/cell-reports/fulltext/S2211-1247(18)30438-8).
# Included were amino acid metabolism (348 genes), carbohydrate metabolism (286 genes),
# integration of energy (110 genes), lipid metabolism (766 genes), nucleotide metabolism (90 genes),
# tricarboxylic acid cycle (TCA cycle, 148 genes) and vitamin & cofactor metabolism (168 genes).
# Table S2 of [Peng et al., 2018] was download as "data-raw/Reactome_metab_TableS2.xlsx"

reac_list <- lapply(1:7, readxl::read_xlsx, path = "data-raw/Reactome_metab_TableS2.xlsx", skip = 1)
reac <- do.call(rbind, reac_list)

trans.table <- getUpdatedSymbols(reac$Genes)
table(trans.table$type)

reacMetabDf <- data.frame(GeneSymbol = trans.table$updated, Pathway = reac$Pathway)
usethis::use_data(reacMetabDf, overwrite = TRUE)
