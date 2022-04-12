## code to prepare `gaudeMetabDf` dataset goes here

# Table S3 of [Gaude et al., 2016] was download as "data-raw/NC_2016_Gaude_TableS3.xlsx"

gaude <- readxl::read_xlsx("data-raw/gaudeMetabDf/NC_2016_Gaude_TableS3_mod.xlsx", skip = 1)

trans.table <- update_symbols(gaude$Genes)
table(trans.table$type)

gaudeMetabDf <- data.frame(GeneSymbol = trans.table$updated, Pathway = gaude$Pathways)

usethis::use_data(gaudeMetabDf, overwrite = TRUE)
