## code to prepare `gaudeMetabDf` dataset goes here

# Table S3 of [Gaude et al., 2016] was download as "data-raw/NC_2016_Gaude_TableS3.xlsx"

# gaude <- readxl::read_xlsx("data-raw/gaudeMetabDf/NC_2016_Gaude_TableS3_mod.xlsx", skip = 1)
gaude <- read.csv("data-raw/gaudeMetabDf/gaude_tolower.csv")
gaude_raw <- readxl::read_excel("data-raw/gaudeMetabDf/NC_2016_Gaude_TableS3_mod.xlsx", skip = 1)
gaude_trans <- data.frame(original = unique(gaude_raw$Pathways), updated = unique(gaude$Pathways))
rownames(gaude_trans) <- gaude_trans$original

trans.table <- update_symbols(gaude$Genes)
table(trans.table$type)

gaudeMetabDf <- data.frame(GeneSymbol = trans.table$updated, Pathway = gaude$Pathways)


usethis::use_data(gaudeMetabDf, overwrite = TRUE)
usethis::use_data(gaude_trans, overwrite = TRUE)
