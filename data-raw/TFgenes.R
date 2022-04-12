## code to prepare `TFgenes` dataset goes here
TFgenes <- read.csv("data-raw/TFgenes/hs_hgnc_tfs_symbols_updated.txt", header = F)$V1

usethis::use_data(TFgenes, overwrite = TRUE)
