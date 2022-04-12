## code to prepare `hgnc` dataset goes here
# hgnc_20210601.txt is download from "https://www.genenames.org/download/custom/" with "Select status"
# only include "Approved" in Jun 1, 2021

hgnc <- read.table("data-raw/internal/hgnc_20210601.txt", sep = "\t", header = T, comment.char = "", quote = "")
hgnc <- hgnc[, c("Approved.symbol", "Previous.symbols")]
usethis::use_data(hgnc, overwrite = TRUE, internal = TRUE)
