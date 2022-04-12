## code to prepare `immuneDf` dataset goes here

immune <- read.csv("data-raw/immunePathway/immunePath.csv")

trans.table <- update_symbols(immune$GeneSymbol)
table(trans.table$type)

immuneDf <- data.frame(GeneSymbol = trans.table$updated, Pathway = immune$Pathway)

usethis::use_data(immuneDf, overwrite = TRUE)
