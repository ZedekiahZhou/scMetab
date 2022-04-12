## code to prepare `cancerHallmarksDf` dataset goes here

# This include six functional cancer hallmarks (i.e., angiogenesis, apoptosis, DNA repair, EMT,
# G2M checkpoint, and inflammatory response) and mTORC1 signaling pathway which has been implicated
# in metabolic dysregulation and cancer development.
# This pathway genes were derived from MsigDB

hallmarks <- read.delim("data-raw/cancerHallmarksDf/h.all.v7.5.symbols.gmt", header = F, fill = T)
rownames(hallmarks) <- gsub("_", "-", sub("HALLMARK_", "", hallmarks$V1))
hallmarks <- hallmarks[, -c(1, 2)]
hallmarks <- as.data.frame(t(hallmarks))

# 6+1 cancer hallmarks
cancerHallmarks <- c("APOPTOSIS", "ANGIOGENESIS", "DNA-REPAIR", "EPITHELIAL-MESENCHYMAL-TRANSITION",
                     "G2M-CHECKPOINT", "INFLAMMATORY-RESPONSE", "MTORC1-SIGNALING")
cancerHallmarks.altname <- c("Apoptosis", "Angiogenesis", "DNA repair", "EMT",
                             "G2M checkpoint", "Inflammatory response", "mTORC1 signaling")
all(cancerHallmarks %in% colnames(hallmarks))

# clean data and convert
Hallmarks.list <- lapply(hallmarks, function(x) x[x != ""])
HallmarksDf <- reshape::melt(Hallmarks.list)
colnames(HallmarksDf) <- c("GeneSymbol", "Pathway")

# update gene symbol
trans.table <- update_symbols(HallmarksDf$GeneSymbol)
table(trans.table$type)


cancerHallmarksDf <- HallmarksDf[HallmarksDf$Pathway %in% cancerHallmarks, ]
cancerHallmarksDf$Pathway <- plyr::mapvalues(cancerHallmarksDf$Pathway,
                                             from = cancerHallmarks, to = cancerHallmarks.altname)

usethis::use_data(HallmarksDf, overwrite = TRUE)
usethis::use_data(cancerHallmarksDf, overwrite = TRUE)
