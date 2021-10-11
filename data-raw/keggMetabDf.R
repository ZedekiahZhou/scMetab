## code to prepare `keggMetabDf` dataset goes here

# "hsa00001.json" were download from KEGG website, which include all kegg pathway info of homo sapiens.

library(rjson)
library(stringr)
kegg <- fromJSON(file = "data-raw/hsa00001.json")

# # the first element if metabolism
metab <- kegg$children[[1]]
metab$name
keggMetabDf <- data.frame()

for (cate in metab$children) {
    print(cate$name)
    for (path in cate$children) {
        cat("\t", path$name, "\n")
        if ("children" %in% names(path)) {
            for (gene in path$children) {
                m1 <- str_match(cate$name, "^([0-9]+) (.+)$")
                m2 <- str_match(path$name, "^([0-9]+) ([^\\[]+)( \\[PATH:hsa[0-9]+\\])?$")
                m3 <- str_match(gene$name, "^([0-9]+) ([A-Z0-9-]+)?;? ?([^\t]+)\t(K[0-9]+) ([^;]+); (.+)$")

                tmpDf <- data.frame('CateID' = m1[2],
                                      'Category' = m1[3],
                                      'PathID' = m2[2],
                                      'Pathway' = m2[3],
                                      'GeneID' = m3[2],
                                      'GeneSymbol' = m3[3],
                                      'GeneDescription' = m3[4],
                                      'ProteinID' = m3[5],
                                      'ProteinName' = m3[6],
                                      'ProteinDescription' = m3[7])
                keggMetabDf <- rbind(keggMetabDf, tmpDf)
            }
        }
    }
}

## multiple modify ==============
# 1. remove "Not included in regular maps"
# 2. remove empty gene symbol
keggMetabDf <- keggMetabDf %>%
    filter(Category != "Not included in regular maps") %>%
    drop_na(GeneSymbol)

# 3. merge category "Metabolism of other amino acids" with "Amino acid metabolism"
keggMetabDf <- keggMetabDf %>%
    mutate(Category = ifelse(test = Category == "Metabolism of other amino acids",
                             yes = "Amino acid metabolism", no = Category))

# 4. some gene to update manually
mtGenes <- data.frame(original = c("ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "CYTB",
                                   "COX1", "COX2", "COX3", "ATP6", "ATP8", "TAZ"),
                      updated = c("MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6", "MT-CYB",
                                  "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6", "MT-ATP8", "TAFAZZIN"))
keggMetabDf <- keggMetabDf %>%
    mutate(GeneSymbol = ifelse(test = GeneSymbol %in% mtGenes$original,
                               yes = mtGenes$updated[match(GeneSymbol, mtGenes$original)],
                               no = GeneSymbol))
all(keggMetabDf$GeneSymbol %in% hgnc$Approved.symbol)

# 5. remove pathway have only 1 gene
keggMetabDf <- keggMetabDf %>%
    group_by(Pathway) %>%
    filter(n() > 1)

usethis::use_data(keggMetabDf, overwrite = TRUE)
