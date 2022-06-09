# maf.R
# ******************************
# takes in large TCGA hg19 dataset and subsets based on 
# project_id (TCGA-CESC, TCGA-HNSC, TCGA-LUSC)

library(tidyverse)

studies <- c("cesc", "hnsc", "lusc")

# Maf - hg19
# ---------------------------
maf <- read_tsv("../data/mafs/mafs-hg19.maf")
for (study in studies) {
  maf %>%
    filter(project_id == paste("TCGA-", toupper(study), sep="")) %>%
    write.table(paste("../data/mafs/TCGA-", toupper(study), ".maf", sep=""),
                quote=F, sep="\t", row.names=F)
}