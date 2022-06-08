# maf.R
# ******************************
# takes in large TCGA hg19 dataset and subsets based on 
# project_id (TCGA-CESC, TCGA-HNSC, TCGA-LUSC)

library(tidyverse)
# library(TCGAbiolinks)

# ----
# # https://bioconductor.riken.jp/packages/3.4/bioc/vignettes/TCGAbiolinks/inst/doc/mutation.html
# query.maf.hg19 <- GDCquery(
#   project = paste("TCGA-", toupper(study), sep=""),
#   data.category = "Simple nucleotide variation", 
#   data.type = "Simple somatic mutation",
#   # access = "open", 
#   legacy = TRUE, # must be hg19 for oncodrive
#   file.type = paste("gsc_", toupper(study), "_pairs.aggregated.capture.tcga.uuid.automated.somatic.maf", sep="")
# )
# 
# GDCdownload(query.maf.hg19, directory = "../data/GDCdata/")
# maf <- GDCprepare(query.maf.hg19, directory = "../data/GDCdata/")

# maf <- getMC3MAF()
# a <- readRDS("../data/maf-hg38") 
# ----

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


# OncodriveFML formatting 
# ---------------------------
# create variant-type files formatted for submission into oncodrivefml
# http://bbglab.irbbarcelona.org/oncodrivefml/home
variants <- c("3'UTR", "5'UTR")
for (study in studies) {
  for (variant in variants) {
    maf <- read_tsv(paste("../data/mafs/TCGA-", toupper(study), ".maf", sep=""))
    maf %>% 
      filter(Variant_Classification %in% variants) %>%
      select(Chromosome, Start_Position, Reference_Allele, 
             Tumor_Seq_Allele2, Tumor_Sample_Barcode) %>%
      rename(chromosome = Chromosome, position = Start_Position,
             ref = Reference_Allele, alt = Tumor_Seq_Allele2,
             sample = Tumor_Sample_Barcode) %>%
      rename_with(toupper) %>%
      write.table(paste("../data/mafs/TCGA-", toupper(study), 
                        "-" ,toupper(variant),".maf", sep=""),
                  quote=F, sep="\t", row.names=F)
  }
}