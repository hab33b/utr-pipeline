# onco.R
# ************************
# input: tcga data
# filters by variant type (3'UTR & 5'UTR)
# outputs files formatted for submission into oncodrivefml
# http://bbglab.irbbarcelona.org/oncodrivefml/home
library(tidyverse)


studies <- c("cesc", "hnsc", "lusc")
coding <- c("Splice_Site", "Nonsense_Mutation") #,"Missense_Mutation", "Silent")
variants <- c(list("3'UTR", "5'UTR", coding))

for (study in studies) {
  for (variant in variants) {
    maf <- read_tsv(paste("../data/mafs/TCGA-", toupper(study), ".maf", sep=""))
    maf %>% 
      filter(Variant_Classification %in% variant) %>%
      select(Chromosome, Start_Position, Reference_Allele, 
             Tumor_Seq_Allele2, Tumor_Sample_Barcode) %>%
      rename(chromosome = Chromosome, position = Start_Position,
             ref = Reference_Allele, alt = Tumor_Seq_Allele2,
             sample = Tumor_Sample_Barcode) %>%
      rename_with(toupper) %>%
      write.table(paste("../data/mafs/TCGA-", toupper(study), 
                        "-" ,toupper(variant[1]),".maf", sep=""),
                  quote=F, sep="\t", row.names=F)
  }
}
