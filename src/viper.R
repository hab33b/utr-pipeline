# viper.R
# ***************************
# input: deseq object & maf file
# 
library(DESeq2)
library(viper)
library(tidyverse) # ggplot2, dplyr left_join() str_squish()
library(mixtools)
library(aracne.networks)
library(biomaRt)

# # test possibility of considering recurrent mutations
df <- read_tsv("../data/mafs/TCGA-HNSC-5'UTR.maf")
CESC <- df %>% dplyr::count(Chromosome, Start_Position, Tumor_Seq_Allele1, Tumor_Seq_Allele2) %>%
  arrange(desc(n))
hist(CESC$n)

study <- "hnsc"

# get VST normalized data from DESeq2
ddsSE <- readRDS(paste("../data/rds/DESeq_", toupper(study), ".rds", sep=""))
vsd <- vst(ddsSE) # apply variance stabilizing transformation on count data

# filter counts for each gene by normal type and tumor patients
ntIdx <- which(colData(vsd)$shortLetterCode == "NT")
tpIdx <- which(colData(vsd)$shortLetterCode == "TP")
# sccIdx <- which(vsd$condition == "SCC")
# nsIdx <- which(vsd$condition == "NS")
normalMat <- assay(vsd)[, ntIdx] # matrix of transformed values
tumorMat <- assay(vsd)[, tpIdx]

# update rows and columns of count matrices
ensgene <- as.character(sub("[.][0-9]*","",rownames(vsd)))
rownames(normalMat) <- ensgene
rownames(tumorMat) <- ensgene
colnames(normalMat) <- as.character(substr(colnames(normalMat), 1, 15))
colnames(tumorMat) <- as.character(substr(colnames(tumorMat), 1, 15))

# VIPER
# ****************
# generate regulon
network <- paste("regulon", tolower(study), sep="")
if (network %in% data(package="aracne.networks")$results[, "Item"]) {
  print("Loading from aracne.networks")
  
  adjfile <- paste("../data/network-", tolower(study),".txt", sep="")
  write.regulon(get(network), file=adjfile)
  
  a <- read_tsv(adjfile)
  
  # match ensemble to entrezID
  entrezIDs <- unique(c(a$Regulator, a$Target))
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")  
  genes <- getBM(filters="entrezgene_id", attributes=c("ensembl_gene_id","entrezgene_id"), values=entrezIDs, mart=ensembl)
  genes <- as_tibble(genes)
  
  # subset and rename ensemble genes in count matrices
  keepGenes <- genes[genes$ensembl_gene_id %in% rownames(normalMat), ] # subset gene dict
  normalMat <- normalMat[keepGenes$ensembl_gene_id, ] # subset count matrix
  tumorMat <- tumorMat[keepGenes$ensembl_gene_id, ]
  rownames(normalMat) <- keepGenes$entrezgene_id # convert names
  rownames(tumorMat) <- keepGenes$entrezgene_id
  
} else {
  print("Loading from files")
  # upload regulon
  regulon <- readRDS(paste("../data/rds/VIPER_regulon_", toupper(study), ".rds", sep=""))
}

# viper network
vpsig <- viperSignature(tumorMat, normalMat, cores = 2)
vpres <- viper(vpsig, get(network), verbose = TRUE)

# for each gene in each sample label as mutant or wt
maf <- read_tsv(paste("../data/mafs/TCGA-", (toupper(study)),".maf", sep=""))
x <- maf %>%
  mutate(Tumor_Sample_Barcode = as.character(substr(Tumor_Sample_Barcode, 1, 15))) %>%
  mutate(Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode),
         Hugo_Symbol = factor(Hugo_Symbol)) %>%
  filter(Variant_Classification %in% c("3'UTR", "5'UTR")) %>%
  dplyr::count(Tumor_Sample_Barcode, Hugo_Symbol, .drop = F) %>%
  mutate(is_mutated = ifelse(n > 0, "Mutated", "WT"))

hgncs <- unique(c(x$Hugo_Symbol))
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
genes2 <- getBM(filters="hgnc_symbol", attributes=c("hgnc_symbol","entrezgene_id"), values=hgncs, mart=ensembl)

keepGenes <- genes2[genes2$entrezgene_id %in% rownames(vpres), ]

# create table of significant genes
genes2$p_value <- ""
for (r in rownames(vpres)) {
  hgnc <- genes2 %>% filter(entrezgene_id == r) %>% pull(hgnc_symbol)
  if (length(hgnc %in% genes2$hgnc_symbol) != 0){
    row <- vpres[rownames(vpres)==r,]
    mut_samples <- x %>%
      filter(Hugo_Symbol == hgnc, is_mutated == "Mutated") %>%
      pull(Tumor_Sample_Barcode)
    wt_samples <- x %>%
      filter(Hugo_Symbol == hgnc, is_mutated == "WT") %>%
      pull(Tumor_Sample_Barcode)
    if (length(mut_samples) > 2 & length(wt_samples) > 0) {
      z = t.test(row[as.character(mut_samples)], row[as.character(wt_samples)])
      genes2$p_value[genes2$entrezgene_id == r] <- z$p.value
    }
  }
}
genes2$p_adj <- ""
genes2$p_adj <- p.adjust(genes2$p_value, method="fdr")




