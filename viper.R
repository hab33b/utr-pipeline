# viper.R
# ***************************
# input: deseq object & maf file
# 

BiocManager::install(c("mixtools","bcellViper","viper"))
library(DESeq2)
library(viper)
library(tidyverse) # ggplot2, dplyr left_join() str_squish()
library(viper)
library(mixtools)
library(bcellViper)
library(aracne.networks)
library(biomaRt)

# test if it's worth doing recurrent mutations
df <- read_tsv("../data/onco/HNSC-3'UTR")
df2 <- df %>% dplyr::count(CHROMOSOME, POSITION, REF, ALT) %>%
  arrange(desc(n))
hist(df2$n)

studies <- c("cesc", "hnsc", "lusc")
study <- "cesc"

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

ensgene <- as.character(sub("[.][0-9]*","",rownames(vsd)))
rownames(normalMat) <- ensgene
rownames(tumorMat) <- ensgene

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
  ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")  
  genes <- getBM(filters="entrezgene_id", attributes=c("ensembl_gene_id","entrezgene_id"), values=entrezIDs, mart=ensembl)
  genes <- as_tibble(genes)
  
  # subset and rename ensemble genes in count matrices
  keepGenes <- genes[genes$ensembl_gene_id %in% rownames(normalMat), ] # subset gene dict
  normalMat <- normalMat[keepGenes$ensembl_gene_id, ] # subset count matrix
  tumorMat <- tumorMat[keepGenes$ensembl_gene_id, ]
  rownames(normalMat) <- keepGenes$entrezgene_id # convert names
  rownames(tumorMat) <- keepGenes$entrezgene_id
  
  # update adjfile to fit 3col format
  b <- a %>%
    dplyr::select(-likelihood) %>%
    write.table(adjfile, quote=F, sep="\t", row.names=F, col.names=F)
  
  # create regulon
  regul <- aracne2regulon(adjfile, eset=cbind(normalMat, tumorMat), format="3col")
} else {
  print("Loading from files")
  
  # upload regulon
  regulon <- readRDS(paste("../data/rds/VIPER_regulon_", toupper(study), ".rds", sep=""))
}

# viper network
vpsig <- viperSignature(tumorMat, normalMat, cores = 2)
vpres <- viper(vpsig, regulon, verbose = TRUE)

maf <- read_tsv("../data/mafs/TCGA-CESC.maf")

# for each gene in each sample label as mutant or wt
x <- maf %>%
  mutate(Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode),
         Hugo_Symbol = factor(Hugo_Symbol)) %>%
  filter(Variant_Classification %in% c("3'UTR", "5'UTR")) %>%
  # group_by(Tumor_Sample_Barcode) %>%
  dplyr::count(Tumor_Sample_Barcode, Hugo_Symbol, .drop = F) %>%
  mutate(is_mutated = ifelse(n > 0, "Mutated", "WT"))

hgncs <- unique(c(x$Hugo_Symbol))
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
genes <- getBM(filters="hgnc_symbol", attributes=c("hgnc_symbol","entrezgene_id"), values=hgncs, mart=ensembl)
genes2 <- getBM(filters="hgnc_symbol", attributes=c("hgnc_symbol","ensembl_gene_id"), values=hgncs, mart=ensembl)

vpres <- tibble::rownames_to_column(as.data.frame(vpres))

# subset and rename ensemble genes in count matrices
keepGenes <- genes2[genes2$entrezgene_id %in% vpres$rowname, ]
vpres <- vpres[keepGenes$ensembl_gene_id, ]
rownames(vpres) <- keepGenes$entrezgene_id

# if vpres is ensemble genes
keepGenes2 <- genes2[genes2$ensembl_gene_id %in% rownames(vpres), ]
vpres <- vpres[keepGenes2$ensembl_gene_id, ]

# create table of significant genes
genes2$p_value <- ""
for (r in rownames(vpres)) {
  hgnc <- genes2 %>% filter(ensembl_gene_id == r) %>% pull(hgnc_symbol)
  if (hgnc %in% genes2$hgnc_symbol){
    row <- vpres[rownames(vpres)==r,]
    mut_samples <- x %>%
      filter(Hugo_Symbol == hgnc, is_mutated == "Mutated") %>%
      pull(Tumor_Sample_Barcode)
    wt_samples <- x %>%
      filter(Hugo_Symbol == hgnc, is_mutated == "WT") %>%
      pull(Tumor_Sample_Barcode)
    z = t.test(row[as.character(mut_samples)], row[as.character(wt_samples)])
    genes2$p_value[genes2$ensembl_gene_id == r] <- z$p.value
  }
}


