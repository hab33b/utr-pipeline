load_libraries <- function(){
  # if (!require("BiocManager", quietly = TRUE))
  #   install.packages("BiocManager")
  # BiocManager::install("TCGAbiolinks")
  # BiocManager::install("DESeq2")
  # BiocManager::install("org.Hs.eg.db")
  # BiocManager::install("AnnotationDbi")
  # BiocManager::install("biomaRt")
  library(TCGAbiolinks)
  library(DESeq2)
  library(org.Hs.eg.db) # hg19?
  library(AnnotationDbi) # mapId
  library(biomaRt)
  
  # install.packages("tidyverse")
  # install.packages("gridExtra")
  # install.packages("ggpubr")
  # install.packages("ggrepel")
  library(tidyverse) # ggplot2, dplyr (left_join())
  library(gridExtra) # grid.arrange()
  library(ggpubr) # ggmaplot()
  library(ggrepel)
  options(ggrepel.max.overlaps = Inf)
}

# 1 Download Raw Counts & Filter Data
download_raw_counts <- function(study) {
  
  query <- GDCquery(
    project = paste("TCGA-", toupper(study), sep=""),
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts" #HTSeq - Counts
  )
  GDCdownload(query, directory = "../GDCdata")
  data <- GDCprepare(query, directory = "../GDCdata")
  
  # filter samples - keep normal/tumor samples
  data <- data[, data$shortLetterCode %in% c("NT", "TP")]
  
  # create new column & choose reference
  colData(data)$tumor_type <- factor(data$shortLetterCode, levels = c("NT", "TP"))
  return(data)
}

# 2 Analyze Data
de_analysis <- function(data) {
  # perform DE analysis
  ddsSE <- DESeqDataSet(data, design = ~ tumor_type)
  ddsSE <- ddsSE[rowSums(counts(ddsSE)) >= 10,] # keep || Q: why 10?
  ddsSE <- DESeq(ddsSE)
  
  res <- results(ddsSE, name = resultsNames(ddsSE)[2])
  summary(res)
  return(res)
}

# 3 Edit Data
edit_dea <- function(res, genes) {
  dea <- as.data.frame(res)
  
  # Method 1: AnnotationDbi ensembl -> hgnc ----
  dea <- tibble::rownames_to_column(dea, "geneID")
  dea$hgnc_symbol <- mapIds(org.Hs.eg.db,
                            keys=as.character(sub("[.][0-9]*","",rownames(res))),
                            column="SYMBOL", # what to convert to
                            keytype="ENSEMBL", # what to convert from
                            multiVals="first")
  # ----
  
  # Method 2: biomaRt  (doesn't output the same length vector) ----
  # mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  # ensembls <- sub("[.][0-9]*","",rownames(res)) # converts ENSG00000000003.15 (version) to ENSG00000000003
  # gene_IDs <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","hgnc_symbol"), values=ensembls, mart=mart)
  # dea <- merge(dea, gene_IDs, by.x=0, by.y="ensembl_gene_id", by="ensembl_gene_id", all.x = T)
  # ----
  
  # remove duplicates
  # Q: duplicate what???
  dea <- dea %>% distinct()
  
  # add expression/regulation
  dea$diffexpressed <- "NS"
  dea$diffexpressed[dea$log2FoldChange > 1 & dea$padj < 0.05] <- "Up"
  dea$diffexpressed[dea$log2FoldChange < -1 & dea$padj < 0.05] <- "Down"
  
  # find top 10 most expressed genes
  # Q: how do you define most expressed? order by padj then pick greatest log2foldchange
  dea$genelabels <- ""
  dea[order(dea$padj),]$genelabels[1:10] <- dea[order(dea$padj),]$hgnc_symbol[1:10]
  
  # find genes of interest
  for (gene in genes) {
    dea$genelabels[which(dea == toupper(gene),  arr.ind = TRUE)[1]] <- paste("*", toupper(gene), sep="")
  }
  
  return(dea)
}

# 4 Label 3'UTR and 5'UTR
label_genes <- function(dea) {
  files <- c("mut", "u3", "u5")
  for (file in files) {
    df <- read.csv(paste("../onco/", file, sep=""), sep="\t")
    genes <- drop_na(df)$GENE_ID[drop_na(df)$Q_VALUE < 0.25]
    dea[file] <- ""
    for (gene in genes) {
      if(gene %in% sub("[.][0-9]*","", dea$geneID)){
        dea[file][sub("[.][0-9]*","", dea$geneID) == gene,] <- df$SYMBOL[df$GENE_ID == gene]
      }
    }
  }
  return(dea)
}

# 5 Plotting
plotting <- function(dea) {
  mycolors <- c("gray70", "coral2", "dodgerblue2")
  names(mycolors) <- unique(dea$diffexpressed)
  
  p <- ggplot(dea, aes(x = log2FoldChange, y = -log10(padj))) + 
    # base hgnc data
    geom_point(data=dea, aes(col = diffexpressed), alpha=0.5, size=0.6) +
    geom_label_repel(data=dea, aes(col = diffexpressed, label = genelabels), show.legend = FALSE, alpha=0.5) +
    geom_vline(xintercept=c(-1, 1), lty="dashed", col="gray30") + 
    geom_hline(yintercept=-log10(0.05), lty="dashed", col="gray30") +
    # theme
    scale_color_manual(values = mycolors) +
    theme_classic() + # removes gray background
    labs(title = "Volcano Plot",
         y = expression(-log[10]("padj")),
         col="Regulation") +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    theme(plot.title = element_text(family="Palatino", hjust = 0.5, face = "bold", size=15),
          plot.subtitle = element_text(hjust = 0.5, size=10),
          legend.title = element_text(face = "bold"),
          legend.position = "bottom",
          legend.background = element_rect(fill="gray97")
    )
    
  p_mut <- p + geom_point(data=dea[dea["mut"] != "",],shape=8, alpha=1, size=1.5, col="black", aes(col="black")) +
    geom_label_repel(data=dea[dea["mut"] != "" & dea$diffexpressed != "NS",], col="black", aes(label=mut), show.legend = FALSE) +
    labs(subtitle = "TOP 10 GENES + CODING MUTATIONS")
  p_u3 <- p + geom_point(data=dea[dea["u3"] != "",], shape=8, alpha=1, size=1.5, col="black") +
    geom_label_repel(data=dea[dea["u3"] != "" & dea$diffexpressed != "NS",], col="black", aes(label=u3), show.legend = FALSE) +
    labs(subtitle = "TOP 10 GENES + UTR3 GENES")
  p_u5 <- p + geom_point(data=dea[dea["u5"] != "",], shape=8, alpha=1, size=1.5, col="black") +
    geom_label_repel(data=dea[dea["u5"] != "" & dea$diffexpressed != "NS",], col="black", aes(label=u5), show.legend = FALSE) +
    labs(subtitle = "TOP 10 GENES + UTR5 GENES")

  # plot all ----
  grid.arrange(p_mut, p_u3, p_u5, layout_matrix=cbind(1,2,3), top = paste("TCGA-",toupper(study),sep=""))
  
  # ma plot ----
  # plotMA(dea, cex=.3, main="Plot MA")
  # ma <- ggmaplot(dea, main = "Plot MA", fdr = 0.05, fc = 2,  # Q: false discovery rate?, fold change threshold?
  #          alpha = 0.6, size = 0.6, top = 10, label.rectangle = TRUE,
  #          genenames = dea$genelabels,
  #          palette = mycolors,
  #          xlab = "mean of normalized counts")
  # ----
  
  
}

# Output Data
data_to_excel <- function(dea, study) {
  file <- paste("../excel/tcga_", tolower(study),".csv", sep="")
  write.csv(dea, file)
}

# Specific Questions From Data
question <- function(dea, genes, study) {
  # what is the fold change for mab21L4 and ret
  print(paste("For ", study, ":", sep=""))
  for (gene in genes) {
    print(paste(gene, "has a log2foldchange of", signif(dea$log2FoldChange[which(dea == toupper(gene),  arr.ind = TRUE)[1]], 3)))
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ALGORITHM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load_libraries()
# Global Variables
study <- str_squish(readline("Enter the TCGA code (cesc, hnsc, lusc, skcm): "))

rawdata = download_raw_counts(study)                # 1   Download Counts & Filter Data
dea = edit_dea(de_analysis(rawdata), genes)         # 2&3 Analyze & Edit Data
dea = label_genes(dea)                              # 4   Label 3'UTR and 5'UTR
plotting(dea)                                       # 5   Plotting
# data_to_excel(dea, study)                         # 5   Output Data
# question(dea, genes, study)                       # 6   Specific Questions From Data

# TCGA - sample_type
