load_libraries <- function(){
  # if (!require("BiocManager", quietly = TRUE))
  #   install.packages("BiocManager")
  # BiocManager::install("TCGAbiolinks")
  # BiocManager::install("DESeq2")
  library(TCGAbiolinks)
  library(DESeq2)

  # install.packages("tidyverse")
  # library(tidyverse) # ggplot2, dplyr (left_join())
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
  
  # output DESeq object
  saveRDS(ddsSE, paste("../data/rds/DESeq_", toupper(study), ".rds", sep=""))
}

load_libraries()
study <- str_squish(readline("Enter the TCGA code (cesc, hnsc, lusc, skcm): "))

rawdata = download_raw_counts(study)      # 1 Download Counts & Filter Data
de_analysis(rawdata)                      # 2 Analyze Data
