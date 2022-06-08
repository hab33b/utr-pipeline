# https://www.biostars.org/p/430015/

library("biomaRt")
kidney <- data.frame(gene_id=c("ENSG00000000003", "ENSG00000000003","ENSG00000000005","ENSG00000000419","ENSG00000000457","ENSG00000000460", "ENSG00000116883"), vals=runif(7
                                                                                                                                                                            ))

# Mehod 1: using peptide
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensembls <- sub("[.][0-9]*","",kidney$gene_id) # converts ENSG00000000003.15 (version) to ENSG00000000003
gene_IDss <- getBM(filters="ensembl_gene_id", 
                  attributes=c("ensembl_gene_id","hgnc_symbol"), 
                  values=ensembls, mart=mart)
# kidney$hgnc_symbol <- gene_IDs$hgnc_symbol

kidneynew <- merge(gene_IDss,kidney, by.x="ensembl_gene_id", by.y="gene_id", by="ensembl_gene_id", all.x = T)
# left_join(kidney, gene_IDs, by = c("gene_id"="ensembl_gene_id"))




kidney$hgnc = mapIds(org.Hs.eg.db,
                  keys=as.character(kidney$gene_id), 
                  column="SYMBOL",
                  keytype="ENSEMBL",
                  multiVals="first")
