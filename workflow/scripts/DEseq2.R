# setup ------------------------------------------------------------------------
## load necessary packages
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))

# import count table -----------------------------------------------------------
count_table <- read_tsv(snakemake@input[[1]]) %>% 
  select(-c(2:6)) %>% 
  column_to_rownames(var = "Geneid")
colnames(count_table) <- gsub(".bam", "", colnames(count_table))

# create colData table ---------------------------------------------------------
coldata <- read_tsv(snakemake@params[["samples"]])

# run DESeq2 -------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(as.matrix(count_table), coldata, design = snakemake@params[["model"]])

# filter out low count genes
keep <- rowSums(counts(dds)) >= snakemake@params[["count_threshold"]]
dds2 <- dds[keep,]

# test for differential gene expression
dds2 <- DESeq(dds2)

# write dds2 object to Rdata file ----------------------------------------------
saveRDS(dds2, file = snakemake@output[[1]])