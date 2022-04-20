# setup ------------------------------------------------------------------------
## load necessary packages
library(DESeq2)
library(tidyverse)

# import dds object ------------------------------------------------------------
## import dds2 object
dds <- readRDS(snakemake@input[[1]])

# get DEseq results for all contrasts ------------------------------------------
contrast <- c("condition", snakemake@params[["contrast"]])
results <- results(dds, contrast = contrast, alpha = snakemake@params[["padj_cutoff"]]) %>% 
  arrange(padj)

# add additional information to results table ----------------------------------

# write output file ------------------------------------------------------------


# add gene_ids and other useful information
#read in GTF file
gtf_fn <- here("/Volumes/TG1/genomics/genomes/fb_dmel-all-r6.26.gtf")
gtf_genes <- import(gtf_fn) %>% 
  mcols() %>% 
  as.data.frame() %>% 
  filter(type == "gene") %>% 
  dplyr::select(gene_id, gene_symbol)

# read in info on nearest gene
tw_bound_genes <- read_tsv(here("data/2022-01_CHIP/results/hc_S2-Tw_aTw_peak_annotations.tsv")) %>% 
  # rename(gene_id = feature, peak_id = name) %>%
  select(seqnames, start, end, geneId, distanceToTSS)

## add useful columns to results table for all genes

for (i in comparisons) {
  ## for each comparison, create a new table containing only significantly differentially expressed genes
  print(i)
  i_name <- paste(i, "_resOrdered", sep = '')
  i_data <- eval(as.symbol(i_name))
  i_data <- as.data.frame(i_data)
  i_data <- tibble::rownames_to_column(i_data, var = "gene_id")
  
  ## add new column to results containing the gene name
  new_table <- dplyr::left_join(i_data, gtf_genes, by = "gene_id")
  
  
  ## add column indicating if gene is differentially expressed with padj < 0.05 and FC > 2
  new_table <- new_table %>%
    dplyr::mutate(is_diff = (padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))) %>%
    replace_na(list(is_diff = FALSE))
  
  #   
  
  ## create a final output table, save it to an object, and write it to a file
  new_name <- paste(i_name, "_all_final.txt", sep = "")
  assign(new_name, new_table)
  fn <- paste("data/2022-01_RNAseq/results/", new_name, sep = "")
  write_tsv(new_table, fn)
  
}


