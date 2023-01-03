# setup ========================================================================
library(tidyverse)

# function definitions =========================================================
# define function to compute RPKM
rpkm <- function(count_table, widths) {
  require(tidyverse)
  
  row_names <- tibble(id = rownames(count_table))
  
  
  # calculate rpkm
  kb_widths <- widths / 10^3
  
  column_rpkm <- function(x, widths = kb_widths) {
    cr <- (x / widths) / (sum(x) / 10^6 )
    return(cr)
  }
  
  all_rpkm <- count_table %>%
    map(column_rpkm) %>%
    as_tibble() %>%
    bind_cols(row_names) %>%
    column_to_rownames(var = "id")
  
  return(all_rpkm)
  
}

# read in count data ===========================================================
counts <- read_tsv(snakemake@input[[1]], comment = "#") |> 
  select(-c(2:5)) |> 
  column_to_rownames(var = "Geneid")

# perform RPKM normalization ===================================================
rpkm_table <- rpkm(count_table = select(counts, -c("Length")), widths = counts$Length)

# write output table ===========================================================
rpkm_table |> 
  write_tsv(snakemake@output[[1]])