# Read cntm sum files and add database column
read_label <- function(dir, columns) {
  file_path <- file.path(dir, "sourmash", "cntm_sum.txt")
  read_tsv(file_path, 
           col_names = columns, show_col_types = FALSE) %>% 
    mutate(Dataset = basename(dir))
}

# Estimate diversity from abundance matrix, uses phyloseq object as intermediate
div_wrapper <- function(abund.mx, meta.mx) {
  phyloseq(otu_table(abund.mx,taxa_are_rows = TRUE),
           sample_data(meta.mx)) %>% 
    rarefy_even_depth2(verbose = FALSE) %>% 
    estimate_diversity(index = 'Shannon')
}

