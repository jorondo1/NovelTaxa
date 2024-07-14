# Read cntm sum files and add database column
read_label <- function(dir, columns) {
  file_path <- file.path(dir, "sourmash", "cntm_sum.txt")
  read_tsv(file_path, 
           col_names = columns, show_col_types = FALSE) %>% 
    mutate(Dataset = basename(dir))
}

# Estimate diversity from abundance matrix, uses phyloseq object as intermediate
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/community_functions.R'))
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/rarefy_even_depth2.R'))

div_wrapper <- function(input, meta.mx = FALSE) {
  if(is(input, "phyloseq")) {
    ps <- input
  } else {
  ps <- phyloseq(otu_table(input,taxa_are_rows = TRUE),
           sample_data(meta.mx)) 
  }
   ps %>% rarefy_even_depth2(verbose = FALSE) %>% 
    estimate_diversity(index = 'Shannon')
}

