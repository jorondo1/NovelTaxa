# Read cntm sum files and add database column
read_label <- function(dir, columns) {
  file_path <- file.path(dir, "sourmash", "cntm_sum.txt")
  read_tsv(file_path, 
           col_names = columns, show_col_types = FALSE) %>% 
    mutate(Dataset = basename(dir))
}
