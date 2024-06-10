library(pacman)
p_load(tidyverse, phyloseq, magrittr)

# Sourmash output parsing function (custom)
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/community_functions.R'))
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/rarefy_even_depth2.R'))
source("scripts/myFunctions.R")

# Import environment arguments
# args <- commandArgs(trailingOnly = TRUE)
# db <- args[1] # db <- 'r214'

# List paths for each project
project_paths <- list.dirs(path='.', recursive=FALSE, full.names=TRUE) %>% 
  # extract 1st level dir name where */sourmash/cntm_sum.txt exists
  purrr::keep(~ file.exists(file.path(.x, "sourmash", "cntm_sum.txt")))

###############################################
### Sample k-mer containment across databases #
###############################################
cntm.df <- map_dfr(project_paths, 
                   ~ read_label(.x , columns = c('Sample', 'db', 'cntm'))) 
dbNames <- cntm.df %$% db %>% unique %>% rev # extract db names
cntm.df %<>% mutate(db = factor(db, levels = dbNames))# Reorder factors 

message('Summary of containment across databases')
cntm.df %>% group_by(Dataset, db) %>% 
  summarise(mean_cntm = mean(cntm),
            sd_cntm = sd(cntm),
            .groups='keep') %>% as.data.frame %>% print

# Containment increase after adding MAGs
cntmDiff <- cntm.df %>% 
  pivot_wider(names_from = db, values_from = cntm) %>% 
  mutate(fold_inc = !!sym(dbNames[2])/!!sym(dbNames[1]))

message('Mean fold increase by project')
cntmDiff %>% group_by(Dataset) %>% # Summarise increase :
  summarise(mean_inc = mean(fold_inc),
            sd_inc = sd(fold_inc)) %>% as.data.frame %>% print

# Confidence interval in mean cntm increase
# within which range do 90 % of sample increases fall with 95% confidence ?
library(boot)

quantile_function <- function(data, indices) {
  data[indices] %>% quantile(probs = c(0.05, 0.95))
}
message('90% of samples fold-increases are within this range (95% confidence):')
# bootstrap
finc <- cntmDiff %>% filter(Dataset == 'Boreal_mosses') %$% fold_inc
bootstrap_results <- boot(finc, statistic = quantile_function, R = 1000)
print(bootstrap_results$t0) # between 2.1 and 9.1 fold increase
#boot.ci(bootstrap_results, type = "perc", index = 1) # CI 95% (default)
#boot.ci(bootstrap_results, type = "perc", index = 2)

###################
# Alpha diversity #
###################

# Parse Sourmash gather output 
SM_out <- list()
for (p in project_paths) {
  for (db in dbNames) {
    name <- basename(p) %>% paste0('_', db)
    SM_out[[name]] <- parse_SM(paste0(p, "/sourmash/", "*", db, "_gather.csv"))
  }
}

# Import moss project metadata
moss.ps <- readRDS(url("https://github.com/jorondo1/borealMoss/raw/main/data/R_out/mossMAGs.RDS",
           method="libcurl"))

# Phyloseq objects 
moss_custom.ps <- phyloseq(otu_table(SM_out$Boreal_mosses_custom, taxa_are_rows = TRUE), 
                           sample_data(moss.ps@sam_data)) 

moss_r214.ps <- phyloseq(otu_table(SM_out$Boreal_mosses_r214, taxa_are_rows = TRUE), 
                          sample_data(moss.ps@sam_data)) 

# Rarefy + compute alpha-diversity for both datasets
moss_div = data.frame(
  V1 = rarefy_even_depth2(moss_r214.ps, verbose = FALSE) %>% # custom fct
    estimate_diversity(index = 'Shannon'), # custom fct
  V2 = rarefy_even_depth2(moss_custom.ps, verbose = FALSE) %>% 
    estimate_diversity(index = 'Shannon')
  ) %>% setNames(dbNames) %>% # var names
  rownames_to_column("Sample")

# pull sample data
sample_data <- moss_custom.ps@sam_data %>% data.frame %>% rownames_to_column("Sample")
# Long version of df for stat tests, add metadata :
moss_div_long <- moss_div %>% 
  pivot_longer(
    values_to = 'Shannon',
    names_to = 'db',
    cols = c(dbNames[2], dbNames[1])) %>% 
  left_join(dplyr::select(sample_data, Sample, Compartment, Host, Location), 
            by = "Sample") 

message('Mean and variance of diversity across db:')
moss_div_long %>% group_by(db) %>% # Summarise diversity 
  summarise(meanDiv = mean(Shannon),
            varDiv = var(Shannon)) %>% as.data.frame %>% print

moss_div %>% mutate(delta_div = (custom-r214)/r214) %>% 
  summarise(mean = mean(delta_div),
            sd = sd(delta_div))

# Linear regression Shannon - Host moss species
lm_r214 <- moss_div_long %>% 
  filter(db == dbNames[1]) %>% 
  lm(Shannon ~ Host, data = .) 

lm_custom <- moss_div_long %>% 
  filter(db == dbNames[2]) %>% 
  lm(Shannon ~ Host, data = .) 

summary(lm_r214) %>%  print # 6% variance explained
aov(lm_r214) %>% summary # barely significant difference
summary(lm_custom) %>% print # 20% variance explained
aov(lm_custom) %>% summary # highly significant difference

residuals(lm_custom) %>% shapiro.test # normal-ish p=0.04
residuals(lm_r214) %>% shapiro.test # normal p=0.50

# non-parametric supports
moss_div_long %>% 
  filter(db == dbNames[1]) %>% 
  kruskal.test(Shannon ~ Host, data = .)

moss_div_long %>% 
  filter(db == dbNames[2]) %>% 
  kruskal.test(Shannon ~ Host, data = .)
