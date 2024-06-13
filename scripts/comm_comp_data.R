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
#dbNames <- cntm.df %$% db %>% unique %>% rev # extract db names
dbNames <- c('r214', 'nMAGs', 'nMAGsBetter')
cntm.df %<>% mutate(db = factor(db, levels = dbNames))# Reorder factors 

message('Summary of containment across databases')
cntm.df %>% group_by(Dataset, db) %>% 
  summarise(mean_cntm = mean(cntm),
            sd_cntm = sd(cntm),
            .groups='keep') %>% as.data.frame %>% print

# Containment increase after adding MAGs
cntmDiff <- cntm.df %>% 
  pivot_wider(names_from = db, values_from = cntm) %>% 
  mutate(fold_inc_12 = !!sym(dbNames[2])/!!sym(dbNames[1]),
         diff_12 = !!sym(dbNames[2]) - !!sym(dbNames[1]),
         fold_inc_23 = !!sym(dbNames[3])/!!sym(dbNames[2]),
         diff_23 = !!sym(dbNames[3]) - !!sym(dbNames[2]))

message('Mean fold increase vs. default')
cntmDiff %>% group_by(Dataset) %>% # Summarise increase :
  summarise(mean_inc_12 = mean(fold_inc_12),
            sd_inc_12 = sd(fold_inc_12),
            mean_inc_23 = mean(fold_inc_23),
            sd_inc_23 = sd(fold_inc_23)) %>% as.data.frame %>% print

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

### COHEN'S D
message("Cohen's d - difference in terms of standard deviations from the mean")
cohenD <- function(df_long, ds) {
  meanSD <- df_long %>% filter(Dataset == ds) %>% 
    group_by(db) %>% 
    summarise(mean = mean(cntm), var = var(cntm), n=n())
  if (meanSD$n[1] != meanSD$n[2]) {message('Warning! Datasets have different number of obs.')}
  n <- meanSD$n[1] # both datasets should have same n !
  pooledSD <- sqrt(((n-1) * (meanSD$var[1]) + (n-1) * (meanSD$var[3]))/(2*n-2))
  (meanSD$mean[3] - meanSD$mean[1])/pooledSD # Cohen's d
}

message(paste("Boreal mosses: ",cohenD(cntm.df, 'Boreal_mosses')))
message(paste("Saliva: ", cohenD(cntm.df, 'Saliva')))

### PAIRED COHEN'S D
cohenD_p <- function(df_wide, ds) {
  df <- df_wide %>% filter(Dataset == ds)
  d_12 <- mean(df$diff_12)/sd(df$diff_12)
  d_23 <- mean(df$diff_23)/sd(df$diff_23)
  message(paste("Cohen's paired d between scenario 1 and 2: ",format(d_12, digits=3)))
  message(paste("Cohen's paired d between scenario 2 and 3: ",format(d_23, digits=3)))
}

cohenD_p(cntmDiff, 'Boreal_mosses')
cohenD_p(cntmDiff, 'Saliva')

###################
# Alpha diversity #
###################

# Import moss project metadata
moss_sam <- readRDS(url("https://github.com/jorondo1/borealMoss/raw/main/data/R_out/mossMAGs.RDS",
                       method="libcurl")) %>% sample_data %>% 
  data.frame(Microbiome='Boreal_mosses') %>% rownames_to_column('Sample') %>% 
  dplyr::select(Sample, Compartment, Host, Location) 

# Parse Sourmash gather output 
SM_out <- list()
for (p in project_paths) {
  for (db in dbNames) {
    SM_out[[basename(p)]][[db]] <- parse_SM(paste0(p, "/sourmash/", "*", db, "_gather.csv"))
  }
}

# Dummy saliva metadata for phyloseq functions
saliva_sam <- tibble(
  Sample=colnames(SM_out$Saliva$r214),
  V1=1) 
  
# Create df with diversity for each dataset
Boreal_moss_div <- lapply(SM_out$Boreal_mosses, div_wrapper, 
       meta.mx = moss_sam %>% column_to_rownames('Sample')) %>% 
  do.call(cbind, .) %>% as.data.frame %>% 
  rownames_to_column('Sample') %>% 
  mutate(Microbiome = 'Boreal_mosses')

Saliva_div <- lapply(SM_out$Saliva, div_wrapper, 
                     meta.mx = saliva_sam %>% column_to_rownames('Sample') %>% as.data.frame) %>% 
  do.call(cbind, .)%>% as.data.frame %>% 
  rownames_to_column('Sample') %>% 
  mutate(Microbiome = 'Saliva')

# lapply(names(SM_out), function(Microbiome) {
#   div <- lapply(SM_out[[Microbiome]])
# })

# # pull moss sample data
# sample_data <- moss_custom.ps@sam_data %>% data.frame %>% rownames_to_column("Sample")

# Long version of df for stat tests, add metadata :
div_long <- rbind(Saliva_div, Boreal_moss_div) %>% 
  pivot_longer(
    values_to = 'Shannon',
    names_to = 'db',
    cols = all_of(dbNames)) %>% 
  #MOss data
  left_join(moss_sam, by = "Sample") %>% 
  left_join(saliva_sam, by = "Sample")

message('Mean and variance of diversity across db:')
div_long %>% group_by(Microbiome, db) %>% # Summarise diversity 
  summarise(meanDiv = mean(Shannon),
            varDiv = var(Shannon)) %>% as.data.frame %>% print

# Linear regression Shannon - Host moss species
lm_r214 <- div_long %>% 
  filter(db == dbNames[1] & Microbiome == 'Boreal_mosses') %>% 
  lm(Shannon ~ Host, data = .) 

lm_nMAGs <- div_long %>% 
  filter(db == dbNames[2] & Microbiome == 'Boreal_mosses') %>% 
  lm(Shannon ~ Host, data = .) 

lm_nMAGsBetter <- div_long %>% 
  filter(db == dbNames[3] & Microbiome == 'Boreal_mosses') %>% 
  lm(Shannon ~ Host, data = .) 

summary(lm_r214) %>%  print # 6% variance explained
aov(lm_r214) %>% summary # barely significant difference
summary(lm_nMAGs) %>% print # 20% variance explained
aov(lm_nMAGs) %>% summary # highly significant difference
summary(lm_nMAGsBetter) %>% print # 20% variance explained
aov(lm_nMAGsBetter) %>% summary # highly significant difference

residuals(lm_nMAGs) %>% shapiro.test # normal-ish p=0.04
residuals(lm_r214) %>% shapiro.test # normal p=0.50

# non-parametric supports
div_long %>% 
  filter(db == dbNames[1]) %>% 
  kruskal.test(Shannon ~ Host, data = .)

div_long %>% 
  filter(db == dbNames[3]) %>% 
  kruskal.test(Shannon ~ Host, data = .)
