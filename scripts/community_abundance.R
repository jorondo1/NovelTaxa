library(tidyverse, phyloseq)

# Sourmash output parsing function (custom)
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/community_functions.R'))
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/rarefy_even_depth2.R'))

# Import environment arguments
args <- commandArgs(trailingOnly = TRUE)
dbName <- args[1] #e.g. r214

# Sample k-mer containment across databases 
cntm.df <- read_delim("sourmash/cntm_sum.txt", 
                      col_names = c("Sample", "db", "cntm"))

cntm.df %>% group_by(db) %>% 
  summarise(mean_cntm = mean(cntm),
            sd_cntm = sd(cntm),
            n_sam = n())

ggplot(cntm.df, aes(x = db, y = cntm, fill=db)) +
    geom_violin() + geom_jitter(alpha = 0.2) + 
    theme_minimal(base_size = 16) + guides(fill = 'none') +
    labs(y = 'Sample k-mer containment', x = 'Reference database')

# Containment increase after adding MAGs 
cntmDiff <- cntm.df %>% 
  pivot_wider(names_from = db, values_from = cntm) %>% 
  mutate(fold_inc = custom/!!sym(dbName))

cntmDiff %>% summarise(
  mean_inc = mean(fold_inc),
  sd_inc = sd(fold_inc)
)

ggplot(cntmDiff, aes(x = "Boreal moss", y = fold_inc)) +
  geom_violin() + geom_jitter(alpha =0.2, width = 0.2) +
  theme_minimal(base_size = 16) + 
  labs(y = 'Fold increase in sample k-mers containment', x = 'Microbiome') +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, max(cntmDiff$fold_inc), by = 2))  # Set minimum y value and increments

# Alpha diversity #

# Parse Sourmash gather output 
custom <- parse_SM('sourmash/*custom_gather.csv')
GTDB <- parse_SM(paste0('sourmash/*',dbName,'_gather.csv'))

# Import metadata from moss project 
moss.ps <- readRDS(url("https://github.com/jorondo1/borealMoss/raw/main/data/R_out/mossMAGs.RDS",
           method="libcurl"))
sample_data <- moss_custom.ps@sam_data %>% data.frame %>% rownames_to_column("Sample")

# Phyloseq objects 
moss_custom.ps <- phyloseq(otu_table(custom, taxa_are_rows = TRUE), 
                           sample_data(moss.ps@sam_data)) 

moss_GTDB.ps <- phyloseq(otu_table(GTDB, taxa_are_rows = TRUE), 
                          sample_data(moss.ps@sam_data)) 

# compute alpha-diversity 
moss_div = data.frame(
  GTDB = rarefy_even_depth2(moss_GTDB.ps) %>% 
    estimate_diversity(index = 'Shannon'),
  custom = rarefy_even_depth2(moss_custom.ps) %>% 
    estimate_diversity(index = 'Shannon')
) %>% rownames_to_column("Sample")

moss_div_long <- moss_div %>% 
  pivot_longer(
    values_to = 'Shannon',
    names_to = 'Database',
    cols = c(dbName, "custom")) %>% 
  left_join(dplyr::select(sample_data, Sample, Compartment, Host), 
            by = "Sample")

moss_div_long %>% group_by(Database) %>% 
  summarise(meanDiv = mean(Shannon),
            sdDiv = sd(Shannon))

ggplot(moss_div_long, aes(x = Database, y = Shannon, group = Sample)) +
  geom_point() + geom_line() + theme_minimal()

# Check how it affects a simple diversity analysis

ggplot(moss_div_long, aes(x = Database, y = Shannon, fill = Host)) +
  facet_grid(~Database, scales = "free") + labs(x="") +
  geom_violin() + theme_light()

lm_GTDB <- moss_div_long %>% 
  filter(Database == dbName) %>% 
  lm(Shannon ~ Host, data = .) 

lm_custom <- moss_div_long %>% 
  filter(Database == 'custom') %>% 
  lm(Shannon ~ Host, data = .) 

aov(lm_GTDB) %>% summary # barely significant difference
aov(lm_custom) %>% summary # highly significant difference
summary(lm_GTDB) %$% r.squared # 6% variance explained
summary(lm_custom) %$% r.squared # 20% variance explained

residuals(lm_custom) %>% shapiro.test # normal-ish p=0.04
residuals(lm_GTDB) %>% shapiro.test # normal p=0.50

# non-parametric supports
moss_div_long %>% 
  filter(Database == 'custom') %>% 
  kruskal.test(Shannon ~ Host, data = .)

moss_div_long %>% 
  filter(Database == dbName) %>% 
  kruskal.test(Shannon ~ Host, data = .)
