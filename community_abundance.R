library(tidyverse)

# Sourmash output parsing function (custom)
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/community_functions.sh'))

custom <- parse_SM('data/sourmash/*custom_gather.csv')
rs214 <- parse_SM('data/sourmash/*rs214_gather.csv')

# Sample k-mer containment across databases 
cntm.df <- read_delim("data/sourmash/cntm_sum.txt", 
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
  mutate(fold_inc = custom/rs214)

cntmDiff %>% summarise(
  mean_inc = mean(fold_inc),
  sd_inc = sd(fold_inc)
)

ggplot(cntmDiff, aes(x = "Boreal moss", y = fold_inc)) +
  geom_violin() + geom_jitter(alpha =0.2, width = 0.2) +
  theme_minimal(base_size = 16) + 
  labs(y = 'Fold increase in sample k-mers containment', x = 'Microbiome') +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, max(cntmDiff$fold_inc), by = 2))  # Set minimum y value and increments

