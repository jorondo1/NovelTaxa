library(pacman)
p_load(tidyverse, phyloseq, magrittr)

source('scripts/comm_comp_data.R')

### Containment across databases and projects
cntm.df %>% # Format db and dataset names
  mutate(Dataset = case_when(Dataset == "Boreal_mosses" ~ "Boreal mosses (4 species)",
                             Dataset == "Provid19" ~ "Human saliva (Sherbrooke)"),
         db = case_when(db == dbNames[1] ~ 'Default GTDB r214',
                        db == dbNames[2] ~ 'r214 with nMAGs')) %>% 
  # Plot !
  ggplot(aes(x = db, y = cntm)) +
  geom_violin(aes(fill=db)) + # so dots don't appear in the legend!
  geom_jitter(alpha = 0.5, width=0.1,size = 2) + 
  facet_wrap(~Dataset, scales= 'free', nrow=2) +
  theme_light(base_size = 24) + 
  labs(y = 'Sample k-mer containment', x = '', fill = '') +
  theme(axis.text.x = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.margin = margin(t = -30),
        axis.text.y = element_text(size = 12)) 

ggsave('figures/cntm.pdf',  bg = 'white', 
       width = 20, height = 40, units = 'cm')

### Containment increase from adding nMAGs (Mosses only)
cntmDiff %>% filter(Dataset == 'Boreal_mosses') %>% 
  ggplot(aes(x = "Boreal moss", y = fold_inc)) +
  geom_violin() + geom_jitter(alpha =0.5, width = 0.2) +
  labs(y = 'Fold increase in sample k-mers containment', x = 'Microbiome') +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, max(cntmDiff$fold_inc), by = 2))  # Set minimum y value and increments

### Sample-wise change in diversity 
ggplot(moss_div_long, aes(x = db, y = Shannon, group = Sample)) +
  geom_point() + geom_line() + theme_minimal()

### Diversity by host across database 
moss_div_long %>% # Format db names
  mutate(db = factor(db, levels = dbNames),
         db = case_when(db == dbNames[1] ~ 'Default GTDB r214',
                        db == dbNames[2] ~ 'r214 with nMAGs')) %>% 
  # Plot !
  ggplot(aes(x = db, y = Shannon, fill = Host)) +
    geom_violin() + theme_light() +
    facet_wrap(~db, scales = "free", nrow = 2) + 
    labs(y = 'Shannon diversity') +
    theme_light(base_size = 24) + 
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank()) +
    scale_fill_manual(values = MetBrewer::met.brewer('VanGogh2', n=4, direction = 1),
                      labels = c("Dicranum undulatum", "Polytrichum commune", "P. juniperinum", "P. piliferum")) +
    theme(axis.text.x = element_blank(),
          legend.position = 'bottom',
          legend.title = element_blank(),
          legend.margin = margin(t = -10),
          axis.text.y = element_text(size = 12),
          legend.key.height = unit(1.2, 'cm')) +
    guides(fill = guide_legend(nrow = 2))

ggsave('figures/div.pdf',  bg = 'white', 
       width = 20, height = 40, units = 'cm')
