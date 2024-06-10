library(pacman)
p_load(tidyverse, phyloseq, magrittr, foreach, iterators, parallel,
       extrafont, showtext)

source('scripts/comm_comp_data.R')
custom_theme <- theme_light(base_size = 28, 
                            base_family = 'Georgia')

### Containment across databases and projects
cntm.df %>% # Format db and dataset names
  mutate(Dataset = case_when(Dataset == "Boreal_mosses" ~ "Boreal mosses",
                             Dataset == "Provid19" ~ "Human saliva"),
         db = case_when(db == dbNames[1] ~ 'Default GTDB r214',
                        db == dbNames[2] ~ 'r214 with nMAGs')) %>% 
  # Plot !
  ggplot(aes(x = db, y = cntm)) +
  geom_violin(aes(fill=db)) + # so dots don't appear in the legend!
  geom_jitter(alpha = 0.5, width=0.1,size = 2) + 
  facet_wrap(~Dataset, scales= 'free', nrow=2) +
  custom_theme +
  labs(y = 'Sample k-mer containment', x = '', fill = '') +
  theme(axis.text.x = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.margin = margin(t = -30),
        axis.text.y = element_text(size = 12)) 

ggsave('figures/cntm.pdf',  bg = 'white', 
       width = 20, height = 45, units = 'cm')

### Containment increase from adding nMAGs (Mosses only)
cntmDiff %>% filter(Dataset == 'Boreal_mosses') %>% 
  ggplot(aes(x = "Boreal moss", y = fold_inc)) +
  geom_violin() + geom_jitter(alpha =0.5, width = 0.2) +
  labs(y = 'Fold increase in sample k-mers containment', x = 'Microbiome') +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, max(cntmDiff$fold_inc), by = 2))  # Set minimum y value and increments


### Diversity by host across database 
moss_div_long_formatted <- moss_div_long %>% # Format db names
  mutate(db = factor(db, levels = dbNames),
         db = case_when(db == dbNames[1] ~ 'Default GTDB r214',
                        db == dbNames[2] ~ 'r214 with nMAGs')) 

### Sample-wise change in diversity 
moss_div_long_formatted %>% 
  mutate(Title = 'Changes in diversity (Boreal mosses)') %>% 
  ggplot(aes(x = db, y = Shannon)) +
    geom_boxplot(aes(fill = db)) +
    geom_point(alpha = 0.5) + geom_line(aes(group = Sample), alpha=0.3) + 
  custom_theme +
  facet_wrap(~Title) +
  labs(y = 'Shannon diversity') +
  #scale_fill_manual(values = MetBrewer::MetPalettes$VanGogh2[[1]][c(5,6)]) +
  #scale_fill_manual(values = MetBrewer::met.brewer('Klimt', n=2, override.order = TRUE)) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.margin = margin(t = -10),
        axis.text.y = element_text(size = 12)) 

ggsave('figures/div_var.pdf',  bg = 'white', 
       width = 20, height = 45, units = 'cm')

# Plot !
moss_div_long_formatted %>% 
  ggplot(aes(x = db, y = Shannon, fill = Host)) +
    geom_boxplot() + theme_light() +
    facet_wrap(~db, scales = "free", nrow = 2) + 
    labs(y = 'Shannon diversity') +
    custom_theme +
    scale_fill_manual(values = MetBrewer::met.brewer('VanGogh2', n=5, direction = 1),
                      labels = c("Dicranum undulatum", "Polytrichum commune", "P. juniperinum", "P. piliferum")) +
    theme(axis.text.x = element_blank(),
          legend.position = 'bottom',
          legend.title = element_blank(),
          legend.margin = margin(t = -10),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          legend.key.height = unit(1.2, 'cm')) +
    guides(fill = guide_legend(nrow = 2))

ggsave('figures/div.pdf',  bg = 'white', 
       width = 20, height = 45, units = 'cm')


# BIOME REPRESENTATION PLOT
genomes <- read_csv('figures/GenomeCatalogue.csv') %>% 
  select(name, unclustered_genome_count) %>% 
  mutate(count = unclustered_genome_count/sum(unclustered_genome_count),
         .keep='unused',
         name = case_when(name == 'Unified Human Gastrointestinal Genome (UHGG) v2.0.2' ~ 'Unified Human Gastrointestinal \nGenome (UHGG) v2.0.2', TRUE ~ name))

genomes %>% 
  ggplot(aes(x = "", y = count, fill = name)) +
  geom_bar(width = 1, stat = 'identity') +
  coord_polar("y", start = 0) +
  theme_void(base_size = 28, 
             base_family = 'Georgia') + 
  theme(legend.position = 'right',
        #   legend.spacing.y = unit(10, 'cm'),
        legend.key.height = unit(1.4, 'cm')
  ) +
  labs(fill = 'MGnify Genome Catalogues') +
  scale_fill_manual(values = MetBrewer::met.brewer('VanGogh2', n=11))

ggsave('figures/biomes.pdf',  bg = 'white', 
       width = 45, height = 20, units = 'cm')

