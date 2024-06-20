library(pacman)
p_load(tidyverse, phyloseq, magrittr, foreach, iterators, parallel,
       extrafont, showtext, patchwork, ggrepel)

source('scripts/comm_comp_data.R')
custom_theme <- theme_light(base_size = 34)
dbCols <- c('#D4FADD','#00A759', '#356946')
mossCols <- c('#00A4F1', '#0070B8')
showtext_auto()

### Containment across databases and projects
p1 <- cntm.df %>% # Format db and dataset names
  mutate(Dataset = case_when(Dataset == "Boreal_mosses" ~ "Boreal mosses",
                             Dataset == "Saliva" ~ "Human saliva"),
         db = case_when(db == dbNames[1] ~ 'Default GTDB r214',
                        db == dbNames[2] ~ 'with nMAGs',
                        db == dbNames[3] ~ 'with nMAGs and better MAGs')) %>% 
  # Plot !
  ggplot(aes(x = db, y = cntm)) +
  geom_violin(aes(fill=db)) + # so dots don't appear in the legend!
  #geom_jitter(alpha = 1, width=0.08, size = 1) + 
  custom_theme +
  facet_wrap(~Dataset, scales= 'free', nrow=2) +
  labs(y = 'Sample k-mer containment', x = '', fill = '') +
  scale_fill_manual(values = dbCols) 

### Diversity by host across database 
div_long_formatted <- div_long %>% # Format db names
  filter(db != dbNames[3]) %>% # bMAGs were not relevant for mosses
  mutate(db = factor(db, levels = dbNames),
         db = case_when(db == dbNames[1] ~ 'Default GTDB r214',
                        db == dbNames[2] ~ 'with nMAGs'),
         Microbiome = case_when(
                        Microbiome == "Boreal_mosses" ~ "Boreal mosses",
                        Microbiome == "Saliva" ~ "Human saliva"))

### Sample-wise change in diversity 
p2 <- div_long_formatted %>% 
  mutate(Title = 'Changes in diversity (Boreal mosses)') %>% # title for facet wrap
  ggplot(aes(x = db, y = Shannon)) +
  geom_boxplot(aes(fill = db)) +
  geom_point(alpha = 0.4, size=2) + geom_line(aes(group = Sample), alpha=0.3) + 
  custom_theme +
  facet_wrap(~Microbiome, scales= 'free', nrow=2) +
  labs(y = 'Shannon diversity') +
  scale_fill_manual(values = dbCols[1:2]) +
  guides(fill = 'none') # remove legend altogether


p1 + theme(plot.margin = unit(c(0,45,0,0), "pt")) + p2 + 
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom',
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.margin = margin(b = -10, t = 25),
        legend.key.height = unit(1.2, 'cm'),
        legend.key.width = unit(1.2, 'cm'),
        axis.text.y = element_text(size = 18, colour = 'black'),
        strip.background = element_rect(fill = "grey50"),
        strip.text = element_text(face = 'bold')) 

ggsave('figures/cntm_div.pdf',  bg = 'white', 
       width = 40, height = 50, units = 'cm')

# Plot diversity changes
div_long_formatted %>% 
  filter(Microbiome=='Boreal mosses') %>% 
  ggplot(aes(x = db, y = Shannon, fill = Host)) +
    geom_boxplot() + theme_light() +
    facet_wrap(~db, scales = "free", nrow = 2) + 
    labs(y = 'Shannon diversity') +
    custom_theme +
    scale_fill_manual(values = c("#D87A00","#FF5552", "#FFE081", "#469aac"),
                      labels = c("D. undulatum", "P. commune", "P. juniperinum", "P. piliferum")) +
    theme(axis.text.x = element_blank(),
          legend.position = 'bottom',
          legend.title = element_blank(),
          legend.margin = margin(b = -10, t=-15),
          axis.text.y = element_text(size = 18, colour = 'black'),
          axis.title.x = element_blank(),
          legend.key.height = unit(1.2, 'cm'),
          strip.background = element_rect(fill = "grey50"), # facet_wrap header
          strip.text = element_text(face = 'bold')) +
    guides(fill = guide_legend(nrow = 2))

ggsave('figures/div.pdf',  bg = 'white', 
       width = 20, height = 50, units = 'cm')

###############################
# BIOME REPRESENTATION PLOT ###
###############################
biomes <- read_csv('figures/Biome.csv')
max_split <- max(sapply(strsplit(biomes$lineage, ':'), length))

# Split lineages and subset
biomes %<>%
  separate(lineage, into = paste0("V", 1:max_split), sep = ":", fill = 'right') %>%
  select(-V1) %>% 
  filter(V2 %in% c('Host-associated', 'Environmental'))

# Format data and compute top biomes
biomes_cat <- biomes %>% 
  group_by(V3) %>% 
  summarise(samples = sum(samples_count)) %>% ungroup %>% 
  mutate(sample_perc = samples/sum(samples), # percentage, then format names
         biome_cat = case_when(V3 == 'Terrestrial' ~ 'Soil',
                               V3 == 'Aquatic' ~ 'Water',
                               V3 == 'Mammals' ~ 'Other mammals',
                               is.na(V3) ~ 'Other biomes',
                               sample_perc < 0.02 ~ 'Other biomes',
                               TRUE ~ V3)) %>% 
  group_by(biome_cat) %>% # summarise Others together
  summarise(samples = sum(samples))

# reorder factors for plot
biome_rank <- c("Human", "Other mammals", "Water", "Plants", "Soil", "Birds", "Other biomes")

biomes_cat %<>% mutate(biome_cat = factor(biome_cat, levels = biome_rank)) 

### Plot biomes piechart
biomes_cat %>% 
  ggplot(aes(x = "", y = samples, fill = biome_cat)) +
  geom_bar(width = 1, stat = 'identity') +
  coord_polar("y", start = 0) +
  theme_void(base_size = 34) + 
  theme(legend.position = 'right',
        legend.key.height = unit(1.4, 'cm'),
        legend.title = element_text(face = 'bold')) +
  labs(fill = 'Biome') +
  scale_fill_manual(values = c("#00a759", "#FFE081", "#388aff", "#D87A00","#FF5552", "#a183b3", "#469aac"))

ggsave('figures/biomes.pdf',  bg = 'white', 
       width = 45, height = 20, units = 'cm')



# ### Containment increase from adding nMAGs (Mosses only)
# cntmDiff %>% 
#   ggplot(aes(x = Dataset, y = fold_inc, fill = Dataset)) +
#   geom_boxplot() + 
#   theme_minimal() +
#   labs(y = 'Fold increase in sample k-mers containment', x = 'Dataset') +
#   scale_y_continuous(limits = c(0, NA), breaks = seq(0, max(cntmDiff$fold_inc), by = 2))  # Set minimum y value and increments
