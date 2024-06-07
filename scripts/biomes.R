library(tidyverse, MetBrewer, magrittr)

genomes <- read_csv('figures/GenomeCatalogue.csv') %>% 
  select(name, unclustered_genome_count) %>% 
  mutate(count = unclustered_genome_count/sum(unclustered_genome_count),
         .keep='unused',
         name = case_when(name == 'Unified Human Gastrointestinal Genome (UHGG) v2.0.2' ~ 'Unified Human Gastrointestinal \nGenome (UHGG) v2.0.2', TRUE ~ name))

genomes %>% 
ggplot(aes(x = "", y = count, fill = name)) +
  geom_bar(width = 1, stat = 'identity') +
  coord_polar("y", start = 0) +
  theme_void(base_size = 24) + 
  theme(legend.position = 'right',
        legend.spacing = unit(5,'cm'),
        legend.key.height = unit(1.2, 'cm')
        ) +
  labs(fill = 'MGnify Genome Catalogues') +
  scale_fill_manual(values = MetBrewer::met.brewer('VanGogh2', n=11))

ggsave('figures/biomes.png', 
       width = 3400, height = 2000, units = 'px',
       bg = 'white')



 biomes <- read_csv('figures/Biome.csv')
max_split <- max(sapply(strsplit(biomes$lineage, ':'), length))

biomes %<>%
  separate(lineage, into = paste0("V", 1:max_split), sep = ":", fill = 'right') %>%
  select(-V1) %>%
  filter(V2 %in% c('Environmental', 'Host-associated'))

### Plot biomes piechart
biomes_V3 <- biomes %>%
#  filter(is.na(V4) & !is.na(V3)) %>%
  group_by(biome_cat) %>%
  summarise(samples = sum(samples)) %>%
  mutate(sample_perc = samples/sum(samples)) %>%
  mutate(
         biome_cat = case_when(sample_perc < 0.02 ~ 'Other biomes',
                               TRUE ~ biome_cat),
         biome_cat = case_when(biome_cat == 'Terrestrial' ~ 'Soil',
                               biome_cat == 'Aquatic' ~ 'Water',
                               biome_cat == 'Mammals' ~ 'Other mammals',
                        TRUE ~ biome_cat))
biomes_V3 %>%
  group_by(biome_cat) %>%
  summarise(sample_perc = sum(sample_perc)) %>%