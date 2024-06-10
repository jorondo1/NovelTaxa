library(tidyverse, MetBrewer, magrittr)





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