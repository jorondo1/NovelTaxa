library(pacman)
p_load(tidyverse, MetBrewer, magrittr)

# Old code : MGnify catalogues

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
