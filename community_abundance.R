library(pacman)
p_load(tidyverse, phyloseq)

# Sourmash parser function
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/community_functions.sh'))

custom <- parse_SM('data/sourmash/*custom_gather.csv')
rs214 <- parse_SM('data/sourmash/*rs214_gather.csv')

