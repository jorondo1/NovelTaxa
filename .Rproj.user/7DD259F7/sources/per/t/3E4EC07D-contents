################################################################################
### AUTHOR : JONATHAN RONDEAU-LECLAIRE #########################################
### Community Abundance Data Processing ########################################
################################################################################

library(pacman)
p_load(tidyverse, magrittr, RColorBrewer, phyloseq, ape)
source("scripts/myFunctions.R")

# Sarah's PS object (based on Kraken taxonomic assignment)
sarah.ps <- readRDS("data/ps_comptype.RDS")

###################
#### Abundance #####
###################

abund_GTDB <- parse_SM("data/SM_abund/*gtdb_gather.csv")
abund_GBNK <- parse_SM("data/SM_abund/*genbank_gather.csv")
abund_MAGs <- parse_SM("data/SM_abund/*custom_gather.csv")

#######################
#### GTDB Taxonomy #####
#######################
taxLvls <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

tax_GTDB <- read_delim("data/gtdbtk_summary.tsv") %>% 
  separate(classification, into = taxLvls, sep = ";", fill = "right") %>% 
  
  # We rename the bins, because we'll add these to the species
  # (since MAG species is not always known, this prevents two species-level bins
  # from the same genus to have the same species name)
  mutate(unique_MAGs = simplify_name(user_genome),
         genome = user_genome) %>% 
  dplyr::select(all_of(c('genome', taxLvls, 'unique_MAGs'))) %>% 
  
  # Adding known GTDB taxa
  bind_rows(
    read_delim("data/gtdb_taxonomy_subset.csv", delim=',', col_names = F) %>% 
      select(-X2) %>% 
      set_names(c('genome',taxLvls))
  ) %>% 
  mutate_at(taxLvls, rename.fun) %>% 
  
  # If the species is unknown, we'll fill the Species field with the nearest known
  # taxonomic rank, and add a unique MAG identifier
  mutate(across(everything(), ~na_if(.x, "")),
         Species = ifelse(
           is.na(Species),
           paste0(coalesce(Genus, Family, Order, Class, Phylum, Domain), 
                  " ", unique_MAGs),
           Species
         )) %>% 
  select(-unique_MAGs)

# Some taxonomic levels have redundancies because higher levels use alternative names,
# herego there can be a duplicate Order whose Class or Phylum is different. Some 
# can be heterotypic synonyms, others outdated taxonomic names.

# Find duplicates; this function will open the viewer with the names to be resolved.
listDupl(tax_GTDB, "Genus")

# Correct the different level duplicates, where "New_name" = "Old_name"
corrPhylum <- c("Pseudomonadota" = "Proteobacteria",
                "Actinobacteriota" = "Actinomycetota",
                "Cyanobacteria" = "Cyanobacteriota")
corrClass <- c("Acidobacteriae" = "Terriglobia")
corrOrder <- c("Acidobacteriales" = "Terriglobales",
               "Enterobacterales_A" = "Enterobacterales")
corrFamily <- c("Enterobacteriaceae_A" = "Enterobacteriaceae",
                "UBA10450" = "SZAS-18",
                "Burkholderiaceae_B" = "Burkholderiaceae")

# correct, but check family too...
tax_GTDB %<>% mutate(Phylum = recode(Phylum, !!!corrPhylum),
                     Class = recode(Class, !!!corrClass),
                     Order = recode(Order, !!!corrOrder),
                     Family = recode(Family, !!!corrFamily))

##########################
#### Genbank Taxonomy #####
##########################

tax_GBNK <- Sys.glob("data/SM_abund/genbank_lineages/*.csv") %>% 
  map_dfr(read_csv, col_types="cccccccccc") %>% # Parse data
  dplyr::select(-taxid) %>%  # drop taxid
  rename_with(~ c("genome", taxLvls, "Strain"), # same colnames as GTDB
              .cols = everything()) %>% 
  rbind(tax_GTDB %>% # Add MAG taxonomy
          mutate(Strain = NA) %>% # rbind needs same # of columns
          filter(str_detect(genome, # take MAG tax only from GTDB tax
                            paste0("^(", paste(c("Brown", "Green"), collapse="|"), ")"))
          )) %>% # then filter to keep only taxa found by sourmash
  right_join(abund_GBNK %>% rownames_to_column("genome") %>% select(genome))

##################
#### Metadata #####
##################

# Species was used to designate Host Species, risk of confusion with Microbiome species
sampleData <- sarah.ps %>% sample_data %>% as("data.frame") %>% 
  mutate(Host = str_c( 
    str_extract(Species, "^[A-Z]"), # Extract the first capital letter
    "_", # add underscore
    str_extract(Species, "(?<= )[a-z]+"), # Extract the second word (species name)
    sep = ""
  ), .keep = 'unused')

# Generate PS objects
mossGTDB.ps <- makePhyloSeq(abund_GTDB, sampleData, tax_GTDB) %>% 
  prune_taxa(taxa_sums(.) > 0,.) # remove taxa abseent from all samples 

mossMAGs.ps <- makePhyloSeq(abund_MAGs, sampleData, tax_GTDB) %>% 
  # TEMPORARY : REMOVE MAG identified as contaminated by GUNC
  prune_taxa(taxa_names(.) !="Green_AO.bin.6",.) %>% 
  prune_taxa(taxa_sums(.) > 0,.) 

mossGBNK.ps <- makePhyloSeq(abund_GBNK, sampleData, tax_GBNK) %>% 
  prune_taxa(taxa_sums(.) > 0,.) # remove taxa abseent from all samples 

write_rds(mossGTDB.ps,"data/R_out/mossGTDB.RDS")
write_rds(mossMAGs.ps,"data/R_out/mossMAGs.RDS")
write_rds(mossGBNK.ps,"data/R_out/mossGBNK.RDS")
