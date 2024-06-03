#####################
### Formatting ####
#####################

# Moss section colours
compColours <- c('darkgoldenrod4', 'darkolivegreen3')
compColoursContour <- c("tan4", "darkolivegreen4")

mossSpecies <- c("D. undulatum", "P. commune", "P. juniperinum", "P. piliferum")
mossNamesFormatted <- c("D_undulatum" = "D. undulatum",
                        "P_commune" = "P. commune",
                        "P_juniperinum" = "P. juniperinum",
                        "P_piliferum" = "P. piliferum")

### Italicize species names
labelsItal <- c(expression(italic("D. undulatum")),
                expression(italic("P. commune")),
                expression(italic("P. juniperinum")),
                expression(italic("P. piliferum")))

labelsReg <- c(expression("D. undulatum"),
                expression("P. commune"),
                expression("P. juniperinum"),
                expression("P. piliferum"))

### Order colours
col_order =  c(
  "Acetobacterales" = "#D4A6C8FF", "Terriglobales" = "#A0CBE8FF", # Acidobacterales, https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=204433&lvl=3&lin=f&keep=1&srchmode=1&unlock
  "Actinomycetales" = "#F28E2BFF", "Armatimonadales" = "#FFBE7DFF",
  "Baltobacterales" = "#FF9D9AFF", "Bryobacterales" = "#FABFD2FF",
  "Burkholderiales" = "#9D7660FF", "Caulobacterales" = "#34B6C6",
  "Chitinophagales" = "#E15749FF", "Chthoniobacterales" = "#D7B5A6FF",
  "Cyanobacteriales" = "#59A14FFF", "Cytophagales" = "#8CD17DFF",
  "Enterobacterales_A" = "#79706EFF", "Ktedonobacterales" = "#BAB0ACFF",
  "Methylacidiphilales" = "#D37295FF", "Mycobacteriales" = "#961F1F",
  "Rhizobiales" = "#499894FF", "Solirubrobacterales" = "#59386C",
  "Sphingobacteriales" = "#B07AA1FF", "Sphingomonadales" = "#4E79A7",
  "Steroidobacterales" = "#132B69", "Streptosporangiales" = "#86BCB6FF",
  "Xanthomonadales" = "#F1CE63FF", "Others" = "grey",
  "Pseudomonadales" = "hotpink", "Reyranellales" = "#B6992DFF"
)

#####################
### Preprocessing ####
#####################

### Custom sourmash output sparsing function
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/main/community_functions.R'))

### Remove .__ from taxa names
rename.fun <- function(x) {str_remove(x,".__")}

# Renaming MAGs in a more compact way
simplify_name <- function(name) {
  pattern = "^([A-Za-z]{2})[^_]*_([^.]*)\\.bin\\.(\\d+)$"
  replacement = "\\1\\2\\3.MAG"
  str_replace_all(name, pattern, replacement)
}

# Make a phyloseq object 
makePhyloSeq <- function(abund, meta, tax, tree = FALSE) {
  ps <- phyloseq(
    otu_table(abund, taxa_are_rows = T),
    sample_data(meta),
    tax_table((tax %>% 
                 filter(genome %in% rownames(abund)) %>% 
                 column_to_rownames("genome")) %>% as.matrix)
    )
  if (tree!=FALSE) {ps <- read.tree(tree) %>% phyloseq %>% merge_phyloseq(ps)}
  return(ps)
}

# Prune
prunePS <- function(ps, keep) {
  ps %>% taxa_names %>% .[(. %in% keep)] %>% prune_taxa(ps)
}

trans.fun <- function(ps, method, keep) {
  ps %>% tax_filter(min_prevalence=0.05) %>%
    purrr::when(method == 'none' ~ ., # if satisfied, go to next pipe, otherwise
                ~ microbiome::transform(.,transform=method)) %>% # that function we defined above
    prunePS(.,keep)
}

# Some taxonomic levels have redundancies because higher levels use alternative names,
# herego there can be a duplicate Order whose Class or Phylum is different. Some 
# can be heterotypic synonyms, others outdated taxonomic names.

listDupl <- function(tax, level) {
  level_sym <- rlang::sym(level)
  
  subset <- tax %>% dplyr::select(Domain:!!level_sym) %>% 
    dplyr::group_by(!!level_sym) %>%
    unique
  
  dupList <- subset %>% 
    dplyr::summarise(n = n()) %>%
    dplyr::filter(n > 1) %>%
    dplyr::pull(!!level_sym)
  
  # Relative abundance from melted ps object 
  relab.fun <- function(df) {
    summarise(df, Abundance = sum(Abundance, na.rm = TRUE), 
              .groups = 'drop') %>% ungroup %>% 
      group_by(Sample) %>% 
      mutate(relAb = Abundance/sum(Abundance)) %>% ungroup
  }
  
  subset %>% arrange(!!level_sym) %>% 
    filter(!!level_sym %in% dupList) %>% print(n = 1000)
}

##################
### DIVERSITY ####
##################

# Variance-stabilizing transformation
vst.fun <- function(ps, var) {
  phyloseq_to_deseq2(
    ps, as.formula(paste0('~', var))) %>% # DESeq2 object
    estimateSizeFactors(., geoMeans = apply(
      counts(.), 1, function(x) exp(sum(log(x[x>0]))/length(x)))) %>% 
    DESeq2::varianceStabilizingTransformation(blind=T) %>% # VST
    SummarizedExperiment::assay(.) %>% t %>% 
    { .[. < 0] <- 0; . } # replace negatives by zeros
}

# PCOA
pcoa.fun <- function(ps, var, vst.mx, dist) {
  # compute distance matrix
  dist.mx <- vegan::vegdist(vst.mx, distance = dist)
  # compute PCoA
  PCoA <- capscale(dist.mx~1, distance = dist)
  # Print first 3 coordinates
  eig <- round(PCoA$CA$eig[1:3]/sum(PCoA$CA$eig),2)
  message("First 3 principal coordinates :")
  message(paste(eig[1], ',', eig[2], ',', eig[3]))
  # create output list
  out <- data.frame(sample_data(ps))
  out$PCoA1 <- scores(PCoA)$sites[,1]
  out$PCoA2 <- scores(PCoA)$sites[,2]
  list(metadata = out, eig = PCoA$CA$eig, dissMx = vst.mx)
}

##########################
### Differential Abundance Analysis ####
########################

# Parse the DAA results for a pairwise or dunnet test
# We want to keep only taxa for which at least one differential test is 
# significant AND which passed the sensitivity analysis. 

# Does NOT work with the global test, as it does not produce LFC values.

parse_DAA_results <- function(DAA, # ANCOMBC output
                              test, # 'pair' or 'dunnet'
                              thr, # p-value threshold
                              gr, # group variable name (string)
                              taxRank, # Taxonomic rank tested
                              ps) { # phyloseq object
  
  # define results list item name 
  listItem <- paste0('res_', test)
  
    # Extract suffixes for every relevant column
    suffixes <- DAA[[listItem]] %>% 
      dplyr::select(starts_with("diff_")) %>% colnames %>% 
      str_replace(paste0("diff_",gr),"")
    
    # Create a character vector of conditions
    conditions <- purrr::map_chr(suffixes, ~ paste0(
      "(`diff_", gr, .x, # passed sensitivity analysis:
      "` == TRUE & `passed_ss_", gr, .x, # adj. p value above certain threshold:
      "` == TRUE & `q_", gr, .x, "` < ", thr, ")" 
    )) %>% paste(collapse = " | ")

  DAA[[listItem]] %>%
    # drop irrelevant columns: 
    dplyr::select(-starts_with('W_'), -starts_with('p_')) %>% 
    # evaluate all these conditions through filter():
    filter(eval(parse(text = conditions))) %>% 
    # pivot to a long dataset, with one line per test per taxa
    pivot_longer(cols = -taxon, 
                 names_to = c(".value", "Group"), 
                 names_pattern = "(lfc|se|q|diff|passed_ss)_(.+)", 
                 values_drop_na = TRUE) %>% 
    mutate( # change LFC values to 0 if not significant (for the plot):
      across(c(lfc,se), ~ case_when(diff == FALSE ~ 0, TRUE ~ .x)),
      # plot colours:
      textcolour = case_when(lfc==0 ~ "white", TRUE ~ "black"),
      # format group names (needs improvement...)
      Group = str_remove(Group, gr)) %>% 
    # add taxonomy to the data (to allow colouring by higher rank)
    left_join(ps@tax_table %>% as.data.frame %>% tibble %>% 
                dplyr::select(Domain:everything()[which(names(.) == taxRank)]) %>% 
                unique,
              by = c("taxon" = taxRank)) %>% 
    # remove all taxa for which there is no difference in any group:
    group_by(taxon) %>% 
    filter(!all(diff == FALSE)) %>% ungroup()
  # add sth similar to remove pairwise comparisons that have 0 DA taxa?
}

################
### PLOTTING ####
################

# Compute top taxa, aggregate residual in others
topTaxa <- function(psmelt, comp, taxLvl, topN) {
  psmelt %>% 
    filter(Compartment == comp) %>% # subset compartment
    group_by(!!sym(taxLvl)) %>% # group by tax level
    summarise(relAb = mean(relAb)) %>% # find top abundant taxa
    arrange(desc(relAb)) %>% 
    mutate(aggTaxo = as.factor(case_when( # aggTaxo will become the plot legend
      row_number() <= topN ~ !!sym(taxLvl), #+++ We'll need to manually order the species!
      row_number() > topN ~ 'Others'))) # +1 to include the Others section!
}

# Abundance summary dataset by compartment using top taxa
df_comm <- function(psmelt, comp, taxLvl, topTaxa) {
  # use topTaxa to join aggTaxo variable
  left_join(psmelt, topTaxa, by=taxLvl) %>% 
    dplyr::filter(Compartment == comp) %>% # subset
    aggregate(Abundance~aggTaxo+Host+Compartment, # Abundance is aggregated...
              data=., FUN = sum) %>% # ...sums "Others" taxa by specified variables
    mutate(aggTaxo = factor(aggTaxo,# reorder 
                            levels = topTaxaLvls)) 
}

# Extract tree labels in the order they appear :
pullTreeLabels <- function(plot, rank) {
  plot$data %>% filter(!is.na(label)) %>% 
    arrange(y) %>% select(!!sym(rank)) %>% 
    # dplyr way of collapsing a tibble into a vector is pull() :
    pull %>% unique
}
