################################################################################
### AUTHOR : JONATHAN RONDEAU-LECLAIRE #########################################
### Manuscript Plots ###########################################################
################################################################################

library(pacman)
p_load(ape, tidyverse, magrittr, RColorBrewer, colorRamp2, patchwork,
       ggtree, ggtreeExtra, treeio, ggnewscale, cowplot, MetBrewer)
source("scripts/myFunctions.R")
moss.ps <- readRDS("data/R_out/mossMAGs.RDS")

# Taxonomic level for trees and groupings
taxLvl <- 'Order'

###################################
#### PLOT 1. Community Overview ####
###################################

topN=10 # Desired # top taxa to display

# Preliminary dataset with variables of interest
MAGs_melt <- moss.ps %>% psmelt %>%
  dplyr::select(Sample, Abundance, Compartment, Host, Domain:Species) %>% 
  # Compute Class level abundance :
  group_by(Sample, across(all_of(taxLvl)), Compartment, Host) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), 
            .groups = 'drop') %>% 
  group_by(Sample) %>% # Convert to relative abundance
  mutate(relAb = Abundance/sum(Abundance)) %>% ungroup

# Find top taxa:
topTaxa_Brown <- topTaxa(MAGs_melt, 'Brown', taxLvl, topN)
topTaxa_Green <- topTaxa(MAGs_melt, 'Green', taxLvl, topN)

# Create ordered list of taxa
topTaxaLvls <- rbind(topTaxa_Brown, topTaxa_Green) %>% 
  group_by(aggTaxo) %>% 
  aggregate(relAb ~ aggTaxo, data = ., FUN = sum) %>% 
  # order taxa by mean relative abundance
  arrange(relAb) %$% aggTaxo %>% as.character %>% 
  # put "Others" first
  setdiff(.,'Others') %>% c('Others',.)

# Build plots dataframe :
df_compart <- rbind(df_comm(MAGs_melt, 'Brown', taxLvl, topTaxa_Brown),
      df_comm(MAGs_melt, 'Green', taxLvl, topTaxa_Green)) %>% 
  mutate(Compartment = factor(Compartment, levels = c('Green', 'Brown')))

# Plot !
comm.p <- ggplot(df_compart, aes(x = Host, y = Abundance, fill = aggTaxo)) +
  geom_bar(stat = "identity", position = "fill",
           colour = 'black', size = 0.2) +
  facet_wrap('Compartment', ncol = 1) +
  labs(fill = taxLvl, 
       y = paste("Mean relative sequence abundance by sample")) +
  scale_fill_manual(values = col_order, breaks = topTaxaLvls) +
  scale_x_discrete(labels = labelsReg) +
  theme_light() + # fix facet headers :
  theme(strip.background =element_rect(fill = 'white'),
    strip.text = element_text(colour = 'black', size = 18),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=14, color="black"), 
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.y  = element_text(size=14, color="black"),
    legend.title = element_text(colour="black", size=14, face="bold"),
    legend.text = element_text(colour="black", size = 12))

#ggplot2::ggsave("out/community.png", bg = 'white', plot = comm.p,
#                width = 1600, height = 2400, units = 'px', dpi = 240)

write_rds(comm.p, 'out/community.RDS')

#####################################
#### PLOT 2. MAGs characteristics ####
#####################################

# extract MAG names
MAG_names <- moss.ps@tax_table %>% rownames %>% .[grep(".bin.", .)] 

# Extract taxonomy
taxLabels <- moss.ps@tax_table %>% as.data.frame %>% 
  rownames_to_column("label") %>% 
  select(label, Domain:Species) %>% 
  tibble

MAG.tree <- read.tree("data/RAxML_bestTree.MAGs_refined.tre")
# Subset tree to MAGs only and add taxonomy
sub.tree <- MAG.tree %>% 
  drop.tip(setdiff(MAG.tree$tip.label, MAG_names)) %>% 
  as_tibble %>% 
  full_join(taxLabels %>% filter(label %in% MAG_names), by = 'label') %>% 
  as.treedata

# Expand a Brewer palette to 12 colours:
nTax <- moss.ps %>% prune_taxa(taxa = MAG_names, .) %>% 
  tax_table %>% as.data.frame %$% Order %>% unique %>% length

# Plot the tree
p <- ggtree(sub.tree, layout="fan", size=0.1) +
  xlim(-0.2, NA) + # prevent the high-level branches from clustering in the middle
  geom_tippoint(mapping = aes(color = Order), size = 2.5) +
  scale_colour_manual(values = col_order) +
  guides(color = guide_legend(ncol = 1)) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

p$data %<>% 
  mutate(!!sym(taxLvl) := # dynamic management of variable name
           factor(!!sym(taxLvl), levels = pullTreeLabels(p, taxLvl)))

# add MAG data:
hm.mx <- read_tsv("data/R_out/MAG_summary.tsv") %>% 
  # only keep MAGs that are in our original dataset
  filter(MAG %in% MAG_names) 

# Add a heatmap of %GC Content :
p2 <- p + guides(colour = "none") +
    geom_fruit(hm.mx, geom_tile, mapping=aes(y=MAG, fill = GC),
               offset = 0.1,width = 0.1) +
  scale_fill_gradient(low = 'seagreen4', high = 'lightgoldenrod1', name="% GC",
                      guide = guide_colourbar(order = 1)) + #enforce legend display order
  # and Length
  new_scale_fill() +
  geom_fruit(hm.mx, geom_tile, mapping = aes(y = MAG, fill = MBP),
             offset = 0.1, width = 0.1) +
  scale_fill_gradient(low = 'ghostwhite', high = 'darkorange1', name="Length (Mbp)",
                      guide = guide_colourbar(order = 2)) +
  # and n_contigs
  new_scale_fill() +
  geom_fruit(hm.mx, geom_tile, mapping = aes(y = MAG, fill =  `Number of contigs`),
             offset = 0.1, width = 0.1) +
  scale_fill_gradient(low = 'ghostwhite', high = 'purple4',
                      guide = guide_colourbar(order = 3)) +
  # and N50
  new_scale_fill() + 
  geom_fruit(hm.mx, geom_tile, mapping = aes(y = MAG, fill = `N50`),
             offset = 0.1, width = 0.1) +
  scale_fill_gradient(low = "ghostwhite", high = 'dodgerblue4',
                      guide = guide_colourbar(order = 4)) +
  # and Quality score
  new_scale_fill() +
  geom_fruit(hm.mx, geom_tile, mapping = aes(y = MAG, fill = QS),
             offset = 0.1, width = 0.1) + 
  scale_fill_gradient(low = "ghostwhite", high = "red4", name="Quality Score",
                      guide = guide_colourbar(order = 5))

#p2 + geom_hilight(node = 121, fill = "NA",size= 5)

# Extract legends as grobs
order_legend <- get_legend(p)
gradient_legends <- get_legend(p2)

(tree_MAGs.plot <- cowplot::ggdraw(
  plot_grid(order_legend, p2 + theme(legend.position = 'none'), gradient_legends, 
            nrow = 1, ncol = 3, scale = c(3,1.2,1), rel_widths=c(3, 5,2))))

ggplot2::ggsave("out/tree_MAGs.png", bg = 'white',
                width = 4200, height = 2600, units = 'px', dpi = 300)

#########################################
### 3. DIFFERENTIAL ABUNDANCE by HOST ####
#########################################
hostDA <- read_rds('data/R_out/DA_host_results.RDS')

# sort Species by taxLvl (descending so species are top-to-bottom on y axis)
speciesLvl <- hostDA %>% arrange(desc(!!sym(taxLvl)), taxon) %$% taxon %>% unique

DA_host.df <- hostDA %>% 
  # reorder taxa by taxLvl
  mutate(taxon = factor(taxon, levels = speciesLvl))

# MAIN DAA PLOT :
DA.p1 <- DA_host.df %>% 
  ggplot(aes(x = Group, y = taxon, fill = lfc)) +
  geom_tile() +
  scale_fill_gradient2(low = met.brewer("Cassatt1")[1], 
                       mid = "white", 
                       high = met.brewer("Cassatt1")[8], 
                       midpoint = 0) +
  geom_text(aes(Group, taxon, label = round(lfc, 2), color=textcolour)) +
  scale_color_identity(guide = FALSE) + 
  theme_void() + 
  theme(axis.text.x = element_text(margin = margin(t = 5, r = 5, b = 5, l = 5))) +
  labs(x = '', y = '', fill = "LFC")

# taxLvl -coloured tile:
DA.p2 <- DA_host.df %>% 
  ggplot(aes(x = '1', y = taxon, fill = !!sym(taxLvl))) +
  geom_tile() + theme_void() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(hjust = 1)) +
  scale_fill_manual(values = col_order) 

(DA_host_tree_05 <- 
  DA.p2 + DA.p1 + plot_layout(
    guides = "collect",
    design = "ABBBBBBBBBBBB") &
  scale_x_discrete(labels = labelsReg[2:4], 
                   position = 'bottom'))

ggplot2::ggsave("out/DA_host_species.png", plot = DA_host_tree_05,
                width = 2700, height = 3600, units = 'px', dpi = 300)

################################################
### 4. DIFFERENTIAL ABUNDANCE by COMPARTMENT ####
################################################

# Add LFC (from DAA) to significant species
speciesLFC <- readRDS("data/R_out/speciesLFC_comp.RDS")
tree <- read.tree("data/RAxML_bestTree.genomes_refined.tre") 

# Subset taxa for tree layer
DA_species <- speciesLFC %$% MAG
DA_sub.tree <- tree %>% 
  drop.tip(setdiff(tree$tip.label, DA_species)) %>% as_tibble %>% 
  # add taxonomy :
  full_join(taxLabels %>% filter(label %in% DA_species), by = 'label') %>% 
  as.treedata # because.

n <- DA_sub.tree@data %>% as.data.frame %>% .[taxLvl] %>% unique %>% dim %>% .[1]

### Taxonomic tree (generated first to establish species factor levels)
tree.p <- ggtree(DA_sub.tree,size = 0.2) 

# Reorder factor levels 
tree.p$data %<>% 
  mutate(!!sym(taxLvl) := # dynamic management of variable name
           factor(!!sym(taxLvl), levels = pullTreeLabels(tree.p, taxLvl) %>% rev))

tree.p <- tree.p +
  geom_tippoint(mapping = aes(color = !!sym(taxLvl)), size = 3) +
  #geom_tiplab(align=TRUE, aes(label = Species)) +
  scale_colour_manual(values = col_order) +
  scale_fill_manual(values = compColours) +
  labs(fill = "Compartment",
       colour = taxLvl) +
  theme(plot.margin = margin(t=20, r=0, b=20,l = 20),
        legend.margin = margin(t=30),
        axis.title.x = element_text(hjust = 0.95)) +
    theme(legend.position = c(.1,.8))

# Extract species levels for alignment of next plots :
orderedSpecies <- tree.p$data %>% select(y, label, Species) %>% 
  unique %>% arrange(y) %>% # sort by plot position
  filter(!is.na(label)) %$% Species # NAs at non-integer positions?!

### Waterfall plot :
speciesLFC %<>% mutate(Species = factor(Species, levels = orderedSpecies))
  
waterfall.p <- speciesLFC %>% 
  #filter(Species %in% subset) %>% 
  ggplot(aes(y = Species, x = LFC, fill = compAss)) + 
  geom_bar(stat = "identity", width = 1, color = "white",
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(xmin = LFC - SE, 
                    xmax = LFC + SE), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(y = NULL, x = "LFC", fill = 'Moss section\nassociation') + 
  scale_fill_manual(values = c("Green" = compColours[2], "Brown" = compColours[1]))+
  scale_color_discrete(name = NULL) +
  theme_minimal() + 
  theme(panel.grid.minor.y = element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 10, l = -5, unit = "pt"),
        legend.position = c(.8,.95),
        legend.background = element_rect(colour='black', fill='white', linewidth=0.2),
        legend.title = element_text(size = rel(.8)))

tree.p + waterfall.p + plot_layout(design = "AAAABBB") #+
#    plot_annotation(caption = 'Significantly abundant at p<0.05 (ajdusted).\nRestricted to species with >10% relative abundance that passed the sensitivity analysis.')

ggplot2::ggsave("out/DA_comp_tree.png", 
                width = 2700, height = 3600, units = 'px', dpi = 300)

#######################################################
### S-X. DIFFERENTIAL ABUNDANCE by HOST with ORDERS ####
#######################################################

read_rds('data/R_out/DA_pw_host_Order.RDS') %>% 
  parse_DAA_results('pair', 0.01, 'Host', 'Order', moss.ps) %>% 
  #  mutate(taxon = factor(taxon, levels = speciesLvl)) %>% 
  ggplot(aes(x = Group, y = taxon, fill = lfc)) +
  geom_tile() +
  scale_fill_gradient2(low = met.brewer("Cassatt1")[1], 
                       mid = "white", 
                       high = met.brewer("Cassatt1")[8], 
                       midpoint = 0) +
  geom_text(aes(Group, taxon, label = round(lfc, 2), color=textcolour)) +
  scale_color_identity(guide = FALSE) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(margin = margin(t = 5, r = 5, b = 5, l = 5))) +
  labs(x = '', y = '', fill = 'LFC') + 
  scale_x_discrete(labels = c('P. commune\n(D. undulatum)',
                              'P. juniperinum\n(D. undulatum)',
                              'P. juniperinum\n(P. commune)',
                              'P. piliferum\n(D. undulatum)',
                              'P. piliferum\n(P. commune)',
                              'P. piliferum\n(P. juniperinum)'
  ))

ggplot2::ggsave("out/DA_host_Orders.png", bg = 'white',
                width = 3600, height = 2800, units = 'px', dpi = 300)

###############################################
### S-X. CONTAINMENT IMPROVEMENT USING MAGS ####
###############################################

dbNames <- c("GTDB", "GTDB + MAGs", "Genbank", "Genbank + MAGs")

raw <- read_delim("data/cntm_sum.txt", col_names = c("Sample", "db", "cntm")) %>% 
  filter(Sample %in% (moss.ps@sam_data %>% rownames))

cntm.df <- raw %>% 
  mutate(db = case_when(db == "gtdb" ~ dbNames[1],
                        db == "custom" ~ dbNames[2],
                        db == "genbank_default" ~ dbNames[3],
                        db == "genbank" ~ dbNames[4]
  ),
  db = factor(db, levels = dbNames))

# Mean cntm ± sd
cntm.df %>% group_by(db) %>% 
  summarise(mean_cntm = mean(cntm),
            sd_cntm = sd(cntm),
            n_sam = n())

# Containment increase by db
cntm.df %>% 
  # one row per sample
  pivot_wider(names_from = 'db', values_from = 'cntm') %>% 
  # compute increases
  mutate(GTDB_increase = `GTDB + MAGs` - GTDB,
         Genbank_increase = `Genbank + MAGs` - Genbank) %>% 
  summarise(GTDB_mean = mean(GTDB_increase),
            GTDB_sd = sd(GTDB_increase))

(ggplot(cntm.df, aes(x = db, y = cntm, fill=db)) +
    geom_violin() + geom_jitter(alpha = 0.2) + 
    theme_minimal(base_size = 16) + guides(fill = 'none') +
    labs(y = 'Sample k-mer containment', x = 'Reference database'))

ggplot2::ggsave("out/cntm_comparison.png", bg = 'white',
                width = 2700, height = 2400, units = 'px', dpi = 300)


(cntm.df %>% filter(db == dbNames[2] | db == dbNames[4]) %>% 
  ggplot(aes(x = db, y = cntm, fill=db, group = Sample)) +
    geom_point() + geom_line() +
    theme_minimal(base_size = 16) + guides(fill = 'none') +
    labs(y = 'Sample k-mer containment', x = 'Reference database'))

# Increase in containment across databases
cont.increase <- raw %>% 
  pivot_wider(names_from = db, values_from = cntm) %>% 
  dplyr::select(Sample, custom, genbank) %>% 
  mutate(increase = 100*(genbank - custom) / custom) %>% 
  left_join(moss.ps@sam_data %>% data.frame %>% rownames_to_column("Sample"), 
            by = "Sample")

cont.increase %>% 
  ggplot(aes(x = Host, y = increase, colour = Host)) + 
  geom_boxplot() +
  theme_light() +
  labs(y = "% increase in containment") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        panel.grid = element_blank())

cont.increase %>% 
  ggplot(aes(x = Compartment, y = custom, colour = Compartment)) + 
  geom_boxplot() +
  theme_light() +
  labs(y = "containment") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        panel.grid = element_blank())


