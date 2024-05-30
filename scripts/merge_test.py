#!/usr/bin/env python3
# coding=utf-8

import pandas as pd
import re

MAGs = pd.read_csv('data/List.txt', sep='\t', header = None)
MAGs.columns = ['MAG_ID']

##################
### Mash output
##################
mash = pd.read_csv('data/mash_dist.tab', sep='\t', header = None).iloc[:,[0,1,2]]
mash.columns = ['Mash_hit', 'MAG_ID', 'Mash_value']

##################
### Parse checkm output and compute quality score
##################
checkm = pd.read_csv('data/quality_report.tsv', sep='\t')[['Name','Completeness', 'Contamination']]
checkm['QS'] = checkm['Completeness'] - 5 * checkm['Contamination']

##################
### skANI output; keep MAGs that have >95% ANI over >60% of aligned fraction with at least 1 ref genome.
##################

# Extract genome reference
def extract_genome_name(path):
    name = 'GC' + path.rsplit('/GC', 1)[-1]
    return name.replace('_genomic.fna.gz', '')

# Parse data
ani_full = (
    pd.read_csv('data/ANI_results.txt', sep = '\t')
    [['Ref_file','Query_file', 'ANI', 'Align_fraction_query']]
    .assign(Ref_genome = lambda x: x['Ref_file'].apply(extract_genome_name))
    .assign(Query_file = lambda x: x['Query_file']
        .str.replace('.fa', '', regex = False)
        .str.replace('MAGs/', '', regex = False)
        )
    .drop('Ref_file', axis=1)
)

ani_95 = ani_full[ani_full['ANI']>95.00][['Query_file']]
##################
### GTDB-Tk output
##################

def extract_species(tax):
    species = tax.split(';')[-1]
    return species.replace('s__', '') if 's__' in species else None

gtdbtk = (
    pd.read_csv("data/gtdbtk.bac120.summary.tsv", sep='\t')
    [['user_genome', 'classification', 'red_value', 'closest_placement_ani', 'closest_placement_af']]
    .assign(Species=lambda x: x['classification'].apply(extract_species))
    )

# Populate table
MAGs['mash_05'] = ~ MAGs['MAG_ID'].isin(mash['MAG_ID']) # True if min Mash dist is >0.05
MAGs['MAG_ID'] = MAGs['MAG_ID'].str.replace('.fa','', regex = False) # reformat MAG ID
MAGs['QS_50'] = MAGs['MAG_ID'].map(checkm.set_index('Name')['QS'] >= 50.00) # True if QS is high enough
MAGs['ANI_95'] = ~ MAGs['MAG_ID'].isin(ani_95['Query_file']) # True if there are no genomes with ANI >95 & AF >60
MAGs['GTDB_s'] = MAGs['MAG_ID'].map(gtdbtk.set_index('user_genome')['Species'] == '') # True if there is no GTDB species assignment

# Extract completeness from gtdb metadata
def genome_finder(s,ID):
    return any(ID in s for ID in ID)

# Extract QS from GTDB metadata 
metadata = pd.read_csv('data/bac120_metadata_r220_short.tsv', sep='\t')
genome_list = ( # Extract MAGs with known species-level genome 
    ani_full.query('ANI > 95.00 & Align_fraction_query > 60')
    ['Ref_genome'].drop_duplicates().tolist()
)
species_hits = metadata[metadata['accession'].apply(lambda x: genome_finder(x, genome_list))]
species_hits.loc[:,'QS_ref'] = species_hits['checkm2_completeness'] - 5 * species_hits['checkm2_contamination']
species_hits.loc[:,'accession'] = species_hits['accession'].apply(lambda x: re.sub(r'.._', '', x, count=1))

MAG_genome_conversion = species_hits.merge(ani_full, left_on='accession', right_on='Ref_genome', how = 'left')
MAGs = MAGs.merge(MAG_genome_conversion[['Query_file', 'QS_ref']], left_on= 'MAG_ID', right_on='Query_file', how = 'left')
print(MAGs)

MAGs.to_csv('data/MAGs.tsv', sep="\t", index=False)


### TO DO !
# show actual QS value for MAG