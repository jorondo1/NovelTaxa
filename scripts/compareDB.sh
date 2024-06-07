#!/usr/bin/env bash

export ANALYSIS=/home/def-ilafores/analysis
export NOVEL=scripts/Novelty_minimal_pipeline.sh
export GTDB_V=r214

cd $ANALYSIS/NovelTaxa

bash $NOVEL -o Provid19 -g ${GTDB_V} \
	-m $ANALYSIS/projet_PROVID19/Saliva/MAG_analysis/drep_genomes/dereplicated_genomes \
	-s $ANALYSIS/projet_PROVID19/Saliva/preproc/preprocessed_reads.sample.tsv

bash $NOVEL -o Boreal_mosses -g ${GTDB_V} \
	-m $ANALYSIS/boreal_moss/MAG_analysis/drep_genomes/dereplicated_genomes \
	-s $ANALYSIS/boreal_moss/preproc/preprocessed_reads.sample.tsv

#####################
### community_abundance.R
#####################

# Rscript my_script.R "${GTDB_V}"

# Containment plot
# Diversity plot

# Produce test stats ?? or just annotate diversity plots

######################
### Eventually :
######################

# DB checks that are common to all jobs (skani, checkm, sourmash refs)

# export common variables from here

# simultaneous job monitoring, export so in child script: 
# declare "${base}"="${value}" # for jobID
# export ${base}

# eval_cont applied to all, again using variables?



# These options could be used : (see https://chatgpt.com/share/a984f09d-60e0-4823-b3ca-6b42b94c2587 )
# if [ $# -lt 1 ]; then
#   echo "Usage: $0 <dataset1> [dataset2] [dataset3] ..."
#   exit 1
# fi
#
# dataset1=$1
# shift
#
# optional_args=()
# while (( "$#" )); do
#   optional_args+=("$1")
#   shift
# done
#
# # Output the parsed arguments
# echo "Required argument: $required_arg"
# if [ ${#optional_args[@]} -gt 0 ]; then
#   echo "Optional arguments: ${optional_args[@]}"
# else
#   echo "No optional arguments provided"
# fi
