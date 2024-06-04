#!/usr/bin/env bash

export ANALYSIS=/home/def-ilafores/analysis
export NOVEL=scripts/Novelty_minimal_pipeline.sh
export GTDB_V=r214

bash $NOVEL -o Provid19 -g ${GTDB_V} \
	-m $ANALYSIS/projet_PROVID19/Saliva/MAG_analysis/drep_genomes/dereplicated_genomes \
	-s $ANALYSIS/projet_PROVID19/Saliva/preproc/preprocessed_reads.sample.tsv

bash $NOVEL -o Boreal_mosses -g ${GTDB_V} \
	-m $ANALYSIS/boreal_moss/MAG_analysis/drep_genomes/dereplicated_genomes \
	-s $ANALYSIS/boreal_moss/preproc/preprocessed_reads.sample.tsv

#
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
