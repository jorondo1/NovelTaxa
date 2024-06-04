#!/usr/bin/env bash

cd $ILAFORES/analysis/NovelTaxa

export OUTDIR=$PWD/Moss_test
export MAIN=${PWD}
export ILAFORES=/home/def-ilafores
export MAG_DIR=$ILAFORES/analysis/boreal_moss/MAG_analysis/drep_genomes/dereplicated_genomes
export SAMPLE_DIR=$ILAFORES/analysis/boreal_moss/preproc/
export ILL_PIPELINES=${ILAFORES}/programs/ILL_pipelines/containers
export DB=${ILAFORES}/ref_dbs
export SKANI=${ILAFORES}/programs/skani/skani
export GTDB_SKANI=${DB}/GTDB/gtdb_skani_${GTDB_V}
export GTDB_V='rs214'
export ANCHOR=/nfs3_ib/nfs-ip34
export SINGULARITY="singularity exec --writable-tmpfs -e -B ${ILAFORES}:${ILAFORES} ${ILL_PIPELINES}/containers"
export SOURMASH="${SINGULARITY}/sourmash.4.7.0.sif sourmash"
export MAGs_IDX=$(find ${OUTDIR}/sourmash/sketches -type f -name 'nMAGs_index*')
export N_SAM=$(wc ${SAMPLE_DIR}/clean_samples.tsv | awk '{print $1}')

