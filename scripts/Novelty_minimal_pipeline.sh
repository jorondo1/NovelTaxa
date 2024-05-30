#!/usr/bin/env bash

#####################
### SETUP
#####################
export ILAFORES=/home/def-ilafores
export MAIN=${ILAFORES}/analysis/NovelSpecies
export ILL_PIPELINES=${ILAFORES}/programs/ILL_pipelines/containers
export DB=${ILAFORES}/ref_dbs
export SKANI=${ILAFORES}/programs/skani/skani
export GTDB_SKANI=${DB}/GTDB/gtdb_skani_database_ani
export GTDB_VERSION=gtdb_genomes_reps_r214
export ANCHOR=/nfs3_ib/nfs-ip34

# Directory containing MAGs: 
export MAG_DIR=$ILAFORES/analysis/boreal_moss/MAG_analysis/drep_genomes/dereplicated_genomes

# Setup project directories
mkdir -p ${MAIN}/out/checkm2 ${MAIN}/sourmash/sketches ${MAIN}/plots ${MAIN}/logs
cd $MAIN
find ${MAG_DIR} -type f -name '*.fa' > out/MAG_list.txt

#####################
### CHECKM
#####################

# Singularity
if [[ -f out/checkm2/quality_report.tsv ]]; then
	module load apptainer
	
	singularity exec --writable-tmpfs -e -B ${ILAFORES}:${ILAFORES} \
		${ILL_PIPELINES}/containers/checkm2.1.0.2.sif \
		checkm2 predict --threads 48 \
		--database_path ${DB}/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd \
		--input ${MAG_DIR} --extension .fa --output-directory ${MAIN}/out/checkm2

	module unload apptainer
fi

#####################
### skANI
#####################

if [[ ! -d ${GTDB_SKANI} ]]; then
	cd $GTDB_SKANI
	find / -name '*.fna.gz' > gtdb_file_names.txt
	${SKANI} sketch -l gtdb_file_names.txt -o gtdb_skani_database_ani -t 72
	cd $MAIN
fi 

if [[ ! -f out/ANI_results.txt ]]; then
	# Compute ANI
	${SKANI} search -d ${GTDB_SKANI} -o out/ANI_results.txt -t 72 --ql out/MAG_list.txt 
	# Format some columns
	cat out/ANI_results.txt | \
	sed -e "s|${MAG_DIR}/||g" \
	-e "s|${GTDB_VERSION}/database/GC./.../.../.../||g" \
	-e "s/_genomic.fna.gz//g" | \
	awk 'BEGIN {FS=OFS="\t"} {gsub(".fa", "", $2); print}' > tmp.txt
	mv tmp.txt out/ANI_results.txt
fi

#####################
### identify_novel_MAGs.py
#####################

if [[ ! -f scripts/novel_MAGs.py ]]; then
	echo "Program is missing!"
	exit 1
fi

python3 scripts/novel_MAGs.py

#####################
### gather_SLURM.sh
#####################

ml apptainer
export SOURMASH="singularity exec --writable-tmpfs -e -B ${ILAFORES}:${ILAFORES} ${ILAFORES}/programs/ILL_pipelines/containers/sourmash.4.7.0.sif sourmash"
# Sketch novel genomes 
$SOURMASH sketch dna -p scaled=1000,k=31,abund \
	--name-from-first --from-file out/nMAG_list.txt \
	--output-dir sourmash/sketches
# if it fails, check if there are empty returns at the end of nMAG_list.txt
#### Eventually, use 'branchwater multisketch'

# Create an index 
$SOURMASH index sourmash/sketches/nMAGs_index sourmash/sketches/*.sig
export MAGs_IDX=$(find ${MAIN}/sourmash/sketches -type f -name 'nMAGs_index*')

# Gather metagenomes
sbatch --array=1-4 \
	--export=ANCHOR,ILAFORES,DB,MAIN,MAGs_IDX,\
	SAMPLE_DIR=${ANCHOR}${ILAFORES}/analysis/boreal_moss/preproc \
	$PWD/scripts/gather_SLURM.sh

#####################
### parse_sourmash.R
#####################




#####################
### Plots.R
#####################











