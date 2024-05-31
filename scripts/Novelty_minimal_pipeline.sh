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
export SINGULARITY="singularity exec --writable-tmpfs -e -B ${ILAFORES}:${ILAFORES} ${ILL_PIPELINES}/containers"

# Directory containing MAGs: 
export MAG_DIR=$ILAFORES/analysis/boreal_moss/MAG_analysis/drep_genomes/dereplicated_genomes

# Setup project directories
mkdir -p ${MAIN}/out/checkm2 ${MAIN}/sourmash/sketches/nMAGs ${MAIN}/plots ${MAIN}/logs
cd $MAIN
find ${MAG_DIR} -type f -name '*.fa' > out/MAG_list.txt

#####################
### CHECKM
#####################

# Singularity
if [[ -f out/checkm2/quality_report.tsv ]]; then
	module load apptainer
	
	"$SINGULARITY"/checkm2.1.0.2.sif \
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

export SOURMASH="${SINGULARITY}/sourmash.4.7.0.sif sourmash"
export SAMPLE_DIR=${ILAFORES}/analysis/boreal_moss/preproc
export MAGs_IDX=$(find ${MAIN}/sourmash/sketches -type f -name 'nMAGs_index*')
export N_SAM=$(wc ${SAMPLE_DIR}/clean_samples.tsv | awk '{print $1}')

ml apptainer

# Sketch novel genomes 
$SOURMASH sketch dna -p scaled=1000,k=31,abund \
	--name-from-first --from-file out/nMAG_list.txt \
	--output-dir sourmash/sketches/nMAGs
# if it fails, check if there are empty returns at the end of nMAG_list.txt
#### Eventually, use 'branchwater multisketch'

#Fix names for novel genomes 
mkdir moss_MAGs_renamed/
for file in $(find sourmash/sketches/nMAGs -type f -name '*.sig'); do
	new_name=$(basename $file)
	$SOURMASH sig rename $file "${new_name%.fa.sig}" -o sourmash/sketches/nMAGs/${new_name}
done

# Create an index 
$SOURMASH index sourmash/sketches/nMAGs_index sourmash/sketches/nMAGs/*.sig

# Gather metagenomes
sbatch --array=1-"${N_SAM}" \
	--export=ANCHOR,ILAFORES,DB,MAIN,MAGs_IDX,SAMPLE_DIR \
	$PWD/scripts/gather_SLURM.sh

wget https://raw.githubusercontent.com/jorondo1/misc_scripts/main/myFunctions.sh
source myFunctions.sh; rm myFunctions.sh

fix_gtdb data/sourmash # there's a comma problem that staggers the columns in gtdb taxonomy
eval_cont data/sourmash


#####################
### community_abundance.R
#####################




#####################
### Plots.R
#####################











