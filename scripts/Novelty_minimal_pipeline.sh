#!/usr/bin/env bash
message="Identify novel species-level MAGs and leverage them for community composition estimation."
echo "$message"

#####################
### PARSING OPTIONS
#####################
function Usage {
    echo "$message"
	echo "Current directory should include scripts/novel_MAGs.py."
	echo "Assumes the GTDB reference genomes are under /home/def-ilafores/ref_dbs/GTDB/gtdb_genomes_reps_*"
	echo "Options:"
    echo "-o	Output directory name."
	echo "-g	GTDB database version (default: r220)"
	echo "-m	directory containing MAGs"
	echo "-s	tab-delimited file with paths to samples (one sample per line, one column per file, first column is the sample ID)"
    exit 1
}

while getopts 'o:g:m:s:' flag; do
    case "${flag}" in
		o) export OUTDIR=${PWD}/${OPTARG} ;;
		g) export GTDB_V=${OPTARG} ;;
		m) export MAG_DIR=${OPTARG} ;;
		s) export SAM_LIST=${OPTARG} ;;
    esac
done

#####################
### SETUP
#####################

if [[ -z "$OUTDIR" ]] || [[ -z "$MAG_DIR" ]] || [[ -z "$SAM_LIST" ]] ; then
	echo "Argument(s) missing."
	Usage
fi
	
if [[ -z "$GTDB_V" ]]; then
	export GTDB_V=r220
fi

if [[ ! -f ${PWD}/scripts/novel_MAGs.py ]]; then
	echo "novel_MAGs.py can't be found!"
	Usage
fi

export MAIN=${PWD}
export ILAFORES=/home/def-ilafores
export SINGULARITY="singularity exec --writable-tmpfs -e -B ${ILAFORES}:${ILAFORES} ${ILAFORES}/programs/ILL_pipelines/containers"
export SOURMASH="${SINGULARITY}/sourmash.4.7.0.sif sourmash"
export DB=${ILAFORES}/ref_dbs
export SM_SK=$OUTDIR/tmp/sourmash/sketches
export SKANI=${ILAFORES}/programs/skani/skani
export GTDB_SKANI=${DB}/GTDB/gtdb_skani_${GTDB_V}
export ANCHOR=/nfs3_ib/nfs-ip34
export N_SAM=$(wc ${SAM_LIST} | awk '{print $1}')
export THREADS=48

if [[ ! -d $GTDB_SKANI ]] && [[ ! -d ${DB}/GTDB/gtdb_genomes_reps_${GTDB_V} ]]; then
	echo 'Reference skANI database unavailable. Download and store under "${DB}/GTDB".'
	Usage
fi
if [[ ! -f ${DB}/sourmash_db/gtdb-rs${GTDB_V/r}-reps.k31.zip ]]; then
	echo 'Reference Sourmash database unavailable. Currently available: '
	ls $DB/sourmash_db/gtdb*k31.zip
	Usage
fi

# Setup project directories
mkdir -p scripts ${SM_SK}/nMAGs ${OUTDIR}/tmp/logs ${OUTDIR}/tmp/checkm2 ${OUTDIR}/sourmash tmp/

export MAGs_IDX=$(find ${SM_SK} -type f -name 'nMAGs_index*')

# gather output post-processing functions
curl -s -o scripts/myFunctions.sh https://raw.githubusercontent.com/jorondo1/misc_scripts/main/myFunctions.sh
source scripts/myFunctions.sh

if [[ ! -f tmp/bac120_metadata_r220_short.tsv ]]; then
# Make a light version of:
curl -s -o tmp/bac120_metadata_r220.tsv.gz https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/bac120_metadata_r220.tsv.gz
zcat tmp/bac120_metadata_r220.tsv.gz | awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $2, $3}' > tmp/bac120_metadata_r220_short.tsv 
fi

cd $OUTDIR

#####################
### CHECKM
#####################

# Singularity
if [[ ! -f $OUTDIR/tmp/checkm2/quality_report.tsv ]]; then
	echo "Running checkm!"
	module load apptainer
	
$SINGULARITY/checkm2.1.0.2.sif \
	checkm2 predict --threads ${THREADS} \
	--database_path ${DB}/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd \
	--input ${MAG_DIR} --extension .fa --output-directory ${OUTDIR}/tmp/checkm2

	module unload apptainer
else echo 'checkM output found! Skipping.'
fi
cp $OUTDIR/tmp/checkm2/quality_report.tsv .

#####################
### skANI
#####################

# List MAGs
find ${MAG_DIR} -type f -name '*.fa' > ${OUTDIR}/tmp/MAG_list.txt

# Sketch GTDB genomes if required:
if [[ ! -d ${GTDB_SKANI} ]]; then
	echo 'sketching GTDB reference genomes...'
	cd ${DB}/GTDB
	find gtdb_genomes_reps_${GTDB_V}/ -name '*.fna.gz' > gtdb_files_${GTDB_V}.txt
	${SKANI} sketch -l gtdb_files_${GTDB_V}.txt -o ${GTDB_SKANI} -t ${THREADS}
	cd $OUTDIR
else echo 'Genome sketch found! Skipping.'
fi

# Compute ANI
if [[ ! -f tmp/ANI_results_raw.txt ]]; then
	echo 'Calculate ANI using skANI...'
	${SKANI} search -d ${GTDB_SKANI} -o tmp/ANI_results_raw.txt -t ${THREADS} --ql tmp/MAG_list.txt 
else echo 'skANI results found! Skipping.'
fi
# Format output: 
cat tmp/ANI_results_raw.txt | sed -e "s|${MAG_DIR}/||g" -e "s|gtdb_genomes_reps_${GTDB_V}/database/GC./.../.../.../||g" \
	-e "s/_genomic.fna.gz//g" | awk 'BEGIN {FS=OFS="\t"} {gsub(".fa", "", $2); print}' > ANI_results.txt

#####################
### identify_novel_MAGs.py
#####################

echo 'Identifying novel MAGs...'
module load python/3.11.5
python3 ${MAIN}/scripts/novel_MAGs.py -a $OUTDIR/ANI_results.txt -g ${MAIN}/tmp/bac120_metadata_r220_short.tsv \
	-m $OUTDIR/tmp/MAG_list.txt -c $OUTDIR/quality_report.tsv -o $OUTDIR
module unload

#####################
### Sketch & index nMAGs
#####################

### Sketch better MAGs (append to list with nMAGs)
better=($(tail -n +2 $OUTDIR/betterMAGs.txt | awk -F'\t' '{print $1}'))


# List the genomes to remove from the index :
# REPS_220=${DB}/GTDB/gtdb_genomes_reps_r220
# for genome in ${better[@]}; do
# 	line=$(grep $genome $REPS_220/genome_paths.tsv)
# if [ ! -z "$line" ]; then
# 	file=${REPS_220}/$(echo $line | awk '{print $2}')$(echo $line | awk '{print $1}')
# 	echo $file >> $OUTDIR/tmp/nMAG_list.txt
# fi
# done

ml apptainer

# Check if all MAGs have a signature sketched
missing_sig=()
while IFS= read -r fa; do 
	SIG="${SM_SK}/nMAGs/$(basename $fa)"
	if [ ! -e "$SIG".sig ]; then 
	missing_sig+=("$SIG")
	fi
done < tmp/nMAG_list.txt

# Sketch novel genomes 
if [ ${#missing_sig[@]} -eq 0 ]; then
	echo 'All MAGs have been sketched. Skipping.'
else
	echo 'Sketching novel genomes...'
	$SOURMASH sketch dna -p scaled=1000,k=31,abund \
	--name-from-first --from-file tmp/nMAG_list.txt \
	--output-dir ${SM_SK}/nMAGs
# it fails, check if there are empty returns at the end of nMAG_list.txt
#### Eventually, use 'branchwater multisketch'

# Rename novel genomes signatures (otherwise the whole name of the first contig is used) 
echo 'Renaming genome sketches...'
for file in $(find ${SM_SK}/nMAGs -type f -name '*.sig'); do
	new_name=$(basename $file)
	new_name=${new_name%.fa*}
	new_name=${new_name%.fna*}
	new_name=${new_name%_genomic*}
	$SOURMASH sig rename $file "$new_name" -o ${SM_SK}/nMAGs/${new_name}
done
fi

# Create an index 
if [[ ! -f ${SM_SK}/nMAGs_index.sbt.zip ]]; then
	echo 'Create index for genome sketches...'
	$SOURMASH index ${SM_SK}/nMAGs_index ${SM_SK}/nMAGs/*.sig
else echo 'nMAGs index found. Skipping.'
fi

module unload apptainer

#####################
### Gather Metagenomes
#####################

# Gather metagenomes :
num_csv=$(find sourmash/ -type f -name '*gather.csv' | wc | awk '{print $1}')
# We expect two gather files per sample (default db + custom db)
if [[ ${num_csv} -lt "$((2 * ${N_SAM}))" ]]; then

echo "Some output files are missing (${num_csv} found, $((2 * ${N_SAM})) expected.)"
# redo=$(grep -vnf <(find tmp/sourmash/ -type f -name '*gather.csv' -print0 | xargs -0 -I {} basename {} | sed 's/_.*//' | sort -u) ${SAM_LIST}/clean_samples.tsv| cut -d':' -f1 |  paste -sd,)

# save jobID
jobID=$(sbatch --array=1-"${N_SAM}" --export=ANCHOR,ILAFORES,DB,OUTDIR,MAGs_IDX,SAM_LIST,GTDB_V,SM_SK \
	$MAIN/scripts/gather_SLURM.sh | awk '{print $4}'); echo "Submitted job array with Job ID: $jobID"

# Periodically check if the job is still running, then whether all expected output files are there
while true; do
	# Check if job runs
    squeue -j "$jobID" > /dev/null 2>&1
    job_status=$?

    if [ $job_status -eq 0 ]; then
        echo "Job $jobID is still running."
		sleep 300
    else
        echo "Job $jobID has finished."
		break
    fi
done

else echo "All "${num_csv}" expected output files were found."
fi

echo "Summarising containment"
fix_gtdb sourmash # there's a comma problem that shifts the columns in gtdb taxonomy
eval_cont sourmash # compute sample containment and show overall stats

echo "Done !"





