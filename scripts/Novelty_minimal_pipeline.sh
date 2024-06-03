#!/usr/bin/env bash

#####################
### PARSING OPTIONS
#####################
function Usage {
    echo "Identify novel species-level MAGs and leverage them for community composition estimation."
	echo "Current directory should include scripts/novel_MAGs.py."
	echo "Assumes the GTDB reference genomes are under /home/def-ilafores/ref_dbs/GTDB/gtdb_genomes_reps_*"
	echo "Options:"
    echo "-o	Output directory name."
	echo "-g	GTDB database version (default: r220)"
	echo "-m	directory containing MAGs"
	echo "-s	directory containing clean samples"
    exit 1
}

while getopts 'o:g:m:s:' flag; do
    case "${flag}" in
		o) export OUTDIR=${PWD}/${OPTARG} ;;
		g) export GTDB_V=${OPTARG} ;;
		m) export MAG_DIR=${OPTARG} ;;
		s) export SAMPLE_DIR=${OPTARG} ;;
    esac
done

#####################
### SETUP
#####################

if [[ -z "$OUTDIR" ]] || [[ -z "$MAG_DIR" ]] || [[ -z "$SAMPLE_DIR" ]] ; then
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
export ILL_PIPELINES=${ILAFORES}/programs/ILL_pipelines/containers
export SINGULARITY="singularity exec --writable-tmpfs -e -B ${ILAFORES}:${ILAFORES} ${ILL_PIPELINES}/containers"
export SOURMASH="${SINGULARITY}/sourmash.4.7.0.sif sourmash"
export DB=${ILAFORES}/ref_dbs
export SM_SK=$OUTDIR/tmp/sourmash/sketches
export MAGs_IDX=$(find ${SM_SK} -type f -name 'nMAGs_index*')
export SKANI=${ILAFORES}/programs/skani/skani
export GTDB_SKANI=${DB}/GTDB/gtdb_skani_${GTDB_V}
export ANCHOR=/nfs3_ib/nfs-ip34
export N_SAM=$(wc ${SAMPLE_DIR}/clean_samples.tsv | awk '{print $1}')

if [[ ! -d $GTDB_SKANI ]] && [[ ! -d ${DB}/GTDB/gtdb_genomes_reps_${GTDB_V} ]]; then
	echo 'Reference database unavailable. Download and store under "${DB}/GTDB".'
	Usage
fi

# Setup project directories
mkdir -p scripts ${SM_SK}/nMAGs ${OUTDIR}/output ${OUTDIR}/tmp/logs ${OUTDIR}/tmp/checkm2

# gather output post-processing
if [[ -f scripts/myFunctions.sh ]]; then
	wget https://raw.githubusercontent.com/jorondo1/misc_scripts/main/myFunctions.sh -P scripts/
fi
source scripts/myFunctions.sh

cd $OUTDIR

# List MAGs
find "${MAG_DIR}" -type f -name '*.fa' > ${OUTDIR}/tmp/MAG_list.txt

exit 1
#####################
### CHECKM
#####################

# Singularity
if [[ -f ${OUTDIR}/tmp/checkm2/quality_report.tsv ]]; then
	echo "Running checkm!"
	module load apptainer
	
	"$SINGULARITY"/checkm2.1.0.2.sif \
	checkm2 predict --threads 48 \
	--database_path ${DB}/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd \
	--input ${MAG_DIR} --extension .fa --output-directory ${OUTDIR}/tmp/checkm2

	module unload apptainer
	cp tmp/checkm2/quality_report.tsv output/
else echo 'checkM output found! Skipping.'
fi

#####################
### skANI
#####################

# Sketch GTDB genomes if required:
if [[ ! -d ${GTDB_SKANI} ]]; then
	echo 'sketching GTDB reference genomes...'
	cd ${DB}/GTDB
	find gtdb_genomes_reps_${GTDB_V}/ -name '*.fna.gz' > gtdb_files_${GTDB_V}.txt
	${SKANI} sketch -l gtdb_files_${GTDB_V}.txt -o ${GTDB_SKANI} -t 72
	cd $OUTDIR
else echo 'Genome sketch found! Skipping.'
fi 

# Compute ANI
if [[ ! -f out/ANI_results.txt ]]; then
	echo 'Calculate ANI using skANI...'
	${SKANI} search -d ${GTDB_SKANI} -o tmp/ANI_results_raw.txt -t 72 --ql tmp/MAG_list.txt 
	# Format output: 
	cat tmp/ANI_results_raw.txt | sed -e "s|${MAG_DIR}/||g" -e "s|${GTDB_SKANI}/database/GC./.../.../.../||g" \
	-e "s/_genomic.fna.gz//g" | awk 'BEGIN {FS=OFS="\t"} {gsub(".fa", "", $2); print}' > output/ANI_results.txt
else echo 'skANI results found! Skipping.'
fi

#####################
### identify_novel_MAGs.py
#####################

echo 'Identifying novel MAGs...'
python3 ${MAIN}/scripts/novel_MAGs.py -a output/ANI_results.txt \
	-m tmp/MAG_list.txt -c output/quality_report.tsv -o tmp

#####################
### gather_SLURM.sh
#####################

ml apptainer

# Check if all MAGs have a signature sketched
missing_sig=()
while IFS= read -r fa; do 
	if [ -e "$(basename ${fa}).sig" ]; then 
	missing_sig+=("$fa")
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
# if it fails, check if there are empty returns at the end of nMAG_list.txt
#### Eventually, use 'branchwater multisketch'

# Rename novel genomes signatures (otherwise the whole name of the first contig is used) 
echo 'Renaming genome sketches...'
for file in $(find ${SM_SK}/nMAGs -type f -name '*.sig'); do
	new_name=$(basename $file)
	$SOURMASH sig rename $file "${new_name%.fa.sig}" -o ${SM_SK}/nMAGs/${new_name}
done
fi

# Create an index 
echo 'Create index for genome sketches...'
$SOURMASH index ${SM_SK}/nMAGs_index ${SM_SK}/nMAGs/*.sig
module unload apptainer

# Gather metagenomes ; savec jobID
jobID=$(sbatch --array=1-"${N_SAM}" --export=ANCHOR,ILAFORES,DB,OUTDIR,MAGs_IDX,SAMPLE_DIR,GTDB_V \
	$MAIN/scripts/gather_SLURM.sh | awk '{print $4}'); echo "Submitted job array with Job ID: $jobID"
sleep 600

# Periodically check if the job is still running, then whether all expected output files are there
while true; do
	# Check if job runs
    squeue -j "$jobID" > /dev/null 2>&1
    job_status=$?

    if [ $job_status -eq 0 ]; then
        echo "Job $jobID is still running."
    else
        echo "Job $jobID has finished."
		num_csv=$(find tmp/sourmash/ -type f -name '*gather.csv' | wc | awk '{print $1}')
		# We expect two gather files per sample (default db + custom db)
		if [[ ${num_csv} -ge "$((2 * ${N_SAM}))" ]]; then
			echo "All "${num_csv}" expected output files have been produced."
		else 
			echo 'Some output files are missing ("${num_csv}" found, "$((2 * ${N_SAM}))" expected.)'
			redo=$(grep -vnf <(find tmp/sourmash/ -type f -name '*gather.csv' -print0 | xargs -0 -I {} basename {} | sed 's/_.*//' | sort -u) ${SAMPLE_DIR}/clean_samples.tsv| cut -d':' -f1 |  paste -sd,)
			# TBC
		fi
    fi
    sleep 30
done

fix_gtdb tmp/sourmash # there's a comma problem that staggers the columns in gtdb taxonomy
eval_cont tmp/sourmash output # compute sample containment and show overall stats

#####################
### community_abundance.R
#####################

# Rscript my_script.R "${GTDB_V}"

# Containment plot
# Diversity plot

# Produce test stats ?? or just annotate diversity plots






